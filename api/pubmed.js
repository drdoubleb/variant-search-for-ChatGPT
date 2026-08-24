// Serverless proxy for NCBI E-utilities PubMed search.
// The browser cannot call eutils.ncbi.nlm.nih.gov directly due to CORS;
// this function runs server-side and forwards the response.
//
// GET /api/pubmed?term=BRAF+V600E&limit=5
//
// Optional env var: NCBI_API_KEY — increases rate limit from 3 to 10 req/s.

import { ncbiFetchJson, ncbiFetchResponse } from './_ncbi.js';

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

async function ncbiFetch(url) {
    return ncbiFetchJson(url);
}

async function ncbiFetchText(url) {
    const res = await ncbiFetchResponse(url, { accept: 'text/xml,application/xml' });
    return res.text();
}

function parseAbstractsFromXml(xmlText) {
    const abstracts = {};
    const articles = xmlText.split(/<PubmedArticle[\s>]/);
    for (const art of articles) {
        const pmidMatch = art.match(/<PMID[^>]*>(\d+)<\/PMID>/);
        if (!pmidMatch) continue;
        const pmid = pmidMatch[1];
        const textParts = [];
        const regex = /<AbstractText[^>]*>([\s\S]*?)<\/AbstractText>/g;
        let m;
        while ((m = regex.exec(art)) !== null) {
            const text = m[1].replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim();
            if (text) textParts.push(text);
        }
        if (textParts.length > 0) {
            abstracts[pmid] = textParts.join(' ');
        }
    }
    return abstracts;
}

import { rejectDisallowedOrigin } from './_origin.js';
import { setEdgeCache, setNoStore } from './_cache.js';

// ── LitVar2 / PubTator3 variant-normalized literature ──────────────────────
//
// A free-text PubMed search for "BRAF V600E" misses papers that write the same
// variant as p.Val600Glu, c.1799T>A, or rs113488022. NCBI's LitVar2 (served by
// the PubTator3 API) normalizes variant mentions across nomenclatures:
//   1. entity/autocomplete resolves "BRAF V600E" to a variant entity id
//      (e.g. "@VARIANT_p.V600E_BRAF_human", backed by litvar/rsID), then
//   2. search/?text=<entity id> returns relevance-RANKED articles with
//      metadata included — optionally refined with a free-text AND clause
//      (used for the tumor-type tab).
// Abstracts still come from efetch, reusing the retrying eutils helper.
const PUBTATOR_BASE = 'https://www.ncbi.nlm.nih.gov/research/pubtator3-api';

function formatAuthors(list) {
    const authList = Array.isArray(list) ? list.filter(Boolean) : [];
    const firstTwo = authList.slice(0, 2);
    return firstTwo.length > 0
        ? firstTwo.join(', ') + (authList.length > 2 ? ' et al.' : '')
        : '';
}

// Returns { total, articles, litvar } or null when no variant entity matched
// (caller falls back to the plain term search). Throws on upstream failure —
// the caller treats that as "fall back" too.
async function pubtatorVariantSearch({ gene, variant, extra, maxResults, apiKey }) {
    const acQuery = [gene, variant].filter(Boolean).join(' ');
    const acUrl = `${PUBTATOR_BASE}/entity/autocomplete/?query=${encodeURIComponent(acQuery)}`;
    const entities = await ncbiFetch(acUrl);
    const geneUpper = String(gene || '').trim().toUpperCase();
    const entity = (Array.isArray(entities) ? entities : []).find((e) => e && e.biotype === 'variant'
        && String(e._id || '').startsWith('@')
        && (!geneUpper || String(e.description || '').toUpperCase().includes(geneUpper)));
    if (!entity) return null;

    const text = extra && String(extra).trim()
        ? `${entity._id} AND "${String(extra).trim().replace(/"/g, '')}"`
        : entity._id;
    const searchData = await ncbiFetch(`${PUBTATOR_BASE}/search/?text=${encodeURIComponent(text)}`);
    const results = (Array.isArray(searchData?.results) ? searchData.results : []).slice(0, maxResults);

    // rsID and a LitVar2 deep link, when the entity is litvar-backed.
    const dbId = String(entity.db_id || '');
    const rsid = /^rs\d+/.test(dbId) ? dbId.replace(/#/g, '') : '';
    // The entity id carries the precise change ("@VARIANT_p.V600E_BRAF_human");
    // entity.name can be a grouped label ("p.V600X") — prefer the precise one.
    const idName = (String(entity._id).match(/^@VARIANT_([^_]+)/) || [])[1];
    const litvar = {
        id: entity._id,
        name: idName || entity.name || String(variant),
        rsid,
        url: entity.db === 'litvar' && dbId
            ? `https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?variant=${encodeURIComponent(`litvar@${dbId}`)}`
            : ''
    };

    const ids = results.map((r) => String(r.pmid)).filter((id) => /^\d+$/.test(id));
    let abstracts = {};
    if (ids.length > 0) {
        try {
            const fetchUrl = `${EUTILS_BASE}/efetch.fcgi?db=pubmed&rettype=xml&retmode=xml&id=${ids.join(',')}${apiKey}`;
            abstracts = parseAbstractsFromXml(await ncbiFetchText(fetchUrl));
        } catch (e) {
            console.warn('PubMed abstract fetch failed (litvar path):', e.message);
        }
    }

    const articles = results.map((r) => ({
        pmid: String(r.pmid),
        title: r.title || '',
        authors: formatAuthors(r.authors),
        journal: r.journal || '',
        year: r.date ? String(r.date).slice(0, 4) : '',
        abstract: abstracts[String(r.pmid)] || ''
    }));
    return { total: Number(searchData?.count) || articles.length, articles, litvar };
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });
    // Block third-party websites from using this deployment as their backend
    // (no-Origin callers pass — see api/_origin.js).
    if (rejectDisallowedOrigin(req, res)) return;

    const { term, limit = '5', gene = '', variant = '', extra = '' } = req.query;
    if ((!term || !String(term).trim()) && !String(variant).trim()) {
        return res.status(400).json({ error: 'Missing required query parameter: term (or variant)' });
    }

    const apiKey = process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';
    const maxResults = Math.min(Math.max(1, parseInt(limit, 10) || 5), 20);

    // Successful lookups are shared across users — let Vercel's CDN serve
    // repeats without re-invoking the function or hitting the upstream.
    setEdgeCache(res, 21600);

    // Variant-normalized path: when the caller names a specific variant, prefer
    // LitVar2/PubTator3 (catches every nomenclature spelling, relevance-ranked).
    // No matching entity, or an upstream failure, falls back to the term search.
    if (String(variant).trim()) {
        try {
            const lv = await pubtatorVariantSearch({ gene, variant, extra, maxResults, apiKey });
            if (lv) return res.status(200).json({ ...lv, backend: 'litvar' });
        } catch (e) {
            console.warn('PubTator3 variant search failed; falling back to term search:', e.message);
        }
    }
    const searchTerm = (term && String(term).trim()) || [gene, variant, extra].filter(Boolean).join(' ');

    try {
        // Step 1: esearch — get PubMed IDs ranked by relevance
        const searchUrl = `${EUTILS_BASE}/esearch.fcgi?db=pubmed&retmode=json&retmax=${maxResults}&sort=relevance&term=${encodeURIComponent(searchTerm)}${apiKey}`;
        const searchData = await ncbiFetch(searchUrl);
        const ids = searchData?.esearchresult?.idlist || [];
        const total = parseInt(searchData?.esearchresult?.count || '0', 10);

        if (ids.length === 0) {
            return res.status(200).json({ total, articles: [], backend: 'term' });
        }

        // Step 2: esummary — get article metadata for each ID
        const sumUrl = `${EUTILS_BASE}/esummary.fcgi?db=pubmed&retmode=json&id=${ids.join(',')}${apiKey}`;
        const sumData = await ncbiFetch(sumUrl);

        // Step 3: efetch — get article abstracts (XML)
        let abstracts = {};
        try {
            const fetchUrl = `${EUTILS_BASE}/efetch.fcgi?db=pubmed&rettype=xml&retmode=xml&id=${ids.join(',')}${apiKey}`;
            const fetchXml = await ncbiFetchText(fetchUrl);
            abstracts = parseAbstractsFromXml(fetchXml);
        } catch (e) {
            console.warn('PubMed abstract fetch failed:', e.message);
        }

        const articles = ids.map((id) => {
            const rec = sumData?.result?.[id] || {};
            const authList = Array.isArray(rec.authors) ? rec.authors : [];
            const firstTwo = authList.slice(0, 2).map((a) => a.name || '').filter(Boolean);
            const authors = firstTwo.length > 0
                ? firstTwo.join(', ') + (authList.length > 2 ? ' et al.' : '')
                : '';
            const year = rec.pubdate ? String(rec.pubdate).slice(0, 4) : '';
            return {
                pmid: id,
                title: rec.title || '',
                authors,
                journal: rec.source || '',
                year,
                abstract: abstracts[id] || ''
            };
        });

        return res.status(200).json({ total, articles, backend: 'term' });
    } catch (err) {
        setNoStore(res); // transient failure — never pin it in the CDN
        console.error('PubMed proxy error:', err);
        return res.status(502).json({ error: 'PubMed lookup failed', detail: err.message });
    }
}
