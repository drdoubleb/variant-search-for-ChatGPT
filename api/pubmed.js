// Serverless proxy for NCBI E-utilities PubMed search.
// The browser cannot call eutils.ncbi.nlm.nih.gov directly due to CORS;
// this function runs server-side and forwards the response.
//
// GET /api/pubmed?term=BRAF+V600E&limit=5 or /api/pubmed?pmids=12826609,30224644
//
// Optional env var: NCBI_API_KEY — increases rate limit from 3 to 10 req/s.

import { ncbiFetchJson, ncbiFetchResponse } from './_ncbi.js';
import { boundedQueryText, positiveQueryInteger } from './_input.js';
import { attribute, child, children, childText, parseXml, text } from './_eutils-xml.js';
import { rejectDisallowedOrigin } from './_origin.js';
import { setEdgeCache, setNoStore } from './_cache.js';

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
const MAX_EXPLICIT_PMIDS = 20;
const PMID_PATTERN = /^[1-9]\d{0,8}$/;

async function ncbiFetch(url) {
    return ncbiFetchJson(url);
}

function parseExplicitPmids(input) {
    const query = boundedQueryText(input, MAX_EXPLICIT_PMIDS * 10 - 1);
    if (!query) return null;
    const ids = query.split(',').map(id => id.trim());
    if (ids.length > MAX_EXPLICIT_PMIDS || ids.some(id => !PMID_PATTERN.test(id))) return null;
    return [...new Set(ids)];
}

export function parsePubmedRecords(xml) {
    const root = parseXml(xml, { allowPubMedDoctype: true });
    if (root.name !== 'PubmedArticleSet' || child(root, 'ERROR') || child(root, 'Error')) {
        throw new Error('PubMed did not return an article set');
    }
    const records = new Map();
    for (const node of root.children) {
        const citation = child(node, 'MedlineCitation') || child(node, 'BookDocument');
        const pmid = childText(citation, 'PMID');
        if (!PMID_PATTERN.test(pmid)) continue;
        const article = child(citation, 'Article') || citation;
        const journal = child(article, 'Journal');
        const book = child(article, 'Book');
        const pubDate = child(child(journal, 'JournalIssue'), 'PubDate') || child(book, 'PubDate');
        const year = childText(pubDate, 'Year') || childText(pubDate, 'MedlineDate').match(/^\d{4}/)?.[0] || '';
        const authorList = child(article, 'AuthorList');
        const authors = children(authorList, 'Author').map(author =>
            childText(author, 'CollectiveName') || [childText(author, 'LastName'), childText(author, 'Initials')].filter(Boolean).join(' '));
        const paragraphs = children(child(article, 'Abstract'), 'AbstractText').flatMap(paragraph => {
            const content = text(paragraph).replace(/\s+/g, ' ');
            const label = attribute(paragraph, 'Label');
            return content ? [label ? `${label}: ${content}` : content] : [];
        });
        records.set(pmid, {
            pmid,
            title: childText(article, 'ArticleTitle') || childText(book, 'BookTitle'),
            authors: formatAuthors(authors),
            journal: childText(journal, 'ISOAbbreviation') || childText(journal, 'Title') || childText(book, 'BookTitle'),
            year,
            abstract: paragraphs.join('\n')
        });
    }
    return records;
}

async function fetchPubmedRecords(ids, apiKey) {
    const url = `${EUTILS_BASE}/efetch.fcgi?db=pubmed&retmode=xml&id=${ids.join(',')}${apiKey}`;
    const res = await ncbiFetchResponse(url, { accept: 'text/xml,application/xml' });
    return parsePubmedRecords(await res.text());
}

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
    let records = new Map();
    let abstractsUnavailable = false;
    if (ids.length > 0) {
        try {
            records = await fetchPubmedRecords(ids, apiKey);
        } catch (e) {
            abstractsUnavailable = true;
            console.warn('PubMed abstract fetch failed (litvar path):', e.message);
        }
    }

    const articles = results.map((r) => ({
        pmid: String(r.pmid),
        title: r.title || '',
        authors: formatAuthors(r.authors),
        journal: r.journal || '',
        year: r.date ? String(r.date).slice(0, 4) : '',
        abstract: records.get(String(r.pmid))?.abstract || ''
    }));
    return { total: Number(searchData?.count) || articles.length, articles, litvar, abstractsUnavailable };
}

// PubMed does not fail on a phrase it has never indexed: it drops the phrase and
// answers the rest of the query. "BRCA2 D362fs" therefore came back as the most
// relevant BRCA2 papers in general — listed under "Gene + Variant" and handed to
// the AI review as literature on the variant. esearch names what it dropped
// (errorlist.phrasesnotfound); when that is the whole variant, PubMed has
// nothing on the variant, and saying so is the answer.
//
// Part of a variant surviving is still a search for it ("p.Glu746_Ala750del"
// loses "p.Glu746" and matches on "Ala750del"), so only a variant with nothing
// left of it counts as dropped — nothing that says which variant, that is. A
// variant is told apart by its position: "c.1026+59T>G" loses "c.1026+59T", and
// what PubMed then searches is PTEN AND "G" — 1,491 papers; "c.489+1G>T" is
// CFTR AND "T", 8,572. A stray base is not the variant, nor is a "del".
export function variantPhraseDropped(variant, phrasesNotFound) {
    const whole = String(variant || '').toLowerCase();
    if (!/[a-z0-9]/.test(whole)) return false;
    let rest = whole;
    for (const phrase of Array.isArray(phrasesNotFound) ? phrasesNotFound : []) {
        const dropped = String(phrase || '').toLowerCase();
        if (dropped) rest = rest.split(dropped).join(' ');
    }
    if (rest === whole) return false; // nothing of the variant was dropped
    return /\d/.test(whole) ? !/\d/.test(rest) : !/[a-z0-9]/.test(rest);
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

    // Bound free-text inputs to prevent oversized public requests to NCBI.
    const query = req.query || {};
    const term = query.term === undefined ? '' : boundedQueryText(query.term, 512);
    const gene = query.gene === undefined ? '' : boundedQueryText(query.gene, 64);
    const variant = query.variant === undefined ? '' : boundedQueryText(query.variant, 128);
    const extra = query.extra === undefined ? '' : boundedQueryText(query.extra, 256);
    const limit = positiveQueryInteger(query.limit === undefined ? '5' : query.limit, 100000);
    if (query.pmids !== undefined) {
        const ids = parseExplicitPmids(query.pmids);
        if (!ids || ['term', 'gene', 'variant', 'extra'].some(key => query[key] !== undefined)) {
            return res.status(400).json({ error: `Supply 1 to ${MAX_EXPLICIT_PMIDS} comma-separated numeric PMIDs without search parameters` });
        }
        try {
            const apiKey = process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';
            const records = await fetchPubmedRecords(ids, apiKey);
            const articles = ids.filter(id => records.has(id)).map(id => records.get(id));
            const missingPmids = ids.filter(id => !records.has(id));
            if (missingPmids.length) setNoStore(res);
            else setEdgeCache(res, 21600);
            return res.status(200).json({ backend: 'pmid', total: articles.length, articles, missingPmids });
        } catch (err) {
            setNoStore(res);
            console.warn('PubMed PMID lookup failed:', err.message);
            return res.status(502).json({ error: 'PubMed PMID lookup failed' });
        }
    }
    if (term === null || gene === null || variant === null || extra === null || limit === null) {
        return res.status(400).json({ error: 'Invalid or oversized PubMed query parameter' });
    }
    if (!term && !variant) {
        return res.status(400).json({ error: 'Missing required query parameter: term (or variant)' });
    }

    const apiKey = process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';
    const maxResults = Math.min(limit, 20);

    // Successful lookups are shared across users — let Vercel's CDN serve
    // repeats without re-invoking the function or hitting the upstream.
    setEdgeCache(res, 21600);

    // Variant-normalized path: when the caller names a specific variant, prefer
    // LitVar2/PubTator3 (catches every nomenclature spelling, relevance-ranked).
    // No matching entity, or an upstream failure, falls back to the term search.
    if (String(variant).trim()) {
        try {
            const lv = await pubtatorVariantSearch({ gene, variant, extra, maxResults, apiKey });
            if (lv) {
                if (lv.abstractsUnavailable) setNoStore(res);
                return res.status(200).json({ ...lv, backend: 'litvar' });
            }
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

        const notFound = searchData?.esearchresult?.errorlist?.phrasesnotfound;
        if (variantPhraseDropped(variant, notFound)) {
            return res.status(200).json({ total: 0, articles: [], backend: 'term', notFound });
        }

        if (ids.length === 0) {
            return res.status(200).json({ total, articles: [], backend: 'term' });
        }

        // Step 2: esummary — get article metadata for each ID
        const sumUrl = `${EUTILS_BASE}/esummary.fcgi?db=pubmed&retmode=json&id=${ids.join(',')}${apiKey}`;
        const sumData = await ncbiFetch(sumUrl);

        // Step 3: efetch — get article abstracts (XML)
        let records = new Map();
        let abstractsUnavailable = false;
        try {
            records = await fetchPubmedRecords(ids, apiKey);
        } catch (e) {
            abstractsUnavailable = true;
            setNoStore(res);
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
                abstract: records.get(String(id))?.abstract || ''
            };
        });

        return res.status(200).json({ total, articles, backend: 'term', abstractsUnavailable });
    } catch (err) {
        setNoStore(res); // transient failure — never pin it in the CDN
        console.error('PubMed proxy error:', err);
        return res.status(502).json({ error: 'PubMed lookup failed', detail: err.message });
    }
}
