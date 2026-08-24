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

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });
    // Block third-party websites from using this deployment as their backend
    // (no-Origin callers pass — see api/_origin.js).
    if (rejectDisallowedOrigin(req, res)) return;

    const { term, limit = '5' } = req.query;
    if (!term || !String(term).trim()) {
        return res.status(400).json({ error: 'Missing required query parameter: term' });
    }

    const apiKey = process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';
    const maxResults = Math.min(Math.max(1, parseInt(limit, 10) || 5), 20);

    try {
        // Step 1: esearch — get PubMed IDs ranked by relevance
        const searchUrl = `${EUTILS_BASE}/esearch.fcgi?db=pubmed&retmode=json&retmax=${maxResults}&sort=relevance&term=${encodeURIComponent(term)}${apiKey}`;
        const searchData = await ncbiFetch(searchUrl);
        const ids = searchData?.esearchresult?.idlist || [];
        const total = parseInt(searchData?.esearchresult?.count || '0', 10);

        if (ids.length === 0) {
            return res.status(200).json({ total, articles: [] });
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

        return res.status(200).json({ total, articles });
    } catch (err) {
        console.error('PubMed proxy error:', err);
        return res.status(502).json({ error: 'PubMed lookup failed', detail: err.message });
    }
}
