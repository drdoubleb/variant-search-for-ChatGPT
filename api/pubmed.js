// Serverless proxy for NCBI E-utilities PubMed search.
// The browser cannot call eutils.ncbi.nlm.nih.gov directly due to CORS;
// this function runs server-side and forwards the response.
//
// GET /api/pubmed?term=BRAF+V600E&limit=5
//
// Optional env var: NCBI_API_KEY — increases rate limit from 3 to 10 req/s.

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

async function ncbiFetch(url) {
    const res = await fetch(url, {
        headers: { 'Accept': 'application/json' },
        signal: AbortSignal.timeout(10000)
    });
    if (!res.ok) throw new Error(`NCBI request failed: ${res.status}`);
    return res.json();
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });

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
                year
            };
        });

        return res.status(200).json({ total, articles });
    } catch (err) {
        console.error('PubMed proxy error:', err);
        return res.status(502).json({ error: 'PubMed lookup failed', detail: err.message });
    }
}
