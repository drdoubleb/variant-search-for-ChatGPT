// Serverless proxy for NCBI ClinVar region queries.
// The browser cannot attach an API key to direct eutils calls safely;
// this function runs server-side so the key stays out of the client bundle.
//
// GET /api/clinvar-region?chrom=7&pos=140453136&window=10
//
// Optional env var: NCBI_API_KEY — increases rate limit from 3 to 10 req/s.

import { ncbiFetchJson } from './_ncbi.js';

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

import { rejectDisallowedOrigin } from './_origin.js';
import { setEdgeCache, setNoStore } from './_cache.js';

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });
    // Block third-party websites from using this deployment as their backend
    // (no-Origin callers pass — see api/_origin.js).
    if (rejectDisallowedOrigin(req, res)) return;

    const { chrom, pos, window: win = '10' } = req.query;
    if (!chrom || !pos) {
        return res.status(400).json({ error: 'Missing required query parameters: chrom, pos' });
    }

    const c = String(chrom).replace(/^chr/i, '').toUpperCase();
    const p = Number(pos);
    if (!c || !Number.isFinite(p)) {
        return res.status(400).json({ error: 'Invalid chrom or pos value' });
    }

    const windowSize = Math.min(Math.max(1, parseInt(win, 10) || 10), 500);
    const start = Math.max(1, p - windowSize);
    const end = p + windowSize;
    const apiKey = process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';

    // Successful lookups are shared across users — let Vercel's CDN serve
    // repeats without re-invoking the function or hitting the upstream.
    setEdgeCache(res, 21600);
    try {
        const term = `${c}[Chromosome] AND ${start}:${end}[Base Position for Assembly GRCh37]`;
        const searchUrl = `${EUTILS_BASE}/esearch.fcgi?db=clinvar&retmode=json&retmax=500&term=${encodeURIComponent(term)}${apiKey}`;
        const searchData = await ncbiFetchJson(searchUrl);
        const ids = searchData?.esearchresult?.idlist || [];
        const totalInClinVar = Number(searchData?.esearchresult?.count || 0);

        if (!Array.isArray(ids) || ids.length === 0) {
            return res.status(200).json({ variants: [], total: 0 });
        }

        // Fetch summaries in batches of 200 to stay within URL limits and
        // supply an explicit retmax (NCBI esummary defaults to 20 without WebEnv).
        const BATCH = 200;
        const resultMap = {};
        for (let i = 0; i < ids.length; i += BATCH) {
            const chunk = ids.slice(i, i + BATCH);
            const sumUrl = `${EUTILS_BASE}/esummary.fcgi?db=clinvar&retmode=json&retmax=${chunk.length}&id=${chunk.join(',')}${apiKey}`;
            const sumData = await ncbiFetchJson(sumUrl);
            Object.assign(resultMap, sumData?.result || {});
        }

        const variants = ids.map((id) => {
            const rec = resultMap[id] || {};
            const varLoc = rec.variation_set?.[0]?.variation_loc?.find?.((l) => String(l.assembly_name || '').toLowerCase().includes('grch37'))
                || rec.location?.find?.((l) => String(l.assembly || '').toLowerCase().includes('grch37'))
                || rec.variation_set?.[0]?.variation_loc?.[0]
                || rec.location?.[0]
                || null;
            const varPos = varLoc ? (varLoc.display_start || varLoc.start || varLoc.chr_start || null) : null;
            return {
                id,
                title: rec.title || '',
                germline: rec.germline_classification?.description || '',
                review: rec.germline_classification?.review_status || '',
                variationId: rec.variation_set?.[0]?.variation_xrefs?.find?.((x) => String(x.db || '').toLowerCase() === 'dbsnp')?.id || rec.variation_set?.[0]?.variation_name || '',
                variationName: rec.variation_set?.[0]?.variation_name || '',
                molecularConsequence: rec.molecular_consequence_list || rec.molecular_consequence || rec.variation_set?.[0]?.molecular_consequence || '',
                pos: varPos !== null ? Number(varPos) : null
            };
        });

        return res.status(200).json({ variants, total: totalInClinVar });
    } catch (err) {
        setNoStore(res); // transient failure — never pin it in the CDN
        console.error('ClinVar region proxy error:', err);
        return res.status(502).json({ error: 'ClinVar region lookup failed', detail: err.message });
    }
}
