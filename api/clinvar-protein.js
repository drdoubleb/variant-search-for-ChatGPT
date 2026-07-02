// Serverless proxy for NCBI ClinVar protein-change queries.
//
// Returns ClinVar variants that share a given amino-acid change, regardless of
// the underlying nucleotide change. This mirrors what the ClinVar website does
// for a search like "CTNNB1 G34R" (which returns both c.100G>C and c.100G>A) and
// catches same-consequence variants that an exact genomic / ±N bp window match
// misses — including ones annotated on a different transcript.
//
// GET /api/clinvar-protein?gene=CTNNB1&change=Gly34Arg
//   - change: the protein change to search for. One or more comma-separated
//     forms may be supplied (e.g. "Gly34Arg,G34R"); they are OR-ed in the query.
//
// The caller is expected to filter the returned variants down to an exact
// protein-change match (the free-text search is recall-oriented); titles are
// returned so that filtering can parse the canonical HGVS p. change.
//
// Optional env var: NCBI_API_KEY — increases rate limit from 3 to 10 req/s.

import { ncbiFetchJson } from './_ncbi.js';

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });

    const { gene, change, retmax: retmaxParam = '200' } = req.query;
    if (!gene || !change) {
        return res.status(400).json({ error: 'Missing required query parameters: gene, change' });
    }

    const safeGene = String(gene).replace(/[^A-Za-z0-9\-_.]/g, '');
    if (!safeGene) {
        return res.status(400).json({ error: 'Invalid gene value' });
    }

    // Accept one or more comma-separated change forms; sanitise each to the
    // alphanumeric/`*` characters that appear in HGVS protein bodies.
    const changeForms = String(change)
        .split(',')
        .map((s) => s.replace(/[^A-Za-z0-9*]/g, ''))
        .filter(Boolean);
    if (changeForms.length === 0) {
        return res.status(400).json({ error: 'Invalid change value' });
    }

    const retmax = Math.min(Math.max(1, parseInt(retmaxParam, 10) || 200), 500);
    const apiKey = process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';

    try {
        // ClinVar's own search help recommends querying same-protein-change
        // variants by gene symbol plus a one- or three-letter protein change
        // (e.g. "BRCA1 S803R"). Keep the proxy to that single website-style
        // query instead of broadening into multiple eUtils searches.
        const preferredChange = changeForms.find((form) => /^[A-Za-z]\d+[A-Za-z*]$/.test(form)) || changeForms[0];
        const term = `${safeGene} ${preferredChange}`;
        const searchUrl = `${EUTILS_BASE}/esearch.fcgi?db=clinvar&retmode=json&retmax=${retmax}&term=${encodeURIComponent(term)}${apiKey}`;
        const searchData = await ncbiFetchJson(searchUrl);
        const ids = searchData?.esearchresult?.idlist || [];
        const totalInClinVar = Number(searchData?.esearchresult?.count || 0);

        if (ids.length === 0) {
            return res.status(200).json({ variants: [], total: 0 });
        }

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
                variationName: rec.variation_set?.[0]?.variation_name || '',
                molecularConsequence: rec.molecular_consequence_list || rec.molecular_consequence || rec.variation_set?.[0]?.molecular_consequence || '',
                pos: varPos !== null ? Number(varPos) : null
            };
        });

        return res.status(200).json({ variants, total: totalInClinVar });
    } catch (err) {
        console.error('ClinVar protein proxy error:', err);
        return res.status(502).json({ error: 'ClinVar protein lookup failed', detail: err.message });
    }
}
