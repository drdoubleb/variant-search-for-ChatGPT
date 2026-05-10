// Serverless proxy for fetching ClinVar esummary for a specific variation ID.
// Returns germline, somatic, and oncogenicity classifications plus per-condition data.
//
// GET /api/clinvar-variant?id=375941
//
// Optional env var: NCBI_API_KEY — increases rate limit from 3 to 10 req/s.

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

async function ncbiFetch(url) {
    const res = await fetch(url, {
        headers: { 'Accept': 'application/json' },
        signal: AbortSignal.timeout(12000)
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

    const { id } = req.query;
    if (!id || !/^\d+$/.test(String(id))) {
        return res.status(400).json({ error: 'Missing or invalid id parameter (must be numeric ClinVar variation ID)' });
    }

    const apiKey = process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';

    try {
        const sumUrl = `${EUTILS_BASE}/esummary.fcgi?db=clinvar&retmode=json&id=${encodeURIComponent(id)}${apiKey}`;
        const sumData = await ncbiFetch(sumUrl);
        const rec = sumData?.result?.[id] || {};

        // Extract per-condition somatic data if present.
        // The esummary clinical_assertion_summary field may list individual SCVs with condition + impact info.
        const somaticConditions = [];
        const assertions = rec.clinical_assertion_summary || [];
        if (Array.isArray(assertions)) {
            assertions.forEach((a) => {
                if (a.submission_type === 'somatic' || a.classification_type === 'somatic_clinical_impact') {
                    somaticConditions.push({
                        condition: a.trait_name || a.condition || '',
                        impact: a.clinical_impact_category || a.somatic_clinical_impact || a.description || '',
                        submitter: a.submitter_name || '',
                        lastEvaluated: a.last_evaluated || '',
                        accession: a.accession || '',
                    });
                }
            });
        }

        return res.status(200).json({
            title: rec.title || '',
            variationName: rec.variation_set?.[0]?.variation_name || '',
            gene: rec.gene_sort || '',
            germline: rec.germline_classification || null,
            somatic: rec.somatic_clinical_impact || null,
            oncogenicity: rec.oncogenicity_classification || null,
            somaticConditions,
        });
    } catch (err) {
        console.error('ClinVar variant proxy error:', err);
        return res.status(502).json({ error: 'ClinVar variant lookup failed', detail: err.message });
    }
}
