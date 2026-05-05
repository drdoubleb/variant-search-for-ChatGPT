// Serverless proxy for the CivicDB GraphQL API.
// The browser may be blocked by CORS when calling civicdb.org directly;
// this function runs server-side and forwards gene + assertion data.
//
// GET /api/civic?gene=BRAF&protein=V600E
//
// No API key required for public CivicDB read access.
// Always returns HTTP 200 — errors are reported in the response body so the
// frontend can degrade gracefully instead of showing a 502.

const CIVIC_GRAPHQL = 'https://civicdb.org/api/graphql';

// Post a GraphQL query and return the parsed JSON body.
// Returns { data, errors } — never throws (catches network + JSON errors).
async function civicPost(query) {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 12000);
    try {
        const res = await fetch(CIVIC_GRAPHQL, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
            body: JSON.stringify({ query }),
            signal: controller.signal
        });
        const text = await res.text();
        let body;
        try { body = JSON.parse(text); } catch { body = {}; }
        if (!res.ok) return { data: null, errors: [{ message: `HTTP ${res.status}: ${text.slice(0, 300)}` }] };
        return body;
    } catch (e) {
        return { data: null, errors: [{ message: e.message || String(e) }] };
    } finally {
        clearTimeout(timer);
    }
}

// Attempt several query shapes to accommodate CivicDB API schema variations.
// CivicDB v2 renamed entities ("genes" may be "features") — we try both.
async function findGene(safeGene) {
    // Attempt 1 — canonical shape used in CivicDB v2 docs
    const r1 = await civicPost(`{
        genes(name: "${safeGene}") {
            nodes {
                id name description
                variants { nodes { id name variantTypes { nodes { name } } } }
            }
        }
    }`);
    const g1 = r1?.data?.genes?.nodes?.[0];
    if (g1) return { gene: g1, queryErrors: null };

    // Attempt 2 — "features" schema (CivicDB may use this name)
    const r2 = await civicPost(`{
        features(name: "${safeGene}") {
            nodes {
                id name
                ... on Gene {
                    description
                    variants { nodes { id name variantTypes { nodes { name } } } }
                }
            }
        }
    }`);
    const g2 = r2?.data?.features?.nodes?.[0];
    if (g2) return { gene: g2, queryErrors: null };

    // Collect errors for debugging
    const errs = [
        ...(r1?.errors || []).map(e => `genes: ${e.message}`),
        ...(r2?.errors || []).map(e => `features: ${e.message}`)
    ];
    return { gene: null, queryErrors: errs.length ? errs : ['No results'] };
}

async function findAssertions(geneId) {
    // Attempt 1 — geneIds filter
    const r1 = await civicPost(`{
        assertions(geneIds: [${geneId}], status: ACCEPTED) {
            nodes {
                id ampLevel clinicalSignificance significance summary
                disease { name }
                therapies { nodes { name } }
            }
        }
    }`);
    if (!r1?.errors?.length && r1?.data?.assertions) {
        return r1.data.assertions.nodes || [];
    }

    // Attempt 2 — featureIds filter
    const r2 = await civicPost(`{
        assertions(featureIds: [${geneId}], status: ACCEPTED) {
            nodes {
                id ampLevel clinicalSignificance significance summary
                disease { name }
                therapies { nodes { name } }
            }
        }
    }`);
    if (!r2?.errors?.length && r2?.data?.assertions) {
        return r2.data.assertions.nodes || [];
    }

    return [];
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });

    const { gene, protein } = req.query;
    if (!gene || !String(gene).trim()) {
        return res.status(400).json({ error: 'Missing required query parameter: gene' });
    }

    const safeGene = String(gene).replace(/[^A-Za-z0-9\-_.]/g, '');
    if (!safeGene) return res.status(400).json({ error: 'Invalid gene name' });

    const { gene: apiGene, queryErrors } = await findGene(safeGene);

    if (!apiGene) {
        console.warn('CivicDB gene lookup failed for', safeGene, queryErrors);
        return res.status(200).json({ gene: null, matchedVariant: null, assertions: [], queryErrors });
    }

    // Match variant by protein change (normalised comparison, supports triple-letter)
    let matchedVariant = null;
    const normInput = String(protein || '').replace(/^p\./i, '').toLowerCase().replace(/[^a-z0-9*_]/g, '');
    if (normInput && Array.isArray(apiGene.variants?.nodes)) {
        for (const v of apiGene.variants.nodes) {
            const vn = String(v.name || '').toLowerCase().replace(/[^a-z0-9*_]/g, '');
            if (vn && (vn === normInput || normInput.includes(vn) || vn.includes(normInput))) {
                matchedVariant = v;
                break;
            }
        }
    }

    const assertions = apiGene.id ? await findAssertions(apiGene.id) : [];

    return res.status(200).json({ gene: apiGene, matchedVariant, assertions });
}
