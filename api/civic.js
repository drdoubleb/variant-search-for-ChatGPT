// Serverless proxy for the CivicDB GraphQL API.
// The browser may be blocked by CORS when calling civicdb.org directly;
// this function runs server-side and forwards gene + assertion data.
//
// GET /api/civic?gene=BRAF&protein=V600E
//
// No API key is required for public read access to CivicDB.

const CIVIC_GRAPHQL = 'https://civicdb.org/api/graphql';

async function civicQuery(query) {
    const res = await fetch(CIVIC_GRAPHQL, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
        body: JSON.stringify({ query }),
        signal: AbortSignal.timeout(12000)
    });
    if (!res.ok) throw new Error(`CivicDB request failed: ${res.status}`);
    const data = await res.json();
    if (data?.errors?.length) throw new Error(`CivicDB GraphQL error: ${data.errors[0]?.message}`);
    return data;
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

    try {
        // Step 1: get gene info + variant list
        const geneData = await civicQuery(`{
            genes(name: "${safeGene}") {
                nodes {
                    id name description
                    variants {
                        nodes {
                            id name
                            variantTypes { nodes { name } }
                        }
                    }
                }
            }
        }`);

        const apiGene = geneData?.data?.genes?.nodes?.[0] || null;
        if (!apiGene) {
            return res.status(200).json({ gene: null, matchedVariant: null, assertions: [] });
        }

        // Match variant by protein change (normalised comparison)
        let matchedVariant = null;
        if (protein && Array.isArray(apiGene.variants?.nodes)) {
            const normProt = String(protein).replace(/^p\./i, '').toLowerCase().replace(/[^a-z0-9*_]/g, '');
            for (const v of apiGene.variants.nodes) {
                const vn = String(v.name || '').toLowerCase().replace(/[^a-z0-9*_]/g, '');
                if (vn && normProt && (vn === normProt || normProt.includes(vn) || vn.includes(normProt))) {
                    matchedVariant = v;
                    break;
                }
            }
        }

        // Step 2: get accepted assertions for this gene
        let assertions = [];
        if (apiGene.id) {
            try {
                const assertData = await civicQuery(`{
                    assertions(geneIds: [${apiGene.id}], status: ACCEPTED) {
                        nodes {
                            id ampLevel clinicalSignificance significance summary
                            disease { name }
                            therapies { nodes { name } }
                        }
                    }
                }`);
                assertions = assertData?.data?.assertions?.nodes || [];
            } catch (assertErr) {
                // Assertions unavailable — return gene data without them
                console.warn('CivicDB assertions fetch failed:', assertErr.message);
            }
        }

        return res.status(200).json({ gene: apiGene, matchedVariant, assertions });
    } catch (err) {
        console.error('CivicDB proxy error:', err);
        return res.status(502).json({ error: 'CivicDB lookup failed', detail: err.message });
    }
}
