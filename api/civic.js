// Serverless proxy for the CivicDB GraphQL API.
// GET /api/civic?gene=BRAF&protein=V600E
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const CIVIC_GRAPHQL = 'https://civicdb.org/api/graphql';

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

// Step 1: find the numeric gene ID by name.
// Tries geneTypeahead (CivicDB v2), then genes(name:), then features(name:).
async function lookupGeneId(safeGene) {
    const errs = [];

    // Attempt A — geneTypeahead (most reliable name search in CivicDB v2)
    const rA = await civicPost(`{ geneTypeahead(queryString: "${safeGene}") { id name } }`);
    if (!rA.errors?.length) {
        const arr = rA.data?.geneTypeahead;
        if (Array.isArray(arr) && arr.length > 0) {
            const match = arr.find(g => String(g.name || '').toUpperCase() === safeGene.toUpperCase()) || arr[0];
            if (match?.id) return { id: match.id, name: match.name, queryErrors: null };
        }
    } else {
        errs.push(...rA.errors.map(e => `typeahead: ${e.message}`));
    }

    // Attempt B — genes(name:) connection
    const rB = await civicPost(`{ genes(name: "${safeGene}") { nodes { id name } } }`);
    if (!rB.errors?.length) {
        const node = rB.data?.genes?.nodes?.[0];
        if (node?.id) return { id: node.id, name: node.name, queryErrors: null };
    } else {
        errs.push(...rB.errors.map(e => `genes: ${e.message}`));
    }

    // Attempt C — features(name:) (CivicDB may have renamed the type)
    const rC = await civicPost(`{ features(name: "${safeGene}") { nodes { id name } } }`);
    if (!rC.errors?.length) {
        const node = rC.data?.features?.nodes?.[0];
        if (node?.id) return { id: node.id, name: node.name, queryErrors: null };
    } else {
        errs.push(...rC.errors.map(e => `features: ${e.message}`));
    }

    return { id: null, name: null, queryErrors: errs.length ? errs : ['Gene not found'] };
}

// Step 2: fetch full gene data by numeric ID.
async function fetchGeneById(geneId) {
    // Try genes(ids: [...]) connection
    const rA = await civicPost(`{
        genes(ids: [${geneId}]) {
            nodes {
                id name description
                variants { nodes { id name variantTypes { nodes { name } } } }
            }
        }
    }`);
    if (!rA.errors?.length) {
        const node = rA.data?.genes?.nodes?.[0];
        if (node) return node;
    }

    // Fallback: features(ids: [...])
    const rB = await civicPost(`{
        features(ids: [${geneId}]) {
            nodes {
                id name
                ... on Gene {
                    description
                    variants { nodes { id name variantTypes { nodes { name } } } }
                }
            }
        }
    }`);
    if (!rB.errors?.length) {
        const node = rB.data?.features?.nodes?.[0];
        if (node) return node;
    }

    // Minimal fallback — just id + name from typeahead step
    return null;
}

// Step 3: fetch accepted assertions for a gene.
async function fetchAssertions(geneId) {
    const rA = await civicPost(`{
        assertions(geneIds: [${geneId}], status: ACCEPTED) {
            nodes {
                id ampLevel clinicalSignificance significance summary
                disease { name }
                therapies { nodes { name } }
            }
        }
    }`);
    if (!rA.errors?.length && rA.data?.assertions) return rA.data.assertions.nodes || [];

    const rB = await civicPost(`{
        assertions(featureIds: [${geneId}], status: ACCEPTED) {
            nodes {
                id ampLevel clinicalSignificance significance summary
                disease { name }
                therapies { nodes { name } }
            }
        }
    }`);
    if (!rB.errors?.length && rB.data?.assertions) return rB.data.assertions.nodes || [];

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

    // Step 1: resolve gene ID by name
    const { id: geneId, name: geneName, queryErrors } = await lookupGeneId(safeGene);
    if (!geneId) {
        console.warn('CivicDB gene not found:', safeGene, queryErrors);
        return res.status(200).json({ gene: null, matchedVariant: null, assertions: [], queryErrors });
    }

    // Steps 2 + 3 in parallel
    const [fullGene, assertions] = await Promise.all([
        fetchGeneById(geneId),
        fetchAssertions(geneId)
    ]);

    // Build gene object (use full data if available, else minimal from typeahead)
    const apiGene = fullGene || { id: geneId, name: geneName, description: '', variants: { nodes: [] } };

    // Match variant by protein change
    let matchedVariant = null;
    const normInput = String(protein || '').replace(/^p\./i, '').toLowerCase().replace(/[^a-z0-9*_]/g, '');
    const variantNodes = apiGene.variants?.nodes || [];
    if (normInput && variantNodes.length > 0) {
        for (const v of variantNodes) {
            const vn = String(v.name || '').toLowerCase().replace(/[^a-z0-9*_]/g, '');
            if (vn && (vn === normInput || normInput.includes(vn) || vn.includes(normInput))) {
                matchedVariant = v;
                break;
            }
        }
    }

    return res.status(200).json({ gene: apiGene, matchedVariant, assertions });
}
