// Serverless proxy for the CivicDB GraphQL API.
// GET /api/civic?gene=BRAF&protein=V600E
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const CIVIC_GRAPHQL = 'https://civicdb.org/api/graphql';

async function civicPost(query, variables = {}) {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 12000);
    try {
        const res = await fetch(CIVIC_GRAPHQL, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
            body: JSON.stringify({ query, variables }),
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

// Step 1: find gene ID using CIViC v2 schema
async function lookupGeneId(safeGene) {
    const errs = [];

    const geneQuery = `
        query CivicGeneBySymbol($symbol: String!) {
            gene(entrezSymbol: $symbol) {
                id
                name
            }
        }
    `;

    const queryGene = safeGene.toUpperCase();
    const rA = await civicPost(geneQuery, { symbol: queryGene });

    if (!rA.errors?.length) {
        const gene = rA.data?.gene;
        if (gene?.id) return { id: Number(gene.id), name: gene.name, queryErrors: null };
    } else {
        errs.push(...rA.errors.map(e => `gene(entrezSymbol): ${e.message}`));
    }

    const genesQuery = `
        query CivicGenesBySymbols($symbols: [String!]!) {
            genes(entrezSymbols: $symbols) {
                nodes {
                    id
                    name
                }
            }
        }
    `;

    const rB = await civicPost(genesQuery, { symbols: [queryGene] });

    if (!rB.errors?.length) {
        const nodes = rB.data?.genes?.nodes || [];
        const match = nodes.find(g => String(g.name || '').toUpperCase() === safeGene.toUpperCase()) || nodes[0];
        if (match?.id) return { id: Number(match.id), name: match.name, queryErrors: null };
    } else {
        errs.push(...rB.errors.map(e => `genes(entrezSymbols): ${e.message}`));
    }

    return { id: null, name: null, queryErrors: errs.length ? errs : ['Gene not found'] };
}

// Step 2: fetch gene data
async function fetchGeneById(geneId) {
    const geneQuery = `
        query CivicGeneById($id: Int!) {
            gene(id: $id) {
                id
                name
                description
                variants(first: 100) {
                    nodes {
                        id
                        name
                        variantTypes { name }
                        singleVariantMolecularProfile {
                            id
                            assertions(first: 10) {
                                nodes {
                                    id
                                    ampLevel
                                    assertionType
                                    significance
                                    status
                                    summary
                                    disease { name }
                                    therapies { name }
                                }
                            }
                        }
                    }
                    totalCount
                }
            }
        }
    `;

    const rA = await civicPost(geneQuery, { id: Number(geneId) });

    if (!rA.errors?.length) {
        const gene = rA.data?.gene;
        if (gene) return normaliseGene(gene);
    }

    const genesQuery = `
        query CivicGenesByIds($ids: [Int!]!) {
            genes(ids: $ids) {
                nodes {
                    id
                    name
                    description
                    variants(first: 100) {
                        nodes {
                            id
                            name
                            variantTypes { name }
                        }
                        totalCount
                    }
                }
            }
        }
    `;

    const rB = await civicPost(genesQuery, { ids: [Number(geneId)] });

    if (!rB.errors?.length) {
        const gene = rB.data?.genes?.nodes?.[0];
        if (gene) return normaliseGene(gene);
    }

    return null;
}

function normaliseGene(gene) {
    if (!gene) return gene;
    const nodes = Array.isArray(gene.variants?.nodes) ? gene.variants.nodes : [];
    return {
        ...gene,
        variants: {
            ...(gene.variants || {}),
            nodes: nodes.map((variant) => ({
                ...variant,
                variantTypes: normaliseVariantTypes(variant.variantTypes)
            }))
        }
    };
}

function normaliseVariantTypes(variantTypes) {
    if (Array.isArray(variantTypes)) return { nodes: variantTypes };
    if (Array.isArray(variantTypes?.nodes)) return variantTypes;
    return { nodes: [] };
}

// Step 3: fetch assertions
async function fetchAssertions(geneId, geneName) {
    const assertionsQuery = `
        query CivicAssertionsForGene($name: String!) {
            assertions(molecularProfileName: $name, status: ACCEPTED, first: 25) {
                nodes {
                    id
                    ampLevel
                    assertionType
                    significance
                    summary
                    disease { name }
                    therapies { name }
                }
            }
        }
    `;

    const rA = await civicPost(assertionsQuery, { name: geneName || '' });

    if (!rA.errors?.length && rA.data?.assertions) {
        return normaliseAssertions(rA.data.assertions.nodes || []);
    }

    return [];
}

function normaliseAssertions(assertions) {
    const seen = new Set();
    return assertions
        .filter((assertion) => {
            const status = String(assertion?.status || 'accepted').toLowerCase();
            return status !== 'rejected' && status !== 'submitted';
        })
        .filter((assertion) => {
            const key = assertion?.id || JSON.stringify(assertion);
            if (seen.has(key)) return false;
            seen.add(key);
            return true;
        })
        .map((assertion) => ({
            ...assertion,
            therapies: Array.isArray(assertion.therapies)
                ? { nodes: assertion.therapies }
                : assertion.therapies
        }));
}

function collectAssertionsFromGene(gene) {
    const variants = Array.isArray(gene?.variants?.nodes) ? gene.variants.nodes : [];
    const nested = [];
    for (const variant of variants) {
        const nodes = variant?.singleVariantMolecularProfile?.assertions?.nodes;
        if (Array.isArray(nodes)) nested.push(...nodes);
    }
    return normaliseAssertions(nested);
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

    const { id: geneId, name: geneName, queryErrors } = await lookupGeneId(safeGene);

    if (!geneId) {
        console.warn('CivicDB gene not found:', safeGene, queryErrors);
        return res.status(200).json({ gene: null, matchedVariant: null, assertions: [], queryErrors });
    }

    const [fullGene, fetchedAssertions] = await Promise.all([
        fetchGeneById(geneId),
        fetchAssertions(geneId, geneName)
    ]);

    const apiGene = fullGene || {
        id: geneId,
        name: geneName,
        description: '',
        variants: { nodes: [] }
    };

    let assertions = fetchedAssertions;
    if (!assertions || assertions.length === 0) {
        assertions = collectAssertionsFromGene(apiGene);
    }

    let matchedVariant = null;

    // Convert three-letter amino acid codes to single-letter so "Gly12Cys" matches CIViC's "G12C"
    const AA3TO1 = {
        ala:'a',arg:'r',asn:'n',asp:'d',cys:'c',gln:'q',glu:'e',gly:'g',
        his:'h',ile:'i',leu:'l',lys:'k',met:'m',phe:'f',pro:'p',ser:'s',
        thr:'t',trp:'w',tyr:'y',val:'v',ter:'*'
    };
    const normProtein = (s) => String(s || '')
        .replace(/^p\./i, '')
        .replace(/[A-Za-z]{3}/g, m => AA3TO1[m.toLowerCase()] || m.toLowerCase())
        .toLowerCase()
        .replace(/[^a-z0-9*_]/g, '');

    const normInput = normProtein(protein);

    const variantNodes = apiGene.variants?.nodes || [];

    if (normInput && variantNodes.length > 0) {
        for (const v of variantNodes) {
            const vn = normProtein(v.name);
            if (vn && (vn === normInput || normInput.includes(vn) || vn.includes(normInput))) {
                matchedVariant = v;
                break;
            }
        }
    }

    return res.status(200).json({ gene: apiGene, matchedVariant, assertions });
}
