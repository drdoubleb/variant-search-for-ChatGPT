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
    const PAGE_SIZE = 100;
    const MAX_VARIANTS = 500;
    const geneQuery = `
        query CivicGeneById($id: Int!, $first: Int!, $after: String) {
            gene(id: $id) {
                id
                name
                description
                variants(first: $first, after: $after) {
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
                    pageInfo {
                        hasNextPage
                        endCursor
                    }
                    totalCount
                }
            }
        }
    `;

    let gene = null;
    const variantNodes = [];
    let after = null;

    while (variantNodes.length < MAX_VARIANTS) {
        const first = Math.min(PAGE_SIZE, MAX_VARIANTS - variantNodes.length);
        const rA = await civicPost(geneQuery, { id: Number(geneId), first, after });

        if (rA.errors?.length) break;
        const pageGene = rA.data?.gene;
        if (!pageGene) break;

        if (!gene) gene = pageGene;
        const pageVariants = Array.isArray(pageGene.variants?.nodes) ? pageGene.variants.nodes : [];
        variantNodes.push(...pageVariants);

        const pageInfo = pageGene.variants?.pageInfo || {};
        after = pageInfo.endCursor || null;
        if (!pageInfo.hasNextPage || !after || pageVariants.length === 0) break;
    }

    if (gene) {
        const seen = new Set();
        const dedupedVariants = variantNodes.filter((variant) => {
            const key = variant?.id || `${variant?.name || ''}:${seen.size}`;
            if (seen.has(key)) return false;
            seen.add(key);
            return true;
        });
        return normaliseGene({
            ...gene,
            variants: {
                ...(gene.variants || {}),
                nodes: dedupedVariants,
                totalCount: gene.variants?.totalCount ?? dedupedVariants.length
            }
        });
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
        const fallbackGene = rB.data?.genes?.nodes?.[0];
        if (fallbackGene) return normaliseGene(fallbackGene);
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

// Assertions for the matched variant specifically, plus gene-wide assertions as a
// fallback/supplement (in case a variant isn't covered per-variant). Variant-specific
// assertions come first and are tagged scope: 'variant'; remaining gene-level ones are
// tagged scope: 'gene'. Deduplicated by id across both sets. (CIViC molecular profiles
// are named per-variant, e.g. "BRAF V600E", so the previous
// assertions(molecularProfileName: <geneSymbol>) query never matched and was dead.)
function buildAssertions(matchedVariant, gene) {
    const variantNodes = matchedVariant?.singleVariantMolecularProfile?.assertions?.nodes || [];
    const variantAssertions = normaliseAssertions(variantNodes).map(a => ({ ...a, scope: 'variant' }));
    const geneAssertions = collectAssertionsFromGene(gene).map(a => ({ ...a, scope: 'gene' }));

    const seen = new Set(variantAssertions.map(a => a.id).filter(id => id != null));
    const merged = [...variantAssertions];
    for (const a of geneAssertions) {
        if (a.id != null && seen.has(a.id)) continue;
        if (a.id != null) seen.add(a.id);
        merged.push(a);
    }
    return merged;
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

    const fullGene = await fetchGeneById(geneId);

    const apiGene = fullGene || {
        id: geneId,
        name: geneName,
        description: '',
        variants: { nodes: [] }
    };

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
        // Prefer exact match; among substring matches take the longest vn (most specific variant).
        let bestSubstring = null;
        for (const v of variantNodes) {
            const vn = normProtein(v.name);
            if (!vn) continue;
            if (vn === normInput) { matchedVariant = v; break; }
            if (normInput.includes(vn) || vn.includes(normInput)) {
                if (!bestSubstring || vn.length > normProtein(bestSubstring.name).length) {
                    bestSubstring = v;
                }
            }
        }
        if (!matchedVariant) matchedVariant = bestSubstring;
    }

    // Assertions scoped to the matched variant, plus gene-wide as supplement/fallback.
    const assertions = buildAssertions(matchedVariant, apiGene);

    return res.status(200).json({ gene: apiGene, matchedVariant, assertions });
}
