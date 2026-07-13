/*
 * Tests for CIViC assertion assembly (helpers copied verbatim from api/civic.js —
 * KEEP IN SYNC). buildAssertions merges the matched variant's assertions (scope
 * 'variant') with gene-wide assertions (scope 'gene'), deduped by id, and drops
 * rejected/submitted assertions. This replaces a dead
 * assertions(molecularProfileName: <geneSymbol>) query.
 *
 * Run with: node tests/civic-assertions.test.js
 */

// --- helpers copied verbatim from api/civic.js ----------------------------

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

// --- fixtures -------------------------------------------------------------

// Matched variant carries assertions 1 and 2 (one accepted, one submitted).
const matchedVariant = {
    name: 'V600E',
    singleVariantMolecularProfile: {
        assertions: { nodes: [
            { id: 1, ampLevel: 'TIER_I_LEVEL_A', status: 'ACCEPTED' },
            { id: 9, ampLevel: 'TIER_II', status: 'SUBMITTED' }
        ] }
    }
};
// Gene has the matched variant (dup ids 1 + 9) plus another variant with assertion 3,
// plus a rejected assertion 4.
const gene = {
    variants: { nodes: [
        matchedVariant,
        { name: 'V600K', singleVariantMolecularProfile: { assertions: { nodes: [
            { id: 3, ampLevel: 'TIER_I_LEVEL_B', status: 'ACCEPTED' },
            { id: 4, ampLevel: 'TIER_III', status: 'REJECTED' }
        ] } } }
    ] }
};

// --- assertions -----------------------------------------------------------

let passed = 0, failed = 0;
function check(name, cond) { if (cond) passed++; else { failed++; console.error(`FAIL: ${name}`); } }

const result = buildAssertions(matchedVariant, gene);
const byId = Object.fromEntries(result.map(a => [a.id, a]));

check('submitted/rejected assertions dropped', !byId[9] && !byId[4]);
check('accepted assertions kept (1 variant, 3 gene)', !!byId[1] && !!byId[3] && result.length === 2);
check('matched-variant assertion tagged scope=variant', byId[1].scope === 'variant');
check('other-variant assertion tagged scope=gene', byId[3].scope === 'gene');
check('matched-variant assertion is not duplicated as gene scope', result.filter(a => a.id === 1).length === 1);
check('variant assertions come first', result[0].scope === 'variant');

// No matched variant → gene-wide only (all accepted assertions across the gene's
// variants: id 1 from V600E and id 3 from V600K), all tagged scope 'gene'.
const geneOnly = buildAssertions(null, gene);
check('null matched variant → gene-scoped only',
    geneOnly.length === 2
    && geneOnly.every(a => a.scope === 'gene')
    && geneOnly.map(a => a.id).sort().join(',') === '1,3');

// Empty gene → empty result.
check('empty gene → no assertions', buildAssertions(null, { variants: { nodes: [] } }).length === 0);

console.log(`\nCIViC assertion tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
