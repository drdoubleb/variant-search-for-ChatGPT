/*
 * Tests for the ClinVar region exact-match fallback.
 *
 * When MyVariant.info carries no `clinvar` block for the queried variant, the
 * ClinVar card recovers the variation ID from the live region pull by matching
 * on exact genomic position + substitution. The helpers live in script.js
 * (complementBase / findExactClinvarRegionMatch). Because the project has no
 * module system, they are re-declared here verbatim — KEEP IN SYNC with
 * script.js when either side changes.
 *
 * Run with: node tests/clinvar-match.test.js
 */

// --- helpers copied verbatim from script.js -------------------------------

function complementBase(b) {
    return { A: 'T', T: 'A', C: 'G', G: 'C' }[String(b || '').toUpperCase()] || '';
}

function findExactClinvarRegionMatch(variants, tuple) {
    if (!Array.isArray(variants) || !tuple || !tuple.pos) return null;
    const pos = Number(tuple.pos);
    const ref = String(tuple.ref || '').toUpperCase();
    const alt = String(tuple.alt || '').toUpperCase();
    if (!/^[ACGT]$/.test(ref) || !/^[ACGT]$/.test(alt)) return null;

    const refC = complementBase(ref);
    const altC = complementBase(alt);
    const candidates = variants.filter((v) => Number(v.pos) === pos);

    for (const v of candidates) {
        const text = [v.title, v.variationName].filter(Boolean).join(' ');
        const m = text.match(/([ACGT])\s*>\s*([ACGT])/i);
        if (!m) continue;
        const tRef = m[1].toUpperCase();
        const tAlt = m[2].toUpperCase();
        if ((tRef === ref && tAlt === alt) || (tRef === refC && tAlt === altC)) {
            return v;
        }
    }
    return null;
}

// --- tiny assertion harness ----------------------------------------------

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; }
    else { failed++; console.error(`  ✗ ${name}`); }
}

// --- fixtures -------------------------------------------------------------

// Two CTNNB1 SNVs sharing chr3:41266103 (the real-world disambiguation case):
// c.100G>C (G34R) and c.100G>A (G34R). MyVariant lacked a block for the queried
// c.100G>C, so it must be recovered from the region pull by allele, not position.
const ctnnb1Region = [
    { id: '375941', pos: 41266103, title: 'NM_001904.4(CTNNB1):c.100G>C (p.Gly34Arg)', variationName: 'NM_001904.4(CTNNB1):c.100G>C (p.Gly34Arg)', germline: '' },
    { id: '17578',  pos: 41266103, title: 'NM_001904.4(CTNNB1):c.100G>A (p.Gly34Arg)', variationName: 'NM_001904.4(CTNNB1):c.100G>A (p.Gly34Arg)', germline: 'Pathogenic' },
    { id: '99999',  pos: 41266110, title: 'NM_001904.4(CTNNB1):c.107T>C (p.Leu36Pro)', variationName: '', germline: 'Likely pathogenic' },
];

// --- tests ----------------------------------------------------------------

// Position alone is ambiguous (two SNVs); the alt allele must disambiguate.
check('matches c.100G>C and not c.100G>A',
    findExactClinvarRegionMatch(ctnnb1Region, { pos: '41266103', ref: 'G', alt: 'C' })?.id === '375941');
check('matches c.100G>A and not c.100G>C',
    findExactClinvarRegionMatch(ctnnb1Region, { pos: '41266103', ref: 'G', alt: 'A' })?.id === '17578');

// No ClinVar entry for the queried allele at that position → no false match.
check('no match when queried allele absent at position',
    findExactClinvarRegionMatch(ctnnb1Region, { pos: '41266103', ref: 'G', alt: 'T' }) === null);

// Minus-strand gene: the ClinVar title's coding-strand bases are the complement
// of the genomic ref/alt. A genomic C>T must match a transcript-strand G>A entry.
const minusStrandRegion = [
    { id: '12345', pos: 7577120, title: 'NM_000546.6(TP53):c.844C>T (p.Arg282Trp)', variationName: '', germline: 'Pathogenic' },
];
check('matches across strand complement (genomic G>A ↔ coding C>T)',
    findExactClinvarRegionMatch(minusStrandRegion, { pos: '7577120', ref: 'G', alt: 'A' })?.id === '12345');

// Indels are not recovered (allele matching is unreliable for them).
check('does not attempt recovery for indels',
    findExactClinvarRegionMatch(ctnnb1Region, { pos: '41266103', ref: 'G', alt: 'GG' }) === null);

// Guard inputs.
check('null variants → null', findExactClinvarRegionMatch(null, { pos: '1', ref: 'A', alt: 'C' }) === null);
check('missing tuple → null', findExactClinvarRegionMatch(ctnnb1Region, null) === null);
check('complementBase round-trips', complementBase('A') === 'T' && complementBase('g') === 'C');

// --- summary --------------------------------------------------------------

console.log(`\nClinVar region match tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
