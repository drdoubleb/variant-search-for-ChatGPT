/*
 * Tests for the nearby-ClinVar lollipop span helpers.
 *
 * Two things went into the plot rework:
 *   - getQueryAffectedSpan: ClinVar's region records place an indel at its first
 *     affected base, while a VCF-anchored tuple points at the untouched anchor
 *     base 5' of it — which put every indel query's diamond one base off from
 *     ClinVar's own record of the same variant. The helper normalises the query
 *     tuple to the ClinVar display convention.
 *   - parseClinvarEventKind: classifies a record's sequence-change kind from the
 *     HGVS tokens in its title so in-frame del/dup/delins draw as span bars and
 *     insertions as carets, without gene symbols like DELE1/INSR misclassifying.
 *
 * Exercises the REAL implementations, imported from script.js via its Node
 * test-export block.
 *
 * Run with: node tests/lollipop-span.test.js
 */

await import('../script.js');
const {
    getQueryAffectedSpan, parseClinvarEventKind, buildVariantCoordinateTuple,
    isTruncatingClinvarVariant, getClinvarEffectiveClassification, getPathogenicityColor,
    mapVariantToResidueSpan
} = globalThis.__variantSearchHelpers;

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; }
    else { failed++; console.error(`  ✗ ${name}`); }
}
function spanEq(span, kind, start, stop) {
    return span && span.kind === kind && span.start === start && span.stop === stop;
}

// --- getQueryAffectedSpan --------------------------------------------------

// SNV: passes through unchanged.
check('SNV span is the substituted base',
    spanEq(getQueryAffectedSpan({ pos: '100', ref: 'G', alt: 'A', start: 100, end: 100, type: 'sub' }), 'sub', 100, 100));

// VCF-anchored deletion (the real off-by-one case): CFTR F508del pasted as a
// VCF row "7 117199644 CTT C" — ClinVar shows the record at 117199645–117199647.
const f508 = buildVariantCoordinateTuple('7 117199644 CFTR CTTT CT', '');
check('VCF-anchored deletion starts at the first deleted base',
    spanEq(getQueryAffectedSpan(f508), 'del', 117199645, 117199646));

check('simple anchored del: CTT>C',
    spanEq(getQueryAffectedSpan({ pos: 100, ref: 'CTT', alt: 'C', type: 'delins' }), 'del', 101, 102));

// VCF-anchored insertion: T>TG at 100 inserts between 100 and 101.
check('VCF-anchored insertion spans its two flanking bases',
    spanEq(getQueryAffectedSpan({ pos: 100, ref: 'T', alt: 'TG', type: 'delins' }), 'ins', 100, 101));

// Shared suffix instead of prefix: A>GA is an insertion before position 100.
check('suffix-anchored insertion resolves left of the anchor',
    spanEq(getQueryAffectedSpan({ pos: 100, ref: 'A', alt: 'GA', type: 'delins' }), 'ins', 99, 100));

// Anchored delins: CAG>CT deletes AG (101–102) and inserts T.
check('anchored delins spans the replaced bases',
    spanEq(getQueryAffectedSpan({ pos: 100, ref: 'CAG', alt: 'CT', type: 'delins' }), 'delins', 101, 102));

// MAF-style rows spell the absent allele as "" — position conventions differ
// from VCF (no anchor base is included).
check('MAF-style deletion keeps its own first-deleted-base convention',
    spanEq(getQueryAffectedSpan({ pos: 100, ref: 'TG', alt: '', type: 'delins' }), 'del', 100, 101));
check('MAF-style insertion lands between pos and pos+1',
    spanEq(getQueryAffectedSpan({ pos: 100, ref: '', alt: 'TG', type: 'delins' }), 'ins', 100, 101));

// HGVS-parsed tuples already carry the affected span.
const hgvsDel = buildVariantCoordinateTuple('', 'chr7:g.117199645_117199647del');
check('HGVS deletion span passes through',
    spanEq(getQueryAffectedSpan(hgvsDel), 'del', 117199645, 117199647));
const hgvsIns = buildVariantCoordinateTuple('', 'chr7:g.117199644_117199645insAAA');
check('HGVS insertion keeps its flanking-base span',
    spanEq(getQueryAffectedSpan(hgvsIns), 'ins', 117199644, 117199645));
const hgvsDup = buildVariantCoordinateTuple('', 'chr17:g.41245466_41245468dup');
check('HGVS duplication span passes through',
    spanEq(getQueryAffectedSpan(hgvsDup), 'dup', 41245466, 41245468));

// Same-length MNV: tuple trimming already moved pos to the first changed base.
const mnv = buildVariantCoordinateTuple('7 117199644 CFTR CT CA', '');
check('MNV trimmed to changed bases is not treated as an indel',
    spanEq(getQueryAffectedSpan(mnv), 'sub', 117199645, 117199645));

check('null/invalid tuples return null',
    getQueryAffectedSpan(null) === null && getQueryAffectedSpan({ pos: 'x' }) === null);

// --- parseClinvarEventKind -------------------------------------------------

check('substitution titles have no event kind',
    parseClinvarEventKind('NM_000546.6(TP53):c.818G>A (p.Arg273His)') === null);
check('deletion title', parseClinvarEventKind('NM_000492.4(CFTR):c.1521_1523del (p.Phe508del)') === 'del');
check('delins wins over del/ins', parseClinvarEventKind('NM_000492.4(CFTR):c.1509_1510delinsT (p.Lys503fs)') === 'delins');
check('duplication title', parseClinvarEventKind('NM_000492.4(CFTR):c.1496dup (p.Gly500fs)') === 'dup');
check('insertion title', parseClinvarEventKind('NM_007294.4(BRCA1):c.5266_5267insC (p.Gln1756fs)') === 'ins');
check('repeat notation classified via its p. token',
    parseClinvarEventKind('NM_000492.4(CFTR):c.1516ATC[1] (p.Ile507del)') === 'del');
check('haplotype titles with bracketed alleles still classify',
    parseClinvarEventKind('NM_000492.3(CFTR):c.[1521_1523delCTT;3080T>C]') === 'del');
check('gene symbols containing del/ins never misclassify',
    parseClinvarEventKind('NM_024316.3(DELE1):c.100A>G (p.Thr34Ala)') === null
    && parseClinvarEventKind('NM_000208.4(INSR):c.100A>G (p.Thr34Ala)') === null);
check('empty/garbage input', parseClinvarEventKind('') === null && parseClinvarEventKind('GRCh37/hg19 7q31.2') === null);

// Frameshift deletions stay circles: the plot only draws span glyphs for
// non-truncating records, so the truncating check must still catch these.
check('frameshift del title is still flagged truncating',
    isTruncatingClinvarVariant({ title: 'NM_000492.4(CFTR):c.1547_1550del (p.Arg516fs)' }) === true);
check('in-frame del title is not flagged truncating',
    isTruncatingClinvarVariant({ title: 'NM_000492.3(CFTR):c.1521_1523del (p.Phe508del)' }) === false);

// --- getClinvarEffectiveClassification --------------------------------------

// BRAF V600E: conflicting germline aggregate, but Oncogenic / Tier I somatic —
// the somatic call should win and color red, with fromSomatic marking it.
const v600e = {
    germline: 'Conflicting classifications of pathogenicity',
    somatic: 'Tier I - Strong', oncogenicity: 'Oncogenic'
};
{
    const eff = getClinvarEffectiveClassification(v600e);
    check('conflicting germline defers to oncogenicity', eff.label === 'Oncogenic' && eff.fromSomatic === true);
    check('oncogenic colors red', getPathogenicityColor(eff.label) === '#dc2626');
}
check('Tier I wins when only clinical impact is present',
    getClinvarEffectiveClassification({ germline: '', somatic: 'Tier I - Strong', oncogenicity: '' }).label === 'Tier I - Strong');
check('Tier II does not count as definitive',
    getClinvarEffectiveClassification({ germline: '', somatic: 'Tier II - Potential', oncogenicity: '' }).fromSomatic === false);
check('a real germline call is never overridden',
    getClinvarEffectiveClassification({ germline: 'Benign', oncogenicity: 'Oncogenic' }).label === 'Benign');
check('likely oncogenic colors like likely pathogenic', getPathogenicityColor('Likely oncogenic') === '#ef4444');
check('benign oncogenicity never masquerades as oncogenic',
    getClinvarEffectiveClassification({ germline: '', oncogenicity: 'Benign', somatic: '' }).fromSomatic === false);
check('records without somatic data keep their germline label',
    getClinvarEffectiveClassification({ germline: 'Uncertain significance' }).label === 'Uncertain significance');

// --- mapVariantToResidueSpan -----------------------------------------------

// Context: BRCA1 C61G — 17:41258504 A>C (GRCh37), minus strand, c.181, p.61.
const brca1Ctx = { queryProteinPos: 61, queryGenomicPos: 41258504, queryC: 181, minusStrand: true };

// The real regression: a gross deletion (g.(?_41258383)_(41258523_?)del) whose
// genomic span covers the query base. Its 3' endpoint maps to residue ~54
// (plausible); its 5' endpoint lands ~40 residues away (implausible), and the
// title has no p. token — this used to drop the record from the protein plot
// entirely. It should now anchor at 54 and extend open toward higher residues.
{
    const grossDel = { title: 'NC_000017.10:g.(?_41258383)_(41258523_?)del', pos: 41258383, stop: 41258523 };
    const r = mapVariantToResidueSpan(grossDel, brca1Ctx);
    check('gross del anchors at its mappable endpoint', r.pos === 54);
    check('gross del extends open toward higher residues', r.openEnd === true && r.openStart === false);
}
// A deletion spanning far past both edges, still covering the query base:
// neither endpoint maps, but coverage of the query residue is certain.
{
    const hugeDel = { title: 'NC_000017.10:g.41250000_41260000del', pos: 41250000, stop: 41260000 };
    const r = mapVariantToResidueSpan(hugeDel, brca1Ctx);
    check('window-spanning del renders open at both ends',
        r.pos === 61 && r.openStart === true && r.openEnd === true);
}
// A deletion entirely elsewhere (not covering the query, both ends implausible)
// with no p. token stays dropped — no smearing across the window.
{
    const farDel = { title: 'NC_000017.10:g.41260000_41262000del', pos: 41260000, stop: 41262000 };
    const r = mapVariantToResidueSpan(farDel, brca1Ctx);
    check('far-away del without p. notation is still dropped', r.pos === null);
}
// Plain SNV mapping unchanged: 41258504 → residue 61; both-plausible span closed.
check('SNV maps to its residue with no open flags', (() => {
    const r = mapVariantToResidueSpan({ title: 'NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)', pos: 41258504, stop: 41258504 }, brca1Ctx);
    return r.pos === 61 && r.stop === null && !r.openStart && !r.openEnd;
})());
check('fully-mappable del keeps a closed residue span', (() => {
    const r = mapVariantToResidueSpan({ title: 'NM_007294.4(BRCA1):c.180_182del', pos: 41258503, stop: 41258505 }, brca1Ctx);
    return r.pos !== null && !r.openStart && !r.openEnd;
})());
// Title p. fallback still works without genomic context.
check('p. range fallback without mapping context', (() => {
    const r = mapVariantToResidueSpan({ title: 'NM_000492.4(CFTR):c.1493_1507del (p.Met498_Ile502del)', pos: null, stop: null }, { queryProteinPos: 500 });
    return r.pos === 498 && r.stop === 502;
})());
// Frameshift spans are not rescued as bars (they render as ✕ from the title).
check('frameshift del is not span-rescued', (() => {
    const r = mapVariantToResidueSpan({ title: 'NM_007294.4(BRCA1):c.68_5000del (p.Glu23fs)', pos: 41200000, stop: 41277500 }, brca1Ctx);
    return r.openStart === false && r.openEnd === false;
})());

// --- summary ---------------------------------------------------------------

console.log(`lollipop-span: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
