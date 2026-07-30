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

// Normalise a cDNA HGVS string for comparison across sources. Strips any accession
// prefix and the optionally spelled-out reference bases, so the user's
// "c.2319_2321delAAT" and ClinVar's "c.2319_2321del" compare equal. Inserted bases
// are kept, because insAA and insTT are genuinely different variants.
function normalizeCdnaForMatch(cdna) {
    let s = String(cdna || '').trim().toLowerCase().replace(/\s+/g, '');
    const colon = s.lastIndexOf(':');
    if (colon !== -1) s = s.slice(colon + 1);
    if (!/^[cn]\./.test(s)) return '';
    s = s.replace(/^[cn]\./, 'c.');
    s = s.replace(/del[acgt]+(?=ins)/, 'del');
    s = s.replace(/(del|dup|inv)[acgt]+$/, '$1');
    return s;
}

// Locate the ClinVar region-pull entry that corresponds to the exact queried

// ClinVar itself has an entry — the per-variant lookup otherwise shows "N/A".
//
// variants: array of region-pull records ({ id, title, variationName, pos, germline }).
// tuple:    queried variant ({ pos, ref, alt }) from buildVariantCoordinateTuple.
//
// Matching is intentionally conservative: only single-nucleotide substitutions
// are recovered, and the substitution bases in the ClinVar title must match the
// queried alleles (allowing for transcript-strand complementation on minus-strand
// genes). When two SNVs share a position (e.g. G>C and G>A) the allele check
// disambiguates; if no entry's alleles match, no fallback is applied.
//
// cdnaForms: optional list of the queried variant's cDNA HGVS strings. Indels have
// no ref/alt to compare, and their genomic coordinates differ between HGVS (3'-shifted)
// and ClinVar's own representation, so they are matched on the c. notation in the
// record title instead — which is exactly how a user would recognise the record.
function findExactClinvarRegionMatch(variants, tuple, cdnaForms) {
    if (!Array.isArray(variants) || !tuple || !tuple.pos) return null;
    const pos = Number(tuple.pos);
    const ref = String(tuple.ref || '').toUpperCase();
    const alt = String(tuple.alt || '').toUpperCase();
    if (!/^[ACGT]$/.test(ref) || !/^[ACGT]$/.test(alt)) {
        // Non-SNV: fall back to cDNA-notation matching when forms were supplied.
        const wanted = new Set(
            (Array.isArray(cdnaForms) ? cdnaForms : [cdnaForms])
                .map(normalizeCdnaForMatch)
                .filter(Boolean)
        );
        if (wanted.size === 0) return null;
        for (const v of variants) {
            const text = [v.title, v.variationName].filter(Boolean).join(' ');
            const tokens = text.match(/c\.[A-Za-z0-9_>+*\-]+/gi) || [];
            if (tokens.some((t) => wanted.has(normalizeCdnaForMatch(t)))) return v;
        }
        return null;
    }

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

function aaThreeToSingle(aa) {
    if (!aa) return null;
    const aaSingle = {
        ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
        HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
        THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
        TER: '*', STOP: '*'
    };
    return aaSingle[String(aa).toUpperCase()] || null;
}

function parseProteinChange(text) {
    const s = String(text || '');
    let m = s.match(/p\.\(?([A-Za-z]{3})(\d+)([A-Za-z]{3})\)?/);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]);
        if (ref && alt) return { ref, pos: Number(m[2]), alt };
    }
    m = s.match(/\b([A-Za-z]{3})(\d+)([A-Za-z]{3})\b/);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]);
        if (ref && alt) return { ref, pos: Number(m[2]), alt };
    }
    m = s.match(/\bp?\.?([A-Z])(\d+)([A-Z*])\b/);
    if (m) return { ref: m[1].toUpperCase(), pos: Number(m[2]), alt: m[3].toUpperCase() };
    return null;
}

function sameProteinChange(a, b) {
    return !!a && !!b && a.ref === b.ref && a.pos === b.pos && a.alt === b.alt;
}

function aaSingleToThree(aa) {
    const map = {
        A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly',
        H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser',
        T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter'
    };
    return map[String(aa || '').toUpperCase()] || null;
}

// Build the eutils change forms the card sends to /api/clinvar-protein.
function buildChangeForms(pc) {
    const forms = [];
    const r3 = aaSingleToThree(pc.ref);
    const a3 = aaSingleToThree(pc.alt);
    if (r3 && a3) forms.push(`${r3}${pc.pos}${a3}`);
    forms.push(`${pc.ref}${pc.pos}${pc.alt}`);
    return forms;
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

// Indels carry no comparable ref/alt, so they are recovered from the c. notation
// in the record title instead. Without cDNA forms there is nothing to match on.
check('indel with no cDNA forms supplied → null',
    findExactClinvarRegionMatch(ctnnb1Region, { pos: '41266103', ref: 'G', alt: 'GG' }) === null);

// TSC2 c.2319_2321delAAT: the recoder normalises this to chr16:g.2122948_2122950del,
// which has no alleles at all. ClinVar titles it "c.2319_2321del" (bases dropped).
const tsc2Region = [
    { id: '821088', pos: 2122948, title: 'NM_000548.5(TSC2):c.2319A>G (p.Leu773=)', variationName: '', germline: 'Benign/Likely benign' },
    { id: '700001', pos: 2122948, title: 'NM_000548.5(TSC2):c.2319_2321del (p.Leu773_Ile774delinsPhe)', variationName: '', germline: 'Likely pathogenic' },
];
check('recovers a deletion by cDNA notation, ignoring spelled-out bases',
    findExactClinvarRegionMatch(tsc2Region, { pos: '2122948', ref: null, alt: null },
        ['c.2319_2321delAAT'])?.id === '700001');
check('accepts an accession-prefixed cDNA form',
    findExactClinvarRegionMatch(tsc2Region, { pos: '2122948', ref: null, alt: null },
        ['NM_000548.5:c.2319_2321del'])?.id === '700001');
check('a different deletion at the same locus does not match',
    findExactClinvarRegionMatch(tsc2Region, { pos: '2122948', ref: null, alt: null },
        ['c.2319_2320del']) === null);
check('insertions are distinguished by their inserted bases',
    findExactClinvarRegionMatch(
        [{ id: '1', pos: 100, title: 'NM_1(X):c.100_101insAA', variationName: '', germline: '' }],
        { pos: '100', ref: null, alt: null }, ['c.100_101insTT']) === null);

// normalizeCdnaForMatch itself.
check('strips deleted bases', normalizeCdnaForMatch('c.2319_2321delAAT') === 'c.2319_2321del');
check('strips duplicated bases', normalizeCdnaForMatch('NM_000548.5:c.100_102dupCTG') === 'c.100_102dup');
check('keeps inserted bases', normalizeCdnaForMatch('c.100_101insAA') === 'c.100_101insaa');
check('strips the deleted side of a delins only',
    normalizeCdnaForMatch('c.100_102delAATinsG') === 'c.100_102delinsg');
check('accepts n. notation from non-coding transcripts',
    normalizeCdnaForMatch('ENST00000463808.1:n.353_355del') === 'c.353_355del');
check('rejects protein notation', normalizeCdnaForMatch('p.Val600Glu') === '');

// Guard inputs.
check('null variants → null', findExactClinvarRegionMatch(null, { pos: '1', ref: 'A', alt: 'C' }) === null);
check('missing tuple → null', findExactClinvarRegionMatch(ctnnb1Region, null) === null);
check('complementBase round-trips', complementBase('A') === 'T' && complementBase('g') === 'C');

// --- protein-change parsing / matching ------------------------------------

// ClinVar titles (three-letter, parenthesised) parse to a single-letter key.
check('parses p.(Gly34Arg) from a ClinVar title',
    sameProteinChange(
        parseProteinChange('NM_001904.4(CTNNB1):c.100G>C (p.Gly34Arg)'),
        { ref: 'G', pos: 34, alt: 'R' }));

// The whole point of this feature: a different codon nucleotide (c.100G>A)
// yielding the SAME amino-acid change must compare equal to the queried one.
check('c.100G>C and c.100G>A both map to G34R and match',
    sameProteinChange(
        parseProteinChange('NM_001904.4(CTNNB1):c.100G>C (p.Gly34Arg)'),
        parseProteinChange('NM_001904.4(CTNNB1):c.100G>A (p.Gly34Arg)')));

// Query-side forms: compact triple (targetProtGlobal) and single-letter.
check('parses compact triple GLY34ARG', sameProteinChange(parseProteinChange('GLY34ARG'), { ref: 'G', pos: 34, alt: 'R' }));
check('parses single-letter G34R', sameProteinChange(parseProteinChange('G34R'), { ref: 'G', pos: 34, alt: 'R' }));
check('parses gene-prefixed p. string', sameProteinChange(parseProteinChange('CTNNB1:p.Gly34Arg'), { ref: 'G', pos: 34, alt: 'R' }));

// A different residue or a different substituted AA must NOT match.
check('different position does not match',
    !sameProteinChange(parseProteinChange('p.Gly34Arg'), parseProteinChange('p.Gly35Arg')));
check('different substituted AA does not match (G34R vs G34V)',
    !sameProteinChange(parseProteinChange('p.Gly34Arg'), parseProteinChange('p.Gly34Val')));

// Nonsense parses (Ter → *); a transcript accession alone yields no change.
check('parses nonsense Arg213Ter → R213*',
    sameProteinChange(parseProteinChange('p.Arg213Ter'), { ref: 'R', pos: 213, alt: '*' }));
check('bare transcript text yields no protein change',
    parseProteinChange('NM_001904.4') === null);

// The eutils search forms the card sends: three-letter (matches ClinVar titles)
// plus single-letter, so the gene+change query has the best recall.
check('builds three-letter + single-letter change forms for G34R',
    JSON.stringify(buildChangeForms({ ref: 'G', pos: 34, alt: 'R' })) === JSON.stringify(['Gly34Arg', 'G34R']));
check('aaSingleToThree maps stop codon * → Ter', aaSingleToThree('*') === 'Ter');
check('three-letter search form parses back to the same change',
    sameProteinChange(parseProteinChange(buildChangeForms({ ref: 'R', pos: 282, alt: 'W' })[0]), { ref: 'R', pos: 282, alt: 'W' }));

// --- review-status stars (copied verbatim from script.js) -----------------

function escapeHtml(text) {
    return String(text ?? '').replace(/[&<>"']/g, (ch) => (
        { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch]
    ));
}

function clinvarReviewStars(status) {
    const s = String(status || '').toLowerCase();
    if (!s) return '';
    let filled = 0;
    if (s.includes('no assertion') || s.includes('no classification')) filled = 0;
    else if (s.includes('practice guideline')) filled = 4;
    else if (s.includes('expert panel')) filled = 3;
    else if (s.includes('multiple submitters')) filled = 2;
    else if (s.includes('criteria provided')) filled = 1;
    else filled = 0;
    const stars = '★'.repeat(filled) + '☆'.repeat(4 - filled);
    return `<span style="color:#f59e0b" title="${escapeHtml(status)}">${stars}</span>`;
}

// Extract just the ★/☆ glyphs from the rendered span for easy assertions.
function starsOf(status) {
    const html = clinvarReviewStars(status);
    const m = html.match(/>([★☆]+)</);
    return m ? m[1] : '';
}

check('empty status renders no stars element', clinvarReviewStars('') === '');
// The regression this guards: "no assertion criteria provided" contains the
// substring "criteria provided" but must score ZERO stars, not one.
check('"no assertion criteria provided" → 0 stars', starsOf('no assertion criteria provided') === '☆☆☆☆');
check('"criteria provided, single submitter" → 1 star', starsOf('criteria provided, single submitter') === '★☆☆☆');
check('"criteria provided, conflicting classifications" → 1 star', starsOf('criteria provided, conflicting classifications') === '★☆☆☆');
check('"criteria provided, multiple submitters, no conflicts" → 2 stars', starsOf('criteria provided, multiple submitters, no conflicts') === '★★☆☆');
check('"reviewed by expert panel" → 3 stars', starsOf('reviewed by expert panel') === '★★★☆');
check('"practice guideline" → 4 stars', starsOf('practice guideline') === '★★★★');
check('review-status tooltip is html-escaped', clinvarReviewStars('a & b').includes('title="a &amp; b"'));
check('escapeHtml neutralises angle brackets', escapeHtml('<x>') === '&lt;x&gt;');

// --- summary --------------------------------------------------------------

console.log(`\nClinVar region match tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
