/*
 * Tests for pasted coordinate-row parsing.
 *
 * Users paste variants out of spreadsheets, VCFs, MAFs and pipeline reports.
 * Only one row shape ("chr pos gene ref alt", whitespace-separated) used to
 * parse; every other shape fell through to the generic "GENE p.Change" branch
 * of normalizeVariantInput and produced a nonsense query that was then
 * reported to the user as "Variant not found. Please verify the genomic
 * coordinate and reference allele." — blaming the input for a parser gap.
 *
 * The regression that prompted these tests:
 *   "X  20148674  EIF1AX  T  TG"  (EIF1AX NM_001412.4:c.388dup)
 * parsed, but the same variant in MAF column order, as CSV, or with a trailing
 * build column did not.
 *
 * Exercises the REAL implementations, imported from script.js via its Node
 * test-export block.
 *
 * Run with: node tests/coordinate-row.test.js
 */

// --- real helpers imported from script.js ---------------------------------

await import('../script.js');
const { looksLikeCoordinateRow, parseCoordinateRow, coordinateRowToGenomicHgvs } = globalThis.__variantSearchHelpers;

// --- test harness ---------------------------------------------------------

let failures = 0;
function check(name, actual, expected) {
    const pass = actual === expected;
    if (!pass) {
        failures++;
        console.error(`FAIL ${name}\n     expected: ${expected}\n     actual:   ${actual}`);
    } else {
        console.log(`ok   ${name}`);
    }
}

const hgvsOf = (row) => coordinateRowToGenomicHgvs(parseCoordinateRow(row));

// ── The reported regression: EIF1AX c.388dup in every paste shape ──────────
// hg19 chrX:20148674 T>TG is NM_001412.4:c.388dup (p.Gln130ProfsTer3).
const EIF1AX = 'chrX:g.20148674_20148675insG';

check('tab-separated row', hgvsOf('X\t20148674\tEIF1AX\tT\tTG'), EIF1AX);
check('space-separated row', hgvsOf('X 20148674 EIF1AX T TG'), EIF1AX);
check('non-breaking spaces', hgvsOf('X 20148674 EIF1AX T TG'), EIF1AX);
check('chr-prefixed chromosome', hgvsOf('chrX 20148674 EIF1AX T TG'), EIF1AX);
check('no gene column', hgvsOf('X 20148674 T TG'), EIF1AX);
check('CSV paste', hgvsOf('X,20148674,EIF1AX,T,TG'), EIF1AX);
check('CSV paste with spaces', hgvsOf('X, 20148674, EIF1AX, T, TG'), EIF1AX);
check('semicolon separated', hgvsOf('X;20148674;EIF1AX;T;TG'), EIF1AX);
check('quoted spreadsheet cells', hgvsOf('"X","20148674","EIF1AX","T","TG"'), EIF1AX);
check('thousands separators in position', hgvsOf('X\t20,148,674\tEIF1AX\tT\tTG'), EIF1AX);
check('MAF column order (gene first)', hgvsOf('EIF1AX X 20148674 T TG'), EIF1AX);
check('trailing build column', hgvsOf('X 20148674 EIF1AX T TG hg19'), EIF1AX);
check('trailing build + zygosity', hgvsOf('X 20148674 EIF1AX T TG GRCh37 somatic'), EIF1AX);
check('MAF-style insertion ("-" ref)', hgvsOf('X 20148674 EIF1AX - G'), EIF1AX);

// ── Other event types ─────────────────────────────────────────────────────
check('SNV with gene column', hgvsOf('7 140453136 BRAF A T'), 'chr7:g.140453136A>T');
check('SNV without gene column', hgvsOf('7 140453136 A T'), 'chr7:g.140453136A>T');
check('VCF-anchored single deletion', hgvsOf('17 7578405 TP53 CC C'), 'chr17:g.7578406del');
check('MAF-style deletion ("-" alt)', hgvsOf('17 7578406 TP53 C -'), 'chr17:g.7578406del');
// ref spans 55242465..55242475; trimming the anchoring T leaves 10 deleted bases.
check('multi-base deletion', hgvsOf('7 55242465 EGFR TAAGAGAAGCA T'), 'chr7:g.55242466_55242475del');
check('delins', hgvsOf('7 140453122 BRAF TCCATCGAGATTTCA TCT'), 'chr7:g.140453124_140453136delinsT');
check('MNV becomes delins', hgvsOf('7 140453136 BRAF AC TG'), 'chr7:g.140453136_140453137delinsTG');
// Same-length pairs carrying shared context bases reduce to the minimal event —
// an untrimmed delins ("CT">"CA" as delinsCA) misses in MyVariant/ClinVar.
check('same-length pair, shared prefix → SNV', hgvsOf('7 140453135 BRAF CT CA'), 'chr7:g.140453136T>A');
check('same-length pair, shared suffix → SNV', hgvsOf('7 140453136 BRAF AG TG'), 'chr7:g.140453136A>T');
check('same-length pair, shared flanks → SNV', hgvsOf('7 140453135 BRAF CAG CTG'), 'chr7:g.140453136A>T');
check('same-length pair, inner MNV keeps delins', hgvsOf('7 140453135 BRAF CACG CTGG'), 'chr7:g.140453136_140453137delinsTG');
check('identical ref/alt is not a variant', hgvsOf('7 140453136 BRAF AC AC'), null);

// ── Gene hint extraction ──────────────────────────────────────────────────
check('gene from middle column', parseCoordinateRow('X 20148674 EIF1AX T TG').gene, 'EIF1AX');
check('gene from MAF order', parseCoordinateRow('EIF1AX X 20148674 T TG').gene, 'EIF1AX');
check('no gene column yields null', parseCoordinateRow('X 20148674 T TG').gene, null);
// A gene symbol that spells DNA must not be mistaken for the reference allele:
// alleles are always the last two DNA-like columns.
check('DNA-spelling gene symbol (TTN)', hgvsOf('2 179000000 TTN T TG'), 'chr2:g.179000000_179000001insG');
check('DNA-spelling gene kept as hint', parseCoordinateRow('2 179000000 TTN T TG').gene, 'TTN');

// ── Rows that must NOT be guessed at ──────────────────────────────────────
check('too few columns', parseCoordinateRow('X 20148674 EIF1AX T'), null);
check('no chromosome/position pair', parseCoordinateRow('foo bar baz qux'), null);
check('non-adjacent alleles', parseCoordinateRow('X 20148674 T EIF1AX TG'), null);
check('CSV row is not welded by comma stripping',
    hgvsOf('1,12345,A,T'), 'chr1:g.12345A>T');

// ── looksLikeCoordinateRow gates the "GENE p.Change" fallback ─────────────
check('detects a full row', looksLikeCoordinateRow('X 20148674 EIF1AX T TG'), true);
check('detects an unparseable row', looksLikeCoordinateRow('X 20148674 EIF1AX T ZZ'), true);
check('ignores HGVS input', looksLikeCoordinateRow('chrX:g.20148675dup'), false);
check('ignores protein input', looksLikeCoordinateRow('BRAF p.Val600Glu'), false);
check('ignores gene + descriptor', looksLikeCoordinateRow('BRAF amplification'), false);

if (failures) {
    console.error(`\n${failures} test(s) failed`);
    process.exit(1);
}
console.log('\nAll coordinate-row tests passed');
