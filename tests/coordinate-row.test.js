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
 * Helpers are copied verbatim from script.js (the project has no module system)
 * — KEEP IN SYNC when either side changes.
 *
 * Run with: node tests/coordinate-row.test.js
 */

// --- helpers copied verbatim from script.js -------------------------------

const COORDINATE_ROW_NOISE = /^(hg18|hg19|hg38|grch3[678]|b3[678]|chr|snp|snv|ins|del|indel|mnp|somatic|germline|het|hom|\+|\.)$/i;

const isDnaAllele = (t) => /^(?:[ACGTN]+|-)$/i.test(String(t));
const isChromToken = (t) => /^(?:chr)?(?:[0-9]{1,2}|X|Y|M|MT)$/i.test(String(t));

function splitCoordinateRow(raw) {
    return String(raw || '')
        .trim()
        // Strip thousands separators inside a coordinate ("20,148,674") before
        // commas are treated as column separators. The grouping pattern is
        // required so a genuine CSV row ("1,12345,A,T") is not welded together.
        .replace(/\b\d{1,3}(?:,\d{3})+(?!\d)/g, (m) => m.replace(/,/g, ''))
        .split(/[\s,;|]+/)
        .map((t) => t.replace(/^["']+|["']+$/g, ''))
        .filter(Boolean);
}

function looksLikeCoordinateRow(raw) {
    const toks = splitCoordinateRow(raw).filter((t) => !COORDINATE_ROW_NOISE.test(t));
    if (toks.length < 4) return false;
    if (/:[gcp]\./i.test(String(raw))) return false;
    for (let i = 0; i < toks.length - 1; i++) {
        if (isChromToken(toks[i]) && /^\d[\d,]*$/.test(toks[i + 1])) return true;
    }
    return false;
}

function parseCoordinateRow(raw) {
    const all = splitCoordinateRow(raw);
    if (all.length < 4) return null;
    const toks = all.filter((t) => !COORDINATE_ROW_NOISE.test(t));
    if (toks.length < 4) return null;

    let chromIdx = -1;
    for (let i = 0; i < toks.length - 1; i++) {
        if (isChromToken(toks[i]) && /^\d[\d,]*$/.test(toks[i + 1])) { chromIdx = i; break; }
    }
    if (chromIdx === -1) return null;

    const chrom = toks[chromIdx].replace(/^chr/i, '').toUpperCase();
    const pos = toks[chromIdx + 1].replace(/,/g, '');
    if (!/^\d+$/.test(pos)) return null;

    const after = toks.slice(chromIdx + 2);
    const dnaIdx = [];
    after.forEach((t, i) => { if (isDnaAllele(t)) dnaIdx.push(i); });
    if (dnaIdx.length < 2) return null;
    const refIdx = dnaIdx[dnaIdx.length - 2];
    const altIdx = dnaIdx[dnaIdx.length - 1];
    if (altIdx !== refIdx + 1) return null;

    const clean = (t) => (t === '-' || t === '.' ? '' : t.toUpperCase());
    const ref = clean(after[refIdx]);
    const alt = clean(after[altIdx]);
    if (!ref && !alt) return null;

    const geneCandidates = toks
        .filter((t, i) => i !== chromIdx && i !== chromIdx + 1)
        .filter((t) => /^[A-Za-z][A-Za-z0-9-]*$/.test(t) && !isChromToken(t));
    const consumed = new Set([after[refIdx], after[altIdx]]);
    const gene = geneCandidates.find((t) => !consumed.has(t) || !isDnaAllele(t)) || null;

    return { chrom, pos, ref, alt, gene: gene ? gene.toUpperCase() : null };
}

function coordinateRowToGenomicHgvs(row) {
    if (!row) return null;
    const { chrom, ref, alt } = row;
    const pos = parseInt(row.pos, 10);
    if (!Number.isFinite(pos)) return null;

    if (ref.length === 1 && alt.length === 1) return `chr${chrom}:g.${pos}${ref}>${alt}`;
    if (!ref) return `chr${chrom}:g.${pos}_${pos + 1}ins${alt}`;
    if (!alt) {
        return ref.length === 1
            ? `chr${chrom}:g.${pos}del`
            : `chr${chrom}:g.${pos}_${pos + ref.length - 1}del`;
    }
    let r = ref;
    let a = alt;
    let start = pos;
    let end = pos + ref.length - 1;
    while (r.length && a.length && r[0] === a[0]) { r = r.slice(1); a = a.slice(1); start += 1; }
    while (r.length && a.length && r[r.length - 1] === a[a.length - 1]) {
        r = r.slice(0, -1); a = a.slice(0, -1); end -= 1;
    }
    if (!r && !a) return null;
    if (!a) {
        return start === end ? `chr${chrom}:g.${start}del` : `chr${chrom}:g.${start}_${end}del`;
    }
    if (!r) {
        return `chr${chrom}:g.${end}_${end + 1}ins${a}`;
    }
    if (r.length === 1 && a.length === 1) return `chr${chrom}:g.${start}${r}>${a}`;
    return start === end
        ? `chr${chrom}:g.${start}delins${a}`
        : `chr${chrom}:g.${start}_${end}delins${a}`;
}

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
