/*
 * Tests for genomic-coordinate parsing and VCF allele resolution.
 *
 * The variant cards need three different things out of a genomic HGVS string:
 *   - a position/span (UCSC link, ClinVar region pull, nearby-variant plot)
 *   - VCF-style anchored alleles (SpliceAI, gnomAD variant pages)
 *   - those alleles left-aligned (gnomAD indexes the 5'-most representation,
 *     while HGVS uses the 3'-most one)
 *
 * Only substitutions used to parse, so anything written as g.START_ENDdel — the
 * form the Ensembl recoder returns for a multi-base deletion such as
 * TSC2 c.2319_2321delAAT — silently disabled all of the above.
 *
 * Helpers are copied verbatim from script.js (the project has no module system)
 * — KEEP IN SYNC when either side changes.
 *
 * Run with: node tests/genomic-coordinates.test.js
 */

// --- helpers copied verbatim from script.js -------------------------------

function parseGenomicHgvs(gVariant) {
    if (!gVariant) return null;
    const m = String(gVariant).trim().match(/^chr([0-9XYMT]+):g\.(\d+)(?:_(\d+))?(.+)$/i);
    if (!m) return null;
    const chrom = `chr${m[1].toUpperCase()}`;
    const start = Number(m[2]);
    const explicitEnd = m[3] ? Number(m[3]) : null;
    const suffix = m[4].trim();
    const mk = (type, end, refSeq, altSeq) => ({
        chrom,
        start,
        end: Number.isFinite(end) ? end : start,
        type,
        refSeq: refSeq || null,
        altSeq: altSeq || null
    });

    // Substitution — including the "REF>-" / "->ALT" forms that the SPDI converter
    // emits for single-base deletions and insertions.
    let s = suffix.match(/^([A-Za-z-]+)>([A-Za-z-]+)$/);
    if (s) {
        const ref = s[1].toUpperCase();
        const alt = s[2].toUpperCase();
        if (ref === '-' && alt !== '-') return mk('ins', start, null, alt);
        if (alt === '-' && ref !== '-') return mk('del', start + ref.length - 1, ref, null);
        if (ref === '-' || alt === '-') return null;
        return mk(ref.length === 1 && alt.length === 1 ? 'sub' : 'delins', explicitEnd ?? (start + ref.length - 1), ref, alt);
    }
    // delins — must be tested before plain "del" so "delinsAA" is not read as a deletion.
    s = suffix.match(/^del([A-Za-z]*)ins([A-Za-z]+)$/i);
    if (s) return mk('delins', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, s[2].toUpperCase());
    s = suffix.match(/^del([A-Za-z]*)$/i);
    if (s) return mk('del', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, null);
    s = suffix.match(/^dup([A-Za-z]*)$/i);
    if (s) return mk('dup', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, null);
    s = suffix.match(/^ins([A-Za-z]+)$/i);
    if (s) return mk('ins', explicitEnd ?? start, null, s[1].toUpperCase());
    s = suffix.match(/^inv([A-Za-z]*)$/i);
    if (s) return mk('inv', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, null);
    return null;
}

function buildVariantCoordinateTuple(rawInput, gVariant) {
    const parseTokenInput = (raw) => {
        if (!raw) return null;
        const toks = String(raw).trim().split(/\s+/).filter(Boolean);
        if (toks.length !== 4 && toks.length !== 5) return null;
        const tokens = toks.slice();
        if (tokens.length === 5) {
            // For 5-token genomic input, always treat token #3 as an optional gene/label
            // and remove it. Example: "17 7573954 TP53 T A".
            tokens.splice(2, 1);
        }
        const [chrTok, posTok, refTok, altTok] = tokens;
        const chrom = String(chrTok).replace(/^chr/i, '').toUpperCase();
        const pos = String(posTok).replace(/,/g, '');
        const ref = String(refTok).toUpperCase();
        const alt = String(altTok).toUpperCase();
        if (!/^[0-9XYMT]+$/.test(chrom)) return null;
        if (!/^\d+$/.test(pos)) return null;
        if (!/^[A-Za-z-]+$/.test(ref) || !/^[A-Za-z-]+$/.test(alt)) return null;
        const startNum = Number(pos);
        return {
            chrom: `chr${chrom}`, pos, ref, alt,
            start: startNum, end: startNum + Math.max(ref.length, 1) - 1,
            type: ref.length === 1 && alt.length === 1 ? 'sub' : 'delins'
        };
    };
    // Parse any genomic HGVS form (sub / del / dup / ins / inv / delins).
    const parseHgvs = (gv) => {
        const desc = parseGenomicHgvs(gv);
        if (!desc) return null;
        return {
            chrom: desc.chrom,
            pos: String(desc.start),
            ref: desc.type === 'sub' ? desc.refSeq : null,
            alt: desc.type === 'sub' ? desc.altSeq : null,
            start: desc.start,
            end: desc.end,
            type: desc.type
        };
    };
    return parseTokenInput(rawInput) || parseHgvs(gVariant) || null;
}

function reverseComplementSeq(seq) {
    const comp = { A: 'T', T: 'A', C: 'G', G: 'C', N: 'N' };
    return String(seq || '').toUpperCase().split('').reverse().map((b) => comp[b] || 'N').join('');
}

function vcfAllelesFromGenomicHgvs(desc, windowStart, seq) {
    if (!desc || !seq) return null;
    const baseAt = (p) => seq[p - windowStart] || '';
    const sliceAt = (a, b) => seq.slice(a - windowStart, b - windowStart + 1);
    const { start, end, type } = desc;
    switch (type) {
        case 'sub': {
            const ref = desc.refSeq || baseAt(start);
            if (!ref || !desc.altSeq) return null;
            return { pos: start, ref, alt: desc.altSeq };
        }
        case 'delins': {
            const ref = sliceAt(start, end) || desc.refSeq;
            if (!ref || !desc.altSeq) return null;
            return { pos: start, ref, alt: desc.altSeq };
        }
        case 'del': {
            const deleted = sliceAt(start, end);
            const anchor = baseAt(start - 1);
            if (!deleted || !anchor) return null;
            return { pos: start - 1, ref: anchor + deleted, alt: anchor };
        }
        case 'dup': {
            // The duplicated copy is inserted immediately after `end`.
            const duplicated = sliceAt(start, end);
            const anchor = baseAt(end);
            if (!duplicated || !anchor) return null;
            return { pos: end, ref: anchor, alt: anchor + duplicated };
        }
        case 'ins': {
            // g.START_ENDins uses START as the base preceding the insertion point.
            const anchor = baseAt(start);
            if (!anchor || !desc.altSeq) return null;
            return { pos: start, ref: anchor, alt: anchor + desc.altSeq };
        }
        case 'inv': {
            const ref = sliceAt(start, end);
            if (!ref) return null;
            return { pos: start, ref, alt: reverseComplementSeq(ref) };
        }
        default:
            return null;
    }
}

function leftNormalizeVcfAlleles(pos, ref, alt, windowStart, seq) {
    let p = Number(pos);
    let r = String(ref || '').toUpperCase();
    let a = String(alt || '').toUpperCase();
    if (!Number.isFinite(p) || !r || !a) return { pos, ref, alt };
    // Bounded so a malformed window can never spin here.
    for (let i = 0; i < 500; i++) {
        if (r.length > 0 && a.length > 0 && r[r.length - 1] === a[a.length - 1]) {
            r = r.slice(0, -1);
            a = a.slice(0, -1);
            continue;
        }
        if (r.length === 0 || a.length === 0) {
            const prevBase = seq && seq[p - 1 - windowStart];
            if (!prevBase) break;
            r = prevBase + r;
            a = prevBase + a;
            p -= 1;
            continue;
        }
        break;
    }
    // A fully consumed allele means we walked off the window; keep the input.
    if (!r || !a) return { pos, ref, alt };
    return { pos: p, ref: r, alt: a };
}

// --- test harness ---------------------------------------------------------

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) {
        passed++;
    } else {
        failed++;
        console.error(`FAIL: ${name}`);
    }
}
const eq = (a, b) => JSON.stringify(a) === JSON.stringify(b);

// --- parseGenomicHgvs -----------------------------------------------------

check('substitution', eq(parseGenomicHgvs('chr7:g.140453136A>T'),
    { chrom: 'chr7', start: 140453136, end: 140453136, type: 'sub', refSeq: 'A', altSeq: 'T' }));
check('multi-base deletion keeps its span', eq(parseGenomicHgvs('chr16:g.2122948_2122950del'),
    { chrom: 'chr16', start: 2122948, end: 2122950, type: 'del', refSeq: null, altSeq: null }));
check('deletion with spelled-out bases', eq(parseGenomicHgvs('chr16:g.2122948_2122950delAAT'),
    { chrom: 'chr16', start: 2122948, end: 2122950, type: 'del', refSeq: 'AAT', altSeq: null }));
check('single-base deletion', eq(parseGenomicHgvs('chr16:g.2122948del'),
    { chrom: 'chr16', start: 2122948, end: 2122948, type: 'del', refSeq: null, altSeq: null }));
check('SPDI-style "REF>-" reads as a deletion, not a substitution',
    eq(parseGenomicHgvs('chr16:g.2122948A>-'),
        { chrom: 'chr16', start: 2122948, end: 2122948, type: 'del', refSeq: 'A', altSeq: null }));
check('SPDI-style "->ALT" reads as an insertion', eq(parseGenomicHgvs('chr16:g.2122948->AT'),
    { chrom: 'chr16', start: 2122948, end: 2122948, type: 'ins', refSeq: null, altSeq: 'AT' }));
check('duplication', eq(parseGenomicHgvs('chr2:g.29443695_29443696dup'),
    { chrom: 'chr2', start: 29443695, end: 29443696, type: 'dup', refSeq: null, altSeq: null }));
check('insertion', eq(parseGenomicHgvs('chr4:g.55593603_55593604insTT'),
    { chrom: 'chr4', start: 55593603, end: 55593604, type: 'ins', refSeq: null, altSeq: 'TT' }));
check('inversion', eq(parseGenomicHgvs('chr3:g.10183694_10183700inv'),
    { chrom: 'chr3', start: 10183694, end: 10183700, type: 'inv', refSeq: null, altSeq: null }));
check('delins is not read as a deletion', eq(parseGenomicHgvs('chr7:g.140453136_140453137delinsAA'),
    { chrom: 'chr7', start: 140453136, end: 140453137, type: 'delins', refSeq: null, altSeq: 'AA' }));
check('delins with spelled-out reference', eq(parseGenomicHgvs('chr7:g.140453136_140453137delTGinsAA'),
    { chrom: 'chr7', start: 140453136, end: 140453137, type: 'delins', refSeq: 'TG', altSeq: 'AA' }));
check('MNV substitution is a delins', eq(parseGenomicHgvs('chr7:g.140453136TG>AA'),
    { chrom: 'chr7', start: 140453136, end: 140453137, type: 'delins', refSeq: 'TG', altSeq: 'AA' }));
check('non-genomic input rejected', parseGenomicHgvs('TSC2:c.2319_2321delAAT') === null);
check('empty input rejected', parseGenomicHgvs('') === null);

// --- buildVariantCoordinateTuple ------------------------------------------

check('raw token input still wins over the g. notation',
    eq(buildVariantCoordinateTuple('17 7573954 TP53 T A', 'chr17:g.7573954T>A'),
        { chrom: 'chr17', pos: '7573954', ref: 'T', alt: 'A', start: 7573954, end: 7573954, type: 'sub' }));
{
    const t = buildVariantCoordinateTuple('tsc2 c.2319_2321delAAT', 'chr16:g.2122948_2122950del');
    check('deletion yields coordinates for position-only consumers',
        t && t.chrom === 'chr16' && t.pos === '2122948' && t.start === 2122948 && t.end === 2122950);
    check('deletion reports no alleles rather than guessing', t && t.ref === null && t.alt === null);
}
check('substitution still yields alleles',
    eq(buildVariantCoordinateTuple('BRAF V600E', 'chr7:g.140453136A>T'),
        { chrom: 'chr7', pos: '140453136', ref: 'A', alt: 'T', start: 140453136, end: 140453136, type: 'sub' }));
check('unparseable input yields null', buildVariantCoordinateTuple('BRAF V600E', 'BRAF:p.Val600Glu') === null);

// --- VCF allele construction + left alignment -----------------------------

// GRCh37 chr16:2122930-2122960 (TSC2), fetched from Ensembl.
const TSC2_START = 2122930;
const TSC2_SEQ = 'TCCAGTGCTGACAGCATTAATCTCTTACCAT';
// GRCh37 chr7:117199630-117199660 (CFTR), fetched from Ensembl.
const CFTR_START = 117199630;
const CFTR_SEQ = 'TTAAAGAAAATATCATCTTTGGTGTTTCCTA';

const resolve = (hgvs, windowStart, seq) => {
    const anchored = vcfAllelesFromGenomicHgvs(parseGenomicHgvs(hgvs), windowStart, seq);
    return leftNormalizeVcfAlleles(anchored.pos, anchored.ref, anchored.alt, windowStart, seq);
};

check('TSC2 c.2319_2321delAAT anchors on the preceding base',
    eq(vcfAllelesFromGenomicHgvs(parseGenomicHgvs('chr16:g.2122948_2122950del'), TSC2_START, TSC2_SEQ),
        { pos: 2122947, ref: 'TAAT', alt: 'T' }));
check('TSC2 c.2319_2321delAAT left-aligns to 16-2122946-TTAA-T',
    eq(resolve('chr16:g.2122948_2122950del', TSC2_START, TSC2_SEQ),
        { pos: 2122946, ref: 'TTAA', alt: 'T' }));
// The representation gnomAD actually indexes for CFTR F508del is 7-117199644-ATCT-A
// (GRCh38 7-117559590-ATCT-A), not the 3'-shifted HGVS coordinates.
check('CFTR F508del left-aligns to the representation gnomAD indexes',
    eq(resolve('chr7:g.117199646_117199648del', CFTR_START, CFTR_SEQ),
        { pos: 117199644, ref: 'ATCT', alt: 'A' }));
check('substitution passes through unchanged',
    eq(resolve('chr16:g.2122948A>C', TSC2_START, TSC2_SEQ), { pos: 2122948, ref: 'A', alt: 'C' }));
check('duplication anchors on the last duplicated base',
    eq(vcfAllelesFromGenomicHgvs(parseGenomicHgvs('chr16:g.2122948_2122949dup'), TSC2_START, TSC2_SEQ),
        { pos: 2122949, ref: 'A', alt: 'AAA' }));
check('insertion anchors on the base before the insertion point',
    eq(vcfAllelesFromGenomicHgvs(parseGenomicHgvs('chr16:g.2122948_2122949insGG'), TSC2_START, TSC2_SEQ),
        { pos: 2122948, ref: 'A', alt: 'AGG' }));
check('delins uses the deleted span as REF',
    eq(vcfAllelesFromGenomicHgvs(parseGenomicHgvs('chr16:g.2122948_2122950delinsG'), TSC2_START, TSC2_SEQ),
        { pos: 2122948, ref: 'AAT', alt: 'G' }));
check('inversion reverse-complements the span',
    eq(vcfAllelesFromGenomicHgvs(parseGenomicHgvs('chr16:g.2122948_2122950inv'), TSC2_START, TSC2_SEQ),
        { pos: 2122948, ref: 'AAT', alt: 'ATT' }));
check('left alignment stops at the edge of the fetched window rather than looping',
    (() => {
        const out = leftNormalizeVcfAlleles(TSC2_START, 'T', '', TSC2_START, TSC2_SEQ);
        return out.pos === TSC2_START && out.ref === 'T' && out.alt === '';
    })());

console.log(`\nGenomic-coordinate tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
