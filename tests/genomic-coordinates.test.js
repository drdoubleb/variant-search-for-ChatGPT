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
    const m = String(gVariant).trim().match(/^chr([0-9XYMT]+):g\.(\d+)(?:_(\d+))?(.*)$/i);
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

    // A bare span with no event ("chr16:g.2136834_2136836") is a region, not a
    // variant. Protein-only queries that cannot be resolved to a nucleotide change
    // are represented this way so position-based cards still work while the
    // allele-dependent ones correctly stay empty.
    if (suffix === '') return mk('region', explicitEnd ?? start, null, null);

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

const COORDINATE_ROW_NOISE = /^(hg18|hg19|hg38|grch3[678]|b3[678]|chr|snp|snv|ins|del|indel|mnp|somatic|germline|het|hom|\+|\.)$/i;

const isDnaAllele = (t) => /^(?:[ACGTN]+|-)$/i.test(String(t));
const isChromToken = (t) => /^(?:chr)?(?:[0-9]{1,2}|X|Y|M|MT)$/i.test(String(t));

function splitCoordinateRow(raw) {
    return String(raw || '')
        .trim()
        .replace(/\b\d{1,3}(?:,\d{3})+(?!\d)/g, (m) => m.replace(/,/g, ''))
        .split(/[\s,;|]+/)
        .map((t) => t.replace(/^["']+|["']+$/g, ''))
        .filter(Boolean);
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

function buildVariantCoordinateTuple(rawInput, gVariant) {
    const parseTokenInput = (raw) => {
        const row = parseCoordinateRow(raw);
        if (!row) return null;
        let startNum = Number(row.pos);
        if (!Number.isFinite(startNum)) return null;
        let ref = row.ref;
        let alt = row.alt;
        // Same-length pairs may carry shared context bases ("CT">"CA" is really a
        // T>A one base along). VCF SNVs/MNVs need no anchor base, so trim to the
        // minimal representation — the one gnomAD and SpliceAI index. Indels keep
        // their anchored VCF form untouched.
        if (ref.length === alt.length && ref.length > 1) {
            while (ref.length > 1 && ref[0] === alt[0]) {
                ref = ref.slice(1); alt = alt.slice(1); startNum += 1;
            }
            while (ref.length > 1 && ref[ref.length - 1] === alt[alt.length - 1]) {
                ref = ref.slice(0, -1); alt = alt.slice(0, -1);
            }
        }
        return {
            chrom: `chr${row.chrom}`, pos: String(startNum), ref, alt,
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

// Build a GRCh37 genomic HGVS string from a VEP result object.
//
// VEP resolves gene-level cDNA HGVS ("TSC2:c.4952delA") that the variant recoder
// can miss, and its response already carries the genomic location — but the app
// used to read only the consequence off it, leaving the g. field showing the
// user's own query and every coordinate-derived card dark.
//
// VEP's coordinate conventions: a deletion's start..end span the deleted bases; an
// insertion is reported with start = end + 1 and a "-" reference.
function buildGenomicHgvsFromVep(vepResult) {
    if (!vepResult) return '';
    const chrom = String(vepResult.seq_region_name || '').replace(/^chr/i, '').toUpperCase();
    const start = Number(vepResult.start);
    const end = Number(vepResult.end);
    const alleles = String(vepResult.allele_string || '').split('/');
    if (!/^[0-9XYMT]+$/.test(chrom) || !Number.isFinite(start) || !Number.isFinite(end)) return '';
    const ref = String(alleles[0] || '').toUpperCase();
    const alt = String(alleles[1] || '').toUpperCase();
    if (!ref || !alt) return '';
    if (ref === '-') {
        // Insertion between `end` and `start` (VEP reports start = end + 1).
        if (!/^[ACGTN]+$/.test(alt)) return '';
        return `chr${chrom}:g.${end}_${start}ins${alt}`;
    }
    if (alt === '-') {
        return start === end
            ? `chr${chrom}:g.${start}del`
            : `chr${chrom}:g.${start}_${end}del`;
    }
    if (!/^[ACGTN]+$/.test(ref) || !/^[ACGTN]+$/.test(alt)) return '';
    if (ref.length === 1 && alt.length === 1) return `chr${chrom}:g.${start}${ref}>${alt}`;
    return start === end
        ? `chr${chrom}:g.${start}delins${alt}`
        : `chr${chrom}:g.${start}_${end}delins${alt}`;
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

// Parse the affected residue range out of a protein HGVS body.
//
// Handles single- and three-letter codes across substitutions, nonsense,
// frameshifts, in-frame del/dup/ins/delins and extensions:
//   p.N1651Mfs*21             -> { start: 1651, end: 1651 }
//   p.Asn1651MetfsTer21       -> { start: 1651, end: 1651 }
//   p.L773_I774delinsF        -> { start: 773,  end: 774  }
//   p.Leu773_Ile774delinsPhe  -> { start: 773,  end: 774  }
//   p.Lys18del                -> { start: 18,   end: 18   }
// Returns null when no residue position can be read.
function parseProteinResidueRange(proteinChange) {
    const body = String(proteinChange || '').trim().replace(/^p\./i, '').replace(/[()]/g, '');
    if (!body) return null;
    // A range is written "<aa><pos>_<aa><pos>"; take the first and last positions.
    const positions = body.match(/(?:[A-Za-z]{3}|[A-Za-z])(\d+)/g);
    if (!positions || positions.length === 0) return null;
    const nums = positions
        .map((tok) => Number((tok.match(/(\d+)$/) || [])[1]))
        .filter((n) => Number.isFinite(n) && n > 0);
    if (nums.length === 0) return null;
    // "fs*21" / "extTer5" trailing counts are lengths, not residue positions, and the
    // regex above cannot match them (they have no preceding amino-acid token), so the
    // remaining numbers are all genuine residue coordinates.
    return { start: Math.min(...nums), end: Math.max(...nums) };
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
check('bare span is a region, not a variant', eq(parseGenomicHgvs('chr16:g.2136834_2136836'),
    { chrom: 'chr16', start: 2136834, end: 2136836, type: 'region', refSeq: null, altSeq: null }));
check('bare position is a region', eq(parseGenomicHgvs('chr16:g.2136835'),
    { chrom: 'chr16', start: 2136835, end: 2136835, type: 'region', refSeq: null, altSeq: null }));
check('a region yields no VCF alleles',
    vcfAllelesFromGenomicHgvs(parseGenomicHgvs('chr16:g.2136834_2136836'), 2136830, 'ACGTACGTACGT') === null);
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
// Same-length pairs with shared context bases reduce to the minimal alleles
// (matching the minimal g. notation the row parser now produces).
check('same-length pair with shared prefix trims to the SNV',
    eq(buildVariantCoordinateTuple('7 140453135 BRAF CT CA', null),
        { chrom: 'chr7', pos: '140453136', ref: 'T', alt: 'A', start: 140453136, end: 140453136, type: 'sub' }));
check('same-length pair with shared flanks trims to the SNV',
    eq(buildVariantCoordinateTuple('7 140453135 BRAF CAG CTG', null),
        { chrom: 'chr7', pos: '140453136', ref: 'A', alt: 'T', start: 140453136, end: 140453136, type: 'sub' }));
{
    const t = buildVariantCoordinateTuple('7 140453136 BRAF AC TG', null);
    check('true MNV keeps both alleles untrimmed',
        t && t.ref === 'AC' && t.alt === 'TG' && t.pos === '140453136' && t.type === 'delins');
}
{
    // VCF-anchored indel: the anchored form is what gnomAD/SpliceAI want — untouched.
    const t = buildVariantCoordinateTuple('X 20148674 EIF1AX T TG', null);
    check('anchored indel alleles are not trimmed',
        t && t.ref === 'T' && t.alt === 'TG' && t.pos === '20148674');
}

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
