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
 * Exercises the REAL implementations, imported from script.js via its Node
 * test-export block.
 *
 * Run with: node tests/genomic-coordinates.test.js
 */

// --- real helpers imported from script.js ---------------------------------

await import('../script.js');
const {
    GRCH37_NC_VERSIONS, assemblyFromNcAccession, splitGenomicHgvsPositions,
    ncGenomicToChr, minimalSpdiForms, parseGenomicHgvs,
    buildVariantCoordinateTuple, vcfAllelesFromGenomicHgvs, leftNormalizeVcfAlleles
} = globalThis.__variantSearchHelpers;

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

// --- assemblyFromNcAccession ----------------------------------------------
// The GRCh37 recoder mirror emits GRCh37 accessions (NC_000007.13); the main
// host emits GRCh38 (NC_000007.14). Candidate conversion must lift ONLY the
// GRCh38 ones — blanket-lifting corrupted already-hg19 substitutions, because
// Ensembl's /map answers wrongly (not with an error) for a coordinate from the
// wrong assembly.

check('GRCh37 hgvsg accession detected', assemblyFromNcAccession('NC_000007.13:g.140453136A>T') === 'GRCh37');
check('GRCh38 hgvsg accession detected', assemblyFromNcAccession('NC_000007.14:g.140753336A>T') === 'GRCh38');
check('GRCh37 SPDI accession detected', assemblyFromNcAccession('NC_000016.9:2122947:AAT:') === 'GRCh37');
check('GRCh38 chrX accession detected', assemblyFromNcAccession('NC_000023.11:g.20148674T>G') === 'GRCh38');
check('GRCh37 chrX accession detected', assemblyFromNcAccession('NC_000023.10:g.20148674T>G') === 'GRCh37');
check('unknown version yields null', assemblyFromNcAccession('NC_000007.99:g.1A>T') === null);
check('mitochondrial accession yields null', assemblyFromNcAccession('NC_012920.1:g.100A>T') === null);
check('non-NC input yields null', assemblyFromNcAccession('chr7:g.140453136A>T') === null);
check('empty input yields null', assemblyFromNcAccession('') === null);

// --- splitGenomicHgvsPositions --------------------------------------------
// The GRCh38 input toggle remaps ONLY the positions of a genomic HGVS string,
// leaving the event description ("rest") byte-for-byte intact — so it must
// split every g. form the app accepts, not just substitutions.

check('split: substitution', eq(splitGenomicHgvsPositions('chr7:g.140753336A>T'),
    { chrom: '7', start: 140753336, end: null, rest: 'A>T' }));
check('split: chr prefix optional', eq(splitGenomicHgvsPositions('7:g.140753336A>T'),
    { chrom: '7', start: 140753336, end: null, rest: 'A>T' }));
check('split: range deletion', eq(splitGenomicHgvsPositions('chr16:g.2110671_2110673del'),
    { chrom: '16', start: 2110671, end: 2110673, rest: 'del' }));
check('split: insertion keeps inserted bases in rest', eq(splitGenomicHgvsPositions('chrX:g.20166791_20166792insG'),
    { chrom: 'X', start: 20166791, end: 20166792, rest: 'insG' }));
check('split: single-position del', eq(splitGenomicHgvsPositions('chr17:g.7674220del'),
    { chrom: '17', start: 7674220, end: null, rest: 'del' }));
check('split: delins', eq(splitGenomicHgvsPositions('chr7:g.140753336_140753337delinsAA'),
    { chrom: '7', start: 140753336, end: 140753337, rest: 'delinsAA' }));
check('split: mitochondrial MT accepted', eq(splitGenomicHgvsPositions('chrMT:g.8993T>G'),
    { chrom: 'MT', start: 8993, end: null, rest: 'T>G' }));
check('split: NC_ form rejected (convert first)', splitGenomicHgvsPositions('NC_000007.14:g.140753336A>T') === null);
check('split: non-genomic rejected', splitGenomicHgvsPositions('BRAF:p.Val600Glu') === null);
check('split: empty rejected', splitGenomicHgvsPositions('') === null);

// --- ncGenomicToChr --------------------------------------------------------

check('NC_ GRCh38 chr7 converts to chr form', ncGenomicToChr('NC_000007.14:g.140753336A>T') === 'chr7:g.140753336A>T');
check('NC_ GRCh37 chr7 converts to chr form', ncGenomicToChr('NC_000007.13:g.140453136A>T') === 'chr7:g.140453136A>T');
check('NC_ chrX converts', ncGenomicToChr('NC_000023.10:g.20148674T>G') === 'chrX:g.20148674T>G');
check('NC_ chrY converts', ncGenomicToChr('NC_000024.9:g.100A>T') === 'chrY:g.100A>T');
check('mitochondrial NC_ left for the recoder', ncGenomicToChr('NC_012920.1:g.100A>T') === null);
check('non-genomic NC_ (SPDI) rejected', ncGenomicToChr('NC_000016.9:2122947:AAT:') === null);
check('chr-form input passes through as null', ncGenomicToChr('chr7:g.140453136A>T') === null);

// --- minimalSpdiForms -------------------------------------------------------
// NCBI Variation Services contextual SPDIs are VCF-anchored; these are real
// responses fetched live. MyVariant indexes only the minimal event, and in
// repeat regions the trim order decides which flank survives — so both minimal
// forms are offered as candidates.

check('SPDI: clean substitution passes through as one form',
    eq(minimalSpdiForms({ seq_id: 'NC_000007.13', position: 140453135, deleted_sequence: 'A', inserted_sequence: 'T' }),
        ['NC_000007.13:140453135:A:T']));
check('SPDI: anchored deletion TAAT→T yields both trim orders',
    eq(minimalSpdiForms({ seq_id: 'NC_000016.9', position: 2122946, deleted_sequence: 'TAAT', inserted_sequence: 'T' }),
        ['NC_000016.9:2122947:AAT:', 'NC_000016.9:2122946:TAA:']));
check('SPDI: anchored duplication G→GG yields both insertion placements',
    eq(minimalSpdiForms({ seq_id: 'NC_000023.10', position: 20148674, deleted_sequence: 'G', inserted_sequence: 'GG' }),
        ['NC_000023.10:20148675::G', 'NC_000023.10:20148674::G']));
check('SPDI: anchored insertion T→TTT',
    eq(minimalSpdiForms({ seq_id: 'NC_000004.11', position: 55593602, deleted_sequence: 'T', inserted_sequence: 'TTT' }),
        ['NC_000004.11:55593603::TT', 'NC_000004.11:55593602::TT']));
check('SPDI: identical sequences describe no variant',
    eq(minimalSpdiForms({ seq_id: 'NC_000007.13', position: 1, deleted_sequence: 'A', inserted_sequence: 'A' }), []));
check('SPDI: MNV with shared flanks trims to the core',
    eq(minimalSpdiForms({ seq_id: 'NC_000007.13', position: 140453134, deleted_sequence: 'CAG', inserted_sequence: 'CTG' }),
        ['NC_000007.13:140453135:A:T']));

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
