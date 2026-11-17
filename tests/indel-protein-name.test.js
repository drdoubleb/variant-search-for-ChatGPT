/*
 * One protein name for an insertion or deletion, however it is spelled.
 *
 * An indel in a repeat can be written anywhere along it. Ensembl's VEP gives
 * every spelling the same coding name, but not the same protein name: it was
 * asked about some 500 spellings of 50 indels, one by one, and
 *
 *   - every frameshift, and every indel in a gene on the forward strand, had
 *     one protein name at all of its spellings;
 *   - an in-frame insertion, duplication or deletion in a gene on the reverse
 *     strand was misnamed at some spellings in 13 of 27 cases. BRAF
 *     c.1794_1796dup — p.Thr599dup — is "p.Thr599_Val600insTer", an inserted
 *     stop, when written at its 3'-most genomic position. The duplication of
 *     NRAS codon 12 is "p.Val14delinsGlyPhe" in the very spelling Ensembl's
 *     recoder gives it, which is the spelling this page asked VEP by: a
 *     healthy lookup of "NRAS c.34_36dup" showed that name. "KRAS c.34_36dup"
 *     came out as p.Val14dup when neither Ensembl's recoder nor NCBI had
 *     placed the name;
 *   - on the reverse strand the left-aligned spelling was named right every
 *     time.
 *
 * So the protein names VEP gives on reverse-strand transcripts are taken from
 * its answer for the left-aligned spelling, whichever spelling the lookup
 * arrived at: a row that answer names replaces the row in hand. Where no name
 * can be had that way (the reference around the indel is unavailable, VEP does
 * not answer), the protein names VEP gives on reverse-strand transcripts are
 * not shown, and the panel says why. Everything else is VEP's for the spelling
 * in hand, as it always was — forward-strand rows (there the left-aligned
 * spelling is the wrong one to ask by: scenario 7) and reverse-strand rows the
 * left-aligned answer has no protein name for (scenario 8).
 *
 * The scenarios replay answers recorded from the live services.
 *
 * Run with: node tests/indel-protein-name.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved, referenceSlice } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const vepAsked = page => page.requests.filter(r => r.includes('/vep/human/hgvs/')).map(r => r.match(/hgvs\/([^?]+)/)[1]);
// The spellings VEP was asked about the allele by, in order, each once. (It is
// asked twice per spelling when its RefSeq set has no row.)
const allelesAsked = page => Array.from(new Set(vepAsked(page).filter(asked => /:g\./.test(asked))));
const withheldRemark = /protein change not shown: Ensembl names an in-frame insertion or deletion on the reverse strand reliably only when it is written left-aligned, and no name could be had for it written that way/;
const askedAsRemark = spelling => new RegExp(`reverse-strand transcripts asked as ${spelling.replace(/[.]/g, '\\.')} — the same allele, left-aligned: `
    + 'Ensembl\'s protein name for an in-frame insertion or deletion there depends on where it is written');
const wait = ms => new Promise(resolve => setTimeout(resolve, ms));

// ── Which alleles this is about ────────────────────────────────────────────
for (const inFrame of ['chr12:g.25398284_25398286dup', 'chr12:g.25398280_25398281insCCA', 'chr7:g.117199646_117199648del',
    'chrX:g.20156694_20156723dup', 'chr7:g.117199646_117199648delCTT', 'chr12:g.25398280_25398281insCCACCA',
    // ...and a replacement that changes the length by whole codons: it may be one of those in another notation.
    'chr7:g.140453142delinsGTAG', 'chr12:g.25398281_25398284delinsC', 'chr12:g.25398284_25398286delinsTTTAAA']) {
    assert.equal(h.inFrameIndel(inFrame), true, inFrame);
}
for (const other of ['chr17:g.41209082dup', 'chr13:g.32906698_32906699insCA', 'chr17:g.41276045_41276046del', 'chr12:g.25398284C>T',
    'chr12:g.25398284_25398286delinsTTT', 'chr12:g.25398284delinsTT', 'chr12:g.25398284_25398286', 'chr12:g.25398280_25398281ins',
    'NM_004985.5:c.34_36dup', 'chr12:g.25398284_25398287dup', '', null, undefined]) {
    assert.equal(h.inFrameIndel(other), false, String(other));
}

// ── The left-aligned spelling ──────────────────────────────────────────────
// KRAS, reverse strand: …GCC ACC AGC… around codon 12. The duplication of
// codon 12 can be written at seven places, or as a duplication.
const kras = cassettes.krasDup.reference.find(w => w.chrom === '12');
const krasSpellings = ['chr12:g.25398280_25398281insCCA', 'chr12:g.25398281_25398282insCAC', 'chr12:g.25398282_25398283insACC',
    'chr12:g.25398283_25398284insCCA', 'chr12:g.25398284_25398285insCAC', 'chr12:g.25398285_25398286insACC',
    'chr12:g.25398286_25398287insCCA', 'chr12:g.25398284_25398286dup', 'chr12:g.25398283_25398285dup'];
assert.equal(h.leftAlignedSpelling(krasSpellings[0], kras.start, kras.seq), '', 'left-aligned already: asked as it stands');
for (const spelling of krasSpellings.slice(1)) {
    assert.equal(h.leftAlignedSpelling(spelling, kras.start, kras.seq), 'chr12:g.25398280_25398281insCCA', spelling);
    assert.equal(h.sameGenomicAllele(spelling, 'chr12:g.25398280_25398281insCCA', kras.start, kras.seq), true, `${spelling}: one allele`);
}
// A deletion: KRAS c.34_36del can be written at four places.
for (const spelling of ['chr12:g.25398282_25398284del', 'chr12:g.25398283_25398285del', 'chr12:g.25398284_25398286del',
    'chr12:g.25398284_25398286delCCA']) {
    assert.equal(h.leftAlignedSpelling(spelling, kras.start, kras.seq), 'chr12:g.25398281_25398283del', spelling);
}
assert.equal(h.leftAlignedSpelling('chr12:g.25398281_25398283del', kras.start, kras.seq), '');
// One that cannot slide is where it is.
assert.equal(h.leftAlignedSpelling('chr12:g.25398285_25398286insGGT', kras.start, kras.seq), '');
// Every indel, not only in-frame ones (which of them matter is inFrameIndel's
// to say).
assert.equal(h.leftAlignedSpelling('chr12:g.25398282dup', kras.start, kras.seq), 'chr12:g.25398280_25398281insC');
assert.equal(h.leftAlignedSpelling('chr12:g.25398282del', kras.start, kras.seq), 'chr12:g.25398281del');
// What the window cannot say is not guessed.
assert.equal(h.leftAlignedSpelling('chr12:g.25398284C>T', kras.start, kras.seq), null, 'not an indel');
assert.equal(h.leftAlignedSpelling('chr12:g.25398284_25398286delAAA', kras.start, kras.seq), null, 'bases the reference does not hold');
assert.equal(h.leftAlignedSpelling('chr12:g.25398284_25398286dup', 25398283, referenceSlice([kras], '12', 25398283, 25398300)), null,
    'slides off the window: the left end is not known');
assert.equal(h.leftAlignedSpelling('chr12:g.99998284_99998286dup', kras.start, kras.seq), null, 'outside the window');
assert.equal(h.leftAlignedSpelling('chr12:g.25398284_25398286dup', kras.start, ''), null);

// An insertion names two neighbouring bases to go between. One that names
// bases five apart, or one base, names no place: Ensembl refuses both, and so
// must this. Read as "after the first", they were re-spelled and shown as an
// allele nobody had named. (Found in review.)
for (const malformed of ['chr12:g.25398280_25398285insCCA', 'chr12:g.25398280insCCA', 'chr12:g.25398280_25398280insCCA',
    'chr12:g.25398281_25398280insCCA']) {
    const desc = h.parseGenomicHgvs(malformed);
    if (desc) assert.equal(h.vcfAllelesFromGenomicHgvs(desc, kras.start, kras.seq), null, malformed);
    assert.equal(h.indelPlacements(malformed, kras.start, kras.seq), null, malformed);
    assert.equal(h.leftAlignedSpelling(malformed, kras.start, kras.seq), null, malformed);
    assert.deepEqual(h.siblingIdCandidates(malformed, kras.start, kras.seq), [], malformed);
    assert.notEqual(h.sameGenomicAllele(malformed, 'chr12:g.25398280_25398281insCCA', kras.start, kras.seq), true, malformed);
}
assert.deepEqual(h.vcfAllelesFromGenomicHgvs(h.parseGenomicHgvs('chr12:g.25398280_25398281insCCA'), kras.start, kras.seq),
    { pos: 25398280, ref: 'G', alt: 'GCCA' });

// An insertion or a deletion written as a replacement. BRAF c.1794_1796dup as
// "g.140453142delinsGTAG" keeps the G it "deletes" and inserts TAG after it —
// the 3'-most spelling, which VEP calls an inserted stop.
const braf = cassettes.brafT599dupRight.reference.find(w => w.chrom === '7' && w.start <= 140452992 && w.start + w.seq.length > 140453292);
assert(braf, 'the recording holds the window around BRAF codon 599');
for (const spelling of ['chr7:g.140453142delinsGTAG', 'chr7:g.140453142_140453143insTAG', 'chr7:g.140453140_140453142dup',
    'chr7:g.140453141_140453142delinsAGTAG']) {
    assert.equal(h.leftAlignedSpelling(spelling, braf.start, braf.seq), 'chr7:g.140453138_140453139insGTA', spelling);
    assert.equal(h.sameGenomicAllele(spelling, 'chr7:g.140453138_140453139insGTA', braf.start, braf.seq), true, `${spelling}: one allele`);
}
// A deletion in the same disguise: KRAS c.34_36del as "CCAC replaced by C".
assert.equal(referenceSlice([kras], '12', 25398281, 25398284), 'CCAC');
assert.equal(h.leftAlignedSpelling('chr12:g.25398281_25398284delinsC', kras.start, kras.seq), 'chr12:g.25398281_25398283del');
assert.equal(h.leftAlignedSpelling('chr12:g.25398281_25398284delCCACinsC', kras.start, kras.seq), 'chr12:g.25398281_25398283del');
// A replacement in earnest has the one spelling: nothing to re-spell.
assert.equal(h.leftAlignedSpelling('chr12:g.25398284_25398286delinsTTTAAA', kras.start, kras.seq), '');
assert.equal(h.leftAlignedSpelling('chr12:g.25398283_25398285delinsTTA', kras.start, kras.seq), '');
// ...and one whose stated bases are not the reference's says nothing.
assert.equal(h.leftAlignedSpelling('chr12:g.25398281_25398284delAAAAinsC', kras.start, kras.seq), null);

// ── What VEP is asked by ───────────────────────────────────────────────────
{
    const calls = [];
    const reference = async (chrom, start, end) => {
        calls.push([chrom, start, end]);
        return referenceSlice([kras, braf], String(chrom).replace(/^chr/, ''), start, end);
    };
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286dup', reference),
        { askedBy: 'chr12:g.25398280_25398281insCCA', checked: true });
    // (The window the lookup reads for the allele's other records: one request.)
    assert.deepEqual(calls, [['chr12', 25398134, 25398436]]);
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398280_25398281insCCA', reference),
        { askedBy: 'chr12:g.25398280_25398281insCCA', checked: true }, 'left-aligned already');
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398283_25398285del', reference),
        { askedBy: 'chr12:g.25398281_25398283del', checked: true });
    assert.deepEqual(await h.vepSpellingFor('chr7:g.140453142delinsGTAG', reference),
        { askedBy: 'chr7:g.140453138_140453139insGTA', checked: true }, 'an insertion written as a replacement');
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286delinsTTTAAA', reference),
        { askedBy: 'chr12:g.25398284_25398286delinsTTTAAA', checked: true }, 'a replacement in earnest');

    // A frameshift, a substitution: asked as spelled, and the reference is not read for it.
    calls.length = 0;
    for (const asIs of ['chr12:g.25398282dup', 'chr12:g.25398282_25398283del', 'chr12:g.25398284C>T', 'chr12:g.25398284_25398286delinsTTT']) {
        assert.deepEqual(await h.vepSpellingFor(asIs, reference), { askedBy: asIs, checked: true }, asIs);
    }
    assert.deepEqual(calls, [], 'nothing read where the spelling does not matter');

    // The reference unavailable: asked as spelled, and said to be unchecked.
    const warn = console.warn;
    console.warn = () => {};
    try {
        assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286dup', async () => { throw new Error('503'); }),
            { askedBy: 'chr12:g.25398284_25398286dup', checked: false });
        // (A replacement is only looked into: with nothing to go by it is left as it was.)
        assert.deepEqual(await h.vepSpellingFor('chr7:g.140453142delinsGTAG', async () => { throw new Error('503'); }),
            { askedBy: 'chr7:g.140453142delinsGTAG', checked: true });
    } finally {
        console.warn = warn;
    }
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286dup', async () => null),
        { askedBy: 'chr12:g.25398284_25398286dup', checked: false });
    // ...as it is for a spelling that does not fit the reference read,
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286delAAA', reference),
        { askedBy: 'chr12:g.25398284_25398286delAAA', checked: false });
    // ...for an insertion that names no place,
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398280_25398285insCCA', reference),
        { askedBy: 'chr12:g.25398280_25398285insCCA', checked: false });
    // ...and for an event wider than the window is read for.
    calls.length = 0;
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25399786del', reference),
        { askedBy: 'chr12:g.25398284_25399786del', checked: false });
    assert.deepEqual(calls, []);
    assert.deepEqual(h.uncheckedVepSpelling('chr12:g.25398284_25398286dup'), { askedBy: 'chr12:g.25398284_25398286dup', checked: false });
    assert.deepEqual(h.uncheckedVepSpelling('chr7:g.140453142delinsGTAG'), { askedBy: 'chr7:g.140453142delinsGTAG', checked: true });
}

// ── Whose rows are taken from which answer ─────────────────────────────────
{
    const rows = h.transcriptsFromVep([
        { transcript_id: 'NM_004985.5', hgvsc: 'NM_004985.5:c.36_38dup', hgvsp: 'NP_004976.2:p.Val14dup', strand: -1, gene_symbol: 'KRAS',
            consequence_terms: ['inframe_insertion'] },
        { transcript_id: 'NM_000492.4', hgvsc: 'NM_000492.4:c.1521_1523del', hgvsp: 'NP_000483.3:p.Phe508del', strand: 1, gene_symbol: 'CFTR' },
        { transcript_id: 'NR_000001.1', hgvsc: 'NR_000001.1:n.5_7dup', strand: -1 }
    ]);
    assert.deepEqual(rows.map(row => row.strand), [-1, 1, -1], 'a row says which strand its transcript is on');
    // (snpEff's name for the allele is its own work, whatever strand it is on.)
    const snpeff = { transcript: 'NM_004985.5', cDNA: 'c.36_38dup', protein: 'p.Gly13dup', source: 'snpeff', strand: -1 };
    const kept = h.withoutUncheckedProteinNames(rows.concat(snpeff));
    assert.deepEqual(kept.map(row => row.protein), ['', 'p.Phe508del', '', 'p.Gly13dup'],
        'reverse-strand names VEP gave go; forward-strand ones and another source\'s stay');
    assert.deepEqual(kept.map(row => row.cDNA), ['c.36_38dup', 'c.1521_1523del', 'n.5_7dup', 'c.36_38dup'], 'the coding name stays');
    assert.deepEqual(kept.map(row => Boolean(row.proteinWithheld)), [true, false, false, false], 'only a name taken away is flagged');
    assert.equal(kept[0].proteinAccession, '');
    assert.equal(rows[0].protein, 'p.Val14dup', 'the rows given are not changed');
    assert.deepEqual(h.withoutUncheckedProteinNames(null), []);

    // VEP's answer for the left-aligned spelling: a reverse-strand row it names
    // is its row; a forward-strand row in it is NOT taken (left-aligned is the
    // wrong spelling to ask about a forward-strand transcript by).
    const leftAligned = h.transcriptsFromVep([
        { transcript_id: 'NM_004985.5', hgvsc: 'NM_004985.5:c.36_38dup', hgvsp: 'NP_004976.2:p.Gly13dup', strand: -1, gene_symbol: 'KRAS',
            consequence_terms: ['inframe_insertion'] },
        { transcript_id: 'NM_000492.4', hgvsc: 'NM_000492.4:c.1521_1523del', strand: 1, gene_symbol: 'CFTR', consequence_terms: ['intron_variant'] },
        { transcript_id: 'NM_033360.4', hgvsc: 'NM_033360.4:c.36_38dup', hgvsp: 'NP_203524.1:p.Gly13dup', strand: -1, gene_symbol: 'KRAS' },
        { transcript_id: 'NR_000001.1', hgvsc: 'NR_000001.1:n.5_7dup', strand: -1, consequence_terms: ['splice_donor_variant'] }
    ]);
    const settled = h.reverseStrandRowsFrom(rows.concat(snpeff), leftAligned);
    assert.deepEqual(settled.rows.map(row => `${row.transcript} ${row.protein || '-'} ${row.source}`),
        ['NM_004985.5 p.Gly13dup vep', 'NM_000492.4 p.Phe508del vep', 'NR_000001.1 - vep', 'NM_004985.5 p.Gly13dup snpeff']);
    assert.equal(settled.replaced, true);
    assert.equal(settled.withheld, false);
    assert.equal(settled.rows[1], rows[1], 'the forward-strand row is the one it was');
    assert.equal(settled.rows[3], snpeff);
    // (...also when the left-aligned answer has a protein name for it.)
    const forwardNamed = h.transcriptsFromVep([{ transcript_id: 'NM_000492.4', hgvsc: 'NM_000492.4:c.1521_1523del',
        hgvsp: 'NP_000483.3:p.Phe508del', strand: 1, gene_symbol: 'CFTR', consequence_terms: ['inframe_deletion', 'splice_region_variant'] }]);
    const forward = h.reverseStrandRowsFrom([rows[1]], forwardNamed);
    assert.equal(forward.rows[0], rows[1]);
    assert.deepEqual([forward.replaced, forward.withheld], [false, false]);
    // (The whole row, its effect with it.)
    assert.equal(settled.rows[0], leftAligned[0]);
    // A row is taken for a protein NAME only. Where the left-aligned answer
    // has none for the transcript, the row in hand stays — its effect is VEP's
    // for the spelling in hand, not for one that may fall outside the exon...
    assert.equal(settled.rows[2], rows[2], 'a row the left-aligned answer has no protein name for is the one it was');
    const unnamed = h.transcriptsFromVep([
        { transcript_id: 'NM_004985.5', hgvsc: 'NM_004985.5:c.110_111+1dup', strand: -1, gene_symbol: 'KRAS',
            consequence_terms: ['inframe_insertion', 'splice_region_variant'] }
    ]);
    const unnamedLeft = h.transcriptsFromVep([
        { transcript_id: 'NM_004985.5', hgvsc: 'NM_004985.5:c.110_111+1dup', strand: -1, gene_symbol: 'KRAS', consequence_terms: ['splice_donor_variant'] }
    ]);
    const asItWas = h.reverseStrandRowsFrom(unnamed, unnamedLeft);
    assert.equal(asItWas.rows[0], unnamed[0]);
    assert.equal(asItWas.rows[0].effect, 'inframe_insertion, splice_region_variant');
    assert.deepEqual([asItWas.replaced, asItWas.withheld], [false, false], 'nothing taken, nothing taken away');
    // ...and a name in hand that the left-aligned answer does not bear out is
    // not shown: the rest of the row stays.
    const unconfirmed = h.reverseStrandRowsFrom(rows, unnamedLeft.map(row => ({ ...row, cDNA: 'c.36_38dup' })));
    assert.deepEqual(unconfirmed.rows.map(row => row.protein), ['', 'p.Phe508del', '']);
    assert.equal(unconfirmed.rows[0].effect, 'inframe_insertion', 'its effect is the one it had');
    assert.deepEqual([unconfirmed.replaced, unconfirmed.withheld], [false, true]);
    // A reverse-strand row the left-aligned answer does not have loses its name...
    const partial = h.reverseStrandRowsFrom(rows, leftAligned.slice(1));
    assert.deepEqual(partial.rows.map(row => row.protein), ['', 'p.Phe508del', '']);
    assert.deepEqual([partial.replaced, partial.withheld], [false, true]);
    // ...as they all do when there is no such answer.
    const none = h.reverseStrandRowsFrom(rows, null);
    assert.deepEqual(none.rows.map(row => row.protein), ['', 'p.Phe508del', '']);
    assert.deepEqual(none.rows.map(row => row.cDNA), ['c.36_38dup', 'c.1521_1523del', 'n.5_7dup']);
    assert.deepEqual([none.replaced, none.withheld], [false, true]);
    // Nothing is withheld where there was no name to withhold.
    assert.equal(h.reverseStrandRowsFrom([rows[1], rows[2]], null).withheld, false);
    // A transcript is matched by its accession as written, version and all.
    const otherVersion = h.reverseStrandRowsFrom(rows, [{ ...leftAligned[0], transcript: 'NM_004985.4' }]);
    assert.equal(otherVersion.rows[0].protein, '', 'another version of the transcript is another row');
    assert.deepEqual(h.reverseStrandRowsFrom(null, null), { rows: [], replaced: false, withheld: false });
}
console.log('Indel protein-name unit checks passed');

const noRecoder = [/variant_recoder/];
const noPlacement = [/variant_recoder/, /api\.ncbi\.nlm\.nih\.gov\/variation/];
const noReference = [/sequence\/region/];

// 1. KRAS c.34_36dup — p.Gly13dup — by every route to it. The spelling the
//    lookup shows depends on who answered; the protein name must not. It was
//    p.Val14dup by the third and fourth routes (VEP, asked by the name or by
//    the spelling it gives back for the name) and for the allele typed so.
{
    const routes = [
        ['KRAS c.34_36dup', {}, 'chr12:g.25398284_25398286dup'],
        ['KRAS c.34_36dup', { without: noRecoder }, 'chr12:g.25398286_25398287insCCA'],
        ['KRAS c.34_36dup', { without: noPlacement }, 'chr12:g.25398285_25398286insACC'],
        ['NM_004985.5:c.34_36dup', { without: noPlacement }, 'chr12:g.25398285_25398286insACC'],
        // Typed as that spelling, and as a row that is not left-aligned.
        ['chr12:g.25398285_25398286insACC', { without: noRecoder }, 'chr12:g.25398285_25398286insACC'],
        ['12 25398285 C CACC', { without: noRecoder }, 'chr12:g.25398285_25398286insACC']
    ];
    for (const [typed, options, genomic] of routes) {
        const label = `${typed}${options.without ? ` without ${options.without.length} service(s)` : ''}`;
        const page = await runLookup(cassettes.krasDup, { input: typed, ...options });
        assertResolved(assert, page, label);
        assert(page.variant.startsWith(`Variantg.: ${genomic}VCF (hg19): 12-25398280-G-GCCA (left-aligned)`), `${label}: ${page.variant.slice(0, 90)}`);
        assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.36_38dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/, label);
        assert.doesNotMatch(page.variant + page.panel, /Val14dup|Gly13_Val14insCys/, `${label}: not Ensembl's name for another spelling`);
        // VEP is asked about the allele as it is spelled, as it always was, and
        // left-aligned for the rows of its reverse-strand transcripts.
        assert.deepEqual(allelesAsked(page), [genomic.replace(/^chr/, ''), '12:g.25398280_25398281insCCA'], label);
        assert.match(page.panel, askedAsRemark('chr12:g.25398280_25398281insCCA'), label);
        assert.doesNotMatch(page.panel, withheldRemark, label);
        // Every transcript row, not only the one on the card.
        assert.match(page.variant, /NM_033360\.4: c\.36_38dup, NP_203524\.1:p\.Gly13dup — inframe_insertion/, label);
        // What the evidence sources are asked about is the same change.
        assert(page.requests.some(r => /\/api\/civic\?gene=KRAS&protein=Gly13dup/.test(r)), `${label}: CIViC is asked about Gly13dup`);
    }
    // The row as pasted when it IS left-aligned: one question, and nothing to say.
    const leftAligned = await runLookup(cassettes.krasDup, { input: '12 25398280 G GCCA', without: noRecoder });
    assertResolved(assert, leftAligned, 'left-aligned row');
    assert.match(leftAligned.variant, /g\.: chr12:g\.25398280_25398281insCCA.*p\.: p\.Gly13dup/);
    assert.deepEqual(allelesAsked(leftAligned), ['12:g.25398280_25398281insCCA']);
    assert.doesNotMatch(leftAligned.panel, /asked as/);
    console.log('ok indel protein name: KRAS c.34_36dup is p.Gly13dup by every route and spelling');
}

// 2. NRAS c.34_36dup — also p.Gly13dup. Here the spelling Ensembl's own
//    recoder gives (g.115258746_115258748dup) is one VEP misnames, so the
//    lookup showed "p.Val14delinsGlyPhe" with every service answering.
{
    const routes = [
        ['NRAS c.34_36dup', {}, 'chr1:g.115258746_115258748dup'],
        ['NRAS c.34_36dup', { without: noRecoder }, 'chr1:g.115258748_115258749insACC'],
        ['NRAS c.34_36dup', { without: noPlacement }, 'chr1:g.115258748_115258749insACC'],
        ['chr1:g.115258746_115258748dup', { without: noRecoder }, 'chr1:g.115258746_115258748dup']
    ];
    for (const [typed, options, genomic] of routes) {
        const label = `${typed}${options.without ? ` without ${options.without.length} service(s)` : ''}`;
        const page = await runLookup(cassettes.nrasDup, { input: typed, ...options });
        assertResolved(assert, page, label);
        assert(page.variant.startsWith(`Variantg.: ${genomic}VCF (hg19): 1-115258741-A-ACAC (left-aligned)`), `${label}: ${page.variant.slice(0, 90)}`);
        assert.match(page.variant, /Transcript: NM_002524\.5 \(MANE Select\)c\.: c\.38_40dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/, label);
        assert.doesNotMatch(page.variant + page.panel, /Val14delinsGlyPhe|Val14dup/, label);
        assert.deepEqual(allelesAsked(page), [genomic.replace(/^chr/, ''), '1:g.115258741_115258742insCAC'], label);
        assert.match(page.panel, askedAsRemark('chr1:g.115258741_115258742insCAC'), label);
    }
    console.log('ok indel protein name: NRAS c.34_36dup is p.Gly13dup, not what Ensembl calls its own spelling of it');
}

// 3. BRAF p.Thr599dup written at its 3'-most genomic position — a spelling
//    VEP calls an inserted stop ("p.Thr599_Val600insTer"). Shown whenever the
//    recoder had not re-spelled the input by the time VEP was asked — and,
//    recoder or not, when the same allele is written as a replacement
//    ("g.140453142delinsGTAG": the G kept, TAG put after it).
{
    for (const [typed, options, genomic] of [
        [undefined, {}, 'chr7:g.140453140_140453142dup'],
        [undefined, { without: noRecoder }, 'chr7:g.140453142_140453143insTAG'],
        ['chr7:g.140453142delinsGTAG', { without: noRecoder }, 'chr7:g.140453142delinsGTAG']
    ]) {
        const label = `${typed || cassettes.brafT599dupRight.input}${options.without ? ' without the recoder' : ''}`;
        const page = await runLookup(cassettes.brafT599dupRight, { ...(typed ? { input: typed } : {}), ...options });
        assertResolved(assert, page, label);
        assert(page.variant.startsWith(`Variantg.: ${genomic}VCF (hg19): 7-140453138-T-TGTA`), `${label}: ${page.variant.slice(0, 90)}`);
        assert.match(page.variant, /Transcript: NM_004333\.6 \(MANE Select\)c\.: c\.1794_1796dupp\.: p\.Thr599dup \(p\.T599dup\)Effect: inframe_insertion/, label);
        assert.doesNotMatch(page.variant + page.panel, /insTer|stop_gained/, `${label}: no stop is inserted`);
        assert.deepEqual(allelesAsked(page), [genomic.replace(/^chr/, ''), '7:g.140453138_140453139insGTA'], label);
        assert.match(page.panel, askedAsRemark('chr7:g.140453138_140453139insGTA'), label);
        if (!typed) assert.match(page.clinvar, /Variation ID: 162793/, `${label}: and ClinVar's record of the duplication`);
    }
    console.log('ok indel protein name: BRAF p.Thr599dup is not an inserted stop, at its 3\'-most spelling or as a replacement');
}

// 4. The reference unavailable: the left-aligned spelling cannot be worked
//    out, and what VEP calls the protein change on a reverse-strand transcript
//    for the spelling in hand is not shown as fact.
{
    const runs = [
        [cassettes.krasDup, 'KRAS c.34_36dup', noReference, '12:g.25398284_25398286dup', /c\.: c\.36_38dupp\.: N\/AEffect: inframe_insertion/],
        [cassettes.krasDup, 'KRAS c.34_36dup', noPlacement.concat(noReference), '12:g.25398285_25398286insACC', /c\.: c\.36_38dupp\.: N\/AEffect: inframe_insertion/],
        [cassettes.nrasDup, 'NRAS c.34_36dup', noReference, '1:g.115258746_115258748dup', /c\.: c\.38_40dupp\.: N\/AEffect: inframe_insertion/],
        [cassettes.brafT599dupRight, undefined, noRecoder.concat(noReference), '7:g.140453142_140453143insTAG', /c\.: c\.1794_1796dupp\.: N\/AEffect: inframe_insertion/]
    ];
    for (const [cassette, typed, without, asked, card] of runs) {
        const label = `${typed || cassette.input}, reference unavailable (${without.length} service(s) away)`;
        const page = await runLookup(cassette, { ...(typed ? { input: typed } : {}), without });
        assertResolved(assert, page, label);
        assert.match(page.variant, card, label);
        assert.doesNotMatch(page.variant + page.panel, /p\.(?:Gly13|Val14|Thr599)/, `${label}: no protein name, right or wrong`);
        assert.match(page.panel, withheldRemark, label);
        assert.match(page.panel, /▲Consequences · Ensembl VEP/, `${label}: flagged, not a routine line`);
        assert.deepEqual(allelesAsked(page), [asked], `${label}: asked by the spelling in hand`);
        assert.doesNotMatch(page.panel, /asked as/, label);
        // Nor is the name sent anywhere as a search term.
        assert(!page.requests.some(r => /Val14|insTer|GlyPhe/.test(r)), `${label}: no source is asked about Ensembl's name`);
    }
    // A forward-strand gene's names do not depend on the spelling: CFTR
    // F508del keeps its name without the reference.
    const cftr = await runLookup(cassettes.cftrF508del, { without: noReference });
    assertResolved(assert, cftr, 'CFTR F508del, reference unavailable');
    assert.match(cftr.variant, /c\.: c\.1521_1523delp\.: p\.Phe508del \(p\.F508del\)Effect: inframe_deletion/);
    assert.doesNotMatch(cftr.panel, withheldRemark);
    // A replacement is only looked into: with no reference to go by, BRAF's is
    // shown as VEP gives it, as it was before any of this.
    const delins = await runLookup(cassettes.brafT599dupRight, { input: 'chr7:g.140453142delinsGTAG', without: noRecoder.concat(noReference) });
    assertResolved(assert, delins, 'the replacement, reference unavailable');
    assert.match(delins.variant, /c\.: c\.1794_1796dupp\.: p\.Thr599_Val600insTer/);
    assert.doesNotMatch(delins.panel, withheldRemark);
    console.log('ok indel protein name: without the reference, a reverse-strand name is withheld and the panel says why');
}

// 5. VEP answers for the name but not for the allele (the last-resort route,
//    then a failure): the rows it gave for the name stay — they carry the
//    coding change — but their protein names were worked out where the name
//    puts the change, which for KRAS c.34_36dup is p.Val14dup.
{
    const page = await runLookup(cassettes.krasDup, { without: noPlacement.concat([/vep\/human\/hgvs\/12:g\./]) });
    assertResolved(assert, page, 'KRAS c.34_36dup, VEP silent about the allele');
    assert.match(page.variant, /g\.: chr12:g\.25398285_25398286insACC/);
    assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.36_38dupp\.: N\/AEffect: inframe_insertion/);
    assert.doesNotMatch(page.variant + page.panel, /Val14dup|Gly13dup/);
    assert.match(page.panel, withheldRemark);
    assert.doesNotMatch(page.panel, /asked as/, 'nothing came of asking');
    assert(vepAsked(page).includes('NM_004985.5:c.34_36dup'), 'VEP was asked by the name');
    assert.deepEqual(allelesAsked(page).sort(), ['12:g.25398280_25398281insCCA', '12:g.25398285_25398286insACC'], 'and, in vain, about the allele');
    // It answers for the allele as spelled but not left-aligned: the same.
    const half = await runLookup(cassettes.krasDup, { without: noPlacement.concat([/vep\/human\/hgvs\/12:g\.25398280_/]) });
    assertResolved(assert, half, 'KRAS c.34_36dup, VEP silent about the left-aligned spelling');
    assert.match(half.variant, /c\.: c\.36_38dupp\.: N\/AEffect: inframe_insertion/);
    assert.doesNotMatch(half.variant + half.panel, /Val14dup|Gly13dup/);
    assert.match(half.panel, withheldRemark);
    console.log('ok indel protein name: names VEP gave for another spelling are not shown when it does not answer for the left-aligned one');
}

// 6. Nothing changes where the spelling does not matter: a frameshift and a
//    substitution are put to VEP as they are spelled, once, and no reference
//    is read for VEP's sake.
{
    const brca1 = await runLookup(cassettes.brca1Dup, {});
    assertResolved(assert, brca1, 'BRCA1 c.5266dup');
    assert.match(brca1.variant, /g\.: chr17:g\.41209082dup.*p\.: p\.Gln1756ProfsTer74/);
    assert.deepEqual(allelesAsked(brca1), ['17:g.41209082dup'], 'a frameshift is asked about as it is spelled');
    assert.doesNotMatch(brca1.panel, /asked as/);

    const tp53 = await runLookup(cassettes.tp53Dup, {});
    assertResolved(assert, tp53, 'TP53 c.455dup row');
    assert.deepEqual(allelesAsked(tp53), ['17:g.7578474_7578475insG']);
    assert.match(tp53.variant, /p\.: p\.Pro153AlafsTer28/);

    const braf = await runLookup(cassettes.brafV600E, {});
    assertResolved(assert, braf, 'BRAF V600E');
    assert.deepEqual(allelesAsked(braf), ['7:g.140453136A>T']);
    assert.match(braf.variant, /p\.: p\.Val600Glu/);
    assert.equal(braf.requests.filter(r => r.includes('/sequence/region/')).length, 0, 'no reference read for a substitution');

    const frameshiftByName = await runLookup(cassettes.brca1InsC, {});
    assertResolved(assert, frameshiftByName, 'BRCA1 c.5265_5266insC');
    assert.match(frameshiftByName.variant, /c\.: c\.5266dupp\.: p\.Gln1756ProfsTer74/);
    assert.doesNotMatch(frameshiftByName.panel, /asked as|protein change not shown/);
    console.log('ok indel protein name: frameshifts and substitutions are asked about as before');
}

// 7. A forward-strand gene's rows are VEP's for the spelling in hand — NOT for
//    the left-aligned one. There the left-aligned spelling is the 5'-most, and
//    VEP works out the effect and the protein name from the bases it is given:
//    BRCA2 c.6842_6844dup sits on the first bases of exon 12, left-aligned it
//    falls in the intron before it, and VEP then gives no protein name at all.
//    Asking about the whole allele left-aligned, as this first did, showed
//    "p.: N/A · Effect: splice_region_variant, …, intron_variant" for it.
//    (Found in review. Both of VEP's answers below are as recorded.)
{
    const asSpelled = 'grch37.rest.ensembl.org/vep/human/hgvs/13:g.32918695_32918697dup?content-type=application/json&hgvs=1&refseq=1&mane=1&canonical=1';
    const leftAligned = 'grch37.rest.ensembl.org/vep/human/hgvs/13:g.32918692_32918693insAGG?content-type=application/json&hgvs=1&refseq=1&mane=1&canonical=1';
    const row = key => cassettes.brca2ExonEdgeDup.responses[key].body[0].transcript_consequences.find(t => t.transcript_id === 'NM_000059.4');
    assert.deepEqual([row(asSpelled).strand, row(asSpelled).hgvsp, row(asSpelled).consequence_terms.join(', ')],
        [1, 'NP_000050.3:p.Gly2281dup', 'inframe_insertion, splice_region_variant']);
    assert.deepEqual([row(leftAligned).hgvsp, row(leftAligned).consequence_terms.join(', ')],
        [undefined, 'splice_region_variant, splice_polypyrimidine_tract_variant, intron_variant'], 'what VEP says of the left-aligned spelling');
    for (const options of [{}, { without: noRecoder }]) {
        const label = `BRCA2 c.6842_6844dup${options.without ? ' without the recoder' : ''}`;
        const page = await runLookup(cassettes.brca2ExonEdgeDup, options);
        assertResolved(assert, page, label);
        assert.match(page.variant, /g\.: chr13:g\.32918695_32918697dupVCF \(hg19\): 13-32918692-T-TAGG \(left-aligned\)Gene: BRCA2/, label);
        assert.match(page.variant, /Transcript: NM_000059\.4 \(MANE Select\)c\.: c\.6842_6844dupp\.: p\.Gly2281dup \(p\.G2281dup\)Effect: inframe_insertion, splice_region_variant/, label);
        assert.doesNotMatch(page.variant + page.panel, /intron_variant|splice_polypyrimidine/, label);
        // Nothing is said of the left-aligned spelling, though it was asked
        // (which strand a transcript is on only VEP's answer says): no row
        // was taken from it.
        assert.deepEqual(allelesAsked(page), ['13:g.32918695_32918697dup', '13:g.32918692_32918693insAGG'], label);
        assert.doesNotMatch(page.panel, /asked as|protein change not shown/, label);
    }
    // Other forward-strand in-frame changes: the EGFR exon 20 duplication is
    // covered in tests/indel-identity.test.js; F508del by name:
    const cftr = await runLookup(cassettes.cftrF508delName, {});
    assertResolved(assert, cftr, 'CFTR F508del');
    assert.match(cftr.variant, /c\.: c\.1521_1523delp\.: p\.Phe508del \(p\.F508del\)Effect: inframe_deletion/);
    assert.deepEqual(allelesAsked(cftr), ['7:g.117199645_117199647del'], 'left-aligned as it is: one question');
    console.log('ok indel protein name: a forward-strand row is VEP\'s for the spelling in hand, not for the left-aligned one');
}

// 8. A reverse-strand row is taken from the left-aligned answer for a protein
//    NAME, and only for one. At the end of an exon the left-aligned spelling
//    can fall on the splice donor: KRAS c.109_111dup — one more Glu37, in front
//    of an untouched "gt" — is "c.110_111+1dup" to VEP at every spelling, with
//    no protein name at any, "inframe_insertion, splice_region_variant" as the
//    lookup spells it and "splice_donor_variant" left-aligned. Taking the whole
//    row from the left-aligned answer would have turned an in-frame insertion
//    into a splice-donor variant. (The mirror image of scenario 7; six of 31
//    exon-edge indels on reverse-strand transcripts were like this when asked.)
{
    const answer = spelling => cassettes.krasExonEndDup.responses[
        `grch37.rest.ensembl.org/vep/human/hgvs/${spelling}?content-type=application/json&hgvs=1&refseq=1&mane=1&canonical=1`]
        .body[0].transcript_consequences.find(t => t.transcript_id === 'NM_004985.5');
    assert.deepEqual([answer('12:g.25398208_25398210dup').strand, answer('12:g.25398208_25398210dup').hgvsp,
        answer('12:g.25398208_25398210dup').consequence_terms.join(', ')], [-1, undefined, 'inframe_insertion, splice_region_variant']);
    assert.deepEqual([answer('12:g.25398206_25398207insCCT').hgvsp, answer('12:g.25398206_25398207insCCT').consequence_terms.join(', ')],
        [undefined, 'splice_donor_variant'], 'what VEP says of the left-aligned spelling');
    for (const [typed, options] of [[undefined, {}], ['chr12:g.25398208_25398210dup', { without: noRecoder }]]) {
        const label = `${typed || 'KRAS c.109_111dup'}${options.without ? ' without the recoder' : ''}`;
        const page = await runLookup(cassettes.krasExonEndDup, { ...(typed ? { input: typed } : {}), ...options });
        assertResolved(assert, page, label);
        assert.match(page.variant, /g\.: chr12:g\.25398208_25398210dupVCF \(hg19\): 12-25398206-A-ACCT \(left-aligned\)Gene: KRAS/, label);
        assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.110_111\+1dupp\.: N\/AEffect: inframe_insertion, splice_region_variant/, label);
        assert.doesNotMatch(page.variant + page.panel, /splice_donor_variant/, label);
        // It was asked both ways; no row came of the second, and nothing was
        // taken away, so there is nothing to say of either.
        assert.deepEqual(allelesAsked(page), ['12:g.25398208_25398210dup', '12:g.25398206_25398207insCCT'], label);
        assert.doesNotMatch(page.panel, /asked as|protein change not shown/, label);
    }
    console.log('ok indel protein name: a reverse-strand row is replaced for a protein name only, not for a splice-site reading of the left-aligned spelling');
}

// 9. A lookup that has been overtaken leaves the next one's rows alone. "KRAS
//    c.34_36dup" (recoder and NCBI silent) is waiting for VEP about the allele
//    when "BRAF V600E" is submitted; VEP then fails. The first lookup, with
//    nothing left to show, went on to take the protein names off "its" rows —
//    which by then were BRAF's: "c.1799T>A · p.: N/A", a warning about
//    in-frame insertions in BRAF's panel, and CIViC asked about BRAF with no
//    protein change. (Found in review.)
{
    const first = cassettes.krasDup;
    const second = cassettes.brafV600E;
    const both = { input: 'KRAS c.34_36dup', reference: [...first.reference, ...second.reference],
        api: { ...(first.api || {}), ...(second.api || {}) }, responses: { ...first.responses, ...second.responses } };
    let submitted = false;
    let secondHasRows = false;
    let failures = 0;
    let lastFailure = 0;
    const firstHasFailed = () => failures >= 4 && Date.now() - lastFailure > 700;
    const page = await runLookup(both, {
        // (Waits at their real length, up to this: the second lookup has to be
        // waiting for its MANE label when the first one's request fails.)
        settleMs: 2500, timerCap: 4000,
        route: async (url, key, respond, options, { document }) => {
            // Neither Ensembl's recoder nor NCBI places the first name: VEP does, by the name.
            if (/variant_recoder\/human\/(?:NM_004985|KRAS)|api\.ncbi\.nlm\.nih\.gov\/variation/.test(key)) {
                return respond({ error: 'Service Unavailable' }, 503);
            }
            // Its question about the allele is out when the next variant is
            // submitted — and fails, late, every attempt of it.
            if (/\/vep\/human\/hgvs\/12:g\./.test(key)) {
                if (!submitted) {
                    submitted = true;
                    await wait(60);
                    document.getElementById('variantInput').value = 'BRAF V600E';
                    document.getElementById('variantForm').dispatchEvent(new document.defaultView.Event('submit', { cancelable: true }));
                    for (let i = 0; i < 600 && !secondHasRows; i++) await wait(5);
                    await wait(80);
                }
                failures += 1;
                lastFailure = Date.now();
                return respond({ error: 'Service Unavailable' }, 503);
            }
            if (/\/vep\/human\/hgvs\/7:g\.140453136A>T/.test(key)) {
                setTimeout(() => { secondHasRows = true; }, 0);
                return undefined;
            }
            // The second lookup has its rows, and waits for the MANE label
            // until the first one's request has failed for good.
            if (submitted && /lookup\/symbol\/homo_sapiens\/BRAF|datasets\/v2\/gene\/symbol\/BRAF\//.test(key)) {
                for (let i = 0; i < 700 && !firstHasFailed(); i++) await wait(5);
                await wait(60);
                return undefined;
            }
            return undefined;
        }
    });
    assert(submitted && secondHasRows && failures >= 4, 'the first lookup\'s VEP request failed while the second waited with its rows');
    assertResolved(assert, page, 'BRAF V600E, submitted while KRAS c.34_36dup waited for VEP');
    assert.match(page.panel, /You entered\s*BRAF V600E/);
    assert.match(page.variant, /g\.: chr7:g\.140453136A>T.*c\.: c\.1799T>Ap\.: p\.Val600Glu \(p\.V600E\)Effect: missense_variant/);
    assert.doesNotMatch(page.panel, /protein change not shown|Consequences · Ensembl VEP/, 'no VEP line: nothing went wrong in this lookup');
    assert.deepEqual(page.api('civic'), ['civic?gene=BRAF&protein=Val600Glu']);
    assert.doesNotMatch(page.panel + page.variant, /KRAS|25398285/, 'nothing of the lookup that was overtaken');
    console.log('ok indel protein name: an overtaken lookup does not take the next one\'s protein names');
}

// 10. An insertion that names no place between two neighbouring bases is not
//     made into one. "chr12:g.25398280_25398285insCCA" (five bases between its
//     flanks: a "delins" that lost its "del", perhaps) and
//     "chr12:g.25398280insCCA" are refused by Ensembl's recoder and by VEP
//     alike, and the lookup fails. Left-aligning read both as "insert after
//     25398280" and showed KRAS p.Gly13dup, "asked as … — the same allele".
//     (Found in review. Ensembl's refusal is HTTP 400; a spelling this
//     recording does not hold is answered as an outage, which fails the lookup
//     the same way.)
{
    for (const malformed of ['chr12:g.25398280_25398285insCCA', 'chr12:g.25398280insCCA']) {
        const page = await runLookup(cassettes.krasDup, { input: malformed });
        assert.match(page.status, /^Error: /, `${malformed}: not resolved`);
        assert.equal(page.hidden, true, `${malformed}: nothing is rendered`);
        assert.deepEqual(allelesAsked(page), [malformed.replace(/^chr/, '')], `${malformed}: VEP is asked about what was typed, and nothing else`);
        assert.doesNotMatch(page.panel + page.variant, /Gly13dup|asked as|the same allele/, malformed);
        assert(!page.requests.some(r => r.includes('25398280_25398281insCCA')), `${malformed}: no allele is made of it`);
    }
    console.log('ok indel protein name: an insertion that names no place is refused, not re-spelled');
}

// 11. VEP's rows for the left-aligned spelling are used whatever became of its
//     rows for the spelling in hand. A row whose gene column is not the symbol
//     VEP uses ("KRAS2", an old name for KRAS) matches none of them, so the
//     rows already there stay: they got the right names from the left-aligned
//     answer all the same. It used to read "p.: N/A", under "asked as … — the
//     same allele, left-aligned · protein change not shown: … it could not be
//     asked that way". (Found in review.)
{
    const page = await runLookup(cassettes.krasDup, { input: '12 25398285 KRAS2 C CACC', without: noRecoder });
    assertResolved(assert, page, 'KRAS duplication, the gene column an alias');
    assert.match(page.variant, /c\.: c\.36_38dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/);
    assert.match(page.panel, askedAsRemark('chr12:g.25398280_25398281insCCA'));
    assert.doesNotMatch(page.panel, withheldRemark);
    console.log('ok indel protein name: the left-aligned answer names the rows that are there, whichever they are');
}

// 12. The reference read that never completes. Its request is timed; an
//     answer that stalls once it has begun is not, and the lookup waited on it
//     for ever. The wait has a deadline now, after which the reverse-strand
//     names are withheld as they are when the reference cannot be read.
//     (Found in review.)
{
    const stalled = runLookup(cassettes.nrasDup, {
        route: (url, key) => (/\/sequence\/region\/human\/1:115258596\.\.115258898/.test(key) ? new Promise(() => {}) : undefined) });
    const page = await Promise.race([stalled, wait(8000).then(() => null)]);
    assert(page, 'the lookup finishes');
    assertResolved(assert, page, 'NRAS c.34_36dup, the reference read stalled');
    assert.match(page.variant, /g\.: chr1:g\.115258746_115258748dup/);
    assert.match(page.variant, /c\.: c\.38_40dupp\.: N\/AEffect: inframe_insertion/);
    assert.match(page.panel, withheldRemark);
    assert.deepEqual(allelesAsked(page), ['1:g.115258746_115258748dup']);
    console.log('ok indel protein name: a reference read that never completes does not hold the result');
}

// 13. The reference is read once for both of its uses (the left-aligned
//     spelling, and the allele's other MyVariant ids).
{
    const page = await runLookup(cassettes.nrasDup, {});
    const windows = page.requests.filter(r => r.includes('/sequence/region/human/1:115258596..115258898'));
    assert.equal(windows.length, 1, `one read of the window: ${page.requests.filter(r => r.includes('/sequence/')).join(' | ')}`);
    console.log('ok indel protein name: one read of the reference serves both steps');
}

// 14. The second question is asked alongside the first, not after it: which
//     strand a transcript is on only VEP's answer says, and waiting for that
//     before asking again would cost every such lookup a round trip. VEP's
//     answer for the allele as spelled is held here until the left-aligned
//     question has been put (or for half a second).
{
    let asked = false;
    let askedWhileWaiting = false;
    const page = await runLookup(cassettes.nrasDup, {
        route: async (url, key) => {
            if (/\/vep\/human\/hgvs\/1:g\.115258741_115258742insCAC/.test(key)) asked = true;
            if (/\/vep\/human\/hgvs\/1:g\.115258746_115258748dup/.test(key)) {
                for (let i = 0; i < 100 && !asked; i++) await wait(5);
                askedWhileWaiting = askedWhileWaiting || asked;
            }
            return undefined;
        }
    });
    assertResolved(assert, page, 'NRAS c.34_36dup, the first answer held');
    assert.equal(askedWhileWaiting, true, 'the left-aligned spelling was asked about while the first answer was awaited');
    assert.match(page.variant, /c\.: c\.38_40dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/);
    console.log('ok indel protein name: the left-aligned question is asked alongside the first');
}

// 15. VEP's answer for the left-aligned spelling only confirms names, and has
//     a deadline of its own: one that never comes does not hold the result,
//     and the names it would have confirmed are not shown.
{
    const pending = runLookup(cassettes.nrasDup, { input: 'chr1:g.115258746_115258748dup', without: noRecoder,
        route: (url, key) => (/\/vep\/human\/hgvs\/1:g\.115258741_115258742insCAC/.test(key) ? new Promise(() => {}) : undefined) });
    const page = await Promise.race([pending, wait(8000).then(() => null)]);
    assert(page, 'the lookup finishes');
    assertResolved(assert, page, 'NRAS duplication, the left-aligned answer never arriving');
    assert.match(page.variant, /g\.: chr1:g\.115258746_115258748dup/);
    assert.match(page.variant, /c\.: c\.38_40dupp\.: N\/AEffect: inframe_insertion/);
    assert.doesNotMatch(page.variant + page.panel, /Val14delinsGlyPhe|Gly13dup/);
    assert.match(page.panel, withheldRemark);
    console.log('ok indel protein name: a left-aligned answer that never comes does not hold the result');
}
