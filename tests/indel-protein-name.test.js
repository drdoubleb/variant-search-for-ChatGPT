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
 *     came out as p.Val14dup whenever the recoder had not answered;
 *   - the left-aligned spelling was named right every time.
 *
 * So VEP is asked about an in-frame indel left-aligned, whichever spelling the
 * lookup arrived at. Where that spelling cannot be worked out (the reference
 * around the indel is unavailable), the protein names VEP gives on reverse-
 * strand transcripts are not shown, and the panel says why.
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
// (VEP is asked twice per spelling when its RefSeq set has no row; once here.)
const allelesAsked = page => Array.from(new Set(vepAsked(page).filter(asked => /:g\./.test(asked))));
const withheldRemark = /protein change not shown: Ensembl names an in-frame insertion or deletion on the reverse strand reliably only when it is written left-aligned, and it could not be asked that way/;
const askedAsRemark = spelling => new RegExp(`asked as ${spelling.replace(/[.]/g, '\\.')} — the same allele, left-aligned: `
    + 'Ensembl\'s protein name for an in-frame insertion or deletion depends on where it is written');

// ── Which alleles this is about ────────────────────────────────────────────
for (const inFrame of ['chr12:g.25398284_25398286dup', 'chr12:g.25398280_25398281insCCA', 'chr7:g.117199646_117199648del',
    'chrX:g.20156694_20156723dup', 'chr7:g.117199646_117199648delCTT', 'chr12:g.25398280_25398281insCCACCA']) {
    assert.equal(h.inFrameIndel(inFrame), true, inFrame);
}
for (const other of ['chr17:g.41209082dup', 'chr13:g.32906698_32906699insCA', 'chr17:g.41276045_41276046del', 'chr12:g.25398284C>T',
    'chr12:g.25398284_25398286delinsTTT', 'chr12:g.25398284_25398286', 'chr12:g.25398280_25398281ins', 'NM_004985.5:c.34_36dup',
    'chr12:g.25398284_25398287dup', '', null, undefined]) {
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
// Every indel, not only in-frame ones (which of them VEP is asked this way is
// inFrameIndel's to say).
assert.equal(h.leftAlignedSpelling('chr12:g.25398282dup', kras.start, kras.seq), 'chr12:g.25398280_25398281insC');
assert.equal(h.leftAlignedSpelling('chr12:g.25398282del', kras.start, kras.seq), 'chr12:g.25398281del');
// What the window cannot say is not guessed.
assert.equal(h.leftAlignedSpelling('chr12:g.25398284C>T', kras.start, kras.seq), null, 'not an indel');
assert.equal(h.leftAlignedSpelling('chr12:g.25398284_25398286delAAA', kras.start, kras.seq), null, 'bases the reference does not hold');
assert.equal(h.leftAlignedSpelling('chr12:g.25398284_25398286dup', 25398283, referenceSlice([kras], '12', 25398283, 25398300)), null,
    'slides off the window: the left end is not known');
assert.equal(h.leftAlignedSpelling('chr12:g.99998284_99998286dup', kras.start, kras.seq), null, 'outside the window');
assert.equal(h.leftAlignedSpelling('chr12:g.25398284_25398286dup', kras.start, ''), null);

// ── What VEP is asked by ───────────────────────────────────────────────────
{
    const calls = [];
    const reference = async (chrom, start, end) => {
        calls.push([chrom, start, end]);
        return referenceSlice([kras], String(chrom).replace(/^chr/, ''), start, end);
    };
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286dup', reference),
        { askedBy: 'chr12:g.25398280_25398281insCCA', checked: true });
    // (The window the lookup reads for the allele's other records: one request.)
    assert.deepEqual(calls, [['chr12', 25398134, 25398436]]);
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398280_25398281insCCA', reference),
        { askedBy: 'chr12:g.25398280_25398281insCCA', checked: true }, 'left-aligned already');
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398283_25398285del', reference),
        { askedBy: 'chr12:g.25398281_25398283del', checked: true });

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
    } finally {
        console.warn = warn;
    }
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286dup', async () => null),
        { askedBy: 'chr12:g.25398284_25398286dup', checked: false });
    // ...as it is for a spelling that does not fit the reference read,
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25398286delAAA', reference),
        { askedBy: 'chr12:g.25398284_25398286delAAA', checked: false });
    // ...and for an event wider than the window is read for.
    calls.length = 0;
    assert.deepEqual(await h.vepSpellingFor('chr12:g.25398284_25399786del', reference),
        { askedBy: 'chr12:g.25398284_25399786del', checked: false });
    assert.deepEqual(calls, []);
}

// ── Names that are not shown ───────────────────────────────────────────────
{
    const rows = h.transcriptsFromVep([
        { transcript_id: 'NM_004985.5', hgvsc: 'NM_004985.5:c.36_38dup', hgvsp: 'NP_004976.2:p.Val14dup', strand: -1, gene_symbol: 'KRAS' },
        { transcript_id: 'NM_000492.4', hgvsc: 'NM_000492.4:c.1521_1523del', hgvsp: 'NP_000483.3:p.Phe508del', strand: 1, gene_symbol: 'CFTR' },
        { transcript_id: 'NR_000001.1', hgvsc: 'NR_000001.1:n.5_7dup', strand: -1 }
    ]);
    assert.deepEqual(rows.map(row => row.strand), [-1, 1, -1], 'a row says which strand its transcript is on');
    // (snpEff's name for the allele is its own work, whatever strand it is on.)
    const kept = h.withoutUncheckedProteinNames(rows.concat(
        { transcript: 'NM_004985.5', cDNA: 'c.36_38dup', protein: 'p.Gly13dup', source: 'snpeff', strand: -1 }));
    assert.deepEqual(kept.map(row => row.protein), ['', 'p.Phe508del', '', 'p.Gly13dup'],
        'reverse-strand names VEP gave go; forward-strand ones and another source\'s stay');
    assert.deepEqual(kept.map(row => row.cDNA), ['c.36_38dup', 'c.1521_1523del', 'n.5_7dup', 'c.36_38dup'], 'the coding name stays');
    assert.deepEqual(kept.map(row => Boolean(row.proteinWithheld)), [true, false, false, false], 'only a name taken away is flagged');
    assert.equal(kept[0].proteinAccession, '');
    assert.equal(rows[0].protein, 'p.Val14dup', 'the rows given are not changed');
    assert.deepEqual(h.withoutUncheckedProteinNames(null), []);
}
console.log('Indel protein-name unit checks passed');

const noRecoder = [/variant_recoder/];
const noPlacement = [/variant_recoder/, /api\.ncbi\.nlm\.nih\.gov\/variation/];
const noReference = [/sequence\/region/];

// 1. KRAS c.34_36dup — p.Gly13dup — by every route to it. The spelling the
//    lookup shows depends on who answered; the protein name must not. It was
//    p.Val14dup by the last two routes (VEP, asked by the name or by the
//    spelling it gives back for the name).
{
    const routes = [
        ['KRAS c.34_36dup', {}, /g\.: chr12:g\.25398284_25398286dup/],
        ['KRAS c.34_36dup', { without: noRecoder }, /g\.: chr12:g\.25398286_25398287insCCA/],
        ['KRAS c.34_36dup', { without: noPlacement }, /g\.: chr12:g\.25398285_25398286insACC/],
        ['NM_004985.5:c.34_36dup', { without: noPlacement }, /g\.: chr12:g\.25398285_25398286insACC/],
        // Typed as that spelling, and as a row that is not left-aligned.
        ['chr12:g.25398285_25398286insACC', { without: noRecoder }, /g\.: chr12:g\.25398285_25398286insACC/],
        ['12 25398285 C CACC', { without: noRecoder }, /g\.: chr12:g\.25398285_25398286insACC/],
        ['chr12:g.25398284_25398285insCAC', { without: noRecoder }, /g\.: chr12:g\.25398284_25398285insCAC/]
    ];
    for (const [typed, options, genomic] of routes) {
        const label = `${typed}${options.without ? ` without ${options.without.length} service(s)` : ''}`;
        const page = await runLookup(cassettes.krasDup, { input: typed, ...options });
        assertResolved(assert, page, label);
        assert.match(page.variant, genomic, label);
        assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.36_38dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/, label);
        assert.doesNotMatch(page.variant + page.panel, /Val14dup|Gly13_Val14insCys/, `${label}: not Ensembl's name for another spelling`);
        // About the allele, VEP is asked once — left-aligned.
        assert.deepEqual(allelesAsked(page), ['12:g.25398280_25398281insCCA'], label);
        assert.match(page.panel, askedAsRemark('chr12:g.25398280_25398281insCCA'), label);
        assert.doesNotMatch(page.panel, withheldRemark, label);
        // Every transcript row, not only the one on the card.
        assert.match(page.variant, /NM_033360\.4: c\.36_38dup, NP_203524\.1:p\.Gly13dup — inframe_insertion/, label);
        // What the evidence sources are asked about is the same change.
        assert(page.requests.some(r => /\/api\/civic\?gene=KRAS&protein=Gly13dup/.test(r)), `${label}: CIViC is asked about Gly13dup`);
    }
    // The row as pasted when it IS left-aligned: nothing to re-spell, nothing said.
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
        ['NRAS c.34_36dup', {}, /g\.: chr1:g\.115258746_115258748dup/],
        ['NRAS c.34_36dup', { without: noRecoder }, /g\.: chr1:g\.115258748_115258749insACC/],
        ['NRAS c.34_36dup', { without: noPlacement }, /g\.: chr1:g\.115258748_115258749insACC/],
        ['chr1:g.115258746_115258748dup', { without: noRecoder }, /g\.: chr1:g\.115258746_115258748dup/]
    ];
    for (const [typed, options, genomic] of routes) {
        const label = `${typed}${options.without ? ` without ${options.without.length} service(s)` : ''}`;
        const page = await runLookup(cassettes.nrasDup, { input: typed, ...options });
        assertResolved(assert, page, label);
        assert.match(page.variant, genomic, label);
        assert.match(page.variant, /VCF \(hg19\): 1-115258741-A-ACAC \(left-aligned\)/, label);
        assert.match(page.variant, /Transcript: NM_002524\.5 \(MANE Select\)c\.: c\.38_40dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/, label);
        assert.doesNotMatch(page.variant + page.panel, /Val14delinsGlyPhe|Val14dup/, label);
        assert.deepEqual(allelesAsked(page), ['1:g.115258741_115258742insCAC'], label);
        assert.match(page.panel, askedAsRemark('chr1:g.115258741_115258742insCAC'), label);
    }
    console.log('ok indel protein name: NRAS c.34_36dup is p.Gly13dup, not what Ensembl calls its own spelling of it');
}

// 3. BRAF p.Thr599dup written at its 3'-most genomic position — a spelling
//    VEP calls an inserted stop ("p.Thr599_Val600insTer"). Shown whenever the
//    recoder had not re-spelled the input by the time VEP was asked.
{
    for (const options of [{}, { without: noRecoder }]) {
        const label = `chr7:g.140453142_140453143insTAG${options.without ? ' without the recoder' : ''}`;
        const page = await runLookup(cassettes.brafT599dupRight, options);
        assertResolved(assert, page, label);
        assert.match(page.variant, options.without ? /g\.: chr7:g\.140453142_140453143insTAG/ : /g\.: chr7:g\.140453140_140453142dup/, label);
        assert.match(page.variant, /Transcript: NM_004333\.6 \(MANE Select\)c\.: c\.1794_1796dupp\.: p\.Thr599dup \(p\.T599dup\)Effect: inframe_insertion/, label);
        assert.doesNotMatch(page.variant + page.panel, /insTer|stop_gained/, `${label}: no stop is inserted`);
        assert.deepEqual(allelesAsked(page), ['7:g.140453138_140453139insGTA'], label);
        assert.match(page.clinvar, /Variation ID: 162793/, `${label}: and ClinVar's record of the duplication`);
    }
    console.log('ok indel protein name: BRAF p.Thr599dup at its 3\'-most spelling is not an inserted stop');
}

// 4. The reference unavailable: the left-aligned spelling cannot be worked
//    out, VEP is asked by the spelling in hand, and what it calls the protein
//    change on a reverse-strand transcript is not shown as fact.
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
    assert(vepAsked(page).includes('NM_004985.5:c.34_36dup'), 'VEP was asked by the name');
    assert.deepEqual(allelesAsked(page), ['12:g.25398280_25398281insCCA'], 'and, in vain, by the left-aligned allele');

    // The same failure for a change whose names do not depend on the spelling
    // (a frameshift) leaves the rows as VEP gave them.
    const frameshift = await runLookup(cassettes.brca1InsC, {});
    assertResolved(assert, frameshift, 'BRCA1 c.5265_5266insC');
    assert.match(frameshift.variant, /c\.: c\.5266dupp\.: p\.Gln1756ProfsTer74/);
    assert.doesNotMatch(frameshift.panel, /asked as|protein change not shown/);
    console.log('ok indel protein name: names VEP gave for the NAME are not shown when it does not answer for the allele');
}

// 6. Nothing changes where the spelling does not matter: a frameshift and a
//    substitution are put to VEP as they are spelled, at once (no reference is
//    read first), and a forward-strand in-frame deletion keeps its name.
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

    const cftr = await runLookup(cassettes.cftrF508delName, {});
    assertResolved(assert, cftr, 'CFTR F508del');
    assert.match(cftr.variant, /c\.: c\.1521_1523delp\.: p\.Phe508del \(p\.F508del\)Effect: inframe_deletion/);
    assert.deepEqual(allelesAsked(cftr), ['7:g.117199645_117199647del'], 'left-aligned: once');
    console.log('ok indel protein name: frameshifts, substitutions and forward-strand names are as before');
}

// 7. The reference is read once for both of its uses (the left-aligned
//    spelling, and the allele's other MyVariant ids).
{
    const page = await runLookup(cassettes.nrasDup, {});
    const windows = page.requests.filter(r => r.includes('/sequence/region/human/1:115258596..115258898'));
    assert.equal(windows.length, 1, `one read of the window: ${page.requests.filter(r => r.includes('/sequence/')).join(' | ')}`);
    console.log('ok indel protein name: one read of the reference serves both steps');
}
