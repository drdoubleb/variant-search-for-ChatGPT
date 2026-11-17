/*
 * The alleles of a protein substitution, without Ensembl.
 *
 * "KRAS Y64*" names a protein change. Which bases changed, only Ensembl's
 * variant recoder said — and it is the slowest thing the page asks: 3 to 17 s
 * per name on the day this was measured, one name in three unanswered after
 * 25 s. The page then waited out four attempts, about 50 s, and showed the
 * codon alone.
 *
 * MyVariant carries dbNSFP, which lists every non-synonymous single-base
 * change in the genome with the residue it changes. Asked for the residue it
 * answers in half a second (proteinSubstitutionAlleles), and what it lists are
 * records. For a substitution named on a gene symbol that the name search did
 * not find, it is now asked first:
 *
 *   - it lists alleles: those are the candidates, and Ensembl is not waited for;
 *   - no single-base change gives the protein change: Ensembl may still name
 *     the larger ones, within the time of one attempt — and where "X" was
 *     typed, that already settles how X is read;
 *   - it does not have the residue under that number, or does not answer:
 *     Ensembl, as before.
 *
 * The scenarios replay answers recorded from the live services.
 *
 * Run with: node tests/protein-substitution-alleles.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));

// ── The pieces ─────────────────────────────────────────────────────────────
const hit = (id, genename, aa, transcriptid) => ({ _id: id, dbnsfp: { genename, aa, ...(transcriptid ? { ensembl: { transcriptid } } : {}) } });
// (What is listed and whether the residue is known; the codon has its own checks below.)
const listed = (hits, gene, change) => {
    const { alleles, residueKnown } = h.proteinSubstitutionAlleles(hits, gene, change);
    return { alleles, residueKnown };
};
{
    // KRAS Tyr64 (TAC): the two changes of the third base make a stop.
    const tyr64 = [
        hit('chr12:g.25380266G>C', ['KRAS', 'KRAS'], { alt: 'X', pos: [64, 64], ref: 'Y' }),
        hit('chr12:g.25380266G>T', ['KRAS', 'KRAS'], { alt: 'X', pos: [64, 64], ref: 'Y' }),
        hit('chr12:g.25380267T>C', ['KRAS', 'KRAS'], { alt: 'C', pos: [64, 64], ref: 'Y' }),
        hit('chr12:g.25380268A>G', ['KRAS', 'KRAS'], { alt: 'H', pos: [64, 64], ref: 'Y' })
    ];
    assert.deepEqual(listed(tyr64, 'KRAS', { ref: 'Y', pos: 64, alt: '*' }),
        { alleles: ['chr12:g.25380266G>C', 'chr12:g.25380266G>T'], residueKnown: true });
    assert.deepEqual(h.proteinSubstitutionAlleles([...tyr64].reverse(), 'kras', { ref: 'Y', pos: 64, alt: '*' }).alleles,
        ['chr12:g.25380266G>C', 'chr12:g.25380266G>T'], 'in one order, whatever MyVariant\'s');
    assert.deepEqual(listed(tyr64, 'KRAS', { ref: 'Y', pos: 64, alt: 'H' }),
        { alleles: ['chr12:g.25380268A>G'], residueKnown: true });
    // No single base makes Tyr64 a tryptophan: the residue is known, the list empty.
    assert.deepEqual(listed(tyr64, 'KRAS', { ref: 'Y', pos: 64, alt: 'W' }), { alleles: [], residueKnown: true });
    // Another residue, another reference amino acid, another gene: not known — nothing is concluded.
    for (const [gene, change] of [['KRAS', { ref: 'Y', pos: 65, alt: '*' }], ['KRAS', { ref: 'F', pos: 64, alt: '*' }],
        ['NRAS', { ref: 'Y', pos: 64, alt: '*' }]]) {
        assert.deepEqual(listed(tyr64, gene, change), { alleles: [], residueKnown: false }, JSON.stringify([gene, change]));
    }
    assert.deepEqual(listed([], 'KRAS', { ref: 'Y', pos: 64, alt: '*' }), { alleles: [], residueKnown: false });
    assert.deepEqual(listed(null, 'KRAS', { ref: 'Y', pos: 64, alt: '*' }), { alleles: [], residueKnown: false });

    // A variant in two reading frames comes as two amino-acid changes (as
    // MyVariant gives MECP2's records). R306C is chrX:g.153296363G>A, which is
    // A301V in the other frame; the residue has to be right within ONE frame.
    const mecp2 = [
        hit('chrX:g.153296363G>A', 'MECP2', [{ alt: 'V', pos: 301, ref: 'A' }, { alt: 'C', pos: [306, 318], ref: 'R' }]),
        hit('chrX:g.153296362C>T', ['MECP2', 'MECP2'], { alt: 'H', pos: [306, 318], ref: 'R' }),
        hit('chrX:g.153296399G>A', 'MECP2', [{ alt: 'L', pos: 289, ref: 'P' }, { alt: 'X', pos: [294, 306], ref: 'R' }]),
        hit('chrX:g.153296398C>T', ['MECP2', 'MECP2'], { alt: 'Q', pos: [294, 306], ref: 'R' })
    ];
    assert.deepEqual(h.proteinSubstitutionAlleles(mecp2, 'MECP2', { ref: 'R', pos: 306, alt: 'C' }).alleles, ['chrX:g.153296363G>A']);
    assert.deepEqual(h.proteinSubstitutionAlleles(mecp2, 'MECP2', { ref: 'A', pos: 301, alt: 'V' }).alleles, ['chrX:g.153296363G>A']);
    // (Residue 306 is an arginine on both transcripts — two different codons.
    // One base makes the one a stop, and none the other.)
    assert.deepEqual(h.proteinSubstitutionAlleles(mecp2, 'MECP2', { ref: 'R', pos: 306, alt: '*' }).alleles, ['chrX:g.153296399G>A']);
    assert.deepEqual(h.proteinSubstitutionAlleles(mecp2, 'MECP2', { ref: 'R', pos: 294, alt: '*' }).alleles, ['chrX:g.153296399G>A']);
    // "R301C" and "A306V" are in no frame, though MyVariant's search, which
    // matches across frames, would return the record for either.
    assert.deepEqual(listed(mecp2, 'MECP2', { ref: 'R', pos: 301, alt: 'C' }), { alleles: [], residueKnown: false });
    assert.deepEqual(listed(mecp2, 'MECP2', { ref: 'A', pos: 306, alt: 'V' }), { alleles: [], residueKnown: false });
    assert.deepEqual(listed(mecp2, 'MECP2', { ref: 'A', pos: 301, alt: 'C' }), { alleles: [], residueKnown: true },
        'the residue is there; that change of it is not');

    // A stop that is lost: dbNSFP writes the stop as X on either side.
    const stopLost = [hit('chr17:g.7572928C>A', ['TP53'], { alt: 'L', pos: [355, 394], ref: 'X' })];
    assert.deepEqual(h.proteinSubstitutionAlleles(stopLost, 'TP53', { ref: '*', pos: 394, alt: 'L' }).alleles, ['chr17:g.7572928C>A']);
    // A symbol written with small letters is the same gene.
    assert.deepEqual(h.proteinSubstitutionAlleles([hit('chr9:g.27567117A>G', ['C9orf72'], { alt: 'T', pos: [1], ref: 'M' })],
        'C9ORF72', { ref: 'M', pos: 1, alt: 'T' }).alleles, ['chr9:g.27567117A>G']);
    // Only single-base changes are taken from it, each once.
    assert.deepEqual(h.proteinSubstitutionAlleles([
        hit('chr1:g.100A>G', 'GENE', { alt: 'D', pos: 5, ref: 'G' }), hit('chr1:g.100A>G', 'GENE', { alt: 'D', pos: 5, ref: 'G' }),
        hit('chr1:g.100_101delinsGG', 'GENE', { alt: 'D', pos: 5, ref: 'G' }), hit(undefined, 'GENE', { alt: 'D', pos: 5, ref: 'G' })
    ], 'GENE', { ref: 'G', pos: 5, alt: 'D' }).alleles, ['chr1:g.100A>G']);

    // Where the residue's codon lies, from which base of it each record changes
    // (as MyVariant gives them). KRAS Gly12 is read backwards: base 1 at
    // 25398285, base 2 at 25398284 — and no record at base 3, where every
    // change of GGT is still a glycine.
    const both = ['ENST00000311936', 'ENST00000256078'];
    const gly12 = [
        hit('chr12:g.25398285C>A', ['KRAS', 'KRAS'], { alt: 'C', codonpos: [1, 1], pos: [12, 12], ref: 'G' }, both),
        hit('chr12:g.25398284C>T', ['KRAS', 'KRAS'], { alt: 'D', codonpos: [2, 2], pos: [12, 12], ref: 'G' }, both)
    ];
    assert.deepEqual(h.proteinSubstitutionAlleles(gly12, 'KRAS', { ref: 'G', pos: 12, alt: '*' }),
        { alleles: [], residueKnown: true, codon: { chrom: 'chr12', start: 25398283, end: 25398285, transcript: 'ENST00000256078' } });
    // EGFR Thr790, read forwards — and numbered 745 on two other transcripts.
    const three = ['ENST00000455089', 'ENST00000454757', 'ENST00000275493'];
    const thr790 = [
        hit('chr7:g.55249070A>G', ['EGFR', 'EGFR', 'EGFR'], { alt: 'A', codonpos: [1, 1, 1], pos: [745, 745, 790], ref: 'T' }, three),
        hit('chr7:g.55249071C>T', ['EGFR', 'EGFR', 'EGFR'], { alt: 'M', codonpos: [2, 2, 2], pos: [745, 745, 790], ref: 'T' }, three)
    ];
    assert.deepEqual(h.proteinSubstitutionAlleles(thr790, 'EGFR', { ref: 'T', pos: 790, alt: 'M' }).codon,
        { chrom: 'chr7', start: 55249070, end: 55249072, transcript: 'ENST00000275493' });
    assert.deepEqual(h.proteinSubstitutionAlleles(thr790, 'EGFR', { ref: 'T', pos: 745, alt: 'M' }).codon,
        { chrom: 'chr7', start: 55249070, end: 55249072, transcript: 'ENST00000454757' }, 'the transcripts that number it so');
    // One base of the codon alone fits either strand: not enough.
    assert.equal(h.proteinSubstitutionAlleles(thr790.slice(0, 1), 'EGFR', { ref: 'T', pos: 790, alt: 'A' }).codon, null);
    // One number, two codons — MECP2's residue 306 on its two transcripts: not guessed at.
    const twoCodons = [
        hit('chrX:g.153296363G>A', 'MECP2', [{ alt: 'V', codonpos: 2, pos: 301, ref: 'A' }, { alt: 'C', codonpos: [1, 1], pos: [306, 318], ref: 'R' }]),
        hit('chrX:g.153296362C>T', ['MECP2', 'MECP2'], { alt: 'H', codonpos: [2, 2], pos: [306, 318], ref: 'R' }),
        hit('chrX:g.153296399G>A', 'MECP2', [{ alt: 'L', codonpos: 2, pos: 289, ref: 'P' }, { alt: 'X', codonpos: [1, 1], pos: [294, 306], ref: 'R' }]),
        hit('chrX:g.153296398C>T', ['MECP2', 'MECP2'], { alt: 'Q', codonpos: [2, 2], pos: [294, 306], ref: 'R' })
    ];
    assert.equal(h.proteinSubstitutionAlleles(twoCodons, 'MECP2', { ref: 'R', pos: 306, alt: 'C' }).codon, null);
    // ...while 318 and 294 are one codon each (the other frame's base is not counted in).
    assert.deepEqual(h.proteinSubstitutionAlleles(twoCodons, 'MECP2', { ref: 'R', pos: 318, alt: 'C' }).codon,
        { chrom: 'chrX', start: 153296361, end: 153296363, transcript: '' });
    assert.deepEqual(h.proteinSubstitutionAlleles(twoCodons, 'MECP2', { ref: 'R', pos: 294, alt: '*' }).codon,
        { chrom: 'chrX', start: 153296397, end: 153296399, transcript: '' });
    // Two reading frames that happen to give a residue the same number: only
    // the frame with the residue asked for says where the codon is.
    const sameNumber = [
        hit('chr1:g.1001C>T', 'GENE', [{ alt: 'V', codonpos: 2, pos: 306, ref: 'A' }, { alt: 'C', codonpos: 1, pos: 306, ref: 'R' }]),
        hit('chr1:g.1000G>A', 'GENE', { alt: 'T', codonpos: 1, pos: 306, ref: 'A' })
    ];
    assert.deepEqual(h.proteinSubstitutionAlleles(sameNumber, 'GENE', { ref: 'A', pos: 306, alt: 'V' }).codon,
        { chrom: 'chr1', start: 1000, end: 1002, transcript: '' });
    // A codon split by an intron (its bases are not next to one another): not guessed at either.
    const split = [
        hit('chr1:g.1000A>G', 'GENE', { alt: 'A', codonpos: 1, pos: 7, ref: 'T' }),
        hit('chr1:g.2500C>T', 'GENE', { alt: 'M', codonpos: 2, pos: 7, ref: 'T' })
    ];
    assert.equal(h.proteinSubstitutionAlleles(split, 'GENE', { ref: 'T', pos: 7, alt: 'M' }).codon, null);
    assert.equal(h.proteinSubstitutionAlleles([], 'GENE', { ref: 'T', pos: 7, alt: 'M' }).codon, null);

    // What MyVariant is asked: the residue, any new amino acid.
    assert.equal(h.proteinResidueQuery('KRAS', { ref: 'Y', pos: 64, alt: '*' }), 'dbnsfp.genename:KRAS AND dbnsfp.aa.ref:Y AND dbnsfp.aa.pos:64');
    assert.equal(h.proteinResidueQuery('TP53', { ref: '*', pos: 394, alt: 'L' }), 'dbnsfp.genename:TP53 AND dbnsfp.aa.ref:X AND dbnsfp.aa.pos:394');
    for (const [gene, change] of [['KRAS AND x', { ref: 'Y', pos: 64, alt: '*' }], ['NP_004976.2', { ref: 'Y', pos: 0, alt: '*' }],
        ['KRAS', { ref: 'Tyr', pos: 64, alt: '*' }], ['KRAS', { ref: 'Y', pos: 6.5, alt: '*' }], ['KRAS', null], ['', { ref: 'Y', pos: 64, alt: '*' }]]) {
        assert.equal(h.proteinResidueQuery(gene, change), '', JSON.stringify([gene, change]));
    }

    // The alleles go the way an answer from Ensembl goes: every one of them is for the name asked.
    const shaped = h.recoderShapedFromAlleles('KRAS:p.Tyr64Ter', ['chr12:g.25380266G>C', 'chr12:g.25380266G>T']);
    assert.deepEqual(h.selectRecoderAlleles(shaped, 'KRAS:p.Tyr64Ter').map(allele => allele.myvariant),
        [['chr12:g.25380266G>C'], ['chr12:g.25380266G>T']]);
    console.log('protein-substitution unit checks passed');
}

const note = page => [...page.document.querySelectorAll('.protein-name-note')].map(el => el.textContent.replace(/\s+/g, ' ').trim());
const allelesNote = page => [...page.document.querySelectorAll('.protein-alleles-note')].map(el => el.textContent.replace(/\s+/g, ' ').trim());
const residueQueries = page => page.requests.filter(r => /dbnsfp\.aa\.pos/.test(r));
// Ensembl's recoder asked and never answering: the request simply stays out.
// (A page that waited for it would wait for ever; that is a failure here, not a hang.)
const never = new Promise(() => {});
const recoderSilent = (url, key) => (/\/variant_recoder\//.test(key) ? never : undefined);
const unwaited = (label, lookup) => Promise.race([lookup, new Promise((resolve, reject) => {
    const timer = setTimeout(() => reject(new Error(`${label}: the lookup is still waiting for Ensembl`)), 8000);
    lookup.then(() => clearTimeout(timer), () => clearTimeout(timer));
})]);

// 1. "KRAS Y64*": two single-base changes make Tyr64 a stop, MyVariant lists
//    both, and the best documented is shown — with Ensembl saying nothing at
//    all. (It used to be 13 s when Ensembl answered, and 54 s for the codon
//    alone when it did not.)
{
    const page = await unwaited('KRAS Y64*', runLookup(cassettes.krasY64Stop, { route: recoderSilent }));
    assertResolved(assert, page, 'KRAS Y64*, Ensembl silent');
    assert.match(page.variant, /g\.: chr12:g\.25380266G>C/);
    assert.match(page.variant, /c\.: c\.192C>Gp\.: p\.Tyr64Ter/);
    assert.deepEqual(allelesNote(page), ['p.Tyr64Ter is 2 alleles on record. This is the best documented one; the other is chr12:g.25380266G>T.']);
    assert.deepEqual(residueQueries(page), ['myvariant.info/v1/query?q=dbnsfp.genename:KRAS AND dbnsfp.aa.ref:Y AND dbnsfp.aa.pos:64'
        + '&fields=dbnsfp.genename,dbnsfp.aa.ref,dbnsfp.aa.alt,dbnsfp.aa.pos,dbnsfp.aa.codonpos,dbnsfp.ensembl.transcriptid&size=60']);
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*not waited for — MyVariant lists the 2 single-base changes that give this protein change/);
    assert.doesNotMatch(page.variant, /Codon region/);
    // With Ensembl answering, the alleles are still MyVariant's — the same two,
    // in the same order, whatever order Ensembl lists them in.
    const key = 'grch37.rest.ensembl.org/variant_recoder/human/KRAS:p.Tyr64Ter?content-type=application/json';
    const recorded = cassettes.krasY64Stop.responses[key].body;
    const alleles = Object.keys(recorded[0]).filter(k => k !== 'warnings');
    for (const order of [alleles, [...alleles].reverse()]) {
        const body = [Object.fromEntries(order.map(a => [a, recorded[0][a]]))];
        const answered = await runLookup(cassettes.krasY64Stop, { route: (url, k, respond) => (k === key ? respond(body) : undefined) });
        assertResolved(assert, answered, `KRAS Y64*, Ensembl listing ${order}`);
        assert.match(answered.variant, /g\.: chr12:g\.25380266G>C/, String(order));
        assert.deepEqual(allelesNote(answered), allelesNote(page), String(order));
        assert.equal(residueQueries(answered).length, 1);
    }
    console.log('ok protein substitution: KRAS Y64* is resolved from MyVariant\'s list, Ensembl or no Ensembl');
}

// 2. "KRAS Q61X": one base makes Gln61 a stop. With Ensembl down this was the
//    codon and a note that the stop could not be looked up.
{
    const page = await runLookup(cassettes.krasQ61X, { without: ['variant_recoder'] });
    assertResolved(assert, page, 'KRAS Q61X, Ensembl down');
    assert.match(page.variant, /g\.: chr12:g\.25380277G>A/);
    assert.match(page.variant, /c\.: c\.181C>T/);
    assert.deepEqual(note(page), ['"X" was read as a stop codon (p.Gln61Ter). If it stands for any amino acid, look up KRAS p.Gln61Xaa.']);
    assert.doesNotMatch(page.variant, /Codon region/);
    console.log('ok protein substitution: KRAS Q61X is the stop with Ensembl down');
}

// 3. "KRAS G12X": no single base makes Gly12 a stop, and that settles how X is
//    read. Ensembl is not waited for: 11 s on a good day, 55 s on a bad one.
{
    const page = await unwaited('KRAS G12X', runLookup(cassettes.krasG12Stop, { input: 'KRAS G12X', route: recoderSilent }));
    assertResolved(assert, page, 'KRAS G12X, Ensembl silent');
    assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.match(page.variant, /p\.: p\.Gly12Xaa/);
    // (What Ensembl's alleles would have added — that no larger change is on
    // record either — is not claimed.)
    assert.deepEqual(note(page), ['"X" was read as any amino acid (p.Gly12Xaa), not as a stop: no single-base change turns Gly12 into a stop codon.']);
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*not waited for — MyVariant shows that no single-base change makes a stop here, which settles how "X" is read/);
    assert(!page.requests.some(r => /\/vep\//.test(r)), 'no consequence is asked for an allele nobody named');
    for (const asked of page.api('pubmed')) assert.match(asked, /term=KRAS G12X&gene=KRAS&variant=G12X(&|$)/);
    // Not even for the one attempt a name without "X" would get: with the
    // page's waits running at their real length, the lookup is over at once.
    const started = Date.now();
    const realTime = await unwaited('KRAS G12X, real timers',
        runLookup(cassettes.krasG12Stop, { input: 'KRAS G12X', route: recoderSilent, timerCap: 60000 }));
    assertResolved(assert, realTime, 'KRAS G12X, real timers');
    assert.deepEqual(note(realTime), note(page));
    assert(Date.now() - started < 5000, `KRAS G12X took ${Date.now() - started} ms with Ensembl silent`);
    console.log('ok protein substitution: KRAS G12X is read as any amino acid without waiting for Ensembl');
}

// 4. "KRAS G12*": a stop that takes more than one base. Ensembl may name the
//    changes, and has the time of one attempt to; after that the codon is shown,
//    with the reason. (It used to be four attempts.)
{
    const silent = await unwaited('KRAS G12*', runLookup(cassettes.krasG12Stop, { input: 'KRAS G12*', route: recoderSilent }));
    assertResolved(assert, silent, 'KRAS G12*, Ensembl silent');
    assert.match(silent.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.match(silent.variant, /p\.: p\.Gly12Ter/);
    assert.deepEqual(note(silent), ['No single-base change gives p.Gly12Ter, and which larger changes do could not be looked up, so only the codon is shown.']);
    assert.match(silent.panel, /Resolve nomenclature · Ensembl\s*not waited for further — no single-base change gives this protein change, and Ensembl has not said which larger ones do/);
    // Answering within that time, it is used as before: three changes, each a link.
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    const slow = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12*', timerCap: 6000,
        route: async (url, key) => { if (/\/variant_recoder\//.test(key)) await wait(40); return undefined; } });
    assertResolved(assert, slow, 'KRAS G12*, Ensembl answering after a while');
    assert.match(note(slow)[0], /^p\.Gly12Ter can be any of 3 nucleotide changes and none of them is on record, so no single one is shown: /);
    console.log('ok protein substitution: KRAS G12* waits one attempt for Ensembl, not four');
}

// 5. Where MyVariant has nothing to say, it is Ensembl as before: the residue
//    is not in dbNSFP under that number, or MyVariant does not answer.
{
    const unknownResidue = await runLookup(cassettes.krasY64Stop, { without: ['dbnsfp.aa.pos'] });
    assertResolved(assert, unknownResidue, 'KRAS Y64*, residue not listed');
    assert.match(unknownResidue.variant, /g\.: chr12:g\.25380266G>C/);
    assert.match(unknownResidue.panel, /Resolve nomenclature · Ensembl\s*\d+ transcripts/);
    const failing = await runLookup(cassettes.krasY64Stop, {
        route: (url, key, respond) => (/dbnsfp\.aa\.pos/.test(key) ? respond({ error: 'Service Unavailable' }, 503) : undefined) });
    assertResolved(assert, failing, 'KRAS Y64*, residue query failing');
    assert.match(failing.variant, /g\.: chr12:g\.25380266G>C/);
    assert.match(failing.panel, /Resolve nomenclature · Ensembl\s*\d+ transcripts/);
    // ...and with neither, the codon, after Ensembl's attempts — as it always was.
    const neither = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12*', without: ['dbnsfp.aa.pos', 'variant_recoder'] });
    assertResolved(assert, neither, 'KRAS G12*, neither source');
    assert.match(neither.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.deepEqual(note(neither), [], 'nothing is known about single-base changes, so nothing is said of them');
    console.log('ok protein substitution: without MyVariant\'s list it is Ensembl, as before');
}

// 5b. With Ensembl not answering anything at all — neither the alleles nor
//     where the codon is — the residue's own records place the codon. "KRAS
//     G12K", a change of two bases, ended in "Could not resolve … to a genomic
//     position" on such a day.
{
    const noEnsembl = [/grch37\.rest\.ensembl\.org/, /\/\/rest\.ensembl\.org|^rest\.ensembl\.org/];
    for (const [typed, protein, said] of [
        ['KRAS G12*', 'p.Gly12Ter', 'No single-base change gives p.Gly12Ter, and which larger changes do could not be looked up, so only the codon is shown.'],
        ['KRAS G12X', 'p.Gly12Xaa', '"X" was read as any amino acid (p.Gly12Xaa), not as a stop: no single-base change turns Gly12 into a stop codon.']]) {
        const page = await runLookup(cassettes.krasG12Stop, { input: typed, without: noEnsembl });
        assertResolved(assert, page, `${typed}, Ensembl down altogether`);
        assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/, typed);
        assert.match(page.variant, new RegExp(`p\\.: ${protein.replace('.', '\\.')}`), typed);
        assert.deepEqual(note(page), [said], typed);
        assert.match(page.variant, /Resolved to residue 12 on ENST00000256078 only\./, typed);
        assert.match(page.panel, /codon placed from MyVariant's records of the residue — Ensembl did not answer/, typed);
        assert(!page.requests.some(r => /map\/translation/.test(r)), `${typed}: nothing was left to ask Ensembl's translation map`);
    }
    // Where Ensembl does place the codon, it is Ensembl's: the same span.
    const answered = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12X', route: recoderSilent });
    assert.match(answered.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.doesNotMatch(answered.panel, /codon placed from MyVariant/);
    console.log('ok protein substitution: with Ensembl down altogether, the residue\'s records place the codon');
}

// 6. Nothing is asked that is not needed: a name the name search finds, a name
//    that is not a substitution, a coding name, a name on an accession.
{
    for (const [cassette, input] of [[cassettes.brafV600E, undefined], [cassettes.cftrF508delName, undefined],
        [cassettes.kras35, undefined], [cassettes.krasDup, 'NM_004985.5:c.34_36dup']]) {
        const page = await runLookup(cassette, input ? { input } : {});
        assertResolved(assert, page, input || cassette.input);
        assert.deepEqual(residueQueries(page), [], `${input || cassette.input}: MyVariant is not asked for a residue`);
    }
    console.log('ok protein substitution: the residue is asked for only where the alleles are not yet known');
}

// 7. A lookup overtaken while MyVariant is being asked for the residue stops
//    there: it does not go on to fetch the alleles it is told of.
{
    const first = cassettes.krasY64Stop;
    const second = cassettes.brafV600E;
    const both = { input: first.input, reference: [...first.reference, ...second.reference],
        api: { ...(first.api || {}), ...(second.api || {}) }, responses: { ...first.responses, ...second.responses } };
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    let submitted = false;
    let residueAnswered = false;
    const page = await runLookup(both, {
        settleMs: 400,
        route: async (url, key, respond, options, { document }) => {
            if (!submitted && /dbnsfp\.aa\.pos/.test(key)) {
                submitted = true;
                document.getElementById('variantInput').value = 'BRAF V600E';
                document.getElementById('variantForm').dispatchEvent(new document.defaultView.Event('submit', { cancelable: true }));
                await wait(30);
                residueAnswered = true;
                return undefined;
            }
            // The second lookup is still under way when that answer arrives.
            if (submitted && /\/vep\/human\/hgvs\/7:g\.140453136A>T/.test(key)) {
                for (let i = 0; i < 200 && !residueAnswered; i++) await wait(2);
                await wait(20);
            }
            return undefined;
        }
    });
    assert(submitted, 'the second lookup was submitted while the first asked for the residue');
    assertResolved(assert, page, 'BRAF V600E, submitted during KRAS Y64*');
    assert.match(page.panel, /You entered\s*BRAF V600E/);
    assert.match(page.variant, /g\.: chr7:g\.140453136A>T/);
    assert.deepEqual(page.requests.filter(r => /myvariant\.info\/v1\/variant\/chr12:/.test(r)), [], 'the first lookup\'s alleles are not fetched');
    assert.doesNotMatch(page.panel + page.variant, /KRAS|Tyr64/, 'nothing of the lookup that was overtaken');
    assert.doesNotMatch(page.panel, /single-base change/, 'nor a line of its panel on this one');
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*\d+ transcripts/);
    console.log('ok protein substitution: a lookup overtaken while asking for the residue stops there');
}

console.log('Protein substitution allele tests passed');
