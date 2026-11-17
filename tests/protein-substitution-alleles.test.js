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
// MyVariant's answers to the residue query, as recorded live ("GENE REF POS").
const residues = JSON.parse(readFileSync(new URL('./fixtures/residue-records.json', import.meta.url)));
// The transcripts a name is counted on, as NCBI lists them for each gene: MANE
// Select (RefSeq Select for CDKN2A), then MANE Plus Clinical.
const STANDARD = {
    KRAS: ['ENST00000311936', 'ENST00000256078'], TP53: ['ENST00000269305'], PDGFRA: ['ENST00000257290'], UGT1A1: ['ENST00000305208'],
    MECP2: ['ENST00000453960', 'ENST00000303391'], CDKN2A: ['ENST00000304494', 'ENST00000579755'], EGFR: ['ENST00000275493'],
    PTEN: ['ENST00000371953'], 'MT-ND4': ['ENST00000361381']
};
const answer = (key, alt, standard) => {
    const [gene, ref, pos] = key.split(' ');
    return h.proteinSubstitutionAlleles(residues[key].hits, gene, { ref: ref === 'X' ? '*' : ref, pos: Number(pos), alt }, standard || STANDARD[gene]);
};
// (What is listed and whether the residue is known; the codon has its own checks below.)
const listed = (key, alt, standard) => {
    const out = answer(key, alt, standard);
    return out && { alleles: out.alleles, residueKnown: out.residueKnown };
};
const hit = (id, genename, aa, transcriptid) => ({ _id: id, dbnsfp: { genename, aa, ...(transcriptid ? { ensembl: { transcriptid } } : {}) } });
{
    // KRAS Tyr64 (TAC, read backwards): the two changes of the third base make a stop.
    assert.deepEqual(listed('KRAS Y 64', '*'), { alleles: ['chr12:g.25380266G>C', 'chr12:g.25380266G>T'], residueKnown: true });
    assert.deepEqual(h.proteinSubstitutionAlleles([...residues['KRAS Y 64'].hits].reverse(), 'kras', { ref: 'Y', pos: 64, alt: '*' }, STANDARD.KRAS).alleles,
        ['chr12:g.25380266G>C', 'chr12:g.25380266G>T'], 'in one order, whatever MyVariant\'s — and whatever case the symbol is in');
    assert.deepEqual(listed('KRAS Y 64', 'H'), { alleles: ['chr12:g.25380268A>G'], residueKnown: true });
    assert.deepEqual(listed('KRAS Y 64', 'F'), { alleles: ['chr12:g.25380267T>A'], residueKnown: true });
    // No single base makes Tyr64 a tryptophan: the residue is known, the list empty.
    assert.deepEqual(listed('KRAS Y 64', 'W'), { alleles: [], residueKnown: true });
    // Which changes give an amino acid follows from the codon, not from what
    // is listed: dbNSFP holds no record of a silent change, and TAT is still a
    // tyrosine.
    assert(!residues['KRAS Y 64'].hits.some(record => record._id === 'chr12:g.25380266G>A'));
    assert.deepEqual(listed('KRAS Y 64', 'Y'), { alleles: ['chr12:g.25380266G>A'], residueKnown: true });
    // Another residue, another reference amino acid, another gene: not known — nothing is concluded.
    const tyr64 = residues['KRAS Y 64'].hits;
    for (const [gene, change, standard] of [['KRAS', { ref: 'Y', pos: 65, alt: '*' }, STANDARD.KRAS], ['KRAS', { ref: 'F', pos: 64, alt: '*' }, STANDARD.KRAS],
        ['NRAS', { ref: 'Y', pos: 64, alt: '*' }, ['ENST00000369535']], ['NRAS', { ref: 'Y', pos: 64, alt: '*' }, STANDARD.KRAS]]) {
        assert.deepEqual(h.proteinSubstitutionAlleles(tyr64, gene, change, standard), { alleles: [], residueKnown: false, codon: null }, JSON.stringify([gene, change]));
    }
    for (const none of [[], null, undefined]) {
        assert.deepEqual(h.proteinSubstitutionAlleles(none, 'KRAS', { ref: 'Y', pos: 64, alt: '*' }, STANDARD.KRAS), { alleles: [], residueKnown: false, codon: null });
    }
    // Which transcripts the name is counted on has to be known: otherwise nothing is said.
    for (const unknown of [null, undefined, [], [''], ['NM_004985.5']]) {
        assert.equal(h.proteinSubstitutionAlleles(tyr64, 'KRAS', { ref: 'Y', pos: 64, alt: '*' }, unknown), null, JSON.stringify(unknown));
    }
    // (An id with its version is the transcript all the same.)
    assert.deepEqual(h.proteinSubstitutionAlleles(tyr64, 'KRAS', { ref: 'Y', pos: 64, alt: '*' }, ['ENST00000311936.8']).alleles,
        ['chr12:g.25380266G>C', 'chr12:g.25380266G>T']);

    // ── The number has to be the residue's on a transcript the name is counted on ──
    // A record lists the residue's number on every transcript it lies in, and
    // which transcript each number is on. Any of them used to count. (Found in
    // review; each of these was shown as the variant.)
    //
    // TP53: the Arg213 codon is numbered 202 on ENST00000615910, a fragment.
    // "TP53 R202X" came out as p.Arg213Ter. Arg202 itself (CGT) cannot become
    // a stop by one base.
    const arg213Stop = residues['TP53 R 202'].hits.find(record => record._id === 'chr17:g.7578212G>A');
    assert.equal(arg213Stop.dbnsfp.aa.pos[arg213Stop.dbnsfp.ensembl.transcriptid.indexOf('ENST00000615910')], 202);
    assert.equal(arg213Stop.dbnsfp.aa.pos[arg213Stop.dbnsfp.ensembl.transcriptid.indexOf('ENST00000269305')], 213);
    assert.deepEqual(listed('TP53 R 202', '*'), { alleles: [], residueKnown: true });
    assert.deepEqual(listed('TP53 R 202', 'C'), { alleles: ['chr17:g.7578245G>A'], residueKnown: true });
    assert.deepEqual(listed('TP53 R 202', '*', ['ENST00000615910']), { alleles: ['chr17:g.7578212G>A'], residueKnown: true },
        'whose residue 202 that allele is');
    // "TP53 W53G" came out as p.Trp146Gly, "the best documented" of two.
    assert.deepEqual(listed('TP53 W 53', 'G'), { alleles: ['chr17:g.7579530A>C'], residueKnown: true });
    assert(residues['TP53 W 53'].hits.some(record => record._id === 'chr17:g.7578494A>C'), 'the other is in MyVariant\'s answer');
    // PDGFRA: Val607's codon is numbered 367 on a transcript of ANOTHER gene
    // the record is shared with. "PDGFRA V367A" came out as p.Val607Ala.
    const val607 = residues['PDGFRA V 367'].hits.find(record => record._id === 'chr4:g.55143588T>C');
    assert.deepEqual([val607.dbnsfp.genename, val607.dbnsfp.aa.pos], [['AC058822.1', 'PDGFRA'], [367, 607]]);
    assert.deepEqual(listed('PDGFRA V 367', 'A'), { alleles: ['chr4:g.55133887T>C'], residueKnown: true });
    // UGT1A1: Lys520's codon is Lys521 of UGT1A3, 4 and 5, which share the
    // exon. "UGT1A1 K521Q" came out as p.Lys520Gln.
    const lys520 = residues['UGT1A1 K 521'].hits.find(record => record._id === 'chr2:g.234681161A>C');
    assert.deepEqual([lys520.dbnsfp.genename.slice(-3), lys520.dbnsfp.aa.pos.slice(-4)], [['UGT1A4', 'UGT1A3', 'UGT1A1'], [521, 521, 521, 520]]);
    assert.deepEqual(listed('UGT1A1 K 521', 'Q'), { alleles: ['chr2:g.234681164A>C'], residueKnown: true });
    // MECP2's T158M is counted on its MANE Plus Clinical transcript, where the
    // clinic counts; on MANE Select the residue is Thr170.
    assert.deepEqual(listed('MECP2 T 158', 'M'), { alleles: ['chrX:g.153296806G>A'], residueKnown: true });
    assert.deepEqual(listed('MECP2 T 158', 'M', ['ENST00000453960']), { alleles: [], residueKnown: false });

    // ── A variant in two reading frames ─────────────────────────────────────
    // ...comes as two amino-acid changes, and the transcripts it names are
    // those of the first. R306C is chrX:g.153296363G>A, which is A301V in
    // another frame; the residue has to be right within ONE frame, and the
    // record places a base of the codon it fits.
    const r306c = residues['MECP2 R 306'].hits.find(record => record._id === 'chrX:g.153296363G>A');
    assert.deepEqual(r306c.dbnsfp.aa.map(aa => `${aa.ref}${[].concat(aa.pos).join('/')}${aa.alt}`), ['A301V', 'R306/318C']);
    assert.deepEqual([].concat(r306c.dbnsfp.ensembl.transcriptid), ['ENST00000622433'], 'the transcript of A301V only');
    assert.deepEqual(listed('MECP2 R 306', 'C'), { alleles: ['chrX:g.153296363G>A'], residueKnown: true });
    assert.deepEqual(listed('MECP2 R 306', 'H'), { alleles: ['chrX:g.153296362C>T'], residueKnown: true });
    // (Residue 306 is an arginine on both transcripts — two different codons.
    // One base makes the one a stop, and none the other.)
    assert.deepEqual(listed('MECP2 R 306', '*'), { alleles: ['chrX:g.153296399G>A'], residueKnown: true });
    // "R301C" and "A306V" are in no frame, though MyVariant's search, which
    // matches across frames, returns the record for either.
    const mecp2 = residues['MECP2 R 306'].hits;
    assert.deepEqual(h.proteinSubstitutionAlleles(mecp2, 'MECP2', { ref: 'R', pos: 301, alt: 'C' }, STANDARD.MECP2), { alleles: [], residueKnown: false, codon: null });
    assert.deepEqual(h.proteinSubstitutionAlleles(mecp2, 'MECP2', { ref: 'A', pos: 306, alt: 'V' }, STANDARD.MECP2), { alleles: [], residueKnown: false, codon: null });
    // CDKN2A's exon 2 is read in two frames throughout (p16 and p14ARF).
    assert.deepEqual(listed('CDKN2A P 114', 'L'), { alleles: ['chr9:g.21971017G>A'], residueKnown: true });
    assert.deepEqual(listed('CDKN2A P 114', 'S'), { alleles: ['chr9:g.21971018G>A'], residueKnown: true });
    // A record that names no transcript and fits no codon placed is not counted in.
    const kras = ['KRAS', 'KRAS'];
    const gly12 = codonpos => ({ codonpos: [codonpos, codonpos], pos: [12, 12], ref: 'G', refcodon: ['GGT', 'GGT'] });
    const placedOnly = [hit('chr12:g.25398285C>A', kras, { alt: 'C', ...gly12(1) }, STANDARD.KRAS)];
    for (const stray of [
        hit('chr12:g.25398290C>T', 'KRAS', [{ alt: 'V', codonpos: 2, pos: 99, ref: 'A', refcodon: 'GCT' }, { alt: 'D', codonpos: 2, pos: 12, ref: 'G', refcodon: 'GGT' }]),
        hit('chr12:g.25398284G>A', 'KRAS', [{ alt: 'V', codonpos: 2, pos: 99, ref: 'A', refcodon: 'GCT' }, { alt: 'D', codonpos: 2, pos: 12, ref: 'G', refcodon: 'GGT' }]),
        hit('chr12:g.25398284C>T', 'KRAS', [{ alt: 'V', codonpos: 2, pos: 99, ref: 'A', refcodon: 'GCT' }, { alt: 'D', codonpos: 2, pos: 12, ref: 'G', refcodon: 'GGC' }])
    ]) {
        assert.equal(h.proteinSubstitutionAlleles(placedOnly.concat(stray), 'KRAS', { ref: 'G', pos: 12, alt: 'D' }, STANDARD.KRAS), null,
            `${stray._id}: the base of the change is not placed, so nothing is said`);
    }
    assert.deepEqual(h.proteinSubstitutionAlleles(placedOnly.concat(
        hit('chr12:g.25398284C>T', 'KRAS', [{ alt: 'V', codonpos: 2, pos: 99, ref: 'A', refcodon: 'GCT' }, { alt: 'D', codonpos: 2, pos: 12, ref: 'G', refcodon: 'GGT' }])),
    'KRAS', { ref: 'G', pos: 12, alt: 'D' }, STANDARD.KRAS).alleles, ['chr12:g.25398284C>T'], 'one that fits places its base');

    // A stop that is lost: dbNSFP writes the stop as X on either side.
    assert.deepEqual(listed('TP53 X 394', 'L'), { alleles: ['chr17:g.7572928C>A'], residueKnown: true });
    // Only single-base changes are read, and only what the record itself bears out.
    const one = ['ENST00000000001'];
    const gene5 = (id, aa) => hit(id, 'GENE', { pos: 5, ref: 'Y', refcodon: 'TAC', ...aa }, one);
    assert.deepEqual(h.proteinSubstitutionAlleles([gene5('chr1:g.1000T>C', { alt: 'H', codonpos: 1 }), gene5('chr1:g.1000T>C', { alt: 'H', codonpos: 1 }),
        gene5('chr1:g.1000_1001delinsGG', { alt: 'H', codonpos: 1 }), gene5(undefined, { alt: 'H', codonpos: 1 })],
    'GENE', { ref: 'Y', pos: 5, alt: 'H' }, one).alleles, ['chr1:g.1000T>C']);
    // A symbol written with small letters is the same gene.
    assert.deepEqual(h.proteinSubstitutionAlleles([hit('chr9:g.27567117A>G', ['C9orf72'], { alt: 'T', codonpos: [2], pos: [1], ref: 'M', refcodon: ['ATG'] }, ['ENST00000380003'])],
        'C9ORF72', { ref: 'M', pos: 1, alt: 'T' }, ['ENST00000380003']).alleles, ['chr9:g.27567117A>G']);
    // Records that contradict the code say nothing: a reference base that is
    // on neither strand of its codon, a codon that is not the amino acid, a
    // change listed that the code does not give.
    for (const [label, records] of [
        ['base on neither strand', [gene5('chr1:g.1000G>C', { alt: 'H', codonpos: 1 })]],
        ['codon of another amino acid', [hit('chr1:g.1000T>C', 'GENE', { alt: 'H', codonpos: 1, pos: 5, ref: 'Y', refcodon: 'TGC' }, one)]],
        ['a listed change the code does not give', [gene5('chr1:g.1000T>C', { alt: 'W', codonpos: 1 })]],
        ['no codon given', [hit('chr1:g.1000T>C', 'GENE', { alt: 'H', pos: 5, ref: 'Y' }, one)]],
        ['one base at two places', [gene5('chr1:g.1000T>C', { alt: 'H', codonpos: 1 }), gene5('chr1:g.1900T>C', { alt: 'H', codonpos: 1 })]]
    ]) {
        const alt = label === 'a listed change the code does not give' ? 'W' : 'H';
        assert.equal(h.proteinSubstitutionAlleles(records, 'GENE', { ref: 'Y', pos: 5, alt }, one), null, label);
    }

    // ── What is said when no single base gives the change ───────────────────
    // KRAS Gly12 is GGT: no single base makes it a stop — the codon says so.
    assert.deepEqual(listed('KRAS G 12', '*'), { alleles: [], residueKnown: true });
    assert.deepEqual(listed('KRAS G 12', 'K'), { alleles: [], residueKnown: true });
    assert.deepEqual(listed('KRAS G 12', 'D'), { alleles: ['chr12:g.25398284C>T'], residueKnown: true });
    // "KRAS G12G": three changes of the third base leave it a glycine. dbNSFP
    // lists none of them (it has no silent changes), and the page concluded
    // "No single-base change gives p.Gly12Gly". No record places that base,
    // so where they lie is not said here: the name is left to Ensembl.
    assert.equal(answer('KRAS G 12', 'G'), null);

    // ── Where the codon lies ────────────────────────────────────────────────
    // Each record places one base of it: its own position, which base of the
    // codon it changes, and the strand (its reference base against the
    // codon's). KRAS Gly12 is read backwards — base 1 at 25398285, base 2 at
    // 25398284 — and no record is at base 3, where every change is silent.
    // That base is not assumed to follow: the span is the two bases placed.
    assert.deepEqual(answer('KRAS G 12', '*').codon, { chrom: 'chr12', start: 25398284, end: 25398285, bases: 2, transcript: 'ENST00000311936' });
    // All three, where each has a record.
    assert.deepEqual(answer('KRAS Y 64', '*').codon, { chrom: 'chr12', start: 25380266, end: 25380268, bases: 3, transcript: 'ENST00000311936' });
    assert.deepEqual(answer('TP53 W 53', 'G').codon, { chrom: 'chr17', start: 7579528, end: 7579530, bases: 3, transcript: 'ENST00000269305' });
    // EGFR Thr790 (ACG), read forwards.
    assert.deepEqual(answer('EGFR T 790', 'M').codon, { chrom: 'chr7', start: 55249070, end: 55249071, bases: 2, transcript: 'ENST00000275493' });
    // One record places one base.
    assert.deepEqual(h.proteinSubstitutionAlleles(residues['EGFR T 790'].hits.filter(record => record._id === 'chr7:g.55249071C>T'),
        'EGFR', { ref: 'T', pos: 790, alt: 'M' }, STANDARD.EGFR),
    { alleles: ['chr7:g.55249071C>T'], residueKnown: true, codon: { chrom: 'chr7', start: 55249071, end: 55249071, bases: 1, transcript: 'ENST00000275493' } });
    // PTEN Leu70 is CT|G with an intron after its second base, and no change
    // of its third is on record. The codon used to be given as
    // chr10:g.89685313_89685315, whose last base is the intron's first. (Found
    // in review.)
    assert.deepEqual(answer('PTEN L 70', '*'), { alleles: [], residueKnown: true,
        codon: { chrom: 'chr10', start: 89685313, end: 89685314, bases: 2, transcript: 'ENST00000371953' } });
    // One number, two codons — MECP2's residue 306 on its two transcripts: not one span.
    assert.equal(answer('MECP2 R 306', 'C').codon, null);
    // A codon an intron splits, with a record at each base: each allele is
    // where its record is, and the three are not one span.
    const split = [gene5('chr1:g.1000T>C', { alt: 'H', codonpos: 1 }), gene5('chr1:g.1001A>G', { alt: 'C', codonpos: 2 }),
        gene5('chr1:g.2500C>A', { alt: 'X', codonpos: 3 })];
    assert.deepEqual(h.proteinSubstitutionAlleles(split, 'GENE', { ref: 'Y', pos: 5, alt: '*' }, one),
        { alleles: ['chr1:g.2500C>A', 'chr1:g.2500C>G'], residueKnown: true, codon: null });

    // ── Mitochondrial genes are left to Ensembl ─────────────────────────────
    // MyVariant files dbNSFP's chrM records by another reference's positions
    // than the rest of this page uses: m.11778G>A, MT-ND4 R340H, is
    // "chrMT:g.11779G>A" there, and was shown as that — one base off, "p.:
    // N/A", VEP refusing it. (Found in review.)
    assert(residues['MT-ND4 R 340'].hits.some(record => record._id === 'chrMT:g.11779G>A'));
    assert.deepEqual(answer('MT-ND4 R 340', 'H'), { alleles: [], residueKnown: false, codon: null });
    assert.equal(h.proteinResidueQuery('MT-ND4', { ref: 'R', pos: 340, alt: 'H' }), '');
    assert.equal(h.proteinResidueQuery('mt-atp6', { ref: 'L', pos: 5, alt: 'P' }), '');

    // What MyVariant is asked: the residue, any new amino acid.
    assert.equal(h.proteinResidueQuery('KRAS', { ref: 'Y', pos: 64, alt: '*' }), 'dbnsfp.genename:KRAS AND dbnsfp.aa.ref:Y AND dbnsfp.aa.pos:64');
    assert.equal(h.proteinResidueQuery('TP53', { ref: '*', pos: 394, alt: 'L' }), 'dbnsfp.genename:TP53 AND dbnsfp.aa.ref:X AND dbnsfp.aa.pos:394');
    for (const [gene, change] of [['KRAS AND x', { ref: 'Y', pos: 64, alt: '*' }], ['NP_004976.2', { ref: 'Y', pos: 0, alt: '*' }],
        ['KRAS', { ref: 'Tyr', pos: 64, alt: '*' }], ['KRAS', { ref: 'Y', pos: 6.5, alt: '*' }], ['KRAS', null], ['', { ref: 'Y', pos: 64, alt: '*' }]]) {
        assert.equal(h.proteinResidueQuery(gene, change), '', JSON.stringify([gene, change]));
    }

    // Which transcripts a gene's names are counted on: NCBI's list, by
    // Ensembl's names. null when it cannot say — nothing is then concluded.
    assert.deepEqual(await h.fetchStandardEnsemblTranscripts('KRAS', async () => ({ ensembl: STANDARD.KRAS })), STANDARD.KRAS);
    assert.deepEqual(await h.fetchStandardEnsemblTranscripts('NOGENE', async () => ({ ensembl: [] })), []);
    assert.equal(await h.fetchStandardEnsemblTranscripts('KRAS', async () => null), null);
    assert.equal(await h.fetchStandardEnsemblTranscripts('KRAS', async () => { throw new Error('down'); }), null);

    // The alleles go the way an answer from Ensembl goes: every one of them is for the name asked.
    const shaped = h.recoderShapedFromAlleles('KRAS:p.Tyr64Ter', ['chr12:g.25380266G>C', 'chr12:g.25380266G>T']);
    assert.deepEqual(h.selectRecoderAlleles(shaped, 'KRAS:p.Tyr64Ter').map(allele => allele.myvariant),
        [['chr12:g.25380266G>C'], ['chr12:g.25380266G>T']]);
    // ...but they are not Ensembl's answer, and are not shown or passed on as
    // one. The page's "Show Ensembl recoder result" held this list, in the
    // shape of an answer, and the AI review was given it as `ensembl_recoder`.
    // (Found in review.)
    const ensembls = [{ T: { input: 'KRAS:p.Tyr64Ter', hgvsg: ['NC_000012.11:g.25380266G>T'] } }];
    const placedByNcbi = Object.assign([{ A: { input: 'NM_004333.4:c.1799T>A', hgvsg: ['NC_000007.13:g.140453136A>T'] } }], { __vsFallback: true });
    assert.equal(h.ensemblRecoderAnswer(shaped), null);
    assert.equal(h.ensemblRecoderAnswer(placedByNcbi), null);
    assert.equal(h.ensemblRecoderAnswer(ensembls), ensembls);
    assert.equal(h.ensemblRecoderAnswer(null), null);
    assert.equal(h.recoderPaneText(shaped, null), 'Ensembl\'s recoder was not waited for. MyVariant lists the single-base changes that give this protein change:\n'
        + JSON.stringify(['chr12:g.25380266G>C', 'chr12:g.25380266G>T'], null, 2));
    assert.equal(h.recoderPaneText(shaped, ensembls), JSON.stringify(ensembls, null, 2), 'Ensembl\'s own answer, where it came as well');
    assert.equal(h.recoderPaneText(ensembls, null), JSON.stringify(ensembls, null, 2));
    assert.match(h.recoderPaneText(placedByNcbi, null), /^Ensembl's recoder had not answered\. The variant was placed by NCBI Variation Services:\n\[/);
    assert.equal(h.recoderPaneText(null, null), '');
    assert.equal(h.recoderPaneText(null, placedByNcbi), '');
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
        + '&fields=dbnsfp.genename,dbnsfp.aa.ref,dbnsfp.aa.alt,dbnsfp.aa.pos,dbnsfp.aa.codonpos,dbnsfp.aa.refcodon,dbnsfp.ensembl.transcriptid&size=60']);
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*not waited for — MyVariant lists the 2 single-base changes that give this protein change/);
    assert.doesNotMatch(page.variant, /Codon region/);
    // What "Show Ensembl recoder result" holds is not passed off as Ensembl's.
    const pane = shown => shown.document.getElementById('ensemblOutput').textContent;
    assert.equal(pane(page), 'Ensembl\'s recoder was not waited for. MyVariant lists the single-base changes that give this protein change:\n'
        + JSON.stringify(['chr12:g.25380266G>C', 'chr12:g.25380266G>T'], null, 2));
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
        // (...and there the pane holds what Ensembl said.)
        assert.deepEqual(JSON.parse(pane(answered)), body, String(order));
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
    // A list MyVariant cut short is no list either: with one of the two
    // alleles missing from it, the other would have stood as "the" change.
    {
        const residueKey = Object.keys(cassettes.krasY64Stop.responses).find(key => /dbnsfp\.aa\.pos/.test(key));
        const whole = cassettes.krasY64Stop.responses[residueKey].body;
        const stopAlleles = h.proteinSubstitutionAlleles(whole.hits, 'KRAS', { ref: 'Y', pos: 64, alt: '*' }, STANDARD.KRAS).alleles;
        assert.deepEqual(stopAlleles, ['chr12:g.25380266G>C', 'chr12:g.25380266G>T']);
        const cut = { ...whole, total: whole.hits.length, hits: whole.hits.filter(hit => hit._id !== 'chr12:g.25380266G>C') };
        assert.equal(cut.hits.length, whole.hits.length - 1);
        const short = await runLookup(cassettes.krasY64Stop, { route: (url, key, respond) => (key === residueKey ? respond(cut) : undefined) });
        assertResolved(assert, short, 'KRAS Y64*, residue list cut short');
        assert.match(short.panel, /Resolve nomenclature · Ensembl\s*\d+ transcripts/, 'Ensembl was waited for');
        assert.doesNotMatch(short.panel, /not waited for/);
        assert.match(short.variant, /p\.Tyr64Ter is 2 alleles on record/, 'and both alleles are known');
    }
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
//
//     As much of it as they place, that is. Gly12 is GGT: every change of its
//     third base is silent, so no record is at that base, and it is not
//     assumed to follow the second — an intron may lie between, as it does in
//     PTEN's Leu70 (see the unit checks). The region is the two bases placed,
//     and the card says so.
{
    const noEnsembl = [/grch37\.rest\.ensembl\.org/, /\/\/rest\.ensembl\.org|^rest\.ensembl\.org/];
    for (const [typed, protein, said] of [
        ['KRAS G12*', 'p.Gly12Ter', 'No single-base change gives p.Gly12Ter, and which larger changes do could not be looked up, so only the codon is shown.'],
        ['KRAS G12X', 'p.Gly12Xaa', '"X" was read as any amino acid (p.Gly12Xaa), not as a stop: no single-base change turns Gly12 into a stop codon.']]) {
        const page = await runLookup(cassettes.krasG12Stop, { input: typed, without: noEnsembl });
        assertResolved(assert, page, `${typed}, Ensembl down altogether`);
        assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398284_25398285/, typed);
        assert.doesNotMatch(page.variant, /25398283/, typed);
        assert.match(page.variant, new RegExp(`p\\.: ${protein.replace('.', '\\.')}`), typed);
        assert.deepEqual(note(page), [said], typed);
        // (On the gene's MANE Select transcript, as everything else is.)
        assert.match(page.variant, /Resolved to residue 12 on ENST00000311936 only\./, typed);
        assert.match(page.variant, /The region is the 2 bases of the codon that MyVariant's records place: Ensembl, which maps the whole codon, did not answer\./, typed);
        assert.match(page.panel, /2 of the codon's three bases placed from MyVariant's records of the residue \(none places the other\) — Ensembl did not answer/, typed);
        assert(!page.requests.some(r => /map\/translation/.test(r)), `${typed}: nothing was left to ask Ensembl's translation map`);
    }
    // A codon whose three bases all have records is shown whole, and nothing is said of it.
    const whole = await runLookup(cassettes.krasY64Stop, { input: 'KRAS Y64W', without: noEnsembl });
    assertResolved(assert, whole, 'KRAS Y64W, Ensembl down altogether');
    assert.match(whole.variant, /Codon region \(hg19\): chr12:g\.25380266_25380268/);
    assert.doesNotMatch(whole.variant, /The region is/);
    assert.match(whole.panel, /codon placed from MyVariant's records of the residue — Ensembl did not answer/);
    // Where Ensembl does place the codon, it is Ensembl's: all three bases.
    const answered = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12X', route: recoderSilent });
    assert.match(answered.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.doesNotMatch(answered.panel, /placed from MyVariant/);
    assert.doesNotMatch(answered.variant, /The region is/);
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

// 8. The residue's number has to be its number on a transcript the name is
//    counted on. MyVariant's answer holds every record in which SOME
//    transcript has the residue under that number, and the first lookups took
//    them all. (Found in review; each of these showed another variant than
//    the one named, and main shows what is asserted here.)
{
    // "TP53 R202X": the Arg213 codon is residue 202 of a fragment transcript,
    // and one base makes it a stop — p.Arg213Ter was shown, as "X read as a
    // stop". Arg202 itself cannot become a stop by one base.
    const r202 = await unwaited('TP53 R202X', runLookup(cassettes.tp53R202X, { route: recoderSilent }));
    assertResolved(assert, r202, 'TP53 R202X');
    assert.match(r202.variant, /Codon region \(hg19\): chr17:g\.7578243_7578245/);
    assert.match(r202.variant, /p\.: p\.Arg202Xaa/);
    assert.deepEqual(note(r202), ['"X" was read as any amino acid (p.Arg202Xaa), not as a stop: no single-base change turns Arg202 into a stop codon.']);
    assert.doesNotMatch(r202.variant + r202.panel, /Arg213|7578212/);
    // "TP53 W53G": p.Trp146Gly was shown, as "the best documented" of two.
    const w53 = await unwaited('TP53 W53G', runLookup(cassettes.tp53W53G, { route: recoderSilent }));
    assertResolved(assert, w53, 'TP53 W53G');
    assert.match(w53.variant, /g\.: chr17:g\.7579530A>C.*c\.: c\.157T>Gp\.: p\.Trp53Gly \(p\.W53G\)Effect: missense_variant/);
    assert.deepEqual(allelesNote(w53), []);
    assert.doesNotMatch(w53.variant + w53.panel, /Trp146|7578494/);
    assert.match(w53.panel, /not waited for — MyVariant lists the single-base change that gives this protein change/);
    // "PDGFRA V367A": p.Val607Ala was shown — 367 on a transcript of another gene.
    const v367 = await unwaited('PDGFRA V367A', runLookup(cassettes.pdgfraV367A, { route: recoderSilent }));
    assertResolved(assert, v367, 'PDGFRA V367A');
    assert.match(v367.variant, /g\.: chr4:g\.55133887T>C.*c\.: c\.1100T>Cp\.: p\.Val367Ala/);
    assert.doesNotMatch(v367.variant + v367.panel, /Val607|55143588/);
    // Which transcripts those are NCBI says. When it cannot — or what it says
    // names none by Ensembl's id — nothing is taken from the residue's
    // records, and the name is Ensembl's to resolve, as before.
    for (const [label, options] of [
        ['NCBI not answering', { without: [/datasets\/v2\/gene\/symbol/] }],
        ['NCBI naming no transcript by Ensembl\'s id', { route: (url, key, respond) => (/datasets\/v2\/gene\/symbol\/KRAS\//.test(key)
            ? respond({ reports: [{ product: { symbol: 'KRAS', transcripts: [{ accession_version: 'NM_004985.5', select_category: 'MANE_SELECT' }] } }] })
            : undefined) }]
    ]) {
        const page = await runLookup(cassettes.krasY64Stop, options);
        assertResolved(assert, page, `KRAS Y64*, ${label}`);
        assert.match(page.variant, /g\.: chr12:g\.25380266G>C/, label);
        assert.match(page.panel, /Resolve nomenclature · Ensembl\s*\d+ transcripts/, `${label}: Ensembl was waited for`);
        assert.doesNotMatch(page.panel, /not waited for/, label);
    }
    console.log('ok protein substitution: the residue is counted on the gene\'s standard transcripts only');
}

// 9. Mitochondrial names are left to Ensembl. MyVariant files dbNSFP's chrM
//    records by another reference's positions: "MT-ND4 R340H" — m.11778G>A,
//    the commonest LHON variant — came out as chrMT:g.11779G>A, "p.: N/A",
//    with VEP refusing the allele. (Found in review.)
{
    const page = await runLookup(cassettes.mtNd4R340H, {});
    assertResolved(assert, page, 'MT-ND4 R340H');
    assert.match(page.variant, /g\.: chrMT:g\.11778G>A.*c\.: c\.1019G>Ap\.: p\.Arg340His \(p\.R340H\)Effect: missense_variant/);
    assert.doesNotMatch(page.variant + page.panel, /11779/);
    assert.deepEqual(residueQueries(page), [], 'MyVariant is not asked for the residue');
    console.log('ok protein substitution: a mitochondrial name is resolved by Ensembl, at its own position');
}

// 10. "No single-base change gives …" is said only where the codon says so.
//     "KRAS G12G": GGT to GGC, GGA or GGG is still a glycine, but dbNSFP has
//     no silent changes, and "the residue is listed, none of its records
//     gives this" was taken for proof — with Ensembl's answer, which names
//     the three, in hand. (Found in review.)
{
    for (const [label, options] of [['Ensembl answering', {}], ['Ensembl\'s recoder down', { without: ['variant_recoder'] }]]) {
        const page = await runLookup(cassettes.krasG12G, options);
        assertResolved(assert, page, `KRAS G12G, ${label}`);
        assert.equal(residueQueries(page).length, 1, label);
        assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/, label);
        assert.match(page.variant, /p\.: p\.Gly12Gly/, label);
        assert.deepEqual(note(page), [], `${label}: nothing is said about single-base changes`);
        assert.doesNotMatch(page.variant + page.panel, /[Nn]o single-base change/, label);
    }
    console.log('ok protein substitution: nothing is concluded about a silent change dbNSFP does not list');
}

// 11. The same name looked up again is resolved as it was the first time.
//     Ensembl's answer to the first lookup arrives after that lookup has what
//     it needs ("answered after … — not needed") and stays on the page. The
//     second lookup took it for its own answer being in hand, did not set its
//     own call aside — and waited for it: 49 s when Ensembl then said nothing,
//     against none the first time. (Found in review.)
{
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    let recoderCalls = 0;
    let cardsBegun = false;
    let lateAnswerGiven = false;
    const lookup = runLookup(cassettes.krasY64Stop, {
        followUp: { input: 'KRAS Y64*' },
        route: async (url, key) => {
            if (/\/variant_recoder\//.test(key)) {
                recoderCalls += 1;
                // The second time, Ensembl says nothing at all.
                if (recoderCalls > 1) return never;
                // The first time it answers late: once the cards are being fetched.
                for (let i = 0; i < 2000 && !cardsBegun; i++) await wait(2);
                setTimeout(() => { lateAnswerGiven = true; }, 0);
                return undefined;
            }
            if (/\/api\/(?:civic|clinvar-region|pubmed)/.test(key) && recoderCalls === 1) {
                cardsBegun = true;
                // ...and the first lookup ends only after that answer is in.
                for (let i = 0; i < 2000 && !lateAnswerGiven; i++) await wait(2);
                await wait(10);
            }
            return undefined;
        }
    });
    const page = await unwaited('KRAS Y64*, looked up a second time', lookup);
    assert.equal(recoderCalls, 2, 'each lookup asked Ensembl');
    assert(lateAnswerGiven, 'the first answer came in late');
    assertResolved(assert, page, 'KRAS Y64*, a second time');
    assert.match(page.variant, /g\.: chr12:g\.25380266G>C/);
    assert.deepEqual(allelesNote(page), ['p.Tyr64Ter is 2 alleles on record. This is the best documented one; the other is chr12:g.25380266G>T.']);
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*not waited for — MyVariant lists the 2 single-base changes that give this protein change/);
    console.log('ok protein substitution: a second lookup of the same name does not wait on the first one\'s late answer');
}

// 12. MyVariant answering for some of what it is asked. (Found in review.)
{
    // "KRAS Y64*" is two alleles, both on record. When the request for one of
    // the two records failed, the other was shown as THE variant, with nothing
    // said — c.192C>A or c.192C>G according to which request it was. It is
    // still shown (it is on record), and the card says what could not be asked.
    for (const [failing, shown] of [['chr12:g.25380266G>C', 'chr12:g.25380266G>T'], ['chr12:g.25380266G>T', 'chr12:g.25380266G>C']]) {
        const page = await unwaited(`KRAS Y64*, no answer for ${failing}`, runLookup(cassettes.krasY64Stop, {
            route: (url, key, respond) => (/\/variant_recoder\//.test(key) ? never
                : (key === `myvariant.info/v1/variant/${failing}` ? respond({ error: 'Service Unavailable' }, 503) : undefined)) }));
        assertResolved(assert, page, `KRAS Y64*, no answer for ${failing}`);
        assert(page.variant.startsWith(`Variantg.: ${shown}`), page.variant.slice(0, 80));
        assert.deepEqual(allelesNote(page), [`p.Tyr64Ter can be 2 nucleotide changes. This one is on record; MyVariant could not be asked about `
            + `the other, which may be better documented: ${failing}.`]);
        assert.deepEqual([...page.document.querySelectorAll('.protein-alleles-note a')].map(a => decodeURIComponent(a.getAttribute('href'))), [`?variant=${failing}`]);
        assert.match(page.panel, new RegExp(`MyVariant could not be asked about another allele that gives this protein change — showing ${shown.replace(/[.>]/g, '\\$&')}, which is on record`));
        assert.doesNotMatch(page.variant, /best documented one/);
    }
    // ...and with both answering, as before: the best documented, and the other named.
    const both = await unwaited('KRAS Y64*', runLookup(cassettes.krasY64Stop, { route: recoderSilent }));
    assert.deepEqual(allelesNote(both), ['p.Tyr64Ter is 2 alleles on record. This is the best documented one; the other is chr12:g.25380266G>T.']);

    // "KRAS Q61X": MyVariant lists the one change that makes Gln61 a stop —
    // and then does not answer for its record, while Ensembl has not answered
    // either. The allele is known all the same. It came out as the codon and
    // "Which nucleotide change a stop here (p.Gln61Ter) would be could not be
    // looked up".
    for (const [label, status, said] of [
        ['MyVariant not answering for the record', 503, /Annotations · MyVariant\.info\s*could not be asked for a record of this allele — it is shown without one/],
        ['MyVariant having no record', 404, null]
    ]) {
        const page = await unwaited(`KRAS Q61X, ${label}`, runLookup(cassettes.krasQ61X, {
            route: (url, key, respond) => (/\/variant_recoder\//.test(key) ? never
                : (/myvariant\.info\/v1\/variant\//.test(key)
                    ? respond(status === 404 ? { success: false, error: 'Not Found.' } : { error: 'Service Unavailable' }, status) : undefined)) }));
        assertResolved(assert, page, `KRAS Q61X, ${label}`);
        assert.match(page.variant, /g\.: chr12:g\.25380277G>A/, label);
        assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.181C>Tp\.: p\.Gln61Ter \(p\.Q61\*\)Effect: stop_gained/, label);
        assert.deepEqual(note(page), ['"X" was read as a stop codon (p.Gln61Ter). If it stands for any amino acid, look up KRAS p.Gln61Xaa.'], label);
        assert.doesNotMatch(page.variant, /Codon region|could not be looked up/, label);
        if (said) assert.match(page.panel, said, label);
        else assert.doesNotMatch(page.panel, /could not be asked for a record/, label);
    }
    console.log('ok protein substitution: an allele MyVariant could not be asked about is named, and one it listed is shown without its record');
}

console.log('Protein substitution allele tests passed');
