/*
 * A coding change named on a gene symbol.
 *
 * "PTEN c.512T>G" does not say which transcript c.512 is counted on. The page
 * left that to chance twice over: among MyVariant's hits, to whichever matched
 * on some transcript; and with no hit, to Ensembl, which picks a transcript of
 * its own and warns that it did. For PTEN c.512T>G Ensembl's pick is a variant
 * that is c.1026+59T>G — an intron — on every RefSeq transcript of the gene,
 * and that is what the page showed, without a word. It also does not answer
 * "BRCA1:c.5266dup" within a minute.
 *
 * A name on a gene symbol is now read on one of the transcripts the gene is
 * reported on — MANE Select, or where the gene has one, MANE Plus Clinical —
 * and from there on the lookup is for "NM_…:c.X".
 *
 * The second transcript matters. MECP2 is counted on NM_004992 in every clinic
 * letter (T158M is c.473C>T there, c.509C>T on MANE Select), NF1 on NM_000267,
 * MUTYH on NM_001128425. Reading such names on MANE Select alone — which this
 * rule first did — refused "MECP2 c.473C>T" and silently turned
 * "NF1 c.6789_6792del", a well-known pathogenic deletion, into four other
 * bases nobody has reported. (Scenarios 9 onwards; found in review.)
 *
 * The scenarios replay answers recorded from the live services.
 *
 * Run with: node tests/gene-symbol-name.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const asked = (page, pattern) => page.requests.filter(r => pattern.test(r));
const recoderQueries = page => asked(page, /\/variant_recoder\/human\//).map(r => r.match(/human\/([^?]+)/)[1]);

// ── What counts as a name on a gene symbol ─────────────────────────────────
assert.deepEqual(h.geneSymbolCodingName('BRCA1:c.5266dup'), { gene: 'BRCA1', change: 'c.5266dup' });
assert.deepEqual(h.geneSymbolCodingName('brca1 c.5266dupC'), { gene: 'BRCA1', change: 'c.5266dupC' });
assert.deepEqual(h.geneSymbolCodingName('NT5C2:c.1A>G'), { gene: 'NT5C2', change: 'c.1A>G' }, 'a symbol that starts like an accession prefix');
for (const other of ['NM_007294.4:c.5266dup', 'ENST00000357654:c.5266dup', 'LRG_292t1:c.5266dup', 'chr17:g.41209082dup',
    'BRAF:p.Val600Glu', 'CHR7:c.1A>G', 'rs80357906', 'BRCA1', '']) {
    assert.equal(h.geneSymbolCodingName(other), null, other);
}

// ── The transcript ClinVar names a gene's variants on ──────────────────────
{
    const named = (accession, gene = 'KRAS') => ({ clinvar: { rcv: [{ preferred_name: `${accession}(${gene}):c.35G>A (p.Gly12Asp)` }] } });
    const many = (n, accession, gene) => Array.from({ length: n }, () => named(accession, gene));
    // As recorded for KRAS: 17 records on the MANE Select, 8 on KRAS4A.
    assert.equal(h.clinvarNamingTranscript([...many(17, 'NM_004985.5'), ...many(8, 'NM_033360.4')], 'KRAS'), 'NM_004985.5');
    assert.equal(h.clinvarNamingTranscript([...many(13, 'NM_004985.5'), ...many(12, 'NM_033360.4')], 'KRAS'), '', 'no clear majority');
    assert.equal(h.clinvarNamingTranscript(many(4, 'NM_004985.5'), 'KRAS'), '', 'too few records to say');
    assert.equal(h.clinvarNamingTranscript(many(9, 'NM_004985.5', 'KRAS-AS1'), 'KRAS'), '', 'named in another gene');
    assert.equal(h.clinvarNamingTranscript([], 'KRAS'), '');
    assert.equal(h.clinvarNamingTranscript(null, 'KRAS'), '');
    // A record with several RCVs counts once; one with none, not at all.
    const twice = { clinvar: { rcv: [{ preferred_name: 'NM_004985.5(KRAS):c.35G>A' }, { preferred_name: 'NM_004985.5(KRAS):c.35G>A' }] } };
    assert.equal(h.clinvarNamingTranscript([twice, twice, twice, twice, {}, { clinvar: {} }], 'KRAS'), '', 'four records, not eight');
    // The recorded answer for PTEN.
    const recorded = cassettes.pten512.responses['myvariant.info/v1/query?q=clinvar.gene.symbol:PTEN AND _exists_:clinvar.rcv.preferred_name&fields=clinvar.rcv.preferred_name&size=25'];
    assert.equal(h.clinvarNamingTranscript(recorded.body.hits, 'PTEN'), 'NM_000314.8');
}

// ── The sources for the standard transcript ────────────────────────────────
{
    const timing = { maneAlone: 30, total: 120, clinical: 40 };
    const after = (ms, value) => () => new Promise(resolve => setTimeout(() => resolve(value), ms));
    const silent = async () => null;                       // NCBI not answering
    const std = (mane, naming, clinical = silent, gene = 'PTEN') => h.standardTranscriptForGene(gene, mane, naming, timing, clinical);
    // Ensembl's MANE record is an authority when it answers in time.
    assert.deepEqual(await std(after(1, ['NM_000314.8', 'ENST00000371953.8']), after(1, 'NM_999999.1')),
        { accession: 'NM_000314.8', basis: 'the MANE Select transcript of PTEN', plusClinical: null });
    // Slow: ClinVar's naming transcript stands in.
    assert.deepEqual(await std(after(500, ['NM_000314.8']), after(1, 'NM_000314.8')),
        { accession: 'NM_000314.8', basis: 'the transcript ClinVar names PTEN variants on', plusClinical: null });
    // ClinVar unclear: MANE gets the rest of the time.
    assert.deepEqual(await std(after(60, ['NM_000314.8']), after(1, '')),
        { accession: 'NM_000314.8', basis: 'the MANE Select transcript of PTEN', plusClinical: null });
    // Neither can say — including by failing.
    assert.equal(await std(after(1, []), after(1, '')), null);
    assert.equal(await std(() => { throw new Error('down'); }, () => Promise.reject(new Error('down')), () => { throw new Error('down'); }), null);
    assert.equal(await std(after(500, ['NM_000314.8']), after(500, 'NM_000314.8')), null, 'both too slow');
    // A gene with no MANE Select (only an Ensembl accession, or none) is not a RefSeq answer.
    assert.equal(await std(after(1, ['ENST00000371953.8']), after(1, '')), null);

    // NCBI's list is the other authority for MANE Select, and the only source
    // for MANE Plus Clinical. Whichever names MANE Select first is taken…
    const mecp2 = { select: 'NM_001110792.2', plusClinical: ['NM_004992.4'] };
    assert.deepEqual(await std(after(500, ['NM_001110792.2']), after(500, ''), after(1, mecp2), 'MECP2'),
        { accession: 'NM_001110792.2', basis: 'the MANE Select transcript of MECP2', plusClinical: ['NM_004992.4'] });
    assert.deepEqual(await std(after(1, ['NM_001110792.2']), after(1, ''), after(10, mecp2), 'MECP2'),
        { accession: 'NM_001110792.2', basis: 'the MANE Select transcript of MECP2', plusClinical: ['NM_004992.4'] },
        'Ensembl first, and NCBI still in time for the Plus Clinical transcript');
    // …a gene with none has an empty list, which is not "unknown"…
    assert.deepEqual((await std(after(1, ['NM_000314.8']), after(1, ''), after(1, { select: 'NM_000314.8', plusClinical: [] }))).plusClinical, []);
    // …and NCBI answering too late leaves it unknown rather than holding the lookup.
    assert.equal((await std(after(1, ['NM_001110792.2']), after(1, ''), after(300, mecp2), 'MECP2')).plusClinical, null);
    // NCBI has no MANE Select for the gene (CDKN2A at the time of writing: its
    // RefSeq record has moved a version on): the other sources name it.
    assert.deepEqual(await std(after(5, ['NM_000077.5']), after(1, ''), after(1, { select: '', plusClinical: ['NM_058195.4'] }), 'CDKN2A'),
        { accession: 'NM_000077.5', basis: 'the MANE Select transcript of CDKN2A', plusClinical: ['NM_058195.4'] });
    // ...and when nobody names a MANE Select transcript at all, RefSeq's own pick.
    assert.deepEqual(await std(after(1, []), after(1, ''), after(1, { select: '', plusClinical: ['NM_058195.4'], refseqSelect: 'NM_000077.6' }), 'CDKN2A'),
        { accession: 'NM_000077.6', basis: 'the RefSeq Select transcript of CDKN2A', plusClinical: ['NM_058195.4'] });
    // MANE Select is never listed among its own alternatives.
    assert.deepEqual((await std(after(1, ['NM_000314.8']), after(1, ''), after(1, { select: 'NM_000314.8', plusClinical: ['NM_000314.7', 'NM_1.1'] }))).plusClinical, ['NM_1.1']);
}

// ── NCBI's list of a gene's MANE transcripts ───────────────────────────────
{
    // (Trimmed from the live answers for MECP2, CDKN2A and "TTN".)
    const report = (symbol, transcripts) => ({ product: { symbol, transcripts: transcripts.map(([accession_version, select_category, ensembl_transcript]) => ({ accession_version, select_category, ensembl_transcript })) } });
    // MANE Select is one transcript with two names; NCBI gives Ensembl's as well.
    const mecp2 = { reports: [report('MECP2', [['NM_001110792.2', 'MANE_SELECT', 'ENST00000453960.7'], ['NM_004992.4', 'MANE_PLUS_CLINICAL', 'ENST00000303391.11'], ['NM_001316337.2'], ['XM_011531166.3']])] };
    assert.deepEqual(h.clinicalTranscriptsFromReport(mecp2, 'mecp2'),
        { select: 'NM_001110792.2', selectEnsembl: 'ENST00000453960.7', plusClinical: ['NM_004992.4'], refseqSelect: '' });
    const cdkn2a = { reports: [report('CDKN2A', [['NM_000077.6', 'REFSEQ_SELECT'], ['NM_058195.4', 'MANE_PLUS_CLINICAL', 'ENST00000579755.2']])] };
    assert.deepEqual(h.clinicalTranscriptsFromReport(cdkn2a, 'CDKN2A'), { select: '', selectEnsembl: '', plusClinical: ['NM_058195.4'], refseqSelect: 'NM_000077.6' },
        'RefSeq Select is not MANE Select, and is kept apart');
    // A symbol search also returns genes that have the symbol as an alias.
    const ttn = { reports: [report('TTR', [['NM_000371.4', 'MANE_SELECT', 'ENST00000237014.8']]), report('TTN', [['NM_001267550.2', 'MANE_SELECT']])] };
    assert.deepEqual(h.clinicalTranscriptsFromReport(ttn, 'TTN'), { select: 'NM_001267550.2', selectEnsembl: '', plusClinical: [], refseqSelect: '' });
    assert.deepEqual(h.clinicalTranscriptsFromReport({}, 'NOSUCHGENE'), { select: '', selectEnsembl: '', plusClinical: [], refseqSelect: '' });
    assert.equal(h.clinicalTranscriptsFromReport({ reports: [report('X1', [['NM_1.1', 'MANE_SELECT', '<img onerror=x>']])] }, 'X1').selectEnsembl, '',
        'only an Ensembl transcript id is taken for one');

    // ── The "(MANE Select)" label: whichever authority names it first ──────
    {
        const after = (ms, value) => () => new Promise(resolve => setTimeout(() => resolve(value), ms));
        const never = () => new Promise(() => {});
        const ncbi = { select: 'NM_004333.6', selectEnsembl: 'ENST00000646891.2', plusClinical: ['NM_001374258.1'], refseqSelect: '' };
        // Ensembl's gene lookup not answering (29 s on the day this was written):
        // NCBI's list names the transcript, under both its names.
        let started = Date.now();
        assert.deepEqual(await h.maneSelectForLabel('BRAF', never, after(5, ncbi)), ['NM_004333.6', 'ENST00000646891.2']);
        assert(Date.now() - started < 500, 'without waiting for Ensembl');
        // Ensembl first: its answer, as before.
        assert.deepEqual(await h.maneSelectForLabel('BRAF', after(5, ['NM_004333.6', 'ENST00000646891.2']), never), ['NM_004333.6', 'ENST00000646891.2']);
        assert.deepEqual(await h.maneSelectForLabel('BRAF', after(1, ['NM_004333.6']), after(40, { ...ncbi, select: 'NM_9.9' })), ['NM_004333.6']);
        // One source having nothing is not the answer while the other may have something.
        assert.deepEqual(await h.maneSelectForLabel('BRAF', after(1, []), after(30, ncbi)), ['NM_004333.6', 'ENST00000646891.2']);
        assert.deepEqual(await h.maneSelectForLabel('BRAF', after(30, ['NM_004333.6']), after(1, null)), ['NM_004333.6'], 'NCBI could not be asked');
        assert.deepEqual(await h.maneSelectForLabel('CDKN2A', after(30, ['NM_000077.5']), after(1, { select: '', selectEnsembl: '', plusClinical: [], refseqSelect: 'NM_000077.6' })),
            ['NM_000077.5'], 'RefSeq Select is not MANE Select');
        // Neither names one: no label, once both have said so.
        assert.deepEqual(await h.maneSelectForLabel('X', after(5, []), after(10, null)), []);
        assert.deepEqual(await h.maneSelectForLabel('X', () => { throw new Error('down'); }, () => Promise.reject(new Error('down'))), []);
    }

    const realFetch = globalThis.fetch;
    const urls = [];
    let status = 200;
    globalThis.fetch = async (address) => {
        urls.push(String(address));
        return new Response(JSON.stringify(mecp2), { status });
    };
    try {
        assert.deepEqual(await h.fetchClinicalTranscripts('mecp2'),
            { select: 'NM_001110792.2', selectEnsembl: 'ENST00000453960.7', plusClinical: ['NM_004992.4'], refseqSelect: '' });
        assert.deepEqual(urls, ['https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/MECP2/taxon/9606/product_report']);
        await h.fetchClinicalTranscripts('MECP2');
        assert.equal(urls.length, 1, 'cached per gene');
        // Not answering is "unknown" — not cached, and not asked again at once.
        status = 503;
        assert.equal(await h.fetchClinicalTranscripts('NF1'), null);
        status = 200;
        assert.equal(await h.fetchClinicalTranscripts('NF1'), null);
        assert.equal(urls.length, 2);
        for (const notAGene of ['chr7', '', 'NM_000546.6', '7']) assert.equal(await h.fetchClinicalTranscripts(notAGene), null);
        assert.equal(urls.length, 2);
    } finally {
        globalThis.fetch = realFetch;
    }
}

// ── The MANE lookup: two small requests ────────────────────────────────────
{
    const realFetch = globalThis.fetch;
    const urls = [];
    let fail = false;
    globalThis.fetch = async (address) => {
        const url = String(address);
        urls.push(url.replace('https://rest.ensembl.org/lookup/', ''));
        if (fail) return new Response('{}', { status: 503 });
        if (url.includes('/lookup/symbol/')) return Response.json({ id: 'ENSG00000171862', canonical_transcript: 'ENST00000371953.8' });
        return Response.json({ id: 'ENST00000371953', MANE: [{ type: 'MANE_Select', refseq_match: 'NM_000314.8', id: 'ENST00000371953', version: 8 }] });
    };
    try {
        assert.deepEqual(await h.fetchManeSelectAccessions('pten'), ['NM_000314.8', 'ENST00000371953.8']);
        assert.deepEqual(urls, ['symbol/homo_sapiens/PTEN?content-type=application/json',
            'id/ENST00000371953?expand=1&mane=1&content-type=application/json'],
            'the gene, then its canonical transcript — never the gene "expanded"');
        await h.fetchManeSelectAccessions('PTEN');
        assert.equal(urls.length, 2, 'cached per gene');
        // A failure is not cached as an answer, and is not retried at once.
        fail = true;
        assert.deepEqual(await h.fetchManeSelectAccessions('TP53'), []);
        fail = false;
        const before = urls.length;
        assert.deepEqual(await h.fetchManeSelectAccessions('TP53'), [], 'not asked again straight away');
        assert.equal(urls.length, before);
        // Not a gene symbol: nothing is asked.
        for (const notAGene of ['chr7', 'CHR17', '', 'NM_000546.6', '7']) assert.deepEqual(await h.fetchManeSelectAccessions(notAGene), []);
        assert.equal(urls.length, before);
    } finally {
        globalThis.fetch = realFetch;
    }
}

// ── A name the services reject ─────────────────────────────────────────────
{
    const ensembl = Object.assign(new Error('Variant recoder request failed (400): {"error":"Unable to parse HGVS notation \'NM_000314.8:c.512T>G\':  : '
        + 'Reference allele extracted from NM_000314.8:89711894-89711894 (A) does not match reference allele given by HGVS notation NM_000314.8:c.512T>G (T)"}'), { status: 400 });
    assert.deepEqual(h.referenceBaseRejection(ensembl), { actual: 'A', stated: 'T' });
    const ncbi = Object.assign(new Error('{"error":{"code":400,"message":"Error while processing \'NM_000314.8:c.512T>G\': The submitted reference allele (T) '
        + 'does not agree with the actual reference sequence (A)."}}'), { status: 400 });
    assert.deepEqual(h.referenceBaseRejection(ncbi), { actual: 'A', stated: 'T' });
    assert.equal(h.referenceBaseRejection(Object.assign(new Error('Service Unavailable'), { status: 503 })), null);
    assert.equal(h.referenceBaseRejection(Object.assign(new Error('Could not get a Transcript object for NM_1.1'), { status: 400 })), null);
    assert.equal(h.referenceBaseRejection(Object.assign(new Error(ensembl.message), { status: 500 })), null, 'only a 400 is an answer about the name');
    assert.equal(h.referenceBaseRejection(null), null);
    // Where the base IS as stated: MyVariant holds c.512T>C on PTEN-L.
    const hits = cassettes.pten512.responses['myvariant.info/v1/query?q=PTEN c.512T>G&size=20'].body.hits;
    assert.deepEqual(h.transcriptsWithReferenceBase(hits, 'c.512T>G', 'NM_000314.8'), ['NM_001304717.5']);
    assert.deepEqual(h.transcriptsWithReferenceBase(hits, 'c.512A>G', 'NM_000314.8'), [], 'a different base is no evidence');
    assert.deepEqual(h.transcriptsWithReferenceBase(hits, 'c.512_513del', 'NM_000314.8'), [], 'only substitutions name a base');
    assert.deepEqual(h.transcriptsWithReferenceBase([], 'c.512T>G'), []);
}

// ── What was typed has to be what is shown ─────────────────────────────────
{
    const rows = [{ transcript: 'NM_000314.8', cDNA: 'c.1026+59T>G' }, { transcript: 'NM_001304717.5', cDNA: 'c.1545+59T>G' }];
    assert.equal(h.codingSubstitutionNotShown('PTEN:c.512T>G', 'c.512t>g', rows), 'c.512T>G');
    assert.equal(h.codingSubstitutionNotShown('NM_000314.8:c.512T>G', 'c.512t>g', rows), 'c.512T>G');
    assert.equal(h.codingSubstitutionNotShown('PTEN:c.1026+59T>G', 'c.1026+59t>g', rows), '', 'it is shown');
    // An indel may come back re-spelled, an Ensembl transcript is not among the
    // RefSeq rows, and a genomic or protein query typed no coding change.
    assert.equal(h.codingSubstitutionNotShown('BRCA1:c.5265_5266insC', 'c.5265_5266insc', [{ cDNA: 'c.5266dup' }]), '');
    assert.equal(h.codingSubstitutionNotShown('ENST00000371953:c.512T>G', 'c.512t>g', rows), '');
    assert.equal(h.codingSubstitutionNotShown('chr10:g.89720934T>G', null, rows), '');
    assert.equal(h.codingSubstitutionNotShown('BRAF:p.Val600Glu', null, rows), '');
    assert.equal(h.codingSubstitutionNotShown('PTEN:c.512T>G', 'c.512t>g', []), '', 'nothing is shown to compare with');
}
// ── Which transcript a name is read on ─────────────────────────────────────
{
    // Records as MyVariant returns them, reduced to what the choice reads.
    const record = (id, names) => ({ _id: id, clinvar: { gene: { symbol: 'MECP2' }, hgvs: { coding: names } }, dbnsfp: { genename: 'MECP2' } });
    const t158m = record('chrX:g.153296806G>A', ['NM_001110792.2:c.509C>T', 'NM_004992.4:c.473C>T', 'NM_001316337.2:c.194C>T']);
    const p251l = record('chrX:g.153296527G>A', ['NM_001110792.2:c.788C>T', 'NM_004992.4:c.752C>T', 'NM_001316337.2:c.473C>T', 'NM_001369391.2:c.473C>T']);
    const select = { accession: 'NM_001110792.2', basis: 'the MANE Select transcript of MECP2' };
    const known = { ...select, plusClinical: ['NM_004992.4'] };
    const choose = (hits, change, standard) => h.chooseNameReading(hits, 'MECP2', change, standard);

    assert.deepEqual(h.transcriptsNamingChange(p251l, 'c.473C>T'), ['NM_001316337.2', 'NM_001369391.2'], 'oldest accession first');
    assert.deepEqual(h.transcriptsNamingChange(t158m, 'c.473C>T'), ['NM_004992.4']);
    assert.deepEqual(h.transcriptsNamingChange(t158m, 'c.474C>T'), []);

    // A record on the Plus Clinical transcript, none on MANE Select: read there,
    // whichever record MyVariant happened to rank first.
    for (const hits of [[t158m, p251l], [p251l, t158m]]) {
        const reading = choose(hits, 'c.473C>T', known);
        assert.equal(reading.accession, 'NM_004992.4');
        assert.equal(reading.basis, 'a MANE Plus Clinical transcript of MECP2');
        assert.equal(reading.hit, t158m);
        // MANE Select is offered; the namesake on a minor isoform is not.
        assert.deepEqual(reading.alternatives, [{ accession: 'NM_001110792.2', kind: 'select', known: false }]);
        assert.deepEqual(reading.untried, []);
    }
    // A record on MANE Select wins, and the other reading is still offered.
    {
        const onSelect = record('chrX:g.1G>A', ['NM_001110792.2:c.473C>T']);
        const reading = choose([t158m, onSelect], 'c.473C>T', known);
        assert.equal(reading.accession, 'NM_001110792.2');
        assert.equal(reading.hit, onSelect);
        assert.deepEqual(reading.alternatives, [{ accession: 'NM_004992.4', kind: 'plus', known: true,
            hitId: 'chrX:g.153296806G>A', inClinvar: true }], 'with the record it rests on, and whether ClinVar has that');
    }
    // One record that is the change on both transcripts (KRAS4A/4B share their
    // numbering up to codon 150): nothing to choose between.
    {
        const both = record('chr12:g.25398284C>T', ['NM_004985.5:c.35G>A', 'NM_033360.4:c.35G>A']);
        const reading = h.chooseNameReading([both], 'KRAS', 'c.35G>A',
            { accession: 'NM_004985.5', basis: 'the MANE Select transcript of KRAS', plusClinical: ['NM_033360.4'] });
        assert.equal(reading.accession, 'NM_004985.5');
        assert.deepEqual(reading.alternatives, []);
    }
    // No record on either: MANE Select, with the Plus Clinical transcript to try
    // if the name does not fit — and offered in any case, since no record
    // decided between them (a deletion names no bases and "fits" either). A
    // namesake elsewhere is now worth a mention too.
    {
        const reading = choose([p251l], 'c.473C>T', known);
        assert.equal(reading.accession, 'NM_001110792.2');
        assert.equal(reading.hit, null);
        assert.deepEqual(reading.untried, ['NM_004992.4']);
        assert.deepEqual(reading.alternatives, [{ accession: 'NM_004992.4', kind: 'plus', known: false },
            { accession: 'NM_001316337.2', kind: 'listed', known: true, hitId: p251l._id, inClinvar: true }]);
        assert.deepEqual(choose([], 'c.473C>T', known).alternatives, [{ accession: 'NM_004992.4', kind: 'plus', known: false }]);
        // (A record that did decide leaves the unrecorded reading unmentioned, as before.)
        assert.deepEqual(choose([t158m], 'c.473C>T', known).alternatives, [{ accession: 'NM_001110792.2', kind: 'select', known: false }]);
    }
    // A gene with no Plus Clinical transcript: MANE Select, as before.
    {
        const reading = choose([p251l], 'c.473C>T', { ...select, plusClinical: [] });
        assert.equal(reading.accession, 'NM_001110792.2');
        assert.deepEqual(reading.untried, []);
    }
    // NCBI did not say (null is not "none"): a record that carries the name on
    // another transcript is all there is to go by — the first-ranked one.
    {
        const reading = choose([t158m, p251l], 'c.473C>T', { ...select, plusClinical: null });
        assert.equal(reading.accession, 'NM_004992.4');
        assert.match(reading.basis, /a transcript ClinVar lists this name on \(whether MECP2 is reported on it could not be learned\)/);
        assert.equal(reading.hit, t158m);
        assert.deepEqual(reading.alternatives, [{ accession: 'NM_001110792.2', kind: 'select', known: false }]);
        // …unless a record carries it on MANE Select itself.
        const onSelect = record('chrX:g.1G>A', ['NM_001110792.2:c.473C>T']);
        assert.equal(choose([t158m, onSelect], 'c.473C>T', { ...select, plusClinical: null }).accession, 'NM_001110792.2');
        assert.equal(choose([], 'c.473C>T', { ...select, plusClinical: null }).accession, 'NM_001110792.2');
    }

    // What the card says.
    const typed = { typed: 'MECP2 c.473C>T', change: 'c.473C>T', select: 'NM_001110792.2' };
    assert.equal(h.nameReadingNote(null), null);
    assert.equal(h.nameReadingNote({ ...typed, accession: 'NM_001110792.2', basis: select.basis, alternatives: [], rejected: [] }), null,
        'MANE Select and nothing else to choose from: nothing to say');
    assert.deepEqual(h.nameReadingNote({ ...typed, accession: 'NM_004992.4', basis: 'a MANE Plus Clinical transcript of MECP2',
        alternatives: [{ accession: 'NM_001110792.2', kind: 'select', known: false }], rejected: [] }),
    { text: 'MECP2 c.473C>T is read on NM_004992.4, a MANE Plus Clinical transcript of MECP2.',
        links: [{ variant: 'NM_001110792.2:c.473C>T', lead: 'If the MANE Select transcript was meant, look up' }] });
    assert.deepEqual(h.nameReadingNote({ ...typed, accession: 'NM_001110792.2', basis: select.basis,
        alternatives: [{ accession: 'NM_004992.4', kind: 'plus', known: true, inClinvar: true }], rejected: [] }).links,
    [{ variant: 'NM_004992.4:c.473C>T', lead: 'Counted on the MANE Plus Clinical transcript the same name is another variant, which ClinVar lists:' }]);
    // ("which ClinVar lists" is said of a record ClinVar has, not of any record.)
    assert.deepEqual(h.nameReadingNote({ ...typed, accession: 'NM_001110792.2', basis: select.basis,
        alternatives: [{ accession: 'NM_004992.4', kind: 'plus', known: true, inClinvar: false }], rejected: [] }).links,
    [{ variant: 'NM_004992.4:c.473C>T', lead: 'Counted on the MANE Plus Clinical transcript the same name is another variant:' }]);
    assert.deepEqual(h.nameReadingNote({ ...typed, accession: 'NM_004992.4', basis: 'a MANE Plus Clinical transcript of MECP2', alternatives: [],
        rejected: [{ accession: 'NM_001110792.2', basis: select.basis, actual: 'G', stated: 'C' }] }),
    { text: 'MECP2 c.473C>T is read on NM_004992.4, a MANE Plus Clinical transcript of MECP2: it is not a change on '
        + 'NM_001110792.2, the MANE Select transcript of MECP2 (G there, not C).', links: [] });

    // ClinVar's own name for a variant, and the transcript it is on.
    const p14 = { clinvar: { rcv: [{ preferred_name: 'NM_058195.4(CDKN2A):c.3G>C (p.Met1?)' }, { preferred_name: 'NM_058195.4(CDKN2A):c.3G>C (p.Met1?)' }] } };
    assert.equal(h.clinvarNamedTranscript(p14, 'c.3G>C'), 'NM_058195.4');
    assert.equal(h.clinvarNamedTranscript(p14, 'c.3G>A'), '');
    assert.equal(h.clinvarNamedTranscript({}, 'c.3G>C'), '');

    // The row shown for an accession is that transcript's, a version behind or not.
    const rows = [
        { transcript: 'NM_001042492.3', cDNA: 'c.6852_6855del', protein: 'p.Tyr2285ThrfsTer5', maneSelect: true },
        { transcript: 'NM_000267.3', cDNA: 'c.6789_6792del', protein: 'p.Tyr2264ThrfsTer5' }
    ];
    assert.equal(h.selectTranscriptForQuery(rows, 'NM_000267.4:c.6789_6792del'), 1, 'the transcript asked for, one version on');
    assert.equal(h.selectTranscriptForQuery(rows, 'NM_000267.3:c.6789_6792del'), 1);
    assert.equal(h.selectTranscriptForQuery(rows, 'NM_000267.4:c.1A>G'), 0, 'not that change on that transcript: the standard row');
    const versions = [{ transcript: 'NM_000267.3', cDNA: 'c.1A>G' }, { transcript: 'NM_000267.4', cDNA: 'c.1A>G' }];
    assert.equal(h.selectTranscriptForQuery(versions, 'NM_000267.4:c.1A>G'), 1, 'the exact version when it is there');
}
console.log('Gene-symbol name unit checks passed');

// 1. ClinVar itself names the variant "BRCA1 c.5266dup" on its standard
//    transcript: settled at once, and Ensembl's recoder — which does not answer
//    this within a minute — is not asked at all.
{
    const page = await runLookup(cassettes.brca1Dup);
    assertResolved(assert, page, 'BRCA1 c.5266dup');
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.match(page.variant, /Transcript: NM_007294\.4 \(MANE Select\)c\.: c\.5266dup/);
    assert.deepEqual(recoderQueries(page), []);
    assert.doesNotMatch(page.panel, /Resolve nomenclature/, 'no step for a service that was not asked');
    // The label comes from the two small lookups, not the 236 KB one.
    assert.deepEqual(asked(page, /\/lookup\//).map(r => r.match(/lookup\/([^?]+\?[^&]*)/)[1]),
        ['symbol/homo_sapiens/BRCA1?content-type=application/json', 'id/ENST00000357654?expand=1']);
    console.log('ok gene-symbol name: BRCA1 c.5266dup, as ClinVar names it, without the recoder');
}

// 2. A name MyVariant has no record for. Was: Ensembl asked "BRCA2:c.1083_1084insCA".
{
    const page = await runLookup(cassettes.brca2GeneInsCA);
    assertResolved(assert, page, 'BRCA2 c.1083_1084insCA');
    assert.match(page.panel, /read as BRCA2:c\.1083_1084insCA · on NM_000059\.4, the MANE Select transcript of BRCA2/);
    assert.deepEqual(recoderQueries(page), ['NM_000059.4:c.1083_1084insCA'], 'asked about the transcript, not the symbol');
    assert.match(page.variant, /g\.: chr13:g\.32906698_32906699insCA/);
    assert.match(page.variant, /Gene: BRCA2Transcript: NM_000059\.4 \(MANE Select\)c\.: c\.1083_1084insCA/);
    assert.match(page.variant, /p\.: p\.Asp362GlnfsTer6/);
    assert(asked(page, /[?&]gene=BRCA2\b/).length > 0, 'the gene cards still know the gene');
    console.log('ok gene-symbol name: read on the MANE Select transcript, and the recoder asked about that');
}

// 3. Neither MANE source answering (Ensembl's record, NCBI's list): ClinVar's
//    naming transcript stands in. With Ensembl alone not answering, NCBI's
//    list names the MANE Select transcript just the same.
{
    const page = await runLookup(cassettes.brca2GeneInsCA, { without: [/\/lookup\//, /\/datasets\//] });
    assertResolved(assert, page, 'BRCA2 c.1083_1084insCA, no MANE lookup');
    assert.match(page.panel, /on NM_000059\.4, the transcript ClinVar names BRCA2 variants on/);
    assert.deepEqual(recoderQueries(page), ['NM_000059.4:c.1083_1084insCA']);
    assert.match(page.variant, /c\.: c\.1083_1084insCA/);
    const ncbiOnly = await runLookup(cassettes.brca2GeneInsCA, { without: [/\/lookup\//] });
    assertResolved(assert, ncbiOnly, 'BRCA2 c.1083_1084insCA, Ensembl\'s MANE lookup down');
    assert.match(ncbiOnly.panel, /on NM_000059\.4, the MANE Select transcript of BRCA2/);
    console.log('ok gene-symbol name: NCBI\'s list, then ClinVar\'s naming transcript, when Ensembl is slow');
}

// 4. PTEN c.512T>G. On the MANE Select transcript the base at c.512 is A, so
//    the name is not a change there. Was: chr10:g.89720934T>G, an intron
//    variant, presented as the answer.
{
    const page = await runLookup(cassettes.pten512);
    assert.equal(page.status, 'Error: PTEN c.512T>G is not a change on NM_000314.8, the MANE Select transcript of PTEN: '
        + 'the reference base at that position is A, not T. The name may be counted on another transcript — enter it with that '
        + 'transcript\'s accession (NM_001304717.5:c.512T>G — the base there is T).');
    assert.equal(page.hidden, true, 'nothing is rendered for a name that does not fit');
    assert.deepEqual(recoderQueries(page), ['NM_000314.8:c.512T>G'], 'Ensembl is never left to pick a transcript');
    assert(!page.requests.some(r => r.includes('89720934')), 'the intron variant is never even looked at');
    assert.match(page.panel, /✕Resolve nomenclature · Ensembl\s*NM_000314\.8 has A there, not T/);
    assert.doesNotMatch(page.status, /retry/i, 'it is the name that is wrong, not the service');
    console.log('ok gene-symbol name: PTEN c.512T>G is refused with the reason and the likely transcript');
}

// 5. The same name typed with the transcript it IS a change on resolves.
//    (Covered end to end in variant-lookup-ui.test.js; here, the suggestion
//    made above is a name the page accepts as an accession query.)
assert.equal(h.geneSymbolCodingName('NM_001304717.5:c.512T>G'), null);

// 6. An accession typed with the wrong base gets the same plain statement.
//    Was: "Lookup could not complete: … failed. Please retry in a moment."
{
    const page = await runLookup(cassettes.pten512, { input: 'NM_000314.8:c.512T>G' });
    assert.equal(page.status, 'Error: NM_000314.8:c.512T>G names a reference base that is not there: NM_000314.8 has A at that '
        + 'position, not T. Check the position, the base and the transcript version.');
    assert.equal(page.hidden, true);
    console.log('ok accession name: a wrong reference base is reported as such');
}

// 7. Neither source can say which transcript is standard. The old behaviour
//    is all that is left — and now it is said, in the panel and beside the
//    coding change that is not the one typed.
{
    const page = await runLookup(cassettes.pten512, { without: [/\/lookup\//, /\/datasets\//, /clinvar\.gene\.symbol/] });
    assertResolved(assert, page, 'PTEN c.512T>G, transcript not settled');
    assert.deepEqual(recoderQueries(page), ['PTEN:c.512T>G']);
    assert.match(page.variant, /g\.: chr10:g\.89720934T>G/);
    assert.match(page.variant, /c\.: c\.1026\+59T>G/);
    assert.match(page.panel, /▲Parse input\s*read as PTEN:c\.512T>G — could not learn which transcript is PTEN's standard one/);
    assert.match(page.variant, /c\.512T>G was entered, but this variant is c\.512T>G on none of the transcripts shown here/);
    console.log('ok gene-symbol name: an unsettled transcript is said, twice');
}

// 8. A substitution ClinVar names on the standard transcript: KRAS G12D, and
//    no warning beside a coding change that is the one typed.
{
    const page = await runLookup(cassettes.kras35);
    assertResolved(assert, page, 'KRAS c.35G>A');
    assert.match(page.variant, /g\.: chr12:g\.25398284C>T/);
    assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.35G>Ap\.: p\.Gly12Asp/);
    assert.equal(page.document.querySelector('.coding-name-note'), null);
    assert.deepEqual(recoderQueries(page), []);
    console.log('ok gene-symbol name: KRAS c.35G>A is G12D');
}

const readingNote = page => page.document.querySelector('.name-reading-note');
// (For "no note": the text, so that a failure prints a sentence and not the page's whole tree.)
const readingText = page => (readingNote(page) ? readingNote(page).textContent : null);
const noteLinks = page => [...readingNote(page).querySelectorAll('a')].map(a => decodeURIComponent(a.getAttribute('href')));

// 9. MECP2 c.473C>T — T158M, the commonest Rett syndrome variant, in the
//    numbering every report uses. It is not a change on MANE Select (the base
//    there is G); MECP2's MANE Plus Clinical transcript is the one it is counted
//    on. Was, read on MANE Select alone: refused, with a hint that led to a
//    benign variant on a third transcript.
{
    const page = await runLookup(cassettes.mecp2T158M);
    assertResolved(assert, page, 'MECP2 c.473C>T');
    assert.match(page.panel, /read as MECP2:c\.473C>T · on NM_004992\.4, a MANE Plus Clinical transcript of MECP2/);
    assert.match(page.variant, /g\.: chrX:g\.153296806G>A/);
    assert.match(page.variant, /Gene: MECP2Transcript: NM_004992\.4c\.: c\.473C>Tp\.: p\.Thr158Met \(p\.T158M\)/);
    assert.match(page.clinvar, /Variation ID: 11811/);
    assert.equal(readingNote(page).textContent, 'MECP2 c.473C>T is read on NM_004992.4, a MANE Plus Clinical transcript of MECP2. '
        + 'If the MANE Select transcript was meant, look up NM_001110792.2:c.473C>T.');
    assert.deepEqual(noteLinks(page), ['?variant=NM_001110792.2:c.473C>T']);
    assert.deepEqual(recoderQueries(page), [], 'MyVariant holds it: Ensembl is not asked');
    assert.equal(asked(page, /\/datasets\/v2\/gene\/symbol\/MECP2\//).length, 1, 'NCBI\'s list is asked once');
    console.log('ok gene-symbol name: MECP2 c.473C>T is T158M, on the Plus Clinical transcript');
}

// 10. NF1 c.6789_6792del. A deletion names no bases, so it "fits" anywhere: on
//     MANE Select it is four bases nobody has reported, and that is what was
//     shown — "Clinical significance: N/A" — for ClinVar 216866, Pathogenic.
{
    const page = await runLookup(cassettes.nf1LegacyDel);
    assertResolved(assert, page, 'NF1 c.6789_6792del');
    assert.match(page.variant, /g\.: chr17:g\.29665752_29665755del/);
    assert.doesNotMatch(page.variant, /29665127/, 'not the MANE Select reading');
    assert.match(page.clinvar, /Variation ID: 216866/);
    // The row shown is the transcript the name was read on (GRCh37's annotation
    // is a version behind NCBI's list), under the name that was typed.
    assert.match(page.variant, /Transcript: NM_000267\.3c\.: c\.6789_6792delp\.: p\.Tyr2264ThrfsTer5/);
    assert.equal(readingNote(page).textContent, 'NF1 c.6789_6792del is read on NM_000267.4, a MANE Plus Clinical transcript of NF1. '
        + 'If the MANE Select transcript was meant, look up NM_001042492.3:c.6789_6792del.');
    console.log('ok gene-symbol name: NF1 c.6789_6792del is the known deletion, not four other bases');
}

// 11. A name no record carries: MECP2 c.474G>C. MANE Select is tried first and
//     says the base there is C; the Plus Clinical transcript is then read, and
//     the card says both.
{
    const page = await runLookup(cassettes.mecp2OnPlusClinical);
    assertResolved(assert, page, 'MECP2 c.474G>C');
    assert.deepEqual(recoderQueries(page), ['NM_001110792.2:c.474G>C', 'NM_004992.4:c.474G>C']);
    assert.match(page.panel, /read as NM_004992\.4:c\.474G>C · not a change on NM_001110792\.2 \(C there, not G\) · on NM_004992\.4, a MANE Plus Clinical transcript of MECP2/);
    assert.match(page.variant, /g\.: chrX:g\.153296805C>G/);
    assert.match(page.variant, /Transcript: NM_004992\.4c\.: c\.474G>Cp\.: p\.Thr158=/);
    assert.equal(readingNote(page).textContent, 'MECP2 c.474G>C is read on NM_004992.4, a MANE Plus Clinical transcript of MECP2: '
        + 'it is not a change on NM_001110792.2, the MANE Select transcript of MECP2 (C there, not G).');
    assert(asked(page, /[?&]gene=MECP2\b/).length > 0, 'the gene cards still know the gene');
    console.log('ok gene-symbol name: a name that does not fit MANE Select is read on the Plus Clinical transcript');
}

// 12. CDKN2A c.3G>C is ClinVar's own name for a variant — on NM_058195.4, the
//     p14ARF transcript. On NM_000077 the same name is p16's start codon. The
//     variant shown is as before; that there was a choice is now said.
{
    const page = await runLookup(cassettes.cdkn2aC3);
    assertResolved(assert, page, 'CDKN2A c.3G>C');
    assert.match(page.variant, /g\.: chr9:g\.21994328C>G/);
    assert.match(page.clinvar, /Variation ID: 1018474/);
    // (NM_000077.6 is RefSeq's pick, which stands in here because Ensembl's MANE
    // record did not answer in this recording and NCBI lists no MANE Select for
    // CDKN2A. It used to be introduced as "the MANE Select transcript".)
    assert.equal(readingNote(page).textContent, 'CDKN2A c.3G>C is read on NM_058195.4, the transcript ClinVar names this variant on '
        + '(a MANE Plus Clinical transcript of CDKN2A). If the RefSeq Select transcript was meant, look up NM_000077.6:c.3G>C.');
    assert.deepEqual(noteLinks(page), ['?variant=NM_000077.6:c.3G>C']);
    // A name ClinVar gives on MANE Select itself has nothing to add (scenario 8).
    assert.equal(readingText(await runLookup(cassettes.kras35)), null);
    assert.equal(readingText(await runLookup(cassettes.brca1Dup)), null);
    console.log('ok gene-symbol name: a name ClinVar gives on another transcript says so');
}

// 13. MUTYH c.536A>G is Y179C on the transcript MUTYH is reported on
//     (NM_001128425). Before any of this it was shown as p.Glu179Gly, ClinVar
//     533304 — the same name on the lowest-numbered transcript that had one.
{
    const page = await runLookup(cassettes.mutyhY179C);
    assertResolved(assert, page, 'MUTYH c.536A>G');
    assert.match(page.variant, /g\.: chr1:g\.45798475T>C/);
    assert.match(page.variant, /Transcript: NM_001128425\.2c\.: c\.536A>Gp\.: p\.Tyr179Cys \(p\.Y179C\)/);
    assert.match(page.clinvar, /Variation ID: 5293/);
    assert.match(readingNote(page).textContent, /^MUTYH c\.536A>G is read on NM_001128425\.2, a MANE Plus Clinical transcript of MUTYH\. /);
    console.log('ok gene-symbol name: MUTYH c.536A>G is Y179C');
}

// 14. NCBI's list not answering. Whether the gene has a second clinical
//     transcript is then unknown — which is not "it has none". A record that
//     carries the name is what there is to go by (as before this rule), and
//     the card says how little that is.
{
    const page = await runLookup(cassettes.mecp2T158M, { without: [/\/datasets\//] });
    assertResolved(assert, page, 'MECP2 c.473C>T, NCBI not answering');
    assert.match(page.clinvar, /Variation ID: 11811/);
    assert.equal(readingNote(page).textContent, 'MECP2 c.473C>T is read on NM_004992.4, a transcript ClinVar lists this name on '
        + '(whether MECP2 is reported on it could not be learned). If the MANE Select transcript was meant, look up NM_001110792.2:c.473C>T.');
    const del = await runLookup(cassettes.nf1LegacyDel, { without: [/\/datasets\//] });
    assert.match(del.clinvar, /Variation ID: 216866/);
    assert.match(readingNote(del).textContent, /^NF1 c\.6789_6792del is read on NM_000267\.3, a transcript ClinVar lists this name on/);
    // With no record to go by, the name is refused on MANE Select as before —
    // and the transcript suggested is the one the gene was long counted on.
    const unrecorded = await runLookup(cassettes.mecp2OnPlusClinical, { without: [/\/datasets\//] });
    assert.equal(unrecorded.status, 'Error: MECP2 c.474G>C is not a change on NM_001110792.2, the MANE Select transcript of MECP2: '
        + 'the reference base at that position is C, not G. The name may be counted on another transcript — enter it with that '
        + 'transcript\'s accession (NM_004992.3:c.474G>C — the base there is G).');
    console.log('ok gene-symbol name: without NCBI\'s list, the record that carries the name — and it says so');
}

// 15. MyVariant's name search failing. The name still has to be read on one of
//     the gene's own transcripts: it used to go to Ensembl as typed, which came
//     back with the intron variant of scenario 4, "Variant resolved".
{
    const page = await runLookup(cassettes.pten512, {
        route: (url, key, respond) => (/myvariant\.info\/v1\/query/.test(key) && /c\.512T>G/.test(key)
            ? respond({ error: 'Service Unavailable' }, 503) : undefined)
    });
    assert.deepEqual(recoderQueries(page), ['NM_000314.8:c.512T>G'], 'never "PTEN:c.512T>G"');
    assert.match(page.panel, /read as PTEN:c\.512T>G · on NM_000314\.8, the MANE Select transcript of PTEN/);
    // MyVariant being down is what the lookup then reports; what it no longer
    // does is resolve to a variant nobody named.
    assert.match(page.status, /^Error: MyVariant\.info query request failed \(503\)/);
    assert.equal(page.hidden, true);
    assert(!page.requests.some(r => r.includes('89720934')));
    console.log('ok gene-symbol name: a failed MyVariant search does not hand the name to Ensembl as typed');
}

// 16. A transcript version newer than GRCh37's annotation. RefSeq moves on and
//     that annotation does not: "NM_000077.6:c.3G>C" — CDKN2A's current record,
//     and what the note in scenario 12 links to — came back "Lookup could not
//     complete … failed. Please retry in a moment." Ensembl reads the accession
//     without its version as the one it has; that is looked up, and said.
{
    const unknown = 'Unable to parse HGVS notation \'NM_000059.9:c.1083_1084insCA\': Could not get a Transcript object for \'NM_000059.9\'';
    assert.equal(h.unknownTranscriptVersion(Object.assign(new Error(`Variant recoder request failed (400): {"error":"${unknown}"}`), { status: 400 })), 'NM_000059.9');
    assert.equal(h.unknownTranscriptVersion(Object.assign(new Error('Could not get a Transcript object for \'NM_000059.9\''), { status: 503 })), '', 'only a rejection counts');
    assert.equal(h.unknownTranscriptVersion(Object.assign(new Error('Unable to parse HGVS notation'), { status: 400 })), '');
    assert.equal(h.unknownTranscriptVersion(null), '');

    const cassette = cassettes.brca2GeneInsCA;
    const known = cassette.responses['grch37.rest.ensembl.org/variant_recoder/human/NM_000059.4:c.1083_1084insCA?content-type=application/json'];
    assert.equal(known.status, 200, 'the recording holds Ensembl\'s answer for the version it has');
    const route = (url, key, respond) => {
        if (/variant_recoder\/human\/NM_000059\.9:/.test(key)) return respond({ error: unknown }, 400);
        if (/api\.ncbi\.nlm\.nih\.gov\/variation/.test(key)) return respond({ error: 'no GRCh37 placement' }, 404);
        if (/\/vep\/human\/hgvs\/NM_000059/.test(key)) return respond({ error: unknown }, 400);
        if (/variant_recoder\/human\/NM_000059:/.test(key)) return respond(known.body);
        return undefined;
    };
    const page = await runLookup(cassette, { input: 'NM_000059.9:c.1083_1084insCA', route });
    assertResolved(assert, page, 'NM_000059.9:c.1083_1084insCA');
    assert.deepEqual(recoderQueries(page).filter((q, i, all) => all.indexOf(q) === i),
        ['NM_000059.9:c.1083_1084insCA', 'NM_000059:c.1083_1084insCA']);
    assert.match(page.panel, /read as NM_000059:c\.1083_1084insCA · NM_000059\.9 is a newer version than GRCh37's annotation has — read on the version it does have/);
    assert.match(page.variant, /g\.: chr13:g\.32906698_32906699insCA/);
    assert.match(page.variant, /Transcript: NM_000059\.4 \(MANE Select\)c\.: c\.1083_1084insCA/);
    assert.match(page.panel, /You entered\s*NM_000059\.9:c\.1083_1084insCA/);
    // If that does not resolve either, the lookup fails once — no third pass.
    const never = await runLookup(cassette, {
        input: 'NM_000059.9:c.1083_1084insCA',
        route: (url, key, respond) => (/variant_recoder|\/vep\/|api\.ncbi\.nlm\.nih\.gov\/variation/.test(key)
            ? respond({ error: unknown.replace(/NM_000059\.9'$/, 'NM_000059.9\'') }, 400) : undefined)
    });
    assert.match(never.status, /^Error: /);
    assert.equal(recoderQueries(never).filter(q => /^NM_000059:/.test(q)).length > 0, true, 'the version-less name was tried');
    assert.equal(never.hidden, true);
    console.log('ok accession name: a version newer than GRCh37\'s annotation is read on the version it has');
}

// 17. The "(MANE Select)" label when Ensembl's gene lookup does not answer. It
//     took 29 s on the day this was written; every result is held up to 4 s
//     for the label, so each lookup sat those out and then came without it.
//     NCBI's list of the gene's MANE transcripts names the same transcript.
{
    const cassette = cassettes.kras35;
    const labelled = /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.35G>A/;
    // Ensembl's lookup failing outright: the label comes from NCBI's list.
    const failing = await runLookup(cassette, { without: [/rest\.ensembl\.org\/lookup/] });
    assertResolved(assert, failing, 'KRAS c.35G>A, Ensembl lookup down');
    assert.match(failing.variant, labelled);
    // Ensembl's lookup not answering at all, with the page's waits at their
    // real length: the result is not held for it.
    let release;
    const held = new Promise((resolve) => { release = resolve; });
    const started = Date.now();
    const pending = runLookup(cassette, {
        timerCap: 6000,
        route: async (url, key, respond) => {
            if (!/rest\.ensembl\.org\/lookup/.test(key)) return undefined;
            await held;
            return respond({ error: 'Service Unavailable' }, 503);
        }
    });
    const page = await Promise.race([pending, new Promise(resolve => setTimeout(() => resolve(null), 3000))]);
    const took = Date.now() - started;
    release();
    assert(page, 'the lookup finished while Ensembl\'s gene lookup had still not answered');
    assert(took < 2500, `without sitting out the 4 s the label is allowed (${took} ms)`);
    assertResolved(assert, page, 'KRAS c.35G>A, Ensembl lookup not answering');
    assert.match(page.variant, labelled);
    // And either source alone is enough; with neither there is no label.
    assert.match((await runLookup(cassette, { without: [/datasets/] })).variant, labelled);
    const neither = await runLookup(cassette, { without: [/rest\.ensembl\.org\/lookup/, /datasets/] });
    assertResolved(assert, neither, 'KRAS c.35G>A, no MANE source');
    assert.match(neither.variant, /Transcript: NM_004985\.5c\.: c\.35G>A/);
    console.log('ok MANE label: from whichever of Ensembl and NCBI names the transcript first');
}

// 18. A lookup that has been overtaken must not run itself again. "MECP2
//     c.474G>C" is not a change on MANE Select, so the lookup runs again on the
//     Plus Clinical transcript — and did so after "BRAF V600E" had been
//     submitted, while its last fallback was still out: the page then showed
//     MECP2's result under a box that said BRAF V600E, and BRAF's own lookup
//     never appeared. (Found in review.)
{
    const first = cassettes.mecp2OnPlusClinical;
    const second = cassettes.brafV600E;
    const both = { input: first.input, reference: [...first.reference, ...second.reference],
        api: { ...(first.api || {}), ...(second.api || {}) }, responses: { ...first.responses, ...second.responses } };
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    let submitted = false;
    const page = await runLookup(both, {
        settleMs: 400,
        route: async (url, key, respond, options, { document }) => {
            // The first lookup's last fallback before it judges the name: the
            // user submits the next variant while that request is out.
            if (!submitted && key.includes('/vep/human/hgvs/NM_001110792.2:c.474G>C') && !key.includes('refseq=1')) {
                submitted = true;
                document.getElementById('variantInput').value = 'BRAF V600E';
                document.getElementById('variantForm').dispatchEvent(new document.defaultView.Event('submit', { cancelable: true }));
                await wait(2);
                return undefined;
            }
            // The second lookup's recoder takes a little while, as Ensembl's does.
            if (key.includes('/variant_recoder/human/BRAF:p.Val600Glu')) { await wait(60); return undefined; }
            return undefined;
        }
    });
    assert(submitted, 'the second lookup was submitted during the first');
    assertResolved(assert, page, 'BRAF V600E, submitted during MECP2 c.474G>C');
    assert.match(page.panel, /You entered\s*BRAF V600E/);
    assert.match(page.variant, /g\.: chr7:g\.140453136A>T/);
    assert.doesNotMatch(page.panel + page.variant, /MECP2|153296805/, 'nothing of the lookup that was overtaken');
    assert(!page.requests.some(r => /variant_recoder\/human\/NM_004992\.4:c\.474G>C/.test(r)), 'which did not go on to its second pass');
    console.log('ok gene-symbol name: a lookup that has been overtaken does not run itself again');
}

// 18b. ...nor does it write onto the page of the lookup that overtook it. The
//      same moment, the next variant being one that asks VEP nothing ("KRAS
//      p.Gly12Xaa" is a codon) and is still under way when the overtaken
//      lookup's VEP request fails: that failure belongs to no step of the
//      lookup now on the page.
{
    const first = cassettes.mecp2OnPlusClinical;
    const second = cassettes.krasG12Stop;
    const both = { input: first.input, reference: [...first.reference, ...second.reference],
        api: { ...(first.api || {}), ...(second.api || {}) }, responses: { ...first.responses, ...second.responses } };
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    let submitted = false;
    let failed = false;
    const page = await runLookup(both, {
        settleMs: 400,
        route: async (url, key, respond, options, { document }) => {
            if (!submitted && key.includes('/vep/human/hgvs/NM_001110792.2:c.474G>C') && !key.includes('refseq=1')) {
                submitted = true;
                document.getElementById('variantInput').value = 'KRAS p.Gly12Xaa';
                document.getElementById('variantForm').dispatchEvent(new document.defaultView.Event('submit', { cancelable: true }));
                await wait(30);
                failed = true;
                return undefined;
            }
            // The second lookup is held at its first Ensembl request until the
            // first one's VEP request has failed, and a little longer.
            if (submitted && key.includes('/variant_recoder/human/KRAS:p.Gly12Xaa')) {
                for (let i = 0; i < 200 && !failed; i++) await wait(2);
                await wait(20);
                return undefined;
            }
            return undefined;
        }
    });
    assert(submitted && failed, 'the second lookup was under way when the first one\'s VEP request failed');
    assertResolved(assert, page, 'KRAS p.Gly12Xaa, submitted during MECP2 c.474G>C');
    assert.match(page.panel, /You entered\s*KRAS p\.Gly12Xaa/);
    assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.doesNotMatch(page.panel, /Consequences · Ensembl VEP/, 'no VEP step: this lookup asked VEP nothing');
    assert.doesNotMatch(page.panel + page.variant, /MECP2|NM_001110792/, 'nothing of the lookup that was overtaken');
    console.log('ok gene-symbol name: an overtaken lookup\'s failure is not written onto the next one');
}

// 18c. And overtaken earlier still — while it waits for Ensembl's answer about
//      the name — it stops when that answer arrives, instead of asking again
//      and going on to look up the alleles in it. ("KRAS G12*" comes back as
//      three; "BRAF V600E" is on the page by then.)
{
    const first = cassettes.krasG12Stop;
    const second = cassettes.brafV600E;
    const both = { input: 'KRAS G12*', reference: [...first.reference, ...second.reference],
        api: { ...(first.api || {}), ...(second.api || {}) }, responses: { ...first.responses, ...second.responses } };
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    let submitted = false;
    const page = await runLookup(both, {
        settleMs: 400,
        // (Ensembl's answer being what the first lookup waits for: MyVariant's
        // list of the residue's single-base changes, which would cut that wait
        // short, is left out of this one.)
        without: ['dbnsfp.aa.pos'],
        route: async (url, key, respond, options, { document }) => {
            if (!submitted && key.includes('/variant_recoder/human/KRAS:p.Gly12Ter')) {
                // Ensembl takes its time. The first lookup gets as far as
                // waiting for it; only then is the next variant submitted.
                await wait(40);
                submitted = true;
                document.getElementById('variantInput').value = 'BRAF V600E';
                document.getElementById('variantForm').dispatchEvent(new document.defaultView.Event('submit', { cancelable: true }));
                await wait(10);
                return undefined;
            }
            return undefined;
        }
    });
    assert(submitted, 'the second lookup was submitted while the first waited for Ensembl');
    assertResolved(assert, page, 'BRAF V600E, submitted during KRAS G12*');
    assert.match(page.panel, /You entered\s*BRAF V600E/);
    assert.match(page.variant, /g\.: chr7:g\.140453136A>T/);
    assert.equal(asked(page, /variant_recoder\/human\/KRAS:p\.Gly12Ter/).length, 1, 'Ensembl is not asked about the first name again');
    assert.deepEqual(asked(page, /myvariant\.info\/v1\/variant\/chr12:/), [], 'and the alleles in its answer are not looked up');
    assert.deepEqual(page.api('clinvar-protein').filter(r => /gene=KRAS/.test(r)), [], 'nor ClinVar asked about it');
    assert.doesNotMatch(page.panel + page.variant, /KRAS|Gly12/, 'nothing of the lookup that was overtaken');
    console.log('ok gene-symbol name: a lookup overtaken while Ensembl worked stops when Ensembl answers');
}

// 19. A name ClinVar itself settles can still be two variants. "MECP2
//     c.916C>T" is ClinVar's name for p.Arg306Ter on MANE Select — and, on the
//     Plus Clinical transcript the clinic counts on, R306C (ClinVar 11824).
//     The first was shown with nothing to say there was a second. (Found in
//     review; as recorded live.)
{
    // The pieces.
    const record = (id, coding) => ({ _id: id, clinvar: { hgvs: { coding } } });
    const stop = record('chrX:g.153296399G>A', ['NM_001110792.2:c.916C>T']);
    const r306c = record('chrX:g.153296363G>A', ['NM_004992.4:c.916C>T']);
    const minor = record('chrX:g.153296084G>A', ['NM_001316337.2:c.916C>T']);
    const settled = { accession: 'NM_001110792.2', change: 'c.916C>T', hitId: stop._id };
    assert.deepEqual(h.settledNameAlternatives([minor, stop, r306c], settled, 'NM_001110792.2', ['NM_004992.4']),
        [{ accession: 'NM_004992.4', kind: 'plus', known: true, hitId: r306c._id, inClinvar: true }],
        'the namesake on the Plus Clinical transcript, and not the one on a minor isoform');
    assert.deepEqual(h.settledNameAlternatives([stop, { _id: r306c._id, snpeff: { ann: [{ feature_id: 'NM_004992.3', hgvs_c: 'c.916C>T', genename: 'MECP2' }] } }],
        settled, 'NM_001110792.2', ['NM_004992.4']),
    [{ accession: 'NM_004992.4', kind: 'plus', known: true, hitId: r306c._id, inClinvar: false }], 'a record ClinVar does not have is not said to be ClinVar\'s');
    assert.deepEqual(h.settledNameAlternatives([stop], settled, 'NM_001110792.2', ['NM_004992.4']), [], 'no record there: nothing to add');
    assert.deepEqual(h.settledNameAlternatives([stop, { ...stop }], settled, 'NM_001110792.2', ['NM_004992.4']), []);
    // The record that settled the name is not its own namesake: KRAS c.35G>A is
    // one record, and carries the name on both of KRAS's clinical transcripts.
    {
        const g12d = record('chr12:g.25398284C>T', ['NM_004985.5:c.35G>A', 'NM_033360.4:c.35G>A']);
        assert.deepEqual(h.settledNameAlternatives([g12d], { accession: 'NM_004985.5', change: 'c.35G>A', hitId: g12d._id },
            'NM_004985.5', ['NM_033360.4']), []);
    }
    // ClinVar's name on the Plus Clinical transcript: MANE Select is offered, known or not.
    const onPlus = { accession: 'NM_004992.4', change: 'c.916C>T', hitId: r306c._id };
    assert.deepEqual(h.settledNameAlternatives([r306c], onPlus, 'NM_001110792.2', ['NM_004992.4']),
        [{ accession: 'NM_001110792.2', kind: 'select', known: false }]);
    assert.deepEqual(h.settledNameAlternatives([r306c, stop], onPlus, 'NM_001110792.3', ['NM_004992.4']),
        [{ accession: 'NM_001110792.3', kind: 'select', known: true, hitId: stop._id, inClinvar: true }], 'in whichever version');
    // How each is put.
    const said = reading => { const n = h.nameReadingNote(reading); return n && [n.text, ...n.links.map(l => `${l.lead} ${l.variant}`)].join(' '); };
    const base = { typed: 'MECP2 c.916C>T', gene: 'MECP2', change: 'c.916C>T', accession: 'NM_001110792.2',
        basis: 'the MANE Select transcript of MECP2', select: 'NM_001110792.2', untried: [], rejected: [] };
    assert.equal(said({ ...base, alternatives: [] }), null);
    assert.equal(said({ ...base, alternatives: [{ accession: 'NM_004992.4', kind: 'plus', known: true, inClinvar: true }] }),
        'MECP2 c.916C>T is read on NM_001110792.2, the MANE Select transcript of MECP2. Counted on the MANE Plus Clinical transcript '
        + 'the same name is another variant, which ClinVar lists: NM_004992.4:c.916C>T');
    assert.equal(said({ ...base, alternatives: [{ accession: 'NM_004992.4', kind: 'plus', known: false }] }),
        'MECP2 c.916C>T is read on NM_001110792.2, the MANE Select transcript of MECP2. If the MANE Plus Clinical transcript was meant, look up NM_004992.4:c.916C>T');
    assert.equal(said({ ...base, searchFailed: true, alternatives: [{ accession: 'NM_004992.4', kind: 'plus', known: false }] }),
        'MECP2 c.916C>T is read on NM_001110792.2, the MANE Select transcript of MECP2. MyVariant could not be searched for the name, '
        + 'so no record decided this. If the MANE Plus Clinical transcript was meant, look up NM_004992.4:c.916C>T');
    assert.equal(said({ ...base, accession: 'NM_004992.4', basis: 'a MANE Plus Clinical transcript of MECP2',
        alternatives: [{ accession: 'NM_001110792.2', kind: 'select', known: true, inClinvar: true }] }),
    'MECP2 c.916C>T is read on NM_004992.4, a MANE Plus Clinical transcript of MECP2. Counted on the MANE Select transcript '
        + 'the same name is another variant, which ClinVar lists: NM_001110792.2:c.916C>T');
    assert.equal(said({ ...base, accession: 'NM_058195.4', basis: 'x', select: 'NM_000077.6', selectLabel: 'RefSeq Select',
        alternatives: [{ accession: 'NM_000077.6', kind: 'select', known: false }] }),
    'MECP2 c.916C>T is read on NM_058195.4, x. If the RefSeq Select transcript was meant, look up NM_000077.6:c.916C>T');
    // (Where it is ClinVar's naming transcript standing in for MANE Select, it is called that.)
    assert.equal(said({ ...base, accession: 'NM_004992.4', basis: 'x', selectLabel: '',
        alternatives: [{ accession: 'NM_001110792.2', kind: 'select', known: false }] }),
    'MECP2 c.916C>T is read on NM_004992.4, x. If the transcript ClinVar names MECP2 variants on was meant, look up NM_001110792.2:c.916C>T');

    // Which of those readings are another variant at all, once the variant's
    // name on each transcript is known. KRAS's two clinical transcripts number
    // exon 2 alike, NF1's everything before the exon only one of them has: the
    // name is then the same variant whichever was meant, and nothing is offered.
    {
        const row = (transcript, cDNA) => ({ transcript, cDNA });
        const kras = { accession: 'NM_004985.5', change: 'c.34_36dup', hitId: '',
            alternatives: [{ accession: 'NM_033360.4', kind: 'plus', known: false }] };
        assert.deepEqual(h.namesakesGivenRows(kras, [row('NM_004985.5', 'c.36_38dup'), row('NM_033360.4', 'c.36_38dup')]), [],
            'the same name on both: no other reading');
        assert.deepEqual(h.namesakesGivenRows(kras, [row('NM_004985.4', 'c.36_38dupTGG'), row('NM_033360.3', 'c.36_38dup')]), [],
            'whatever the version, and with the bases written out');
        // MECP2's two count 36 bases apart: the name on the other one is elsewhere.
        const mecp2 = { accession: 'NM_001110792.2', change: 'c.916C>T', hitId: 'chrX:g.153296399G>A',
            alternatives: [{ accession: 'NM_004992.4', kind: 'plus', known: true, hitId: 'chrX:g.153296363G>A', inClinvar: true }] };
        assert.deepEqual(h.namesakesGivenRows(mecp2, [row('NM_001110792.2', 'c.916C>T'), row('NM_004992.4', 'c.880C>T')]), mecp2.alternatives);
        // No row for one of them: two substitutions under two ids are two variants...
        assert.deepEqual(h.namesakesGivenRows(mecp2, [row('NM_001110792.2', 'c.916C>T')]), mecp2.alternatives);
        assert.deepEqual(h.namesakesGivenRows(mecp2, []), mecp2.alternatives);
        // ...two records of an indel need not be, and the claim is then not made.
        const indel = { accession: 'NM_004985.5', change: 'c.34_36dup', hitId: 'chr12:g.25398284_25398286dup',
            alternatives: [{ accession: 'NM_033360.4', kind: 'plus', known: true, hitId: 'chr12:g.25398285_25398286insACC', inClinvar: true }] };
        assert.deepEqual(h.namesakesGivenRows(indel, []), [{ ...indel.alternatives[0], known: false }]);
        assert.deepEqual(h.namesakesGivenRows(indel, [row('NM_004985.5', 'c.36_38dup'), row('NM_033360.4', 'c.36_38dup')]), []);
        // An unrecorded reading is left as it is.
        assert.deepEqual(h.namesakesGivenRows({ ...kras }, []), kras.alternatives);
        // A namesake on some further transcript is worth a mention only while
        // the name is not a known variant where it was read.
        const brca1 = { accession: 'NM_007294.4', change: 'c.5265_5266insC', hitId: '',
            alternatives: [{ accession: 'NM_001407610.1', kind: 'listed', known: true, hitId: 'chr17:g.41209139_41209140insG', inClinvar: true }] };
        assert.deepEqual(h.namesakesGivenRows(brca1, [row('NM_007294.4', 'c.5266dup')], { onRecord: true }), []);
        assert.deepEqual(h.namesakesGivenRows(brca1, [row('NM_007294.4', 'c.5266dup')]), [{ ...brca1.alternatives[0], known: false }]);
        // (The gene's other clinical transcript is another matter: that reading stays.)
        assert.deepEqual(h.namesakesGivenRows(mecp2, [row('NM_001110792.2', 'c.916C>T'), row('NM_004992.4', 'c.880C>T')], { onRecord: true }), mecp2.alternatives);
        assert.deepEqual(h.namesakesGivenRows(null, []), []);
    }

    const page = await runLookup(cassettes.mecp2R306);
    assertResolved(assert, page, 'MECP2 c.916C>T');
    assert.match(page.variant, /g\.: chrX:g\.153296399G>A/);
    assert.match(page.variant, /Transcript: NM_001110792\.2 \(MANE Select\)c\.: c\.916C>Tp\.: p\.Arg306Ter/);
    assert.match(page.clinvar, /Variation ID: 11819/);
    assert.equal(readingNote(page).textContent, 'MECP2 c.916C>T is read on NM_001110792.2, the MANE Select transcript of MECP2. '
        + 'Counted on the MANE Plus Clinical transcript the same name is another variant, which ClinVar lists: NM_004992.4:c.916C>T.');
    assert.deepEqual(noteLinks(page), ['?variant=NM_004992.4:c.916C>T']);
    console.log('ok gene-symbol name: a name ClinVar settles on MANE Select still names its namesake on the Plus Clinical transcript');
}

// 20. No record to decide the reading. "NF1 c.6789_6792del" is counted on
//     NM_000267 (ClinVar 216866); a deletion names no bases, so on MANE Select
//     it "fits" as well — four other bases. With MyVariant's name search not
//     answering that is what was shown, "Variant resolved", without a word
//     (the README's own "before" row). It is still what is read; the card now
//     says no record decided it, and offers the other transcript. (Found in
//     review; the MANE Select answers as recorded live.)
{
    const cassette = cassettes.nf1LegacyDel;
    const searchDown = (url, key, respond) => (key.startsWith('myvariant.info/v1/query?q=NF1') ? respond({ error: 'Service Unavailable' }, 503) : undefined);
    const down = await runLookup(cassette, { route: searchDown });
    assertResolved(assert, down, 'NF1 c.6789_6792del, name search down');
    assert.match(down.variant, /g\.: chr17:g\.29665127_29665130del/, 'read on MANE Select, there being nothing to say otherwise');
    assert.equal(readingNote(down).textContent, 'NF1 c.6789_6792del is read on NM_001042492.3, the MANE Select transcript of NF1. '
        + 'MyVariant could not be searched for the name, so no record decided this. '
        + 'If the MANE Plus Clinical transcript was meant, look up NM_000267.4:c.6789_6792del.');
    assert.deepEqual(noteLinks(down), ['?variant=NM_000267.4:c.6789_6792del']);
    assert.match(down.panel, /Annotations · MyVariant\.info/, 'and the panel has a MyVariant step');
    // The search answering that it holds no record under the name: the same
    // reading and the same offer, without the sentence about MyVariant.
    const none = await runLookup(cassette, { route: (url, key, respond) => (key.startsWith('myvariant.info/v1/query?q=NF1') ? respond({ total: 0, hits: [] }) : undefined) });
    assertResolved(assert, none, 'NF1 c.6789_6792del, no record under the name');
    assert.match(none.variant, /g\.: chr17:g\.29665127_29665130del/);
    assert.equal(readingNote(none).textContent, 'NF1 c.6789_6792del is read on NM_001042492.3, the MANE Select transcript of NF1. '
        + 'If the MANE Plus Clinical transcript was meant, look up NM_000267.4:c.6789_6792del.');
    // With the search answering as recorded, the record decides — as before.
    const found = await runLookup(cassette);
    assert.match(found.variant, /g\.: chr17:g\.29665752_29665755del/);
    assert.match(found.clinvar, /Variation ID: 216866/);
    console.log('ok gene-symbol name: with no record to decide, the other clinical transcript is offered');
}

// 21. "Read on the version it does have" is said of the transcript it happened
//     to. Here MANE Select's version is newer than GRCh37's annotation, the
//     name does not fit that transcript in any version, and it is read on the
//     Plus Clinical one: the final pass must not say its version was dropped.
{
    const cassette = cassettes.mecp2OnPlusClinical;
    const unknown = 'Unable to parse HGVS notation \'NM_001110792.2:c.474G>C\':  : Could not get a Transcript object for \'NM_001110792.2\'';
    const recordedFor = (key) => cassette.responses[key.replace('NM_001110792:', 'NM_001110792.2:')];
    const page = await runLookup(cassette, {
        route: (url, key, respond) => {
            // The versioned accession is unknown to GRCh37's annotation and to NCBI's placement...
            if (/NM_001110792\.2:c\.474G>C/.test(key) && /variant_recoder|\/vep\//.test(key)) return respond({ error: unknown }, 400);
            if (/api\.ncbi\.nlm\.nih\.gov\/variation\/v0\/hgvs\/NM_001110792\.2:/.test(key)) return respond({ error: 'no GRCh37 placement' }, 404);
            // ...and without its version it is read, and rejected for its base, as recorded.
            if (/NM_001110792:c\.474G>C/.test(key)) {
                const recorded = recordedFor(key);
                return recorded ? respond(recorded.body, recorded.status) : undefined;
            }
            return undefined;
        }
    });
    assertResolved(assert, page, 'MECP2 c.474G>C, MANE Select version unknown to GRCh37');
    assert.deepEqual(recoderQueries(page).filter((q, i, all) => all.indexOf(q) === i),
        ['NM_001110792.2:c.474G>C', 'NM_001110792:c.474G>C', 'NM_004992.4:c.474G>C']);
    assert.match(page.variant, /Transcript: NM_004992\.4c\.: c\.474G>Cp\.: p\.Thr158=/);
    assert.match(page.panel, /read as NM_004992\.4:c\.474G>C/);
    assert.doesNotMatch(page.panel, /is a newer version than GRCh37's annotation has/, 'that was another transcript');
    console.log('ok gene-symbol name: the version note stays with the transcript it is about');
}

// 22. A name Ensembl re-spells, on a gene that lies on the minus strand.
//     "KRAS c.34_36dup" duplicates a glycine codon; HGVS places it two bases
//     along, and Ensembl answers c.36_38dup. That answer was discarded for not
//     carrying the name as typed, the name went to VEP as a last resort, and
//     VEP's answer — "-/GGT", strand -1: bases of the transcript's strand — was
//     taken as genomic: chr12:g.25398285_25398286insGGT, shown as
//     c.33_34insACC, p.Ala11_Gly12insThr. Another variant, on `main` as well.
{
    // The pieces. VEP's coordinates are the genome's; its bases are the strand's.
    const vep = (allele_string, start, end, strand = -1) => ({ seq_region_name: '12', start, end, strand, allele_string });
    assert.equal(h.buildGenomicHgvsFromVep(vep('-/GGT', 25398286, 25398285)), 'chr12:g.25398285_25398286insACC');
    assert.equal(h.buildGenomicHgvsFromVep(vep('G/A', 25398284, 25398284)), 'chr12:g.25398284C>T', 'c.35G>A is C>T on the genome');
    assert.equal(h.buildGenomicHgvsFromVep(vep('GG/TT', 25398284, 25398285)), 'chr12:g.25398284_25398285delinsAA');
    assert.equal(h.buildGenomicHgvsFromVep(vep('GG/-', 25398284, 25398285)), 'chr12:g.25398284_25398285del');
    // (Asked by a genomic name, or of a gene on the plus strand, VEP answers in the genome's strand.)
    assert.equal(h.buildGenomicHgvsFromVep(vep('-/CCA', 25398284, 25398283, 1)), 'chr12:g.25398283_25398284insCCA');
    assert.equal(h.buildGenomicHgvsFromVep({ seq_region_name: '7', start: 55259515, end: 55259515, strand: 1, allele_string: 'T/G' }), 'chr7:g.55259515T>G');
    assert.equal(h.buildGenomicHgvsFromVep({ seq_region_name: '7', start: 55259515, end: 55259515, allele_string: 'T/G' }), 'chr7:g.55259515T>G');

    // Ensembl's own spelling of the name it was asked.
    const answer = { input: 'NM_004985.5:c.34_36dup', hgvsg: ['NC_000012.11:g.25398284_25398286dup'],
        hgvsc: ['ENST00000256078.4:c.36_38dup', 'NM_004985.5:c.36_38dup', 'NM_033360.4:c.36_38dup'] };
    assert.equal(h.recoderRespelling(answer, 'NM_004985.5:c.34_36dup'), 'NM_004985.5:c.36_38dup');
    assert.equal(h.selectRecoderAlleles([{ GGT: answer }], 'NM_004985.5:c.34_36dup').length, 1, 'the answer to the name is the variant');
    assert.equal(h.selectRecoderAlleles([{ GGT: answer }], 'NM_004985:c.34_36dup').length, 0, 'not the answer to this query');
    // Nothing is left to chance: the answer must be to this query, on this
    // transcript, and an indel of the same size.
    assert.equal(h.recoderRespelling(answer, 'NM_004985.5:c.31_33dup'), '', 'the answer to another name');
    assert.equal(h.recoderRespelling({ ...answer, input: 'KRAS:c.34_36dup' }, 'KRAS:c.34_36dup'), '', 'a gene symbol leaves Ensembl to pick the transcript');
    assert.equal(h.recoderRespelling({ ...answer, hgvsc: ['NM_033360.4:c.36_38dup'] }, 'NM_004985.5:c.34_36dup'), '', 'no name on the transcript asked');
    assert.equal(h.recoderRespelling({ ...answer, hgvsc: ['NM_004985.5:c.36_39dup'] }, 'NM_004985.5:c.34_36dup'), '', 'another size');
    assert.equal(h.recoderRespelling({ ...answer, hgvsc: ['NM_004985.5:c.34_36dup'] }, 'NM_004985.5:c.34_36dup'), '', 'spelled as asked: nothing to say');
    assert.equal(h.recoderRespelling({ input: 'NM_004985.5:c.35G>A', hgvsc: ['NM_004985.5:c.35G>T'] }, 'NM_004985.5:c.35G>A'), '',
        'a substitution is never re-spelled');
    // (An insertion that repeats what precedes it is a duplication to HGVS.)
    assert.equal(h.recoderRespelling({ input: 'NM_004985.5:c.36_37insTGG', hgvsc: ['NM_004985.5:c.36_38dup'] }, 'NM_004985.5:c.36_37insTGG'),
        'NM_004985.5:c.36_38dup');
    assert.equal(h.codingIndelSize('c.34_36dup'), 3);
    assert.equal(h.codingIndelSize('c.5266dupC'), 1);
    assert.equal(h.codingIndelSize('c.1521_1523del'), -3);
    assert.equal(h.codingIndelSize('c.35_36insAAA'), 3);
    assert.equal(h.codingIndelSize('c.34_35delinsTT'), 0);
    assert.equal(h.codingIndelSize('c.489+1_489+4del'), -4);
    assert.equal(h.codingIndelSize('c.*5_*7dup'), 3);
    assert(Number.isNaN(h.codingIndelSize('c.-2_1del')), 'across the start codon: not read off the name');
    assert(Number.isNaN(h.codingIndelSize('c.489_490-1del')), 'across an exon boundary');
    assert.equal(h.codingIndelSize('c.35G>A'), null);

    // As recorded live. Whichever services answer, it is the same allele.
    const cassette = cassettes.krasDup;
    const noRecoder = [/variant_recoder/];
    const noPlacement = [/variant_recoder/, /api\.ncbi\.nlm\.nih\.gov\/variation/];
    const outcomes = [
        ['KRAS c.34_36dup', {}, /g\.: chr12:g\.25398284_25398286dup/],
        ['NM_004985.5:c.34_36dup', {}, /g\.: chr12:g\.25398284_25398286dup/],
        ['KRAS c.34_36dup', { without: noRecoder }, /g\.: chr12:g\.25398286_25398287insCCA/],
        // Ensembl's recoder and NCBI both silent: VEP, by the name, as a last resort.
        ['KRAS c.34_36dup', { without: noPlacement }, /g\.: chr12:g\.25398285_25398286insACC/],
        ['NM_004985.5:c.34_36dup', { without: noPlacement }, /g\.: chr12:g\.25398285_25398286insACC/]
    ];
    for (const [typed, options, genomic] of outcomes) {
        const label = `${typed}${options.without ? ` without ${options.without.length} service(s)` : ''}`;
        const page = await runLookup(cassette, { input: typed, ...options });
        assertResolved(assert, page, label);
        assert.match(page.variant, genomic, label);
        assert.match(page.variant, /VCF \(hg19\): 12-25398280-G-GCCA \(left-aligned\)/, `${label}: one allele, however it was reached`);
        assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.36_38dup/, label);
        assert.doesNotMatch(page.variant + page.panel, /insGGT|c\.33_34insACC|Ala11_Gly12insThr/, `${label}: not the insertion of ACC`);
        // The name is the same variant on both of KRAS's clinical transcripts:
        // there is no other reading to offer.
        assert.equal(readingText(page), null, `${label}: no other reading`);
        if (options.without === noPlacement) {
            assert(asked(page, /\/vep\/human\/hgvs\/NM_004985\.5:c\.34_36dup/).length > 0, `${label}: VEP was asked by the name`);
        } else {
            assert.match(page.variant, /p\.: p\.Gly13dup/, label);
        }
        if (!options.without) {
            assert.match(page.panel, /Ensembl gives this change as c\.36_38dup — the same variant, written where HGVS places it/, label);
            assert.match(page.panel, /Resolve nomenclature · Ensembl\s*8 transcripts/, label);
            assert(!asked(page, /\/vep\/human\/hgvs\/NM_/).length, `${label}: not left to VEP`);
        }
    }
    console.log('ok gene-symbol name: KRAS c.34_36dup is the duplication, whichever service answers');
}

// 22b. The same where ClinVar's own name for the variant is on the gene's
//      second clinical transcript: if both number the variant alike there is
//      nothing to say about which was read. (KRAS c.35G>A as recorded, with
//      ClinVar's preferred name moved to NM_033360.4 — the transcript it used
//      to name KRAS variants on.)
{
    const cassette = cassettes.kras35;
    const searchKey = 'myvariant.info/v1/query?q=KRAS c.35G>A&size=20';
    const search = JSON.parse(JSON.stringify(cassette.responses[searchKey].body).replace(/NM_004985\.5\(KRAS\):c\.35G>A/g, 'NM_033360.4(KRAS):c.35G>A'));
    assert.notDeepEqual(search, cassette.responses[searchKey].body, 'ClinVar\'s preferred name was on NM_004985.5');
    const page = await runLookup(cassette, { route: (url, key, respond) => (key === searchKey ? respond(search) : undefined) });
    assertResolved(assert, page, 'KRAS c.35G>A, named by ClinVar on NM_033360.4');
    assert.match(page.variant, /g\.: chr12:g\.25398284C>T/);
    assert.match(page.clinvar, /Variation ID: 12582/);
    assert.equal(readingText(page), null, 'c.35G>A is the same variant on both transcripts');
    console.log('ok gene-symbol name: a name that is the same variant on both clinical transcripts needs no note');
}

// 22c. And a name that is on no record as typed, because HGVS writes it
//      otherwise: "BRCA1 c.5265_5266insC" inserts a C after a C, which is
//      c.5266dup — 5382insC, ClinVar 17677. (Ensembl did not answer for it in
//      this recording; NCBI placed it.) A record that carries the name as typed
//      on a minor isoform was offered as what might have been meant:
//      "Counted on another transcript the same name is another variant, which
//      ClinVar lists: NM_001407610.1:c.5265_5266insC". Read on MANE Select the
//      name is a known variant; that is no longer mentioned.
{
    const page = await runLookup(cassettes.brca1InsC);
    assertResolved(assert, page, 'BRCA1 c.5265_5266insC');
    assert.match(page.variant, /VCF \(hg19\): 17-41209079-T-TG/);
    assert.match(page.variant, /Transcript: NM_007294\.4 \(MANE Select\)c\.: c\.5266dupp\.: p\.Gln1756ProfsTer74/);
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.equal(readingText(page), null, 'no other reading is offered for a known variant on MANE Select');
    console.log('ok gene-symbol name: BRCA1 c.5265_5266insC is c.5266dup, and nothing else is suggested');
}

// 23. What the other transcript is called. MANE Select has two authorities; a
//     gene neither names one for is read on RefSeq's own pick, and when
//     neither answers in time ClinVar's naming transcript stands in. Offered
//     as the other reading, each is called what it is — the first of the two
//     used to be introduced as "the MANE Select transcript".
{
    const reportKey = 'api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/MECP2/taxon/9606/product_report';
    // (NCBI's list with the gene's MANE Select tag turned into RefSeq Select:
    // the shape of its answer for a gene that has no MANE Select.)
    const report = structuredClone(cassettes.mecp2T158M.responses[reportKey].body);
    for (const t of report.reports[0].product.transcripts) if (t.select_category === 'MANE_SELECT') t.select_category = 'REFSEQ_SELECT';
    const refseq = await runLookup(cassettes.mecp2T158M, { without: [/rest\.ensembl\.org\/lookup/, /preferred_name/],
        route: (url, key, respond) => (key === reportKey ? respond(report) : undefined) });
    assertResolved(assert, refseq, 'MECP2 c.473C>T, RefSeq Select standing in');
    assert.match(refseq.clinvar, /Variation ID: 11811/);
    assert.equal(readingText(refseq), 'MECP2 c.473C>T is read on NM_004992.4, a MANE Plus Clinical transcript of MECP2. '
        + 'If the RefSeq Select transcript was meant, look up NM_001110792.2:c.473C>T.');

    const standIn = await runLookup(cassettes.mecp2T158M, { without: [/rest\.ensembl\.org\/lookup/, /\/datasets\//] });
    assertResolved(assert, standIn, 'MECP2 c.473C>T, ClinVar\'s naming transcript standing in');
    assert.match(standIn.clinvar, /Variation ID: 11811/);
    assert.equal(readingText(standIn), 'MECP2 c.473C>T is read on NM_004992.4, a transcript ClinVar lists this name on '
        + '(whether MECP2 is reported on it could not be learned). If the transcript ClinVar names MECP2 variants on was meant, '
        + 'look up NM_001110792.2:c.473C>T.');
    console.log('ok gene-symbol name: the other transcript is called what it is');
}

// 24. The same mix-up of strands, for a substitution. "NM_004333.4:c.1799T>A"
//     is BRAF V600E on an older version of its transcript. Ensembl answers the
//     name with rows on the current version only, so the name goes to VEP —
//     whose "T/A", strand -1, was shown as chr7:g.140453136T>A: an allele that
//     does not exist (the genome has A there), with no ClinVar record, no rsID
//     and no frequency. On `main` this is what the page showed whenever
//     Ensembl answered before NCBI did.
{
    const cassette = cassettes.brafOlderVersion;
    const situations = [
        ['every service answering', {}],
        ['NCBI not answering', { without: [/api\.ncbi\.nlm\.nih\.gov\/variation/] }],
        ['Ensembl\'s recoder not answering', { without: [/variant_recoder/] }],
        ['neither answering', { without: [/variant_recoder/, /api\.ncbi\.nlm\.nih\.gov\/variation/] }]
    ];
    for (const [label, options] of situations) {
        const page = await runLookup(cassette, options);
        assertResolved(assert, page, `NM_004333.4:c.1799T>A, ${label}`);
        assert.match(page.variant, /g\.: chr7:g\.140453136A>T/, label);
        assert.match(page.variant, /VCF \(hg19\): 7-140453136-A-T/, label);
        assert.match(page.variant, /c\.: c\.1799T>Ap\.: p\.Val600Glu/, label);
        // The allele being known, its record is asked for: V600E is not "N/A".
        assert.match(page.clinvar, /Variation ID: 13961/, label);
        assert.doesNotMatch(page.variant + page.panel, /140453136T>A/, `${label}: not the other strand's bases`);
    }
    // (Where VEP was the route, it was asked by the name, and MyVariant then by the allele.)
    const page = await runLookup(cassette, { without: [/api\.ncbi\.nlm\.nih\.gov\/variation/] });
    assert(asked(page, /\/vep\/human\/hgvs\/NM_004333\.4:c\.1799T>A/).length > 0);
    assert.match(page.panel, /Annotations · MyVariant\.info\s*annotation found/);
    console.log('ok accession name: BRAF V600E on an older transcript version is A>T, with its record, whichever service answers');
}

console.log('Gene-symbol name tests passed');
