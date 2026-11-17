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
    const report = (symbol, transcripts) => ({ product: { symbol, transcripts: transcripts.map(([accession_version, select_category]) => ({ accession_version, select_category })) } });
    const mecp2 = { reports: [report('MECP2', [['NM_001110792.2', 'MANE_SELECT'], ['NM_004992.4', 'MANE_PLUS_CLINICAL'], ['NM_001316337.2'], ['XM_011531166.3']])] };
    assert.deepEqual(h.clinicalTranscriptsFromReport(mecp2, 'mecp2'), { select: 'NM_001110792.2', plusClinical: ['NM_004992.4'], refseqSelect: '' });
    const cdkn2a = { reports: [report('CDKN2A', [['NM_000077.6', 'REFSEQ_SELECT'], ['NM_058195.4', 'MANE_PLUS_CLINICAL']])] };
    assert.deepEqual(h.clinicalTranscriptsFromReport(cdkn2a, 'CDKN2A'), { select: '', plusClinical: ['NM_058195.4'], refseqSelect: 'NM_000077.6' },
        'RefSeq Select is not MANE Select, and is kept apart');
    // A symbol search also returns genes that have the symbol as an alias.
    const ttn = { reports: [report('TTR', [['NM_000371.4', 'MANE_SELECT']]), report('TTN', [['NM_001267550.2', 'MANE_SELECT']])] };
    assert.deepEqual(h.clinicalTranscriptsFromReport(ttn, 'TTN'), { select: 'NM_001267550.2', plusClinical: [], refseqSelect: '' });
    assert.deepEqual(h.clinicalTranscriptsFromReport({}, 'NOSUCHGENE'), { select: '', plusClinical: [], refseqSelect: '' });

    const realFetch = globalThis.fetch;
    const urls = [];
    let status = 200;
    globalThis.fetch = async (address) => {
        urls.push(String(address));
        return new Response(JSON.stringify(mecp2), { status });
    };
    try {
        assert.deepEqual(await h.fetchClinicalTranscripts('mecp2'), { select: 'NM_001110792.2', plusClinical: ['NM_004992.4'], refseqSelect: '' });
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
        assert.deepEqual(reading.alternatives, [{ accession: 'NM_004992.4', kind: 'plus', known: true }]);
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
    // if the name does not fit — and a namesake elsewhere is now worth a mention.
    {
        const reading = choose([p251l], 'c.473C>T', known);
        assert.equal(reading.accession, 'NM_001110792.2');
        assert.equal(reading.hit, null);
        assert.deepEqual(reading.untried, ['NM_004992.4']);
        assert.deepEqual(reading.alternatives, [{ accession: 'NM_001316337.2', kind: 'listed', known: true }]);
        assert.deepEqual(choose([], 'c.473C>T', known).alternatives, []);
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
        alternatives: [{ accession: 'NM_004992.4', kind: 'plus', known: true }], rejected: [] }).links,
    [{ variant: 'NM_004992.4:c.473C>T', lead: 'Counted on the MANE Plus Clinical transcript the same name is another variant, which ClinVar lists:' }]);
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
    assert.equal(readingNote(page).textContent, 'CDKN2A c.3G>C is read on NM_058195.4, the transcript ClinVar names this variant on '
        + '(a MANE Plus Clinical transcript of CDKN2A). If the MANE Select transcript was meant, look up NM_000077.6:c.3G>C.');
    assert.deepEqual(noteLinks(page), ['?variant=NM_000077.6:c.3G>C']);
    // A name ClinVar gives on MANE Select itself has nothing to add (scenario 8).
    assert.equal(readingNote(await runLookup(cassettes.kras35)), null);
    assert.equal(readingNote(await runLookup(cassettes.brca1Dup)), null);
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

console.log('Gene-symbol name tests passed');
