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
 * A name on a gene symbol is now read on the gene's standard transcript, MANE
 * Select, and from there on the lookup is for "NM_…:c.X".
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

// ── Two sources for the standard transcript ────────────────────────────────
{
    const timing = { maneAlone: 30, total: 120 };
    const after = (ms, value) => () => new Promise(resolve => setTimeout(() => resolve(value), ms));
    const std = (mane, naming) => h.standardTranscriptForGene('PTEN', mane, naming, timing);
    // Ensembl's MANE record is the authority when it answers in time.
    assert.deepEqual(await std(after(1, ['NM_000314.8', 'ENST00000371953.8']), after(1, 'NM_999999.1')),
        { accession: 'NM_000314.8', basis: 'the MANE Select transcript of PTEN' });
    // Slow: ClinVar's naming transcript stands in.
    assert.deepEqual(await std(after(500, ['NM_000314.8']), after(1, 'NM_000314.8')),
        { accession: 'NM_000314.8', basis: 'the transcript ClinVar names PTEN variants on' });
    // ClinVar unclear: MANE gets the rest of the time.
    assert.deepEqual(await std(after(60, ['NM_000314.8']), after(1, '')),
        { accession: 'NM_000314.8', basis: 'the MANE Select transcript of PTEN' });
    // Neither can say — including by failing.
    assert.equal(await std(after(1, []), after(1, '')), null);
    assert.equal(await std(() => { throw new Error('down'); }, () => Promise.reject(new Error('down'))), null);
    assert.equal(await std(after(500, ['NM_000314.8']), after(500, 'NM_000314.8')), null, 'both too slow');
    // A gene with no MANE Select (only an Ensembl accession, or none) is not a RefSeq answer.
    assert.equal(await std(after(1, ['ENST00000371953.8']), after(1, '')), null);
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

// 3. Ensembl's MANE lookup not answering: ClinVar's naming transcript stands in.
{
    const page = await runLookup(cassettes.brca2GeneInsCA, { without: [/\/lookup\//] });
    assertResolved(assert, page, 'BRCA2 c.1083_1084insCA, no MANE lookup');
    assert.match(page.panel, /on NM_000059\.4, the transcript ClinVar names BRCA2 variants on/);
    assert.deepEqual(recoderQueries(page), ['NM_000059.4:c.1083_1084insCA']);
    assert.match(page.variant, /c\.: c\.1083_1084insCA/);
    console.log('ok gene-symbol name: ClinVar\'s naming transcript when Ensembl is slow');
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
    const page = await runLookup(cassettes.pten512, { without: [/\/lookup\//, /clinvar\.gene\.symbol/] });
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

console.log('Gene-symbol name tests passed');
