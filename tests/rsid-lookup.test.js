/*
 * rsID lookups.
 *
 * An rsID used to be handed to Ensembl's variant recoder and waited for. For
 * rs80357906 (BRCA1 5382insC) Ensembl needs about 45 seconds; this page allows
 * 12, tries four times, and then falls back to the MyVariant search it could
 * have started with. The record that search finds is one of MyVariant's
 * dbSNP-only ones — an rsID and a gene, no transcript — so after 49.5 seconds
 * the page read "c.: N/A, p.: N/A, Effect: delins".
 *
 * And an rsID is not a variant when the site has several alleles. Which one was
 * shown depended on the order Ensembl listed them in: rs113488022, the rsID
 * everyone uses for BRAF V600E, came back as V600A.
 *
 * The scenarios replay answers recorded from the live services with the page's
 * real timeouts, so a request Ensembl did not answer in time is a timeout here.
 *
 * Run with: node tests/rsid-lookup.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const asked = (page, pattern) => page.requests.filter(r => pattern.test(r));

// ── Which allele an rsID means ─────────────────────────────────────────────
{
    const rcvs = n => ({ rcv: Array.from({ length: n }, (_, i) => ({ accession: `RCV${i}` })) });
    const hits = [
        { _id: 'chr7:g.140453136A>C', clinvar: rcvs(2) },
        { _id: 'chr7:g.140453136A>G', clinvar: rcvs(3), gnomad_exome: { af: { af: 0.5 } } },
        { _id: 'chr7:g.140453136A>T', clinvar: rcvs(40) }
    ];
    assert.deepEqual(h.rankRsidAlleles(hits).map(hit => hit._id),
        ['chr7:g.140453136A>T', 'chr7:g.140453136A>G', 'chr7:g.140453136A>C'], 'most ClinVar submissions first');
    // Without ClinVar, the common allele; with nothing to go on, MyVariant's order.
    assert.deepEqual(h.rankRsidAlleles([
        { _id: 'rare', gnomad_genome: { af: { af: 0.0001 } } },
        { _id: 'common', gnomad_exome: { af: { af: 0.3 } } },
        { _id: 'unseen' }
    ]).map(hit => hit._id), ['common', 'rare', 'unseen']);
    assert.deepEqual(h.rankRsidAlleles([{ _id: 'a' }, { _id: 'b' }, { _id: 'c' }]).map(hit => hit._id), ['a', 'b', 'c']);
    assert.deepEqual(h.rankRsidAlleles([{ _id: 'only' }]).map(hit => hit._id), ['only']);
    assert.deepEqual(h.rankRsidAlleles(null), []);
    // A single ClinVar object (not an array) and a single RCV both count.
    assert.deepEqual(h.rankRsidAlleles([{ _id: 'none' }, { _id: 'one', clinvar: { rcv: { accession: 'RCV1' } } }]).map(hit => hit._id),
        ['one', 'none']);
}

// ── A deadline that does not cancel ────────────────────────────────────────
{
    assert.equal(await h.settledWithin(Promise.resolve('answer'), 50, 'late'), 'answer');
    await assert.rejects(h.settledWithin(Promise.reject(new Error('down')), 50, 'late'), /down/);
    let finished = false;
    const slow = new Promise(resolve => setTimeout(() => { finished = true; resolve('eventually'); }, 40));
    await assert.rejects(h.settledWithin(slow, 5, 'did not answer in time'), err => err.timedOut === true && /did not answer in time/.test(err.message));
    assert.equal(await slow, 'eventually', 'the request itself is left to finish');
    assert.equal(finished, true);
    // A request that fails after its deadline must not surface as an unhandled rejection.
    const failsLate = new Promise((_, reject) => setTimeout(() => reject(new Error('late failure')), 20));
    await assert.rejects(h.settledWithin(failsLate, 5, 'gave up'), /gave up/);
    await new Promise(resolve => setTimeout(resolve, 30));
}

// ── A bare c. does not displace a row that says what the change does ───────
{
    // As MyVariant serves TP53 R175H: snpEff on NM_000546.5, ClinVar on .6.
    const annotation = {
        snpeff: { ann: [{ feature_id: 'NM_000546.5', feature_type: 'transcript', hgvs_c: 'c.524G>A', hgvs_p: 'p.Arg175His', effect: 'missense_variant', genename: 'TP53' }] },
        clinvar: { hgvs: { coding: ['NM_000546.6:c.524G>A', 'NM_001126112.3:c.524G>A'] } }
    };
    const rows = h.buildCanonicalFromAnnotation(annotation, null).transcriptsList;
    assert.deepEqual(rows.map(t => `${t.transcript}:${t.cDNA} ${t.protein}`.trim()),
        ['NM_000546.5:c.524G>A p.Arg175His', 'NM_001126112.3:c.524G>A'],
        'the same transcript\'s bare c. is dropped; another transcript\'s is kept');
    // ...and so the MANE flag on the bare row cannot win the selection.
    const marked = h.markManeTranscripts(rows, ['NM_000546.6']);
    assert.equal(marked[h.selectCanonicalTranscript(marked)].protein, 'p.Arg175His');
    // "c.5266dupC" (snpEff) and "c.5266dup" (ClinVar) are one change.
    const dup = h.buildCanonicalFromAnnotation({
        snpeff: { ann: [{ feature_id: 'NM_007294.3', feature_type: 'transcript', hgvs_c: 'c.5266dupC', hgvs_p: 'p.Gln1756fs', effect: 'frameshift_variant' }] },
        clinvar: { hgvs: { coding: ['NM_007294.4:c.5266dup'] } }
    }, null).transcriptsList;
    assert.equal(dup.length, 1);
}
console.log('rsID unit checks passed');

// 1. rs80357906. Ensembl's recoder never answers within its deadline (as
//    recorded). Was: 49.5 s, then no transcript, no coding change, no protein.
{
    const page = await runLookup(cassettes.rs80357906);
    assertResolved(assert, page, 'rs80357906');
    assert.match(page.variant, /g\.: chr17:g\.41209082dup/);
    assert.match(page.variant, /Gene: BRCA1Transcript: NM_007294\.4/);
    assert.match(page.variant, /c\.: c\.5266dupp\.: p\.Gln1756ProfsTer74 \(p\.Q1756Pfs\*74\)Effect: frameshift_variant/);
    assert.doesNotMatch(page.variant, /names \d+ alleles/, 'one allele: nothing to choose');
    assert.match(page.panel, /Annotations · MyVariant\.info\s*rs80357906 is one allele here: chr17:g\.41209082dup/);
    assert.match(page.panel, /Consequences · Ensembl VEP\s*6 transcripts/);
    // MyVariant is asked once, at the start — not after the recoder gives up.
    assert.deepEqual(asked(page, /myvariant\.info\/v1\/query/), ['myvariant.info/v1/query?q=rs80357906&size=20']);
    const order = page.requests.map(r => (/\/v1\/query/.test(r) ? 'search' : /variant_recoder/.test(r) ? 'recoder' : null)).filter(Boolean);
    assert.equal(order.indexOf('search') <= 1, true, 'the search does not wait for the recoder');
    assert.deepEqual(asked(page, /\/vep\/human\/hgvs\//).map(r => r.match(/hgvs\/([^?]+)/)[1]), ['17:g.41209082dup']);
    console.log('ok rsID lookup: rs80357906 resolves without the recoder, with its transcripts');
}

// 2. rs113488022 names three alleles. Was: A>G, p.Val600Ala — the first one
//    Ensembl lists.
{
    const page = await runLookup(cassettes.rs113488022);
    assertResolved(assert, page, 'rs113488022');
    assert.match(page.variant, /g\.: chr7:g\.140453136A>T/);
    assert.match(page.variant, /c\.: c\.1799T>Ap\.: p\.Val600Glu \(p\.V600E\)Effect: missense_variant/);
    assert.match(page.clinvar, /Variation ID: 13961/);
    // The choice is stated where the result is read, with the others a click away.
    assert.match(page.variant,
        /rs113488022 names 3 alleles\. This is the best documented one; the others are chr7:g\.140453136A>C, chr7:g\.140453136A>G\./);
    const note = page.document.querySelector('.rsid-alleles-note');
    assert.deepEqual([...note.querySelectorAll('a')].map(a => decodeURIComponent(a.getAttribute('href'))),
        ['?variant=chr7:g.140453136A>C', '?variant=chr7:g.140453136A>G']);
    assert.match(page.panel, /rs113488022 names 3 alleles here — showing the best documented, chr7:g\.140453136A>T/);
    // Only the chosen allele's record is fetched.
    assert.deepEqual(asked(page, /myvariant\.info\/v1\/variant\//), ['myvariant.info/v1/variant/chr7:g.140453136A>T']);
    console.log('ok rsID lookup: rs113488022 is shown as V600E, and says it chose');
}

// 3. rs28934578 (TP53 R175H). Was: 51.8 s, then "p.: N/A, Effect: N/A".
{
    const page = await runLookup(cassettes.rs28934578);
    assertResolved(assert, page, 'rs28934578');
    assert.match(page.variant, /g\.: chr17:g\.7578406C>T/);
    assert.match(page.variant, /Transcript: NM_000546\.6/);
    assert.match(page.variant, /c\.: c\.524G>Ap\.: p\.Arg175His \(p\.R175H\)Effect: missense_variant/);
    assert.match(page.variant, /rs28934578 names 3 alleles/);
    console.log('ok rsID lookup: rs28934578 is R175H, with its protein change');
}

// 4. The same with VEP not answering either. MyVariant's own rows are then
//    all there is — and the one that says what the change does is the one shown.
{
    const page = await runLookup(cassettes.rs28934578, { without: [/\/vep\//] });
    assertResolved(assert, page, 'rs28934578, no VEP');
    assert.match(page.variant, /Transcript: NM_000546\.5c\.: c\.524G>Ap\.: p\.Arg175His \(p\.R175H\)Effect: missense_variant/);
    assert.doesNotMatch(page.variant, /p\.: N\/A/);
    assert.match(page.panel, /▲Consequences · Ensembl VEP\s*(?:HTTP 503|Ensembl VEP did not answer in time)/);
    console.log('ok rsID lookup: without VEP, the described row is shown, not the bare one');
}

// 5. A direct genomic lookup of a dbSNP-only record takes the same road: the
//    recoder is not waited for once VEP has answered.
{
    const page = await runLookup(cassettes.rs80357906, { input: 'chr17:g.41209082dup' });
    assertResolved(assert, page, 'chr17:g.41209082dup');
    assert.match(page.variant, /Gene: BRCA1Transcript: NM_007294\.4/);
    assert.match(page.variant, /p\.: p\.Gln1756ProfsTer74/);
    console.log('ok genomic lookup: a record with no transcript rows gets them from VEP');
}

// 6. Not waiting for the recoder means it answers later — typically a second
//    or two after VEP, in the middle of what remains of the lookup. Its rows
//    carry no protein change, and arriving then they replaced VEP's: F508del by
//    its g. name came out "c.1521_1523del · p.: N/A · Effect: N/A" in a real
//    browser, whenever Ensembl's recoder took 3 s rather than 1. Here the
//    recoder comes back only once the panel shows VEP's rows taken up, while
//    the gene's MANE lookup keeps the label step that follows waiting.
{
    const cassette = cassettes.cftrF508del;
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    const late = async (outcome) => {
        const order = [];
        const page = await runLookup(cassette, {
            input: 'chr7:g.117199645_117199647del',
            timerCap: 150,
            settleMs: 500,
            route: async (url, key, respond, options, { document }) => {
                if (/variant_recoder/.test(key)) {
                    const vepStep = () => document.querySelector('[data-step="vep"]')?.textContent || '';
                    for (let i = 0; i < 2000 && !/\d transcripts?/.test(vepStep()); i++) await wait(1);
                    order.push(`VEP: ${vepStep().match(/\d transcripts?/)?.[0]}`, `recoder ${outcome}`);
                    const recorded = cassette.responses[key];
                    return outcome === 'answers' ? respond(recorded.body, recorded.status) : respond({ error: 'Bad Gateway' }, 400);
                }
                if (/rest\.ensembl\.org\/lookup\//.test(key)) { await wait(450); return respond({ error: 'late' }, 503); }
                return undefined;
            }
        });
        return { page, order };
    };
    for (const outcome of ['answers', 'fails']) {
        const { page, order } = await late(outcome);
        assertResolved(assert, page, `F508del, recoder ${outcome} late`);
        assert.deepEqual(order, ['VEP: 2 transcripts', `recoder ${outcome}`], 'the recoder came back after the rows were decided');
        assert.match(page.variant, /Transcript: NM_000492\.4c\.: c\.1521_1523delp\.: p\.Phe508del \(p\.F508del\)Effect: inframe_deletion/, outcome);
        assert.doesNotMatch(page.variant, /p\.: N\/A|Effect: N\/A/, outcome);
        assert.match(page.panel, /CodingNM_000492\.4:c\.1521_1523delProteinp\.Phe508delEffectinframe_deletion/, outcome);
        assert.match(page.panel, outcome === 'answers'
            ? /◍Resolve nomenclature · Ensemblanswered after the result was ready — not needed: the variant was resolved without it/
            : /✕Resolve nomenclature · EnsemblHTTP 400/);
    }
    // Answering in good time, it is used as before (its rows, paired with VEP's).
    const prompt = await runLookup(cassette, { input: 'chr7:g.117199645_117199647del' });
    assert.match(prompt.panel, /●Resolve nomenclature · Ensembl5 transcripts/);
    assert.match(prompt.variant, /Transcript: NM_000492\.4 \(MANE Select\)c\.: c\.1521_1523delp\.: p\.Phe508del/);
    console.log('ok genomic lookup: a recoder answer that arrives late does not replace VEP\'s rows');
}

console.log('rsID lookup tests passed');
