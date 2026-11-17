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
//    recoder is not waited for once VEP has answered. Shown with the page's
//    waits at their real length (it gives the recoder 6 s where it has to) and
//    the recoder's answer held back until the lookup is over: the lookup must
//    finish regardless, and long before that wait would have run out.
{
    let release;
    const held = new Promise((resolve) => { release = resolve; });
    let asked = 0;
    const started = Date.now();
    const lookup = runLookup(cassettes.rs80357906, {
        input: 'chr17:g.41209082dup',
        timerCap: 6000,
        route: async (url, key, respond) => {
            if (!/variant_recoder/.test(key)) return undefined;
            asked += 1;
            await held;
            return respond({ error: 'Service Unavailable' }, 503);
        }
    });
    const page = await Promise.race([lookup, new Promise(resolve => setTimeout(() => resolve(null), 3500))]);
    const took = Date.now() - started;
    release();
    assert(page, 'the lookup finished while the recoder had still not answered');
    assert(asked >= 1, 'which it had been asked');
    assert(took < 3000, `and without sitting out the 6 s it allows the recoder (${took} ms)`);
    assertResolved(assert, page, 'chr17:g.41209082dup');
    assert.match(page.variant, /Gene: BRCA1Transcript: NM_007294\.4/);
    assert.match(page.variant, /p\.: p\.Gln1756ProfsTer74/);
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*still answering — not needed: the variant was resolved without it/);
    console.log('ok genomic lookup: a record with no transcript rows gets them from VEP, without waiting for the recoder');
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

// 7. The same lateness, with an rsID — where the recoder's answer holds several
//    alleles and its rows are at first the FIRST one's. They are put right for
//    the allele the lookup chose, but that used to happen once, before the wait
//    for VEP: answering during that wait, with VEP then giving up, the recoder
//    left V600A's rows on V600E's page —
//
//        g.: chr7:g.140453136A>T   ClinVar 13961   c.: c.1799T>C
//
//    and PubMed was asked about "BRAF c.1799T>C". (Found in review.)
{
    const cassette = cassettes.rs113488022;
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    const first = Object.keys(cassette.responses['grch37.rest.ensembl.org/variant_recoder/human/rs113488022?content-type=application/json'].body[0])[0];
    assert.equal(first, 'G', 'the recorded answer lists A>G (V600A) first');
    const lateRecoder = (record = null) => runLookup(cassette, {
        timerCap: 150,
        settleMs: 700,
        route: async (url, key, respond, options, { document }) => {
            if (/variant_recoder/.test(key)) {
                // Answers once VEP is being waited for: the allele is chosen by then.
                const vep = () => document.querySelector('[data-step="vep"]');
                for (let i = 0; i < 3000 && !(vep() && /lp-running/.test(vep().className)); i++) await wait(1);
                const recorded = cassette.responses[key];
                return respond(recorded.body, recorded.status);
            }
            if (/\/vep\//.test(key)) { await wait(900); return respond({ error: 'late' }, 503); }
            if (record && /myvariant\.info\/v1\/variant\/chr7:g\.140453136A>T/.test(key)) return respond(record);
            return undefined;
        }
    });
    {
        const page = await lateRecoder();
        assertResolved(assert, page, 'rs113488022, recoder late, VEP gives up');
        assert.match(page.variant, /g\.: chr7:g\.140453136A>T/);
        assert.match(page.panel, /Consequences · Ensembl VEPEnsembl VEP did not answer in time/);
        // The record's own rows say what the change does; the recoder's do not.
        assert.match(page.variant, /Transcript: NM_004333\.4c\.: c\.1799T>Ap\.: p\.Val600Glu \(p\.V600E\)Effect: missense_variant/);
        assert.doesNotMatch(page.variant, /c\.1799T>C/);
        assert.deepEqual(page.api('pubmed').map(r => (r.match(/variant=([^&]+)/) || [])[1]).filter(Boolean), ['V600E']);
    }
    // A record with no transcript rows of its own, as a dbSNP-only one is: no
    // snpEff or dbNSFP rows, and no ClinVar names to make rows from. The
    // recoder's rows are then all there is to show — for the allele chosen.
    const bare = structuredClone(cassette.responses['myvariant.info/v1/variant/chr7:g.140453136A>T'].body);
    for (const block of ['snpeff', 'dbnsfp', 'cadd', 'clinvar', 'hgvsc', 'hgvsp', 'civic', 'cosmic', 'docm', 'cgi']) delete bare[block];
    const noOtherAllele = (page, label) => {
        assertResolved(assert, page, label);
        assert.match(page.variant, /g\.: chr7:g\.140453136A>T/, label);
        assert.match(page.variant, /Transcript: NM_004333\.\d+(?: \(MANE Select\))?c\.: c\.1799T>A/, `${label}: the recoder's row for A>T`);
        assert.doesNotMatch(page.variant, /c\.1799T>[CG]/, `${label}: no row of another allele`);
        assert(!page.requests.some(r => /c\.1799T>[CG]/.test(r)), `${label}: and no source asked about one`);
    };
    // Arriving late, during the wait for VEP: taken for the chosen allele then.
    noOtherAllele(await lateRecoder(bare), 'rs113488022, recoder late, record without rows');
    // Arriving at once, before the allele is chosen: put right when it is.
    noOtherAllele(await runLookup(cassette, {
        route: (url, key, respond) => {
            if (/\/vep\//.test(key)) return respond({ error: 'Service Unavailable' }, 503);
            if (/myvariant\.info\/v1\/variant\/chr7:g\.140453136A>T/.test(key)) return respond(bare);
            return undefined;
        }
    }), 'rs113488022, recoder prompt, VEP down, record without rows');
    // The same window, with VEP answering after all — from Ensembl's own
    // transcript set, its RefSeq one having failed. That answer is not allowed
    // to displace the record's RefSeq row, and the panel says "the record's own
    // are shown": but the recoder's rows, arrived meanwhile, were what was
    // shown, "p.: N/A · Effect: N/A". (Found in review.)
    {
        const refseqKey = 'grch37.rest.ensembl.org/vep/human/hgvs/7:g.140453136A>T?content-type=application/json&hgvs=1&refseq=1&mane=1&canonical=1';
        const recorded = cassette.responses[refseqKey].body[0];
        const ensemblOnly = [{ ...recorded, transcript_consequences: recorded.transcript_consequences.slice(0, 1).map(row => ({ ...row,
            gene_id: 'ENSG00000157764', transcript_id: 'ENST00000288602', hgvsc: 'ENST00000288602.6:c.1799T>A', hgvsp: 'ENSP00000288602.6:p.Val600Glu' })) }];
        let recoderAnswered = false;
        const page = await runLookup(cassette, {
            timerCap: 600,
            settleMs: 700,
            route: async (url, key, respond, options, { document }) => {
                if (/variant_recoder/.test(key)) {
                    const vep = () => document.querySelector('[data-step="vep"]');
                    for (let i = 0; i < 3000 && !(vep() && /lp-running/.test(vep().className)); i++) await wait(1);
                    setTimeout(() => { recoderAnswered = true; }, 20);
                    const answer = cassette.responses[key];
                    return respond(answer.body, answer.status);
                }
                if (/\/vep\//.test(key)) {
                    if (/refseq=1/.test(key)) return respond({ error: 'Internal Server Error' }, 500);
                    for (let i = 0; i < 3000 && !recoderAnswered; i++) await wait(1);
                    return respond(ensemblOnly);
                }
                return undefined;
            }
        });
        assertResolved(assert, page, 'rs113488022, recoder late, VEP RefSeq set down');
        assert(recoderAnswered, 'the recoder answered while VEP was being waited for');
        assert.match(page.panel, /Consequences · Ensembl VEPRefSeq transcripts unavailable — the record's own are shown/);
        assert.match(page.variant, /Transcript: NM_004333\.4c\.: c\.1799T>Ap\.: p\.Val600Glu \(p\.V600E\)Effect: missense_variant/, 'and they are');
        assert.doesNotMatch(page.variant, /p\.: N\/A/);
    }
    console.log('ok rsID lookup: rows that arrive late are the chosen allele\'s');
}

// 8. A lookup's recoder outliving it. The recoder takes Ensembl up to a minute
//    over an rsID and is no longer waited for, so it gives up while the NEXT
//    lookup is running — and its failure used to be written onto that one:
//    "rs80357906", then "BRCA1 c.5382insC", ended in "Lookup could not
//    complete … upstream outage … retry", and ClinVar's other names, which
//    resolve that name, were never asked. (Found in review.)
{
    const first = cassettes.rs80357906;
    const second = cassettes.brca1OtherNames;
    const both = { input: first.input, reference: [...first.reference, ...second.reference], api: { ...(second.api || {}) },
        responses: { ...first.responses, ...second.responses } };
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    let staleAttempts = 0;
    const page = await runLookup(both, {
        followUp: { input: 'BRCA1 c.5382insC' },
        timerCap: 60,
        settleMs: 400,
        route: async (url, key, respond, options, { document }) => {
            if (/variant_recoder\/human\/rs80357906/.test(key)) {
                // The first lookup's recoder fails only once the second is under way…
                for (let i = 0; i < 4000 && !/BRCA1 c\.5382insC/.test(document.getElementById('lookupProgress')?.textContent || ''); i++) await wait(1);
                staleAttempts += 1;
                return respond({ error: 'Service Unavailable' }, 503);
            }
            if (/variant_recoder\/human\/NM_007294\.4:c\.5382insC/.test(key)) {
                // …and the second lookup's own question is answered after all of
                // those attempts, and the failure that follows them, have happened.
                let last = -1;
                let quiet = 0;
                for (let i = 0; i < 6000 && quiet < 250; i++) {
                    await wait(1);
                    if (staleAttempts !== last) { last = staleAttempts; quiet = 0; } else if (staleAttempts > 0) quiet += 1;
                }
            }
            return undefined;
        }
    });
    assert(staleAttempts >= 2, 'the first lookup\'s recoder was still trying during the second');
    assertResolved(assert, page, 'BRCA1 c.5382insC, after rs80357906');
    assert.match(page.panel, /You entered\s*BRCA1 c\.5382insC/);
    assert.match(page.variant, /c\.: c\.5266dup/);
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.equal(page.api('clinvar-protein').filter(r => /retmax=200/.test(r)).length, 1, 'ClinVar\'s other names were asked');
    assert.doesNotMatch(page.panel, /Resolve nomenclature · Ensembl\s*HTTP 503/, 'the first lookup\'s failure is not on this panel');
    console.log('ok rsID lookup: a finished lookup\'s recoder does not fail the next one');
}

// 9. ...nor does it take the next lookup's rows away by answering. The page
//    keeps the recoder's last answer to read it again later; a finished
//    lookup's answer landing there displaced the current one's, whose
//    transcript rows — all it had, with no record and no VEP — then went.
{
    const first = cassettes.rs113488022;
    const second = cassettes.cftrF508del;
    const both = { input: first.input, reference: [...first.reference, ...second.reference], responses: { ...first.responses, ...second.responses } };
    const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
    let vepAsked = false;
    let staleAnswered = false;
    const page = await runLookup(both, {
        followUp: { input: 'chr7:g.117199645_117199647del' },
        settleMs: 300,
        route: async (url, key, respond, options) => {
            if (/variant_recoder\/human\/rs113488022/.test(key)) {
                // The first lookup's recoder answers — properly, with its three
                // alleles — only once the second lookup is waiting for VEP.
                for (let i = 0; i < 6000 && !vepAsked; i++) await wait(1);
                const recorded = first.responses['grch37.rest.ensembl.org/variant_recoder/human/rs113488022?content-type=application/json'];
                setTimeout(() => { staleAnswered = true; }, 30);
                return respond(recorded.body, recorded.status);
            }
            // The second lookup has only the recoder to go by: MyVariant holds
            // nothing for the allele under any of its ids, and VEP fails once
            // the stale answer is in.
            if (/myvariant\.info\/v1\/variant\/chr7:g\.1171996/.test(key)) return respond({ success: false, error: 'Not Found.' }, 404);
            if (url.hostname === 'myvariant.info' && String(options.method).toUpperCase() === 'POST') {
                return respond(new URLSearchParams(String(options.body || '')).get('ids').split(',').map(id => ({ query: id, notfound: true })));
            }
            if (/\/vep\/human\/hgvs\/7:g\.1171996/.test(key)) {
                vepAsked = true;
                for (let i = 0; i < 6000 && !staleAnswered; i++) await wait(1);
                return respond({ error: 'Service Unavailable' }, 503);
            }
            return undefined;
        }
    });
    assert(staleAnswered, 'the first lookup\'s recoder answered during the second');
    assertResolved(assert, page, 'chr7:g.117199645_117199647del, after rs113488022');
    assert.match(page.panel, /You entered\s*chr7:g\.117199645_117199647del/);
    assert.match(page.variant, /Transcript: NM_000492\.\d+(?: \(MANE Select\))?c\.: c\.1521_1523del/, 'its own recoder rows are still there');
    assert.doesNotMatch(page.variant, /1799T|BRAF/);
    assert.doesNotMatch(page.variant, /p\.Phe508del/, '(the recoder\'s rows: MyVariant and VEP contributed none)');
    console.log('ok rsID lookup: a finished lookup\'s recoder answer does not displace the next one\'s');
}

console.log('rsID lookup tests passed');
