/*
 * Where the GRCh37 reference bases come from.
 *
 * The page read them from one place, Ensembl's sequence service. That service
 * fails on its own at times — HTTP 500 for every region, or answers that take
 * many seconds, while Ensembl's recoder and VEP go on answering — and
 * everything that reads the reference then went without it:
 *
 *   X 20156693 EIF1AX C CAGATTCATTCTCATTCTTACCCCTGCGTCT
 *     Confirm allele · GRCh37 reference   not confirmed — reference sequence unavailable
 *     Annotations · MyVariant.info        its other ids for this allele were not checked
 *     c.: c.34_63dup · p.: N/A            "protein change not shown: … no name could be
 *                                          had for it written that way"
 *
 * The protein name was the casualty that showed: EIF1AX is on the reverse
 * strand, and Ensembl's protein name for an in-frame insertion there is shown
 * only for a spelling known to be left-aligned — which takes the reference to
 * know. With it the row is p.Arg12_Ser21dup.
 *
 * UCSC serves the same assembly (hg19), base for base. Where Ensembl fails, or
 * has not answered after a short head start, UCSC is asked as well and the
 * first of the two to give the bases is taken.
 *
 * The scenarios replay answers recorded from the live services.
 *
 * Run with: node tests/reference-source.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const wait = ms => new Promise(resolve => setTimeout(resolve, ms));

// ── The pieces ─────────────────────────────────────────────────────────────
{
    // UCSC's names for GRCh37's chromosomes. Its "chrM" is an older
    // mitochondrial sequence; the one GRCh37 uses (the rCRS) is its "chrMT".
    assert.deepEqual(['1', '9', '10', '22', 'X', 'Y', 'chrX', 'y', 'chr7'].map(h.ucscHg19Chromosome),
        ['chr1', 'chr9', 'chr10', 'chr22', 'chrX', 'chrY', 'chrX', 'chrY', 'chr7']);
    assert.deepEqual(['MT', 'chrMT', 'mt'].map(h.ucscHg19Chromosome), ['chrMT', 'chrMT', 'chrMT']);
    assert.deepEqual(['M', 'chrM', '0', '23', '220', 'GL000207.1', 'X_random', '', null, undefined].map(h.ucscHg19Chromosome),
        ['', '', '', '', '', '', '', '', '', ''], 'nothing else has a name there');

    const realFetch = globalThis.fetch;
    const realWarn = console.warn;
    const real = { head: h.API_TIMEOUT_MS.referenceHeadStart };
    // (The page notes each failed read on the console.)
    console.warn = () => {};
    // (GRCh37 X:20156693-20156725, as both services gave it. UCSC writes
    // repeats in lower case.)
    const bases = 'CAGATTCATTCTCATTCTTACCCCTGCGTCTGT';
    const ensemblUrl = 'https://grch37.rest.ensembl.org/sequence/region/human/X:20156693..20156725?content-type=text/plain';
    const ucscUrl = 'https://api.genome.ucsc.edu/getData/sequence?genome=hg19&chrom=chrX&start=20156692&end=20156725';
    const ucscAnswer = (over = {}) => ({ downloadTime: '2026:11:18T14:34:35Z', downloadTimeStamp: 1795012475,
        genome: 'hg19', chrom: 'chrX', start: 20156692, end: 20156725, dna: 'cagattcattctcattcttaccccTGCGTCTGT', ...over });
    const text = (body, status = 200) => () => new Response(body, { status });
    const json = (body, status = 200) => () => new Response(JSON.stringify(body), { status });
    const down = text('<!doctype html><html lang="en"><title>Internal Server Error</title></html>', 500);
    // A service that takes `ms` to answer — or never does, until it is given up on.
    const after = (ms, answer) => (url, options) => new Promise((resolve, reject) => {
        const timer = setTimeout(() => resolve(answer()), ms);
        options.signal?.addEventListener('abort', () => {
            clearTimeout(timer);
            reject(Object.assign(new Error('This operation was aborted'), { name: 'AbortError' }));
        });
    });
    const never = after(3600000, () => new Response(''));
    const serve = ({ ensembl, ucsc }) => {
        const asked = [];
        globalThis.fetch = async (address, options = {}) => {
            const url = String(address);
            asked.push(url);
            return (url.includes('ensembl.org') ? ensembl : ucsc)(url, options);
        };
        return asked;
    };
    const read = (ms = 5000) => h.fetchGrch37ReferenceSequence('chrX', 20156693, 20156725, ms);
    try {
        // Ensembl answers: its bases, and nobody else is asked.
        let asked = serve({ ensembl: text(`${bases}\n`), ucsc: json(ucscAnswer()) });
        assert.equal(await read(), bases);
        assert.deepEqual(asked, [ensemblUrl]);

        // Ensembl fails — as it did, for every region: UCSC is asked at once.
        asked = serve({ ensembl: down, ucsc: json(ucscAnswer()) });
        let began = Date.now();
        assert.equal(await read(), bases, 'UCSC\'s bases, in upper case');
        assert.deepEqual(asked, [ensemblUrl, ucscUrl], 'the same stretch, its start counted from 0');
        assert(Date.now() - began < real.head, 'without waiting out the head start');
        // ...and when what Ensembl sends is not the bases asked for.
        for (const notTheBases of [text('<!doctype html><html></html>'), text(bases.slice(1)), text(`${bases}A`), text(''), text('CAGATTCATTCTCATTCTTACCCCTGCGTCTGR')]) {
            asked = serve({ ensembl: notTheBases, ucsc: json(ucscAnswer()) });
            assert.equal(await read(), bases);
            assert.equal(asked.length, 2);
        }

        // Neither has them: the error is Ensembl's, the service the page reads from.
        asked = serve({ ensembl: down, ucsc: json({ error: 'server overloaded' }, 500) });
        await assert.rejects(read(), /^Error: Ensembl sequence request failed \(500\)$/);
        assert.equal(asked.length, 2);
        // What UCSC sends has to be the stretch asked for, of this assembly, and bases.
        for (const [why, answer] of [
            ['another stretch', ucscAnswer({ start: 20156693 })],
            ['another end', ucscAnswer({ end: 20156726 })],
            ['another chromosome', ucscAnswer({ chrom: 'chr12' })],
            ['another assembly', ucscAnswer({ genome: 'hg38' })],
            ['too few bases', ucscAnswer({ dna: bases.slice(1) })],
            ['too many', ucscAnswer({ dna: `${bases}a` })],
            ['not bases', ucscAnswer({ dna: 'cagattcattctcattcttaccccTGCGTCTG-' })],
            ['no bases', ucscAnswer({ dna: undefined })],
            ['an error', { error: 'can not find specified chrom', statusCode: 400 }],
            ['nothing', null]
        ]) {
            serve({ ensembl: down, ucsc: json(answer) });
            await assert.rejects(read(), /Ensembl sequence request failed \(500\)/, why);
        }
        serve({ ensembl: down, ucsc: text('not json') });
        await assert.rejects(read(), /Ensembl sequence request failed \(500\)/);
        // (A failure is a failure, whatever comes with it.)
        serve({ ensembl: down, ucsc: json(ucscAnswer(), 503) });
        await assert.rejects(read(), /Ensembl sequence request failed \(500\)/);
        // (Bases the assembly does not know are N in both.)
        serve({ ensembl: down, ucsc: json(ucscAnswer({ dna: 'nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn' })) });
        assert.equal(await read(), 'N'.repeat(33));

        // The mitochondrial genome is UCSC's "chrMT".
        asked = serve({ ensembl: down, ucsc: json({ genome: 'hg19', chrom: 'chrMT', start: 8989, end: 8996, dna: 'CCCTGGC' }) });
        assert.equal(await h.fetchGrch37ReferenceSequence('chrMT', 8990, 8996, 5000), 'CCCTGGC');
        assert.equal(asked[1], 'https://api.genome.ucsc.edu/getData/sequence?genome=hg19&chrom=chrMT&start=8989&end=8996');
        // A sequence UCSC has no name for is not asked for there.
        asked = serve({ ensembl: down, ucsc: json(ucscAnswer()) });
        await assert.rejects(h.fetchGrch37ReferenceSequence('GL000207.1', 100, 110, 5000), /Ensembl sequence request failed \(500\)/);
        assert.equal(asked.length, 1);
        // Nothing to read: nobody is asked.
        asked = serve({ ensembl: text(bases), ucsc: json(ucscAnswer()) });
        for (const [chrom, start, end] of [['', 1, 2], ['chrX', 10, 9], ['chrX', NaN, 9], ['chrX', 1, Infinity], [null, 1, 2]]) {
            assert.equal(await h.fetchGrch37ReferenceSequence(chrom, start, end), null);
        }
        assert.deepEqual(asked, []);

        // Ensembl is given a head start, and UCSC is asked when it is over.
        h.API_TIMEOUT_MS.referenceHeadStart = 60;
        // (So that it can be told whose answer was taken, UCSC's differs here.)
        const other = 'ACGTACGTACGTACGTACGTACGTACGTACGTA';
        const ucscOther = json(ucscAnswer({ dna: other.toLowerCase() }));
        // Ensembl in good time: UCSC never hears of it.
        asked = serve({ ensembl: after(15, text(bases)), ucsc: ucscOther });
        assert.equal(await read(), bases);
        await wait(90);
        assert.deepEqual(asked, [ensemblUrl]);
        // Ensembl slow — its answers took 4 and 8 s that day: UCSC's is taken.
        asked = serve({ ensembl: after(600, text(bases)), ucsc: after(20, ucscOther) });
        began = Date.now();
        assert.equal(await read(), other);
        assert(Date.now() - began < 400, 'and Ensembl is not waited for');
        assert.deepEqual(asked, [ensemblUrl, ucscUrl]);
        // Ensembl slow, and still the first of the two.
        asked = serve({ ensembl: after(120, text(bases)), ucsc: after(600, ucscOther) });
        began = Date.now();
        assert.equal(await read(), bases);
        assert(Date.now() - began < 400);
        assert.equal(asked.length, 2);
        // Ensembl slow and then failing: UCSC's answer is waited for.
        asked = serve({ ensembl: after(100, down), ucsc: after(250, ucscOther) });
        assert.equal(await read(), other);
        // ...as Ensembl's is, when UCSC is the one that fails.
        asked = serve({ ensembl: after(250, text(bases)), ucsc: after(80, json({ error: 'x' }, 500)) });
        assert.equal(await read(), bases);
        // Neither answers: the read is given up when its time is up — the one
        // stretch of time it was given, not one for each of them.
        h.API_TIMEOUT_MS.referenceHeadStart = 400;
        asked = serve({ ensembl: never, ucsc: never });
        began = Date.now();
        await assert.rejects(read(1000), /aborted|timed out/i);
        const waited = Date.now() - began;
        assert(waited >= 950 && waited < 1300, `given up after ${waited} ms`);
        assert.equal(asked.length, 2);
    } finally {
        globalThis.fetch = realFetch;
        console.warn = realWarn;
        h.API_TIMEOUT_MS.referenceHeadStart = real.head;
    }
    console.log('reference-source unit checks passed');
}

const ucscRequests = page => page.requests.filter(r => r.startsWith('api.genome.ucsc.edu/'));
const ensemblSequenceRequests = page => page.requests.filter(r => /ensembl\.org\/sequence\/region\//.test(r));
// (Ensembl's sequence service down: every request to it fails, as on the day.)
const sequenceDown = [/ensembl\.org\/sequence\/region\//];

// 1. The row as recorded, on a day Ensembl's sequence service answered — after
//    4 to 8 seconds, or with HTTP 500. With its answers, the result needs
//    nobody else.
{
    const page = await runLookup(cassettes.eif1axDup30, {});
    assertResolved(assert, page, 'EIF1AX row');
    assert.match(page.variant, /g\.: chrX:g\.20156694_20156723dupVCF \(hg19\): X-20156693-C-CAGATTCATTCTCATTCTTACCCCTGCGTCTGene: EIF1AX/);
    assert.match(page.variant, /Transcript: NM_001412\.4 \(MANE Select\)c\.: c\.34_63dupp\.: p\.Arg12_Ser21dup \(p\.R12_S21dup\)Effect: inframe_insertion/);
    assert.match(page.panel, /Confirm allele · GRCh37 reference\s*chrX:g\.20156694_20156723dup is the same allele as the input/);
    assert.equal(ensemblSequenceRequests(page).length, 3);
    assert.deepEqual(ucscRequests(page), [], 'UCSC is not asked while Ensembl answers');
    console.log('ok reference source: with Ensembl\'s sequence service answering, nobody else is asked');
}

// 2. Ensembl's sequence service down, its recoder and VEP answering — the
//    state the row was looked up in. Was: "not confirmed — reference sequence
//    unavailable", "its other ids for this allele were not checked", and
//    "p.: N/A". UCSC has the same bases, and the result is the same result.
{
    const usual = await runLookup(cassettes.eif1axDup30, {});
    const page = await runLookup(cassettes.eif1axDup30, { without: sequenceDown, ucsc: true });
    assertResolved(assert, page, 'EIF1AX row, Ensembl\'s sequence service down');
    assert.match(page.variant, /Transcript: NM_001412\.4 \(MANE Select\)c\.: c\.34_63dupp\.: p\.Arg12_Ser21dup \(p\.R12_S21dup\)Effect: inframe_insertion/);
    assert.match(page.panel, /Protein\s*p\.Arg12_Ser21dup/);
    assert.match(page.panel, /Confirm allele · GRCh37 reference\s*chrX:g\.20156694_20156723dup is the same allele as the input/);
    assert.match(page.panel, /Annotations · MyVariant\.info\s*no record for this variant\s*0\.0s/, 'its other ids were looked for');
    assert.doesNotMatch(page.panel, /reference sequence unavailable|not confirmed|protein change not shown/);
    assert.equal(page.variant, usual.variant, 'the card is the one Ensembl\'s bases give');
    assert.equal(page.panel, usual.panel);
    // UCSC is asked for each stretch Ensembl failed to give — the first of
    // them answered as recorded from UCSC itself.
    assert.deepEqual(ucscRequests(page), [
        'api.genome.ucsc.edu/getData/sequence?genome=hg19&chrom=chrX&start=20156692&end=20156693',
        'api.genome.ucsc.edu/getData/sequence?genome=hg19&chrom=chrX&start=20156690&end=20156725',
        'api.genome.ucsc.edu/getData/sequence?genome=hg19&chrom=chrX&start=20156543&end=20156873'
    ]);
    assert.deepEqual(ensemblSequenceRequests(page).map(r => r.match(/X:\d+\.\.\d+/)[0]),
        ['X:20156693..20156693', 'X:20156691..20156725', 'X:20156544..20156873'], 'the same three stretches');
    console.log('ok reference source: with Ensembl\'s sequence service down, UCSC\'s bases give the same result — p.Arg12_Ser21dup');
}

// 3. ...and with Ensembl's recoder failing as well (it answered HTTP 500 that
//    day too): VEP names the allele as it was entered, which the reference
//    shows to be the left-aligned spelling already.
{
    const page = await runLookup(cassettes.eif1axDup30, { without: [...sequenceDown, /variant_recoder/], ucsc: true });
    assertResolved(assert, page, 'EIF1AX row, Ensembl\'s recoder and sequence service down');
    assert.match(page.variant, /g\.: chrX:g\.20156693_20156694insAGATTCATTCTCATTCTTACCCCTGCGTCT/);
    assert.match(page.variant, /c\.: c\.34_63dupp\.: p\.Arg12_Ser21dup \(p\.R12_S21dup\)Effect: inframe_insertion/);
    assert.doesNotMatch(page.panel, /reference sequence unavailable|protein change not shown/);
    console.log('ok reference source: the same with Ensembl\'s recoder down as well');
}

// 4. No reference from anyone. Then the spelling cannot be shown to be
//    left-aligned, and a protein name Ensembl gives for a spelling that may not
//    be is still not shown: that name is wrong often enough (KRAS c.34_36dup
//    came out as another variant). The rest of the result stands.
{
    const page = await runLookup(cassettes.eif1axDup30, { without: [...sequenceDown, /ucsc\.edu/] });
    assertResolved(assert, page, 'EIF1AX row, no reference to be had');
    assert.match(page.variant, /Transcript: NM_001412\.4 \(MANE Select\)c\.: c\.34_63dupp\.: N\/AEffect: inframe_insertion/);
    assert.match(page.panel, /Confirm allele · GRCh37 reference\s*not confirmed — reference sequence unavailable/);
    assert.match(page.panel, /protein change not shown: Ensembl names an in-frame insertion or deletion on the reverse strand reliably only when it is written left-aligned/);
    assert.equal(ucscRequests(page).length, 3, 'UCSC was asked');
    // (UCSC answering with something that is not the stretch asked for is no reference either.)
    const wrong = await runLookup(cassettes.eif1axDup30, { without: sequenceDown,
        route: (url, key, respond) => (/ucsc\.edu/.test(key) ? respond({ genome: 'hg19', chrom: 'chrX', start: 1, end: 2, dna: 'A' }) : undefined) });
    assert.match(wrong.variant, /c\.: c\.34_63dupp\.: N\/A/);
    assert.match(wrong.panel, /not confirmed — reference sequence unavailable/);
    console.log('ok reference source: with no reference from anyone, an unchecked protein name is still not shown');
}

// 5. Ensembl's sequence service not failing but not answering — its answers
//    took up to 8 s that day, longer than an identity check waits. UCSC is
//    asked when the head start is over, and the result does not wait for Ensembl.
{
    let held = 0;
    const page = await runLookup(cassettes.eif1axDup30, { ucsc: true, timerCap: 40, settleMs: 500,
        route: (url, key) => {
            if (!/ensembl\.org\/sequence\/region\//.test(key)) return undefined;
            held += 1;
            return new Promise(() => {});
        } });
    assertResolved(assert, page, 'EIF1AX row, Ensembl\'s sequence service silent');
    assert.match(page.variant, /c\.: c\.34_63dupp\.: p\.Arg12_Ser21dup \(p\.R12_S21dup\)Effect: inframe_insertion/);
    assert.match(page.panel, /Confirm allele · GRCh37 reference\s*chrX:g\.20156694_20156723dup is the same allele as the input/);
    assert.equal(held, 3);
    assert.equal(ucscRequests(page).length, 3);
    console.log('ok reference source: a sequence service that does not answer is not waited out');
}

// 6. Whatever reads the reference reads the same bases from UCSC: rows that
//    are left-aligned against it, alleles found under their other ids, the
//    reverse-strand names — each scenario with Ensembl's sequence service down
//    comes out as it does with it up.
{
    const same = ['brca2Neighbour', 'cftrF508del', 'tp53Dup', 'brca1Dup', 'krasDup', 'brafT599dupRight', 'brca2ExonEdgeDup', 'rs80357906', 'geneAliasRow'];
    const shown = page => ({ status: page.status, variant: page.variant, clinvar: page.clinvar, panel: page.panel });
    for (const name of same) {
        const usual = await runLookup(cassettes[name], {});
        const page = await runLookup(cassettes[name], { without: sequenceDown, ucsc: true });
        assert(ensemblSequenceRequests(usual).length > 0 && ucscRequests(usual).length === 0, `${name}: reads the reference, from Ensembl`);
        assert.equal(ucscRequests(page).length, ensemblSequenceRequests(page).length, `${name}: every stretch asked of UCSC`);
        assert.deepEqual(shown(page), shown(usual), `${name}: the same result from UCSC's bases`);
        // (Without either, the result is not the same: the scenario does depend on the reference.)
        const without = await runLookup(cassettes[name], { without: sequenceDown });
        assert.notDeepEqual(shown(without), shown(usual), `${name}: differs without the reference`);
    }
    console.log(`ok reference source: ${same.length} recorded scenarios come out the same from UCSC's bases`);
}

console.log('Reference source tests passed');
