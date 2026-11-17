/*
 * What a protein-level name comes to.
 *
 * A protein change does not say which bases changed; Ensembl's recoder answers
 * with every way it knows, and the page took the first.
 *
 *   "KRAS G12X"  X is the old code for a stop (TP53 R213X) and, in HGVS (Xaa)
 *                and in oncology shorthand, any amino acid. Read as a stop it
 *                came back as chr12:g.25398283_25398285delinsTTA,
 *                c.34_36delinsTAA, stop_gained — an allele nobody has
 *                reported, shown as the variant.
 *   "KRAS G12*"  a stop written as one: three alleles (GGT -> TAA, TAG, TGA),
 *                none on record, in an order Ensembl changes between calls.
 *                Two live runs showed two different variants.
 *
 * The scenarios replay answers recorded from the live services
 * (tests/fixtures/lookup-cassettes.json) through the real form handler; each
 * fails on the code before the fix.
 *
 * Run with: node tests/protein-name-reading.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));

// ── The pieces ─────────────────────────────────────────────────────────────
for (const id of ['chr12:g.25380277G>A', 'chrX:g.5T>C', '17:g.7578212G>A']) assert(h.isSingleBaseChange(id), id);
for (const id of ['chr12:g.25398283_25398285delinsTTA', 'chr12:g.25398283ACC>TTA', 'chr7:g.140453136_140453137delinsTA',
    'chr17:g.41209082dup', 'chr7:g.117199646_117199648del', 'chr12:g.25398283_25398285', 'KRAS:p.Gly12Ter', '']) {
    assert(!h.isSingleBaseChange(id), id);
}

// "X" typed as the new residue and read as a stop — and nothing else.
for (const [typed, query] of [['KRAS G12X', 'KRAS:p.Gly12Ter'], ['kras g12x', 'KRAS:p.Gly12Ter'], ['KRAS p.G12X', 'KRAS:p.Gly12Ter'],
    ['KRAS:p.Gly12X', 'KRAS:p.Gly12Ter'], ['KRAS p.(Gly12X)', 'KRAS:p.Gly12Ter'], ['  TP53 R213X ', 'TP53:p.Arg213Ter']]) {
    assert.equal(h.typedResidueX(typed, query), true, typed);
}
for (const [typed, query] of [['KRAS G12*', 'KRAS:p.Gly12Ter'], ['KRAS G12Ter', 'KRAS:p.Gly12Ter'], ['KRAS p.Gly12Xaa', 'KRAS:p.Gly12Xaa'],
    ['BRCA1 Q1756PfsX74', 'BRCA1:p.Q1756PfsX74'], ['KRAS G12D', 'KRAS:p.Gly12Asp'], ['X 20148674 T TG', 'chrX:g.20148674_20148675insG'],
    ['EIF1AX', 'EIF1AX']]) {
    assert.equal(h.typedResidueX(typed, query), false, typed);
}

// One allele, written twice by Ensembl, is one allele; the HGVS spelling is kept.
assert.deepEqual(h.distinctAlleleSpellings([
    'chr12:g.25398283_25398285delinsTTA', 'chr12:g.25398283ACC>TTA',
    'chr12:g.25398283ACC>CTA', 'chr12:g.25398283_25398285delinsCTA',
    'chr12:g.25398283_25398285delinsTCA'
]), ['chr12:g.25398283_25398285delinsTTA', 'chr12:g.25398283_25398285delinsCTA', 'chr12:g.25398283_25398285delinsTCA']);
assert.deepEqual(h.distinctAlleleSpellings(['chr12:g.25380277G>A', 'chr12:g.25380277G>A', 'chr12:g.25380277G>T']),
    ['chr12:g.25380277G>A', 'chr12:g.25380277G>T']);
assert.deepEqual(h.distinctAlleleSpellings([]), []);

// The rule.
const three = ['chr12:g.25398283_25398285delinsTTA', 'chr12:g.25398283_25398285delinsCTA', 'chr12:g.25398283_25398285delinsTCA'];
// A record for one of them settles it, however the name was typed.
assert.equal(h.proteinNameReading({ typedX: true, alleles: 3, candidates: three, onRecord: true }), 'allele');
assert.equal(h.proteinNameReading({ typedX: false, alleles: 3, candidates: three, onRecord: true }), 'allele');
// "G12X": a stop needs the whole codon replaced — X was not a stop.
assert.equal(h.proteinNameReading({ typedX: true, alleles: 3, candidates: three }), 'any-residue');
assert.equal(h.proteinNameReading({ typedX: true, alleles: 1, candidates: ['chr7:g.140453136_140453137delinsTA'] }), 'any-residue');
// "Q61X": one base makes the stop (CAA -> TAA), so a stop is what was meant.
assert.equal(h.proteinNameReading({ typedX: true, alleles: 1, candidates: ['chr12:g.25380277G>A'] }), 'allele');
// "G12*": a stop, but which of three?
assert.equal(h.proteinNameReading({ typedX: false, alleles: 3, candidates: three }), 'undetermined');
// "V600*": Ensembl gives one way; that is the allele, as before.
assert.equal(h.proteinNameReading({ typedX: false, alleles: 1, candidates: ['chr7:g.140453136_140453137delinsTA'] }), 'allele');
assert.equal(h.proteinNameReading(), 'allele');
// MyVariant not answering is not "nothing on record": that cannot be said, so
// the old behaviour stands — but a stop that needs the codon replaced is still
// not what "X" meant.
assert.equal(h.proteinNameReading({ typedX: false, alleles: 3, candidates: three, recordsKnown: false }), 'allele');
assert.equal(h.proteinNameReading({ typedX: true, alleles: 3, candidates: three, recordsKnown: false }), 'any-residue');
// The change is the second thing typed, whatever follows it.
assert.equal(h.typedResidueX('KRAS G12X somatic', 'KRAS:p.Gly12Ter'), true);
assert.equal(h.typedResidueX('KRAS G12* X', 'KRAS:p.Gly12Ter'), false);

// "Any amino acid" is searched the way it is written about.
assert.equal(h.proteinSearchTerm('p.Gly12Xaa'), 'G12X');
assert.equal(h.proteinSearchTerm('GLY12XAA'), 'G12X');
assert.equal(h.proteinSearchTerm('p.G12Xaa'), 'G12X');
assert.equal(h.proteinSearchTerm('p.Gly12Ter'), 'G12X');
assert.equal(h.proteinSearchTerm('p.Gly12Asp'), 'G12D');
console.log('protein-name unit checks passed');

const note = page => [...page.document.querySelectorAll('.protein-name-note')].map(el => el.textContent.replace(/\s+/g, ' ').trim());
const noteLinks = page => [...page.document.querySelectorAll('.protein-name-note a')].map(a => a.getAttribute('href'));

// 1. "KRAS G12X": the codon, read as any amino acid. It used to be
//    "g.: chr12:g.25398283_25398285delinsTTA … c.34_36delinsTAA … stop_gained".
for (const typed of ['KRAS G12X', 'kras g12x', 'KRAS p.Gly12X', 'KRAS:p.G12X', 'KRAS p.(G12X)', 'KRAS p.(Gly12X)']) {
    const page = await runLookup(cassettes.krasG12Stop, { input: typed });
    assertResolved(assert, page, typed);
    assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/, typed);
    assert.doesNotMatch(page.variant, /delins|stop_gained|Gly12Ter|Gly12\*/, typed);
    assert.match(page.variant, /p\.: p\.Gly12Xaa/, typed);
    assert.deepEqual(note(page), ['"X" was read as any amino acid (p.Gly12Xaa), not as a stop: no single-base change turns Gly12 '
        + 'into a stop codon, and no such variant is on record.'], typed);
    assert.match(page.panel, /read as KRAS:p\.Gly12Ter · in the end as any amino acid \(Xaa\), not as a stop/, typed);
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*a stop here takes more than one base change — "X" read as any amino acid/, typed);
    // Nothing is asked about a stop at Gly12, or about an allele nobody named.
    assert.deepEqual(page.api('clinvar-protein').filter(r => !/retmax=/.test(r)), [], typed);
    // ClinVar's names are consulted for it, once.
    assert.equal(page.api('clinvar-protein').filter(r => /retmax=/.test(r)).length, 1, typed);
    assert(!page.requests.some(r => /\/vep\//.test(r)), `${typed}: no consequence is asked for an allele nobody named`);
    for (const asked of page.api('pubmed')) assert.match(asked, /term=KRAS G12X&gene=KRAS&variant=G12X(&|$)/, typed);
    assert.deepEqual(page.searches('google.com/search').slice(-2),
        ['q=pathogenicity of KRAS G12X', 'q=clinical significance of KRAS G12X'], typed);
}
console.log('ok protein name: KRAS G12X is the codon, X read as any amino acid');

// 2. "KRAS G12*": a stop, but three alleles and none on record. Whichever
//    order Ensembl lists them in, none is shown as the variant; all are offered.
{
    const key = 'grch37.rest.ensembl.org/variant_recoder/human/KRAS:p.Gly12Ter?content-type=application/json';
    const recorded = cassettes.krasG12Stop.responses[key].body;
    const alleles = Object.keys(recorded[0]).filter(k => k !== 'warnings');
    assert.equal(alleles.length, 3, 'Ensembl gives three alleles for p.Gly12Ter');
    const shown = new Set();
    for (const order of [[0, 1, 2], [2, 0, 1], [1, 2, 0]]) {
        const body = [Object.fromEntries(order.map(i => [alleles[i], recorded[0][alleles[i]]]))];
        const page = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12*',
            route: (url, k, respond) => (k === key ? respond(body) : undefined) });
        assertResolved(assert, page, 'KRAS G12*');
        assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
        assert.doesNotMatch(page.variant, /g\.: chr12|c\.34_36delins|stop_gained/);
        assert.match(page.variant, /p\.: p\.Gly12Ter/);
        assert.match(note(page)[0], /^p\.Gly12Ter can be any of 3 nucleotide changes and none of them is on record, so no single one is shown: /);
        // One link per allele, each a lookup of its own; the same three whatever the order.
        assert.deepEqual(noteLinks(page).map(decodeURIComponent).sort(), [
            '?variant=chr12:g.25398283_25398285delinsCTA',
            '?variant=chr12:g.25398283_25398285delinsTCA',
            '?variant=chr12:g.25398283_25398285delinsTTA']);
        assert.match(page.panel, /Resolve nomenclature · Ensembl\s*3 different nucleotide changes give this protein change, none on record/);
        // It IS a stop that was asked about: the same-protein search stands.
        assert.deepEqual(page.api('clinvar-protein').filter(r => !/retmax=/.test(r)), ['clinvar-protein?gene=KRAS&change=Gly12Ter,G12*']);
        // (...for the stop, in the form a record of it is titled — not for "G12".)
        assert.deepEqual(page.api('clinvar-protein').filter(r => /retmax=/.test(r)), ['clinvar-protein?gene=KRAS&change=Gly12Ter&retmax=200'],
            'ClinVar\'s names are consulted once');
        shown.add(page.variant.match(/Codon region \(hg19\): (\S+?)(?=p\.|Gene|Resolved|$)/)?.[1] || page.variant.slice(0, 80));
    }
    assert.equal(shown.size, 1, 'the result does not depend on the order of Ensembl\'s answer');
    console.log('ok protein name: KRAS G12* is three alleles, none shown as the variant');
}

// 3. "KRAS Q61X": one base makes Gln61 a stop, so that is the reading — and
//    because X could have meant any amino acid, the card says which was taken.
{
    const page = await runLookup(cassettes.krasQ61X);
    assertResolved(assert, page, 'KRAS Q61X');
    assert.match(page.variant, /g\.: chr12:g\.25380277G>A/);
    assert.match(page.variant, /c\.: c\.181C>T/);
    assert.match(page.variant, /stop_gained/);
    assert.deepEqual(note(page), ['"X" was read as a stop codon (p.Gln61Ter). If it stands for any amino acid, look up KRAS p.Gln61Xaa.']);
    assert.deepEqual(noteLinks(page).map(decodeURIComponent), ['?variant=KRAS p.Gln61Xaa']);
    console.log('ok protein name: KRAS Q61X is the stop one base makes, and says so');
}

// 4. The same for the stop everyone means by it, and no note where the input
//    said "stop" itself.
{
    const typedX = await runLookup(cassettes.tp53Stop, { input: 'TP53 R213X' });
    assertResolved(assert, typedX, 'TP53 R213X');
    assert.match(typedX.variant, /g\.: chr17:g\.7578212G>A/);
    assert.deepEqual(note(typedX), ['"X" was read as a stop codon (p.Arg213Ter). If it stands for any amino acid, look up TP53 p.Arg213Xaa.']);
    // (Written as a prediction, in parentheses, it is the same name.)
    const predicted = await runLookup(cassettes.tp53Stop, { input: 'TP53 p.(R213X)' });
    assertResolved(assert, predicted, 'TP53 p.(R213X)');
    assert.match(predicted.variant, /g\.: chr17:g\.7578212G>A/);
    assert.deepEqual(note(predicted), ['"X" was read as a stop codon (p.Arg213Ter). If it stands for any amino acid, look up TP53 p.Arg213Xaa.']);
    for (const typed of ['TP53 R213*', 'TP53 p.Arg213Ter', 'TP53 p.(Arg213Ter)', 'TP53 p.(R213*)']) {
        const page = await runLookup(cassettes.tp53Stop, { input: typed });
        assertResolved(assert, page, typed);
        assert.match(page.variant, /g\.: chr17:g\.7578212G>A/, typed);
        assert.deepEqual(note(page), [], typed);
    }
    console.log('ok protein name: R213X says how X was read; R213* has nothing to say');
}

// 5. The link the note offers: "KRAS p.Gly12Xaa" is the codon, with no stop
//    asked about anywhere. (Ensembl answers 400 for Xaa; nothing is recorded
//    for it here, which the harness answers with a failure as well.)
{
    const page = await runLookup(cassettes.krasG12Stop, { input: 'KRAS p.Gly12Xaa' });
    assertResolved(assert, page, 'KRAS p.Gly12Xaa');
    assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.match(page.variant, /p\.: p\.Gly12Xaa/);
    assert.deepEqual(note(page), []);
    assert(!page.requests.some(r => /Gly12Ter/.test(r)), 'no source is asked about a stop');
    for (const asked of page.api('pubmed')) assert.match(asked, /term=KRAS G12X&gene=KRAS&variant=G12X(&|$)/);
    console.log('ok protein name: p.Gly12Xaa is the codon');
}

// 6. Ensembl not answering at all: nothing says whether one base would do, and
//    no variant is on record for the stop — X is shown for what it is.
{
    const page = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12X', without: ['variant_recoder'] });
    assertResolved(assert, page, 'KRAS G12X, recoder down');
    assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.match(page.variant, /p\.: p\.Gly12Xaa/);
    assert.deepEqual(note(page), ['"X" can mean a stop codon or any amino acid. Which nucleotide change a stop here (p.Gly12Ter) would be '
        + 'could not be looked up, so only the codon is shown, as any amino acid (p.Gly12Xaa).']);
    // Nothing is said about what is on record. For "KRAS Q61X" the same
    // sentence used to read "No variant is on record for the stop
    // (p.Gln61Ter)" — a variant MyVariant holds, and the one the page shows
    // when Ensembl does answer; its name search simply does not find it.
    // (Found in review. The codon's position is supplied here: the recording
    // of Q61X, which resolves, never needed it.)
    const q61 = cassettes.krasQ61X;
    assert.equal(q61.responses['myvariant.info/v1/variant/chr12:g.25380277G>A'].status, 200, 'MyVariant holds the stop at Gln61');
    const codon61 = { mappings: [{ assembly_name: 'GRCh37', seq_region_name: '12', start: 25380275, end: 25380277, strand: -1, rank: 0, gap: 0, coord_system: 'chromosome' }] };
    const q61WithCodon = { ...q61, responses: { ...q61.responses,
        'grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/KRAS?expand=1&content-type=application/json':
            cassettes.krasG12Stop.responses['grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/KRAS?expand=1&content-type=application/json'],
        'grch37.rest.ensembl.org/map/translation/ENSP00000256078/61..61?content-type=application/json': { status: 200, body: codon61 } } };
    const unanswered = await runLookup(q61WithCodon, { without: ['variant_recoder'] });
    assertResolved(assert, unanswered, 'KRAS Q61X, recoder down');
    assert.match(unanswered.variant, /Codon region \(hg19\): chr12:g\.25380275_25380277/);
    assert.deepEqual(note(unanswered), ['"X" can mean a stop codon or any amino acid. Which nucleotide change a stop here (p.Gln61Ter) would be '
        + 'could not be looked up, so only the codon is shown, as any amino acid (p.Gln61Xaa).']);
    assert.doesNotMatch(unanswered.variant, /on record/);
    // Where one base does make the stop but Ensembl gives a second way too, and
    // MyVariant holds neither: X was a stop, and it is the stop that has two
    // alleles — not "p.Gln61Xaa can be any of 2 nucleotide changes".
    const recoderKey = Object.keys(q61.responses).find(key => key.includes('/variant_recoder/'));
    const one = Object.values(q61.responses[recoderKey].body[0]).find(value => value && value.hgvsg);
    const two = [{ T: one, TT: { ...one, hgvsg: ['NC_000012.11:g.25380276_25380277delinsAA'], spdi: ['NC_000012.11:25380275:TG:AA'] } }];
    const twoWays = await runLookup(q61WithCodon, {
        route: (url, key, respond) => {
            if (key === recoderKey) return respond(two);
            if (/myvariant\.info\/v1\/variant\/chr12:g\.2538027/.test(key)) return respond({ success: false, error: 'Not Found.' }, 404);
            return undefined;
        }
    });
    assertResolved(assert, twoWays, 'KRAS Q61X, two alleles, none on record');
    assert.match(twoWays.variant, /p\.: p\.Gln61Ter/);
    assert.doesNotMatch(twoWays.variant, /Xaa/);
    assert.match(note(twoWays)[0], /^"X" was read as a stop codon\. p\.Gln61Ter can be any of 2 nucleotide changes and none of them is on record, so no single one is shown: /);
    console.log('ok protein name: G12X with the recoder down is still not shown as a stop');
}

// 7. MyVariant down. A 404 says "no record"; an outage says nothing, so
//    "none of them is on record" is not claimed — while a stop that takes the
//    whole codon is still not what X meant, record or no record.
{
    const down = (url, key, respond) => (url.hostname === 'myvariant.info' ? respond({ error: 'Service Unavailable' }, 503) : undefined);
    const typedX = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12X', route: down });
    assertResolved(assert, typedX, 'KRAS G12X, MyVariant down');
    assert.match(typedX.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.deepEqual(note(typedX), ['"X" was read as any amino acid (p.Gly12Xaa), not as a stop: no single-base change turns Gly12 into a stop codon.']);
    const stop = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12*', route: down });
    assert.doesNotMatch(note(stop).join(' '), /on record/, 'nothing is said about records MyVariant was not asked for');
    console.log('ok protein name: a MyVariant outage is not "nothing on record"');
}

// 8. A missense with a record is what it always was.
{
    const page = await runLookup(cassettes.brafV600E);
    assertResolved(assert, page, 'BRAF V600E');
    assert.match(page.variant, /g\.: chr7:g\.140453136A>T/);
    assert.deepEqual(note(page), []);
    console.log('ok protein name: BRAF V600E unchanged');
}

// 9. Several alleles ON record. "KRAS Y64*" is TAC -> TAG or TAA, MyVariant
//    holds both, and its name search finds neither — so which was shown was
//    again the order of Ensembl's answer (c.192C>G in one call, c.192C>A in
//    another, live). The best documented one is shown, whatever the order, and
//    the card names the other. (Found in review.)
{
    const rec = (id, extra) => ({ id, record: { _id: id, ...extra } });
    const rcv = n => ({ rcv: Array.from({ length: n }, (_, i) => ({ accession: `RCV${i}` })) });
    // ClinVar submissions, then frequency, then whether dbSNP has seen it...
    let out = h.rankAllelesOnRecord([rec('chr1:g.5A>C', {}), rec('chr1:g.5A>G', { clinvar: rcv(2) }), rec('chr1:g.5A>T', { clinvar: rcv(5) })]);
    assert.deepEqual(out.ranked.map(e => e.id), ['chr1:g.5A>T', 'chr1:g.5A>G', 'chr1:g.5A>C']);
    assert.equal(out.tied, false);
    out = h.rankAllelesOnRecord([rec('chr1:g.5A>C', { gnomad_exome: { af: { af: 0.001 } } }), rec('chr1:g.5A>G', { gnomad_genome: { af: { af: 0.2 } } })]);
    assert.deepEqual(out.ranked.map(e => e.id), ['chr1:g.5A>G', 'chr1:g.5A>C']);
    out = h.rankAllelesOnRecord([rec('chr1:g.5A>T', {}), rec('chr1:g.5A>C', { dbsnp: { rsid: 'rs1' } })]);
    assert.deepEqual(out.ranked.map(e => e.id), ['chr1:g.5A>C', 'chr1:g.5A>T']);
    assert.equal(out.tied, false);
    // ...and where nothing tells them apart, the same one every time, and it says so.
    for (const order of [['chr1:g.5A>T', 'chr1:g.5A>C'], ['chr1:g.5A>C', 'chr1:g.5A>T']]) {
        out = h.rankAllelesOnRecord(order.map(id => rec(id, {})));
        assert.deepEqual(out.ranked.map(e => e.id), ['chr1:g.5A>C', 'chr1:g.5A>T']);
        assert.equal(out.tied, true);
    }
    assert.deepEqual(h.rankAllelesOnRecord([]), { ranked: [], tied: false });

    const cassette = cassettes.krasY64Stop;
    const key = 'grch37.rest.ensembl.org/variant_recoder/human/KRAS:p.Tyr64Ter?content-type=application/json';
    const recorded = cassette.responses[key].body;
    const alleles = Object.keys(recorded[0]).filter(k => k !== 'warnings');
    assert.equal(alleles.length, 2, 'Ensembl gives two alleles for p.Tyr64Ter');
    for (const id of ['chr12:g.25380266G>C', 'chr12:g.25380266G>T']) {
        assert.equal(cassette.responses[`myvariant.info/v1/variant/${id}`].status, 200, `MyVariant holds ${id}`);
    }
    // (Of the two, dbSNP has seen G>C; neither is in ClinVar or gnomAD.)
    for (const typed of ['KRAS Y64*', 'KRAS Y64X']) {
        for (const order of [alleles, [...alleles].reverse()]) {
            const body = [Object.fromEntries(order.map(a => [a, recorded[0][a]]))];
            const page = await runLookup(cassette, { input: typed, route: (url, k, respond) => (k === key ? respond(body) : undefined) });
            assertResolved(assert, page, `${typed}, alleles listed ${order.join(',')}`);
            assert.match(page.variant, /g\.: chr12:g\.25380266G>C/, `${typed} ${order}`);
            assert.match(page.variant, /c\.: c\.192C>G/, `${typed} ${order}`);
            const alleleNote = page.document.querySelector('.protein-alleles-note');
            assert.equal(alleleNote.textContent, 'p.Tyr64Ter is 2 alleles on record. This is the best documented one; the other is chr12:g.25380266G>T.', `${typed} ${order}`);
            assert.deepEqual([...alleleNote.querySelectorAll('a')].map(a => decodeURIComponent(a.getAttribute('href'))), ['?variant=chr12:g.25380266G>T']);
            assert.match(page.panel, /2 alleles on record for this protein change — showing the best documented, chr12:g\.25380266G>C/, `${typed} ${order}`);
            // Nothing is asked about the allele that is not shown.
            assert(!page.requests.some(r => /vep\/human\/hgvs\/12:g\.25380266G>T/.test(r)), `${typed} ${order}`);
            assert.deepEqual(note(page), typed.endsWith('X')
                ? ['"X" was read as a stop codon (p.Tyr64Ter). If it stands for any amino acid, look up KRAS p.Tyr64Xaa.'] : [], `${typed} ${order}`);
        }
    }
    // One allele on record is simply the variant, as before: no such note.
    assert.equal((await runLookup(cassettes.krasQ61X)).document.querySelector('.protein-alleles-note'), null);
    console.log('ok protein name: of two alleles on record the best documented is shown, whatever Ensembl\'s order');
}

console.log('Protein name reading tests passed');
