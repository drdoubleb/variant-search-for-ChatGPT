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
// MyVariant not answering changes what can be SAID — "none on record" cannot —
// and not what is shown: there is still no record to choose one of three by.
// (The first one Ensembl listed used to be shown then.)
assert.equal(h.proteinNameReading({ typedX: false, alleles: 3, candidates: three, onRecord: false }), 'undetermined');
assert.equal(h.proteinNameReading({ typedX: true, alleles: 3, candidates: three, onRecord: false }), 'any-residue');
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
for (const typed of ['KRAS G12X', 'kras g12x', 'KRAS p.Gly12X', 'KRAS:p.G12X', 'KRAS p.(G12X)', 'KRAS p.(Gly12X)', 'KRAS:p.(G12X)']) {
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
        // ...and listed in the same order: the sentence does not change from
        // one call to the next either.
        assert.deepEqual(noteLinks(page).map(decodeURIComponent), [
            '?variant=chr12:g.25398283_25398285delinsCTA',
            '?variant=chr12:g.25398283_25398285delinsTCA',
            '?variant=chr12:g.25398283_25398285delinsTTA'], String(order));
        assert.equal(note(page)[0], 'p.Gly12Ter can be any of 3 nucleotide changes and none of them is on record, so no single one is shown: '
            + 'chr12:g.25398283_25398285delinsCTA, chr12:g.25398283_25398285delinsTCA, chr12:g.25398283_25398285delinsTTA.', String(order));
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

// 2b. One allele, and no record of it: that allele is what the name means —
//     but ClinVar's names are asked first, because Ensembl's allele is its way
//     of producing the protein change and the variant known by the name can be
//     another one. (A substitution is not asked about ahead of Ensembl; this is
//     where it is asked.)
{
    const key = 'grch37.rest.ensembl.org/variant_recoder/human/KRAS:p.Gly12Ter?content-type=application/json';
    const recorded = cassettes.krasG12Stop.responses[key].body;
    const [first] = Object.keys(recorded[0]).filter(k => k !== 'warnings');
    const page = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12*',
        route: (url, k, respond) => (k === key ? respond([{ [first]: recorded[0][first] }]) : undefined) });
    assertResolved(assert, page, 'KRAS G12* (one allele)');
    assert.match(page.variant, /g\.: chr12:g\.25398283_25398285delins[ACGT]{3}/);
    assert.doesNotMatch(page.variant, /Codon region/);
    assert.deepEqual(note(page), []);
    assert.deepEqual(page.api('clinvar-protein').filter(r => /retmax=/.test(r)), ['clinvar-protein?gene=KRAS&change=Gly12Ter&retmax=200'],
        'ClinVar\'s names are consulted before the allele is shown');
    assert.match(page.panel, /Other names · ClinVar\s*no KRAS variant listed under "G12\*"/);
    console.log('ok protein name: a single unrecorded allele is shown, after ClinVar\'s names were asked');
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
    for (const typed of ['TP53 R213*', 'TP53 p.Arg213Ter', 'TP53 p.(Arg213Ter)', 'TP53 p.(R213*)', 'TP53:p.(Arg213Ter)']) {
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

// 6. Neither source of the alleles answering — Ensembl, nor MyVariant's list
//    of the single-base changes at the residue (which settles these names on
//    its own: tests/protein-substitution-alleles.test.js). Nothing then says
//    whether one base would do, and X is shown for what it is.
{
    const neither = ['variant_recoder', 'dbnsfp.aa.pos'];
    const page = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12X', without: neither });
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
    const unanswered = await runLookup(q61WithCodon, { without: neither });
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
        without: ['dbnsfp.aa.pos'],
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
    // A stop written as one, three ways to make it, and no telling which is on
    // record: none is shown as the variant. The first one Ensembl listed used
    // to be — a different allele from one run to the next.
    const key = 'grch37.rest.ensembl.org/variant_recoder/human/KRAS:p.Gly12Ter?content-type=application/json';
    const recorded = cassettes.krasG12Stop.responses[key].body;
    const alleles = Object.keys(recorded[0]).filter(k => k !== 'warnings');
    for (const order of [[0, 1, 2], [2, 0, 1]]) {
        const body = [Object.fromEntries(order.map(i => [alleles[i], recorded[0][alleles[i]]]))];
        const stop = await runLookup(cassettes.krasG12Stop, { input: 'KRAS G12*',
            route: (url, k, respond) => (k === key ? respond(body) : down(url, k, respond)) });
        assertResolved(assert, stop, 'KRAS G12*, MyVariant down');
        assert.match(stop.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/, String(order));
        assert.doesNotMatch(stop.variant, /g\.: chr12|c\.34_36delins|stop_gained/, String(order));
        assert.equal(note(stop)[0].replace(/: chr12.*$/, ':'), 'p.Gly12Ter can be any of 3 nucleotide changes. MyVariant could not be asked which '
            + 'of them is on record, so no single one is shown:', String(order));
        assert.deepEqual(noteLinks(stop).map(decodeURIComponent).sort(), [
            '?variant=chr12:g.25398283_25398285delinsCTA',
            '?variant=chr12:g.25398283_25398285delinsTCA',
            '?variant=chr12:g.25398283_25398285delinsTTA'], String(order));
        assert.doesNotMatch(note(stop).join(' '), /none of them is on record/, 'nothing is said about records MyVariant was not asked for');
        assert.match(stop.panel, /3 different nucleotide changes give this protein change, and MyVariant could not be asked which is on record/);
    }
    console.log('ok protein name: with MyVariant down, several alleles are the codon and "none on record" is not said');
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
    // (The one dbSNP has seen, though its id sorts last.)
    out = h.rankAllelesOnRecord([rec('chr1:g.5A>C', {}), rec('chr1:g.5A>T', { dbsnp: { rsid: 'rs1' } })]);
    assert.deepEqual(out.ranked.map(e => e.id), ['chr1:g.5A>T', 'chr1:g.5A>C']);
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
            // (Ensembl's answer being all there is: MyVariant's own list of the
            // residue's single-base changes, which has no order to depend on, is
            // withheld here and tested in protein-substitution-alleles.test.js.)
            const page = await runLookup(cassette, { input: typed, without: ['dbnsfp.aa.pos'],
                route: (url, k, respond) => (k === key ? respond(body) : undefined) });
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
    assert.equal((await runLookup(cassettes.krasQ61X)).document.querySelectorAll('.protein-alleles-note').length, 0);
    console.log('ok protein name: of two alleles on record the best documented is shown, whatever Ensembl\'s order');
}

// 10. "X" typed, and ClinVar lists several variants with the stop's protein
//     change. Then a stop is what the name was: the codon is shown as the
//     stop, with ClinVar's variants, and is not turned into "any amino acid"
//     on its way there. (ClinVar holds no stop at KRAS Gly12; the two records
//     are made up, in the shape of its answer.)
{
    const answer = { variants: [
        { id: '900001', title: 'NM_004985.5(KRAS):c.34_35delinsTA (p.Gly12Ter)', aliases: [], variationCount: 1 },
        { id: '900002', title: 'NM_004985.5(KRAS):c.34_36delinsTAG (p.Gly12Ter)', aliases: [], variationCount: 1 }], total: 2, notFound: [] };
    const shared = page => [...page.document.querySelectorAll('.protein-change-shared-note')].map(el => el.textContent.replace(/\s+/g, ' ').trim());
    for (const typed of ['KRAS G12X', 'KRAS G12*']) {
        const page = await runLookup(cassettes.krasG12Stop, { input: typed, backend: { 'clinvar-protein': () => answer } });
        assertResolved(assert, page, typed);
        assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/, typed);
        assert.match(page.variant, /p\.: p\.Gly12Ter/, typed);
        assert.doesNotMatch(page.variant + page.panel, /Xaa\b(?! *\.)|any amino acid \(/, `${typed}: not re-read as any amino acid`);
        assert.deepEqual(shared(page), [`2 KRAS variants in ClinVar have the protein change "${typed.split(' ')[1]}", and the name does not say which: `
            + 'NM_004985.5:c.34_35delinsTA, NM_004985.5:c.34_36delinsTAG.'], typed);
        // How X was read is said where it was typed, with the other reading as a link.
        assert.deepEqual(note(page), typed.endsWith('X')
            ? ['"X" was read as a stop codon (p.Gly12Ter). If it stands for any amino acid, look up KRAS p.Gly12Xaa.'] : [], typed);
        // The stop is what ClinVar is asked about.
        assert.deepEqual(page.api('clinvar-protein').filter(r => !/retmax=/.test(r)), ['clinvar-protein?gene=KRAS&change=Gly12Ter,G12*'], typed);
        assert.match(page.panel, /Resolve nomenclature · Ensembl\s*its allele is one of several ClinVar has for this name — not shown as the variant/, typed);
    }
    console.log('ok protein name: X typed where ClinVar lists stop variants stays a stop');
}

// 11. The same name on a protein accession. All of the above went through the
//     gene symbol, which is how the codon is found — so "NP_004976.2:p.Gly12Ter"
//     still took the first of Ensembl's three alleles. The codon is found on
//     the accession's own protein instead: Ensembl's answer names the change on
//     each protein it is the same change on, its own among them.
{
    const answer = [
        { input: 'NP_004976.2:p.Gly12Ter', hgvsg: ['NC_000012.11:g.25398283_25398285delinsCTA'],
            hgvsp: ['ENSP00000256078.4:p.Gly12Ter', 'ENSP00000308495.3:p.Gly12Ter', 'NP_004976.2:p.Gly12Ter'] }
    ];
    assert.equal(h.translationForProteinName(answer, 'NP_004976.2:p.Gly12Ter'), 'ENSP00000256078');
    assert.equal(h.translationForProteinName(answer, 'NP_004976.2:p.(Gly12*)'), 'ENSP00000256078', 'however the change is written');
    assert.equal(h.translationForProteinName(answer, 'NP_004976.2:p.Gly13Ter'), '', 'another residue: none of them');
    assert.equal(h.translationForProteinName([{ hgvsp: ['NP_004976.2:p.Gly12Ter'] }], 'NP_004976.2:p.Gly12Ter'), '', 'no Ensembl protein in the answer');
    assert.equal(h.translationForProteinName([], 'ENSP00000256078.4:p.Gly12Ter'), 'ENSP00000256078', 'an Ensembl accession is its own');
    assert.equal(h.translationForProteinName([{ hgvsp: ['ENSP00000275493.2:p.Glu746_Ala750del'] }], 'NP_005219.2:p.Glu746_Ala750del'),
        'ENSP00000275493', 'a change that is not a substitution, by its text');
    assert.equal(h.translationForProteinName(null, ''), '');

    const cassette = cassettes.npGly12Ter;
    const key = 'grch37.rest.ensembl.org/variant_recoder/human/NP_004976.2:p.Gly12Ter?content-type=application/json';
    const recorded = cassette.responses[key].body;
    const alleles = Object.keys(recorded[0]).filter(k => k !== 'warnings');
    assert.equal(alleles.length, 3, 'Ensembl gives three alleles for NP_004976.2:p.Gly12Ter');
    for (const order of [[0, 1, 2], [2, 0, 1], [1, 2, 0]]) {
        const body = [Object.fromEntries(order.map(i => [alleles[i], recorded[0][alleles[i]]]))];
        const page = await runLookup(cassette, { route: (url, k, respond) => (k === key ? respond(body) : undefined) });
        assertResolved(assert, page, 'NP_004976.2:p.Gly12Ter');
        assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/, String(order));
        assert.doesNotMatch(page.variant, /g\.: chr12|c\.34_36delins|stop_gained/, String(order));
        assert.match(page.variant, /Gene: KRAS/, String(order));
        assert.match(page.variant, /p\.: p\.Gly12Ter/, String(order));
        assert.match(note(page)[0], /^p\.Gly12Ter can be any of 3 nucleotide changes and none of them is on record, so no single one is shown: /);
        assert.deepEqual(noteLinks(page).map(decodeURIComponent).sort(), [
            '?variant=chr12:g.25398283_25398285delinsCTA',
            '?variant=chr12:g.25398283_25398285delinsTCA',
            '?variant=chr12:g.25398283_25398285delinsTTA'], String(order));
        // On the protein that was named, found through Ensembl's own id for it.
        assert.match(page.variant, /Resolved to residue 12 on NP_004976\.2 only\./, String(order));
        assert(page.requests.some(r => /map\/translation\/ENSP00000256078\/12\.\.12/.test(r)), 'the codon is mapped on Ensembl\'s protein for it');
        assert(page.requests.some(r => /datasets\/v2\/gene\/accession\/NP_004976\.2$/.test(r)), 'and the gene is NCBI\'s for the accession');
        // The map can read the archive after resolution. Check the lookup used to resolve the codon.
        assert(!page.requests.some(r => r.startsWith('grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/KRAS?expand=1')),
            'codon resolution skips the gene canonical lookup');
    }
    // MyVariant not answering: the codon just the same, and nothing said of records.
    const down = (url, k, respond) => (url.hostname === 'myvariant.info' ? respond({ error: 'Service Unavailable' }, 503) : undefined);
    const unasked = await runLookup(cassette, { route: down });
    assertResolved(assert, unasked, 'NP_004976.2:p.Gly12Ter, MyVariant down');
    assert.match(unasked.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.match(note(unasked)[0], /^p\.Gly12Ter can be any of 3 nucleotide changes\. MyVariant could not be asked which of them is on record/);
    // NCBI not saying whose accession it is: the codon is still the codon.
    const nameless = await runLookup(cassette, { without: ['gene/accession'] });
    assertResolved(assert, nameless, 'NP_004976.2:p.Gly12Ter, gene not learned');
    assert.match(nameless.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/);
    assert.equal(noteLinks(nameless).length, 3);

    // Two alleles ON record under an accession: the best documented, whatever
    // order Ensembl lists them in. (The first listed used to be shown here as
    // well: no record "fitted" the name, because outside ClinVar MyVariant's
    // protein notations carry no accession. It is Ensembl's answer that says
    // each allele is p.Tyr64Ter on NP_004976.2.)
    const named = { input: 'NP_004976.2:p.Tyr64Ter', hgvsp: ['ENSP00000256078.4:p.Tyr64Ter', 'NP_004976.2:p.Tyr64Ter', 'NP_203524.1:p.Tyr64Ter'] };
    assert.equal(h.recoderAlleleNamesProtein(named, 'NP_004976.2:p.Tyr64Ter'), true);
    assert.equal(h.recoderAlleleNamesProtein(named, 'NP_004976.2:p.(Y64*)'), true);
    assert.equal(h.recoderAlleleNamesProtein(named, 'NP_004976.1:p.Tyr64Ter'), false, 'another version of the accession');
    assert.equal(h.recoderAlleleNamesProtein(named, 'NP_004976.2:p.Tyr64Cys'), false, 'another change');
    assert.equal(h.recoderAlleleNamesProtein(named, 'KRAS:p.Tyr64Ter'), false, 'a gene symbol is matched by the record, not by this');
    assert.equal(h.recoderAlleleNamesProtein(null, 'NP_004976.2:p.Tyr64Ter'), false);
    const twoKey = 'grch37.rest.ensembl.org/variant_recoder/human/NP_004976.2:p.Tyr64Ter?content-type=application/json';
    const twoRecorded = cassettes.npTyr64Ter.responses[twoKey].body;
    const twoAlleles = Object.keys(twoRecorded[0]).filter(k => k !== 'warnings');
    assert.equal(twoAlleles.length, 2);
    for (const order of [twoAlleles, [...twoAlleles].reverse()]) {
        const body = [Object.fromEntries(order.map(a => [a, twoRecorded[0][a]]))];
        const two = await runLookup(cassettes.npTyr64Ter, { route: (url, k, respond) => (k === twoKey ? respond(body) : undefined) });
        assertResolved(assert, two, `NP_004976.2:p.Tyr64Ter, alleles listed ${order}`);
        assert.match(two.variant, /g\.: chr12:g\.25380266G>C/, String(order));
        assert.deepEqual([...two.document.querySelectorAll('.protein-alleles-note')].map(el => el.textContent),
            ['p.Tyr64Ter is 2 alleles on record. This is the best documented one; the other is chr12:g.25380266G>T.'], String(order));
    }
    // The same on Ensembl's own id for the protein — which, written without
    // its version, reads like a gene symbol. It is not one: the card said
    // "Gene: ENSP00000256078", and ClinVar's names, CIViC, PubMed and the MANE
    // lookup were asked about a gene of that name. (Found in review. Ensembl's
    // answers below are the recorded one for the RefSeq accession, and what it
    // says of its protein and transcript.)
    for (const typed of ['ENSP00000256078:p.Gly12Ter', 'ENSP00000256078.4:p.Gly12Ter', 'ensp00000256078:p.G12*']) {
        const page = await runLookup(cassette, { input: typed,
            route: (url, k, respond) => {
                if (/\/variant_recoder\/human\/ENSP00000256078(?:\.4)?:p\.Gly12Ter/i.test(k)) return respond(recorded);
                if (/lookup\/id\/ENSP00000256078\?/.test(k)) return respond({ id: 'ENSP00000256078', Parent: 'ENST00000256078', object_type: 'Translation' });
                if (/lookup\/id\/ENST00000256078\?content-type/.test(k)) return respond({ id: 'ENST00000256078', display_name: 'KRAS-004', Parent: 'ENSG00000133703' });
                return undefined;
            } });
        assertResolved(assert, page, typed);
        assert.match(page.variant, /Codon region \(hg19\): chr12:g\.25398283_25398285/, typed);
        assert.match(page.variant, /Gene: KRAS/, typed);
        assert.doesNotMatch(page.variant, /Gene: ENSP/i, typed);
        assert.match(page.variant, /Resolved to residue 12 on ENSP00000256078(?:\.4)? only\./i, typed);
        assert.deepEqual(page.requests.filter(r => /ENSP00000256078/i.test(r) && !/variant_recoder|map\/translation|myvariant\.info|lookup\/id\/ENSP/.test(r)), [],
            `${typed}: nothing else is asked about the accession, as if it were a gene`);
        assert(page.requests.some(r => /[?&]gene=KRAS\b/.test(r)), `${typed}: gene-level cards ask about KRAS`);
    }
    console.log('ok protein name: on a protein accession, several unrecorded alleles are the codon too');
}

console.log('Protein name reading tests passed');
