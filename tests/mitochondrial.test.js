/*
 * Mitochondrial variants.
 *
 * m.8993T>G (MT-ATP6 p.Leu156Arg — NARP / Leigh syndrome) was presented as
 *
 *     Gene: MT-ATP6, MT-ATP8, RNR2-END   Transcript: RNR2-END
 *     c.: n.8993T>G   p.: N/A   Effect: intergenic_region
 *
 * because every source the page read first is wrong about the mitochondrial
 * genome in its own way:
 *
 *   - snpEff (hg19) has no mitochondrial gene models beyond the two rRNAs and
 *     calls everything else intergenic, against a feature named "RNR2-END";
 *   - the one gene MyVariant keeps from ClinVar is a neighbour (MT-ATP8), or for
 *     m.1555A>G not near the variant at all (MT-ND1);
 *   - Ensembl's RefSeq transcript set has no accessions for these genes and
 *     makes up its own ("ATP6", "ATP6.1");
 *   - Ensembl answers "NC_012920.1:m.8993T>G", which was not recognised as the
 *     allele asked about, so its transcripts were dropped and VEP never asked;
 *   - and Ensembl cannot parse "chrMT:g.…", MyVariant's spelling.
 *
 * The scenarios replay answers recorded from the live services.
 *
 * Run with: node tests/mitochondrial.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const recordOf = (cassette, id) => cassette.responses[`myvariant.info/v1/variant/${id}`].body;
const recoderOf = cassette => Object.entries(cassette.responses).find(([key]) => key.includes('/variant_recoder/'))[1].body;

// ── One sequence, however it is named ──────────────────────────────────────
assert.equal(h.isMitochondrialAccession('NC_012920.1'), true);
assert.equal(h.isMitochondrialAccession('NC_012920.1:8992:T:G'), true);
assert.equal(h.isMitochondrialAccession('NC_000017.10'), false);
assert.equal(h.isMitochondrialAccession('NC_001807.4'), false, 'the old hg19 chrM is a different sequence');
assert.equal(h.mitochondrialToChr('NC_012920.1:m.8993T>G'), 'chrMT:g.8993T>G');
assert.equal(h.mitochondrialToChr('NC_012920.1:g.8993T>G'), 'chrMT:g.8993T>G');
assert.equal(h.mitochondrialToChr('NC_000017.10:g.41209082dup'), null);
for (const spelling of ['chrMT:g.8993T>G', 'chrM:g.8993T>G', 'MT:g.8993T>G', 'M:g.8993T>G', 'NC_012920.1:m.8993T>G']) {
    assert.equal(h.normalizedGenomicId(spelling), 'chrmt:g.8993t>g', spelling);
}
// A nuclear accession still declares its build.
assert.equal(h.normalizedGenomicId('NC_000007.13:g.140453136A>T'), 'chr7:g.140453136a>t');
assert.equal(h.normalizedGenomicId('NC_000007.14:g.140753336A>T'), '', 'GRCh38 is never compared as hg19');

// ── Ensembl's answer belongs to the allele asked about ─────────────────────
{
    // Asked "MT:g.8993T>G", it echoes that and answers in "m." notation.
    const recoder = recoderOf(cassettes.mt8993);
    assert.equal(h.selectRecoderAlleles(recoder, 'chrMT:g.8993T>G').length, 1, 'the echoed "MT:" input is the query');
    const rows = h.recoderTranscriptsForGenomic(recoder, 'chrMT:g.8993T>G', 'chrMT:g.8993T>G');
    assert.deepEqual(rows.map(t => `${t.transcript}:${t.cDNA}`), ['ENST00000361899.2:c.467T>G'],
        'the real transcript, not the recoder\'s own "ATP6.1"');
}
{
    // A deletion in a run of T's: the record is filed as g.3271del, Ensembl
    // answers g.3273del. Same allele; its rows were being discarded.
    const recoder = recoderOf(cassettes.mt3271del);
    assert.deepEqual(Object.values(recoder[0])[0].hgvsg, ['NC_012920.1:m.3273del']);
    const rows = h.recoderTranscriptsForGenomic(recoder, 'chrMT:g.3271del', 'chrMT:g.3271del');
    assert.deepEqual(rows.map(t => `${t.transcript}:${t.cDNA}`), ['ENST00000386347.1:n.44del']);
    // Only a genomic query has a single, known allele. For a name the id decides.
    assert.deepEqual(h.recoderTranscriptsForGenomic(recoder, 'chrMT:g.3271del', 'MT-TL1:n.44del'), []);
    assert.deepEqual(h.recoderTranscriptsForGenomic(recoder, 'chrMT:g.3271del'), []);
}

// ── What the record says about the gene ────────────────────────────────────
{
    const raw = recordOf(cassettes.mt8993, 'chrMT:g.8993T>G');
    assert.equal(raw.snpeff.ann.effect, 'intergenic_region', 'as MyVariant serves it');
    assert.equal(raw.clinvar.gene.symbol, 'MT-ATP8', 'as MyVariant serves it');
    const record = h.withoutUnusableBlocks(raw);
    assert.equal(record.snpeff, undefined);
    assert.match(record._not_used.snpeff, /no mitochondrial gene models/);
    assert.equal(record.clinvar, raw.clinvar, 'nothing else is touched');
    const geneRow = rec => h.buildSummary(rec, rec._id).find(row => row.name === 'Gene(s)').value;
    assert.equal(geneRow(record), 'MT-ATP6');
    // m.1555A>G is in MT-RNR1. MyVariant's ClinVar block says MT-ND1 twice over.
    const rnr1 = recordOf(cassettes.mt1555, 'chrMT:g.1555A>G');
    assert.equal(rnr1.clinvar.gene.symbol, 'MT-ND1');
    assert.equal(geneRow(h.withoutUnusableBlocks(rnr1)), 'MT-RNR1');
    // No CADD row for an indel: ClinVar's name is then the only one there is.
    assert.equal(geneRow(h.withoutUnusableBlocks(recordOf(cassettes.mt3271del, 'chrMT:g.3271del'))), 'MT-TL1');
    // A nuclear record keeps its snpEff block, and both of ClinVar's labels.
    const braf = recordOf(cassettes.brafV600E, 'chr7:g.140453136A>T');
    assert.equal(h.withoutUnusableBlocks(braf), braf);
    assert(h.getHitGenesUpper(raw).has('MT-ATP6'), 'CADD\'s gene and ClinVar\'s named gene both count as labels');
}
console.log('Mitochondrial unit checks passed');

// ── Full lookups ───────────────────────────────────────────────────────────
const asked = (page, pattern) => page.requests.filter(r => pattern.test(r));

// 1. The variant, in every way it gets written. The first is how it comes out
//    of an hg19 VCF; the fourth is ClinVar's own name for it.
for (const typed of [
    'M 8993 T G', 'chrM\t8993\tT\tG', 'MT 8993 T G', 'NC_012920.1:m.8993T>G', 'm.8993T>G', 'm.8993t>g',
    'chrM:g.8993T>G', 'chrMT:g.8993T>G', 'MT:g.8993T>G', 'MT-ATP6 m.8993T>G', 'mt-atp6:m.8993T>G'
]) {
    const page = await runLookup(cassettes.mt8993, { input: typed });
    assertResolved(assert, page, typed);
    assert.match(page.panel, typed === 'chrMT:g.8993T>G' ? /HGVS input/ : /read as chrMT:g\.8993T>G/, typed);
    assert.match(page.variant, /g\.: chrMT:g\.8993T>G/, typed);
    assert.match(page.variant, /Gene: MT-ATP6Transcript: ENST00000361899\.2/, `${typed}: one gene, the right one`);
    assert.match(page.variant, /c\.: c\.467T>Gp\.: p\.Leu156Arg \(p\.L156R\)Effect: missense_variant/, typed);
    assert.doesNotMatch(page.variant + page.panel, /RNR2|intergenic|MT-ATP8|ATP6\.1/, typed);
    assert.match(page.clinvar, /Variation ID: 9641/, typed);
    // Ensembl is asked in its own spelling, for its own transcript set, once.
    assert.deepEqual(asked(page, /\/variant_recoder\//).map(r => r.match(/human\/([^?]+)/)[1]), ['MT:g.8993T>G'], typed);
    assert.deepEqual(asked(page, /\/vep\/human\/hgvs\//).map(r => r.match(/hgvs\/(.+)$/)[1]),
        ['MT:g.8993T>G?content-type=application/json&hgvs=1&mane=1&canonical=1'], `${typed}: no refseq=1`);
    assert.deepEqual(asked(page, /\/map\/human\//), [], `${typed}: nothing to lift`);
    assert.deepEqual(asked(page, /myvariant\.info/), ['myvariant.info/v1/variant/chrMT:g.8993T>G'], `${typed}: found at once`);
    // Gene-level cards are asked about the gene, and only that.
    const genes = [...new Set(page.requests.map(r => (r.match(/[?&]gene=([^&]+)/) || [])[1]).filter(Boolean))];
    assert.deepEqual(genes, ['MT-ATP6'], typed);
    assert(!page.requests.some(r => /RNR2-END/.test(r)), typed);
}
console.log('ok mitochondrial lookup: m.8993T>G is MT-ATP6 p.Leu156Arg, in eleven spellings');

// 2. An rRNA variant. Was: "Gene: MT-RNR1, MT-ND1, RNR2, RNR1".
{
    const page = await runLookup(cassettes.mt1555);
    assertResolved(assert, page, 'm.1555A>G');
    assert.match(page.variant, /Gene: MT-RNR1Transcript: ENST00000389680\.2/);
    assert.match(page.variant, /c\.: n\.908A>Gp\.: N\/AEffect: non_coding_transcript_exon_variant/);
    assert.doesNotMatch(page.variant, /MT-ND1|RNR2|n\.906/);
    console.log('ok mitochondrial lookup: m.1555A>G is in MT-RNR1');
}

// 3. A tRNA variant, typed as ClinVar names it.
{
    const page = await runLookup(cassettes.mt3243);
    assertResolved(assert, page, 'm.3243A>G');
    assert.match(page.variant, /g\.: chrMT:g\.3243A>G/);
    assert.match(page.variant, /Gene: MT-TL1Transcript: ENST00000386347\.1/);
    assert.match(page.variant, /c\.: n\.14A>Gp\.: N\/AEffect: non_coding_transcript_exon_variant/);
    assert.doesNotMatch(page.variant, /RNR|downstream|intergenic/);
    console.log('ok mitochondrial lookup: m.3243A>G is in MT-TL1');
}

// 4. A deletion Ensembl spells three bases further along.
{
    const page = await runLookup(cassettes.mt3271del);
    assertResolved(assert, page, 'm.3271del');
    assert.match(page.variant, /g\.: chrMT:g\.3271del/);
    assert.match(page.variant, /Gene: MT-TL1Transcript: ENST00000386347\.1/);
    assert.match(page.variant, /c\.: n\.44delp\.: N\/AEffect: non_coding_transcript_exon_variant/);
    console.log('ok mitochondrial lookup: m.3271del keeps its transcript');
}

// 5. Ensembl not answering. Nothing computed the coding change, so none is
//    shown — rather than snpEff's "intergenic" — and CADD still names the gene.
{
    const page = await runLookup(cassettes.mt8993, { without: [/ensembl\.org/] });
    assertResolved(assert, page, 'm.8993T>G, Ensembl down');
    assert.match(page.variant, /Gene: MT-ATP6c\.: N\/Ap\.: N\/AEffect: NON_SYNONYMOUS/);
    assert.doesNotMatch(page.variant + page.panel, /RNR2|intergenic|MT-ATP8/);
    assert.match(page.clinvar, /Variation ID: 9641/);
    console.log('ok mitochondrial lookup: without Ensembl, no invented consequence');
}

// 6. Not a mitochondrial matter, but the same list: snpEff's protein-structure
//    contact rows are not transcripts. TP53 R213* carries 762 of them (three are
//    kept in the recording) and each was listed under "Other transcripts".
{
    const record = recordOf(cassettes.tp53Stop, 'chr17:g.7578212G>A');
    const contacts = record.snpeff.ann.filter(a => a.feature_type === 'interaction');
    assert.equal(contacts.length, 3);
    assert.equal(contacts[0].effect, 'structural_interaction_variant');
    const page = await runLookup(cassettes.tp53Stop);
    assertResolved(assert, page, 'TP53 R213*');
    assert.match(page.variant, /Transcript: NM_000546\.5c\.: c\.637C>T/);
    assert.doesNotMatch(page.variant, /structural_interaction_variant|4QO1|:B_/);
    console.log('ok snpEff rows: protein-structure contacts are not listed as transcripts');
}

// 6b. The GRCh38 box ticked. The mitochondrial genome is one sequence in both
//     assemblies, so there is nothing to convert; but every spelling is read as
//     "chrMT:g.…" by now, and that was sent to Ensembl's assembly map all the
//     same — failing the lookup ("could not map chrMT:g.8993T>G (GRCh38) to
//     GRCh37") whenever the map was down, as it is in this recording.
for (const typed of ['M 8993 T G', 'chrMT:g.8993T>G', 'NC_012920.1:m.8993T>G', 'm.8993T>G']) {
    const page = await runLookup(cassettes.mt8993, { input: typed, grch38: true });
    assertResolved(assert, page, `${typed}, GRCh38 box ticked`);
    assert.deepEqual(asked(page, /\/map\//), [], `${typed}: the assembly map is not asked`);
    assert.match(page.variant, /g\.: chrMT:g\.8993T>G/, typed);
    assert.match(page.variant, /c\.: c\.467T>Gp\.: p\.Leu156Arg/, typed);
    assert.match(page.panel, /Map GRCh38 → GRCh37 · Ensemblnot needed — the mitochondrial sequence is the same in both assemblies/, typed);
}
console.log('ok mitochondrial lookup: the GRCh38 box changes nothing, and asks nothing');

// 7. gnomAD. Mitochondrial variants are in a call set of their own, counted by
//    heteroplasmy, which neither half of the gnomAD card read: v2 has none, and
//    the v4 proxy — nuclear variants only — answered every one of them with a
//    liftover failure that the card reported as
//
//        gnomAD v4.1 unavailable: the Ensembl liftover service is not responding,
//        so M:8993 could not be converted to GRCh38. This is usually temporary —
//        try again shortly.
//
//    beside a "View on gnomAD (v2.1.1)" link to a page saying v2 has no
//    mitochondrial variants. (The answers below are gnomAD's own, recorded.)
{
    const recorded = {
        'M-1555-A-G': { data: { mitochondrial_variant: { variant_id: 'M-1555-A-G', an: 56401, ac_hom: 63, ac_het: 11, max_heteroplasmy: 1 } } },
        'M-3243-A-G': { data: { mitochondrial_variant: { variant_id: 'M-3243-A-G', an: 56383, ac_hom: 0, ac_het: 6, max_heteroplasmy: 0.464 } } },
        'M-8993-T-G': { data: { mitochondrial_variant: { variant_id: 'M-8993-T-G', an: 56433, ac_hom: 0, ac_het: 0, max_heteroplasmy: 0 } } }
    };
    const notFound = { errors: [{ message: 'Variant not found' }], data: { mitochondrial_variant: null } };
    // gnomAD's API as a route for the replay harness; `asked` collects the ids.
    const gnomad = (answer = (id, respond) => respond(recorded[id] || notFound)) => {
        const asked = [];
        const route = (url, key, respond, options) => {
            if (url.hostname !== 'gnomad.broadinstitute.org') return undefined;
            const body = JSON.parse(options.body);
            if (!/mitochondrial_variant/.test(body.query)) return undefined;     // (gene constraint)
            asked.push(body.variables.id);
            return answer(body.variables.id, respond);
        };
        return { route, asked };
    };
    const gnomadLinks = page => [...page.document.querySelectorAll('.card')]
        .filter(card => card.querySelector('h3')?.textContent === 'gnomAD')
        .flatMap(card => [...card.querySelectorAll('a')])
        .filter(a => !/display:\s*none/.test(a.getAttribute('style') || ''))
        .map(a => a.getAttribute('href'));
    const never = /liftover|v4\.1|try again|No v2\.1\.1 data/;

    // m.1555A>G: 63 homoplasmic and 11 heteroplasmic carriers among 56,401.
    {
        const api = gnomad();
        const page = await runLookup(cassettes.mt1555, { route: api.route });
        assertResolved(assert, page, 'm.1555A>G');
        assert.deepEqual(api.asked, ['M-1555-A-G'], 'asked once, under the position as it is: no liftover');
        assert.deepEqual(page.api('gnomad-v4'), [], 'the nuclear-variant proxy is not asked');
        assert.match(page.cards.gnomAD, /gnomAD v2\.1\.1 has no mitochondrial variants\./);
        assert.match(page.cards.gnomAD, /Mitochondrial call set \(gnomAD API\)/);
        assert.match(page.cards.gnomAD, /Homoplasmic AF: 0\.1117% AC\/AN: 63\/56401/);
        assert.match(page.cards.gnomAD, /Heteroplasmic AF: 0\.0195% AC\/AN: 11\/56401/);
        assert.match(page.cards.gnomAD, /Highest level seen: 100%\./);
        assert.doesNotMatch(page.cards.gnomAD, never);
        assert.deepEqual(gnomadLinks(page), ['https://gnomad.broadinstitute.org/variant/M-1555-A-G?dataset=gnomad_r4'],
            'one link, to the page that exists');
    }
    // m.3243A>G: heteroplasmic only, at most 46.4 % of a carrier's copies.
    {
        const page = await runLookup(cassettes.mt3243, { route: gnomad().route });
        assert.match(page.cards.gnomAD, /Homoplasmic AF: 0\.0000% AC\/AN: 0\/56383/);
        assert.match(page.cards.gnomAD, /Heteroplasmic AF: 0\.0106% AC\/AN: 6\/56383/);
        assert.match(page.cards.gnomAD, /Highest level seen: 46\.4%\./);
    }
    // m.8993T>G: listed, with no carrier at 10 % or more.
    {
        const page = await runLookup(cassettes.mt8993, { route: gnomad().route });
        assert.match(page.cards.gnomAD, /Homoplasmic AF: 0\.0000% AC\/AN: 0\/56433Heteroplasmic AF: 0\.0000% AC\/AN: 0\/56433/);
        assert.doesNotMatch(page.cards.gnomAD, /Highest level seen/);
        assert.doesNotMatch(page.cards.gnomAD, never);
    }
    // m.3271del: gnomAD does not hold it. Asked under VCF alleles read off the reference.
    {
        const api = gnomad();
        const page = await runLookup(cassettes.mt3271del, { route: api.route });
        assert.deepEqual(api.asked, ['M-3270-CT-C']);
        assert.match(page.cards.gnomAD, /Not in gnomAD's mitochondrial call set \(queried: M-3270-CT-C\)\./);
        assert.deepEqual(gnomadLinks(page), ['https://gnomad.broadinstitute.org/variant/M-3270-CT-C?dataset=gnomad_r4']);
    }
    // gnomAD not answering is not "gnomAD does not have it".
    for (const down of [(id, respond) => respond({ error: 'unavailable' }, 503),
        (id, respond) => respond({ errors: [{ message: 'Internal server error' }], data: null })]) {
        const page = await runLookup(cassettes.mt1555, { route: gnomad(down).route });
        assertResolved(assert, page, 'm.1555A>G, gnomAD down');
        assert.match(page.cards.gnomAD, /gnomAD's mitochondrial call set could not be read\./);
        assert.doesNotMatch(page.cards.gnomAD, /Not in gnomAD|Homoplasmic AF/);
        assert.deepEqual(gnomadLinks(page), ['https://gnomad.broadinstitute.org/variant/M-1555-A-G?dataset=gnomad_r4']);
    }
    // A nuclear variant's card is as it was: v2 link, and the proxy asked.
    {
        const api = gnomad();
        const page = await runLookup(cassettes.brafV600E, { route: api.route });
        assert.deepEqual(api.asked, []);
        assert.equal(page.api('gnomad-v4').length, 1);
        assert.match(page.cards.gnomAD, /v2\.1\.1 · GRCh37.*v4\.1 · GRCh38 \(gnomAD API\)/);
        assert(gnomadLinks(page).some(href => /\/variant\/7-140453136-A-T\?dataset=gnomad_r2_1$/.test(href)));
    }
    console.log('ok mitochondrial lookup: gnomAD\'s mitochondrial call set, not a liftover error');

    // The pieces.
    for (const name of ['chrMT', 'MT', 'chrM', 'M', 'mt', 'chrm']) assert.equal(h.isMitochondrialChromosome(name), true, name);
    for (const name of ['1', 'X', 'chrX', 'MTOR', 'chr17', '', null]) assert.equal(h.isMitochondrialChromosome(name), false, String(name));
    assert.equal(h.gnomadMitochondrialId('8993', 't', 'g'), 'M-8993-T-G');
    assert.equal(h.gnomadMitochondrialId(3270, 'CT', 'C'), 'M-3270-CT-C');
    assert.equal(h.gnomadMitochondrialId(8993, 'T', '-'), '', 'no id without plain bases');
    assert.equal(h.gnomadMitochondrialId('x', 'T', 'G'), '');
    assert.equal(h.buildGnomadVariantUrl('', 'chrMT:g.8993T>G', '', null),
        'https://gnomad.broadinstitute.org/variant/M-8993-T-G?dataset=gnomad_r4');
    assert.equal(h.buildGnomadVariantUrl('', 'chrMT:g.3271del', '', null), '', 'no v2 region view for the mitochondrial genome');
    assert.equal(h.buildGnomadVariantUrl('', 'chrMT:g.3271del', '', null, { chrom: 'MT', pos: '3270', ref: 'CT', alt: 'C' }),
        'https://gnomad.broadinstitute.org/variant/M-3270-CT-C?dataset=gnomad_r4');
    assert.equal(h.buildGnomadVariantUrl('', 'chr7:g.140453136A>T', '', null),
        'https://gnomad.broadinstitute.org/variant/7-140453136-A-T?dataset=gnomad_r2_1');
    const realFetch = globalThis.fetch;
    try {
        const sent = [];
        const answerWith = (data, status = 200) => {
            globalThis.fetch = async (url, options) => {
                sent.push({ url: String(url), body: JSON.parse(options.body) });
                return new Response(JSON.stringify(data), { status });
            };
        };
        answerWith(recorded['M-1555-A-G']);
        const found = await h.fetchGnomadMitochondrial(1555, 'A', 'G');
        assert.equal(sent[0].url, 'https://gnomad.broadinstitute.org/api');
        assert.deepEqual(sent[0].body.variables, { id: 'M-1555-A-G' });
        assert.match(sent[0].body.query, /mitochondrial_variant\(variant_id: \$id, dataset: gnomad_r4\)/);
        assert.equal(found.status, 'found');
        assert.deepEqual([found.variant_id, found.an, found.ac_hom, found.ac_het, found.max_heteroplasmy], ['M-1555-A-G', 56401, 63, 11, 1]);
        assert.equal(found.af_hom.toFixed(6), '0.001117');
        assert.equal(found.af_het.toFixed(6), '0.000195');
        assert.match(found.note, /95%.*homoplasmic.*10-95%.*heteroplasmic/, 'the AI review is told what the counts mean');
        answerWith(notFound);
        assert.deepEqual(await h.fetchGnomadMitochondrial(8993, 'T', 'A'),
            { status: 'not_found', dataset: 'gnomAD mitochondrial call set', variant_id: 'M-8993-T-A' });
        answerWith({ data: { mitochondrial_variant: null } });
        assert.equal((await h.fetchGnomadMitochondrial(8993, 'T', 'A')).status, 'not_found');
        answerWith({ errors: [{ message: 'Mitochondrial variants are not available for gnomAD v2' }], data: { mitochondrial_variant: null } });
        await assert.rejects(h.fetchGnomadMitochondrial(8993, 'T', 'A'), /not available/, 'an error that is not "not found" is a failure');
        answerWith({}, 502);
        await assert.rejects(h.fetchGnomadMitochondrial(8993, 'T', 'A'), /502/);
        sent.length = 0;
        assert.equal(await h.fetchGnomadMitochondrial(8993, 'T', '-'), null);
        assert.equal(sent.length, 0, 'nothing is asked without alleles');
    } finally {
        globalThis.fetch = realFetch;
    }
    console.log('ok gnomAD mitochondrial call set: found, not found, and could not be asked');
}

// 8. A variant between genes. m.16300A>C is in the control region, 277 bp
//    upstream of MT-TP (VEP's rows below, as recorded). The gene shown was the
//    first row's — MT-ND4, 4 kb away, which the variant lies *downstream* of —
//    beside the effect "upstream_gene_variant".
{
    const row = (gene_symbol, term, distance) => ({ gene_symbol, consequence_terms: [term], distance });
    const controlRegion = [
        row('MT-ND4', 'downstream_gene_variant', 4163), row('MT-ND5', 'downstream_gene_variant', 2152),
        row('MT-ND6', 'upstream_gene_variant', 1627), row('MT-CYB', 'downstream_gene_variant', 413),
        row('MT-TH', 'downstream_gene_variant', 4094), row('MT-TS2', 'downstream_gene_variant', 4035),
        row('MT-TL2', 'downstream_gene_variant', 3964), row('MT-TE', 'upstream_gene_variant', 1558),
        row('MT-TT', 'downstream_gene_variant', 347), row('MT-TP', 'upstream_gene_variant', 277)
    ];
    assert.equal(await h.resolveGeneSymbolFromVep(controlRegion, null, 'upstream_gene_variant'), 'MT-TP',
        'the nearest gene the reported effect is about');
    assert.equal(await h.resolveGeneSymbolFromVep([...controlRegion].reverse(), null, 'upstream_gene_variant'), 'MT-TP',
        'whatever order VEP lists them in');
    // Without the effect to go by, the nearest gene of all.
    assert.equal(await h.resolveGeneSymbolFromVep(controlRegion, null), 'MT-TP');
    assert.equal(await h.resolveGeneSymbolFromVep(controlRegion.slice(0, 9), null), 'MT-TT');
    // A gene the variant is in still comes before any neighbour, however near…
    const within = { gene_symbol: 'PARP3', consequence_terms: ['inframe_deletion'], hgvsc: 'NM_005485.6:c.1_3del', hgvsp: 'NP_005476.4:p.Met1del' };
    assert.equal(await h.resolveGeneSymbolFromVep([row('RRP9', 'upstream_gene_variant', 12), within], null, 'inframe_deletion'), 'PARP3');
    // …and of two genes it is in, the one the reported effect belongs to.
    const intronic = { gene_symbol: 'GENEB', consequence_terms: ['intron_variant'], hgvsc: 'NM_000002.1:c.100+20A>G' };
    const missense = { gene_symbol: 'GENEA', consequence_terms: ['missense_variant'], hgvsc: 'NM_000001.1:c.5A>G', hgvsp: 'NP_000001.1:p.Lys2Arg' };
    assert.equal(await h.resolveGeneSymbolFromVep([intronic, missense], null, 'missense_variant'), 'GENEA');
    assert.equal(await h.resolveGeneSymbolFromVep([intronic, missense], null), 'GENEB', 'unchanged when the effect is not given');
    console.log('ok gene from VEP: the gene the effect is about, nearest first');
}

// 8b. ...and the page hands VEP's reported effect over, at both places where it
//     names the gene from VEP's rows (the function's own tests above cannot
//     tell). F508del by its g. name, with nothing from MyVariant: VEP gives a
//     row for CFTR (inframe_deletion) and one for CFTR-AS1, the antisense gene
//     the deletion is merely intronic in. Listed in the other order than they
//     were recorded, the two stand equal until the effect is taken into account
//     — and the variant was labelled CFTR-AS1.
{
    const cassette = cassettes.cftrF508del;
    const vepKey = Object.keys(cassette.responses).find(key => key.includes('/vep/human/hgvs/'));
    const recorded = cassette.responses[vepKey].body[0];
    assert.deepEqual(recorded.transcript_consequences.map(row => row.gene_symbol), ['CFTR', 'CFTR-AS1']);
    assert.equal(recorded.most_severe_consequence, 'inframe_deletion');
    const antisenseFirst = [{ ...recorded, transcript_consequences: [...recorded.transcript_consequences].reverse() }];
    const route = (url, key, respond) => (/\/vep\/human\/hgvs\//.test(key) ? respond(antisenseFirst) : undefined);
    // The recoder resolves the allele and MyVariant holds no record of it...
    const viaRecoder = await runLookup(cassette, { input: 'chr7:g.117199645_117199647del', route, without: [/myvariant\.info/] });
    assertResolved(assert, viaRecoder, 'F508del by g. name, no MyVariant record');
    assert.match(viaRecoder.variant, /Gene: CFTR(?!-AS1)/);
    assert.match(viaRecoder.panel, /Resolve nomenclature · Ensembl\s*\d+ transcripts?/, 'the recoder answered');
    assert(viaRecoder.api('civic').every(asked => /gene=CFTR(&|$)/.test(asked)), 'and the gene cards are CFTR\'s');
    // ...and with the recoder not answering either, where VEP is the last resort.
    const lastResort = await runLookup(cassette, { input: 'chr7:g.117199645_117199647del', route, without: [/myvariant\.info/, /variant_recoder/] });
    assertResolved(assert, lastResort, 'F508del by g. name, VEP alone');
    assert.match(lastResort.variant, /Gene: CFTR(?!-AS1)/);
    assert.doesNotMatch(lastResort.panel, /Resolve nomenclature · Ensembl\s*\d+ transcripts?/, 'the recoder did not answer');
    console.log('ok gene from VEP: the reported effect reaches the ranking, by either route');
}

console.log('Mitochondrial tests passed');
