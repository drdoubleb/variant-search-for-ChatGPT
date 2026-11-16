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

console.log('Mitochondrial tests passed');
