import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import '../script.js';

const {
    parseGenomicHgvs, geneMapTranscript, geneMapVariantLocus, geneMapProteinSpan,
    geneMapDomains, geneMapGeometry
} = globalThis.__variantSearchHelpers;
const fixtures = JSON.parse(readFileSync(new URL('./fixtures/ensembl-gene-maps.json', import.meta.url), 'utf8'));

const tp53 = geneMapTranscript(fixtures.TP53.transcript, 'GRCh37');
assert.equal(tp53.exons.length, 11);
assert.equal(tp53.strand, -1);
assert.equal(tp53.exons[0].end, tp53.end);
assert.equal(geneMapTranscript(fixtures.TP53.transcript, 'GRCh38'), null);

for (const [hgvs, start, end] of [
    ['chr17:g.7577572T>A', 237, 237],    // M237L
    ['chr17:g.7577120C>T', 273, 273],    // R273H
    ['chr17:g.7577571_7577574del', 236, 237],
    ['chr17:g.7577572_7577573insA', 236, 237]
]) {
    assert.deepEqual(geneMapProteinSpan(tp53, parseGenomicHgvs(hgvs)), { start, end }, hgvs);
}
const intron = parseGenomicHgvs('chr17:g.7578000');
assert.deepEqual(geneMapVariantLocus(tp53, intron), { start: 7578000, end: 7578000 });
assert.equal(geneMapProteinSpan(tp53, intron), null);
assert.equal(geneMapProteinSpan(tp53, parseGenomicHgvs('chr17:g.7572927')), null, 'terminal stop codon');
assert.equal(geneMapVariantLocus(tp53, parseGenomicHgvs('chr7:g.7577572T>A')), null);
assert.equal(geneMapVariantLocus(tp53, parseGenomicHgvs('chr17:g.7571700T>A')), null, 'coordinate outside the transcript');
assert.equal(geneMapProteinSpan({ ...tp53, coding: tp53.coding.slice(1) }, parseGenomicHgvs('chr17:g.7577572T>A')), null,
    'incomplete CDS');

const egfr = geneMapTranscript(fixtures.EGFR.transcript, 'GRCh37');
assert.equal(egfr.strand, 1);
assert.deepEqual(geneMapProteinSpan(egfr, parseGenomicHgvs('chr7:g.55259515T>G')), { start: 858, end: 858 });
const malat1 = geneMapTranscript(fixtures.MALAT1.transcript, 'GRCh38');
const noncoding = parseGenomicHgvs('chr11:g.65500000');
assert.deepEqual(geneMapVariantLocus(malat1, noncoding), { start: 65500000, end: 65500000 });
assert.equal(malat1.protein, null);
assert.equal(geneMapProteinSpan(malat1, noncoding), null);

const domains = geneMapDomains(fixtures.TP53.domains, tp53.protein.length);
assert.deepEqual(domains.map(({ id, start, end }) => [id, start, end]), [
    ['PF08563', 6, 29], ['PF00870', 95, 288], ['PF07710', 318, 358]
]);
assert.equal(geneMapDomains([...fixtures.TP53.domains, { type: 'Pfam', id: 'PF00000', start: 1, end: 394 }], 393).length, 3);
assert.deepEqual(geneMapGeometry({ start: 100, end: 199 }, { start: 189, end: 198 }, true), { left: 1, width: 10 });
console.log('Gene map assembly, strand, coding span and domain tests passed');
