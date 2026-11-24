import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import handler from '../api/tp53.js';
import { captureResponse } from './support/api-response.js';

const csv = readFileSync(new URL('./fixtures/tp53-mutationview-r21.csv', import.meta.url), 'utf8');
const datasetUrl = 'https://tp53.cancer.gov/static/data/MutationView_r21.csv';
const previousUrl = process.env.TP53_MUTATION_DATASET_URL;
const originalFetch = globalThis.fetch;
process.env.TP53_MUTATION_DATASET_URL = datasetUrl;
let download = csv;
let requests = 0;
globalThis.fetch = async url => {
    assert.equal(url, datasetUrl);
    requests++;
    return new Response(download, { headers: { 'Content-Type': 'text/csv' } });
};
try {
    async function query(body, endpoint = handler) {
        const response = await endpoint({ method: 'POST', headers: {}, body: { gene: 'TP53', ...body } }, captureResponse());
        assert.equal(response.statusCode, 200);
        return response.body;
    }
    const protein = 'p.Met237Leu';
    const cdna = 'NM_000546.6:c.709A>T';
    const g37 = { assembly: 'GRCh37', chrom: '17', pos: 7577572, ref: 'T', alt: 'A' };
    const exact = await query({ protein, cdna, genomic_allele: g37 });
    assert.deepEqual(exact.selection, { status: 'exact', match_type: 'exact_genomic', candidate_count: 1 });
    assert.deepEqual(exact.matches[0].matched_genomic_allele, g37);
    assert.equal(exact.matches[0].mut_id, '3082');
    assert.equal(exact.matches[0].transcript_hgvs, 'NM_000546.5:c.709A>T');
    assert.equal(exact.pathogenicity.somatic_count, '3');
    assert.equal(exact.pathogenicity.splice_ai.DS_AG, '0.08');
    assert.equal(exact.pathogenicity.splice_ai.DS_AL, '0');
    assert.equal(exact.pathogenicity.ta_targets.NOXA, '78.2');
    assert.match(exact.functional_sources.structure_prediction.method, /Delaunay.*trained on the Kato/);
    const [kato, giacomelli, kotler] = exact.assay_comparison.studies;
    assert.equal(exact.assay_comparison.includes_mutation_view, true);
    assert.equal(kato.measurements[0].value, 40.6);
    assert.deepEqual(kato.measurements.filter(m => m.label === 'WAF1').map(m => m.value), [6.6]);
    assert.equal(giacomelli.scope, 'protein');
    assert.equal(giacomelli.classifications[0].value, 'Unclassified');
    assert.equal(giacomelli.measurements.find(m => m.label === 'p53 WT · nutlin-3').value, 0.7432024599950522);
    assert.equal(kotler.scope, 'coding_construct');
    assert.match(kotler.variant, /c\.709A>T.*302465/);
    assert.equal(kotler.measurements[0].value, -1.666323792);
    const otherLeu = exact.matches.find(m => m.mut_id === '3080');
    assert.equal(otherLeu.match_type, 'same_protein');
    assert.equal(otherLeu.somatic_count, '2');
    assert.equal(otherLeu.splice_ai.DS_AL, '0.01');
    assert.ok(exact.matches.some(m => m.match_type === 'same_codon'));
    const exact38 = await query({ protein, genomic: 'NC_000017.11:g.7674254T>A' });
    assert.equal(exact38.matches[0].matched_genomic_allele.assembly, 'GRCh38');
    const exactCdna = await query({ protein, cdna });
    assert.equal(exactCdna.selection.match_type, 'exact_cdna');
    assert.equal(exactCdna.matches[0].match_basis, 'identical_refseq_cds_for_coding_snv');
    const exactC = await query({ protein, genomic_allele: { ...g37, alt: 'G' } });
    assert.equal(exactC.pathogenicity.somatic_count, '2');
    assert.equal(exactC.matches[0].mut_id, '3080');
    assert.match(exactC.assay_comparison.studies.at(-1).variant, /c\.709A>C.*302463/);
    assert.equal(exactC.assay_comparison.studies.at(-1).measurements[0].value, -1.779194794);
    assert.ok((await query({ protein })).assay_comparison.studies.every(study => study.scope === 'protein'));
    const exactG = await query({ protein, genomic_allele: { ...g37, alt: 'C' } });
    assert.equal(exactG.assay_comparison.protein, 'p.M237V');
    assert.equal(exactG.assay_comparison.studies.length, 2);
    assert.equal(exactG.assay_comparison.studies[0].measurements.find(m => m.label === 'WAF1').value, 22.8);
    assert.deepEqual(exactG.assay_comparison.studies[1].classifications.map(m => m.value), ['DNE_LOF']);
    const conflicting = await query({ protein, cdna, genomic_allele: { ...g37, ref: 'C' } });
    assert.equal(conflicting.selection.status, 'related_only');
    assert.equal(conflicting.pathogenicity, null);
    const coordinate = await query({ genomic: 'chr17:g.7577572', assembly: 'GRCh37' });
    assert.ok(coordinate.matches.length > 1);
    assert.ok(coordinate.matches.every(m => m.match_type === 'coordinate_only'));
    assert.equal(coordinate.selection.status, 'related_only');
    assert.equal(coordinate.pathogenicity, null);
    assert.equal(coordinate.assay_comparison, null);

    const r175hGenomic = await query({ genomic: 'NC_000017.10:g.7578406C>T' });
    const r175hProtein = await query({ protein: 'p.Arg175His' });
    assert.equal(r175hGenomic.selection.status, 'exact');
    assert.equal(r175hProtein.selection.status, 'related_only');
    const studies = r175hProtein.assay_comparison.studies;
    assert.deepEqual(studies.map(s => s.study), ['Kato et al. (2003)', 'Kakudo et al. (2005)', 'Giacomelli et al. (2018)', 'Dominant-negative transactivation']);
    const results = ({ measurements, classifications }) => ({ measurements, classifications });
    assert.deepEqual(studies.map(results), r175hGenomic.assay_comparison.studies.map(results));
    const saos2 = studies[1].measurements;
    assert.equal(saos2.find(m => m.label === 'WAF1').value, -2.2);
    assert.equal(saos2.find(m => m.label === 'BAX').value, 0);
    assert.equal(saos2.find(m => m.label === 'Sub-G1 apoptosis').value, 30);
    assert.deepEqual(studies[2].classifications.map(m => m.value), ['DNE_LOF']);
    assert.ok(studies.every(study => study.measurements.every(result => Number.isFinite(result.value))));
    for (const input of [
        { genomic: 'chr17:g.7577572T>A' },
        { genomic: 'chr17:g.7577572T>A', assembly: 'GRCh38' },
        { genomic_allele: { ...g37, chrom: '16' } },
        { protein, cdna: 'c.709A>T' },
        { protein, cdna: 'NM_000546:c.709A>T' },
        { protein, cdna: 'NM_001126112.3:c.709A>T' },
        { protein, cdna: 'NM_000546.6:c.709A>G' , genomic_allele: { ...g37, pos: 7674254 } },
        { protein }
    ]) {
        assert.equal((await query(input)).pathogenicity, null, JSON.stringify(input));
    }
    assert.equal(requests, 1, 'the downloaded CSV is reused');

    // A duplicate exact dataset row (e.g. sample-level data) stays ambiguous.
    const alternate = (await import('../api/tp53.js?duplicate-row')).default;
    download = csv + csv.split('\n').find(line => line.startsWith('3082,')) + '\n';
    const ambiguous = await query({ protein, cdna, genomic_allele: g37 }, alternate);
    assert.equal(ambiguous.selection.status, 'ambiguous');
    assert.equal(ambiguous.match_count, 2);
    assert.equal(ambiguous.pathogenicity, null);
    assert.equal(ambiguous.prevalence, '');
    assert.ok(ambiguous.assay_comparison.studies.every(study => study.scope === 'protein'));
    assert.equal(ambiguous.matches.filter(m => m.match_type === 'exact_genomic').length, 2);

    // A hypothetical discrepancy between dataset rows keeps both measurements.
    const conflictingRows = (await import('../api/tp53.js?conflicting-values')).default;
    download = csv.split('\n').map(line => line.startsWith('3082,') ? line.replace(',6.6,23.4,', ',7.1,23.4,') : line).join('\n');
    const disagreement = await query({ protein }, conflictingRows);
    const measurements = disagreement.assay_comparison.studies[0].measurements;
    assert.deepEqual(measurements.filter(m => m.label === 'WAF1').map(m => [m.value, m.source_records]), [
        [6.6, ['MutationView 3080 · NM_000546.5:c.709A>C']], [7.1, ['MutationView 3082 · NM_000546.5:c.709A>T']]
    ]);
    assert.ok(!measurements.some(m => m.label.startsWith('Median')));
} finally {
    globalThis.fetch = originalFetch;
    if (previousUrl === undefined) delete process.env.TP53_MUTATION_DATASET_URL;
    else process.env.TP53_MUTATION_DATASET_URL = previousUrl;
}
console.log('TP53 exact and related identity tests passed');
