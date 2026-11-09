// Regression fixtures reproduce both wrong ClinVar IDs reported for PTEN c.-8T>G.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const fixture = JSON.parse(readFileSync(new URL('./fixtures/pten-utr.json', import.meta.url)));
const genomic = 'chr10:g.89624219T>G';

for (const [query, key] of [
    ['PTEN c.-8T>G', 'mv-gene-space'],
    ['NM_000314.8 c.-8T>G', 'mv-nm-short-space'],
    ['NM_001304717.5 c.512T>G', 'mv-nm-long-space']
]) {
    assert.equal(h.selectMyVariantHit(fixture.searchHits[key], query), null, query);
    assert.equal(h.selectMyVariantHit(fixture.searchHits[key], query, genomic), null);
}
const exact = {
    _id: genomic,
    snpeff: { ann: { feature_id: 'NM_000314.8', hgvs_c: 'c.-8T>G', genename: 'PTEN' } },
    clinvar: { variant_id: 4757098, gene: { symbol: 'PTEN' }, hgvs: { coding: ['NM_000314.8:c.-8T>G'] } }
};
assert.equal(h.selectMyVariantHit([...fixture.searchHits['mv-gene-space'], exact], 'PTEN c.-8T>G'), exact);
assert.equal(h.selectMyVariantHit([exact], 'NM_000314.8:c.-8T>G'), exact);
assert.equal(h.selectMyVariantHit([exact], 'NM_000314:c.-8T>G'), exact);
assert.equal(h.selectMyVariantHit([exact], 'NM_000314.7:c.-8T>G'), null, 'respect explicit transcript versions');
assert.equal(h.selectMyVariantHit([exact], 'NM_001304717.5:c.-8T>G'), null, 'respect transcript identity');
assert.equal(h.selectMyVariantHit([exact], 'OTHER:c.-8T>G'), null);
assert.equal(h.selectMyVariantHit([exact], genomic), exact);
assert.equal(h.selectMyVariantHit([{ ...exact, _id: 'chr10:g.89624219T>C' }], 'PTEN:c.-8T>G', genomic), null);
assert.equal(h.selectMyVariantHit([exact], 'PTEN:c.8T>G'), null, 'preserve UTR minus sign');
assert.equal(h.selectMyVariantHit([exact], 'PTEN:c.*8T>G'), null, 'distinguish 3-prime UTR');
const proteinHit = { _id: 'chr7:g.140453136A>T', dbnsfp: { genename: ['BRAF','BRAF'], hgvsp: 'p.Val600Glu' } };
assert.equal(h.selectMyVariantHit([proteinHit], 'BRAF:p.V600E'), proteinHit);
assert.equal(h.selectMyVariantHit([proteinHit], 'BRAF:p.V600K'), null);

for (const wrong of ['chr10:g.89711981T>C', 'chr10:g.89711966T>C', 'chr10:g.89624219T>C']) {
    assert.throws(() => h.assertResolvedVariantIdentity(genomic, wrong, { _id: wrong }), /identity mismatch/);
    assert.throws(() => h.assertResolvedVariantIdentity('PTEN:c.-8T>G', genomic, { _id: wrong }), /identity mismatch/);
}
h.assertResolvedVariantIdentity(genomic, genomic, { _id: genomic });

const alleles = h.selectRecoderAlleles(fixture.recoder, genomic);
assert.equal(alleles.length, 1);
const recoderRows = h.transcriptsFromRecoderAllele(alleles[0]);
assert.equal(recoderRows.length, 6);
assert(recoderRows.every(t => t.protein === ''), 'independent arrays must not be zipped');
assert(h.transcriptsFromRecoderAllele({hgvsc: ['NM_1:c.1A>T'],hgvsp:['NP_2:p.Lys100Asn']}).every(t => !t.protein));
assert.equal(h.selectRecoderAlleles(fixture.recoder, 'NM_000314.8:c.-8T>C').length, 0);
assert.equal(h.selectRecoderAlleles(fixture.recoder, 'NM_001304717.5:c.-8T>G').length, 0);
const rows = h.transcriptsFromVep(fixture.vep[0].transcript_consequences);
for (const query of [genomic, 'PTEN:c.-8T>G', 'NM_000314.8:c.-8T>G']) {
    const display = h.buildTranscriptDisplay(rows, query);
    assert.equal(display.selected.transcript, 'NM_000314.8');
    assert.equal(display.selected.cDNA, 'c.-8T>G');
    assert.equal(display.canonicalProtein, '');
    assert.equal(h.transcriptEffect(display.selected), '5_prime_UTR_variant');
    assert.equal(display.transcriptsList.filter(t => t.canonical).length, 1);
    const extended = display.transcriptsList.find(t => t.transcript === 'NM_001304717.5');
    assert.equal(extended.cDNA, 'c.512T>G');
    assert.equal(extended.protein, 'p.Leu171Arg');
    assert.equal(extended.proteinAccession, 'NP_001291646.4');
}
assert.equal(h.buildTranscriptDisplay(rows,'NM_001304717.5:c.512T>G').selected.protein, 'p.Leu171Arg');
assert.equal(h.buildTranscriptDisplay(recoderRows,genomic).selected.transcript,'NM_000314.8');
assert.equal(h.buildTranscriptDisplay([{transcript:'NM_1', cDNA:'c.-8T>G',protein:'p.Leu171Arg'}]).canonicalProtein, '');
assert.equal(h.buildTranscriptDisplay([
    {transcript:'NM_1',cDNA:'c.1A>T'}, {transcript:'NM_2',cDNA:'c.-8T>G',maneSelect:true}
]).selected.transcript,'NM_2', 'MANE outranks accession heuristics, even for UTR variants');
assert.equal(h.buildCanonicalFromAnnotation({hgvsc:['NM_1:c.-8T>G'],hgvsp:['NP_2:p.Leu171Arg']}).canonicalProtein,'');
const region = [fixture.clinvar, {id:'631126',pos:89624219,title:'NM_000314.8(PTEN):c.-8T>C',germline:'Uncertain significance'}];
assert.equal(h.findExactClinvarRegionMatch(region, { pos:89624219, ref:'T', alt:'G' }).id, '4757098');
console.log('Variant identity and transcript regression checks passed');

const mane = h.maneSelectAccessions(fixture.mane);
assert(mane.includes('NM_000314.8'));
const marked = h.markManeTranscripts(rows, mane);
assert.equal(h.buildTranscriptDisplay(marked, genomic).selected.transcript, 'NM_000314.8');
assert.equal(h.buildTranscriptDisplay(marked, genomic).selected.maneSelect, true);
assert.equal(h.buildTranscriptDisplay(marked, 'NM_001304717.5:c.512T>G').selected.transcript, 'NM_001304717.5');
assert.equal(h.recoderTranscriptsForGenomic(fixture.recoder, 'chr10:g.89711981T>C').length, 0);
assert.equal(h.recoderTranscriptsForGenomic(fixture.recoder, genomic).length, 6);
assert.equal(h.buildTranscriptDisplay([{transcript:'NM_1',cDNA:'c.-10_20del',protein:'p.Met1?'}]).canonicalProtein,'p.Met1?');
