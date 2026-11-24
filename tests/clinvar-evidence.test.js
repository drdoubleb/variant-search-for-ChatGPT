import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import handler, { parseClinvarVariant } from '../api/clinvar-variant.js';
import { captureResponse } from './support/api-response.js';

const fixture = id => readFileSync(new URL(`./fixtures/clinvar-${id}.xml`, import.meta.url), 'utf8');
const retrievedAt = '2026-11-24T16:13:00Z';
const m237l = parseClinvarVariant(fixture(2774649), '2774649', retrievedAt);
assert.equal(m237l.germline.description, 'Uncertain significance');
assert.equal(m237l.somatic, null);
assert.equal(m237l.record.variationId, '2774649');
assert.equal(m237l.record.accession, 'VCV002774649');
assert.equal(m237l.record.retrievedAt, retrievedAt);
assert.match(m237l.record.name, /c\.709A>T/);
assert.equal(m237l.submissions.length, 1);
const [color] = m237l.submissions;
assert.equal(color.accession, 'SCV004359989');
assert.equal(color.version, '2');
assert.equal(color.lastEvaluated, '2023-06-20');
assert.match(color.submitter, /Color Health/);
assert.equal(color.conditions[0].name, 'Hereditary cancer-predisposing syndrome');
assert.deepEqual(color.conditions[0].identifiers, [{ db: 'MedGen', id: 'C0027672' }]);
assert.deepEqual(color.assertionCriteria, ['ACMG Guidelines, 2015']);
assert.deepEqual(color.classifications.map(c => c.type), ['germline']);
for (const phrase of [/partially reduced transactivation activity/, /inconsistent results in human cell growth suppression/, /functional in human cell proliferation/, /threshold >= 0\.7/]) {
    assert.match(color.classificationComments[0], phrase);
}
assert.deepEqual(color.citations.filter(c => c.scope === 'variant_evidence').map(c => c.pmid).sort(), ['12826609', '29979965', '30224644']);
assert.deepEqual(color.citations.filter(c => c.scope === 'assertion_method').map(c => c.pmid), ['25741868']);

const m237i = parseClinvarVariant(fixture(142714), 142714);
assert.equal(m237i.germline.description, 'Pathogenic');
assert.equal(m237i.oncogenicity.description, 'Oncogenic');
assert.equal(m237i.somatic.description, 'Tier II - Potential');
assert.equal(m237i.submissions.length, 5);
assert.ok(m237i.record.reportedSubmissionCount > m237i.submissions.length);
assert.ok(m237i.submissions.some(s => s.reviewStatus === 'reviewed by expert panel'));
assert.ok(m237i.submissions.some(s => s.classifications[0].description === 'Uncertain significance' && !s.classificationComments.length));
assert.equal(m237i.submissions.find(s => s.classifications[0].type === 'somatic').classifications[0].assertionType, 'Diagnostic');
assert.ok(m237i.submissions.some(s => s.classifications[0].type === 'oncogenicity'));
const repeated = m237i.submissions.find(s => s.accession === 'SCV000545260');
assert.equal(repeated.citations.filter(c => c.pmid === '12826609').length, 1);
const braf = parseClinvarVariant(fixture(375941), 375941);
assert.deepEqual(braf.somaticConditions.map(c => c.condition).sort(), ['Melanoma', 'Spitz Melanocytoma']);
assert.equal(braf.submissions.find(s => s.submitter.includes('CIViC')).classifications[0].drugs, 'Dabrafenib;Trametinib');

// Edge cases constructed in the schema's structure exercise scoping and decoding.
const edge = `<VariationArchive VariationID="7" Accession="VCV000000007">
<ClassifiedRecord><Classifications/><ClinicalAssertionList><ClinicalAssertion ID="22">
  <ClinVarAccession Accession="SCV000000001" Version="1" SubmitterName="A &amp; B"/>
  <Classification><NoClassification>not provided</NoClassification>
    <ExplanationOfClassification>Limited evidence</ExplanationOfClassification>
    <Comment>&lt;strong&gt;literal&lt;/strong&gt; &amp; 5 &lt; &#55;</Comment>
  </Classification>
  <Comment>Record wording</Comment>
  <ObservedInList>
    <ObservedIn><Sample><Origin>somatic</Origin></Sample><ObservedData>
      <Attribute Type="Description">Recorded &quot;signal&quot;</Attribute><Comment>Observation wording</Comment>
      <Citation><ID Source="PubMed">765</ID></Citation>
    </ObservedData><TraitSet><Trait><Name><ElementValue Type="Preferred">Sample condition</ElementValue></Name></Trait></TraitSet></ObservedIn>
    <ObservedIn><Sample/><FunctionalData><FunctionalEffect Source="Laboratory">Reduced activity</FunctionalEffect>
      <FunctionalConsequence>Observed</FunctionalConsequence><Result>Partial</Result>
    </FunctionalData></ObservedIn>
  </ObservedInList>
  <TraitSet Type="Disease"><Trait><Name><ElementValue Type="Preferred">First &amp; related</ElementValue></Name></Trait>
    <Trait><Name><ElementValue Type="Preferred"><Emphasis>Second</Emphasis> &#8211; other</ElementValue></Name></Trait>
    <Comment>Condition note</Comment>
  </TraitSet>
</ClinicalAssertion></ClinicalAssertionList></ClassifiedRecord></VariationArchive>`;
const [parsed] = parseClinvarVariant(edge, 7).submissions;
assert.equal(parsed.submitter, 'A & B');
assert.equal(parsed.classifications[0].type, 'unclassified');
assert.deepEqual(parsed.classificationComments, ['Limited evidence', '<strong>literal</strong> & 5 < 7']);
assert.deepEqual(parsed.submitterComments, ['Record wording']);
assert.deepEqual(parsed.conditions.map(c => c.name), ['First & related', 'Second – other']);
assert.equal(parsed.observations[0].conditions[0].name, 'Sample condition');
assert.deepEqual(parsed.observations[0].comments, ['Observation wording']);
assert.equal(parsed.observations[0].descriptions[0].text, 'Recorded "signal"');
assert.deepEqual(parsed.observations[1].descriptions[0], { type: 'FunctionalEffect', text: 'Reduced activity', source: 'Laboratory' });
assert.deepEqual(parsed.citations, [{ pmid: '765', scope: 'variant_evidence', locations: ['observed_data'] }]);

const realFetch = globalThis.fetch;
try {
    let source;
    let fetches = 0;
    globalThis.fetch = async address => {
        fetches++;
        const query = new URL(address).searchParams;
        assert.equal(query.get('id'), '2774649');
        assert.ok(query.has('is_variationid'));
        return new Response(source);
    };
    async function call(id) {
        return handler({ method: 'GET', headers: {}, query: { id } }, captureResponse());
    }
    assert.equal((await call(['2774649'])).statusCode, 400);
    assert.equal(fetches, 0);
    source = fixture(2774649);
    assert.equal((await call('2774649')).body.submissions[0].accession, color.accession);
    for (source of [fixture(142714), '<ClinVarResult-Set><ERROR>Service error</ERROR></ClinVarResult-Set>',
        '<ClinVarResult-Set><VariationArchive></ClinVarResult-Set>', '<html><body>Gateway failure</body></html>',
        '<VariationArchive VariationID="2774649" Accession="VCV002774649"><ClassifiedRecord/></VariationArchive>',
        '<!DOCTYPE VariationArchive SYSTEM "https://invalid.test/clinvar.dtd">' + edge,
        '<!DOCTYPE VariationArchive [<!ENTITY external SYSTEM "https://example.invalid/ignored-entity">]>' + edge]) {
        const res = await call('2774649');
        assert.equal(res.statusCode, 502);
        assert.match(res.headers['Cache-Control'], /no-store/);
    }
} finally {
    globalThis.fetch = realFetch;
}
console.log('ClinVar evidence and upstream failure tests passed');
