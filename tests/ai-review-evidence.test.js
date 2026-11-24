import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import handler from '../api/ai-review.js';
import { parseClinvarVariant } from '../api/clinvar-variant.js';
import { captureResponse } from './support/api-response.js';

const xml = readFileSync(new URL('./fixtures/clinvar-2774649.xml', import.meta.url), 'utf8');
const evidence = parseClinvarVariant(xml, 2774649);
const context = {
    submitted_variant: 'NM_000546.6:c.709A>T',
    user_notes: 'A note containing {{CONTEXT_JSON}} and {{SCHEMA_JSON}}',
    supplemental_card_data: {
        clinvar_variant_record: evidence,
        openfda_drug_labels: { results: [{ label: 'a'.repeat(740000) }, { label: 'b'.repeat(800000) }] }
    }
};
function promptContext(prompt) {
    const [, json] = prompt.split('Context JSON:\n\n');
    return JSON.parse(json.split('\n\nEnd of supplied context.')[0]);
}
const originalFetch = globalThis.fetch;
try {
    const prompts = [];
    globalThis.fetch = async (url, options) => {
        assert.equal(new URL(url).hostname, 'openrouter.ai');
        prompts.push(JSON.parse(options.body).messages[1].content);
        return Response.json({ choices: [{ message: { content: '{"pathogenicity":"VUS"}' } }] });
    };
    const call = body => handler({ method: 'POST', headers: {}, body }, captureResponse());
    const copy = await call({ context, mode: 'prompt' });
    assert.equal(prompts.length, 0);
    const review = await call({ context, openrouter_api_key: 'sk-local-prompt-test-0000' });
    assert.equal(copy.statusCode, 200);
    assert.equal(review.statusCode, 200);
    for (const [response, prompt] of [[copy, copy.body.prompt], [review, prompts[0]]]) {
        const payload = promptContext(prompt);
        assert.equal(payload.supplemental_card_data.openfda_drug_labels.results.length, 1);
        assert.deepEqual(payload.supplemental_card_data.clinvar_variant_record, evidence);
        assert.deepEqual(payload.context_omissions, [{ source: 'openfda_drug_labels', records_omitted: 1 }]);
        assert.deepEqual(response.body.context_omissions, payload.context_omissions);
        assert.ok(prompt.includes(context.user_notes));
    }
    const tooLarge = await call({ mode: 'prompt', context: { supplemental_card_data: { clinvar_variant_record: 'a'.repeat(1600000) } } });
    assert.equal(tooLarge.statusCode, 413);
    assert.match(tooLarge.body.error, /exceeds the prompt size limit/);
} finally {
    globalThis.fetch = originalFetch;
}
console.log('Both AI prompt paths preserve ClinVar evidence and report size limits');
