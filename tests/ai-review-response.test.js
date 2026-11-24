/*
 * Tests for how api/ai-review.js handles OpenRouter responses that are not clean
 * JSON. These all used to surface to the user as a raw parser message —
 * "Expected ',' or ']' after array element in JSON at position 4348" — which reads
 * like a prompt bug but is actually the model running out of output tokens.
 *
 * The budget matters because max_tokens is shared with REASONING tokens on
 * reasoning models: a model that thinks for 1500 tokens had only 1500 left for the
 * answer under the old 3000 cap, so the JSON stopped mid-array (or came back empty
 * when reasoning consumed all of it).
 *
 * Run with: node tests/ai-review-response.test.js
 */

import handler from '../api/ai-review.js';

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; } else { failed++; console.error(`FAIL: ${name}`); }
}

const mkRes = () => ({
    code: 0, body: null, headers: {},
    setHeader(k, v) { this.headers[k] = v; },
    status(c) { this.code = c; return this; },
    json(b) { this.body = b; return this; },
    end() { return this; }
});

// BYO-key so the call skips Turnstile and the rate limiter.
const mkReq = () => ({
    method: 'POST',
    headers: { origin: 'https://drdoubleb.com' },
    body: {
        model: 'openai/gpt-6-luna',
        openrouter_api_key: 'sk-or-v1-0123456789abcdef0123',
        context: { gene: 'PTEN', submitted_variant: 'p.His123Asp' }
    }
});

async function callWith(payload, tumorType) {
    const realFetch = global.fetch;
    global.fetch = async () => ({ ok: true, status: 200, text: async () => JSON.stringify(payload) });
    try {
        const res = mkRes();
        const req = mkReq();
        if (tumorType) req.body.context.tumor_type = tumorType;
        await handler(req, res);
        return res;
    } finally {
        global.fetch = realFetch;
    }
}

// Cut off right after an object inside an array — the exact shape that produced
// "Expected ',' or ']' after array element" in the wild.
const CUTOFF = '{"pathogenicity":"Likely Pathogenic","amp_tier":"Tier IA","fda_approved_therapies":[{"drug":"capivasertib"}';

// A truncated generation is reported as a token-budget problem, not a parse error.
const truncated = await callWith({ choices: [{ finish_reason: 'length', message: { content: CUTOFF } }] });
check('truncated response returns 502', truncated.code === 502);
check('truncated response blames the token budget', /ran out of output tokens/.test(truncated.body.error));
check('truncated response does not leak a parser message', !/Expected ',' or '\]'/.test(truncated.body.error));
check('truncated response reports the cap in detail', /max_tokens=/.test(truncated.body.detail || ''));

// OpenRouter reports some models' cutoff only in native_finish_reason.
const nativeTrunc = await callWith({ choices: [{ native_finish_reason: 'length', message: { content: CUTOFF } }] });
check('native_finish_reason=length is also caught', /ran out of output tokens/.test(nativeTrunc.body.error));

// Reasoning consumed the whole budget, leaving no content.
const empty = await callWith({ choices: [{ finish_reason: 'stop', message: { content: '' } }] });
check('empty content returns 502', empty.code === 502);
check('empty content mentions internal reasoning', /internal reasoning/.test(empty.body.error));

// Model ignored the format and wrote prose.
const prose = await callWith({ choices: [{ finish_reason: 'stop', message: { content: 'PTEN H123D is likely pathogenic.' } }] });
check('prose response returns 502', prose.code === 502);
check('prose response is actionable', /text instead of JSON/.test(prose.body.error));

// Malformed JSON without a truncation flag still gets a readable message.
const malformed = await callWith({ choices: [{ finish_reason: 'stop', message: { content: CUTOFF } }] });
check('malformed response returns 502', malformed.code === 502);
check('malformed response is readable', /incomplete or malformed/.test(malformed.body.error));
check('malformed response does not leak a parser message', !/position \d+/.test(malformed.body.error));

// The happy path still works, including JSON wrapped in stray prose.
const mockGermline = { germline_pathogenicity: { classification: 'Pathogenic', condition: 'Example syndrome', rationale: 'Germline fixture' } };
const mockReview = {
    ...mockGermline,
    somatic_oncogenicity: { classification: 'Likely oncogenic', rationale: 'Somatic fixture' },
    clinical_significance: { status: 'assessed', tumor_type: 'Example tumor', tier: 'Tier III', evidence_level: null, rationale: 'Clinical evidence fixture' }
};
const wrapped = review => ({ choices: [{ finish_reason: 'stop', message: { content: JSON.stringify(review) } }] });
const ok = await callWith(wrapped(mockReview), 'Example tumor');
check('valid JSON returns 200', ok.code === 200);
check('the germline and somatic assessments stay independent', ok.body.review.pathogenicity === 'Pathogenic'
    && ok.body.review.germline_pathogenicity.condition === 'Example syndrome'
    && ok.body.review.somatic_oncogenicity.classification === 'Likely oncogenic');
check('an oncogenic finding can have Tier III clinical significance', ok.body.review.clinical_significance.tier === 'Tier III'
    && ok.body.review.amp_tier === 'Tier III');

const fenced = await callWith({ choices: [{ finish_reason: 'stop', message: { content: 'Here you go:\n```json\n{"germline_pathogenicity":{"classification":"Benign"}}\n```' } }] });
check('JSON wrapped in prose/fences is recovered', fenced.code === 200 && fenced.body.review.pathogenicity === 'Benign');

const noTumor = await callWith(wrapped({ ...mockReview, amp_tier: 'Tier IA' }));
check('a missing tumor type suppresses even a supplied model tier', noTumor.body.review.clinical_significance.status === 'insufficient_disease_context'
    && noTumor.body.review.clinical_significance.tier === null && noTumor.body.review.amp_tier === 'Insufficient disease context');
const otherTumor = await callWith(wrapped(mockReview), 'Different tumor');
check('a model assessment for a different tumor stays unassigned', otherTumor.body.review.clinical_significance.tier === null);
const one = { ...mockReview, clinical_significance: { ...mockReview.clinical_significance, tier: 'Tier I', evidence_level: 'A' } };
const supported = await callWith(wrapped(one), 'Example tumor');
check('a supported tier keeps its evidence level and compatibility field', supported.body.review.clinical_significance.tier === 'Tier I'
    && supported.body.review.amp_tier === 'Tier IA');
one.clinical_significance.evidence_level = 'D';
const incompatible = await callWith(wrapped(one), 'Example tumor');
check('incompatible tier and evidence level stay unassigned', incompatible.body.review.clinical_significance.status === 'insufficient_evidence'
    && incompatible.body.review.clinical_significance.tier === null);

console.log(`\nAI-review response tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
