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
        model: 'openai/gpt-5-mini',
        openrouter_api_key: 'sk-or-v1-0123456789abcdef0123',
        context: { gene: 'PTEN', submitted_variant: 'p.His123Asp' }
    }
});

async function callWith(payload) {
    const realFetch = global.fetch;
    global.fetch = async () => ({ ok: true, status: 200, text: async () => JSON.stringify(payload) });
    try {
        const res = mkRes();
        await handler(mkReq(), res);
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
const ok = await callWith({ choices: [{ finish_reason: 'stop', message: { content: '{"pathogenicity":"Pathogenic","amp_tier":"Tier IA"}' } }] });
check('valid JSON returns 200', ok.code === 200);
check('valid JSON is parsed', ok.body.review.pathogenicity === 'Pathogenic');

const fenced = await callWith({ choices: [{ finish_reason: 'stop', message: { content: 'Here you go:\n```json\n{"pathogenicity":"Benign"}\n```' } }] });
check('JSON wrapped in prose/fences is recovered', fenced.code === 200 && fenced.body.review.pathogenicity === 'Benign');

console.log(`\nAI-review response tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
