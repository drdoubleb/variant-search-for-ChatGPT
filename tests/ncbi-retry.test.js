/*
 * Tests for the shared NCBI E-utilities fetch helper (api/_ncbi.js).
 *
 * The ClinVar proxy endpoints share this helper. Its job is to survive the
 * transient 429 "Too Many Requests" responses that eUtils returns under the
 * 3 req/s anonymous rate limit — a single UI lookup fans out to several ClinVar
 * calls at once, so bursts routinely trip the limit. Before retry was added,
 * that 429 surfaced as a 502 and the "Same protein change" section silently
 * vanished.
 *
 * global.fetch is stubbed so no network is touched. Run with:
 *   node tests/ncbi-retry.test.js
 */

import { ncbiFetchResponse, ncbiFetchJson, ncbiApiKeyParam } from '../api/_ncbi.js';

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; }
    else { failed++; console.error(`  ✗ ${name}`); }
}

function makeResponse({ status = 200, body = '', headers = {} } = {}) {
    const h = new Map(Object.entries(headers).map(([k, v]) => [k.toLowerCase(), String(v)]));
    return {
        ok: status >= 200 && status < 300,
        status,
        headers: { get: (k) => (h.has(k.toLowerCase()) ? h.get(k.toLowerCase()) : null) },
        json: async () => JSON.parse(body),
        text: async () => body,
    };
}

// Speed the tests up: neutralise the backoff sleeps (they use setTimeout).
const realSetTimeout = global.setTimeout;
global.setTimeout = (fn) => realSetTimeout(fn, 0);

async function run() {
    // 1. Retries through two 429s, then succeeds on the third attempt.
    {
        let calls = 0;
        global.fetch = async () => {
            calls++;
            if (calls < 3) return makeResponse({ status: 429 });
            return makeResponse({ status: 200, body: JSON.stringify({ ok: true }) });
        };
        const data = await ncbiFetchJson('http://x');
        check('retries past 429s and returns parsed JSON on eventual success', data.ok === true);
        check('made exactly 3 attempts (2 retries)', calls === 3);
    }

    // 2. Persistent 429 exhausts the budget and throws (4 attempts total).
    {
        let calls = 0;
        global.fetch = async () => { calls++; return makeResponse({ status: 429 }); };
        let threw = false;
        try { await ncbiFetchResponse('http://x'); } catch { threw = true; }
        check('throws after exhausting retries on persistent 429', threw);
        check('caps at 4 attempts', calls === 4);
    }

    // 3. A non-retryable client error (400) fails fast without retrying.
    {
        let calls = 0;
        global.fetch = async () => { calls++; return makeResponse({ status: 400 }); };
        let threw = false;
        try { await ncbiFetchResponse('http://x'); } catch { threw = true; }
        check('throws on non-retryable 400', threw);
        check('does NOT retry a 400 (single attempt)', calls === 1);
    }

    // 4. Network/timeout errors are retried, then can succeed.
    {
        let calls = 0;
        global.fetch = async () => {
            calls++;
            if (calls === 1) throw new Error('network down');
            return makeResponse({ status: 200, body: JSON.stringify({ recovered: true }) });
        };
        const data = await ncbiFetchJson('http://x');
        check('retries after a network error and recovers', data.recovered === true && calls === 2);
    }

    // 5. Retryable 5xx (503) is retried like a 429.
    {
        let calls = 0;
        global.fetch = async () => {
            calls++;
            return calls === 1
                ? makeResponse({ status: 503 })
                : makeResponse({ status: 200, body: JSON.stringify({ ok: 1 }) });
        };
        const data = await ncbiFetchJson('http://x');
        check('retries a 503 and recovers', data.ok === 1 && calls === 2);
    }

    // 6. Immediate success makes exactly one call.
    {
        let calls = 0;
        global.fetch = async () => { calls++; return makeResponse({ status: 200, body: '{}' }); };
        await ncbiFetchJson('http://x');
        check('single attempt on immediate success', calls === 1);
    }

    // 7. API-key param helper reflects the env var.
    {
        const saved = process.env.NCBI_API_KEY;
        delete process.env.NCBI_API_KEY;
        check('no api_key param when env unset', ncbiApiKeyParam() === '');
        process.env.NCBI_API_KEY = 'abc 123';
        check('api_key param is url-encoded when env set', ncbiApiKeyParam() === '&api_key=abc%20123');
        if (saved === undefined) delete process.env.NCBI_API_KEY; else process.env.NCBI_API_KEY = saved;
    }

    console.log(`\nNCBI retry tests: ${passed} passed, ${failed} failed`);
    if (failed > 0) process.exit(1);
}

run();
