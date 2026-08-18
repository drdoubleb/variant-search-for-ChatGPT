/*
 * Tests for the retry layer in front of the annotation APIs.
 *
 * grch37.rest.ensembl.org — which every nomenclature lookup depends on —
 * intermittently returns 500/503 and has been measured taking 20-40 s to answer
 * a cold variant_recoder or VEP query. The app used to make exactly one attempt
 * behind a 6-7 s deadline and translate any failure into
 * "Variant not found. Please verify the genomic coordinate and reference
 * allele.", which blames the user for an outage.
 *
 * What matters here:
 *   - transport errors and 5xx are retried, 4xx are not (a 404 from
 *     MyVariant.info is a real answer: it holds no record for this indel)
 *   - the failure that escapes is tagged `retryable` so the caller can tell an
 *     outage apart from a variant that genuinely does not exist
 *
 * Helpers are copied verbatim from script.js (the project has no module system)
 * — KEEP IN SYNC when either side changes.
 *
 * Run with: node tests/fetch-retry.test.js
 */

// --- helpers copied verbatim from script.js -------------------------------

const RETRYABLE_HTTP_STATUS = new Set([408, 425, 429, 500, 502, 503, 504]);

function isRetryableFetchError(err) {
    const msg = String((err && err.message) || err || '');
    return /timed out after|Failed to fetch|NetworkError|network error|load failed|ECONNRESET|terminated/i.test(msg);
}

function describeUpstreamFailure(err) {
    if (!err) return 'failed';
    if (err.status) return `HTTP ${err.status}`;
    const msg = String(err.message || err);
    if (/timed out after/i.test(msg)) return 'timed out';
    if (/Failed to fetch|NetworkError|network error|load failed/i.test(msg)) return 'network error';
    const m = msg.match(/\((\d{3})\)/);
    if (m) return `HTTP ${m[1]}`;
    return msg.length > 60 ? `${msg.slice(0, 57)}…` : msg;
}

// fetchWithRetry, copied verbatim except that fetchWithTimeout is injected so
// the test can drive it without a network or a real AbortController.
function makeFetchWithRetry(fetchWithTimeout) {
    return async function fetchWithRetry(url, options = {}, timeoutMs = 6000, retryOpts = {}) {
        const attempts = Math.max(1, retryOpts.attempts ?? 3);
        const baseDelayMs = retryOpts.baseDelayMs ?? 400;
        const onAttempt = typeof retryOpts.onAttempt === 'function' ? retryOpts.onAttempt : null;
        let lastErr = null;
        for (let attempt = 1; attempt <= attempts; attempt++) {
            if (onAttempt) onAttempt({ attempt, attempts });
            try {
                const res = await fetchWithTimeout(url, options, timeoutMs);
                if (RETRYABLE_HTTP_STATUS.has(res.status) && attempt < attempts) {
                    lastErr = new Error(`Upstream returned ${res.status}`);
                    lastErr.retryable = true;
                    lastErr.status = res.status;
                } else {
                    return res;
                }
            } catch (err) {
                if (!isRetryableFetchError(err) || attempt === attempts) {
                    if (isRetryableFetchError(err)) err.retryable = true;
                    throw err;
                }
                lastErr = err;
                lastErr.retryable = true;
            }
            const delay = baseDelayMs * Math.pow(2, attempt - 1) * (0.7 + Math.random() * 0.6);
            await new Promise((resolve) => setTimeout(resolve, delay));
        }
        throw lastErr || new Error('Request failed');
    };
}

// --- test harness ---------------------------------------------------------

let failures = 0;
function check(name, actual, expected) {
    const pass = actual === expected;
    if (!pass) {
        failures++;
        console.error(`FAIL ${name}\n     expected: ${expected}\n     actual:   ${actual}`);
    } else {
        console.log(`ok   ${name}`);
    }
}

// Sequence-driven stub: each entry is a status number or an Error to throw.
function stubFetch(sequence) {
    const calls = { count: 0 };
    const fn = async () => {
        const next = sequence[Math.min(calls.count, sequence.length - 1)];
        calls.count++;
        if (next instanceof Error) throw next;
        return { status: next, ok: next >= 200 && next < 300 };
    };
    return { fn, calls };
}

const timeoutError = () => new Error('Request timed out after 6000ms');
const networkError = () => new Error('Failed to fetch');

(async () => {
    // ── 5xx is retried, and a later success is returned ────────────────────
    {
        const { fn, calls } = stubFetch([503, 503, 200]);
        const res = await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 3, baseDelayMs: 1 });
        check('503 then 200: returns the success', res.status, 200);
        check('503 then 200: made three attempts', calls.count, 3);
    }

    // The exact failure observed against grch37.rest.ensembl.org: a cold query
    // blows the deadline, the retry lands warm.
    {
        const { fn, calls } = stubFetch([timeoutError(), 200]);
        const res = await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 3, baseDelayMs: 1 });
        check('timeout then 200: returns the success', res.status, 200);
        check('timeout then 200: made two attempts', calls.count, 2);
    }

    {
        const { fn, calls } = stubFetch([networkError(), networkError(), 200]);
        const res = await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 3, baseDelayMs: 1 });
        check('network error recovers', res.status, 200);
        check('network error used all three attempts', calls.count, 3);
    }

    // ── 4xx is a real answer and must not be retried ───────────────────────
    {
        const { fn, calls } = stubFetch([404]);
        const res = await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 3, baseDelayMs: 1 });
        check('404 is returned as-is', res.status, 404);
        check('404 is not retried', calls.count, 1);
    }
    {
        const { fn, calls } = stubFetch([400]);
        await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 3, baseDelayMs: 1 });
        check('400 is not retried', calls.count, 1);
    }

    // ── Exhausted retries surface a retryable error ────────────────────────
    {
        const { fn, calls } = stubFetch([503, 503, 503]);
        const res = await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 3, baseDelayMs: 1 });
        // The last attempt returns the response rather than throwing, so the
        // caller can read the status and classify it.
        check('persistent 503 returns the last response', res.status, 503);
        check('persistent 503 used all attempts', calls.count, 3);
    }
    {
        const { fn, calls } = stubFetch([timeoutError()]);
        let caught = null;
        try {
            await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 2, baseDelayMs: 1 });
        } catch (e) { caught = e; }
        check('persistent timeout throws', Boolean(caught), true);
        check('persistent timeout is tagged retryable', caught && caught.retryable, true);
        check('persistent timeout used all attempts', calls.count, 2);
    }

    // ── Non-retryable throws escape immediately ────────────────────────────
    {
        const { fn, calls } = stubFetch([new TypeError('bad url')]);
        let caught = null;
        try {
            await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 3, baseDelayMs: 1 });
        } catch (e) { caught = e; }
        check('programming error is not retried', calls.count, 1);
        check('programming error is not tagged retryable', caught && caught.retryable, undefined);
    }

    // ── attempts:1 preserves single-shot behaviour ─────────────────────────
    {
        const { fn, calls } = stubFetch([503, 200]);
        const res = await makeFetchWithRetry(fn)('u', {}, 100, { attempts: 1, baseDelayMs: 1 });
        check('attempts:1 does not retry', calls.count, 1);
        check('attempts:1 returns the failure', res.status, 503);
    }

    // ── onAttempt reports progress for the panel ───────────────────────────
    {
        const seen = [];
        const { fn } = stubFetch([503, 503, 200]);
        await makeFetchWithRetry(fn)('u', {}, 100, {
            attempts: 3, baseDelayMs: 1,
            onAttempt: ({ attempt, attempts }) => seen.push(`${attempt}/${attempts}`)
        });
        check('onAttempt fires once per attempt', seen.join(' '), '1/3 2/3 3/3');
    }

    // ── describeUpstreamFailure produces short, honest labels ──────────────
    check('describes an HTTP status', describeUpstreamFailure(Object.assign(new Error('x'), { status: 503 })), 'HTTP 503');
    check('describes a timeout', describeUpstreamFailure(timeoutError()), 'timed out');
    check('describes a network error', describeUpstreamFailure(networkError()), 'network error');
    check('extracts a status from a message',
        describeUpstreamFailure(new Error('Variant recoder request failed (500): boom')), 'HTTP 500');
    check('falls back to the message', describeUpstreamFailure(new Error('No VEP HGVS data found')), 'No VEP HGVS data found');

    if (failures) {
        console.error(`\n${failures} test(s) failed`);
        process.exit(1);
    }
    console.log('\nAll fetch-retry tests passed');
})();
