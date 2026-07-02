// Shared NCBI E-utilities fetch helper for the ClinVar proxy endpoints.
//
// The `_` prefix keeps Vercel from treating this file as a routable serverless
// function — it is imported by the sibling clinvar-*.js handlers only.
//
// Why this exists: E-utilities rate-limits anonymous callers to 3 requests/sec
// (10/sec with an NCBI_API_KEY) and answers 429 "Too Many Requests" once the
// limit is tripped. A single variant lookup in the UI fans out to several
// ClinVar endpoints at once (region + variant + protein, each doing an
// esearch/esummary pair), so bursts routinely exceed the limit — and because
// they share one Vercel egress IP, so do concurrent users. Without retry the
// 429 surfaces to the browser as a 502 and the affected card/section silently
// disappears (e.g. the "Same protein change" list intermittently not showing).
//
// Retry transient failures (429 + retryable 5xx + network/timeout errors) with
// exponential backoff and jitter, honouring a Retry-After header when present.

const RETRYABLE_STATUS = new Set([429, 500, 502, 503, 504]);
const MAX_ATTEMPTS = 4;

function sleep(ms) {
    return new Promise((resolve) => setTimeout(resolve, ms));
}

// Query-string fragment for the optional NCBI API key (raises the rate limit).
export function ncbiApiKeyParam() {
    return process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';
}

// Fetch an E-utilities URL, retrying transient failures. Returns the Response
// (only ever an ok one); throws on a non-retryable status or after exhausting
// the retry budget.
export async function ncbiFetchResponse(url, { timeout = 12000, accept = 'application/json' } = {}) {
    let lastErr = null;
    for (let attempt = 1; attempt <= MAX_ATTEMPTS; attempt++) {
        let res = null;
        try {
            res = await fetch(url, { headers: { Accept: accept }, signal: AbortSignal.timeout(timeout) });
        } catch (err) {
            lastErr = err; // network error / timeout — always retryable
        }
        if (res) {
            if (res.ok) return res;
            // Fail fast on client errors that a retry cannot fix (400/404/...).
            if (!RETRYABLE_STATUS.has(res.status)) {
                throw new Error(`NCBI request failed: ${res.status}`);
            }
            lastErr = new Error(`NCBI request failed: ${res.status}`);
        }
        if (attempt === MAX_ATTEMPTS) break;
        // Exponential backoff with jitter; never wait below the Retry-After hint.
        let waitMs = 400 * 2 ** (attempt - 1) + Math.floor(Math.random() * 200);
        const retryAfter = res ? Number(res.headers.get('retry-after')) : NaN;
        if (Number.isFinite(retryAfter) && retryAfter > 0) {
            waitMs = Math.max(waitMs, Math.min(retryAfter * 1000, 5000));
        }
        await sleep(waitMs);
    }
    throw lastErr || new Error('NCBI request failed');
}

// Convenience wrapper: fetch (with retry) and parse a JSON E-utilities response.
export async function ncbiFetchJson(url, opts = {}) {
    const res = await ncbiFetchResponse(url, { accept: 'application/json', ...opts });
    return res.json();
}
