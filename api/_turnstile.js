// Cloudflare Turnstile verification for the AI-review proxy.
//
// Graceful by design: if TURNSTILE_SECRET_KEY is not set, verification is skipped
// (returns { ok: true, skipped: true }) so the endpoint works before Turnstile is
// provisioned. Once the secret is set, a missing/invalid token is rejected.
//
// Set in the deployment env:
//   TURNSTILE_SECRET_KEY  (Cloudflare Turnstile secret)
// and put the matching site key in the frontend widget.

const SITEVERIFY_URL = 'https://challenges.cloudflare.com/turnstile/v0/siteverify';

// Returns { ok: true } / { ok: true, skipped: true } when allowed,
// or { ok: false, reason } when the token is missing or fails verification.
export async function verifyTurnstile(token, ip) {
    const secret = process.env.TURNSTILE_SECRET_KEY;
    if (!secret) return { ok: true, skipped: true }; // not configured → skip

    if (!token || typeof token !== 'string') return { ok: false, reason: 'missing-token' };

    try {
        const form = new URLSearchParams();
        form.append('secret', secret);
        form.append('response', token);
        if (ip) form.append('remoteip', ip);

        const resp = await fetch(SITEVERIFY_URL, { method: 'POST', body: form });
        const data = await resp.json().catch(() => ({}));
        if (data && data.success) return { ok: true };
        return { ok: false, reason: (data && data['error-codes'] ? data['error-codes'].join(',') : '') || 'verification-failed' };
    } catch (err) {
        // Fail open on a Cloudflare outage — the rate limiter still bounds spend, and
        // blocking all reviews because siteverify is unreachable is worse for users.
        console.warn('[ai-review] turnstile verify error, allowing:', err?.message);
        return { ok: true, skipped: true };
    }
}
