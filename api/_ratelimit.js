// Distributed rate limiting for the AI-review proxy, backed by Upstash Redis.
//
// This module is intentionally graceful: if the Upstash environment variables are
// not set (or the optional @upstash packages are not installed), every check
// returns { ok: true, skipped: true } so the endpoint keeps working with no rate
// limiting rather than failing. Enable it by setting, in the deployment env:
//   UPSTASH_REDIS_REST_URL, UPSTASH_REDIS_REST_TOKEN
// Optional overrides (defaults in parens):
//   AI_RL_IP_PER_MIN (12), AI_RL_IP_PER_DAY (120), AI_RL_GLOBAL_PER_DAY (1500)
//
// Serverless instances are ephemeral, so per-instance counters cannot bound total
// spend across a fleet — a shared store (Redis) is required for the global daily
// circuit breaker to be meaningful.

let _limitersPromise = null;

function retryAfterSeconds(result) {
    const reset = Number(result?.reset || 0);
    if (!reset) return 60;
    // Upstash returns reset as an epoch-ms timestamp.
    return Math.max(1, Math.ceil((reset - Date.now()) / 1000));
}

async function buildLimiters() {
    const url = process.env.UPSTASH_REDIS_REST_URL;
    const token = process.env.UPSTASH_REDIS_REST_TOKEN;
    if (!url || !token) return null; // not configured → disabled

    // Import lazily so the endpoint runs even when these optional deps are absent.
    const { Ratelimit } = await import('@upstash/ratelimit');
    const { Redis } = await import('@upstash/redis');
    const redis = new Redis({ url, token });

    const ipPerMin = Math.max(1, Number(process.env.AI_RL_IP_PER_MIN || 12));
    const ipPerDay = Math.max(1, Number(process.env.AI_RL_IP_PER_DAY || 120));
    const globalPerDay = Math.max(1, Number(process.env.AI_RL_GLOBAL_PER_DAY || 1500));

    return {
        ipMinute: new Ratelimit({ redis, limiter: Ratelimit.slidingWindow(ipPerMin, '60 s'), prefix: 'airl:ipm', analytics: false }),
        ipDay: new Ratelimit({ redis, limiter: Ratelimit.slidingWindow(ipPerDay, '86400 s'), prefix: 'airl:ipd', analytics: false }),
        globalDay: new Ratelimit({ redis, limiter: Ratelimit.slidingWindow(globalPerDay, '86400 s'), prefix: 'airl:gd', analytics: false })
    };
}

function getLimiters() {
    if (_limitersPromise === null) {
        _limitersPromise = buildLimiters().catch((err) => {
            console.warn('[ai-review] rate limiter init failed, disabling:', err?.message);
            return null;
        });
    }
    return _limitersPromise;
}

// Returns { ok: true } / { ok: true, skipped: true } when allowed,
// or { ok: false, scope, retryAfter } when a limit is exceeded.
export async function checkRateLimits(ip) {
    const limiters = await getLimiters();
    if (!limiters) return { ok: true, skipped: true };

    const id = ip || 'unknown';
    try {
        // Per-IP checks first so an abusive client is stopped before it consumes the
        // shared global budget.
        const minute = await limiters.ipMinute.limit(id);
        if (!minute.success) return { ok: false, scope: 'ip-minute', retryAfter: retryAfterSeconds(minute) };

        const day = await limiters.ipDay.limit(id);
        if (!day.success) return { ok: false, scope: 'ip-day', retryAfter: retryAfterSeconds(day) };

        const global = await limiters.globalDay.limit('global');
        if (!global.success) return { ok: false, scope: 'global', retryAfter: retryAfterSeconds(global) };

        return { ok: true };
    } catch (err) {
        // Never let a rate-store outage take down the endpoint; log and allow.
        console.warn('[ai-review] rate limit check failed, allowing:', err?.message);
        return { ok: true, skipped: true };
    }
}
