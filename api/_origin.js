// Shared browser-Origin allowlist for every serverless proxy.
//
// The `_` prefix keeps Vercel from treating this file as a routable serverless
// function — it is imported by the sibling handlers only.
//
// Why this exists: these endpoints were open CORS proxies — any third-party
// website's JavaScript could use this deployment (and its NCBI_API_KEY quota,
// Vercel invocation budget, and for ai-review the owner's OpenRouter key) as
// its backend. Browser requests are now restricted to the hosts that actually
// serve this frontend; requests with NO Origin header (curl, server-to-server
// integrations, GPT actions) are deliberately allowed through, so only
// browser-embedded third-party use is blocked.
//
// Entries may be full origins or host patterns with "*" wildcards (a wildcard
// matches within a single DNS label, so it cannot be stretched across
// domains); a bare "*" disables the check entirely.
//
// Env overrides (comma-separated):
//   AI_ALLOWED_ORIGINS     — used by api/ai-review.js
//   PROXY_ALLOWED_ORIGINS  — used by the data proxies; falls back to
//                            AI_ALLOWED_ORIGINS, then the defaults below.

export const DEFAULT_ALLOWED_ORIGINS = [
    'https://drdoubleb.com',
    'https://www.drdoubleb.com',
    'https://variant-search-for-chat-gpt.vercel.app',
    'variant-search-for-chat-gpt-*.vercel.app'
];

function parseAllowedOrigins(envValue) {
    if (envValue && envValue.trim()) return envValue.split(',').map((s) => s.trim()).filter(Boolean);
    return DEFAULT_ALLOWED_ORIGINS;
}

export function isOriginAllowed(origin, envValue) {
    if (!origin) return true; // no browser Origin → not a third-party website embedding
    const allowed = parseAllowedOrigins(envValue);
    if (allowed.includes('*')) return true;
    let parsedOrigin;
    try {
        parsedOrigin = new URL(origin);
        // Browser Origin headers contain exactly a scheme + host (+ optional port).
        // Reject path/userinfo variants and non-HTTP(S) protocols.
        if (parsedOrigin.origin !== origin || !['https:', 'http:'].includes(parsedOrigin.protocol)) return false;
    } catch { return false; }
    const host = parsedOrigin.host;
    const isHttps = parsedOrigin.protocol === 'https:';
    return allowed.some((entry) => {
        if (entry === origin) return true;
        if (entry.startsWith('*.')) return isHttps && (host === entry.slice(2) || host.endsWith(entry.slice(1)));
        // Host pattern with an embedded wildcard, e.g.
        // "variant-search-for-chat-gpt-*.vercel.app". The wildcard matches within
        // a single DNS label (no dots), so it cannot be stretched across domains.
        if (entry.includes('*')) {
            const escaped = entry.split('*')
                .map((part) => part.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'))
                .join('[a-z0-9-]*');
            return isHttps && new RegExp(`^${escaped}// Shared browser-Origin allowlist for every serverless proxy.
//
// The `_` prefix keeps Vercel from treating this file as a routable serverless
// function — it is imported by the sibling handlers only.
//
// Why this exists: these endpoints were open CORS proxies — any third-party
// website's JavaScript could use this deployment (and its NCBI_API_KEY quota,
// Vercel invocation budget, and for ai-review the owner's OpenRouter key) as
// its backend. Browser requests are now restricted to the hosts that actually
// serve this frontend; requests with NO Origin header (curl, server-to-server
// integrations, GPT actions) are deliberately allowed through, so only
// browser-embedded third-party use is blocked.
//
// Entries may be full origins or host patterns with "*" wildcards (a wildcard
// matches within a single DNS label, so it cannot be stretched across
// domains); a bare "*" disables the check entirely.
//
// Env overrides (comma-separated):
//   AI_ALLOWED_ORIGINS     — used by api/ai-review.js
//   PROXY_ALLOWED_ORIGINS  — used by the data proxies; falls back to
//                            AI_ALLOWED_ORIGINS, then the defaults below.

export const DEFAULT_ALLOWED_ORIGINS = [
    'https://drdoubleb.com',
    'https://www.drdoubleb.com',
    'https://variant-search-for-chat-gpt.vercel.app',
    'variant-search-for-chat-gpt-*.vercel.app'
];

function parseAllowedOrigins(envValue) {
    if (envValue && envValue.trim()) return envValue.split(',').map((s) => s.trim()).filter(Boolean);
    return DEFAULT_ALLOWED_ORIGINS;
}

export function isOriginAllowed(origin, envValue) {
    if (!origin) return true; // no browser Origin → not a third-party website embedding
    const allowed = parseAllowedOrigins(envValue);
    if (allowed.includes('*')) return true;
    let parsedOrigin;
    try {
        parsedOrigin = new URL(origin);
        // Browser Origin headers contain exactly a scheme + host (+ optional port).
        // Reject path/userinfo variants and non-HTTP(S) protocols.
        if (parsedOrigin.origin !== origin || !['https:', 'http:'].includes(parsedOrigin.protocol)) return false;
    } catch { return false; }
    const host = parsedOrigin.host;
    const isHttps = parsedOrigin.protocol === 'https:';
    return allowed.some((entry) => {
        if (entry === origin) return true;
        if (entry.startsWith('*.')) return isHttps && (host === entry.slice(2) || host.endsWith(entry.slice(1)));
        // Host pattern with an embedded wildcard, e.g.
        // "variant-search-for-chat-gpt-*.vercel.app". The wildcard matches within
        // a single DNS label (no dots), so it cannot be stretched across domains.
        if (entry.includes('*')) {
            const escaped = entry.split('*')
                .map((part) => part.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'))
                .join('[a-z0-9-]*');
, 'i').test(host);
        }
        if (entry === 'localhost') return host === 'localhost' || host.startsWith('localhost:');
        // For explicit allowlist URLs, compare the full origin, not just its host.
        // In particular, http://drdoubleb.com must not match https://drdoubleb.com.
        try { return new URL(entry).origin === parsedOrigin.origin; } catch { return false; }
    });
}

// Guard for the data proxies: sends the 403 itself and returns false when the
// request must be rejected. Call after the OPTIONS/method checks.
export function rejectDisallowedOrigin(req, res) {
    const envValue = process.env.PROXY_ALLOWED_ORIGINS ?? process.env.AI_ALLOWED_ORIGINS;
    const origin = req && req.headers ? req.headers.origin : undefined;
    if (isOriginAllowed(origin, envValue)) return false;
    res.status(403).json({ error: 'Origin not allowed' });
    return true;
}
