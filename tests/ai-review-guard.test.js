/*
 * Tests for the AI-review proxy guards (helpers copied verbatim from
 * api/ai-review.js — KEEP IN SYNC). Covers the Origin allowlist that protects the
 * owner's OpenRouter key and the BYO-key format check.
 *
 * Run with: node tests/ai-review-guard.test.js
 */

// --- helpers copied verbatim from api/ai-review.js ------------------------

const DEFAULT_ALLOWED_ORIGINS = [
    'https://drdoubleb.com',
    'https://www.drdoubleb.com',
    'https://variant-search-for-chat-gpt.vercel.app',
    '*.vercel.app'
];

function parseAllowedOrigins() {
    const env = process.env.AI_ALLOWED_ORIGINS;
    if (env && env.trim()) return env.split(',').map((s) => s.trim()).filter(Boolean);
    return DEFAULT_ALLOWED_ORIGINS;
}

function isOriginAllowed(origin) {
    if (!origin) return true; // no browser Origin → handled by rate-limit/Turnstile
    const allowed = parseAllowedOrigins();
    if (allowed.includes('*')) return true;
    let host;
    try { host = new URL(origin).host; } catch { return false; }
    return allowed.some((entry) => {
        if (entry === origin) return true;
        if (entry.startsWith('*.')) return host === entry.slice(2) || host.endsWith(entry.slice(1));
        if (entry === 'localhost') return host === 'localhost' || host.startsWith('localhost:');
        try { return new URL(entry).host === host; } catch { return false; }
    });
}

function parseUserKey(value) {
    if (value === undefined || value === null || value === '') return '';
    const key = String(value).trim();
    if (!key) return '';
    return /^sk-[A-Za-z0-9_-]{16,}$/.test(key) ? key : null;
}

// --- assertions -----------------------------------------------------------

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; } else { failed++; console.error(`FAIL: ${name}`); }
}

// Origin allowlist — allowed
check('no Origin header is allowed', isOriginAllowed(undefined) === true);
check('empty Origin is allowed', isOriginAllowed('') === true);
check('live Hostinger origin allowed', isOriginAllowed('https://drdoubleb.com') === true);
check('www Hostinger origin allowed', isOriginAllowed('https://www.drdoubleb.com') === true);
check('production Vercel origin allowed', isOriginAllowed('https://variant-search-for-chat-gpt.vercel.app') === true);
check('vercel preview wildcard allowed', isOriginAllowed('https://variant-search-git-branch-abc.vercel.app') === true);

// Origin allowlist — blocked
check('unrelated origin blocked', isOriginAllowed('https://evil.com') === false);
check('suffix-spoof origin blocked', isOriginAllowed('https://drdoubleb.com.evil.com') === false);
check('non-dot vercel lookalike blocked', isOriginAllowed('https://evilvercel.app') === false);
check('malformed origin blocked', isOriginAllowed('not-a-url') === false);

// Origin allowlist — env override
process.env.AI_ALLOWED_ORIGINS = 'https://example.org';
check('env override allows configured origin', isOriginAllowed('https://example.org') === true);
check('env override blocks default origin', isOriginAllowed('https://drdoubleb.com') === false);
process.env.AI_ALLOWED_ORIGINS = '*';
check('wildcard * allows any origin', isOriginAllowed('https://anything.example') === true);
delete process.env.AI_ALLOWED_ORIGINS;

// BYO-key format
check('missing key → empty string (use owner key)', parseUserKey(undefined) === '');
check('empty key → empty string', parseUserKey('') === '');
check('whitespace key → empty string', parseUserKey('   ') === '');
check('valid OpenRouter key accepted', parseUserKey('sk-or-v1-0123456789abcdef0123') === 'sk-or-v1-0123456789abcdef0123');
check('valid key is trimmed', parseUserKey('  sk-or-v1-0123456789abcdef0123  ') === 'sk-or-v1-0123456789abcdef0123');
check('too-short key → null (malformed)', parseUserKey('sk-short') === null);
check('non-sk key → null (malformed)', parseUserKey('bearer-abcdefghijklmnop') === null);

console.log(`\nAI-review guard tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
