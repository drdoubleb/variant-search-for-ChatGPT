/*
 * Tests for the proxy guards: the shared Origin allowlist (api/_origin.js,
 * imported directly — no copy to drift) that protects the owner's OpenRouter
 * key and keeps third-party websites from using the data proxies as their
 * backend, plus the BYO-key format check (imported from
 * api/ai-review.js).
 *
 * Run with: node tests/ai-review-guard.test.js
 */

import { isOriginAllowed } from '../api/_origin.js';

// --- real helper imported from api/ai-review.js ---------------------------

import { parseUserKey } from '../api/ai-review.js';

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
check('project deploy-hash preview allowed', isOriginAllowed('https://variant-search-for-chat-gpt-a1b2c3d4e-drdoubleb.vercel.app') === true);
check('project git-branch preview allowed', isOriginAllowed('https://variant-search-for-chat-gpt-git-fix-liftover-drdoubleb.vercel.app') === true);

// Origin allowlist — blocked
check('unrelated origin blocked', isOriginAllowed('https://evil.com') === false);
check('suffix-spoof origin blocked', isOriginAllowed('https://drdoubleb.com.evil.com') === false);
check('non-dot vercel lookalike blocked', isOriginAllowed('https://evilvercel.app') === false);
check('someone else\'s Vercel deployment blocked', isOriginAllowed('https://attacker-frontend.vercel.app') === false);
check('project-prefix on another domain blocked', isOriginAllowed('https://variant-search-for-chat-gpt-x.evil.com') === false);
check('wildcard cannot span a dot', isOriginAllowed('https://variant-search-for-chat-gpt-x.y.vercel.app') === false);
check('malformed origin blocked', isOriginAllowed('not-a-url') === false);

// Origin allowlist — env override (passed as the second argument)
check('env override allows configured origin', isOriginAllowed('https://example.org', 'https://example.org') === true);
check('env override blocks default origin', isOriginAllowed('https://drdoubleb.com', 'https://example.org') === false);
check('wildcard * allows any origin', isOriginAllowed('https://anything.example', '*') === true);

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
