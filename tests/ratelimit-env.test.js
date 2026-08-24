/*
 * Tests for resolveRedisRestEnv (imported from api/_ratelimit.js). It must accept both the canonical UPSTASH_REDIS_REST_* names and the prefixed
 * <PREFIX>_REST_API_URL / <PREFIX>_REST_API_TOKEN pair that the Vercel/Upstash
 * Marketplace integration injects (any prefix), and must ignore a redis:// *_URL.
 *
 * Run with: node tests/ratelimit-env.test.js
 */

// --- real helper imported from api/_ratelimit.js --------------------------

import { resolveRedisRestEnv } from '../api/_ratelimit.js';

// --- assertions -----------------------------------------------------------

let passed = 0, failed = 0;
function check(name, cond) { if (cond) passed++; else { failed++; console.error(`FAIL: ${name}`); } }

const KEYS = [
    'UPSTASH_REDIS_REST_URL', 'UPSTASH_REDIS_REST_TOKEN',
    'STORAGE_REST_API_URL', 'STORAGE_REST_API_TOKEN', 'STORAGE_URL',
    'KV_REST_API_URL', 'KV_REST_API_TOKEN', 'KV_URL'
];
function clearEnv() { KEYS.forEach(k => { delete process.env[k]; }); }

// 1. Canonical names.
clearEnv();
process.env.UPSTASH_REDIS_REST_URL = 'https://canon.upstash.io';
process.env.UPSTASH_REDIS_REST_TOKEN = 'canon-token';
let r = resolveRedisRestEnv();
check('canonical UPSTASH_* resolved', r && r.url === 'https://canon.upstash.io' && r.token === 'canon-token');

// 2. STORAGE prefix (default Vercel Marketplace prefix).
clearEnv();
process.env.STORAGE_URL = 'redis://ignore-me';       // connection string, must be ignored
process.env.STORAGE_REST_API_URL = 'https://storage.upstash.io';
process.env.STORAGE_REST_API_TOKEN = 'storage-token';
r = resolveRedisRestEnv();
check('STORAGE_REST_API_* resolved', r && r.url === 'https://storage.upstash.io' && r.token === 'storage-token');
check('redis:// *_URL is not used as REST url', !r || r.url !== 'redis://ignore-me');

// 3. KV prefix.
clearEnv();
process.env.KV_REST_API_URL = 'https://kv.upstash.io';
process.env.KV_REST_API_TOKEN = 'kv-token';
r = resolveRedisRestEnv();
check('KV_REST_API_* resolved', r && r.url === 'https://kv.upstash.io' && r.token === 'kv-token');

// 4. Canonical takes priority over a prefixed pair.
clearEnv();
process.env.UPSTASH_REDIS_REST_URL = 'https://canon.upstash.io';
process.env.UPSTASH_REDIS_REST_TOKEN = 'canon-token';
process.env.STORAGE_REST_API_URL = 'https://storage.upstash.io';
process.env.STORAGE_REST_API_TOKEN = 'storage-token';
r = resolveRedisRestEnv();
check('canonical wins over prefixed', r && r.url === 'https://canon.upstash.io');

// 5. URL present but no matching token → not resolved.
clearEnv();
process.env.STORAGE_REST_API_URL = 'https://storage.upstash.io';
r = resolveRedisRestEnv();
check('url without matching token → null', r === null);

// 6. Nothing set → null (rate limiting stays disabled/graceful).
clearEnv();
r = resolveRedisRestEnv();
check('unset → null', r === null);

clearEnv();
console.log(`\nRate-limit env tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
