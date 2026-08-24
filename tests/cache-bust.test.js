/*
 * Cache-bust guard: index.html's ?v= query strings must change whenever
 * script.js or style.css change, or Hostinger/browser caches keep serving the
 * stale copy. Remembering to bump the string by hand fails eventually, so this
 * test derives the expected value from the file content itself: each ?v= must
 * END with the first 10 hex chars of the file's SHA-256. Any prefix (a date,
 * for readability) is allowed before it.
 *
 * When this test fails it prints the exact value to paste into index.html.
 *
 * Run with: node tests/cache-bust.test.js
 */

import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const indexHtml = readFileSync(join(root, 'index.html'), 'utf8');

let passed = 0;
let failed = 0;

function checkAsset(file) {
    const hash = createHash('sha256').update(readFileSync(join(root, file))).digest('hex').slice(0, 10);
    const m = indexHtml.match(new RegExp(`${file.replace('.', '\\.')}\\?v=([A-Za-z0-9.-]+)`));
    const current = m ? m[1] : '(no ?v= found)';
    const today = new Date().toISOString().slice(0, 10).replace(/-/g, '');
    if (m && m[1].endsWith(hash)) {
        passed++;
        console.log(`ok   ${file} ?v= matches content hash (${current})`);
    } else {
        failed++;
        console.error(`FAIL ${file} changed but index.html still loads ?v=${current}`);
        console.error(`     set it to:  ${file}?v=${today}-${hash}`);
    }
}

checkAsset('script.js');
checkAsset('style.css');

console.log(`\nCache-bust tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
