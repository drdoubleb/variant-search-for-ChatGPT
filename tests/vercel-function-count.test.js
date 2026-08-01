/*
 * Guards the Vercel Serverless Function budget.
 *
 * Vercel turns every .js file in api/ into a Serverless Function, EXCEPT files
 * whose name starts with an underscore — those are treated as plain imported
 * modules. The Hobby plan caps a deployment at 12 functions, and exceeding it
 * fails the deployment with no local signal: `npm test` stays green and the
 * error only shows up in the Vercel build log.
 *
 * That is exactly how it broke once: adding two shared prompt modules to api/
 * as ai-review-prompt-core.js / -human.js pushed the count from 12 to 14 and
 * failed the deploy. They are now _ai-review-prompt-core.js / _-human.js.
 *
 * So: if you add a file to api/ that is NOT an HTTP endpoint, prefix it with an
 * underscore. If you add a genuine new endpoint and this test fails, you need to
 * remove an endpoint or move the project off the Hobby plan — bumping the limit
 * below without doing one of those will just move the failure to deploy time.
 *
 * Run with: node tests/vercel-function-count.test.js
 */

import { readdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const HOBBY_FUNCTION_LIMIT = 12;

const apiDir = join(dirname(fileURLToPath(import.meta.url)), '..', 'api');
const functions = readdirSync(apiDir)
    .filter((name) => name.endsWith('.js') && !name.startsWith('_'))
    .sort();

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; } else { failed++; console.error(`FAIL: ${name}`); }
}

check(
    `api/ has ${functions.length} Serverless Functions, within the Hobby limit of ${HOBBY_FUNCTION_LIMIT}`,
    functions.length <= HOBBY_FUNCTION_LIMIT
);

// Every counted file must actually be an endpoint: a default-exported handler.
// A shared module that forgot its underscore both burns a function slot and
// would 500 if anything requested its path.
for (const name of functions) {
    const mod = await import(join(apiDir, name));
    check(
        `api/${name} default-exports a handler (else prefix it with "_")`,
        typeof mod.default === 'function'
    );
}

if (failed > 0) {
    console.error(`\nCounted functions:\n  ${functions.join('\n  ')}`);
}
console.log(`\nVercel function-count tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
