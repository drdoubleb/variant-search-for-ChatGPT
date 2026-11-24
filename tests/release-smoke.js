// Run after deploying the backend and copying the static assets.
// Browser tooling is the same as in tests/browser-smoke.js.
import { createRequire } from 'node:module';
import { checkM237lRelease, checkReleaseAssets, M237L_QUERY } from './support/m237l-release.js';

const require = createRequire(import.meta.url);
const { chromium } = require('playwright-core');
const address = new URL(process.argv[2] || 'https://www.drdoubleb.com/variant/');
address.searchParams.set('variant', M237L_QUERY);
const browser = await chromium.launch({
    executablePath: process.env.CHROMIUM_PATH,
    headless: true,
    args: ['--no-sandbox']
});
try {
    const page = await browser.newPage({ viewport: { width: 390, height: 844 } });
    const response = await page.goto(address.href, { waitUntil: 'domcontentloaded', timeout: 30000 });
    if (!response?.ok()) throw new Error(`The site returned HTTP ${response?.status() ?? 'unknown'}.`);
    await checkReleaseAssets(page);
    await checkM237lRelease(page);
    console.log('Release smoke passed. Both assets match this checkout and M237L resolves to MutationView record 3082 with the new evidence sections.');
} catch (err) {
    console.error(`Release smoke failed. ${err.message}`);
    process.exitCode = 1;
} finally {
    await browser.close();
}
