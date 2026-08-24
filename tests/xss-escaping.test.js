/*
 * Tests for the HTML/URL escaping helpers used to neutralise untrusted third-party
 * text (ClinVar / CIViC / COSMIC / model output) before it is written via innerHTML.
 * Helpers imported from script.js.
 *
 * Run with: node tests/xss-escaping.test.js
 */

// safeUrl resolves relative URLs against window.location.href — provide a stub.
globalThis.window = { location: { href: 'https://drdoubleb.com/variant/' } };

// --- real helpers imported from script.js ---------------------------------

await import('../script.js');
const { escapeHtml, safeUrl } = globalThis.__variantSearchHelpers;

// --- assertions -----------------------------------------------------------

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; } else { failed++; console.error(`FAIL: ${name}`); }
}

// escapeHtml neutralises the characters that break out of text/attribute context.
check('script tag escaped', escapeHtml('<script>alert(1)</script>') === '&lt;script&gt;alert(1)&lt;/script&gt;');
check('img onerror escaped', escapeHtml('<img src=x onerror=alert(1)>') === '&lt;img src=x onerror=alert(1)&gt;');
check('attribute breakout escaped', escapeHtml('" onmouseover="alert(1)') === '&quot; onmouseover=&quot;alert(1)');
check('ampersand escaped', escapeHtml('A & B') === 'A &amp; B');
check('single quote escaped', escapeHtml("it's") === 'it&#39;s');
check('null/undefined → empty', escapeHtml(null) === '' && escapeHtml(undefined) === '');
check('numbers coerced', escapeHtml(42) === '42');
check('plain text unchanged', escapeHtml('Pathogenic') === 'Pathogenic');

// safeUrl allows only http(s)/mailto; everything else becomes ''.
check('https allowed', safeUrl('https://civicdb.org/x') === 'https://civicdb.org/x');
check('http allowed', safeUrl('http://example.org') === 'http://example.org');
check('mailto allowed', safeUrl('mailto:a@b.com') === 'mailto:a@b.com');
check('javascript: blocked', safeUrl('javascript:alert(1)') === '');
check('data: blocked', safeUrl('data:text/html,<script>alert(1)</script>') === '');
check('vbscript: blocked', safeUrl('vbscript:msgbox(1)') === '');
check('whitespace/case javascript blocked', safeUrl('  JavaScript:alert(1)  ') === '');
check('relative URL resolves to allowed scheme, returns raw', safeUrl('/variants/5/summary') === '/variants/5/summary');
check('empty → empty', safeUrl('') === '' && safeUrl(null) === '');

console.log(`\nXSS escaping tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
