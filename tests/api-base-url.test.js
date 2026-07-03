/*
 * Tests for backend API base-URL resolution (copied verbatim from script.js —
 * KEEP IN SYNC). The key behaviour: a Vercel preview/branch deployment should
 * call its OWN same-origin /api so a branch can be tested end-to-end before it
 * is merged to production, while Hostinger / custom domains / the production
 * apex keep using the pinned production backend.
 *
 * Run with: node tests/api-base-url.test.js
 */

const DEFAULT_BACKEND_API_BASE_URL = 'https://variant-search-for-chat-gpt.vercel.app';
const PRODUCTION_VERCEL_HOST = 'variant-search-for-chat-gpt.vercel.app';

function trimTrailingSlash(value) {
    return String(value || '').trim().replace(/\/+$/, '');
}

function isVercelPreviewHost() {
    const host = (typeof window !== 'undefined' && window.location && window.location.hostname) || '';
    return /\.vercel\.app$/i.test(host) && host !== PRODUCTION_VERCEL_HOST;
}

function getBackendApiBaseUrl() {
    if (typeof window.BACKEND_API_BASE_URL === 'string') {
        return trimTrailingSlash(window.BACKEND_API_BASE_URL);
    }
    if (isVercelPreviewHost()) return '';
    return DEFAULT_BACKEND_API_BASE_URL;
}

function getConfiguredApiEndpoint(globalName, apiPath) {
    const configured = String(window[globalName] || '').trim();
    if (configured) return configured;
    const base = getBackendApiBaseUrl();
    return base ? `${base}${apiPath}` : apiPath;
}

// --- harness --------------------------------------------------------------

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; } else { failed++; console.error(`  ✗ ${name}`); }
}
function setWindow(hostname, extra = {}) {
    global.window = { location: { hostname }, ...extra };
}

// Also mirror index.html's default assignment so we test the two halves together.
function applyIndexHtmlDefault() {
    if (typeof window.BACKEND_API_BASE_URL === 'undefined') {
        const h = (window.location && window.location.hostname) || '';
        const isPreview = /\.vercel\.app$/i.test(h) && h !== 'variant-search-for-chat-gpt.vercel.app';
        window.BACKEND_API_BASE_URL = isPreview ? '' : 'https://variant-search-for-chat-gpt.vercel.app';
    }
}

const PREVIEW_HOST = 'variant-search-for-chat-gpt-git-clau-847bf8-drdoublebs-projects.vercel.app';
const PROTEIN = '/api/clinvar-protein';

// 1. Preview host, no explicit config (script.js standalone) → same-origin.
setWindow(PREVIEW_HOST);
check('preview host → empty base (same origin)', getBackendApiBaseUrl() === '');
check('preview host → relative clinvar-protein endpoint', getConfiguredApiEndpoint('CLINVAR_PROTEIN_API_ENDPOINT', PROTEIN) === PROTEIN);

// 2. Production apex host, no explicit config → pinned production backend.
setWindow(PRODUCTION_VERCEL_HOST);
check('production apex → default backend', getBackendApiBaseUrl() === DEFAULT_BACKEND_API_BASE_URL);
check('production apex → absolute endpoint', getConfiguredApiEndpoint('CLINVAR_PROTEIN_API_ENDPOINT', PROTEIN) === DEFAULT_BACKEND_API_BASE_URL + PROTEIN);

// 3. Hostinger / custom domain → pinned production backend.
setWindow('drdoubleb.com');
check('custom domain → default backend', getBackendApiBaseUrl() === DEFAULT_BACKEND_API_BASE_URL);
check('custom domain → absolute endpoint', getConfiguredApiEndpoint('CLINVAR_PROTEIN_API_ENDPOINT', PROTEIN) === DEFAULT_BACKEND_API_BASE_URL + PROTEIN);

// 4. Explicit empty string means "same origin" and must not fall through to the default.
setWindow('drdoubleb.com', { BACKEND_API_BASE_URL: '' });
check('explicit empty base → same origin even off-Vercel', getBackendApiBaseUrl() === '');
check('explicit empty base → relative endpoint', getConfiguredApiEndpoint('CLINVAR_PROTEIN_API_ENDPOINT', PROTEIN) === PROTEIN);

// 5. Explicit custom base is respected and trailing-slash trimmed.
setWindow('example.org', { BACKEND_API_BASE_URL: 'https://my-backend.example/' });
check('explicit base is used and trimmed', getBackendApiBaseUrl() === 'https://my-backend.example');

// 6. index.html default assignment: preview host → '' → relative endpoints.
setWindow(PREVIEW_HOST);
applyIndexHtmlDefault();
check('index.html sets empty base on preview', window.BACKEND_API_BASE_URL === '');
check('index.html-derived civic endpoint is same-origin on preview',
    `${window.BACKEND_API_BASE_URL}/api/civic` === '/api/civic');

// 7. index.html default assignment: custom domain → pinned backend.
setWindow('drdoubleb.com');
applyIndexHtmlDefault();
check('index.html pins production backend off-Vercel', window.BACKEND_API_BASE_URL === DEFAULT_BACKEND_API_BASE_URL);

console.log(`\nAPI base-URL tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
