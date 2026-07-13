/*
 * Tests for ClinicalTrials.gov query building + filter helpers (copied verbatim from
 * api/clinicaltrials.js — KEEP IN SYNC). The proxy now pushes recruiting /
 * interventional / phase-2+ / US filters into the server query (percent-encoded), so
 * the page budget is spent on real candidates instead of being silently truncated at
 * 1000.
 *
 * Run with: node tests/clinicaltrials-query.test.js
 */

// --- helpers copied verbatim from api/clinicaltrials.js -------------------

const PHASES_PHASE2_PLUS = new Set(['PHASE2', 'PHASE3', 'PHASE4']);

function buildCtQuery(params) {
    return Object.entries(params)
        .map(([k, v]) => `${k}=${encodeURIComponent(v)}`)
        .join('&');
}

function isPhase2Plus(phases) {
    if (!Array.isArray(phases) || phases.length === 0) return false;
    return phases.some(p => PHASES_PHASE2_PLUS.has(String(p).toUpperCase()));
}

function hasUsLocation(locations) {
    if (!Array.isArray(locations)) return false;
    return locations.some(loc => String(loc.country || '').toLowerCase() === 'united states');
}

// --- assertions -----------------------------------------------------------

let passed = 0, failed = 0;
function check(name, cond) { if (cond) passed++; else { failed++; console.error(`FAIL: ${name}`); } }

// buildCtQuery percent-encodes values (spaces -> %20, brackets/colons/commas encoded)
// — the form the live v2 API was verified to accept for aggFilters/filter.advanced.
const q = buildCtQuery({
    'query.term': 'BRAF melanoma',
    'filter.overallStatus': 'RECRUITING',
    'aggFilters': 'studyType:int,phase:2 3 4',
    'filter.advanced': 'AREA[LocationCountry]United States'
});
check('spaces encoded as %20 (no + or raw space)', !q.includes(' ') && !q.includes('+') && q.includes('%20'));
check('aggFilters colon/comma encoded', q.includes('aggFilters=studyType%3Aint%2Cphase%3A2%203%204'));
check('filter.advanced brackets encoded', q.includes('filter.advanced=AREA%5BLocationCountry%5DUnited%20States'));
check('overallStatus present', q.includes('filter.overallStatus=RECRUITING'));
check('query.term present and encoded', q.includes('query.term=BRAF%20melanoma'));

// Phase filter (client-side backstop)
check('phase 2 passes', isPhase2Plus(['PHASE2']) === true);
check('phase 3 passes', isPhase2Plus(['PHASE1', 'PHASE3']) === true);
check('phase 1 only fails', isPhase2Plus(['PHASE1']) === false);
check('early phase 1 fails', isPhase2Plus(['EARLY_PHASE1']) === false);
check('no phases fails', isPhase2Plus([]) === false && isPhase2Plus(null) === false);

// US location (client-side backstop)
check('US location detected', hasUsLocation([{ country: 'United States' }]) === true);
check('mixed with US passes', hasUsLocation([{ country: 'France' }, { country: 'United States' }]) === true);
check('non-US fails', hasUsLocation([{ country: 'Canada' }]) === false);
check('empty/invalid fails', hasUsLocation([]) === false && hasUsLocation(null) === false);

console.log(`\nClinicalTrials query tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
