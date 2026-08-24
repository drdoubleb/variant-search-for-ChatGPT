/*
 * Tests for TP53 backend helpers (imported from api/tp53.js):
 *   - parseCsv: single-pass parser that preserves quoted newlines
 *   - proteinCodonKey / sameProteinCodon: same-codon relative matching
 *   - coordinateHasPosition: exact integer-token match (no substring false positives)
 *   - the discovery URL extraction (anchored so InducedMutationView isn't mistaken
 *     for MutationView)
 *
 * Run with: node tests/tp53-parse.test.js
 */

// --- real helpers imported from api/tp53.js -------------------------------

import { parseCsv, sameProteinCodon, coordinateHasPosition, extractDatasetUrls } from '../api/tp53.js';

// --- assertions -----------------------------------------------------------

let passed = 0, failed = 0;
function check(name, cond) { if (cond) passed++; else { failed++; console.error(`FAIL: ${name}`); } }

// parseCsv — a quoted field containing a newline must NOT split its row.
const csv = 'a,b,c\n1,"line one\nline two",3\n4,"has ""quotes"" in it",6\n';
const parsed = parseCsv(csv);
check('two data rows parsed', parsed.length === 2);
check('embedded newline kept in field', parsed[0].b === 'line one\nline two');
check('column alignment preserved after embedded newline', parsed[0].a === '1' && parsed[0].c === '3');
check('doubled quotes unescaped', parsed[1].b === 'has "quotes" in it');
check('row after quoted-newline row intact', parsed[1].a === '4' && parsed[1].c === '6');

// parseCsv — plain rows and trailing newline handling.
check('plain csv parses', parseCsv('h1,h2\nx,y').length === 1);
check('blank trailing rows skipped', parseCsv('h1,h2\nx,y\n\n').length === 1);

// coordinateHasPosition — exact integer token, not substring.
check('exact coordinate matches', coordinateHasPosition('7578406', 7578406) === true);
check('substring does NOT match (7578406 vs 17578406)', coordinateHasPosition('17578406', 7578406) === false);
check('range start matches', coordinateHasPosition('7578406-7578490', 7578406) === true);
check('range end matches', coordinateHasPosition('7578406-7578490', 7578490) === true);
check('unrelated coordinate does not match', coordinateHasPosition('7578400', 7578406) === false);

// sameProteinCodon — same ref+position, different alt.
check('R175H vs R175C same codon', sameProteinCodon('R175H', 'R175C') === true);
check('R175H vs R175H same codon (also exact elsewhere)', sameProteinCodon('R175H', 'R175H') === true);
check('R175H vs R273H different codon', sameProteinCodon('R175H', 'R273H') === false);
check('R175fs vs R175H same codon', sameProteinCodon('R175fs', 'R175H') === true);
check('missing ref/number → not same codon', sameProteinCodon('', 'R175H') === false);

// discovery extraction — pick MutationView (not InducedMutationView) and resolve.
const html = `
  <a href="/static/data/InducedMutationView_r22.csv">induced</a>
  <a href="/static/data/MutationView_r22.csv">mutation</a>
  <a href="/static/data/TumorVariantDownload_r22.csv">tumor</a>
  <a href="/static/data/GermlineDownload_r22.csv">germline</a>
`;
const urls = extractDatasetUrls(html);
check('discovery picks MutationView, not InducedMutationView',
  urls[0] === 'https://tp53.cancer.gov/static/data/MutationView_r22.csv');
check('discovery resolves all three needed datasets', urls.length === 3);
check('discovery returns absolute URLs', urls.every(u => u.startsWith('https://tp53.cancer.gov/static/data/')));
check('discovery empty when no CSV links', extractDatasetUrls('<p>no data here</p>').length === 0);

console.log(`\nTP53 parse/match tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
