/*
 * Tests for the two AI-review prompt templates.
 *
 * The site's "Run AI review" needs strict JSON; the "Copy prompt" button needs a
 * human-readable Markdown write-up. These are two separate templates that share
 * their clinical guidance from api/ai-review-prompt-core.js. The tests below guard
 * the two ways that arrangement can break:
 *   1. Drift — someone edits the tiering rules in one prompt but not the other.
 *   2. Contradiction — the human prompt regains a JSON-output instruction (the
 *      earlier design appended "disregard the JSON instructions" to the JSON
 *      prompt, which confused weaker models).
 *
 * Run with: node tests/ai-review-prompt.test.js
 */

import JSON_PROMPT from '../api/ai-review-prompt.js';
import HUMAN_PROMPT from '../api/ai-review-prompt-human.js';
import { CORE_GUIDANCE, SELF_CHECK, CONTEXT_BLOCK } from '../api/ai-review-prompt-core.js';

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; } else { failed++; console.error(`FAIL: ${name}`); }
}

const count = (haystack, needle) => haystack.split(needle).length - 1;

// --- shared clinical core: identical in both prompts -----------------------

check('JSON prompt opens with the shared core', JSON_PROMPT.startsWith(CORE_GUIDANCE));
check('human prompt opens with the shared core', HUMAN_PROMPT.startsWith(CORE_GUIDANCE));
check('JSON prompt includes the shared self-check', JSON_PROMPT.includes(SELF_CHECK));
check('human prompt includes the shared self-check', HUMAN_PROMPT.includes(SELF_CHECK));
check('JSON prompt ends with the context block', JSON_PROMPT.endsWith(CONTEXT_BLOCK));
check('human prompt ends with the context block', HUMAN_PROMPT.endsWith(CONTEXT_BLOCK));

// Spot-check that the tiering substance really is in both, not just byte-equal
// blobs — these are the rules most likely to be edited in one place only.
for (const [label, prompt] of [['JSON', JSON_PROMPT], ['human', HUMAN_PROMPT]]) {
    check(`${label} prompt defines Tier IIE`, prompt.includes('Tier IIE (tentative):'));
    check(`${label} prompt keeps the tissue-agnostic rule`, prompt.includes('Tissue-agnostic / all-solid-tumor approvals count as same-tumor evidence'));
    check(`${label} prompt keeps the KRAS G12C colorectal rule`, prompt.includes('KRAS G12C inhibitor + anti-EGFR combination therapy'));
    check(`${label} prompt keeps the Pathogenic floor of Tier IIE`, prompt.includes('must receive a minimum of Tier IIE'));
}

// --- placeholders ---------------------------------------------------------

check('JSON prompt has exactly one user-notes placeholder', count(JSON_PROMPT, '{{USER_NOTES_SECTION}}') === 1);
check('human prompt has exactly one user-notes placeholder', count(HUMAN_PROMPT, '{{USER_NOTES_SECTION}}') === 1);
check('JSON prompt has exactly one context placeholder', count(JSON_PROMPT, '{{CONTEXT_JSON}}') === 1);
check('human prompt has exactly one context placeholder', count(HUMAN_PROMPT, '{{CONTEXT_JSON}}') === 1);
check('JSON prompt has exactly one schema placeholder', count(JSON_PROMPT, '{{SCHEMA_JSON}}') === 1);
check('human prompt has no schema placeholder', !HUMAN_PROMPT.includes('{{SCHEMA_JSON}}'));

// No placeholder outside the set api/ai-review.js substitutes, or it would ship
// to the model as literal braces.
const KNOWN = new Set(['{{USER_NOTES_SECTION}}', '{{SCHEMA_JSON}}', '{{CONTEXT_JSON}}']);
for (const [label, prompt] of [['JSON', JSON_PROMPT], ['human', HUMAN_PROMPT]]) {
    const unknown = (prompt.match(/\{\{[^}]*\}\}/g) || []).filter((tok) => !KNOWN.has(tok));
    check(`${label} prompt has no unknown placeholders`, unknown.length === 0);
}

// --- output-format separation ---------------------------------------------

check('JSON prompt asks for strict JSON', JSON_PROMPT.includes('Return ONLY valid JSON matching this schema'));

// The human prompt must not mention JSON output at all. Its only legitimate use of
// the word is the "Context JSON:" header on the supplied data payload, so strip the
// context block before checking.
const humanWithoutContext = HUMAN_PROMPT.slice(0, HUMAN_PROMPT.length - CONTEXT_BLOCK.length);
check('human prompt never mentions JSON outside the context block', !humanWithoutContext.includes('JSON'));
check('human prompt never mentions a schema', !/schema/i.test(humanWithoutContext));
check('human prompt asks for Markdown', HUMAN_PROMPT.includes('Use Markdown'));

// Guard against reintroducing the self-contradicting override.
for (const phrase of ['Disregard every instruction above', 'supersedes any earlier formatting instruction', 'Do NOT output JSON']) {
    check(`human prompt does not contain override phrase: ${phrase}`, !HUMAN_PROMPT.includes(phrase));
}

// Each prompt states exactly one output format.
check('JSON prompt has one output-format section', count(JSON_PROMPT, 'Output format:') === 1);
check('human prompt has one output-format section', count(HUMAN_PROMPT, 'Output format:') === 1);

// --- template-literal hazards ---------------------------------------------
// The prompt files are template literals; unescaped backticks or ${...} would have
// thrown at import, but a stray backslash silently mangles the text.
for (const [label, prompt] of [['JSON', JSON_PROMPT], ['human', HUMAN_PROMPT]]) {
    check(`${label} prompt has no stray backslashes`, !prompt.includes('\\'));
}

console.log(`\nAI-review prompt tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
