// JSON prompt — used by the site's own "Run AI review", whose response is parsed
// into the result cards. The human-readable counterpart used by the "Copy prompt"
// button lives in ./_ai-review-prompt-human.js.
//
// Both prompts share their clinical guidance verbatim from
// ./_ai-review-prompt-core.js. Only the OUTPUT FORMAT section below is specific to
// this prompt, so tiering/evidence edits belong in the core file, not here.
//
// Placeholders substituted at request time by api/ai-review.js:
//   {{USER_NOTES_SECTION}} — optional user notes block (or empty)   [from core]
//   {{SCHEMA_JSON}}        — JSON schema the model must conform to  [below]
//   {{CONTEXT_JSON}}       — clamped variant context payload        [from core]
//
// Avoid using backticks, dollar-brace interpolation, or backslashes inside the
// string below, or escape them if you need to — they are special inside template
// literals.

import { CORE_GUIDANCE, SELF_CHECK, CONTEXT_BLOCK } from './_ai-review-prompt-core.js';

const OUTPUT_FORMAT = `Output format:

Return ONLY valid JSON matching this schema; do not wrap it in markdown:

{{SCHEMA_JSON}}

- Use one of the exact pathogenicity strings listed above as the value of "pathogenicity", and one of the exact AMP tier strings as the value of "amp_tier".
- "amp_tier_rationale" is a brief rationale covering the tumor-type considerations and the evidence level supporting the selected tier.
- "fda_approved_therapies", "resistance_or_lack_of_benefit", and "clinical_trials" are arrays. Where the interpretation rules above say to state plainly that nothing was identified, return an empty array for that field instead and explain in "summary" or "limitations".
- "summary" is a brief string; "limitations" is an array of brief strings.`;

export default [CORE_GUIDANCE, OUTPUT_FORMAT, SELF_CHECK, CONTEXT_BLOCK].join('\n\n');
