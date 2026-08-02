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

// Repeated AFTER the context payload, because that payload can run to hundreds of
// thousands of tokens — without this, the output-format instruction is buried far
// from the end and weak models drift into prose or truncate mid-array. The old
// single-template prompt got this signal for free: its self-check was headed
// "Self-check before returning JSON" and named schema fields, so JSON was the last
// thing the model read. Making the self-check format-neutral removed that, so the
// JSON reminder is restated here explicitly.
const FINAL_REMINDER = `Reminder — output format for this response:
Return the entire response as a single valid JSON object matching the schema given above, with the keys pathogenicity, amp_tier, amp_tier_rationale, fda_approved_therapies, resistance_or_lack_of_benefit, clinical_trials, summary, and limitations.
Output nothing before or after the JSON object: no preamble, no explanation, no markdown code fences.
Close every array and object. Keep the content concise so the response completes within the token limit rather than being cut off mid-array.`;

export default [CORE_GUIDANCE, OUTPUT_FORMAT, SELF_CHECK, CONTEXT_BLOCK, FINAL_REMINDER].join('\n\n');
