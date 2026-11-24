// Human-readable prompt — used by the "Copy prompt" button, so that what the user
// pastes into their own LLM asks for a readable clinical write-up rather than the
// strict JSON the site parses into cards. The JSON counterpart lives in
// ./_ai-review-prompt.js.
//
// This prompt is fully self-contained and never mentions JSON output, schemas, or
// arrays: earlier versions appended a "disregard the JSON instructions" override to
// the JSON prompt, and the resulting contradiction confused weaker models. Keep it
// that way — if you add a format instruction here, make sure nothing in
// ./_ai-review-prompt-core.js contradicts it.
//
// Both prompts share their clinical guidance verbatim from
// ./_ai-review-prompt-core.js. Only the OUTPUT FORMAT section below is specific to
// this prompt, so tiering/evidence edits belong in the core file, not here.
//
// Placeholders substituted at request time by api/ai-review.js:
//   {{USER_NOTES_SECTION}} — optional user notes block (or empty)   [from core]
//   {{CONTEXT_JSON}}       — clamped variant context payload        [from core]
//
// Avoid using backticks, dollar-brace interpolation, or backslashes inside the
// string below, or escape them if you need to — they are special inside template
// literals.

import { CORE_GUIDANCE, SELF_CHECK, CONTEXT_BLOCK } from './_ai-review-prompt-core.js';

const OUTPUT_FORMAT = `Output format:

Write a clear, human-readable clinical interpretation for a knowledgeable clinical or scientific reader. Use Markdown, with a short heading per section and bullet points within sections. Do not return a machine-readable data structure — this answer is for a person to read.

Use these sections in order.

1. Germline pathogenicity — give the classification, the evaluated condition, and the independent rationale.
2. Somatic oncogenicity — give the classification and the supporting rationale.
3. Tumor-specific clinical significance — name the submitted tumor type and the assigned AMP/ASCO/CAP tier. Give the evidence level for an assessed Tier I or II. Use Insufficient disease context or Insufficient evidence when applicable.
4. FDA-approved therapies — give the drug, indication, biomarker, and supporting evidence for each finding.
5. Resistance or lack of benefit — give the therapy, tumor type, biomarker, and evidence with a source for each finding.
6. Clinical trials — give the NCT ID, title, phase, intervention, relevance, and URL for each lead.
7. Summary — briefly explain the variant and the most relevant evidence.
8. Limitations and verification needs — state what is uncertain and what sources or criteria need checking.

State plainly when the supplied evidence supports no findings for a section. Present the interpretation as a research and education summary that requires verification before clinical use.`;

// Repeated AFTER the context payload for the same reason as the JSON prompt: the
// payload can run to hundreds of thousands of tokens, leaving the output-format
// instruction far from the end. Phrased positively (write Markdown prose) rather
// than as a prohibition, so this prompt still never mentions JSON at all.
const FINAL_REMINDER = `Reminder — output format for this response:
Write the answer as Markdown prose for a human reader, using the eight sections in the order given above.
Do not return a machine-readable data structure. Keep each section concise.`;

export default [CORE_GUIDANCE, OUTPUT_FORMAT, SELF_CHECK, CONTEXT_BLOCK, FINAL_REMINDER].join('\n\n');
