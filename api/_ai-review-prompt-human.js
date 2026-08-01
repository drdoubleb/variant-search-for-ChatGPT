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

Cover these sections, in this order:

1. Pathogenicity — the classification (Pathogenic, Likely Pathogenic, VUS, Likely Benign, or Benign) and a one-line justification.
2. AMP/ASCO/CAP tier — the tier (Tier IA, Tier IB, Tier IIC, Tier IID, Tier IIE (tentative), Tier III, or Tier IV) and a brief rationale referencing the evidence level and the submitted tumor type.
3. FDA-approved therapies — one bullet per therapy, giving the drug, indication, biomarker context, and evidence. If none are supported, say so explicitly.
4. Resistance or lack of benefit — one bullet per therapy, giving the therapy, tumor type, biomarker context, and evidence. If none are supported, say so explicitly.
5. Clinical trials — one bullet per trial, giving the NCT ID, title, phase, intervention, relevance, and URL. If none are relevant, say so explicitly.
6. Summary — a short paragraph on the variant interpretation and its clinical significance.
7. Limitations and verification needs — brief bullets on caveats and what must be confirmed.

Keep the clinical caution and disclaimers throughout: this is for research and education only, it is not medical advice, and it must be verified against curated databases, current FDA labeling, clinical guidelines, and trial eligibility criteria before any clinical use.`;

export default [CORE_GUIDANCE, OUTPUT_FORMAT, SELF_CHECK, CONTEXT_BLOCK].join('\n\n');
