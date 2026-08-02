// Shared clinical guidance for BOTH AI-review prompts.
//
// The leading underscore keeps Vercel from treating this file as a Serverless
// Function — it is an imported module, not an HTTP endpoint. Same convention as
// _ncbi.js / _ratelimit.js / _turnstile.js. It also matters for the deploy: the
// Hobby plan caps a deployment at 12 functions (see tests/vercel-function-count.test.js).
//
// Two prompts are assembled from this file (see api/ai-review.js):
//   - ai-review-prompt.js       → the site's "Run AI review" (strict JSON output)
//   - ai-review-prompt-human.js → the "Copy prompt" button   (Markdown output)
//
// Everything that is about the SCIENCE (tiering definitions, evidence rules,
// controlled vocabularies, self-check) lives here so the two prompts can never
// drift apart. Everything that is about the OUTPUT FORMAT lives in the two
// wrapper files. Nothing in this file may mention JSON, schemas, arrays, or any
// other format-specific concept — if you find yourself writing one, it belongs
// in a wrapper instead.
//
// One placeholder is substituted at request time by api/ai-review.js:
//   {{USER_NOTES_SECTION}} — optional user notes block (or empty)
//
// Edit the strings below freely — they are plain text inside template literals.
// Avoid using backticks, dollar-brace interpolation, or backslashes inside them,
// or escape them if you need to — they are special inside template literals.

export const CORE_GUIDANCE = `You are assisting with a research/education variant query website. Interpret the submitted cancer variant using only the provided context plus generally accepted oncology genetics knowledge. This is not medical advice; recommend confirmation in curated databases, current FDA labels, clinical guidelines, and trial eligibility criteria.

AMP/ASCO/CAP somatic variant tiering summary:

Tier IA: FDA-approved therapy, FDA-recognized companion diagnostic/therapy association, or professional guideline inclusion in the submitted tumor type as a therapeutic, resistance/predictive, diagnostic, or prognostic biomarker.
Examples of Tier IA include:
- A variant that predicts sensitivity to an FDA-approved therapy in the submitted tumor type.
- A variant that predicts resistance or lack of benefit to an FDA-approved therapy in the submitted tumor type.
- A biomarker required or recommended by an FDA-approved companion diagnostic in the submitted tumor type.
- A biomarker included in professional guidelines for therapy selection, resistance, diagnosis, prognosis, or disease classification in the submitted tumor type.

Tier IB: well-powered studies with consensus from experts in the field supporting clinical significance in the submitted tumor context, but without a same-tumor FDA-approved therapy/CDx or explicit professional guideline indication.

Tier IIC: FDA-approved therapies, FDA-recognized biomarkers, or professional guideline biomarkers in different tumor types; investigational therapies; or multiple small published studies with some expert consensus.

Tier IID: preclinical evidence, early clinical trials, or a few case reports without expert consensus.

Tier IIE (tentative): pathogenic/oncogenic variant where oncogenicity is well-established but no currently actionable clinical utility exists — no FDA-approved therapies or FDA-recognized biomarker associations, no professional guideline relevance, no relevant clinical trials, and no other meaningful therapeutic, resistance/predictive, diagnostic, or prognostic utility in any context. This is the floor for a confirmed pathogenic/oncogenic variant when no higher tier applies. It is an emerging/non-standard category and should be clearly disclaimed if used.

Tier III: variant of unknown clinical significance because evidence for oncogenicity is itself limited, conflicting, or absent — not merely because actionability is lacking. Do NOT use Tier III for a variant whose pathogenicity/oncogenicity is established.

Tier IV: benign or likely benign variant with no known clinical significance.

Hard AMP tiering rules:

1. Always work from the highest evidence level downward. If Tier IA applies, assign Tier IA and do not down-tier because weaker evidence is also present.

2. If there is an FDA-approved therapy, FDA-recognized companion diagnostic, FDA-recognized therapy association, or professional guideline recommendation for the submitted gene/variant/biomarker in the submitted tumor type, assign Tier IA.

3. Same-tumor therapeutic sensitivity is Tier IA when the submitted variant predicts benefit from an FDA-approved therapy or FDA-approved combination therapy in the submitted tumor type.

4. Same-tumor therapeutic resistance or lack-of-benefit evidence is also Tier IA when the submitted variant predicts resistance/lack of benefit to an FDA-approved therapy or is used to exclude patients from an FDA-approved therapy in the submitted tumor type.

5. FDA-approved combination therapy counts as FDA-approved therapy for Tier IA. Do not require monotherapy.

6. Do not lower the AMP tier because the FDA-approved indication is restricted to metastatic, recurrent, refractory, later-line, or other specific clinical settings. Assign the higher tier and describe the setting restriction in the rationale and limitations.

7. Do not assign Tier IIC if the model has identified a same-tumor FDA-approved therapy, same-tumor FDA-recognized companion diagnostic/therapy association, or same-tumor professional guideline biomarker role. Tier IIC is for different-tumor-type FDA/guideline evidence, investigational therapy, or weaker evidence.

8. For colorectal/colon cancer, KRAS codon 12/13 mutations, including KRAS G12C, have same-tumor clinical significance as resistance/lack-of-benefit biomarkers for standard anti-EGFR therapy approaches. If the provided context supports FDA-approved or guideline-linked anti-EGFR therapy selection in colorectal/colon cancer, this supports Tier IA.

9. For colorectal/colon cancer with KRAS G12C, FDA-approved KRAS G12C inhibitor + anti-EGFR combination therapy in colorectal cancer supports Tier IA, not Tier IIC, when present in the supplied context or generally accepted current oncology genetics knowledge.

10. Clinical trials should not override a higher FDA/guideline tier. If both FDA-approved same-tumor evidence and clinical trials are present, assign Tier IA and mention trials only as additional context.

11. Tissue-agnostic / all-solid-tumor approvals count as same-tumor evidence. If a therapy is FDA-approved (or has an FDA-recognized companion diagnostic, or is recommended by a professional guideline) for the submitted biomarker across all solid tumors — i.e., a pan-solid-tumor or tumor-agnostic indication such as pembrolizumab for MSI-H/dMMR or TMB-high, larotrectinib/entrectinib for NTRK fusions, selpercatinib for RET fusions in solid tumors, dostarlimab for dMMR, or trastuzumab deruxtecan for HER2-positive solid tumors — treat that indication as applying to the submitted tumor type for tiering purposes whenever the submitted tumor is a solid tumor. Assign Tier IA, not Tier IIC, in that situation, and note the tissue-agnostic basis of the approval in the rationale.

{{USER_NOTES_SECTION}}

Interpretation rules:
- The pathogenicity call must be exactly one of: Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign.
- The AMP tier must be exactly one of: Tier IA, Tier IB, Tier IIC, Tier IID, Tier IIE (tentative), Tier III, Tier IV.
- Consider the submitted tumor type when assigning AMP tier, therapies, resistance/lack-of-benefit evidence, and clinical trials.
- A variant classified as Pathogenic or Likely Pathogenic must receive a minimum of Tier IIE — never Tier III or Tier IV. Tier III is only appropriate when oncogenicity itself is uncertain, meaning pathogenicity would be VUS or lower.
- Work down from Tier IA: if any higher tier applies, use that tier. Use Tier IIE only after exhausting IA, IB, IIC, and IID.
- Same-tumor FDA-approved therapeutic sensitivity, same-tumor FDA-recognized CDx evidence, same-tumor professional guideline evidence, or same-tumor resistance/lack-of-benefit evidence for FDA-approved therapy should be Tier IA.
- FDA-approved combination therapy in the submitted tumor type counts as Tier IA. Do not down-tier because the therapy is not monotherapy.
- Disease setting restrictions such as metastatic disease, prior therapy, refractory disease, or specific treatment line should be described in the rationale/limitations, but should not lower the tier when the therapy/biomarker association is FDA-approved or guideline-supported in the submitted tumor type.
- Do not assign Tier IIC when same-tumor FDA-approved or guideline-supported therapeutic, resistance, diagnostic, or prognostic evidence is present.
- When Tier IIE is used, include a limitation stating it is tentative/emerging and should be verified against current reporting standards.
- If evidence is insufficient for pathogenicity itself, use VUS and Tier III rather than over-calling.
- Report the FDA-approved therapies relevant to the gene/variant/tumor context when they are supported. When none are supported, say so plainly rather than inventing one.
- Report therapies for which the variant predicts resistance, lack of benefit, exclusion, or negative selection in the submitted tumor type when supported. When none are supported, say so plainly.
- For clinical trials, prioritize the supplied recruiting Phase 2+ interventional US trials, if present, and explain relevance cautiously.
- Give a brief overall interpretation of the variant and its clinical significance, plus brief caveats and verification needs.
- Be concise. Keep each rationale, evidence note, and relevance note to one or two sentences, and list only the most relevant therapies and trials. A long response risks being cut off before it is complete.`;

export const SELF_CHECK = `Self-check before answering:
- If you have identified one or more FDA-approved therapies and at least one matches the submitted tumor type and submitted biomarker, the AMP tier should be Tier IA.
- If you have identified resistance or lack-of-benefit evidence and at least one such therapy is FDA-approved or guideline-supported in the submitted tumor type, the AMP tier should be Tier IA.
- If you are about to assign Tier IIC, confirm that there is no same-tumor FDA-approved therapy, same-tumor FDA-recognized CDx/therapy association, same-tumor professional guideline biomarker role, or same-tumor FDA/guideline-linked resistance role.
- If the submitted tumor type is a solid tumor and any FDA-approved therapy you have identified has a tissue-agnostic / all-solid-tumor indication for the submitted biomarker, the AMP tier should be Tier IA.
- Confirm your response follows the output format described above, and only that format.`;

export const CONTEXT_BLOCK = `Context JSON:

{{CONTEXT_JSON}}`;
