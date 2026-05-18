// Serverless proxy for optional OpenRouter AI variant interpretation.
// POST /api/ai-review
// Requires OPENROUTER_API_KEY in the Vercel environment.

const OPENROUTER_API_URL = 'https://openrouter.ai/api/v1/chat/completions';
const DEFAULT_MODEL = 'openai/gpt-4.1-mini';
const ALLOWED_MODELS = new Set([
    'openai/gpt-4.1-mini',
    'openai/gpt-oss-120b',
    'openai/gpt-4o-mini',
    'openai/gpt-5-mini',
    'openai/gpt-5.4-nano',
    'openai/gpt-5-nano',
    'google/gemini-2.5-flash-lite',
    'google/gemini-3-flash-preview',
    'google/gemini-3.1-flash-lite',
    'anthropic/claude-3-haiku',
    'deepseek/deepseek-v4-flash',
    'deepseek/deepseek-v4-pro',
    'x-ai/grok-4.3',
    'minimax/minimax-m2.7',
    'nvidia/nemotron-3-super-120b-a12b:free'
]);

function clampString(value, maxLength = 20000) {
    const text = typeof value === 'string' ? value : JSON.stringify(value ?? null);
    return text.length > maxLength ? `${text.slice(0, maxLength)}\n...[truncated]` : text;
}

function buildPrompt(context) {
    const ampGuidelines = [
        'AMP/ASCO/CAP somatic variant tiering summary:',
        '',
        'Tier IA: FDA-approved therapy, FDA-recognized companion diagnostic/therapy association, or professional guideline inclusion in the submitted tumor type as a therapeutic, resistance/predictive, diagnostic, or prognostic biomarker.',
        'Examples of Tier IA include:',
        '- A variant that predicts sensitivity to an FDA-approved therapy in the submitted tumor type.',
        '- A variant that predicts resistance or lack of benefit to an FDA-approved therapy in the submitted tumor type.',
        '- A biomarker required or recommended by an FDA-approved companion diagnostic in the submitted tumor type.',
        '- A biomarker included in professional guidelines for therapy selection, resistance, diagnosis, prognosis, or disease classification in the submitted tumor type.',
        '',
        'Tier IB: well-powered studies with consensus from experts in the field supporting clinical significance in the submitted tumor context, but without a same-tumor FDA-approved therapy/CDx or explicit professional guideline indication.',
        '',
        'Tier IIC: FDA-approved therapies, FDA-recognized biomarkers, or professional guideline biomarkers in different tumor types; investigational therapies; or multiple small published studies with some expert consensus.',
        '',
        'Tier IID: preclinical evidence, early clinical trials, or a few case reports without expert consensus.',
        '',
        'Tier IIE (tentative): pathogenic/oncogenic variant where oncogenicity is well-established but no currently actionable clinical utility exists — no FDA-approved therapies or FDA-recognized biomarker associations, no professional guideline relevance, no relevant clinical trials, and no other meaningful therapeutic, resistance/predictive, diagnostic, or prognostic utility in any context. This is the floor for a confirmed pathogenic/oncogenic variant when no higher tier applies. It is an emerging/non-standard category and should be clearly disclaimed if used.',
        '',
        'Tier III: variant of unknown clinical significance because evidence for oncogenicity is itself limited, conflicting, or absent — not merely because actionability is lacking. Do NOT use Tier III for a variant whose pathogenicity/oncogenicity is established.',
        '',
        'Tier IV: benign or likely benign variant with no known clinical significance.'
    ].join('\n');

    const hardTierRules = [
        'Hard AMP tiering rules:',
        '',
        '1. Always work from the highest evidence level downward. If Tier IA applies, assign Tier IA and do not down-tier because weaker evidence is also present.',
        '',
        '2. If there is an FDA-approved therapy, FDA-recognized companion diagnostic, FDA-recognized therapy association, or professional guideline recommendation for the submitted gene/variant/biomarker in the submitted tumor type, assign Tier IA.',
        '',
        '3. Same-tumor therapeutic sensitivity is Tier IA when the submitted variant predicts benefit from an FDA-approved therapy or FDA-approved combination therapy in the submitted tumor type.',
        '',
        '4. Same-tumor therapeutic resistance or lack-of-benefit evidence is also Tier IA when the submitted variant predicts resistance/lack of benefit to an FDA-approved therapy or is used to exclude patients from an FDA-approved therapy in the submitted tumor type.',
        '',
        '5. FDA-approved combination therapy counts as FDA-approved therapy for Tier IA. Do not require monotherapy.',
        '',
        '6. Do not lower the AMP tier because the FDA-approved indication is restricted to metastatic, recurrent, refractory, later-line, or other specific clinical settings. Assign the higher tier and describe the setting restriction in the rationale and limitations.',
        '',
        '7. Do not assign Tier IIC if the model has identified a same-tumor FDA-approved therapy, same-tumor FDA-recognized companion diagnostic/therapy association, or same-tumor professional guideline biomarker role. Tier IIC is for different-tumor-type FDA/guideline evidence, investigational therapy, or weaker evidence.',
        '',
        '8. For colorectal/colon cancer, KRAS codon 12/13 mutations, including KRAS G12C, have same-tumor clinical significance as resistance/lack-of-benefit biomarkers for standard anti-EGFR therapy approaches. If the provided context supports FDA-approved or guideline-linked anti-EGFR therapy selection in colorectal/colon cancer, this supports Tier IA.',
        '',
        '9. For colorectal/colon cancer with KRAS G12C, FDA-approved KRAS G12C inhibitor + anti-EGFR combination therapy in colorectal cancer supports Tier IA, not Tier IIC, when present in the supplied context or generally accepted current oncology genetics knowledge.',
        '',
        '10. Clinical trials should not override a higher FDA/guideline tier. If both FDA-approved same-tumor evidence and clinical trials are present, assign Tier IA and mention trials only as additional context.'
    ].join('\n');

    const schema = {
        pathogenicity: 'Pathogenic | Likely Pathogenic | VUS | Likely Benign | Benign',
        amp_tier: 'Tier IA | Tier IB | Tier IIC | Tier IID | Tier IIE (tentative) | Tier III | Tier IV',
        amp_tier_rationale: 'brief rationale with tumor-type considerations and the evidence level supporting the selected AMP tier',
        fda_approved_therapies: [
            {
                drug: 'string',
                indication: 'string',
                biomarker_context: 'string',
                evidence: 'string'
            }
        ],
        resistance_or_lack_of_benefit: [
            {
                therapy: 'string',
                tumor_type: 'string',
                biomarker_context: 'string',
                evidence: 'string'
            }
        ],
        clinical_trials: [
            {
                nct_id: 'string',
                title: 'string',
                phase: 'string',
                intervention: 'string',
                relevance: 'string',
                url: 'string'
            }
        ],
        summary: 'brief variant interpretation and clinical significance',
        limitations: ['brief caveats and verification needs']
    };

    return [
        'You are assisting with a research/education variant query website. Interpret the submitted cancer variant using only the provided context plus generally accepted oncology genetics knowledge. This is not medical advice; recommend confirmation in curated databases, current FDA labels, clinical guidelines, and trial eligibility criteria.',
        ampGuidelines,
        hardTierRules,
        'Return ONLY valid JSON matching this schema; do not wrap it in markdown:',
        JSON.stringify(schema),
        'Rules:',
        '- Pathogenicity must be exactly one of: Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign.',
        '- AMP tier must be exactly one of: Tier IA, Tier IB, Tier IIC, Tier IID, Tier IIE (tentative), Tier III, Tier IV.',
        '- Consider tumor_type when assigning AMP tier, therapies, resistance/lack-of-benefit evidence, and clinical trials.',
        '- A variant classified as Pathogenic or Likely Pathogenic must receive a minimum of Tier IIE — never Tier III or Tier IV. Tier III is only appropriate when oncogenicity itself is uncertain, meaning pathogenicity would be VUS or lower.',
        '- Work down from Tier IA: if any higher tier applies, use that tier. Use Tier IIE only after exhausting IA, IB, IIC, and IID.',
        '- Same-tumor FDA-approved therapeutic sensitivity, same-tumor FDA-recognized CDx evidence, same-tumor professional guideline evidence, or same-tumor resistance/lack-of-benefit evidence for FDA-approved therapy should be Tier IA.',
        '- FDA-approved combination therapy in the submitted tumor type counts as Tier IA. Do not down-tier because the therapy is not monotherapy.',
        '- Disease setting restrictions such as metastatic disease, prior therapy, refractory disease, or specific treatment line should be described in the rationale/limitations, but should not lower the tier when the therapy/biomarker association is FDA-approved or guideline-supported in the submitted tumor type.',
        '- Do not assign Tier IIC when same-tumor FDA-approved or guideline-supported therapeutic, resistance, diagnostic, or prognostic evidence is present.',
        '- When Tier IIE is used, include a limitation stating it is tentative/emerging and should be verified against current reporting standards.',
        '- If evidence is insufficient for pathogenicity itself, use VUS and Tier III rather than over-calling.',
        '- fda_approved_therapies should include FDA-approved therapies relevant to the gene/variant/tumor context when supported; otherwise return an empty array and explain in summary/limitations.',
        '- resistance_or_lack_of_benefit should include therapies for which the variant predicts resistance, lack of benefit, exclusion, or negative selection in the submitted tumor type when supported; otherwise return an empty array.',
        '- Clinical trials should prioritize supplied recruiting Phase 2+ interventional US trials, if present, and explain relevance cautiously.',
        'Self-check before returning JSON:',
        '- If fda_approved_therapies is non-empty and at least one therapy matches the submitted tumor type and submitted biomarker, amp_tier should be Tier IA.',
        '- If resistance_or_lack_of_benefit is non-empty and at least one therapy is FDA-approved or guideline-supported in the submitted tumor type, amp_tier should be Tier IA.',
        '- If amp_tier is Tier IIC, confirm that there is no same-tumor FDA-approved therapy, same-tumor FDA-recognized CDx/therapy association, same-tumor professional guideline biomarker role, or same-tumor FDA/guideline-linked resistance role.',
        'Context JSON:',
        clampString(context, 120000)
    ].join('\n\n');
}

function parseModel(value) {
    const model = String(value || '').trim();
    return ALLOWED_MODELS.has(model) ? model : DEFAULT_MODEL;
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'POST,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'POST') return res.status(405).json({ error: 'Method not allowed' });

    const apiKey = process.env.OPENROUTER_API_KEY;
    if (!apiKey) {
        return res.status(500).json({ error: 'OpenRouter is not configured. Set OPENROUTER_API_KEY in the deployment environment.' });
    }

    const body = req.body || {};
    const model = parseModel(body.model);
    const context = body.context || {};

    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 45000);

    try {
        const response = await fetch(OPENROUTER_API_URL, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${apiKey}`,
                'HTTP-Referer': req.headers.origin || 'https://variant-search-for-chat-gpt.vercel.app',
                'X-Title': 'Variant Search AI Review'
            },
            body: JSON.stringify({
                model,
                messages: [
                    { role: 'system', content: 'You are a careful molecular oncology variant interpretation assistant. Return strictly valid JSON.' },
                    { role: 'user', content: buildPrompt(context) }
                ],
                temperature: 0.1,
                response_format: { type: 'json_object' }
            }),
            signal: controller.signal
        });

        const text = await response.text();
        let data;
        try { data = JSON.parse(text); } catch { data = { raw: text }; }

        if (!response.ok) {
            const detail = data?.error?.message || data?.error || text.slice(0, 500) || `OpenRouter returned ${response.status}`;
            return res.status(502).json({ error: 'OpenRouter request failed', detail });
        }

        const content = data?.choices?.[0]?.message?.content || '';
        let review;
        try {
            review = JSON.parse(content);
        } catch {
            const match = String(content).match(/\{[\s\S]*\}/);
            if (!match) throw new Error('OpenRouter response did not contain JSON');
            review = JSON.parse(match[0]);
        }

        return res.status(200).json({ model, review, usage: data.usage || null });
    } catch (err) {
        const message = err?.name === 'AbortError' ? 'OpenRouter request timed out' : (err.message || 'AI review failed');
        return res.status(502).json({ error: message });
    } finally {
        clearTimeout(timer);
    }
}
