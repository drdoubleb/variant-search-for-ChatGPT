// Shared schema and output normalization for the site's AI review.
const GERMLINE_CLASSIFICATIONS = ['Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign', 'Insufficient evidence'];
const SOMATIC_CLASSIFICATIONS = ['Oncogenic', 'Likely oncogenic', 'Uncertain significance', 'Likely benign', 'Benign', 'Insufficient evidence'];
const TIER_LEVELS = new Map([
    ['Tier I', ['A', 'B']], ['Tier II', ['C', 'D']], ['Tier III', [null]], ['Tier IV', [null]]
]);

export const REVIEW_SCHEMA = {
    germline_pathogenicity: {
        classification: GERMLINE_CLASSIFICATIONS.join(' | '),
        condition: 'disease or predisposition syndrome evaluated, or null',
        rationale: 'brief explanation with attributed evidence'
    },
    somatic_oncogenicity: {
        classification: SOMATIC_CLASSIFICATIONS.join(' | '),
        rationale: 'brief explanation with attributed evidence'
    },
    clinical_significance: {
        status: 'assessed | insufficient_disease_context | insufficient_evidence',
        tumor_type: 'repeat the exact submitted tumor type, or null if absent',
        tier: `${[...TIER_LEVELS.keys()].join(' | ')} | null`,
        evidence_level: 'A | B | C | D | null',
        rationale: 'brief explanation with the indication, alteration and evidence for the tier when assigned'
    },
    fda_approved_therapies: [{ drug: 'string', indication: 'string', biomarker_context: 'string', evidence: 'string' }],
    resistance_or_lack_of_benefit: [{ therapy: 'string', tumor_type: 'string', biomarker_context: 'string', evidence: 'string', url: 'source URL or null' }],
    clinical_trials: [{ nct_id: 'string', title: 'string', phase: 'string', intervention: 'string', relevance: 'string', url: 'string' }],
    summary: 'brief variant interpretation and clinical significance',
    limitations: ['brief caveats and verification needs']
};

const cleanText = value => typeof value === 'string' ? value.trim() : '';
const canonicalText = value => cleanText(value).replace(/\s+/g, ' ').toLowerCase();

function normalizeAssessment(assessment, labels) {
    const submitted = canonicalText(assessment?.classification);
    const classification = labels.find(label => label.toLowerCase() === submitted);
    return { classification: classification || 'Insufficient evidence', rationale: classification ? cleanText(assessment.rationale) : '' };
}

function normalizeClinicalSignificance(assessment, tumorType) {
    const result = { tumor_type: tumorType || null, tier: null, evidence_level: null };
    if (!tumorType) {
        return { ...result, status: 'insufficient_disease_context', rationale: 'Enter a tumor type to assess tumor-specific clinical significance.' };
    }
    const sameTumor = canonicalText(assessment?.tumor_type) === canonicalText(tumorType);
    if (!sameTumor) {
        return { ...result, status: 'insufficient_evidence', rationale: 'The review did not provide an assessment for the submitted tumor type.' };
    }
    const rationale = cleanText(assessment.rationale);
    const { tier, status } = assessment;
    const allowedLevels = TIER_LEVELS.get(tier);
    const evidenceLevel = cleanText(assessment.evidence_level).toUpperCase() || null;
    if (status === 'assessed' && allowedLevels?.includes(evidenceLevel) && rationale) {
        return { ...result, status, tier, evidence_level: evidenceLevel, rationale };
    }
    if (status === 'insufficient_disease_context' || status === 'insufficient_evidence') {
        return { ...result, status, rationale };
    }
    return { ...result, status: 'insufficient_evidence', rationale: 'The review did not provide a supported tier and evidence level for this tumor type.' };
}

export function normalizeReview(review, context) {
    if (!review || typeof review !== 'object' || Array.isArray(review)) {
        throw new Error('The model did not return a structured variant review. Try running the review again, or pick a different model.');
    }
    const germline = normalizeAssessment(review.germline_pathogenicity, GERMLINE_CLASSIFICATIONS);
    germline.condition = germline.classification !== 'Insufficient evidence' ? cleanText(review.germline_pathogenicity?.condition) || null : null;
    const somatic = normalizeAssessment(review.somatic_oncogenicity, SOMATIC_CLASSIFICATIONS);
    const clinical = normalizeClinicalSignificance(review.clinical_significance, cleanText(context.tumor_type));
    const legacyTier = clinical.tier ? `${clinical.tier}${clinical.evidence_level || ''}`
        : clinical.status === 'insufficient_disease_context' ? 'Insufficient disease context' : 'Insufficient evidence';
    return {
        ...review,
        germline_pathogenicity: germline,
        somatic_oncogenicity: somatic,
        clinical_significance: clinical,
        // Existing clients display these flat fields. Derive them from the same
        // validated conclusions so clients cannot revive a discarded tier.
        pathogenicity: germline.classification,
        amp_tier: legacyTier,
        amp_tier_rationale: clinical.rationale
    };
}
