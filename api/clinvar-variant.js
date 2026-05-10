// Serverless proxy for fetching ClinVar VCV data for a specific variation ID.
// Uses efetch (VCV XML) for reliable per-condition somatic data, supplemented
// by esummary for aggregate germline/somatic/oncogenicity classifications.
//
// GET /api/clinvar-variant?id=375941
//
// Optional env var: NCBI_API_KEY — increases rate limit from 3 to 10 req/s.

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

function apiKey() {
    return process.env.NCBI_API_KEY ? `&api_key=${encodeURIComponent(process.env.NCBI_API_KEY)}` : '';
}

// Extract text content of the first occurrence of a tag (non-nested).
function tagText(xml, tag) {
    const m = xml.match(new RegExp(`<${tag}(?:\\s[^>]*)?>([^<]*)</${tag}>`, 'i'));
    return m ? m[1].trim() : '';
}

// Extract the value of a named attribute from an opening tag string.
function attr(tagStr, name) {
    const m = tagStr.match(new RegExp(`\\b${name}="([^"]*)"`));
    return m ? m[1] : '';
}

// Capitalise first letter of a string.
function cap(s) { return s ? s[0].toUpperCase() + s.slice(1) : s; }

// Format a date string: keep only the YYYY-MM-DD portion if a time component is present.
function fmtDate(s) {
    if (!s) return '';
    return s.includes('T') ? s.split('T')[0] : s;
}

// Parse per-condition somatic entries from VCV XML RCVAccession blocks.
// Returns an array of { condition, tier, assertionType, clinSig, lastEvaluated, accession }.
function parseSomaticConditions(xml) {
    const results = [];

    // Split by <RCVAccession opening tag; index 0 is pre-content, skip it.
    const parts = xml.split(/<RCVAccession\b/);
    for (let i = 1; i < parts.length; i++) {
        const chunk = parts[i];

        // Only process RCV blocks that mention somatic impact.
        if (!chunk.includes('SomaticClinicalImpact')) continue;

        const accession = attr(chunk, 'Accession');

        // Condition name: first ClassifiedCondition text.
        const condMatch = chunk.match(/<ClassifiedCondition[^>]*>([^<]*)<\/ClassifiedCondition>/i);
        const condition = condMatch ? condMatch[1].trim() : '';

        // SomaticClinicalImpact block within this RCV chunk.
        const sciMatch = chunk.match(/<SomaticClinicalImpact\b([^>]*)>([\s\S]*?)<\/SomaticClinicalImpact>/i);
        if (!sciMatch) continue;

        const sciAttrs = sciMatch[1];
        const sciBlock = sciMatch[2];

        const lastEvaluated = fmtDate(attr(sciAttrs, 'DateLastEvaluated'));

        // Description element with assertion type and clinical significance.
        const descMatch = sciBlock.match(/<Description\b([^>]*)>([^<]*)<\/Description>/i);
        const tier = descMatch ? descMatch[2].trim() : '';
        const descAttrs = descMatch ? descMatch[1] : '';
        const assertionType = cap(attr(descAttrs, 'ClinicalImpactAssertionType'));
        const clinSig = cap(attr(descAttrs, 'ClinicalImpactClinicalSignificance').replace(/\//g, '/'));

        results.push({ condition, tier, assertionType, clinSig, lastEvaluated, accession });
    }

    return results;
}

// Parse aggregate classifications from the top-level <Classifications> block.
// Returns { germline, somatic, oncogenicity } each as { description, lastEvaluated } or null.
function parseAggregateClassifications(xml) {
    // The top-level <Classifications> appears inside <ClassifiedRecord>, after <RCVList>.
    // We find it by looking for the block that is NOT inside an RCVAccession.
    // Strategy: strip everything between <RCVList> and </RCVList> first.
    const stripped = xml.replace(/<RCVList\b[\s\S]*?<\/RCVList>/i, '');

    const classMatch = stripped.match(/<Classifications\b[^>]*>([\s\S]*?)<\/Classifications>/i);
    const classBlock = classMatch ? classMatch[1] : stripped;

    function extractClass(block, tag) {
        const m = block.match(new RegExp(`<${tag}\\b([^>]*)>([\\s\\S]*?)</${tag}>`, 'i'));
        if (!m) return null;
        const tagAttrs = m[1];
        const inner = m[2];
        // Description may carry DateLastEvaluated as attribute or as a child attribute.
        const descMatch = inner.match(/<Description\b([^>]*)>([^<]*)<\/Description>/i);
        const description = descMatch ? descMatch[2].trim() : tagText(inner, 'Description');
        const dateFromDesc = descMatch ? fmtDate(attr(descMatch[1], 'DateLastEvaluated')) : '';
        const dateFromTag = fmtDate(attr(tagAttrs, 'DateLastEvaluated'));
        const lastEvaluated = dateFromDesc || dateFromTag;
        if (!description) return null;
        return { description, lastEvaluated };
    }

    return {
        germline: extractClass(classBlock, 'GermlineClassification'),
        somatic: extractClass(classBlock, 'SomaticClinicalImpact'),
        oncogenicity: extractClass(classBlock, 'OncogenicityClassification'),
    };
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });

    const { id } = req.query;
    if (!id || !/^\d+$/.test(String(id))) {
        return res.status(400).json({ error: 'Missing or invalid id parameter (must be numeric ClinVar variation ID)' });
    }

    try {
        // VCV efetch gives the full classification tree including per-condition somatic data.
        const vcvUrl = `${EUTILS_BASE}/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=${encodeURIComponent(id)}${apiKey()}`;
        const vcvRes = await fetch(vcvUrl, {
            headers: { 'Accept': 'application/xml, text/xml' },
            signal: AbortSignal.timeout(15000),
        });
        if (!vcvRes.ok) throw new Error(`VCV efetch failed: ${vcvRes.status}`);
        const xml = await vcvRes.text();

        const aggregate = parseAggregateClassifications(xml);
        const somaticConditions = parseSomaticConditions(xml);

        return res.status(200).json({
            germline: aggregate.germline,
            somatic: aggregate.somatic,
            oncogenicity: aggregate.oncogenicity,
            somaticConditions,
        });
    } catch (err) {
        console.error('ClinVar variant proxy error:', err);
        return res.status(502).json({ error: 'ClinVar variant lookup failed', detail: err.message });
    }
}
