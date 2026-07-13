// Serverless proxy for the ClinicalTrials.gov v2 REST API.
// GET /api/clinicaltrials?gene=BRAF&tumorType=melanoma
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const CT_API_BASE = 'https://clinicaltrials.gov/api/v2/studies';

// Pagination bounds. Server-side filtering (below) shrinks the candidate set enough
// that one page almost always suffices; the extra page is a backstop for broad genes.
const CT_PAGE_SIZE = 1000;
const CT_MAX_PAGES = 2;

const PHASES_PHASE2_PLUS = new Set(['PHASE2', 'PHASE3', 'PHASE4']);

function buildCtQuery(params) {
    // Percent-encode values (space -> %20). The v2 API's aggFilters/filter.advanced
    // expressions were verified to accept this encoding.
    return Object.entries(params)
        .map(([k, v]) => `${k}=${encodeURIComponent(v)}`)
        .join('&');
}

// Fetch up to CT_MAX_PAGES pages, following nextPageToken. Returns the accumulated
// studies, the server-reported total (countTotal), and whether more pages remained
// beyond the cap. On a non-200 FIRST page, returns { ok: false, status } so the caller
// can fall back. Later-page failures just stop pagination with what was collected.
async function fetchCtPages(baseParams, signal) {
    const studies = [];
    let total = null;
    let pageToken = null;
    let morePages = false;
    for (let page = 0; page < CT_MAX_PAGES; page++) {
        const params = { ...baseParams, pageSize: String(CT_PAGE_SIZE), countTotal: 'true' };
        if (pageToken) params.pageToken = pageToken;
        const resp = await fetch(`${CT_API_BASE}?${buildCtQuery(params)}`, {
            headers: { 'Accept': 'application/json', 'User-Agent': 'variant-search-tool/1.0' },
            signal
        });
        if (!resp.ok) {
            if (page === 0) return { ok: false, status: resp.status, studies, total };
            break;
        }
        const data = await resp.json();
        if (total === null && typeof data.totalCount === 'number') total = data.totalCount;
        const batch = data.studies || [];
        studies.push(...batch);
        pageToken = data.nextPageToken || null;
        if (!pageToken || batch.length === 0) break;
        if (page === CT_MAX_PAGES - 1) morePages = true;
    }
    return { ok: true, studies, total, morePages };
}

function isPhase2Plus(phases) {
    if (!Array.isArray(phases) || phases.length === 0) return false;
    return phases.some(p => PHASES_PHASE2_PLUS.has(String(p).toUpperCase()));
}

function hasUsLocation(locations) {
    if (!Array.isArray(locations)) return false;
    return locations.some(loc => String(loc.country || '').toLowerCase() === 'united states');
}

function extractUsLocations(locations) {
    if (!Array.isArray(locations)) return [];
    return locations
        .filter(loc => String(loc.country || '').toLowerCase() === 'united states')
        .map(loc => [loc.facility, loc.city, loc.state].filter(Boolean).join(', '))
        .filter(Boolean);
}

function extractInclusionCriteria(text) {
    if (!text) return '';
    const inclMatch = text.match(/inclusion criteria[:\s]*([\s\S]*?)(?:\s*exclusion criteria|$)/i);
    if (inclMatch) return inclMatch[1].trim().slice(0, 1000);
    return text.slice(0, 1000);
}

function mapStudy(study) {
    const proto = study.protocolSection || {};
    const idMod = proto.identificationModule || {};
    const statusMod = proto.statusModule || {};
    const designMod = proto.designModule || {};
    const descMod = proto.descriptionModule || {};
    const armsMod = proto.armsInterventionsModule || {};
    const locMod = proto.contactsLocationsModule || {};
    const condMod = proto.conditionsModule || {};
    const eligMod = proto.eligibilityModule || {};

    const nctId = idMod.nctId || '';
    const title = idMod.briefTitle || '';
    const phases = designMod.phases || [];
    const studyType = designMod.studyType || '';
    const primaryPurpose = (designMod.designInfo || {}).primaryPurpose || '';
    const overallStatus = statusMod.overallStatus || '';
    const locations = locMod.locations || [];

    const interventions = (armsMod.interventions || [])
        .map(i => ({ type: i.type || '', name: i.name || '' }))
        .filter(i => i.name);

    const usLocations = extractUsLocations(locations);
    const usLocationCount = usLocations.length;

    const briefSummary = String(descMod.briefSummary || '').trim().slice(0, 1500);
    const conditions = condMod.conditions || [];
    const inclusionCriteria = extractInclusionCriteria(String(eligMod.eligibilityCriteria || ''));

    return {
        nctId,
        title,
        phases,
        studyType,
        primaryPurpose,
        overallStatus,
        interventions,
        usLocationCount,
        usLocationSample: usLocations.slice(0, 3),
        briefSummary,
        conditions,
        inclusionCriteria,
        url: nctId ? `https://clinicaltrials.gov/study/${nctId}` : ''
    };
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });

    const { gene, tumorType } = req.query;

    if (!gene || !String(gene).trim()) {
        return res.status(400).json({ error: 'Missing required query parameter: gene' });
    }

    const safeGene = String(gene).replace(/[^A-Za-z0-9\-_.]/g, '').slice(0, 20);
    if (!safeGene) return res.status(400).json({ error: 'Invalid gene name' });

    const safeTumorType = tumorType
        ? String(tumorType).replace(/[^A-Za-z0-9\s\-,.'()]/g, '').trim().slice(0, 100)
        : '';

    const searchTerm = safeTumorType ? `${safeGene} ${safeTumorType}` : safeGene;

    // Push the cheap, high-selectivity filters into the query so the page budget is
    // spent on candidates that can actually pass: recruiting, interventional, phase
    // 2/3/4, and a US location. (Verified against the live v2 API.) The remaining
    // filters — treatment purpose and the gene word-boundary text check — run
    // client-side below on the much smaller returned set.
    const filteredParams = {
        'query.term': searchTerm,
        'filter.overallStatus': 'RECRUITING',
        'aggFilters': 'studyType:int,phase:2 3 4',
        'filter.advanced': 'AREA[LocationCountry]United States'
    };

    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 14000);

    try {
        let usedServerFilters = true;
        let result = await fetchCtPages(filteredParams, controller.signal);
        if (!result.ok) {
            // A filter combination the API rejects should degrade, not break the card:
            // fall back to the plain term query and rely on client-side filtering.
            console.warn('ClinicalTrials.gov filtered query failed:', result.status, '— falling back to term-only');
            usedServerFilters = false;
            result = await fetchCtPages({ 'query.term': searchTerm }, controller.signal);
            if (!result.ok) {
                return res.status(200).json({ error: `ClinicalTrials.gov returned ${result.status}`, studies: [], total: 0 });
            }
        }

        const allStudies = result.studies;
        const serverMatchedTotal = result.total;
        // True when more studies matched than we retrieved (only realistic in the
        // unfiltered fallback path). Surfaced so the UI/AI never implies completeness.
        const truncated = !!result.morePages;

        // Word-boundary check so "BRAF" doesn't match "CRAFT" or eGFR contexts
        const geneRegex = new RegExp(`\\b${safeGene}\\b`, 'i');

        // Filter: gene present in text, interventional, treatment purpose, phase 2+, recruiting, has US location
        const filtered = allStudies.filter(study => {
            const proto = study.protocolSection || {};
            const designMod = proto.designModule || {};
            const studyType = String(designMod.studyType || '').toUpperCase();
            const primaryPurpose = String((designMod.designInfo || {}).primaryPurpose || '').toUpperCase();
            const phases = designMod.phases || [];
            const overallStatus = String((proto.statusModule || {}).overallStatus || '').toUpperCase();
            const locations = (proto.contactsLocationsModule || {}).locations || [];

            if (overallStatus !== 'RECRUITING') return false;
            if (studyType !== 'INTERVENTIONAL') return false;
            if (primaryPurpose && primaryPurpose !== 'TREATMENT') return false;
            if (!isPhase2Plus(phases)) return false;
            if (!hasUsLocation(locations)) return false;

            // Verify gene appears in the study text (title, summary, eligibility, keywords)
            const idMod = proto.identificationModule || {};
            const descMod = proto.descriptionModule || {};
            const eligMod = proto.eligibilityModule || {};
            const condMod = proto.conditionsModule || {};
            const haystack = [
                idMod.briefTitle,
                idMod.officialTitle,
                descMod.briefSummary,
                descMod.detailedDescription,
                eligMod.eligibilityCriteria,
                ...(condMod.conditions || []),
                ...(condMod.keywords || [])
            ].filter(Boolean).join(' ');
            return geneRegex.test(haystack);
        });

        const studies = filtered.map(mapStudy);

        return res.status(200).json({
            total: studies.length,
            studies,
            searchTerm,
            server_filtered: usedServerFilters,
            total_matched_before_client_filter: serverMatchedTotal,
            truncated
        });
    } catch (e) {
        if (e.name === 'AbortError') {
            return res.status(200).json({ error: 'Request timed out', studies: [], total: 0 });
        }
        console.error('ClinicalTrials.gov proxy error:', e.message);
        return res.status(200).json({ error: e.message || String(e), studies: [], total: 0 });
    } finally {
        clearTimeout(timer);
    }
}
