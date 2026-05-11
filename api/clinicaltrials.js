// Serverless proxy for the ClinicalTrials.gov v2 REST API.
// GET /api/clinicaltrials?gene=BRAF&tumorType=melanoma
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const CT_API_BASE = 'https://clinicaltrials.gov/api/v2/studies';

const PHASES_PHASE2_PLUS = new Set(['PHASE2', 'PHASE3', 'PHASE4']);

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

function mapStudy(study) {
    const proto = study.protocolSection || {};
    const idMod = proto.identificationModule || {};
    const statusMod = proto.statusModule || {};
    const designMod = proto.designModule || {};
    const descMod = proto.descriptionModule || {};
    const armsMod = proto.armsInterventionsModule || {};
    const locMod = proto.contactsLocationsModule || {};

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

    const briefSummary = String(descMod.briefSummary || '').trim();

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
        briefSummary: briefSummary.slice(0, 400),
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

    const params = new URLSearchParams({
        'query.term': searchTerm,
        'filter.overallStatus': 'RECRUITING',
        'aggFilters': 'phase:phase2,phase3,phase4,studyType:int',
        'pageSize': '100',
        'format': 'json',
    });

    const url = `${CT_API_BASE}?${params}`;

    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 14000);

    try {
        const response = await fetch(url, {
            headers: { 'Accept': 'application/json', 'User-Agent': 'variant-search-tool/1.0' },
            signal: controller.signal
        });

        if (!response.ok) {
            const text = await response.text();
            console.warn('ClinicalTrials.gov API error:', response.status, text.slice(0, 200));
            return res.status(200).json({ error: `ClinicalTrials.gov returned ${response.status}`, studies: [], total: 0 });
        }

        const data = await response.json();
        const allStudies = data.studies || [];

        // Filter: interventional, treatment purpose, phase 2+, has US location
        const filtered = allStudies.filter(study => {
            const proto = study.protocolSection || {};
            const designMod = proto.designModule || {};
            const studyType = String(designMod.studyType || '').toUpperCase();
            const primaryPurpose = String((designMod.designInfo || {}).primaryPurpose || '').toUpperCase();
            const phases = designMod.phases || [];
            const locations = (proto.contactsLocationsModule || {}).locations || [];

            if (studyType !== 'INTERVENTIONAL') return false;
            if (primaryPurpose && primaryPurpose !== 'TREATMENT') return false;
            if (!isPhase2Plus(phases)) return false;
            if (!hasUsLocation(locations)) return false;
            return true;
        });

        const studies = filtered.map(mapStudy);

        return res.status(200).json({
            total: studies.length,
            studies,
            searchTerm
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
