// GET /api/clinvar-variant?id=375941 — ClinVar VCV XML via NCBI EFetch.
import { ncbiFetchResponse, ncbiApiKeyParam } from './_ncbi.js';
import { attribute, child, children, childText, parseXml, text } from './_eutils-xml.js';
import { positiveQueryInteger } from './_input.js';
import { rejectDisallowedOrigin } from './_origin.js';
import { setEdgeCache, setNoStore } from './_cache.js';

const EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
const CLASSIFICATIONS = {
    germline: 'GermlineClassification',
    somatic: 'SomaticClinicalImpact',
    oncogenicity: 'OncogenicityClassification'
};

const date = value => value.split('T')[0];
const capitalize = value => value ? value[0].toUpperCase() + value.slice(1) : '';

function aggregateClassification(node) {
    const description = child(node, 'Description');
    if (!text(description)) return null;
    return {
        description: text(description),
        lastEvaluated: date(attribute(description, 'DateLastEvaluated') || attribute(node, 'DateLastEvaluated'))
    };
}

function clinicalImpact(node) {
    return {
        assertionType: capitalize(attribute(node, 'ClinicalImpactAssertionType')),
        clinSig: capitalize(attribute(node, 'ClinicalImpactClinicalSignificance')),
        drugs: attribute(node, 'DrugForTherapeuticAssertion')
    };
}

function parseSomaticConditions(rcvRecords) {
    return rcvRecords.flatMap(rcv => {
        const classification = child(child(rcv, 'RCVClassifications'), CLASSIFICATIONS.somatic);
        if (!classification) return [];
        const descriptions = children(classification, 'Description');
        const condition = children(child(rcv, 'ClassifiedConditionList'), 'ClassifiedCondition').map(text).filter(Boolean).join('; ');
        return descriptions.map(description => {
            const { assertionType, clinSig } = clinicalImpact(description);
            return {
                condition, tier: text(description), assertionType, clinSig,
                lastEvaluated: date(attribute(description, 'DateLastEvaluated') || attribute(classification, 'DateLastEvaluated')),
                accession: attribute(rcv, 'Accession')
            };
        });
    });
}

function parseConditions(assertion, mappings) {
    // VCV often gives the submitted trait as an identifier only. ClinVar's
    // mapping is keyed to the assertion, so a co-submitted condition stays separate.
    const sets = children(assertion, 'TraitSet');
    const mapped = mappings.map(mapping => {
        const medgen = child(mapping, 'MedGen');
        const identifiers = attribute(medgen, 'CUI') ? [{ db: 'MedGen', id: attribute(medgen, 'CUI') }] : [];
        if (attribute(mapping, 'MappingType') === 'XRef' && attribute(mapping, 'MappingRef') !== 'MedGen') {
            identifiers.push({ db: attribute(mapping, 'MappingRef'), id: attribute(mapping, 'MappingValue') });
        }
        return {
            name: attribute(medgen, 'Name'),
            type: attribute(mapping, 'TraitType'),
            identifiers: identifiers.filter(ref => ref.db && ref.id)
        };
    }).filter(condition => condition.name);
    const submitted = sets.flatMap(set => children(set, 'Trait').map(trait => {
        const names = children(trait, 'Name').flatMap(n => children(n, 'ElementValue'));
        const preferred = names.find(n => attribute(n, 'Type') === 'Preferred');
        return {
            name: text(preferred || names[0]),
            type: attribute(trait, 'Type') || attribute(set, 'Type'),
            identifiers: children(trait, 'XRef').map(ref => ({ db: attribute(ref, 'DB'), id: attribute(ref, 'ID') }))
                .filter(ref => ref.db && ref.id)
        };
    }));
    const conditions = [...mapped];
    for (const condition of submitted) {
        const mappedCondition = mapped.find(item => item.name === condition.name || item.identifiers.some(ref =>
            condition.identifiers.some(submittedRef => ref.db === submittedRef.db && ref.id === submittedRef.id)));
        if (mappedCondition) {
            for (const ref of condition.identifiers) {
                if (!mappedCondition.identifiers.some(item => item.db === ref.db && item.id === ref.id)) {
                    mappedCondition.identifiers.push(ref);
                }
            }
        } else if (condition.name || condition.identifiers.length) {
            conditions.push(condition);
        }
    }
    return conditions.filter((condition, index) => conditions.findIndex(other =>
        other.name === condition.name && other.type === condition.type
        && JSON.stringify(other.identifiers) === JSON.stringify(condition.identifiers)) === index);
}

function comments(node) {
    return children(node, 'Comment').map(text).filter(Boolean);
}

function evidenceDescriptions(node) {
    const descriptions = children(node, 'Description').map(element => ({ type: 'Description', text: text(element) }));
    const attributes = children(node, 'Attribute')
        .filter(element => /description|evidence/i.test(attribute(element, 'Type')))
        .map(element => ({ type: attribute(element, 'Type'), text: text(element) }));
    return [...descriptions, ...attributes].filter(entry => entry.text);
}

function parseSubmission(assertion, mappings) {
    const accessionNode = child(assertion, 'ClinVarAccession');
    const classification = child(assertion, 'Classification');
    const contributes = attribute(assertion, 'ContributesToAggregateClassification');
    const citationMap = new Map();
    const addCitations = (node, location, scope) => {
        for (const citation of children(node, 'Citation')) {
            const citationScope = scope === 'assertion_method' || attribute(citation, 'Type') !== 'general' ? scope : 'general';
            for (const id of children(citation, 'ID')) {
                const pmid = text(id);
                if (attribute(id, 'Source') !== 'PubMed' || !/^[1-9]\d{0,8}$/.test(pmid)) continue;
                const key = `${citationScope}:${pmid}`;
                const entry = citationMap.get(key) ?? { pmid, scope: citationScope, locations: [] };
                if (!entry.locations.includes(location)) entry.locations.push(location);
                citationMap.set(key, entry);
            }
        }
    };
    addCitations(classification, 'classification', 'variant_evidence');
    addCitations(assertion, 'submission', 'general');
    const classificationComments = [childText(classification, 'ExplanationOfClassification'), ...comments(classification)].filter(Boolean);
    const classifications = Object.entries(CLASSIFICATIONS).flatMap(([type, name]) => {
        const value = child(classification, name);
        if (!value) return [];
        addCitations(value, 'classification', 'variant_evidence');
        classificationComments.push(...comments(value));
        const description = child(value, 'Description');
        return [{
            type,
            description: text(description || value),
            lastEvaluated: date(attribute(description, 'DateLastEvaluated') || attribute(value, 'DateLastEvaluated')
                || attribute(classification, 'DateLastEvaluated')),
            ...(type === 'somatic' ? clinicalImpact(description || value) : {})
        }];
    });
    const noClassification = childText(classification, 'NoClassification');
    if (noClassification) classifications.push({ type: 'unclassified', description: noClassification,
        lastEvaluated: date(attribute(classification, 'DateLastEvaluated')) });
    const assertionCriteria = [];
    for (const set of children(assertion, 'AttributeSet')) {
        const criteria = children(set, 'Attribute').filter(node => attribute(node, 'Type') === 'AssertionMethod');
        assertionCriteria.push(...criteria.map(text).filter(Boolean));
        addCitations(set, 'submission_attribute', criteria.length ? 'assertion_method' : 'general');
    }
    const observations = children(child(assertion, 'ObservedInList'), 'ObservedIn').map(observation => {
        const sample = child(observation, 'Sample');
        const methodTypes = [];
        addCitations(observation, 'observation', 'variant_evidence');
        const descriptions = evidenceDescriptions(observation);
        const evidenceComments = comments(observation);
        for (const method of children(observation, 'Method')) {
            methodTypes.push(...children(method, 'MethodType').map(text).filter(Boolean));
            addCitations(method, 'observation_method', 'general');
        }
        for (const evidence of children(observation, 'ObservedData')) {
            descriptions.push(...evidenceDescriptions(evidence));
            evidenceComments.push(...comments(evidence));
            addCitations(evidence, 'observed_data', 'variant_evidence');
        }
        const functional = child(observation, 'FunctionalData');
        for (const type of ['FunctionalEffect', 'FunctionalConsequence', 'Result']) {
            descriptions.push(...children(functional, type).map(node => ({
                type, text: text(node), source: attribute(node, 'Source')
            })).filter(item => item.text));
        }
        const functionalComment = childText(functional, 'FunctionalConsequenceComment');
        if (functionalComment) evidenceComments.push(functionalComment);
        addCitations(functional, 'functional_data', 'variant_evidence');
        return {
            origin: childText(sample, 'Origin'), methodTypes, descriptions, comments: evidenceComments,
            conditions: parseConditions(observation, [])
        };
    });
    return {
        accession: attribute(accessionNode, 'Accession'),
        version: attribute(accessionNode, 'Version'),
        submitter: attribute(accessionNode, 'SubmitterName'),
        recordStatus: childText(assertion, 'RecordStatus'),
        contributesToAggregate: contributes === 'true' || contributes === '1' ? true
            : contributes === 'false' || contributes === '0' ? false : null,
        assertion: childText(assertion, 'Assertion'),
        reviewStatus: childText(classification, 'ReviewStatus'),
        lastEvaluated: date(attribute(classification, 'DateLastEvaluated')),
        classifications,
        conditions: parseConditions(assertion, mappings),
        assertionCriteria,
        classificationComments,
        submitterComments: comments(assertion),
        studyName: childText(assertion, 'StudyName'),
        studyDescription: childText(assertion, 'StudyDescription'),
        observations,
        citations: [...citationMap.values()]
    };
}

export function parseClinvarVariant(xml, id, retrievedAt = new Date().toISOString()) {
    const root = parseXml(xml);
    const archives = root.name === 'VariationArchive' ? [root]
        : root.name === 'ClinVarResult-Set' && !child(root, 'ERROR') && !child(root, 'Error') ? children(root, 'VariationArchive') : [];
    const archive = archives.find(node => attribute(node, 'VariationID') === String(id));
    const accession = attribute(archive, 'Accession');
    const classified = child(archive, 'ClassifiedRecord');
    const included = child(archive, 'IncludedRecord');
    const classifications = child(classified || included, 'Classifications');
    const assertionList = child(classified, 'ClinicalAssertionList');
    const includedAllele = child(included, 'SimpleAllele') || child(included, 'Haplotype');
    if (!/^VCV\d+$/.test(accession) || !classifications
        || (classified ? !children(assertionList, 'ClinicalAssertion').length : !includedAllele)) {
        throw new Error('ClinVar did not return a VCV record for the requested variation');
    }
    const rcvRecords = children(child(classified, 'RCVList'), 'RCVAccession');
    const mappingByAssertion = new Map();
    for (const mapping of children(child(classified, 'TraitMappingList'), 'TraitMapping')) {
        const assertionId = attribute(mapping, 'ClinicalAssertionID');
        if (!mappingByAssertion.has(assertionId)) mappingByAssertion.set(assertionId, []);
        mappingByAssertion.get(assertionId).push(mapping);
    }
    const submissions = children(assertionList, 'ClinicalAssertion')
        .map(assertion => parseSubmission(assertion, mappingByAssertion.get(attribute(assertion, 'ID')) ?? []));
    const reportedCount = attribute(archive, 'NumberOfSubmissions');
    return {
        ...Object.fromEntries(Object.entries(CLASSIFICATIONS).map(([type, name]) => [type, aggregateClassification(child(classifications, name))])),
        somaticConditions: parseSomaticConditions(rcvRecords),
        record: {
            variationId: String(id), accession, version: attribute(archive, 'Version'),
            name: attribute(archive, 'VariationName'), recordStatus: childText(archive, 'RecordStatus'),
            lastUpdated: date(attribute(archive, 'DateLastUpdated')),
            reportedSubmissionCount: /^\d+$/.test(reportedCount) ? Number(reportedCount) : null,
            returnedSubmissionCount: submissions.length,
            source: 'NCBI ClinVar VCV (EFetch)', url: `https://www.ncbi.nlm.nih.gov/clinvar/variation/${id}/`, retrievedAt
        },
        submissions
    };
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });
    if (rejectDisallowedOrigin(req, res)) return;

    const id = positiveQueryInteger(req.query?.id);
    if (id === null) {
        return res.status(400).json({ error: 'Missing or invalid id parameter (must be numeric ClinVar variation ID)' });
    }
    try {
        const vcvUrl = `${EUTILS_BASE}/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=${id}${ncbiApiKeyParam()}`;
        const response = await ncbiFetchResponse(vcvUrl, { accept: 'application/xml, text/xml', timeout: 15000 });
        const data = parseClinvarVariant(await response.text(), id);
        setEdgeCache(res, 21600);
        return res.status(200).json(data);
    } catch (err) {
        setNoStore(res);
        console.error('ClinVar variant proxy error:', err.message);
        return res.status(502).json({ error: 'ClinVar variant lookup failed', detail: 'Could not read the requested VCV record from NCBI' });
    }
}
