const PERCENT_WT = '% of wild type';
const NCI_DATA = { dataLabel: 'NCI TP53 data dictionary', dataUrl: 'https://tp53.cancer.gov/help' };
const CLINGEN_DATA = {
    cutoffLabel: 'ClinGen TP53 rules v2.4',
    cutoffUrl: 'https://erepo.clinicalgenome.org/cspec/ui/svi/doc/135637581?version=2.4.0'
};

// The NCI rules apply to the eight-promoter median and to two specific
// Giacomelli screens. The ClinGen rule uses the corrected H1299 RFS scale.
const SCORE_RANGES = {
    katoMedian: {
        label: 'NCI eight-promoter median (% of wild type)', inclusive: true,
        bands: [
            { upper: 20, label: 'Non-functional', tone: 'red' },
            { upper: 75, label: 'Partially functional', tone: 'orange' },
            { upper: 140, label: 'Functional', tone: 'green' },
            { label: 'Supertransactivation', tone: 'green' }
        ]
    },
    wtNutlin: {
        label: 'NCI WT nutlin-3 Z-score', inclusive: false,
        bands: [
            { upper: 0.61, label: 'Below DNE cutoff', tone: 'green' },
            { label: 'Meets DNE cutoff', tone: 'orange' }
        ]
    },
    nullEtoposide: {
        label: 'NCI null etoposide Z-score', inclusive: true,
        bands: [
            { upper: -0.21, label: 'LOF in this assay', tone: 'red' },
            { label: 'No LOF in this assay', tone: 'green' }
        ]
    },
    h1299Rfs: {
        label: 'ClinGen H1299 RFS', inclusive: false,
        bands: [
            { upper: -1, label: 'No LOF in this assay', tone: 'green' },
            { label: 'LOF in this assay', tone: 'red' }
        ]
    }
};

const DNE_LOF_LABELS = new Map([
    ['DNE_LOF', { label: 'DNE and LOF', tone: 'red' }],
    ['notDNE_LOF', { label: 'LOF without DNE', tone: 'red' }],
    ['notDNE_notLOF', { label: 'No DNE or LOF', tone: 'green' }]
]);
const DNE_LABELS = new Map([
    ['Yes', { label: 'DNE observed', tone: 'red' }],
    ['Moderate', { label: 'DNE on some promoters', tone: 'orange' }],
    ['No', { label: 'No DNE observed', tone: 'green' }]
]);

export function assayScoreInterpretation(value, assay) {
    if (!Number.isFinite(value) || !Object.hasOwn(SCORE_RANGES, assay)) return null;
    const { inclusive, bands } = SCORE_RANGES[assay];
    const { label, tone } = bands.find(({ upper }) => upper === undefined || (inclusive ? value <= upper : value < upper));
    return { label, tone };
}

function scoreRangeText(assay) {
    const { label, inclusive, bands } = SCORE_RANGES[assay];
    const format = value => String(value).replace(/^-/, '−');
    const descriptions = bands.map((band, index) => {
        const previous = bands[index - 1];
        const lower = previous ? `${inclusive ? '>' : '≥'} ${format(previous.upper)}` : '';
        const upper = band.upper === undefined ? '' : `${inclusive ? '≤' : '<'} ${format(band.upper)}`;
        return `${band.label} ${[lower, upper].filter(Boolean).join(' and ')}`;
    });
    return `${label} — ${descriptions.join(' · ')}`;
}

function classifiedResults(results, labels) {
    return results.map(result => {
        const interpretation = labels.get(result.value);
        return interpretation ? { ...result, interpretation } : result;
    });
}

const KATO = {
    ...NCI_DATA,
    study: 'Kato et al. (2003)',
    article: 'https://pubmed.ncbi.nlm.nih.gov/12826609/',
    system: 'Yeast reporter assay with eight p53-responsive promoters',
    endpoint: 'Promoter transactivation relative to wild-type p53',
    comparator: 'Each promoter’s wild-type p53 activity is 100%.'
};
const KAKUDO = {
    ...NCI_DATA,
    study: 'Kakudo et al. (2005)',
    article: 'https://pubmed.ncbi.nlm.nih.gov/15781620/',
    system: 'Human Saos-2 cells',
    endpoint: 'Promoter transactivation and sub-G1 apoptosis readout',
    comparator: 'The reporter results use p53-null vector normalization. Values are expressed relative to wild-type p53.'
};
const GIACOMELLI = {
    ...NCI_DATA,
    study: 'Giacomelli et al. (2018)',
    article: 'https://pubmed.ncbi.nlm.nih.gov/30224644/',
    system: 'Isogenic human A549 cells with wild-type or null p53',
    endpoint: 'Variant enrichment or depletion after treatment relative to the early sample',
    comparator: 'NCI combines the WT nutlin-3 and null etoposide screens. When the WT nutlin-3 score meets the DNE cutoff and etoposide falls in the no-LOF range, the combined class is Unclassified. DNE means dominant-negative effect. LOF means loss of function.'
};
const CURATED_DNE = {
    ...NCI_DATA,
    study: 'Dominant-negative transactivation',
    system: 'Multiple transactivation studies curated by the NCI TP53 Database',
    endpoint: 'Dominant-negative effect on wild-type p53 transactivation'
};

// These supplemental results were checked independently of MutationView. The
// general panel uses the measurements and classifications in that dataset.
// Kotler’s values apply to the specified coding construct and backbone only.
const M237L_SCORES = [
    { label: 'p53 WT · nutlin-3', value: 0.7432024599950522, assay: 'wtNutlin' },
    { label: 'p53 null · nutlin-3', value: -0.1555150626148093 },
    { label: 'p53 null · etoposide', value: 0.2662258591870247, assay: 'nullEtoposide' }
].map(({ label, value, assay }) => ({
    label, value, units: 'Z-score',
    ...(assay ? { interpretation: assayScoreInterpretation(value, assay) } : {})
}));
const GIACOMELLI_TABLE = 'https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-018-0204-y/MediaObjects/41588_2018_204_MOESM6_ESM.xlsx';
const KOTLER_CONSTRUCTS = new Map([
    ['C.709A>C', { id: '302463', score: -1.779194794 }],
    ['C.709A>T', { id: '302465', score: -1.666323792 }]
]);

function numericValue(text) {
    if (text == null || String(text).trim() === '') return null;
    const value = Number(text);
    return Number.isFinite(value) ? value : null;
}

function classificationValue(text) {
    const value = String(text || '').trim();
    if (!value || /^(?:N\/?A|ND|-+)$/i.test(value)) return null;
    return /^unclass\.?$/i.test(value) ? 'Unclassified' : value;
}

function studyResults(matches, readFields, readValue, units = '') {
    const targets = new Map();
    for (const match of matches) {
        const accession = match.transcript_hgvs || match.cdna_or_genomic;
        const source = [match.mut_id && `MutationView ${match.mut_id}`, accession].filter(Boolean).join(' · ');
        for (const [label, raw] of Object.entries(readFields(match))) {
            const value = readValue(raw);
            if (value === null) continue;
            if (!targets.has(label)) targets.set(label, new Map());
            const results = targets.get(label);
            if (!results.has(value)) results.set(value, new Set());
            if (source) results.get(value).add(source);
        }
    }
    return [...targets].flatMap(([label, results]) => [...results].map(([value, references]) => ({
        label, value, units,
        ...(results.size > 1 && references.size ? { source_records: [...references] } : {})
    })));
}

export function tp53AssayComparison({ protein, proteinKey, matches, exactTranscriptHgvs }) {
    const studies = [];
    const addProteinStudy = (study, measurements, classifications = []) => {
        if (measurements.length || classifications.length) {
            studies.push({
                ...study, scope: 'protein', variant: protein, measurements,
                ...(classifications.length ? { classifications } : {})
            });
        }
    };

    const yeast = studyResults(matches, match => match.ta_targets, numericValue, PERCENT_WT);
    let kato = KATO;
    if (yeast.length === 8 && new Set(yeast.map(m => m.label)).size === 8) {
        const sorted = yeast.map(m => m.value).sort((a, b) => a - b);
        const value = (sorted[3] + sorted[4]) / 2;
        const median = { label: 'Median across eight promoters', value, units: PERCENT_WT };
        if (/^[A-Z]\d+[A-Z]$/.test(proteinKey) && proteinKey[0] !== proteinKey.at(-1)) {
            median.interpretation = assayScoreInterpretation(value, 'katoMedian');
            kato = { ...KATO, scoreRanges: [scoreRangeText('katoMedian')] };
        }
        yeast.unshift(median);
    }
    addProteinStudy(kato, yeast);
    addProteinStudy(KAKUDO, studyResults(matches, match => match.saos2_targets, numericValue, PERCENT_WT));

    const growthClass = studyResults(matches, match => ({ 'NCI combined DNE/LOF class': match.lof_class }), classificationValue);
    const isM237L = proteinKey === 'M237L';
    const giacomelli = isM237L ? {
        ...GIACOMELLI,
        supplementUrl: GIACOMELLI_TABLE,
        scoreRanges: [scoreRangeText('wtNutlin'), scoreRangeText('nullEtoposide')]
    } : GIACOMELLI;
    addProteinStudy(giacomelli, isM237L ? M237L_SCORES : [], classifiedResults(growthClass, DNE_LOF_LABELS));
    const dneClass = studyResults(matches, match => ({ 'NCI dominant-negative class': match.dne_class }), classificationValue);
    addProteinStudy(CURATED_DNE, [], classifiedResults(dneClass, DNE_LABELS));

    const cdna = String(exactTranscriptHgvs || '').toUpperCase().match(/^NM_000546\.[56]:(C\.709A>[CT])$/)?.[1];
    const construct = isM237L && KOTLER_CONSTRUCTS.get(cdna);
    if (construct) {
        studies.push({
            study: 'Kotler et al. (2018)',
            article: 'https://pubmed.ncbi.nlm.nih.gov/29979965/',
            dataLabel: 'Published Table S3 · GSE115072',
            dataUrl: 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115072/suppl/GSE115072_H1299_RFS_sequence_variants.xlsx',
            ...CLINGEN_DATA,
            scoreRanges: [scoreRangeText('h1299Rfs')],
            correctionUrl: 'https://pubmed.ncbi.nlm.nih.gov/30193102/',
            scope: 'coding_construct',
            variant: `${exactTranscriptHgvs} · wild-type backbone · construct ${construct.id}`,
            system: 'Human p53-null H1299 cells expressing a TP53 construct',
            endpoint: 'Change in relative cell abundance during the proliferation assay',
            comparator: 'Wild-type-like synonymous constructs were used to normalize the sub-libraries. Lower scores indicate greater depletion.',
            measurements: [{
                label: 'H1299 relative fitness score', value: construct.score, units: 'log₂ scale',
                interpretation: assayScoreInterpretation(construct.score, 'h1299Rfs')
            }]
        });
    }
    return studies.length ? { protein, includes_mutation_view: true, studies } : null;
}
