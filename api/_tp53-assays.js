// A small comparison checked against the source tables. See README.md for row
// identifiers and links. The H1299 result is specific to the reported construct.
const KOTLER_CONSTRUCTS = new Map([
    ['C.709A>C', { id: '302463', score: -1.779194794 }],
    ['C.709A>T', { id: '302465', score: -1.666323792 }]
]);

const M237L_PROTEIN_ASSAYS = [
    {
        study: 'Kato et al. (2003)',
        article: 'https://pubmed.ncbi.nlm.nih.gov/12826609/',
        dataLabel: 'NCI TP53 Database r21',
        dataUrl: 'https://tp53.cancer.gov/static/data/FunctionIshiokaDownload_r21.csv',
        scope: 'protein',
        variant: 'p.Met237Leu',
        system: 'Yeast reporter assay with eight p53-responsive promoters',
        endpoint: 'Promoter transactivation relative to wild-type p53',
        comparator: 'Each promoter’s wild-type p53 activity is 100%.',
        measurements: [
            { label: 'Median across eight promoters', value: 40.6, units: '% of wild type' },
            { label: 'Lowest promoter activity', value: 6.6, units: '% of wild type' },
            { label: 'Highest promoter activity', value: 78.2, units: '% of wild type' }
        ]
    },
    {
        study: 'Giacomelli et al. (2018)',
        article: 'https://pubmed.ncbi.nlm.nih.gov/30224644/',
        dataLabel: 'Supplementary Table 3',
        dataUrl: 'https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-018-0204-y/MediaObjects/41588_2018_204_MOESM6_ESM.xlsx',
        scope: 'protein',
        variant: 'p.Met237Leu',
        system: 'Isogenic human A549 cells with wild-type or null p53',
        endpoint: 'Variant enrichment or depletion after treatment relative to the early sample',
        comparator: 'Synonymous variants were included for comparison. Table 3 reports scores by amino-acid change.',
        measurements: [
            { label: 'p53 WT · nutlin-3', value: 0.7432024599950522, units: 'Z-score' },
            { label: 'p53 null · nutlin-3', value: -0.1555150626148093, units: 'Z-score' },
            { label: 'p53 null · etoposide', value: 0.2662258591870247, units: 'Z-score' }
        ]
    }
];

export function tp53AssayComparison(protein, exactTranscriptHgvs) {
    if (protein !== 'M237L') return null;
    const studies = [...M237L_PROTEIN_ASSAYS];
    const cdna = String(exactTranscriptHgvs || '').toUpperCase().match(/^NM_000546\.[56]:(C\.709A>[CT])$/)?.[1];
    const construct = KOTLER_CONSTRUCTS.get(cdna);
    if (construct) {
        studies.push({
            study: 'Kotler et al. (2018)',
            article: 'https://pubmed.ncbi.nlm.nih.gov/29979965/',
            dataLabel: 'Published Table S3 · GSE115072',
            dataUrl: 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115072/suppl/GSE115072_H1299_RFS_sequence_variants.xlsx',
            correctionUrl: 'https://pubmed.ncbi.nlm.nih.gov/30193102/',
            scope: 'coding_construct',
            variant: `${exactTranscriptHgvs} · wild-type backbone · construct ${construct.id}`,
            system: 'Human p53-null H1299 cells expressing a TP53 construct',
            endpoint: 'Change in relative cell abundance during the proliferation assay',
            comparator: 'Wild-type-like synonymous constructs were used to normalize the sub-libraries. Lower scores indicate greater depletion.',
            measurements: [{ label: 'H1299 relative fitness score', value: construct.score, units: 'log₂ scale' }]
        });
    }
    return { protein: 'p.Met237Leu', studies };
}
