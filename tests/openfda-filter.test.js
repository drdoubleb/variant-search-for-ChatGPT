/*
 * Tests for the openFDA gene-negation filter.
 *
 * The browser-side helpers live in script.js (see "openFDA gene-negation
 * filter" section). Because the project has no module system, the regexes
 * are re-declared here verbatim — KEEP IN SYNC with script.js when either
 * side changes.
 *
 * Run with: node tests/openfda-filter.test.js
 */

const OPENFDA_GENE_TOKEN = String.raw`[A-Z][A-Z0-9]{1,7}(?:[-/][A-Z0-9]{1,7})?`;
const OPENFDA_NEGATION_LIST_RE = new RegExp(
    String.raw`\b(?:[Nn]o|[Ww]ithout|[Ww]hose\s+tumors?\s+(?:have\s+no|do\s+not\s+have|lack)|[Ll]acking|[Aa]bsence\s+of)\s+` +
    `(?:${OPENFDA_GENE_TOKEN}(?:\\s*,\\s*|\\s+(?:or|and)\\s+|\\s+))*${OPENFDA_GENE_TOKEN}` +
    String.raw`\s+(?:genomic\s+(?:tumor\s+)?)?(?:[Aa]berration|[Aa]lteration|[Mm]utation|[Rr]earrangement|[Ff]usion|[Dd]river|[Aa]mplification)s?\b`,
    'g'
);
const OPENFDA_WILDTYPE_TRAIL_RE = new RegExp(
    String.raw`\b${OPENFDA_GENE_TOKEN}[- ]wild[- ]?type\b`,
    'g'
);
const OPENFDA_WILDTYPE_LEAD_RE = new RegExp(
    String.raw`\bwild[- ]?type\s+${OPENFDA_GENE_TOKEN}\b`,
    'g'
);

function findOpenFdaNegationSpans(text) {
    if (!text) return [];
    const spans = [];
    for (const re of [OPENFDA_NEGATION_LIST_RE, OPENFDA_WILDTYPE_TRAIL_RE, OPENFDA_WILDTYPE_LEAD_RE]) {
        re.lastIndex = 0;
        let m;
        while ((m = re.exec(text)) !== null) {
            spans.push({ start: m.index, end: m.index + m[0].length });
        }
    }
    return spans;
}

function openFdaGeneOnlyInNegativeContext(text, gene) {
    if (!text || !gene) return false;
    const positions = [];
    let i = 0;
    while ((i = text.indexOf(gene, i)) !== -1) {
        positions.push(i);
        i += gene.length;
    }
    if (positions.length === 0) return false;
    const spans = findOpenFdaNegationSpans(text);
    if (spans.length === 0) return false;
    return positions.every(pos => spans.some(s => pos >= s.start && pos < s.end));
}

// ── Test cases ────────────────────────────────────────────────────────────
// Each case lists the queried gene, an excerpt that mimics real FDA
// indications_and_usage text, and the expected outcome (true = excluded
// because every gene mention is in a negation; false = kept).

const cases = [
    // ── The boilerplate the user called out ──
    {
        name: 'pembrolizumab NSCLC: "no EGFR or ALK genomic tumor aberrations" — query EGFR',
        gene: 'EGFR',
        text: 'KEYTRUDA, in combination with pemetrexed and platinum chemotherapy, is indicated for the first-line treatment of patients with metastatic nonsquamous non-small cell lung cancer (NSCLC), with no EGFR or ALK genomic tumor aberrations.',
        expectExcluded: true,
    },
    {
        name: 'same indication, query ALK',
        gene: 'ALK',
        text: 'KEYTRUDA, in combination with pemetrexed and platinum chemotherapy, is indicated for the first-line treatment of patients with metastatic nonsquamous non-small cell lung cancer (NSCLC), with no EGFR or ALK genomic tumor aberrations.',
        expectExcluded: true,
    },
    {
        name: 'three-gene list "no EGFR, ALK or ROS1 aberrations" — query ROS1',
        gene: 'ROS1',
        text: 'Indicated for patients with metastatic NSCLC with no EGFR, ALK or ROS1 aberrations and PD-L1 expression ≥1%.',
        expectExcluded: true,
    },
    {
        name: 'three-gene list — query ALK',
        gene: 'ALK',
        text: 'Indicated for patients with metastatic NSCLC with no EGFR, ALK or ROS1 aberrations and PD-L1 expression ≥1%.',
        expectExcluded: true,
    },

    // ── Other negation phrasings ──
    {
        name: '"without EGFR or ALK alterations"',
        gene: 'EGFR',
        text: 'Treatment of metastatic NSCLC without EGFR or ALK alterations.',
        expectExcluded: true,
    },
    {
        name: '"whose tumors have no EGFR or ALK mutations"',
        gene: 'EGFR',
        text: 'For adults whose tumors have no EGFR or ALK mutations.',
        expectExcluded: true,
    },
    {
        name: '"BRAF wild-type" melanoma — query BRAF',
        gene: 'BRAF',
        text: 'Indicated for the treatment of unresectable or metastatic BRAF wild-type melanoma.',
        expectExcluded: true,
    },
    {
        name: '"wild-type EGFR" leading form',
        gene: 'EGFR',
        text: 'Indicated for wild-type EGFR colorectal cancer.',
        expectExcluded: true,
    },

    // ── True positives that should NOT be excluded ──
    {
        name: 'osimertinib-style positive indication',
        gene: 'EGFR',
        text: 'Indicated for the first-line treatment of patients with metastatic NSCLC whose tumors have EGFR exon 19 deletions or exon 21 L858R substitution mutations.',
        expectExcluded: false,
    },
    {
        name: 'mixed: positive EGFR mention + negation list mention — keep',
        gene: 'EGFR',
        text: 'Indicated for EGFR T790M-mutant NSCLC. Not for patients with no EGFR or ALK genomic tumor aberrations.',
        expectExcluded: false,
    },
    {
        name: 'BRAF V600E mutation indication',
        gene: 'BRAF',
        text: 'Indicated for unresectable or metastatic melanoma with BRAF V600E mutation as detected by an FDA-approved test.',
        expectExcluded: false,
    },
    {
        name: 'ALK-positive NSCLC indication',
        gene: 'ALK',
        text: 'Indicated for the treatment of patients with metastatic ALK-positive non-small cell lung cancer.',
        expectExcluded: false,
    },
    {
        name: '"no significant EGFR mutations" — soft negation, should NOT exclude',
        gene: 'EGFR',
        text: 'Studies showed no significant EGFR pathway activation in responders.',
        expectExcluded: false,
    },

    // ── Edge cases ──
    {
        name: 'EGFRvIII positive context with separate negation list',
        gene: 'EGFR',
        text: 'Indicated for EGFRvIII-amplified glioblastoma. Not indicated for patients with no EGFR or ALK genomic tumor aberrations.',
        // EGFR appears at EGFRvIII (positive) and inside the negation span.
        // Not all positions are negated → keep.
        expectExcluded: false,
    },
    {
        name: 'eGFR (lowercase, kidney function) — case-sensitive substring filter would already drop this upstream',
        gene: 'EGFR',
        text: 'Adjust dose in patients with eGFR < 30 mL/min/1.73m².',
        // This case is handled by the upstream `ind.includes(gene)` check
        // (case-sensitive); negation filter is irrelevant. We assert that
        // the gene literally does not appear at all.
        expectExcluded: false,
        expectGeneAbsent: true,
    },
    {
        name: '"no EGFR, ALK, ROS1 or BRAF mutations" — Oxford-comma variant, query BRAF',
        gene: 'BRAF',
        text: 'For patients with NSCLC and no EGFR, ALK, ROS1 or BRAF mutations.',
        expectExcluded: true,
    },
];

let passed = 0;
let failed = 0;
for (const tc of cases) {
    const geneAppears = tc.text.includes(tc.gene);
    if (tc.expectGeneAbsent) {
        if (geneAppears) {
            console.error(`FAIL: ${tc.name} — expected gene "${tc.gene}" to be absent from text`);
            failed++;
        } else {
            console.log(`PASS: ${tc.name}`);
            passed++;
        }
        continue;
    }
    if (!geneAppears) {
        console.error(`FAIL: ${tc.name} — gene "${tc.gene}" does not appear in test text (case-sensitive)`);
        failed++;
        continue;
    }
    const got = openFdaGeneOnlyInNegativeContext(tc.text, tc.gene);
    if (got === tc.expectExcluded) {
        console.log(`PASS: ${tc.name}`);
        passed++;
    } else {
        console.error(`FAIL: ${tc.name}`);
        console.error(`  expected excluded=${tc.expectExcluded}, got ${got}`);
        console.error(`  spans: ${JSON.stringify(findOpenFdaNegationSpans(tc.text).map(s => tc.text.slice(s.start, s.end)))}`);
        failed++;
    }
}

console.log(`\n${passed} passed, ${failed} failed`);
process.exit(failed === 0 ? 0 : 1);
