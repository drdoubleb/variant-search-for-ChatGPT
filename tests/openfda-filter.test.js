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
const OPENFDA_NOUN = String.raw`(?:[Aa]berration|[Aa]lteration|[Mm]utation|[Rr]earrangement|[Ff]usion|[Dd]river|[Aa]mplification)s?\b`;
const OPENFDA_NEGATION_LIST_RE = new RegExp(
    String.raw`\b(?:[Nn]o|[Ww]ithout|[Ww]hose\s+tumors?\s+(?:have\s+no|do\s+not\s+have|lack)|[Ll]acking|[Aa]bsence\s+of|[Nn]egative\s+for)\s+` +
    String.raw`[^.;:]{0,80}?\b${OPENFDA_GENE_TOKEN}\b` +
    String.raw`[^.;:]{0,80}?${OPENFDA_NOUN}` +
    String.raw`(?:\s+(?:or|and)\s+[^.;:]{0,80}?\b${OPENFDA_GENE_TOKEN}\b[^.;:]{0,80}?${OPENFDA_NOUN})*`,
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
const OPENFDA_PRIOR_THERAPY_RE = new RegExp(
    String.raw`\b[Pp]atients\s+with\s+[^.;:]{0,150}?\b${OPENFDA_GENE_TOKEN}\b` +
    String.raw`[^.;:]{0,200}?${OPENFDA_NOUN}` +
    String.raw`[^.;:]{0,150}?should\s+have\s+(?:disease\s+)?progression\s+on\s+FDA[- ]approved\s+therapy`,
    'g'
);
// "anti-<GENE>" / "anti‑<GENE>" / "anti <GENE>" — almost always refers to a
// prior therapy class (e.g. "previously treated with ... an anti-EGFR
// therapy" in Fruzaqla, Stivarga, Lonsurf for chemo-refractory mCRC).
// EGFR/HER2/VEGF-targeted labels themselves describe their drug as an
// "EGFR antagonist" / "HER2-directed antibody" / etc., not "anti-X", so
// this pattern doesn't suppress true positives. Includes regular hyphen
// (U+002D) and non-breaking hyphen (U+2011) seen in FDA label text.
const OPENFDA_ANTI_GENE_RE = new RegExp(
    String.raw`\banti[\-‑ ]${OPENFDA_GENE_TOKEN}\b`,
    'g'
);

function findOpenFdaNegationSpans(text) {
    if (!text) return [];
    const spans = [];
    const regexes = [
        OPENFDA_NEGATION_LIST_RE,
        OPENFDA_WILDTYPE_TRAIL_RE,
        OPENFDA_WILDTYPE_LEAD_RE,
        OPENFDA_PRIOR_THERAPY_RE,
        OPENFDA_ANTI_GENE_RE,
    ];
    for (const re of regexes) {
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

// Word-boundary filter — see script.js "openFDA word-boundary filter".
function openFdaEscapeRegExp(s) {
    return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
function openFdaGeneAppearsAsWord(text, gene) {
    if (!text || !gene) return false;
    return new RegExp(String.raw`(?<![A-Za-z])${openFdaEscapeRegExp(gene)}(?![A-Za-z])`).test(text);
}

// Gene synonyms & false-positive exceptions — see script.js
// "openFDA gene synonyms & false-positive exceptions". KEEP IN SYNC.
const OPENFDA_GENE_SYNONYMS = {
    KIT: ['c-Kit', 'CD117'],
};
const OPENFDA_FALSE_POSITIVE_PHRASES = {
    KIT: ['bowel prep kit'],
};
function openFdaSynonymsFor(gene) {
    return OPENFDA_GENE_SYNONYMS[gene] || [];
}
function openFdaWordSpans(text, term) {
    if (!text || !term) return [];
    const re = new RegExp(String.raw`(?<![A-Za-z])${openFdaEscapeRegExp(term)}(?![A-Za-z])`, 'g');
    const spans = [];
    let m;
    while ((m = re.exec(text)) !== null) {
        spans.push({ start: m.index, end: m.index + m[0].length });
        if (re.lastIndex === m.index) re.lastIndex++;
    }
    return spans;
}
function openFdaFalsePositiveSpans(text, phrases) {
    if (!text || !phrases || !phrases.length) return [];
    const spans = [];
    for (const p of phrases) {
        const re = new RegExp(openFdaEscapeRegExp(p), 'gi');
        let m;
        while ((m = re.exec(text)) !== null) {
            spans.push({ start: m.index, end: m.index + m[0].length });
            if (re.lastIndex === m.index) re.lastIndex++;
        }
    }
    return spans;
}
function openFdaRecordExclusionReason(ind, gene) {
    if (!ind || !gene) return 'case';
    const synonyms = openFdaSynonymsFor(gene);
    const fpSpans = openFdaFalsePositiveSpans(ind, OPENFDA_FALSE_POSITIVE_PHRASES[gene] || []);
    const outsideFp = (s) => !fpSpans.some(f => s.start >= f.start && s.end <= f.end);

    if (!ind.includes(gene) && !synonyms.some(t => ind.includes(t))) return 'case';

    const geneSpans = openFdaWordSpans(ind, gene).filter(outsideFp);
    const synSpans = synonyms.flatMap(t => openFdaWordSpans(ind, t)).filter(outsideFp);

    if (geneSpans.length === 0 && synSpans.length === 0) {
        const hadRawWord = openFdaGeneAppearsAsWord(ind, gene)
            || synonyms.some(t => openFdaGeneAppearsAsWord(ind, t));
        return (hadRawWord && fpSpans.length) ? 'falsePositive' : 'boundary';
    }

    if (synSpans.length === 0 && openFdaGeneOnlyInNegativeContext(ind, gene)) return 'negation';
    return '';
}

// ── Test cases ────────────────────────────────────────────────────────────
// Each case lists the queried gene, an excerpt that mimics real FDA
// indications_and_usage text, and the expected outcome (true = excluded
// because every gene mention is in a negation; false = kept).

const cases = [
    // ── The boilerplate the user originally called out ──
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

    // ── Adjective intervening between "no" and gene (Round 2 fix) ──
    {
        name: 'IMJUDO: "no sensitizing EGFR (...) mutation or ALK (...) aberrations"',
        gene: 'EGFR',
        text: 'IMJUDO, in combination with durvalumab and platinum-based chemotherapy, is indicated for the treatment of adult patients with metastatic NSCLC with no sensitizing epidermal growth factor receptor (EGFR) mutations or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.',
        expectExcluded: true,
    },
    {
        name: 'IMJUDO same text — query ALK',
        gene: 'ALK',
        text: 'IMJUDO, in combination with durvalumab and platinum-based chemotherapy, is indicated for the treatment of adult patients with metastatic NSCLC with no sensitizing epidermal growth factor receptor (EGFR) mutations or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.',
        expectExcluded: true,
    },
    {
        name: 'OPDIVO: "no known EGFR mutations or ALK rearrangements"',
        gene: 'EGFR',
        text: 'OPDIVO, in combination with platinum-doublet chemotherapy, is indicated for the neoadjuvant treatment of adult patients with resectable NSCLC and no known EGFR mutations or ALK rearrangements, followed by single-agent OPDIVO as adjuvant treatment after surgery.',
        expectExcluded: true,
    },
    {
        name: 'IMFINZI: "no sensitizing EGFR mutations or ALK genomic tumor aberrations"',
        gene: 'EGFR',
        text: 'IMFINZI, in combination with tremelimumab-actl and platinum-based chemotherapy, is indicated for the treatment of adult patients with metastatic NSCLC with no sensitizing EGFR mutations or ALK genomic tumor aberrations.',
        expectExcluded: true,
    },

    // ── Parenthetical expansion of gene name (Round 2 fix) ──
    {
        name: 'Pemetrexed: parenthetical expansion for both genes',
        gene: 'EGFR',
        text: 'Pemetrexed Injection is indicated in combination with pembrolizumab and platinum chemotherapy, for the initial treatment of patients with metastatic non-squamous non-small cell lung cancer (NSCLC), with no epidermal growth factor receptor (EGFR) or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.',
        expectExcluded: true,
    },
    {
        name: 'Pemetrexed parenthetical — query ALK',
        gene: 'ALK',
        text: 'Pemetrexed Injection is indicated in combination with pembrolizumab and platinum chemotherapy, for the initial treatment of patients with metastatic non-squamous non-small cell lung cancer (NSCLC), with no epidermal growth factor receptor (EGFR) or anaplastic lymphoma kinase (ALK) genomic tumor aberrations.',
        expectExcluded: true,
    },

    // ── "Patients with X aberrations should have progression on FDA-approved therapy" (Round 2 fix) ──
    {
        name: 'KEYTRUDA prior-therapy boilerplate — query EGFR',
        gene: 'EGFR',
        text: 'Patients with EGFR or ALK genomic tumor aberrations should have disease progression on FDA-approved therapy for these aberrations prior to receiving KEYTRUDA.',
        expectExcluded: true,
    },
    {
        name: 'TECENTRIQ prior-therapy boilerplate — query ALK',
        gene: 'ALK',
        text: 'Patients with EGFR or ALK genomic tumor aberrations should have disease progression on FDA-approved therapy for NSCLC harboring these aberrations prior to receiving TECENTRIQ.',
        expectExcluded: true,
    },
    {
        name: 'OPDIVO QVANTIG: combined "no known" + prior-therapy boilerplate',
        gene: 'EGFR',
        text: 'OPDIVO QVANTIG, in combination with platinum-doublet chemotherapy, is indicated for the neoadjuvant treatment of adult patients with resectable NSCLC and no known epidermal growth factor receptor (EGFR) mutations or anaplastic lymphoma kinase (ALK) rearrangements. Patients with EGFR or ALK genomic tumor aberrations should have disease progression on FDA-approved therapy for these aberrations prior to receiving OPDIVO QVANTIG.',
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
        name: 'CYRAMZA: erlotinib indication for EGFR-mutant NSCLC + sequencing carve-out for docetaxel arm — must KEEP',
        gene: 'EGFR',
        text: 'CYRAMZA, in combination with erlotinib, is indicated for the first-line treatment of adults with metastatic non-small cell lung cancer (NSCLC) whose tumors have epidermal growth factor receptor (EGFR) exon 19 deletions or exon 21 (L858R) substitution mutations. CYRAMZA, in combination with docetaxel, is indicated for the treatment of adults with metastatic non-small cell lung cancer (NSCLC) with disease progression on or after platinum-based chemotherapy. Patients with epidermal growth factor receptor (EGFR) or anaplastic lymphoma kinase (ALK) genomic tumor aberrations should have disease progression on FDA-approved therapy for these aberrations prior to receiving CYRAMZA.',
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
        name: '"no significant EGFR pathway activation" — soft negation, should NOT exclude',
        // The widened trigger→gene window does let "no significant EGFR" reach
        // the gene token, but the trailing noun ("activation") isn't in our
        // aberration|alteration|mutation|… list, so no negation span fires and
        // the record is correctly kept.
        gene: 'EGFR',
        text: 'Studies showed no significant EGFR pathway activation in responders.',
        expectExcluded: false,
    },

    // ── Edge cases ──
    {
        name: 'EGFRvIII positive context with separate negation list',
        gene: 'EGFR',
        text: 'Indicated for EGFRvIII-amplified glioblastoma. Not indicated for patients with no EGFR or ALK genomic tumor aberrations.',
        expectExcluded: false,
    },
    {
        name: 'eGFR (lowercase, kidney function) — case-sensitive substring filter would already drop this upstream',
        gene: 'EGFR',
        text: 'Adjust dose in patients with eGFR < 30 mL/min/1.73m².',
        expectExcluded: false,
        expectGeneAbsent: true,
    },
    {
        name: '"no EGFR, ALK, ROS1 or BRAF mutations" — Oxford-comma variant, query BRAF',
        gene: 'BRAF',
        text: 'For patients with NSCLC and no EGFR, ALK, ROS1 or BRAF mutations.',
        expectExcluded: true,
    },
    {
        name: 'sentence boundary: positive EGFR mention in next sentence after negation',
        gene: 'EGFR',
        text: 'Not indicated for patients with no EGFR or ALK genomic tumor aberrations. Separately, EGFR T790M-mutant patients should receive osimertinib.',
        expectExcluded: false,
    },
    {
        name: '"non-small cell lung cancer" must NOT trip the "no" trigger',
        gene: 'EGFR',
        text: 'For non-small cell lung cancer patients with EGFR exon 19 deletions.',
        expectExcluded: false,
    },
    // ── "anti-<GENE>" prior-therapy mentions in mCRC labels ──
    {
        name: 'Fruzaqla mCRC: "previously treated with ... an anti-EGFR therapy" — query EGFR',
        gene: 'EGFR',
        text: 'FRUZAQLA is indicated for the treatment of adult patients with metastatic colorectal cancer (mCRC) who have been previously treated with fluoropyrimidine-, oxaliplatin- and irinotecan-based chemotherapy, an anti-VEGF therapy, and, if RAS wild-type and medically appropriate, an anti-EGFR therapy.',
        expectExcluded: true,
    },
    {
        name: 'Stivarga mCRC: "previously treated with ... anti-EGFR therapy" — query EGFR',
        gene: 'EGFR',
        text: 'STIVARGA is indicated for the treatment of patients with metastatic colorectal cancer (CRC) who have been previously treated with fluoropyrimidine-, oxaliplatin- and irinotecan-based chemotherapy, an anti-VEGF therapy, and, if RAS wild-type, an anti-EGFR therapy.',
        expectExcluded: true,
    },
    {
        name: 'Lonsurf mCRC: "previously treated with ... an anti-EGFR therapy" — query EGFR',
        gene: 'EGFR',
        text: 'LONSURF is indicated as a single agent or in combination with bevacizumab for the treatment of adult patients with metastatic colorectal cancer previously treated with fluoropyrimidine-, oxaliplatin- and irinotecan-based chemotherapy, an anti-VEGF biological therapy, and if RAS wild-type, an anti-EGFR therapy.',
        expectExcluded: true,
    },
    {
        name: 'non-breaking hyphen variant: "anti‑EGFR therapy" — query EGFR',
        gene: 'EGFR',
        text: 'Indicated for patients previously treated with an anti‑EGFR therapy.',
        expectExcluded: true,
    },
    {
        name: 'Erbitux-style true positive: "EGFR antagonist" must NOT match anti-<GENE>',
        gene: 'EGFR',
        text: 'ERBITUX is an epidermal growth factor receptor (EGFR) antagonist indicated for the treatment of patients with KRAS wild-type, EGFR-expressing metastatic colorectal cancer.',
        expectExcluded: false,
    },

    // ── MET word-boundary cases (the original complaint) ──
    // "MET" is the worst offender because three uppercase letters land inside
    // many brand names and homeopathic ingredients. The word-boundary filter
    // must drop those while keeping genuine MET-oncogene labels.
    {
        name: 'AVANDAMET brand: "MET" embedded in brand name — query MET, boundary-excluded',
        gene: 'MET',
        text: 'AVANDAMET is indicated as an adjunct to diet and exercise to improve glycemic control in adults with type 2 diabetes mellitus.',
        expectBoundaryExcluded: true,
    },
    {
        name: 'homeopathic ARGENTUM METALICUM — query MET, boundary-excluded',
        gene: 'MET',
        text: 'ARGENTUM METALICUM is indicated for the temporary relief of hoarseness and laryngitis.',
        expectBoundaryExcluded: true,
    },
    {
        name: 'METHOTREXATE: gene letters lead a larger word — query MET, boundary-excluded',
        gene: 'MET',
        text: 'METHOTREXATE is indicated in the treatment of patients with rheumatoid arthritis.',
        expectBoundaryExcluded: true,
    },
    {
        name: 'capmatinib (TABRECTA) true positive: "(MET) exon 14 skipping" — keep',
        gene: 'MET',
        text: 'TABRECTA is indicated for the treatment of adult patients with metastatic non-small cell lung cancer (NSCLC) whose tumors have a mutation that leads to mesenchymal epithelial transition (MET) exon 14 skipping as detected by an FDA-approved test.',
        expectExcluded: false,
    },
    {
        name: 'MET amplification positive indication — keep',
        gene: 'MET',
        text: 'Indicated for adult patients with advanced solid tumors that have MET amplification.',
        expectExcluded: false,
    },
    {
        name: 'hyphenated "MET-amplified" still counts as a whole word — keep',
        gene: 'MET',
        text: 'Indicated for patients with MET-amplified non-small cell lung cancer.',
        expectExcluded: false,
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
    const appearsAsWord = openFdaGeneAppearsAsWord(tc.text, tc.gene);
    if (tc.expectBoundaryExcluded) {
        // Gene present only as a substring of a larger word (AVANDAMET,
        // METALICUM, …) — the word-boundary filter must drop it.
        if (appearsAsWord) {
            console.error(`FAIL: ${tc.name}`);
            console.error(`  expected boundary-excluded but gene "${tc.gene}" matched as a whole word`);
            failed++;
        } else {
            console.log(`PASS: ${tc.name}`);
            passed++;
        }
        continue;
    }
    // Every kept/negation case must survive the word-boundary filter, i.e.
    // the gene appears at least once as a standalone token.
    if (!appearsAsWord) {
        console.error(`FAIL: ${tc.name} — word-boundary filter wrongly dropped a real-word "${tc.gene}" mention`);
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

// ── Synonym & false-positive classification (openFdaRecordExclusionReason) ──
// gene is always passed upper-cased by the caller, mirroring fetchOpenFdaDrugLabels.
const reasonCases = [
    // ── KIT synonyms: the Gleevec/imatinib label never says all-caps "KIT" ──
    {
        name: 'Gleevec: "Kit (CD117)-positive" — kept via CD117 synonym',
        gene: 'KIT',
        text: 'Gleevec is indicated for the treatment of adult patients with Kit (CD117)-positive unresectable and/or metastatic malignant gastrointestinal stromal tumors (GIST).',
        expect: '',
    },
    {
        name: 'imatinib: "c-Kit mutation" — kept via c-Kit synonym',
        gene: 'KIT',
        text: 'Indicated for adult patients with c-Kit mutation-positive gastrointestinal stromal tumors.',
        expect: '',
    },
    {
        name: 'genuine all-caps KIT positive indication — kept',
        gene: 'KIT',
        text: 'Indicated for patients with KIT exon 11 mutations.',
        expect: '',
    },
    // ── KIT false positive: "BOWEL PREP KIT" is a colonoscopy prep ──
    {
        name: 'BOWEL PREP KIT — all-caps KIT only inside false-positive phrase',
        gene: 'KIT',
        text: 'SUTAB is a BOWEL PREP KIT indicated for cleansing of the colon in preparation for colonoscopy in adults.',
        expect: 'falsePositive',
    },
    {
        name: 'false-positive phrase plus a real KIT mention elsewhere — kept',
        gene: 'KIT',
        text: 'Supplied as a BOWEL PREP KIT. Also indicated for KIT-mutant GIST.',
        expect: '',
    },
    // ── Bare "Kit"/"kit" must NOT match — the "Kit" alias was removed for
    //    false positives; only all-caps KIT, c-Kit and CD117 are synonyms. ──
    {
        name: 'lowercase "kit" (e.g. dosing kit) — not a match',
        gene: 'KIT',
        text: 'Each carton contains a dosing kit with a syringe and adapter.',
        expect: 'case',
    },
    {
        name: 'capitalized standalone "Kit" (no CD117/c-Kit) — not a match',
        gene: 'KIT',
        text: 'Supplied as a Kit containing two vials and a diluent.',
        expect: 'case',
    },
    // ── Synonyms never trip negation boilerplate (all-caps token only) ──
    {
        name: 'c-Kit positive mention survives even with a negated all-caps KIT elsewhere',
        gene: 'KIT',
        text: 'Indicated for c-Kit mutation-positive GIST. Not for tumors with no KIT mutations.',
        expect: '',
    },
    // ── Genes without synonyms keep original behaviour exactly ──
    {
        name: 'no-synonym gene: EGFR negation boilerplate still excluded',
        gene: 'EGFR',
        text: 'Indicated for metastatic NSCLC with no EGFR or ALK genomic tumor aberrations.',
        expect: 'negation',
    },
    {
        name: 'no-synonym gene: MET substring of AVANDAMET still boundary-excluded',
        gene: 'MET',
        text: 'AVANDAMET is indicated as an adjunct to diet and exercise in type 2 diabetes mellitus.',
        expect: 'boundary',
    },
    {
        name: 'no-synonym gene: MET positive indication still kept',
        gene: 'MET',
        text: 'Indicated for adult patients with advanced solid tumors that have MET amplification.',
        expect: '',
    },
];

for (const tc of reasonCases) {
    const got = openFdaRecordExclusionReason(tc.text, tc.gene);
    if (got === tc.expect) {
        console.log(`PASS: ${tc.name}`);
        passed++;
    } else {
        console.error(`FAIL: ${tc.name}`);
        console.error(`  expected reason="${tc.expect}", got "${got}"`);
        failed++;
    }
}

console.log(`\n${passed} passed, ${failed} failed`);
process.exit(failed === 0 ? 0 : 1);
