// Test comment added to verify GitHub PR integration.
// Helper mapping from NCBI reference sequence accessions to chromosome names.
// According to NCBI, accessions NC_000001.11 ... NC_000024.10 correspond to chromosomes 1-22, X (23) and Y (24).
// NC_012920.1 corresponds to mitochondrial DNA (MT).
const accessionToChr = (accession) => {
    // Strip version suffix (e.g. ".11")
    const match = accession.match(/^NC_(\d{6})/);
    if (!match) return null;
    const numStr = match[1];
    const num = parseInt(numStr, 10);
    if (num >= 1 && num <= 22) return String(num);
    if (num === 23) return 'X';
    if (num === 24) return 'Y';
    // mitochondria (NC_012920)
    if (num === 12920) return 'MT';
    return null;
};


const DEFAULT_BACKEND_API_BASE_URL = 'https://variant-search-for-chat-gpt.vercel.app';

function trimTrailingSlash(value) {
    return String(value || '').trim().replace(/\/+$/, '');
}

function getBackendApiBaseUrl() {
    return trimTrailingSlash(window.BACKEND_API_BASE_URL || DEFAULT_BACKEND_API_BASE_URL);
}

function getConfiguredApiEndpoint(globalName, apiPath) {
    const configured = String(window[globalName] || '').trim();
    if (configured) return configured;
    const base = getBackendApiBaseUrl();
    return base ? `${base}${apiPath}` : apiPath;
}

function appendQueryParams(endpoint, params) {
    const separator = endpoint.includes('?') ? '&' : '?';
    return `${endpoint}${separator}${params}`;
}

// When the user provides a gene symbol in a five‑token genomic variant input
// (e.g. "7 140453122 BRAF TCCATCGAGATTTCA TCT"), we temporarily store that
// gene here during normalisation so that it can be used later when
// constructing a minimal annotation. Without this hint the code attempts to
// infer a gene symbol from the normalised genomic HGVS string, which may
// incorrectly yield "chr7" or similar. See normalizeVariantInput() below.
let geneHintGlobal = null;
let alterationTypeGlobal = null;
let isGeneOnlyMode = false;

// Keep third-party API latency from blocking the UI for long periods.
const API_TIMEOUT_MS = {
    myvariant: 7000,
    recoder: 6000,
    vep: 7000,
    lookup: 4000,
    liftover: 4000,
    cosmic: 2500,
    cosmicMeta: 1500,
    tp53: 15000,
    clinvar: 15000,
    civic: 8000,
    pubmed: 25000,
    fda: 8000,
    gnomadV4: 10000,
    clinicalTrials: 20000,
    spliceai: 25000,
    aiReview: 60000,
    openFda: 12000
};

async function fetchWithTimeout(url, options = {}, timeoutMs = 6000) {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), timeoutMs);
    try {
        return await fetch(url, { ...options, signal: controller.signal });
    } catch (err) {
        if (err && err.name === 'AbortError') {
            throw new Error(`Request timed out after ${timeoutMs}ms`);
        }
        throw err;
    } finally {
        clearTimeout(timeoutId);
    }
}

// Chromosome-like placeholders (e.g. "CHR12") can appear in fallback annotations
// and are not useful for user-facing search links.
function isChromosomeLikeGeneSymbol(symbol) {
    if (!symbol) return false;
    const s = String(symbol).trim().toUpperCase();
    return /^CHR(?:[0-9]+|X|Y|M|MT)$/.test(s);
}

// Determine whether the variant string is already in genomic (g.) format.
const isGenomicVariant = (variant) => {
    // Accept patterns like 'chr7:g.123456A>T' or 'NC_000007.13:g.123456A>T'.
    return /:g\./i.test(variant);
};


// Build a SpliceAI lookup tuple from either raw user input (preferred) or a normalised g. variant.
// Returns an object { chrom, pos, ref, alt } when enough information is available, else null.
function buildSpliceAiLookupTuple(rawInput, gVariant) {
    const parseTokenInput = (raw) => {
        if (!raw) return null;
        const toks = String(raw).trim().split(/\s+/).filter(Boolean);
        if (toks.length !== 4 && toks.length !== 5) return null;
        const tokens = toks.slice();
        if (tokens.length === 5) {
            // For 5-token genomic input, always treat token #3 as an optional gene/label
            // and remove it. Example: "17 7573954 TP53 T A".
            tokens.splice(2, 1);
        }
        const [chrTok, posTok, refTok, altTok] = tokens;
        const chrom = String(chrTok).replace(/^chr/i, '').toUpperCase();
        const pos = String(posTok).replace(/,/g, '');
        const ref = String(refTok).toUpperCase();
        const alt = String(altTok).toUpperCase();
        if (!/^[0-9XYMT]+$/.test(chrom)) return null;
        if (!/^\d+$/.test(pos)) return null;
        if (!/^[A-Za-z-]+$/.test(ref) || !/^[A-Za-z-]+$/.test(alt)) return null;
        return { chrom: `chr${chrom}`, pos, ref, alt };
    };
    const parseSimpleGenomic = (gv) => {
        if (!gv) return null;
        const m = String(gv).match(/^chr([0-9XYMT]+):g\.(\d+)([A-Za-z-]+)>([A-Za-z-]+)$/i);
        if (!m) return null;
        return { chrom: `chr${m[1].toUpperCase()}`, pos: m[2], ref: m[3].toUpperCase(), alt: m[4].toUpperCase() };
    };
    // Parses delins notation (MNVs/complex indels) to extract position; ref is unavailable.
    const parseDelins = (gv) => {
        if (!gv) return null;
        const m = String(gv).match(/^chr([0-9XYMT]+):g\.(\d+)(?:_(\d+))?delins([A-Za-z]+)$/i);
        if (!m) return null;
        return { chrom: `chr${m[1].toUpperCase()}`, pos: m[2], ref: null, alt: m[4].toUpperCase() };
    };
    return parseTokenInput(rawInput) || parseSimpleGenomic(gVariant) || parseDelins(gVariant) || null;
}

// Build a gnomAD variant-browser URL from either raw token input (preferred) or
// a simple genomic substitution HGVS string.
// vcfData: optional annotation.vcf object {ref, alt, position} from myvariant.info —
// used to recover REF for MNVs where the delins gVariant notation drops the reference.
function buildGnomadVariantUrl(rawInput, gVariant, annotationId, vcfData) {
    // Raw token input (e.g. "10 89692913 PTEN G GG") is preferred because it
    // preserves REF/ALT for indels even when gVariant is normalised to ins/delins.
    const tuple = buildSpliceAiLookupTuple(rawInput, gVariant);
    if (tuple) {
        const chrom = tuple.chrom.replace(/^chr/i, '');
        const ref = tuple.ref || (vcfData && String(vcfData.ref || '').toUpperCase()) || null;
        const alt = tuple.alt || (vcfData && String(vcfData.alt || '').toUpperCase()) || null;
        if (ref && alt) {
            return `https://gnomad.broadinstitute.org/variant/${chrom}-${tuple.pos}-${ref}-${alt}?dataset=gnomad_r2_1`;
        }
        // No REF available even from vcfData — fall back to a region view so the button still appears.
        const posNum = Number(tuple.pos);
        const altLen = (alt || '').length || 1;
        return `https://gnomad.broadinstitute.org/region/${chrom}-${posNum}-${posNum + altLen - 1}?dataset=gnomad_r2_1`;
    }
    // Fallback: use annotation/gVariant only for simple substitutions.
    const source = annotationId || gVariant || '';
    const m = String(source).match(/^chr([0-9XYMT]+):g\.(\d+)([A-Za-z-]+)>([A-Za-z-]+)$/i);
    if (!m) return '';
    return `https://gnomad.broadinstitute.org/variant/${m[1]}-${m[2]}-${m[3].toUpperCase()}-${m[4].toUpperCase()}?dataset=gnomad_r2_1`;
}

// Build a UCSC Genome Browser link (hg19/GRCh37) centered on a 20-nt window
// around the variant location (10 upstream, 9 downstream).
function buildUcscHg19Url(rawInput, gVariant, annotation) {
    const toUcscChrom = (chrom) => {
        if (!chrom) return '';
        const bare = String(chrom).trim().replace(/^chr/i, '').toUpperCase();
        if (!bare) return '';
        // UCSC uses chrM (not chrMT) for mitochondrial chromosome labels.
        const ucscBare = bare === 'MT' ? 'M' : bare;
        return `chr${ucscBare}`;
    };

    const toTwentyNtWindow = (start, end = start) => {
        const s = Number(start);
        const e = Number(end);
        if (!Number.isFinite(s) || !Number.isFinite(e)) return null;
        const center = Math.round((Math.min(s, e) + Math.max(s, e)) / 2);
        const winStart = Math.max(1, center - 10);
        const winEnd = winStart + 19;
        return { start: winStart, end: winEnd };
    };

    const buildUcscUrl = (region, highlightRegion) => {
        const qs = new URLSearchParams({ db: 'hg19', position: region });
        if (highlightRegion) {
            // UCSC highlight syntax: <db>.<chr>:<start>-<end>
            qs.set('highlight', `hg19.${highlightRegion}`);
        }
        return `https://genome.ucsc.edu/cgi-bin/hgTracks?${qs.toString()}`;
    };

    const tuple = buildSpliceAiLookupTuple(rawInput, gVariant);
    if (tuple && tuple.chrom && tuple.pos) {
        const ucscChrom = toUcscChrom(tuple.chrom);
        const pos = Number(tuple.pos);
        const win = toTwentyNtWindow(pos);
        if (ucscChrom && Number.isFinite(pos) && win) {
            const region = `${ucscChrom}:${win.start}-${win.end}`;
            const highlightRegion = `${ucscChrom}:${pos}-${pos}`;
            return buildUcscUrl(region, highlightRegion);
        }
    }

    const hg19 = annotation?.hg19 || annotation?.dbsnp?.hg19;
    const chrom = annotation?.chrom || annotation?.cadd?.chrom || annotation?.dbsnp?.chrom;
    if (hg19?.start !== undefined && chrom) {
        const ucscChrom = toUcscChrom(chrom);
        const hStart = Number(hg19.start);
        const hEnd = Number(hg19.end ?? hg19.start);
        const win = toTwentyNtWindow(hStart, hEnd);
        if (ucscChrom && Number.isFinite(hStart) && Number.isFinite(hEnd) && win) {
            const region = `${ucscChrom}:${win.start}-${win.end}`;
            const hMin = Math.min(hStart, hEnd);
            const hMax = Math.max(hStart, hEnd);
            const highlightRegion = `${ucscChrom}:${hMin}-${hMax}`;
            return buildUcscUrl(region, highlightRegion);
        }
    }
    return '';
}

// Convert a triple‑letter amino acid change (e.g. VAL600GLU) to a single‑letter code (V600E).
// Accepts uppercase three‑letter codes and returns uppercase single‑letter code if mapping exists.
function tripleToSingle(prot) {
    if (!prot) return null;
    const m = prot.match(/^([A-Z]{3})(\d+)([A-Z]{3})$/);
    if (!m) return null;
    const aaSingle = {
        ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
        HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
        THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V'
    };
    const ref = aaSingle[m[1]];
    const pos = m[2];
    const alt = aaSingle[m[3]];
    if (ref && alt) return `${ref}${pos}${alt}`;
    return null;
}

function aaThreeToSingle(aa) {
    if (!aa) return null;
    const aaSingle = {
        ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
        HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
        THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
        TER: '*', STOP: '*'
    };
    return aaSingle[String(aa).toUpperCase()] || null;
}

// Convert a 3-letter protein HGVS body (without "p.") to single-letter HGVS body.
// Handles substitutions, nonsense and common frameshift forms.
function convertProteinBodyToSingle(proteinBody) {
    if (!proteinBody) return null;
    const body = String(proteinBody).trim();

    // Frameshift with explicit downstream stop: Val133GlyfsTer47 -> V133Gfs*47
    let m = body.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3})fs(?:Ter|Stop|\*)(\d+)$/i);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]);
        if (ref && alt) return `${ref}${m[2]}${alt}fs*${m[4]}`;
    }

    // Frameshift without explicit alt amino acid: Arg97fsTer12 -> R97fs*12
    m = body.match(/^([A-Za-z]{3})(\d+)fs(?:Ter|Stop|\*)(\d+)$/i);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        if (ref) return `${ref}${m[2]}fs*${m[3]}`;
    }

    // Simple substitution / nonsense: Arg714Ter -> R714*
    m = body.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|Stop|\*)$/i);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]) || (m[3] === '*' ? '*' : null);
        if (ref && alt) return `${ref}${m[2]}${alt}`;
    }

    // Fallback to strict triple substitution converter.
    const compact = body.toUpperCase();
    return tripleToSingle(compact);
}

// Format a protein HGVS string to include the single-letter amino acid code in
// parentheses when a three-letter protein change is present.
// Example: "p.Val600Glu" -> "p.Val600Glu (p.V600E)".
function formatProteinDisplayWithSingleLetter(proteinHgvs) {
    if (!proteinHgvs) return '';
    const proteinText = String(proteinHgvs).trim();
    const proteinIndex = proteinText.toLowerCase().lastIndexOf('p.');
    if (proteinIndex === -1) return proteinText;
    const proteinBody = proteinText.slice(proteinIndex + 2).trim();
    const single = convertProteinBodyToSingle(proteinBody);
    if (!single) return proteinText;
    return `${proteinText} (p.${single})`;
}

function normaliseProteinForCivicMatch(proteinChange) {
    const proteinText = String(proteinChange || '').trim();
    if (!proteinText) return '';
    const proteinIndex = proteinText.toLowerCase().lastIndexOf('p.');
    const proteinBody = proteinIndex === -1 ? proteinText : proteinText.slice(proteinIndex + 2);
    const single = convertProteinBodyToSingle(proteinBody);
    return String(single || proteinBody)
        .replace(/^p\./i, '')
        .toUpperCase()
        .replace(/[^A-Z0-9*_]/g, '');
}

function findBestCivicEntryForProtein(entries, proteinChange) {
    if (!Array.isArray(entries) || entries.length === 0) return null;
    const normalisedProtein = normaliseProteinForCivicMatch(proteinChange);
    if (!normalisedProtein) return null;
    return entries.find((entry) => {
        const entryProtein = normaliseProteinForCivicMatch(entry?.protein_change);
        return entryProtein && (entryProtein === normalisedProtein
            || entryProtein.includes(normalisedProtein)
            || normalisedProtein.includes(entryProtein));
    }) || null;
}

function isTp53Gene(geneNames) {
    if (!geneNames) return false;
    return String(geneNames)
        .split(',')
        .map(g => g.trim().toUpperCase())
        .some(g => g === 'TP53');
}

async function fetchTp53MutationDatabase(payload) {
    const endpoint = getConfiguredApiEndpoint('TP53_API_ENDPOINT', '/api/tp53');
    if (!endpoint) return null;
    const response = await fetchWithTimeout(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload || {})
    }, API_TIMEOUT_MS.tp53);
    if (!response.ok) {
        const text = await response.text();
        throw new Error(`TP53 endpoint failed (${response.status}): ${text}`);
    }
    return response.json();
}

// Extract the numeric coordinate from a cDNA string (e.g. "c.1799T>A" -> 1799). If no
// numeric coordinate can be parsed, returns null. This helper ignores intronic
// suffixes (e.g. "+43", "-12") and simply extracts the first integer following
// the "c." prefix. When ranges are present (e.g. "c.178_186del"), the start
// coordinate is returned.
function parseCdnaCoordinate(cdna) {
    if (!cdna) return null;
    const m = String(cdna).match(/c\.\s*(-?\d+)/i);
    if (m) {
        const num = parseInt(m[1], 10);
        return isNaN(num) ? null : num;
    }
    return null;
}

// Extract the numeric residue position from a protein change string (e.g. "p.Val600Glu"
// or "Val600Glu" -> 600). Returns null when no position is found. The function
// ignores prefixes such as "p." and accepts both single and triple letter
// amino acid codes preceding and following the numeric portion.
function parseProteinCoordinate(prot) {
    if (!prot) return null;
    const p = String(prot).replace(/^p\./i, '');
    const m = p.match(/\D(\d+)/);
    if (m) {
        const num = parseInt(m[1], 10);
        return isNaN(num) ? null : num;
    }
    return null;
}

// Given an array of candidate transcripts and an optional target protein (triple‑coded
// uppercase string from the user's query), compute the canonical transcript index.
// Each candidate must have properties: transcript (string), cDNA (string), protein
// (string), and source (one of 'snpeff', 'dbnsfp', 'root'). The function returns
// the index of the canonical candidate using a weighted scoring system:
//   * Prefer RefSeq NM_ transcripts (assign high base score and rank by numeric accession
//     with smaller numbers ranked higher).
//   * Prefer candidates whose cDNA coordinate and protein residue position are
//     closest to the median of all candidate coordinates (to reflect widely used
//     isoforms when multiple variants exist).
//   * Prefer candidates that include a protein annotation over those that do not.
//   * Apply a minor source priority: snpEff > dbNSFP > root‑level.
//   * When targetProtGlobal is provided, strongly favour candidates whose protein
//     change matches the target (either triple‑coded or converted to single letter).
function selectCanonicalTranscript(candidates, targetProtGlobal) {
    if (!candidates || candidates.length === 0) return 0;
    // Precompute numeric coordinates for cDNA and protein to calculate medians.
    const cdnaCoords = [];
    const protCoords = [];
    for (const c of candidates) {
        const cc = parseCdnaCoordinate(c.cDNA);
        if (cc !== null) cdnaCoords.push(cc);
        const pc = parseProteinCoordinate(c.protein);
        if (pc !== null) protCoords.push(pc);
    }
    // Compute median values when possible.
    const median = (arr) => {
        if (arr.length === 0) return null;
        const sorted = [...arr].sort((a, b) => a - b);
        const idx = Math.floor(sorted.length / 2);
        return sorted[idx];
    };
    const medianCdna = median(cdnaCoords);
    const medianProt = median(protCoords);
    // Normalise NM accession numbers to numeric part. Returns null if not an NM transcript.
    const getNmNumber = (tx) => {
        if (!tx || !/^NM_/i.test(tx)) return null;
        const m = tx.match(/^NM_0*([0-9]+)(?:\.|$)/i);
        if (m) {
            const num = parseInt(m[1], 10);
            return isNaN(num) ? null : num;
        }
        return null;
    };
    // Source priority mapping (higher value means preferred). snpEff (3) > dbNSFP (2) > root (1)
    const sourcePriority = { snpeff: 3, dbnsfp: 2, root: 1 };
    // Compute scores for each candidate.
    let bestIndex = 0;
    let bestScore = -Infinity;
    candidates.forEach((c, idx) => {
        let score = 0;
        // High base score for NM transcripts. Smaller accession numbers yield higher score.
        const nmNum = getNmNumber(c.transcript);
        if (nmNum !== null) {
            // Use a million base to emphasise NM transcripts. Subtract nmNum so that smaller numbers are preferred.
            score += 1_000_000 - nmNum;
        }
        // Add moderate bonus when a protein change exists.
        if (c.protein) {
            score += 10_000;
        }
        // Subtract distance from median cDNA coordinate (if both exist). Smaller distance increases score.
        const cc = parseCdnaCoordinate(c.cDNA);
        if (medianCdna !== null && cc !== null) {
            const dist = Math.abs(cc - medianCdna);
            score -= dist;
        }
        // Subtract distance from median protein coordinate (if both exist).
        const pc = parseProteinCoordinate(c.protein);
        if (medianProt !== null && pc !== null) {
            const dist = Math.abs(pc - medianProt);
            score -= dist;
        }
        // Apply source priority multiplier. Higher priority sources add more to the score.
        const pri = sourcePriority[c.source] || 0;
        score += pri * 100;
        // If target protein is specified, greatly boost candidates matching it.
        if (targetProtGlobal && c.protein) {
            // Compare triple‑coded and single‑letter representations.
            // Remove the leading p. and non-alphanumeric characters from the protein string.
            const prot = String(c.protein)
                .replace(/^p\./i, '')
                .replace(/[^A-Za-z0-9]/g, '')
                .toUpperCase();
            const triple = prot;            // e.g. VAL600GLU
            const single = tripleToSingle(prot); // e.g. V600E
            // A match occurs when either the triple‑coded candidate or the single‑letter
            // form equals the target protein (also in triple‑coded form).
            if (triple === targetProtGlobal || (single && single === targetProtGlobal)) {
                score += 100_000_000; // huge bonus to ensure match wins
            }
        }
        // Keep track of the highest scoring candidate.
        if (score > bestScore) {
            bestScore = score;
            bestIndex = idx;
        }
    });
    return bestIndex;
}

// Build a list of transcript candidates from the MyVariant.info annotation and select
// the canonical transcript using selectCanonicalTranscript(). Returns an object
// containing a transcriptsList (with canonical flag), formatted cDNAHTML,
// formatted proteinHTML, and the canonical protein string. Only used for
// genomic variants (when Ensembl recoder transcripts are not available).
function buildCanonicalFromAnnotation(annotation, targetProtGlobal) {
    const candidates = [];
    if (!annotation) return { transcriptsList: [], cDNAHTML: '', proteinHTML: '', canonicalProtein: '' };
    // Gather candidates from snpEff annotations (preferred source).
    if (annotation.snpeff && annotation.snpeff.ann) {
        const annList = Array.isArray(annotation.snpeff.ann) ? annotation.snpeff.ann : [annotation.snpeff.ann];
        for (const ann of annList) {
            if (ann && ann.hgvs_c) {
                const txId = ann.feature_id || '';
                const cDNAVal = ann.hgvs_c;
                const protVal = ann.hgvs_p || '';
                candidates.push({ transcript: txId, cDNA: cDNAVal, protein: protVal, source: 'snpeff' });
            }
        }
    }
    // Gather candidates from dbNSFP hgvsc/hgvsp arrays.
    if (annotation.dbnsfp && annotation.dbnsfp.hgvsc) {
        const hgvsc = Array.isArray(annotation.dbnsfp.hgvsc) ? annotation.dbnsfp.hgvsc : [annotation.dbnsfp.hgvsc];
        const hgvsp = annotation.dbnsfp.hgvsp ? (Array.isArray(annotation.dbnsfp.hgvsp) ? annotation.dbnsfp.hgvsp : [annotation.dbnsfp.hgvsp]) : [];
        for (let i = 0; i < hgvsc.length; i++) {
            const sc = hgvsc[i];
            if (!sc) continue;
            const parts = String(sc).split(':');
            const txId = parts[0];
            const cpart = parts.slice(1).join(':');
            let ppart = '';
            if (hgvsp[i]) {
                const pparts = String(hgvsp[i]).split(':');
                ppart = pparts.slice(1).join(':');
            }
            candidates.push({ transcript: txId, cDNA: cpart, protein: ppart, source: 'dbnsfp' });
        }
    }
    // Gather candidates from root-level hgvsc/hgvsp fields. Only include those not already captured above.
    if (annotation.hgvsc) {
        const hgvscList = Array.isArray(annotation.hgvsc) ? annotation.hgvsc : [annotation.hgvsc];
        const hgvspList = annotation.hgvsp ? (Array.isArray(annotation.hgvsp) ? annotation.hgvsp : [annotation.hgvsp]) : [];
        for (let i = 0; i < hgvscList.length; i++) {
            const h = hgvscList[i];
            if (!h) continue;
            const parts = String(h).split(':');
            const txId = parts[0] || '';
            const cpart = parts.slice(1).join(':');
            let ppart = '';
            if (hgvspList[i]) {
                const pparts = String(hgvspList[i]).split(':');
                ppart = pparts.slice(1).join(':');
            }
            candidates.push({ transcript: txId, cDNA: cpart, protein: ppart, source: 'root' });
        }
    }
    if (candidates.length === 0) {
        return { transcriptsList: [], cDNAHTML: '', proteinHTML: '', canonicalProtein: '' };
    }
    // Deduplicate candidates by transcript ID + cDNA + protein. When duplicates exist from
    // multiple sources, retain only the candidate from the highest priority source.
    const unique = [];
    const seen = {};
    const sourceRank = { snpeff: 3, dbnsfp: 2, root: 1 };
    for (const cand of candidates) {
        const key = `${cand.transcript}|${cand.cDNA}|${cand.protein}`;
        if (!seen[key]) {
            seen[key] = cand;
        } else {
            // Replace existing candidate only if new one has higher source rank.
            const existing = seen[key];
            if ((sourceRank[cand.source] || 0) > (sourceRank[existing.source] || 0)) {
                seen[key] = cand;
            }
        }
    }
    for (const key in seen) {
        unique.push(seen[key]);
    }
    // Determine the canonical index among unique candidates.
    const canonicalIdx = selectCanonicalTranscript(unique, targetProtGlobal);
    // Mark canonical flag and build formatted strings.
    let canonicalProtein = '';
    const cDNAHTML = unique
        .map((c, idx) => {
            if (idx === canonicalIdx) {
                return `<strong>${c.cDNA}</strong>`;
            }
            return c.cDNA;
        })
        .join(', ');
    const proteinHTML = unique
        .map((c, idx) => {
            if (!c.protein) return '';
            if (idx === canonicalIdx) {
                return `<strong>${c.protein}</strong>`;
            }
            return c.protein;
        })
        .filter(Boolean)
        .join(', ');
    // Determine canonical protein string (if any). Prefer the protein from the canonical candidate
    // if present; otherwise use the first available protein among all candidates.
    if (unique[canonicalIdx] && unique[canonicalIdx].protein) {
        canonicalProtein = unique[canonicalIdx].protein;
    } else {
        const firstProt = unique.find(c => c.protein);
        if (firstProt) canonicalProtein = firstProt.protein;
    }
    // Set canonical property on unique candidates.
    unique.forEach((c, idx) => { c.canonical = (idx === canonicalIdx); });
    return { transcriptsList: unique, cDNAHTML, proteinHTML, canonicalProtein };
}

// Convert hgvsg notation (e.g. "NC_000001.11:g.230710048A>C" or "NC_000001.11:g.230710047_230710048delinsGA")
// to MyVariant.info notation ("chr1:g.230710048A>C" etc.).
function convertHgvsgToMyVariant(hgvsg) {
    // Split accession and variant part
    const [accession, remainder] = hgvsg.split(':g.');
    const chr = accessionToChr(accession);
    if (!chr) throw new Error(`Cannot map accession ${accession} to chromosome`);
    return `chr${chr}:g.${remainder}`;
}

// Convert SPDI notation (e.g. "NC_000001.11:230710047:A:C") to MyVariant.info notation.
function convertSpdiToMyVariant(spdi) {
    // Format: accession:position:ref:alt, 0-based coordinate
    const parts = spdi.split(':');
    if (parts.length !== 4) throw new Error(`Invalid SPDI: ${spdi}`);
    const [accession, posStr, ref, alt] = parts;
    const chr = accessionToChr(accession);
    if (!chr) throw new Error(`Cannot map accession ${accession} to chromosome`);
    const pos0 = parseInt(posStr, 10);
    // Convert to 1‑based coordinates. SPDI positions are 0‑based and refer to the start of the reference sequence.
    const start = pos0 + 1;
    // The SPDI format expresses simple substitutions, insertions and deletions. See:
    // https://www.ncbi.nlm.nih.gov/variation/notation/ for specification.
    // When the reference (ref) string is non‑empty and the alternate (alt) string is non‑empty, this is a substitution/multi‑nucleotide change.
    // When the ref is non‑empty and alt is empty, this is a deletion. The deleted region spans ref.length bases starting at the 0‑based position.
    // When ref is empty and alt is non‑empty, this is an insertion. The inserted sequence is inserted between the base at pos0 and pos0+1.
    // Construct MyVariant.info HGVS strings accordingly:
    if (ref && alt) {
        // Substitution or multi‑nucleotide change: chrN:g.startREF>ALT
        const start1 = start;
        return `chr${chr}:g.${start1}${ref}>${alt}`;
    }
    if (ref && !alt) {
        // Deletion: deletion of ref.length bases starting at start
        const end = start + ref.length - 1;
        // If a single base is deleted, use g.posdel; for multi‑base deletion use g.start_enddel
        if (ref.length === 1) {
            return `chr${chr}:g.${start}${ref}>-`;
        }
        return `chr${chr}:g.${start}_${end}del`;
    }
    if (!ref && alt) {
        // Insertion: insertion of alt after base at pos0. In HGVS this is start_(start+1)insALT
        const end = start; // end position for insertion is start position
        return `chr${chr}:g.${start}_${start + 1}ins${alt}`;
    }
    // If both ref and alt are empty (unlikely), fall back to simple representation
    return `chr${chr}:g.${start}>${alt}`;
}

async function fetchVariantRecoder(query) {
    const encoded = encodeURIComponent(query);
    // Use the Ensembl GRCh37 server for variant_recoder requests.  Complex
    // clinical variants are often catalogued relative to the hg19/GRCh37
    // assembly, and the default GRCh38 server may return 3' UTR or other
    // non‑canonical annotations.  Switching to grch37.rest.ensembl.org
    // improves consistency and allows detection of variants such as
    // BRAF c.1799_1811delinsA.
    const url = `https://grch37.rest.ensembl.org/variant_recoder/human/${encoded}?content-type=application/json`;
    const response = await fetchWithTimeout(url, {
        headers: {
            'Accept': 'application/json'
        }
    }, API_TIMEOUT_MS.recoder);
    if (!response.ok) {
        const text = await response.text();
        throw new Error(`Variant recoder request failed (${response.status}): ${text}`);
    }
    return response.json();
}

// Liftover a genomic variant from hg38 to hg19 using Ensembl REST API. If the
// conversion fails or the variant is already hg19, returns the original
// variant. Expects variant in form 'chr7:g.140753336A>T'.
async function liftoverHg38ToHg19(variant) {
    const m = variant.match(/^chr([\w]+):g\.(\d+)([A-Za-z-]+)>([A-Za-z-]+)/);
    if (!m) return variant;
    const chrom = m[1];
    const pos = parseInt(m[2], 10);
    const ref = m[3];
    const alt = m[4];
    // Build URL for liftover using Ensembl map endpoint
    const url = `https://rest.ensembl.org/map/human/GRCh38/${chrom}:${pos}..${pos}:1/GRCh37?content-type=application/json`;
    try {
        const res = await fetchWithTimeout(url, { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.liftover);
        if (!res.ok) {
            // If error (e.g. 400), just return original
            return variant;
        }
        const data = await res.json();
        if (data && data.mappings && data.mappings.length > 0) {
            const mapped = data.mappings[0].mapped;
            if (mapped && mapped.start) {
                const newPos = mapped.start;
                return `chr${chrom}:g.${newPos}${ref}>${alt}`;
            }
        }
    } catch (err) {
        // console.warn('Liftover error', err);
    }
    return variant;
}



// Fetch gnomAD v4.1 (GRCh38) data for a variant.
// Liftover GRCh37→GRCh38 is done server-side by the proxy.
// Returns { status, data?, grch38Id?, message? } or null on hard failure.
async function fetchGnomadV4(chrom, pos37, ref, alt) {
    const c = String(chrom).replace(/^chr/i, '');
    if (!c || !pos37 || !ref || !alt) return null;
    const endpoint = getConfiguredApiEndpoint('GNOMAD_V4_API_ENDPOINT', '/api/gnomad-v4');
    const params = new URLSearchParams({ chrom: c, pos37: String(pos37), ref: ref.toUpperCase(), alt: alt.toUpperCase() });
    const url = `${endpoint}?${params}`;
    try {
        const res = await fetchWithTimeout(url, {}, API_TIMEOUT_MS.gnomadV4);
        if (!res.ok) return null;
        return await res.json();
    } catch (err) {
        return null;
    }
}

async function fetchClinvarRegionVariants(chrom, pos, windowSize = 10) {
    const c = String(chrom || '').replace(/^chr/i, '').toUpperCase();
    const p = Number(pos);
    if (!c || !Number.isFinite(p)) return [];
    const endpoint = getConfiguredApiEndpoint('CLINVAR_REGION_API_ENDPOINT', '/api/clinvar-region');
    const params = new URLSearchParams({ chrom: c, pos: String(p), window: String(windowSize) });
    const res = await fetchWithTimeout(`${endpoint}?${params}`, {}, API_TIMEOUT_MS.clinvar);
    if (!res.ok) throw new Error(`ClinVar region proxy error: ${res.status}`);
    const data = await res.json();
    if (data?.error) throw new Error(data.error);
    return { variants: data.variants || [], total: data.total ?? (data.variants || []).length };
}

async function fetchClinvarGeneVariants(gene, retmax = 500) {
    const safeGene = String(gene || '').replace(/[^A-Za-z0-9\-_.]/g, '');
    if (!safeGene) return [];
    const endpoint = getConfiguredApiEndpoint('CLINVAR_GENE_API_ENDPOINT', '/api/clinvar-gene');
    const params = new URLSearchParams({ gene: safeGene, retmax: String(retmax) });
    const res = await fetchWithTimeout(`${endpoint}?${params}`, {}, API_TIMEOUT_MS.clinvar);
    if (!res.ok) throw new Error(`ClinVar gene proxy error: ${res.status}`);
    const data = await res.json();
    if (data?.error) throw new Error(data.error);
    return data.variants || [];
}

async function fetchClinvarVariant(variationId) {
    if (!variationId || !/^\d+$/.test(String(variationId))) return null;
    const endpoint = getConfiguredApiEndpoint('CLINVAR_VARIANT_API_ENDPOINT', '/api/clinvar-variant');
    const params = new URLSearchParams({ id: String(variationId) });
    const res = await fetchWithTimeout(`${endpoint}?${params}`, {}, API_TIMEOUT_MS.clinvar);
    if (!res.ok) throw new Error(`ClinVar variant proxy error: ${res.status}`);
    const data = await res.json();
    if (data?.error) throw new Error(data.error);
    return data;
}
function getClinvarVariantText(variant) {
    return [
        variant?.title,
        variant?.molecularConsequence,
        variant?.variationName
    ].filter(Boolean).join(' ');
}

function isSynonymousClinvarVariant(variant) {
    const text = getClinvarVariantText(variant).toLowerCase();
    return /\bsynonymous\b/.test(text) || /\bp\.\s*\(?=\)?/.test(text) || /\bp\.\s*[a-z]{1,3}\d+=/i.test(text);
}

function isTruncatingClinvarVariant(variant) {
    const text = getClinvarVariantText(variant).toLowerCase();
    return /\b(?:nonsense|frameshift|frame[- ]shift|stop[- ]gained|stop gained|stop_gained|protein truncating|truncating)\b/.test(text)
        || /\bp\.\s*\(?[a-z]{1,3}\d+(?:ter|\*|x)\)?/i.test(text)
        || /\bp\.\s*\(?[a-z]{1,3}\d+[a-z]{1,3}fs/i.test(text)
        || /\b(?:fs|ter|\*)\d*\b/i.test(text);
}

function isBelowAxisVariant(variant) {
    if (isSynonymousClinvarVariant(variant)) return true;
    const c = String(variant?.germline || '').toLowerCase();
    return c.includes('benign') && !c.includes('pathogenic');
}

function parseProteinPos(text) {
    // Handle both p.Arg273His and p.(Arg273His) (predicted-effect parenthesised form)
    const m = String(text || '').match(/\bp\.\(?(?:[A-Za-z*?]{1,5})?(\d+)/i);
    return m ? Number(m[1]) : null;
}

// Returns a color hex string for a ClinVar/CIViC pathogenicity classification.
function getPathogenicityColor(classification, variant = null) {
    if (variant && isTruncatingClinvarVariant(variant)) return '#111827';
    if (variant && isSynonymousClinvarVariant(variant)) return '#9ca3af';
    const c = String(classification || '').toLowerCase();
    if (c === 'pathogenic') return '#dc2626';
    if (c.includes('likely pathogenic')) return '#ef4444';
    if (c === 'benign') return '#16a34a';
    if (c.includes('likely benign')) return '#22c55e';
    if (c.includes('drug response')) return '#dc2626';
    if (c.includes('conflicting')) return '#f59e0b';
    if (c.includes('pathogenic')) return '#ef4444';
    if (c.includes('benign')) return '#22c55e';
    return '#f59e0b';
}

// Build an SVG lollipop plot for nearby ClinVar variants.
// variants: [{id, title, germline, pos}], queryPos: integer position
// minusStrand: if true, x-axis is oriented 3'→5' genomically (i.e. c./p. increases left-to-right)
// opts.range: half-width of the displayed window (default 30)
// opts.axisLabel: coordinate label shown at left of axis ('g.' or 'p.')
function buildLollipopPlot(variants, queryPos, minusStrand = false, { range = 30, axisLabel = 'g.' } = {}) {
    const NS = 'http://www.w3.org/2000/svg';
    const W = 280, ML = 18, MR = 18, PW = W - ML - MR;
    const RANGE = range;
    const STACK_SPACING = 13;
    const LOLLIPOP_OFFSET = 16;
    const MIN_AXIS_Y = 72;
    const TOP_PADDING = 20;
    const BOTTOM_PADDING = 23;
    const posToX = minusStrand
        ? (p) => ML + ((queryPos + RANGE - p) / (2 * RANGE)) * PW
        : (p) => ML + ((p - queryPos + RANGE) / (2 * RANGE)) * PW;
    const bucketForX = (x) => Math.round(x / 4) * 4;

    const plottableVariants = variants
        .map((v) => ({ ...v, x: Number.isFinite(Number(v.pos)) ? posToX(Number(v.pos)) : null, belowAxis: isBelowAxisVariant(v) }))
        .filter((v) => v.x !== null && v.x >= ML - 8 && v.x <= W - MR + 8);

    const aboveVars = plottableVariants.filter(v => !v.belowAxis);
    const belowVars = plottableVariants.filter(v => v.belowAxis);

    const countBuckets = (vs) => vs.reduce((acc, v) => {
        const b = bucketForX(v.x); acc[b] = (acc[b] || 0) + 1; return acc;
    }, {});
    const maxStackAbove = Math.max(1, ...Object.values(countBuckets(aboveVars)));
    const maxStackBelow = Math.max(0, ...Object.values(countBuckets(belowVars)));

    const AY = Math.max(MIN_AXIS_Y, TOP_PADDING + LOLLIPOP_OFFSET + (maxStackAbove - 1) * STACK_SPACING);
    const belowH = maxStackBelow > 0 ? LOLLIPOP_OFFSET + (maxStackBelow - 1) * STACK_SPACING + 14 : 0;
    const H = AY + Math.max(BOTTOM_PADDING, belowH);

    const bucketsAbove = {}, bucketsBelow = {};
    const stackHeightAbove = (x) => {
        const b = bucketForX(x);
        if (!bucketsAbove[b]) bucketsAbove[b] = 0;
        const h = bucketsAbove[b]; bucketsAbove[b] += STACK_SPACING; return h;
    };
    const stackHeightBelow = (x) => {
        const b = bucketForX(x);
        if (!bucketsBelow[b]) bucketsBelow[b] = 0;
        const h = bucketsBelow[b]; bucketsBelow[b] += STACK_SPACING; return h;
    };

    const svg = document.createElementNS(NS, 'svg');
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.setAttribute('width', '100%');
    svg.setAttribute('height', String(H));
    svg.style.display = 'block';
    svg.style.overflow = 'visible';

    // X axis line
    const axis = document.createElementNS(NS, 'line');
    axis.setAttribute('x1', String(ML)); axis.setAttribute('y1', String(AY));
    axis.setAttribute('x2', String(W - MR)); axis.setAttribute('y2', String(AY));
    axis.setAttribute('stroke', '#cbd5e1'); axis.setAttribute('stroke-width', '1');
    svg.appendChild(axis);

    // Axis coordinate label (g. or p.) at left of axis
    const axisLabelEl = document.createElementNS(NS, 'text');
    axisLabelEl.setAttribute('x', String(ML - 2)); axisLabelEl.setAttribute('y', String(AY + 4));
    axisLabelEl.setAttribute('text-anchor', 'end'); axisLabelEl.setAttribute('font-size', '7');
    axisLabelEl.setAttribute('fill', '#94a3b8'); axisLabelEl.textContent = axisLabel;
    svg.appendChild(axisLabelEl);

    // Subtle directional hints flanking the axis
    const upHint = document.createElementNS(NS, 'text');
    upHint.setAttribute('x', String(ML + 1)); upHint.setAttribute('y', String(AY - 3));
    upHint.setAttribute('font-size', '6'); upHint.setAttribute('fill', '#e2e8f0');
    upHint.textContent = '↑ path';
    svg.appendChild(upHint);
    const dnHint = document.createElementNS(NS, 'text');
    dnHint.setAttribute('x', String(ML + 1)); dnHint.setAttribute('y', String(AY + 22));
    dnHint.setAttribute('font-size', '6'); dnHint.setAttribute('fill', '#e2e8f0');
    dnHint.textContent = '↓ benign';
    svg.appendChild(dnHint);

    // Ticks and labels at -range, -range/2, 0, +range/2, +range
    [-RANGE, -RANGE / 2, 0, RANGE / 2, RANGE].forEach((o) => {
        const x = posToX(queryPos + (minusStrand ? -o : o));
        const tick = document.createElementNS(NS, 'line');
        tick.setAttribute('x1', String(x)); tick.setAttribute('y1', String(AY));
        tick.setAttribute('x2', String(x)); tick.setAttribute('y2', String(AY + 4));
        tick.setAttribute('stroke', '#94a3b8'); tick.setAttribute('stroke-width', '1');
        svg.appendChild(tick);
        const lbl = document.createElementNS(NS, 'text');
        lbl.setAttribute('x', String(x)); lbl.setAttribute('y', String(AY + 13));
        lbl.setAttribute('text-anchor', 'middle'); lbl.setAttribute('font-size', '7');
        lbl.setAttribute('fill', '#64748b');
        lbl.textContent = o === 0 ? '0' : (o > 0 ? `+${o}` : String(o));
        svg.appendChild(lbl);
    });

    // Query position marker: blue diamond on axis
    const qx = posToX(queryPos);
    const dia = document.createElementNS(NS, 'polygon');
    dia.setAttribute('points', `${qx},${AY - 5} ${qx + 4},${AY} ${qx},${AY + 5} ${qx - 4},${AY}`);
    dia.setAttribute('fill', '#3b82f6');
    svg.appendChild(dia);

    // Legend
    [['#111827', 'Truncating'], ['#dc2626', 'Pathogenic/LP'], ['#16a34a', 'Benign/LB'], ['#f59e0b', 'VUS/Other'], ['#9ca3af', 'Synonymous']].forEach(([col, lbl], i) => {
        const lx = 5 + i * 54;
        const lc = document.createElementNS(NS, 'circle');
        lc.setAttribute('cx', String(lx)); lc.setAttribute('cy', '7');
        lc.setAttribute('r', '4'); lc.setAttribute('fill', col);
        svg.appendChild(lc);
        const lt = document.createElementNS(NS, 'text');
        lt.setAttribute('x', String(lx + 7)); lt.setAttribute('y', '11');
        lt.setAttribute('font-size', '7.5'); lt.setAttribute('fill', '#374151');
        lt.textContent = lbl;
        svg.appendChild(lt);
    });

    // Plot each variant, stacking above or below the axis by pathogenicity
    plottableVariants.forEach((v) => {
        const x = v.x;
        const color = getPathogenicityColor(v.germline, v);
        let cy, stemY1, stemY2;
        if (v.belowAxis) {
            const sh = stackHeightBelow(x);
            cy = AY + LOLLIPOP_OFFSET + sh;
            stemY1 = AY;
            stemY2 = cy - 5;
        } else {
            const sh = stackHeightAbove(x);
            cy = AY - LOLLIPOP_OFFSET - sh;
            stemY1 = Math.max(cy + 5, 18);
            stemY2 = AY;
        }

        const stem = document.createElementNS(NS, 'line');
        stem.setAttribute('x1', String(x)); stem.setAttribute('y1', String(stemY1));
        stem.setAttribute('x2', String(x)); stem.setAttribute('y2', String(stemY2));
        stem.setAttribute('stroke', color); stem.setAttribute('stroke-width', '1.5');
        stem.setAttribute('opacity', '0.65');
        svg.appendChild(stem);

        const circ = document.createElementNS(NS, 'circle');
        circ.setAttribute('cx', String(x)); circ.setAttribute('cy', String(cy));
        circ.setAttribute('r', '4.5'); circ.setAttribute('fill', color);
        circ.setAttribute('opacity', '0.88');
        circ.setAttribute('style', 'cursor:pointer');
        const tip = document.createElementNS(NS, 'title');
        const consequenceLabel = isTruncatingClinvarVariant(v) ? 'Truncating; ' : (isSynonymousClinvarVariant(v) ? 'Synonymous; ' : '');
        tip.textContent = `${consequenceLabel}${v.germline || 'Unknown'} (pos ${v.pos}): ${v.title || v.id}`;
        circ.appendChild(tip);
        svg.appendChild(circ);
    });

    // Note for variants without position
    const noPos = variants.filter((v) => v.pos === null || v.pos === undefined || !Number.isFinite(Number(v.pos)));
    if (noPos.length > 0) {
        const note = document.createElementNS(NS, 'text');
        note.setAttribute('x', String(W - MR)); note.setAttribute('y', String(H - 2));
        note.setAttribute('text-anchor', 'end'); note.setAttribute('font-size', '6.5');
        note.setAttribute('fill', '#94a3b8');
        note.textContent = `${noPos.length} variant${noPos.length > 1 ? 's' : ''} without position not shown`;
        svg.appendChild(note);
    }

    return svg;
}

// Build a protein-position lollipop plot from the same nearby ClinVar variants.
// Protein positions are derived from the canonical c. number + genomic offset when
// possible (transcript-agnostic), falling back to parsing p. from the variant title.
// queryGenomicPos / queryC / minusStrand enable the genomic-offset path.
// The axis range auto-scales to the actual spread of the data.
// Returns null if no variants have parseable protein positions.
function buildProteinLollipopPlot(variants, queryProteinPos, { queryGenomicPos = null, queryC = null, minusStrand = false } = {}) {
    const strandFactor = minusStrand ? -1 : 1;
    const proteinVariants = variants.map(v => {
        let pos = null;
        if (v.pos != null && queryGenomicPos != null && queryC != null) {
            const cCanonical = queryC + strandFactor * (v.pos - queryGenomicPos);
            if (cCanonical > 0) {
                const p = Math.ceil(cCanonical / 3);
                // Only trust the genomic offset when it gives a plausible same-exon result.
                // A large deviation signals an intron or exon boundary in between; fall through
                // to title-parsed p. instead (which, though transcript-specific, is better than
                // placing the variant hundreds of residues off).
                if (Math.abs(p - queryProteinPos) <= 20) pos = p;
            }
        }
        if (pos === null) pos = parseProteinPos(getClinvarVariantText(v));
        return { ...v, pos };
    });
    const withPos = proteinVariants.filter(v => v.pos !== null);
    if (withPos.length === 0) return null;

    // Round up to a tidy number so tick labels land on clean values, capped at ±10.
    const niceCeil = (n) => {
        if (n <= 5) return 5;
        const step = n <= 20 ? 5 : n <= 50 ? 10 : 25;
        return Math.min(10, Math.ceil(n / step) * step);
    };
    const maxDev = Math.max(...withPos.map(v => Math.abs(v.pos - queryProteinPos)));
    const range = niceCeil(Math.max(5, maxDev + 2));

    return buildLollipopPlot(proteinVariants, queryProteinPos, false, { range, axisLabel: 'p.' });
}

// Query CivicDB via the /api/civic serverless proxy (avoids browser CORS).
async function fetchCivicApiData(geneName, proteinChange) {
    if (!geneName) return null;
    const safeGene = String(geneName).replace(/[^A-Za-z0-9\-_.]/g, '');
    if (!safeGene) return null;
    try {
        const params = new URLSearchParams({ gene: safeGene });
        if (proteinChange) params.set('protein', String(proteinChange).replace(/^p\./i, ''));
        const endpoint = getConfiguredApiEndpoint('CIVIC_API_ENDPOINT', '/api/civic');
        const res = await fetchWithTimeout(appendQueryParams(endpoint, params), {}, API_TIMEOUT_MS.civic);
        if (!res.ok) return null;
        const data = await res.json();
        if (data?.error) return null;
        return data; // { gene, matchedVariant, assertions }
    } catch (e) {
        console.warn('CivicDB API fetch failed', e);
        return null;
    }
}

/**
 * Query the FDA companion diagnostics API for structured drug/disease/biomarker rows
 * for the given gene. Returns the raw records array from the API response.
 */
async function fetchFdaCompanionDiagnostics(gene) {
    if (!gene) return [];
    const params = new URLSearchParams({ gene });
    const url = `https://drdoubleb.com/BBKB/gene_api.php?${params}`;
    let resp;
    try {
        resp = await fetchWithTimeout(url, {}, API_TIMEOUT_MS.fda);
    } catch (e) {
        console.warn('FDA companion diagnostics fetch failed', e);
        return [];
    }
    if (!resp.ok) throw new Error(`FDA companion diagnostics API error: ${resp.status}`);
    const data = await resp.json();
    return data.records || [];
}

// ── openFDA gene-negation filter ─────────────────────────────────────────
// Many FDA oncology labels (esp. checkpoint inhibitors) carve out genes in
// the *exclusion* criteria of their indication, e.g.
//   "...with no EGFR or ALK genomic tumor aberrations"
//   "...with no sensitizing epidermal growth factor receptor (EGFR) mutation
//      or anaplastic lymphoma kinase (ALK) genomic tumor aberrations"
//   "...with no known EGFR mutations or ALK rearrangements"
//   "...for BRAF wild-type melanoma"
//   "Patients with EGFR or ALK genomic tumor aberrations should have disease
//      progression on FDA-approved therapy for these aberrations prior to
//      receiving X"   (sequencing carve-out used by pembrolizumab,
//                      atezolizumab, nivolumab, ramucirumab, etc.)
// A naive substring match treats these as positive mentions of the gene
// and floods AI review with non-matches. findOpenFdaNegationSpans() returns
// character ranges covering such phrases; openFdaGeneOnlyInNegativeContext()
// excludes a record only when *every* occurrence of the queried gene falls
// inside one of those spans, so a label that also mentions the gene
// positively elsewhere (e.g. CYRAMZA + erlotinib for EGFR exon 19) is kept.
//
// Spans are bounded by sentence terminators (. ; :) and short lazy windows
// so a trigger word in one sentence can't accidentally negate a positive
// mention several sentences later.
//
// KEEP IN SYNC with tests/openfda-filter.test.js — no module system in the
// browser, so the test file re-declares the same regexes.
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

// ── openFDA word-boundary filter ─────────────────────────────────────────
// openFDA's search engine is a free-text substring match, so a query for the
// gene "MET" also returns labels where the three letters merely sit inside a
// larger word — e.g. the diabetes brand AVANDAMET, the drug METHOTREXATE, or
// the homeopathic ingredient METALICUM ("Argentum metalicum"). None of those
// are the MET oncogene. The API has no regex support, so we can't require a
// word boundary server-side; instead we keep a record only when at least one
// occurrence of the gene is flanked by non-letters. Punctuation/space
// boundaries such as "(MET)", "MET exon 14" and "MET-amplified" still qualify,
// so true targeted-therapy labels (capmatinib, amivantamab, …) are retained.
// KEEP IN SYNC with tests/openfda-filter.test.js.
function openFdaEscapeRegExp(s) {
    return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
function openFdaGeneAppearsAsWord(text, gene) {
    if (!text || !gene) return false;
    return new RegExp(String.raw`(?<![A-Za-z])${openFdaEscapeRegExp(gene)}(?![A-Za-z])`).test(text);
}

async function fetchOpenFdaDrugLabels(gene) {
    if (!gene) return { total: 0, fetched: 0, excluded: 0, excludedCase: 0, excludedBoundary: 0, excludedNegation: 0, results: [] };
    const LIMIT = 100;
    const MAX_RECORDS = 1000;
    // Restrict to prescription drugs server-side. The vast majority of false
    // positives (homeopathic remedies, OTC monographs whose ingredient names
    // contain the gene's letters) are HUMAN OTC DRUG / unapproved homeopathic
    // products, so this clause removes them before they consume our paging
    // budget — important when a gene like "MET" otherwise matches 600+ labels
    // and real targeted therapies could fall past the MAX_RECORDS cap. All FDA
    // targeted oncology therapies and biologics carry product_type
    // "HUMAN PRESCRIPTION DRUG", so this does not drop labels we want.
    const PRESCRIPTION_FILTER = '+AND+openfda.product_type:%22HUMAN+PRESCRIPTION+DRUG%22';
    const baseSearch = `indications_and_usage:%22${encodeURIComponent(gene)}%22${PRESCRIPTION_FILTER}`;

    const fetchPage = async (skip) => {
        const url = `https://api.fda.gov/drug/label.json?search=${baseSearch}&limit=${LIMIT}${skip ? `&skip=${skip}` : ''}`;
        const resp = await fetchWithTimeout(url, {}, API_TIMEOUT_MS.openFda);
        if (!resp.ok) {
            if (resp.status === 404) return null;
            throw new Error(`openFDA API error: ${resp.status}`);
        }
        return resp.json();
    };

    try {
        const firstData = await fetchPage(0);
        if (!firstData) return { total: 0, fetched: 0, excluded: 0, excludedCase: 0, excludedBoundary: 0, excludedNegation: 0, results: [] };

        const total = firstData?.meta?.results?.total || 0;
        let raw = firstData?.results || [];

        if (total > LIMIT) {
            const pagesToFetch = Math.ceil(Math.min(total, MAX_RECORDS) / LIMIT);
            for (let page = 1; page < pagesToFetch; page++) {
                const data = await fetchPage(page * LIMIT);
                if (data?.results) raw = raw.concat(data.results);
            }
        }

        let excludedCase = 0;
        let excludedBoundary = 0;
        let excludedNegation = 0;
        const results = [];
        for (const item of raw) {
            const ind = item.indications_and_usage?.[0] || '';
            if (!ind.includes(gene)) { excludedCase++; continue; }
            if (!openFdaGeneAppearsAsWord(ind, gene)) { excludedBoundary++; continue; }
            if (openFdaGeneOnlyInNegativeContext(ind, gene)) { excludedNegation++; continue; }
            results.push({
                brand_name: item.openfda?.brand_name?.[0] || '',
                generic_name: item.openfda?.generic_name?.[0] || '',
                manufacturer: item.openfda?.manufacturer_name?.[0] || '',
                route: (item.openfda?.route || []).join(', '),
                application_number: item.openfda?.application_number?.[0] || '',
                indications_and_usage: ind,
                purpose: item.purpose?.[0] || '',
            });
        }
        return {
            total,
            fetched: raw.length,
            excluded: excludedCase + excludedBoundary + excludedNegation,
            excludedCase,
            excludedBoundary,
            excludedNegation,
            results,
        };
    } catch (e) {
        console.warn('openFDA fetch failed', e);
        throw e;
    }
}

async function fetchPubmedArticles(searchTerm, limit = 5) {
    if (!searchTerm) return { total: 0, articles: [] };
    const params = new URLSearchParams({ term: searchTerm, limit: String(limit) });
    const endpoint = getConfiguredApiEndpoint('PUBMED_API_ENDPOINT', '/api/pubmed');
    const res = await fetchWithTimeout(appendQueryParams(endpoint, params), {}, API_TIMEOUT_MS.pubmed);
    if (!res.ok) throw new Error(`PubMed proxy error: ${res.status}`);
    const data = await res.json();
    if (data?.error) throw new Error(data.error);
    return { total: data.total || 0, articles: data.articles || [] };
}

async function fetchClinicalTrials(gene, tumorType) {
    if (!gene) return { total: 0, studies: [] };
    const params = new URLSearchParams({ gene });
    if (tumorType) params.set('tumorType', tumorType);
    const endpoint = getConfiguredApiEndpoint('CLINICAL_TRIALS_API_ENDPOINT', '/api/clinicaltrials');
    try {
        const res = await fetchWithTimeout(appendQueryParams(endpoint, params), {}, API_TIMEOUT_MS.clinicalTrials);
        if (!res.ok) return { total: 0, studies: [] };
        const data = await res.json();
        if (data?.error) { console.warn('ClinicalTrials.gov proxy error:', data.error); return { total: 0, studies: [] }; }
        return { total: data.total || 0, studies: data.studies || [], searchTerm: data.searchTerm || gene };
    } catch (e) {
        console.warn('fetchClinicalTrials failed', e);
        return { total: 0, studies: [] };
    }
}

const OPENROUTER_MODEL_OPTIONS = [
    'openai/gpt-5-mini',
    'openai/gpt-5.4-nano',
    'openai/gpt-oss-120b',
    'openai/gpt-4o-mini',
    'openai/gpt-5-nano',
    'google/gemini-2.5-flash-lite',
    'google/gemini-3-flash-preview',
    'google/gemini-3.1-flash-lite',
    'anthropic/claude-3-haiku',
    'deepseek/deepseek-v4-flash',
    'deepseek/deepseek-v4-pro',
    'x-ai/grok-4.3',
    'minimax/minimax-m2.7',
    'nvidia/nemotron-3-super-120b-a12b:free',
    'qwen/qwen3-32b'
];

async function fetchAiReview(context, model) {
    const endpoint = getConfiguredApiEndpoint('AI_REVIEW_API_ENDPOINT', '/api/ai-review');
    const res = await fetchWithTimeout(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ model, context })
    }, API_TIMEOUT_MS.aiReview);
    const data = await res.json().catch(() => ({}));
    if (!res.ok || data?.error) {
        const detail = data?.detail ? `: ${data.detail}` : '';
        throw new Error(`${data?.error || `AI review request failed (${res.status})`}${detail}`);
    }
    return data;
}

function appendListOrEmpty(parent, values, emptyText) {
    if (!Array.isArray(values) || values.length === 0) {
        const empty = document.createElement('div');
        empty.className = 'ai-review-empty';
        empty.textContent = emptyText;
        parent.appendChild(empty);
        return;
    }
    const ul = document.createElement('ul');
    ul.className = 'ai-review-list';
    values.forEach((value) => {
        const li = document.createElement('li');
        if (value && typeof value === 'object') {
            const title = value.drug || value.nct_id || value.nctId || value.title || value.intervention || 'Item';
            const bits = [
                value.indication,
                value.biomarker_context,
                value.phase,
                value.intervention,
                value.relevance,
                value.evidence
            ].filter(Boolean);
            if (value.url) {
                const a = document.createElement('a');
                a.href = value.url;
                a.target = '_blank';
                a.rel = 'noopener noreferrer';
                a.textContent = title;
                li.appendChild(a);
            } else {
                const strong = document.createElement('strong');
                strong.textContent = title;
                li.appendChild(strong);
            }
            if (bits.length) li.appendChild(document.createTextNode(` — ${bits.join(' · ')}`));
        } else {
            li.textContent = String(value);
        }
        ul.appendChild(li);
    });
    parent.appendChild(ul);
}

function renderAiReview(review, targetEl) {
    targetEl.innerHTML = '';
    const disclaimer = document.createElement('div');
    disclaimer.className = 'ai-review-disclaimer';
    disclaimer.textContent = 'AI-generated research summary only. Verify with current FDA labeling, guidelines, curated databases, and trial eligibility criteria before any clinical use.';
    targetEl.appendChild(disclaimer);

    const badges = document.createElement('div');
    badges.className = 'ai-review-badges';
    [['Pathogenicity', review?.pathogenicity], ['AMP Tier', review?.amp_tier]].forEach(([label, value]) => {
        const badge = document.createElement('span');
        badge.className = 'ai-review-badge';
        badge.textContent = `${label}: ${value || 'Not provided'}`;
        badges.appendChild(badge);
    });
    targetEl.appendChild(badges);

    const summary = document.createElement('p');
    summary.className = 'ai-review-summary';
    summary.textContent = review?.summary || 'No AI summary returned.';
    targetEl.appendChild(summary);

    if (review?.amp_tier_rationale) {
        const rationale = document.createElement('p');
        rationale.className = 'ai-review-rationale';
        rationale.textContent = `AMP tier rationale: ${review.amp_tier_rationale}`;
        targetEl.appendChild(rationale);
    }

    const therapiesTitle = document.createElement('h4');
    therapiesTitle.textContent = 'FDA-approved therapies';
    targetEl.appendChild(therapiesTitle);
    appendListOrEmpty(targetEl, review?.fda_approved_therapies, 'No FDA-approved therapies were returned by the AI review for this context.');

    const trialsTitle = document.createElement('h4');
    trialsTitle.textContent = 'Clinical trials';
    targetEl.appendChild(trialsTitle);
    appendListOrEmpty(targetEl, review?.clinical_trials, 'No clinical trials were returned by the AI review for this context.');

    if (Array.isArray(review?.limitations) && review.limitations.length > 0) {
        const limitationsTitle = document.createElement('h4');
        limitationsTitle.textContent = 'Limitations';
        targetEl.appendChild(limitationsTitle);
        appendListOrEmpty(targetEl, review.limitations, '');
    }
}


function buildSpliceAiApiVariant(rawInput, gVariant, annotation) {
    const tuple = buildSpliceAiLookupTuple(rawInput, gVariant);
    const vcfData = annotation && annotation.vcf;
    if (!tuple) return '';
    const chrom = tuple.chrom.replace(/^chr/i, 'chr');
    const ref = tuple.ref || (vcfData && String(vcfData.ref || '').toUpperCase()) || null;
    const alt = tuple.alt || (vcfData && String(vcfData.alt || '').toUpperCase()) || null;
    if (!chrom || !tuple.pos || !ref || !alt) return '';
    return `${chrom}-${tuple.pos}-${ref}-${alt}`;
}

async function fetchSpliceAiPrediction(variant, options = {}) {
    if (!variant) return null;
    const params = new URLSearchParams({
        variant,
        hg: String(options.hg || '37'),
        distance: String(options.distance || 500),
        mask: String(options.mask ?? 0),
        bc: String(options.bc || 'basic')
    });
    const endpoint = getConfiguredApiEndpoint('SPLICEAI_API_ENDPOINT', '/api/spliceai');
    const res = await fetchWithTimeout(appendQueryParams(endpoint, params), {}, API_TIMEOUT_MS.spliceai);
    if (!res.ok) throw new Error(`SpliceAI proxy error: ${res.status}`);
    const payload = await res.json();
    if (payload?.error) throw new Error(payload.error);
    return payload;
}

function getSpliceAiScoreSummary(payload) {
    const scores = Array.isArray(payload?.data?.scores) ? payload.data.scores : [];
    const deltaKeys = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'];
    const positionKeys = { DS_AG: 'DP_AG', DS_AL: 'DP_AL', DS_DG: 'DP_DG', DS_DL: 'DP_DL' };
    let best = null;
    const transcripts = scores.map((score) => {
        const deltas = deltaKeys.map((key) => {
            const value = Number(score[key]);
            return {
                key,
                label: key.replace('DS_', ''),
                value: Number.isFinite(value) ? value : null,
                position: score[positionKeys[key]] ?? null
            };
        });
        const transcriptBest = deltas.reduce((acc, item) => {
            if (item.value === null) return acc;
            if (!acc || item.value > acc.value) return item;
            return acc;
        }, null);
        const row = {
            transcript: score.t_id || score.NAME || '',
            gene: score.gene_name || score.g_name || score.SYMBOL || '',
            strand: score.t_strand || score.STRAND || '',
            priority: score.t_priority || '',
            best: transcriptBest,
            deltas
        };
        if (transcriptBest && (!best || transcriptBest.value > best.value)) {
            best = { ...transcriptBest, transcript: row.transcript, gene: row.gene, priority: row.priority };
        }
        return row;
    });
    return { best, transcripts };
}

async function fetchMyVariant(variant) {
    const encoded = encodeURIComponent(variant);
    const url = `https://myvariant.info/v1/variant/${encoded}`;
    const response = await fetchWithTimeout(url, {}, API_TIMEOUT_MS.myvariant);
    if (!response.ok) {
        const text = await response.text();
        throw new Error(`MyVariant.info request failed (${response.status}): ${text}`);
    }
    return response.json();
}

// Fetch variant consequence information from the Ensembl VEP HGVS endpoint on the GRCh37 server.
// Given a genomic HGVS string (e.g. "chr14:g.95583003del"), this helper removes the
// leading "chr" and calls the GRCh37 VEP HGVS endpoint. It returns an array of
// transcript consequences (if any) from the response. If the request fails or no
// consequences are found, it throws an error. We use GRCh37 here because many
// clinical variant coordinates (e.g. hg19) are based on this assembly.


// Resolve a reliable gene symbol from VEP consequences. Prefer explicit non-chromosome
// symbols from VEP, then resolve Ensembl gene IDs (ENSG...) via lookup, and finally
// fall back to transcript-gene mapping from the variant_recoder payload.
async function resolveGeneSymbolFromVep(consequences, recoderData) {
    if (!Array.isArray(consequences) || consequences.length === 0) return '';

    // 1) Prefer any direct gene_symbol that is not chromosome-like.
    for (const c of consequences) {
        const sym = c?.gene_symbol ? String(c.gene_symbol).trim() : '';
        if (sym && !isChromosomeLikeGeneSymbol(sym)) return sym;
    }

    // 2) Resolve via Ensembl gene IDs when available.
    const geneIds = Array.from(new Set(
        consequences
            .map(c => (c?.gene_id ? String(c.gene_id).trim() : ''))
            .filter(id => /^ENSG\d+/i.test(id))
    ));
    for (const geneId of geneIds) {
        try {
            const url = `https://grch37.rest.ensembl.org/lookup/id/${encodeURIComponent(geneId)}?content-type=application/json`;
            const resp = await fetchWithTimeout(url, { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.lookup);
            if (!resp.ok) continue;
            const data = await resp.json();
            const sym = data?.display_name ? String(data.display_name).trim() : '';
            if (sym && !isChromosomeLikeGeneSymbol(sym)) return sym;
        } catch {
            // ignore and continue other IDs
        }
    }

    // 3) Derive gene from transcript IDs via variant_recoder, then resolve gene symbol.
    const txIds = Array.from(new Set(
        consequences
            .map(c => (c?.transcript_id ? String(c.transcript_id).trim() : ''))
            .filter(Boolean)
    ));
    if (Array.isArray(recoderData) && recoderData.length > 0 && txIds.length > 0) {
        const recEntry = recoderData[0];
        for (const subVal of Object.values(recEntry || {})) {
            if (!subVal || typeof subVal !== 'object') continue;
            const hgvscList = Array.isArray(subVal.hgvsc) ? subVal.hgvsc : (subVal.hgvsc ? [subVal.hgvsc] : []);
            for (const sc of hgvscList) {
                const tx = String(sc).split(':')[0].trim();
                if (!txIds.includes(tx)) continue;
                if (/^ENST\d+/i.test(tx)) {
                    try {
                        const txUrl = `https://grch37.rest.ensembl.org/lookup/id/${encodeURIComponent(tx)}?content-type=application/json`;
                        const txResp = await fetchWithTimeout(txUrl, { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.lookup);
                        if (!txResp.ok) continue;
                        const txData = await txResp.json();
                        const parentGeneId = txData?.Parent ? String(txData.Parent).trim() : '';
                        if (!/^ENSG\d+/i.test(parentGeneId)) continue;
                        const gUrl = `https://grch37.rest.ensembl.org/lookup/id/${encodeURIComponent(parentGeneId)}?content-type=application/json`;
                        const gResp = await fetchWithTimeout(gUrl, { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.lookup);
                        if (!gResp.ok) continue;
                        const gData = await gResp.json();
                        const sym = gData?.display_name ? String(gData.display_name).trim() : '';
                        if (sym && !isChromosomeLikeGeneSymbol(sym)) return sym;
                    } catch {
                        // keep trying
                    }
                }
            }
        }
    }

    // 4) Last resort: keep first available non-empty symbol, even if chromosome-like.
    for (const c of consequences) {
        const sym = c?.gene_symbol ? String(c.gene_symbol).trim() : '';
        if (sym) return sym;
    }
    return '';
}

async function fetchVepHgvsHg19(variant) {
    if (!variant) throw new Error('No variant provided');
    // Normalize: strip leading "chr" prefix (case-insensitive) before the colon
    let hgvs = variant;
    const m = variant.match(/^chr([\w]+):g\.(.+)$/i);
    if (m) {
        hgvs = `${m[1]}:g.${m[2]}`;
    }
    const url = `https://grch37.rest.ensembl.org/vep/human/hgvs/${encodeURIComponent(hgvs)}?content-type=application/json`;
    const response = await fetchWithTimeout(url, {
        headers: {
            'Accept': 'application/json'
        }
    }, API_TIMEOUT_MS.vep);
    if (!response.ok) {
        const text = await response.text();
        throw new Error(`VEP HGVS request failed (${response.status}): ${text}`);
    }
    const data = await response.json();
    if (!Array.isArray(data) || data.length === 0) {
        throw new Error('No VEP HGVS data found');
    }
    // Extract transcript consequences list from the first result
    const first = data[0];
    const consequences = first.transcript_consequences || [];
    return { vepData: data, consequences };
}

function parseProteinTargetFromQueryVariant(variantPart) {
    if (!variantPart) return null;
    const aaMap = {
        ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
        HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
        THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V'
    };
    const cleaned = String(variantPart).replace(/^p\./i, '').replace(/[^A-Za-z0-9]/g, '');
    const mTriple = cleaned.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3})$/);
    if (mTriple) {
        const triple = `${mTriple[1]}${mTriple[2]}${mTriple[3]}`.toUpperCase();
        const refSingle = aaMap[mTriple[1].toUpperCase()];
        const altSingle = aaMap[mTriple[3].toUpperCase()];
        const single = (refSingle && altSingle) ? `${refSingle}${mTriple[2]}${altSingle}` : '';
        return { triple, single };
    }
    const mSingle = cleaned.match(/^([A-Za-z])(\d+)([A-Za-z])$/);
    if (mSingle) {
        return { triple: '', single: `${mSingle[1].toUpperCase()}${mSingle[2]}${mSingle[3].toUpperCase()}` };
    }
    return null;
}

function parseMyVariantQueryIntent(identifier) {
    const out = { gene: '', proteinTarget: null, cdnaTarget: '' };
    if (!identifier) return out;
    const s = String(identifier).trim();
    const m = s.match(/^([A-Za-z0-9_.-]+)\s*[:\s]\s*(.+)$/);
    if (!m) return out;
    out.gene = m[1].toUpperCase();
    const variantPart = m[2].trim();
    if (/^p\./i.test(variantPart)) {
        out.proteinTarget = parseProteinTargetFromQueryVariant(variantPart);
    } else if (/^c\./i.test(variantPart)) {
        out.cdnaTarget = variantPart.toLowerCase();
    }
    return out;
}

function getHitGenesUpper(hit) {
    const genes = new Set();
    const addGene = (g) => {
        if (!g) return;
        genes.add(String(g).trim().toUpperCase());
    };
    addGene(hit?.dbnsfp?.genename);
    addGene(hit?.clinvar?.gene?.symbol);
    if (hit?.snpeff?.ann) {
        const annList = Array.isArray(hit.snpeff.ann) ? hit.snpeff.ann : [hit.snpeff.ann];
        for (const ann of annList) {
            addGene(ann?.genename || ann?.gene_name);
        }
    }
    return genes;
}

function extractProteinNotationsFromHit(hit) {
    const vals = [];
    const pushVal = (v) => {
        if (!v) return;
        vals.push(String(v).trim());
    };
    if (hit?.dbnsfp?.hgvsp) {
        const list = Array.isArray(hit.dbnsfp.hgvsp) ? hit.dbnsfp.hgvsp : [hit.dbnsfp.hgvsp];
        for (const v of list) pushVal(v);
    }
    if (hit?.hgvsp) {
        const list = Array.isArray(hit.hgvsp) ? hit.hgvsp : [hit.hgvsp];
        for (const v of list) pushVal(v);
    }
    if (hit?.snpeff?.ann) {
        const annList = Array.isArray(hit.snpeff.ann) ? hit.snpeff.ann : [hit.snpeff.ann];
        for (const ann of annList) pushVal(ann?.hgvs_p);
    }
    return vals;
}

function extractCdnaNotationsFromHit(hit) {
    const vals = [];
    const pushVal = (v) => {
        if (!v) return;
        vals.push(String(v).trim().toLowerCase());
    };
    if (hit?.dbnsfp?.hgvsc) {
        const list = Array.isArray(hit.dbnsfp.hgvsc) ? hit.dbnsfp.hgvsc : [hit.dbnsfp.hgvsc];
        for (const v of list) pushVal(v);
    }
    if (hit?.hgvsc) {
        const list = Array.isArray(hit.hgvsc) ? hit.hgvsc : [hit.hgvsc];
        for (const v of list) pushVal(v);
    }
    if (hit?.snpeff?.ann) {
        const annList = Array.isArray(hit.snpeff.ann) ? hit.snpeff.ann : [hit.snpeff.ann];
        for (const ann of annList) pushVal(ann?.hgvs_c);
    }
    return vals;
}

// Query MyVariant.info by rsID or other variant identifier. Returns the best matching hit when possible.
async function queryMyVariantById(identifier) {
    const encoded = encodeURIComponent(identifier);
    const url = `https://myvariant.info/v1/query?q=${encoded}&size=20`;
    const response = await fetchWithTimeout(url, {}, API_TIMEOUT_MS.myvariant);
    if (!response.ok) {
        const text = await response.text();
        throw new Error(`MyVariant.info query request failed (${response.status}): ${text}`);
    }
    const data = await response.json();
    if (data && data.hits && data.hits.length > 0) {
        const intent = parseMyVariantQueryIntent(identifier);
        const scored = data.hits.map((hit, idx) => {
            let score = 0;
            const id = String(hit?._id || '').trim().toLowerCase();
            const rawId = String(identifier || '').trim().toLowerCase();
            if (id && rawId && id === rawId) score += 1_000_000;
            const hitGenes = getHitGenesUpper(hit);
            if (intent.gene && hitGenes.has(intent.gene)) score += 50_000;
            if (intent.proteinTarget) {
                const protVals = extractProteinNotationsFromHit(hit);
                for (const p of protVals) {
                    const parsed = parseProteinTargetFromQueryVariant(p);
                    if (!parsed) continue;
                    if (intent.proteinTarget.triple && parsed.triple && parsed.triple === intent.proteinTarget.triple) {
                        score += 500_000;
                    } else if (intent.proteinTarget.single && parsed.single && parsed.single === intent.proteinTarget.single) {
                        score += 350_000;
                    }
                }
            }
            if (intent.cdnaTarget) {
                const cdnaVals = extractCdnaNotationsFromHit(hit);
                if (cdnaVals.some(v => v.endsWith(`:${intent.cdnaTarget}`) || v === intent.cdnaTarget)) {
                    score += 300_000;
                }
            }
            // Prefer richer records on ties.
            if (hit?.clinvar) score += 10;
            if (hit?.dbnsfp) score += 10;
            return { idx, score, hit };
        });
        scored.sort((a, b) => b.score - a.score || a.idx - b.idx);
        return scored[0].hit;
    }
    return null;
}

// Build summary table rows from annotation object.
function buildSummary(annotation, variant) {
    const rows = [];
    // Variant ID
    rows.push({ name: 'Variant ID', value: annotation._id || variant });
    // Genes: collect from various sources
    const geneNames = new Set();
    // cadd.gene.genename may be string or object
    if (annotation.cadd?.gene) {
        const gname = annotation.cadd.gene.genename || annotation.cadd.gene.symbol;
        if (gname) geneNames.add(gname);
    }
    if (annotation.dbnsfp?.genename) {
        geneNames.add(annotation.dbnsfp.genename);
    }
    if (annotation.clinvar?.gene?.symbol) {
        geneNames.add(annotation.clinvar.gene.symbol);
    }
    if (annotation.civic?.gene?.name) {
        geneNames.add(annotation.civic.gene.name);
    }
    // Extract gene symbols from dbsnp gene list when available
    if (annotation.dbsnp?.gene) {
        const dbsnpGenes = Array.isArray(annotation.dbsnp.gene) ? annotation.dbsnp.gene : [annotation.dbsnp.gene];
        for (const g of dbsnpGenes) {
            if (g?.symbol) geneNames.add(g.symbol);
            else if (g?.genename) geneNames.add(g.genename);
        }
    }
    // Extract gene names from snpEff annotations
    if (annotation.snpeff?.ann) {
        const annList = Array.isArray(annotation.snpeff.ann) ? annotation.snpeff.ann : [annotation.snpeff.ann];
        for (const ann of annList) {
            if (ann?.genename) geneNames.add(ann.genename);
            else if (ann?.gene_name) geneNames.add(ann.gene_name);
        }
    }
    // unify genes, join with comma
    if (geneNames.size > 0) {
        rows.push({ name: 'Gene(s)', value: Array.from(geneNames).join(', ') });
    }
    // Type
    const type = annotation.cadd?.annotype || annotation.dbsnp?.vartype || annotation.type;
    if (type) rows.push({ name: 'Variant Type', value: type });
    // Consequence or effect
    const consequence = annotation.cadd?.consequence || annotation.exac?.functional_annotation;
    if (consequence) rows.push({ name: 'Consequence', value: consequence });
    // ClinVar significance
    let clinSig;
    if (annotation.clinvar?.rcv) {
        // rcv can be array or single object. Extract the clinical_significance field.
        const rcv = annotation.clinvar.rcv;
        const sigs = [];
        const extractSig = (rc) => {
            if (!rc) return;
            // rcv entries may define clinical_significance as a string or object. Prefer description when available.
            let cs = rc.clinical_significance;
            if (cs) {
                if (typeof cs === 'object' && cs.description) {
                    sigs.push(cs.description);
                } else {
                    sigs.push(cs);
                }
            }
        };
        if (Array.isArray(rcv)) {
            rcv.forEach(extractSig);
        } else if (typeof rcv === 'object') {
            extractSig(rcv);
        }
        if (sigs.length > 0) clinSig = Array.from(new Set(sigs.filter(Boolean))).join(', ');
    }
    // If no significance extracted from RCV, fall back to top-level clinical_significance
    if (!clinSig && annotation.clinvar) {
        let top = annotation.clinvar.clinical_significance;
        if (top) {
            // clinical_significance may be an object with description, a string, or an array
            if (typeof top === 'object' && top.description) {
                top = top.description;
            }
            if (Array.isArray(top)) {
                const unique = Array.from(new Set(top.filter(Boolean)));
                if (unique.length > 0) clinSig = unique.join(', ');
            } else {
                clinSig = String(top);
            }
        }
    }
    if (clinSig) rows.push({ name: 'Clinical Significance', value: clinSig });
    // PolyPhen prediction
    const polyphenPred = annotation.dbnsfp?.polyphen2?.hdiv?.pred;
    if (polyphenPred) rows.push({ name: 'PolyPhen2 Prediction', value: polyphenPred });
    // SIFT prediction
    const siftPred = annotation.dbnsfp?.sift?.pred;
    if (siftPred) rows.push({ name: 'SIFT Prediction', value: siftPred });
    // CADD phred score
    const caddPhred = annotation.cadd?.phred;
    if (caddPhred !== undefined) rows.push({ name: 'CADD Phred', value: caddPhred });
    return rows;
}

// Build a structured details view from the annotation object. Returns an array of
// category objects, each with a title and an array of {name, value} pairs.
function buildDetailsData(annotation, rawInput, gVariant) {
    const details = [];
    // Basic genomic details
    const basic = {};
    // Chromosome
    basic['Chromosome'] = annotation.chrom || annotation.cadd?.chrom || annotation.dbsnp?.chrom || '';
    // Genomic positions (hg19, hg38)
    if (annotation.hg19 || annotation.dbsnp?.hg19 || annotation.cadd?.pos) {
        if (annotation.hg19) {
            basic['hg19 Position'] = `${annotation.hg19.start}${annotation.hg19.end && annotation.hg19.end !== annotation.hg19.start ? '-' + annotation.hg19.end : ''}`;
        } else if (annotation.dbsnp?.hg19) {
            const start = annotation.dbsnp.hg19.start;
            const end = annotation.dbsnp.hg19.end;
            basic['hg19 Position'] = `${start}${end && end !== start ? '-' + end : ''}`;
        }
        if (annotation.hg38) {
            basic['hg38 Position'] = `${annotation.hg38.start}${annotation.hg38.end && annotation.hg38.end !== annotation.hg38.start ? '-' + annotation.hg38.end : ''}`;
        } else if (annotation.dbsnp?.hg38) {
            const start = annotation.dbsnp.hg38.start;
            const end = annotation.dbsnp.hg38.end;
            basic['hg38 Position'] = `${start}${end && end !== start ? '-' + end : ''}`;
        }
    }
    // Exon information
    if (annotation.cadd?.exon) basic['Exon'] = annotation.cadd.exon;
    // HGVS cDNA/protein from dbNSFP if available
    if (annotation.dbnsfp?.hgvsc) {
        basic['HGVS cDNA'] = Array.isArray(annotation.dbnsfp.hgvsc) ? annotation.dbnsfp.hgvsc.join(', ') : annotation.dbnsfp.hgvsc;
    }
    if (annotation.dbnsfp?.hgvsp) {
        basic['HGVS protein'] = Array.isArray(annotation.dbnsfp.hgvsp) ? annotation.dbnsfp.hgvsp.join(', ') : annotation.dbnsfp.hgvsp;
    }
    if (Object.keys(basic).length > 0) {
        details.push({ title: 'Genomic Details', items: basic });
    }

    // ClinVar section
    if (annotation.clinvar) {
        const clin = {};
        // Overall significance if present
        if (annotation.clinvar.clinical_significance?.description) {
            clin['Clinical Significance'] = annotation.clinvar.clinical_significance.description;
        }
        // Variation ID
        if (annotation.clinvar.variant_id) {
            clin['ClinVar Variant ID'] = annotation.clinvar.variant_id;
        }
        // RCV accessions and conditions
        const conditions = new Set();
        const sigs = new Set();
        const reviewStatuses = new Set();
        if (annotation.clinvar.rcv) {
            const rcv = annotation.clinvar.rcv;
            const entries = Array.isArray(rcv) ? rcv : [rcv];
            entries.forEach(rc => {
                // Condition names
                const conds = rc.conditions?.name || rc.conditions?.name_normalized || rc.conditions?.trait;
                if (conds) {
                    if (Array.isArray(conds)) conds.forEach(c => conditions.add(c));
                    else conditions.add(conds);
                }
                // Clinical significance
                const desc = rc.clinical_significance?.description || rc.clinical_significance;
                if (desc) {
                    if (Array.isArray(desc)) desc.forEach(d => sigs.add(d)); else sigs.add(desc);
                }
                // Review status
                if (rc.review_status) reviewStatuses.add(rc.review_status);
            });
        }
        if (conditions.size > 0) clin['Condition(s)'] = Array.from(conditions).join(', ');
        if (sigs.size > 0) clin['Significance'] = Array.from(sigs).join(', ');
        if (reviewStatuses.size > 0) clin['Review Status'] = Array.from(reviewStatuses).join(', ');
        if (Object.keys(clin).length > 0) {
            details.push({ title: 'ClinVar', items: clin });
        }
    }

    // CADD details
    if (annotation.cadd) {
        const cadd = {};
        if (annotation.cadd.phred !== undefined) cadd['CADD Phred'] = annotation.cadd.phred;
        if (annotation.cadd.rawscore !== undefined) cadd['CADD Raw Score'] = annotation.cadd.rawscore;
        if (annotation.cadd.consequence) cadd['Consequence'] = annotation.cadd.consequence;
        if (annotation.cadd.consdetail) cadd['Detailed Consequence'] = annotation.cadd.consdetail;
        if (Object.keys(cadd).length > 0) details.push({ title: 'CADD', items: cadd });
    }

    // Functional prediction scores from dbNSFP
    if (annotation.dbnsfp) {
        const preds = {};
        const addPred = (label, obj, predKey = 'pred', scoreKey = 'score', catKey = 'cat', valKey = 'val') => {
            if (!obj) return;
            let pred = '';
            // Some fields like polyphen2 have hdiv and hvar sub-objects
            if (obj.hdiv) {
                const hdiv = obj.hdiv;
                pred += `HDIV: ${hdiv.pred || hdiv.cat || ''}`;
                if (hdiv.score || hdiv.val) pred += ` (${hdiv.score !== undefined ? hdiv.score : hdiv.val})`;
                const hvar = obj.hvar;
                if (hvar) {
                    pred += `; HVAR: ${hvar.pred || hvar.cat || ''}`;
                    if (hvar.score || hvar.val) pred += ` (${hvar.score !== undefined ? hvar.score : hvar.val})`;
                }
            } else {
                const p = obj[predKey] || obj[catKey];
                const s = obj[scoreKey] !== undefined ? obj[scoreKey] : obj[valKey];
                if (p !== undefined) {
                    pred += String(p);
                    if (s !== undefined) pred += ` (${s})`;
                }
            }
            if (pred) preds[label] = pred;
        };
        addPred('SIFT', annotation.dbnsfp.sift);
        addPred('PolyPhen2', annotation.dbnsfp.polyphen2);
        addPred('ClinPred', annotation.dbnsfp.clinpred);
        addPred('MutationTaster', annotation.dbnsfp.mutationtaster);
        addPred('MutationAssessor', annotation.dbnsfp.mutationassessor);
        addPred('FATHMM', annotation.dbnsfp.fathmm);
        addPred('LRT', annotation.dbnsfp.lrt);
        addPred('PROVEAN', annotation.dbnsfp.provean);
        addPred('PrimateAI', annotation.dbnsfp.primateai);
        addPred('DANN', annotation.dbnsfp.dann);
        addPred('Deleteriousness (BayesDel)', annotation.dbnsfp.bayesdel?.no_af);

        // Attempt to add AlphaMissense prediction if available. Some MyVariant.info annotations store
        // AlphaMissense under different keys. Try to detect common patterns.
        if (annotation.dbnsfp) {
            // Some fields may appear as alphamissense_pred (single prediction) or nested objects
            if (annotation.dbnsfp.alphamissense_pred) {
                const val = annotation.dbnsfp.alphamissense_pred;
                preds['AlphaMissense'] = Array.isArray(val) ? Array.from(new Set(val)).join(',') : String(val);
            } else if (annotation.dbnsfp.alphamissense) {
                const am = annotation.dbnsfp.alphamissense;
                let val;
                if (typeof am === 'string') {
                    val = am;
                } else {
                    // Try to read .pred or .cat
                    val = am.pred || am.cat || am.value;
                }
                if (val !== undefined) {
                    if (Array.isArray(val)) {
                        val = val.filter(Boolean);
                        val = Array.from(new Set(val.map(String)));
                        preds['AlphaMissense'] = val.join(',');
                    } else {
                        // If comma-separated string, deduplicate
                        const parts = String(val).split(/[,;\s]+/).filter(Boolean);
                        preds['AlphaMissense'] = Array.from(new Set(parts)).join(',');
                    }
                }
            } else if (annotation.dbnsfp['alpha_missense']) {
                const am = annotation.dbnsfp['alpha_missense'];
                let val = am.pred || am.cat || am.value || am;
                if (val !== undefined) {
                    if (Array.isArray(val)) {
                        val = val.filter(Boolean);
                        val = Array.from(new Set(val.map(String)));
                        preds['AlphaMissense'] = val.join(',');
                    } else {
                        const parts = String(val).split(/[,;\s]+/).filter(Boolean);
                        preds['AlphaMissense'] = Array.from(new Set(parts)).join(',');
                    }
                }
            }
        }
        if (Object.keys(preds).length > 0) details.push({ title: 'Functional Predictions', items: preds });
    }

    // Allele frequencies from dbSNP (exac/gnomad)
    if (annotation.dbsnp?.alleles) {
        const freqs = {};
        annotation.dbsnp.alleles.forEach(a => {
            const allele = a.allele;
            if (a.freq) {
                const parts = [];
                if (a.freq.exac !== undefined) parts.push(`ExAC: ${a.freq.exac}`);
                if (a.freq.gnomad_exomes !== undefined) parts.push(`gnomAD exomes: ${a.freq.gnomad_exomes}`);
                if (a.freq.gnomad_genomes !== undefined) parts.push(`gnomAD genomes: ${a.freq.gnomad_genomes}`);
                if (parts.length > 0) freqs[`Allele ${allele}`] = parts.join('; ');
            }
        });
        if (Object.keys(freqs).length > 0) details.push({ title: 'Allele Frequencies', items: freqs });
    }

    // COSMIC mutation
    if (annotation.cosmic) {
        const cosmic = {};
        if (annotation.cosmic.cosmic_id) cosmic['COSMIC ID'] = annotation.cosmic.cosmic_id;
        if (annotation.cosmic.tumor_site) cosmic['Tumor Site'] = annotation.cosmic.tumor_site;
        if (annotation.cosmic.mut_freq !== undefined) cosmic['Mutation Frequency'] = annotation.cosmic.mut_freq;
        if (Object.keys(cosmic).length > 0) details.push({ title: 'COSMIC', items: cosmic });
    }

    // UniProt entries
    if (annotation.dbnsfp?.uniprot) {
        const upEntries = annotation.dbnsfp.uniprot;
        const uniprot = {};
        if (Array.isArray(upEntries)) {
            const accs = upEntries.map(u => u.acc || u.entry).filter(Boolean);
            if (accs.length > 0) uniprot['UniProt'] = accs.join(', ');
        }
        if (Object.keys(uniprot).length > 0) details.push({ title: 'UniProt', items: uniprot });
    }

    // gnomAD frequencies and link
    {
        const gnomad = {};
        // Helper to convert values that may be nested objects into a readable string
        const formatGnomadValue = (val) => {
            if (val === null || val === undefined) return val;
            // If the value is an object (e.g. allele-specific values), join key:value pairs
            if (typeof val === 'object') {
                const parts = [];
                for (const [allele, v] of Object.entries(val)) {
                    parts.push(`${allele}: ${v}`);
                }
                return parts.join(', ');
            }
            return val;
        };
        // Helper to add allele frequency and counts
        const addGnomadFields = (obj, prefix) => {
            if (!obj) return;
            if (obj.af !== undefined) gnomad[`${prefix} AF`] = formatGnomadValue(obj.af);
            if (obj.ac !== undefined) gnomad[`${prefix} AC`] = formatGnomadValue(obj.ac);
            if (obj.an !== undefined) gnomad[`${prefix} AN`] = formatGnomadValue(obj.an);
        };
        // Common gnomAD field names seen in MyVariant
        addGnomadFields(annotation.gnomad_exome, 'Exome');
        addGnomadFields(annotation.gnomad_genome, 'Genome');
        addGnomadFields(annotation.gnomad_exomes, 'Exome');
        addGnomadFields(annotation.gnomad_genomes, 'Genome');
        // Some annotations may use a top-level gnomad object
        addGnomadFields(annotation.gnomad, 'gnomAD');
        // Build link to gnomAD browser. Prefer raw token input so indels
        // still link to a useful region/variant view even without MyVariant hits.
        const url = buildGnomadVariantUrl(rawInput, gVariant, annotation._id, annotation.vcf);
        if (url) {
            gnomad['gnomAD Link'] = { html: `<a href="${url}" target="_blank" rel="noopener noreferrer">View on gnomAD</a>` };
        }
        if (Object.keys(gnomad).length > 0) details.push({ title: 'gnomAD', items: gnomad });
    }

    // CIViC annotation
    // MyVariant previously provided a `civic` object, but recent releases use the `cgi` array for CIViC evidence.
    if (annotation.civic || annotation.cgi) {
        // If a legacy `civic` object exists and is not an array, preserve the existing behaviour.
        if (annotation.civic && typeof annotation.civic === 'object' && !Array.isArray(annotation.civic)) {
            const civic = {};
            const c = annotation.civic;
            if (c.gene && c.gene.name) civic['Gene'] = c.gene.name;
            if (c.variant && c.variant.name) civic['Variant'] = c.variant.name;
            if (c.variant_id !== undefined) civic['CIViC Variant ID'] = c.variant_id;
            if (c.entrez_id !== undefined) civic['Entrez ID'] = c.entrez_id;
            if (c.evidence_items && Array.isArray(c.evidence_items)) civic['Evidence Items'] = c.evidence_items.length;
            if (c.evidence_level) civic['Evidence Level'] = c.evidence_level;
            if (c.source && c.source.url) civic['CIViC Source'] = { html: `<a href="${c.source.url}" target="_blank" rel="noopener noreferrer">${c.source.url}</a>` };
            if (Object.keys(civic).length > 0) details.push({ title: 'CIViC', items: civic });
        } else {
            // Use whichever array of evidence entries is available (annotation.cgi or annotation.civic if array)
            const entries = Array.isArray(annotation.cgi) ? annotation.cgi : (Array.isArray(annotation.civic) ? annotation.civic : []);
            if (entries.length > 0) {
                const geneNames = new Set();
                const proteinChanges = new Set();
                const evidenceLevels = new Set();
                const drugs = new Set();
                entries.forEach(item => {
                    if (item.gene) geneNames.add(item.gene);
                    if (item.protein_change) proteinChanges.add(item.protein_change);
                    if (item.evidence_level) evidenceLevels.add(item.evidence_level);
                    if (item.drug) drugs.add(item.drug);
                });
                const civic = {};
                if (geneNames.size > 0) civic['Gene(s)'] = Array.from(geneNames).join(', ');
                if (proteinChanges.size > 0) civic['Protein Change(s)'] = Array.from(proteinChanges).join(', ');
                if (drugs.size > 0) civic['Drug(s)'] = Array.from(drugs).join(', ');
                if (evidenceLevels.size > 0) civic['Evidence Level(s)'] = Array.from(evidenceLevels).join(', ');
                civic['Evidence Items'] = entries.length;
                if (geneNames.size > 0) {
                    const encodedGene = encodeURIComponent(Array.from(geneNames)[0]);
                    const civicUrl = `https://civicdb.org/genes/${encodedGene}`;
                    civic['CIViC Gene Page'] = { html: `<a href="${civicUrl}" target="_blank" rel="noopener noreferrer">View gene in CIViC</a>` };
                }
                if (Object.keys(civic).length > 0) details.push({ title: 'CIViC', items: civic });
            }
        }
    }
    return details;
}

document.addEventListener('DOMContentLoaded', () => {
    // Holds a triple‑coded protein change parsed from the user's query (e.g. VAL600GLU).
    // This will be used to help select the canonical variant and to build search queries.
    let targetProtGlobal = null;
    // Holds the cDNA change parsed from the user's query (e.g. c.178_186del). Used to prioritise
    // candidate genomic variants returned by the variant recoder.
    let targetCdnaGlobal = null;
    // Holds a list of transcript annotations (transcript ID, cDNA, protein) fetched via the Ensembl variant recoder.
    // This will be populated when the user submits a variant and used to build the transcripts list in the variant card.
    let transcriptsFromRecoder = [];

    // Apply a stable thematic class to each card so CSS can render distinct colors per section.
    const applyCardTheme = (cardEl, cardTitle) => {
        if (!cardEl || !cardTitle) return;
        const key = String(cardTitle).toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
        if (key) cardEl.setAttribute('data-card', key);
    };

    /**
     * Fetch transcripts for a given variant using the Ensembl variant recoder. Returns an array of
     * objects with transcript, cDNA and protein fields. If no transcripts can be obtained, returns an empty array.
     * @param {string} q Variant query (e.g. protein HGVS, cDNA or genomic variant).
     */
    async function getTranscriptsList(q) {
        try {
            const recResults = await fetchVariantRecoder(q);
            // The Ensembl variant_recoder response is an array with a single object.
            // Historically the transcripts were stored under recResults[0].A, but newer
            // versions use lettered keys (A, B, C, ... or even G) depending on the
            // mapping category. Each of these objects can contain hgvsc/hgvsp arrays.
            if (Array.isArray(recResults) && recResults.length > 0) {
                const recEntry = recResults[0];
                /*
                 * Ensembl's variant_recoder may return multiple sub-objects keyed by different alternate
                 * sequences (e.g. duplication or insertion sequences) when a query maps to more than one
                 * genomic event. In these cases the first object is not guaranteed to correspond to the
                 * user's intended variant. To improve accuracy, attempt to select the sub-object that
                 * contains a transcript whose cDNA exactly matches the user's query (for c. inputs) or,
                 * failing that, that contains a protein change matching a p. input. If neither match,
                 * fall back to the legacy logic of picking the first object with hgvsc/hgvsp data.
                 */
                let objWithTranscripts = null;
                // Extract the variant part from the query (e.g. "c.181_189dup" from "MSH3:c.181_189dup").
                let queryVariantPart = null;
                const qParts = String(q).split(':');
                if (qParts.length > 1) {
                    // Use everything after the first colon as the variant part
                    queryVariantPart = qParts.slice(1).join(':').trim().toLowerCase();
                }
                // If this is a genomic variant (g.) with an insertion (ins) or duplication, attempt to
                // extract the inserted sequence and match it against the keys of the recoder response. For
                // example, "chr5:g.79950735_79950736insGCAGCGCCC" should match the "GCAGCGCCC" key.
                if (!objWithTranscripts && queryVariantPart && /^g\./i.test(queryVariantPart)) {
                    // Attempt to capture the inserted sequence following "ins".
                    const mIns = queryVariantPart.match(/ins([A-Za-z]+)/i);
                    if (mIns) {
                        const altSeq = mIns[1].toLowerCase();
                        for (const [key, subVal] of Object.entries(recEntry)) {
                            if (typeof key === 'string' && key.toLowerCase() === altSeq && subVal && typeof subVal === 'object') {
                                objWithTranscripts = subVal;
                                break;
                            }
                        }
                    }
                }
                // First pass: look for a sub-object whose hgvsc or hgvsp array contains the query variant part.
                if (!objWithTranscripts && queryVariantPart) {
                    for (const subVal of Object.values(recEntry)) {
                        if (!subVal || typeof subVal !== 'object') continue;
                        const scs = Array.isArray(subVal.hgvsc) ? subVal.hgvsc : (subVal.hgvsc ? [subVal.hgvsc] : []);
                        const sps = Array.isArray(subVal.hgvsp) ? subVal.hgvsp : (subVal.hgvsp ? [subVal.hgvsp] : []);
                        // Check if any cDNA matches the variant part exactly (case-insensitive)
                        let matchFound = false;
                        for (const sc of scs) {
                            const cPart = String(sc).split(':').slice(1).join(':').trim().toLowerCase();
                            if (cPart === queryVariantPart) {
                                matchFound = true;
                                break;
                            }
                        }
                        // Check for protein match if no cDNA match yet and query looks like a protein variant
                        if (!matchFound && /^p\./i.test(queryVariantPart)) {
                            for (const sp of sps) {
                                const pPart = String(sp).split(':').slice(1).join(':').trim().toLowerCase();
                                if (pPart === queryVariantPart) {
                                    matchFound = true;
                                    break;
                                }
                            }
                        }
                        if (matchFound) {
                            objWithTranscripts = subVal;
                            break;
                        }
                    }
                }
                // If no direct match found, prefer the traditional "A" key if it contains transcript data.
                if (!objWithTranscripts && recEntry.A && typeof recEntry.A === 'object' && (recEntry.A.hgvsc || recEntry.A.hgvsp)) {
                    objWithTranscripts = recEntry.A;
                }
                // Otherwise pick the first sub-object that has hgvsc/hgvsp arrays.
                if (!objWithTranscripts) {
                    for (const subVal of Object.values(recEntry)) {
                        if (subVal && typeof subVal === 'object' && (subVal.hgvsc || subVal.hgvsp)) {
                            objWithTranscripts = subVal;
                            break;
                        }
                    }
                }
                if (objWithTranscripts) {
                    const hgvscs = Array.isArray(objWithTranscripts.hgvsc) ? objWithTranscripts.hgvsc : (objWithTranscripts.hgvsc ? [objWithTranscripts.hgvsc] : []);
                    const hgvsp = Array.isArray(objWithTranscripts.hgvsp) ? objWithTranscripts.hgvsp : (objWithTranscripts.hgvsp ? [objWithTranscripts.hgvsp] : []);
                    const len = Math.max(hgvscs.length, hgvsp.length);
                    const list = [];
                    // Limit the number of transcripts we process to avoid freezing the UI on variants with
                    // extremely large transcript sets (e.g. TP53). We keep the first 200 entries.
                    const maxEntries = 200;
                    for (let i = 0; i < Math.min(len, maxEntries); i++) {
                        const sc = hgvscs[i];
                        const sp = hgvsp[i];
                        if (sc) {
                            const parts = String(sc).split(':');
                            const transcriptId = parts[0];
                            // Join remaining parts for cDNA and trim whitespace. Some variant_recoder
                            // responses include a space after the colon (e.g. "NM_003620.4: c.1518del"), so
                            // trimming ensures consistent matching.
                            const cpart = parts.slice(1).join(':').trim();
                            let prot = '';
                            if (sp) {
                                const pparts = String(sp).split(':');
                                prot = pparts.slice(1).join(':').trim();
                            }
                            list.push({ transcript: transcriptId, cDNA: cpart, protein: prot });
                        }
                    }
                    return list;
                }
            }
        } catch (err) {
            // If the recoder request fails, fall through and return an empty array.
        }
        return [];
    }
    const form = document.getElementById('variantForm');
    const input = document.getElementById('variantInput');
    const tumorTypeInput = document.getElementById('tumorTypeInput');
    const statusEl = document.getElementById('status');
    const resultSection = document.getElementById('resultSection');
    const summaryTable = document.getElementById('summaryTable');
    const rawOutput = document.getElementById('rawOutput');
    // Pre element to display the Ensembl variant_recoder JSON response
    const ensemblOutput = document.getElementById('ensemblOutput');
    // Keep a shareable URL in sync with the current search value so external sites can link
    // directly to this page with a query string like `/?variant=<any-supported-input>`.
    // Parse query params from window.location.search while preserving literal "+" characters.
    // This helps with intronic HGVS such as c.76+1G>A when a referring site provides an
    // unescaped plus sign in the URL.
    const getVariantFromUrl = () => {
        const search = String(window.location.search || '');
        const m = search.match(/[?&]variant=([^&]*)/);
        if (!m) return '';
        try {
            const rawValue = m[1];
            const decodedLiteralPlus = decodeURIComponent(rawValue);
            const decodedFormStyle = decodeURIComponent(rawValue.replace(/\+/g, '%20'));
            // Heuristic:
            // - For generic free-text inputs, treat "+" as space (standard query-string behavior).
            // - Preserve literal "+" for intronic HGVS patterns (e.g. c.76+1G>A) and for
            //   explicitly encoded plus signs (%2B), which decode identically in both modes.
            const looksLikeIntronicHgvs = /(?:^|[\s:])c\.[^\s]*[+-]\d+/i.test(decodedLiteralPlus);
            if (looksLikeIntronicHgvs) return decodedLiteralPlus;
            return decodedFormStyle;
        } catch {
            return m[1];
        }
    };

    const syncVariantInUrl = (variantValue) => {
        try {
            const url = new URL(window.location.href);
            const value = String(variantValue || '').trim();
            if (value) {
                url.searchParams.set('variant', value);
            } else {
                url.searchParams.delete('variant');
            }
            const search = url.searchParams.toString();
            const nextUrl = `${url.pathname}${search ? `?${search}` : ''}${url.hash || ''}`;
            window.history.replaceState({}, '', nextUrl);
        } catch {
            // Ignore URL update errors and continue normal search flow.
        }
    };

    const getTumorTypeFromUrl = () => {
        const search = String(window.location.search || '');
        const m = search.match(/[?&]tumorType=([^&]*)/);
        if (!m) return '';
        try {
            return decodeURIComponent(m[1].replace(/\+/g, '%20'));
        } catch {
            return m[1];
        }
    };

    const syncTumorTypeInUrl = (tumorTypeValue) => {
        try {
            const url = new URL(window.location.href);
            const value = String(tumorTypeValue || '').trim();
            if (value) {
                url.searchParams.set('tumorType', value);
            } else {
                url.searchParams.delete('tumorType');
            }
            const search = url.searchParams.toString();
            const nextUrl = `${url.pathname}${search ? `?${search}` : ''}${url.hash || ''}`;
            window.history.replaceState({}, '', nextUrl);
        } catch {
            // Ignore URL update errors and continue normal search flow.
        }
    };

    form.addEventListener('submit', async (e) => {
        e.preventDefault();
        // Reset any previous gene hint.  Without resetting, a prior search that
        // included a gene symbol could incorrectly influence the next query.
        geneHintGlobal = null;
        alterationTypeGlobal = null;
        isGeneOnlyMode = false;

        // Normalise user input to improve variant parsing. Accept inputs like
        // "BRAF V600E", "braf:p.v600e", "BRAF:V600E" etc. by converting them
        // to standard HGVS-like strings. This heuristic uppercases gene symbols
        // and prepends "p." to protein variants when missing.
        const rawInput = input.value.trim();
        const tumorType = tumorTypeInput ? tumorTypeInput.value.trim() : '';
        // Update the URL with the raw submitted value so inbound/outbound links can
        // preserve flexible user-entered formats (HGVS g./c./p., space-separated tokens, etc.).
        syncVariantInUrl(rawInput);
        syncTumorTypeInUrl(tumorType);
        let query = rawInput;

        // Detect gene-only or gene+descriptor input before normalisation.
        {
            const ALTERATION_DESCRIPTORS = /^(mutation|fusion|rearrangement|loss|amplification|gain|deletion|overexpression|splice|truncation|frameshift|alteration)s?$/i;
            const toks = rawInput.split(/[\s:]+/).filter(Boolean);
            const looksLikeGene = (s) => /^[A-Za-z][A-Za-z0-9-]{1,9}$/.test(s) && !/^(?:chr|NM_|NP_|rs\d)/i.test(s);
            if (toks.length === 1 && looksLikeGene(toks[0])) {
                geneHintGlobal = toks[0].toUpperCase();
                alterationTypeGlobal = null;
                isGeneOnlyMode = true;
            } else if (toks.length === 2 && looksLikeGene(toks[0]) && ALTERATION_DESCRIPTORS.test(toks[1])) {
                geneHintGlobal = toks[0].toUpperCase();
                alterationTypeGlobal = toks[1].toLowerCase().replace(/s$/, '');
                isGeneOnlyMode = true;
            }
        }

        const normalizeVariantInput = (raw) => {
            let s = raw.trim();
            // If the input contains a '>' (looks like genomic substitution), attempt to parse flexible genomic formats.
            if (s.includes('>')) {
                try {
            // Split around '>' to get left (chr/pos/ref) and right (alt)
            const partsGT = s.split('>');
            if (partsGT.length === 2) {
                let left = partsGT[0].trim();
                let right = partsGT[1].trim();
                // Clean alt by removing whitespace and punctuation. Uppercase nucleotides for consistency.
                const alt = right.replace(/\s+/g, '').replace(/[^A-Za-z-]/g, '').toUpperCase();
                // Normalise the left portion by removing commas (thousands separators) and trimming whitespace.
                const leftNorm = left.replace(/,/g, '').trim();
                /*
                 * Attempt to parse the chromosome, position and reference allele from the left side. Support
                 * flexible formats including:
                 *   chr7:g.140453136A
                 *   chr7:140453136A
                 *   7:140453136A
                 *   chr7 g.140453136A
                 *   7 140453136A
                 * The regular expression captures:
                 *   - An optional "chr" prefix
                 *   - The chromosome (digits or X/Y/MT)
                 *   - Optional separator consisting of whitespace or a colon
                 *   - An optional "g" followed by an optional dot
                 *   - The numeric coordinate (group 2)
                 *   - The reference allele (letters or hyphen) (group 3)
                 */
                const m2 = leftNorm.match(/^(?:chr)?([0-9XYMT]+)[\s:]*g?\.?([0-9]+)([A-Za-z-]+)$/i);
                if (m2) {
                    let chrom = m2[1].toUpperCase();
                    const pos = m2[2];
                    const ref = m2[3].toUpperCase();
                    return `chr${chrom}:g.${pos}${ref}>${alt}`;
                }
            }
                } catch {
                    // fall through to other normalizations
                }
            }
            // Handle cases where genomic variant is specified as separate tokens: "chr7 140453136 A T" or "7 140453136 A T".
            {
                const toks = s.split(/\s+/).filter(Boolean);
                // Support both four-token inputs (chr pos ref alt) and five-token inputs where a gene symbol appears between position and ref.
                if (toks.length === 4 || toks.length === 5) {
                    // Create a copy of tokens and remove a potential gene symbol at index 2 when length is 5.
                    let tokens = toks.slice();
                    if (tokens.length === 5) {
                        const token3 = tokens[2];
                        // For 5-token genomic input, always ignore token #3 as an optional gene/label
                        // so values like TP53, ENSG IDs or other labels do not block normalisation.
                        // Preserve a simple uppercase token as a gene hint for fallback displays.
                        if (/^[A-Za-z][A-Za-z0-9-]*$/.test(token3)) {
                            geneHintGlobal = token3.toUpperCase();
                        }
                        tokens.splice(2, 1);
                    }
                    if (tokens.length === 4) {
                        const [chrTok, posTok, refTok, altTok] = tokens;
                        // Ensure chromosome token, position, ref and alt are valid DNA-like sequences
                        if (/^[0-9XYMT]+$/.test(chrTok.replace(/^chr/i, '')) && /^\d+$/.test(posTok.replace(/,/g, '')) && /^[A-Za-z-]+$/.test(refTok) && /^[A-Za-z-]+$/.test(altTok)) {
                            let chrom = chrTok.replace(/^chr/i, '').toUpperCase();
                            // Strip commas from the coordinate before parsing (e.g. "140,453,136" -> "140453136").
                            const pos = parseInt(posTok.replace(/,/g, ''), 10);
                            const ref = refTok.toUpperCase();
                            const alt = altTok.toUpperCase();
                            // Handle length differences between ref and alt: deletions, insertions and complex delins.
                            if (ref.length !== alt.length) {
                                // Deletion or delins: alt shorter than ref indicates a deletion or a complex substitution
                                if (alt.length < ref.length) {
                                    /*
                                     * Trim common prefix and suffix between ref and alt to obtain the
                                     * minimal representation for delins variants.  For example, the
                                     * input "7 140453122 TCCATCGAGATTTCA TCT" should normalise to
                                     * chr7:g.140453124_140453136delinsT, because the shared "TC" prefix
                                     * can be removed and the genomic start coordinate shifted.  If the
                                     * resulting alt sequence is empty after trimming, this represents a
                                     * pure deletion.  See documentation on variant normalisation.
                                     */
                                    let refTrim = ref;
                                    let altTrim = alt;
                                    let startPos = pos;
                                    let endPos = pos + ref.length - 1;
                                    // Trim from the left while both strings share the same leading base
                                    while (altTrim.length > 0 && refTrim.length > 0 && altTrim[0] === refTrim[0]) {
                                        altTrim = altTrim.slice(1);
                                        refTrim = refTrim.slice(1);
                                        startPos += 1;
                                    }
                                    // Trim from the right while both strings share the same trailing base
                                    while (altTrim.length > 0 && refTrim.length > 0 && altTrim[altTrim.length - 1] === refTrim[refTrim.length - 1]) {
                                        altTrim = altTrim.slice(0, -1);
                                        refTrim = refTrim.slice(0, -1);
                                        endPos -= 1;
                                    }
                                    // After trimming, if altTrim is empty, it's a pure deletion.
                                    if (altTrim.length === 0) {
                                        if (startPos === endPos) {
                                            return `chr${chrom}:g.${startPos}del`;
                                        }
                                        return `chr${chrom}:g.${startPos}_${endPos}del`;
                                    }
                                    // Otherwise it's a delins with the trimmed alt sequence.
                                    return `chr${chrom}:g.${startPos}_${endPos}delins${altTrim}`;
                                }
                                // Insertion or delins: alt longer than ref
                                else if (alt.length > ref.length) {
                                    // If ref is empty, simple insertion between pos and pos+1
                                    if (ref === '') {
                                        return `chr${chrom}:g.${pos}_${pos + 1}ins${alt}`;
                                    }
                                    // If ref is a prefix of alt, treat as insertion of the suffix
                                    if (alt.startsWith(ref)) {
                                        const insSeq = alt.slice(ref.length);
                                        const insPosStart = pos + ref.length - 1;
                                        const insPosEnd = pos + ref.length;
                                        return `chr${chrom}:g.${insPosStart}_${insPosEnd}ins${insSeq}`;
                                    }
                                    // Otherwise treat as delins substitution
                                    const delStart2 = pos;
                                    const delEnd2 = pos + ref.length - 1;
                                    return `chr${chrom}:g.${delStart2}_${delEnd2}delins${alt}`;
                                }
                            }
                            // Lengths equal: represent SNVs as ref>alt, but use delins for MNVs.
                            if (ref.length === 1) {
                                return `chr${chrom}:g.${pos}${ref}>${alt}`;
                            }
                            const endPos = pos + ref.length - 1;
                            return `chr${chrom}:g.${pos}_${endPos}delins${alt}`;
                        }
                    }
                }
            }
            // If no special cases above matched, continue with other normalisations
            const sOrig = s;
            // If already looks like genomic HGVS (chrX:g.) or contains canonical prefixes, return as-is after trimming.
            // If it looks like a genomic variant (contains ':g.') handle separately: just normalise chr prefix
            if (isGenomicVariant(s)) {
                // Make 'chr' prefix lowercase and leave the rest untouched
                const parts = s.split(':');
                const gene = parts[0];
                const rest = parts.slice(1).join(':');
                // Ensure 'chr' is lowercase and chromosome symbol is unchanged
                return gene.replace(/^chr/i, 'chr') + (rest ? ':' + rest : '');
            }
            // If variant contains a colon with p. or c., still normalise variant part
            if (/:[pc]\./i.test(s)) {
                const [genePart, variantPartRaw] = s.split(/:/);
                let gene = genePart.toUpperCase();
                let variantPart = variantPartRaw.trim();
                // Normalize the variant prefix to lowercase if it's p. or c. (e.g., "P." -> "p." or "C." -> "c.")
                if (/^[pPcC]\./.test(variantPart)) {
                    variantPart = variantPart.charAt(0).toLowerCase() + variantPart.slice(1);
                }
                // Normalize case for three-letter amino acid codes if p. variant
                if (/^p\./i.test(variantPart)) {
                    // remove leading 'p.' for processing then add back later
                    let vp = variantPart.replace(/^p\./i, '');
                    // If triple-letter format (Val600Glu etc.)
                    const tripleMatch = vp.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3})$/);
                    if (tripleMatch) {
                        const cap = (str) => str.charAt(0).toUpperCase() + str.slice(1).toLowerCase();
                        vp = `${cap(tripleMatch[1])}${tripleMatch[2]}${cap(tripleMatch[3])}`;
                    } else {
                        // If single-letter format (V600E)
                        const singleMatch = vp.match(/^([A-Za-z])(\d+)([A-Za-z])$/);
                        if (singleMatch) {
                            const aaMap = {
                                A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly',
                                H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser',
                                T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val'
                            };
                            const ref = aaMap[singleMatch[1].toUpperCase()];
                            const pos = singleMatch[2];
                            const alt = aaMap[singleMatch[3].toUpperCase()];
                            if (ref && alt) {
                                vp = `${ref}${pos}${alt}`;
                            }
                        }
                    }
                    variantPart = 'p.' + vp;
                }
                // For c. variant, just ensure lower-case c and uppercase gene part
                if (/^c\./i.test(variantPart)) {
                    variantPart = 'c.' + variantPart.slice(2);
                }
                return `${gene}:${variantPart}`;
            }
            // Split on whitespace or colon. This catches inputs like "BRAF V600E" or "BRAF:V600E".
            const parts = s.split(/[:\s]+/).filter(Boolean);
            if (parts.length >= 2) {
                let gene = parts[0].toUpperCase();
                let variantPart = parts[1].trim();
                // If the variant part begins with an explicit p. or c. prefix (case-insensitive),
                // normalise the prefix to lowercase and defer any case adjustments until later.
                if (/^[pPcC]\./.test(variantPart)) {
                    variantPart = variantPart.charAt(0).toLowerCase() + variantPart.slice(1);
                } else {
                    // For variants without an explicit prefix, capitalise the first letter
                    // and ensure the trailing amino-acid code (if present) is uppercase.
                    variantPart = variantPart
                        .replace(/^[a-z]/, (m) => m.toUpperCase())
                        .replace(/([0-9])([a-z])$/, (_, num, aa) => num + aa.toUpperCase());
                }
                // Prepend "p." if variant part does not already start with p. or c. or g. (case-insensitive)
                if (!/^p\./i.test(variantPart) && !/^c\./i.test(variantPart) && !/^g\./i.test(variantPart)) {
                    variantPart = 'p.' + variantPart;
                }
                // Normalize case for three-letter amino acid codes. If the variant is already
                // in a triple-coded format (e.g. p.val600glu), capitalise the first letter of
                // each amino acid and lowercase the remaining letters.
                {
                    const m = variantPart.match(/^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$/);
                    if (m) {
                        const ref3 = m[1];
                        const pos = m[2];
                        const alt3 = m[3];
                        const cap = (s) => s.charAt(0).toUpperCase() + s.slice(1).toLowerCase();
                        variantPart = `p.${cap(ref3)}${pos}${cap(alt3)}`;
                    }
                }
                // Convert single-letter amino acid codes to three-letter codes for better HGVS recognition.
                const aaMap = {
                    A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly',
                    H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser',
                    T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val'
                };
                const convertSingleToTriple = (v) => {
                    const m = v.match(/^p\.([A-Za-z])([0-9]+)([A-Za-z])$/);
                    if (m) {
                        const ref = aaMap[m[1].toUpperCase()];
                        const pos = m[2];
                        const alt = aaMap[m[3].toUpperCase()];
                        if (ref && alt) {
                            return `p.${ref}${pos}${alt}`;
                        }
                    }
                    return v;
                };
                variantPart = convertSingleToTriple(variantPart);
                return `${gene}:${variantPart}`;
            }
            return s;
        };
        query = normalizeVariantInput(query);
        // Log normalized query for debugging
        console.log('[DEBUG] Normalized query:', query);

        // ── Gene-only / gene+descriptor mode ────────────────────────────────
        // When the input is just a gene symbol or a gene + alteration descriptor
        // (e.g. "BRAF" or "BRAF amplification"), skip the variant annotation
        // pipeline entirely and render a focused set of gene-level cards instead.
        if (isGeneOnlyMode && geneHintGlobal) {
            const gene = geneHintGlobal;
            const altType = alterationTypeGlobal; // e.g. "amplification" or null

            statusEl.textContent = '';
            resultSection.classList.remove('hidden');

            const detailsContainer = document.getElementById('detailsContainer');
            if (detailsContainer) {
                detailsContainer.innerHTML = '';
                detailsContainer.classList.add('hidden');
            }
            summaryTable.innerHTML = '';
            summaryTable.classList.add('hidden');

            const cardsContainer = document.getElementById('cardsContainer');
            if (cardsContainer) cardsContainer.innerHTML = '';
            if (cardsContainer) cardsContainer.classList.remove('hidden');

            // Local aiReviewExtras for gene-only mode (populated by async card callbacks)
            const geneOnlyAiExtras = {};

            // ── Card: CIViC ────────────────────────────────────────────────
            {
                const civicCard = document.createElement('div');
                civicCard.className = 'card';
                const civicTitle = document.createElement('h3');
                civicTitle.textContent = 'CIViC';
                applyCardTheme(civicCard, 'CIViC');
                civicCard.appendChild(civicTitle);
                const content = document.createElement('div');
                content.className = 'card-content';

                // Gene link
                const linksDiv = document.createElement('div');
                linksDiv.style.marginBottom = '0.4rem';
                const encodedGene = encodeURIComponent(gene);
                const geneLinkEl = document.createElement('a');
                geneLinkEl.href = `https://civicdb.org/genes/${gene}/summary`;
                geneLinkEl.target = '_blank';
                geneLinkEl.rel = 'noopener noreferrer';
                geneLinkEl.textContent = `View ${gene} on CIViC`;
                linksDiv.appendChild(geneLinkEl);
                content.appendChild(linksDiv);

                if (altType) {
                    const noteEl = document.createElement('div');
                    noteEl.style.cssText = 'font-size:0.85rem;color:#6b7280;margin-bottom:0.4rem;';
                    noteEl.textContent = `Showing gene-level data for ${gene}. Alteration filter: ${altType}.`;
                    content.appendChild(noteEl);
                }

                // CIViC API section (async)
                const civicApiDiv = document.createElement('div');
                civicApiDiv.style.marginTop = '0.5rem';
                const civicApiSpinner = document.createElement('div');
                civicApiSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                civicApiSpinner.textContent = 'Loading CIViC API data…';
                civicApiDiv.appendChild(civicApiSpinner);
                content.appendChild(civicApiDiv);
                civicCard.appendChild(content);
                if (cardsContainer) cardsContainer.appendChild(civicCard);

                fetchCivicApiData(gene, '').then((civicApiData) => {
                    geneOnlyAiExtras.civic_api = civicApiData;
                    civicApiDiv.innerHTML = '';
                    if (!civicApiData) {
                        civicApiDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">CIViC API not accessible from browser — use the link above to search CIViC directly.</div>';
                        return;
                    }
                    const { gene: apiGene, matchedVariant, assertions } = civicApiData;
                    if (!apiGene) {
                        civicApiDiv.innerHTML = `<div style="font-size:0.82rem;color:#9ca3af;">Gene "${gene}" not found in CIViC.</div>`;
                        return;
                    }

                    // Upgrade gene link
                    if (apiGene.id) {
                        geneLinkEl.href = `https://civicdb.org/features/${apiGene.id}/summary`;
                    }

                    const geneLevelHeading = document.createElement('div');
                    geneLevelHeading.style.cssText = 'font-size:0.9rem;font-weight:600;margin:0.45rem 0 0.25rem;';
                    geneLevelHeading.textContent = 'Gene-level data from CIViC:';
                    civicApiDiv.appendChild(geneLevelHeading);

                    // AMP assertion level
                    if (assertions && assertions.length > 0) {
                        const topAssert = assertions[0];
                        const ampLevel = topAssert.ampLevel || topAssert.significance || '';
                        if (ampLevel) {
                            const ampEl = document.createElement('div');
                            ampEl.style.cssText = 'font-size:0.9rem;font-weight:600;margin-bottom:4px;';
                            ampEl.innerHTML = `<strong>AMP/ACMG tier (CIViC):</strong> ${ampLevel}`;
                            civicApiDiv.appendChild(ampEl);
                        }
                    }

                    // Gene description
                    if (apiGene.description) {
                        const descDet = document.createElement('details');
                        const descSum = document.createElement('summary');
                        descSum.style.fontSize = '0.85rem';
                        descSum.textContent = 'CIViC gene description';
                        descDet.appendChild(descSum);
                        const descText = document.createElement('div');
                        descText.style.cssText = 'font-size:0.82rem;padding:4px 0;line-height:1.45;color:#374151;';
                        const desc = String(apiGene.description);
                        descText.textContent = desc.length > 600 ? desc.slice(0, 600) + '…' : desc;
                        descDet.appendChild(descText);
                        civicApiDiv.appendChild(descDet);
                    }

                    // Variants list (filter by altType if set, but show all with highlight)
                    const allVariants = Array.isArray(apiGene.variants?.nodes) ? apiGene.variants.nodes : [];
                    if (allVariants.length > 0) {
                        const varDet = document.createElement('details');
                        const varSum = document.createElement('summary');
                        varSum.style.fontSize = '0.85rem';
                        varSum.textContent = `CIViC variants for ${apiGene.name} (${allVariants.length})`;
                        varDet.appendChild(varSum);
                        const varTable = document.createElement('table');
                        varTable.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.79rem;margin-top:4px;';
                        const varThead = varTable.createTHead();
                        const varHrow = varThead.insertRow();
                        ['Variant', 'Types'].forEach((h) => {
                            const th = document.createElement('th');
                            th.textContent = h;
                            th.style.cssText = 'text-align:left;padding:2px 5px;background:#f3f4f6;font-size:0.77rem;';
                            varHrow.appendChild(th);
                        });
                        const varTbody = varTable.createTBody();
                        allVariants.slice(0, 30).forEach((v) => {
                            const tr = varTbody.insertRow();
                            const nameMatch = altType && v.name && v.name.toLowerCase().includes(altType.toLowerCase());
                            if (nameMatch) tr.style.background = '#fffbeb';
                            const vTypes = Array.isArray(v.variantTypes?.nodes)
                                ? v.variantTypes.nodes.map((vt) => vt.name).filter(Boolean).join(', ')
                                : '';
                            [v.name || 'N/A', vTypes || 'N/A'].forEach((val, ci) => {
                                const td = tr.insertCell();
                                if (ci === 0 && v.id) {
                                    const a = document.createElement('a');
                                    a.href = `https://civicdb.org/variants/${v.id}/summary`;
                                    a.target = '_blank';
                                    a.rel = 'noopener noreferrer';
                                    a.textContent = val;
                                    td.appendChild(a);
                                } else {
                                    td.textContent = val;
                                }
                                td.style.cssText = 'padding:2px 5px;border-bottom:1px solid #f0f0f0;vertical-align:top;';
                            });
                        });
                        varDet.appendChild(varTable);
                        if (allVariants.length > 30) {
                            const moreEl = document.createElement('div');
                            moreEl.style.cssText = 'font-size:0.78rem;color:#666;padding:3px 5px;';
                            moreEl.textContent = `Showing 30 of ${allVariants.length} variants. Visit CIViC for the full list.`;
                            varDet.appendChild(moreEl);
                        }
                        civicApiDiv.appendChild(varDet);
                        if (altType) {
                            const filtNote = document.createElement('div');
                            filtNote.style.cssText = 'font-size:0.78rem;color:#b45309;margin-top:3px;';
                            filtNote.textContent = `Rows highlighted in yellow contain "${altType}" in the variant name.`;
                            civicApiDiv.appendChild(filtNote);
                        }
                    }

                    // Assertions table
                    if (assertions && assertions.length > 0) {
                        const assertDet = document.createElement('details');
                        const assertSum = document.createElement('summary');
                        assertSum.style.fontSize = '0.85rem';
                        assertSum.textContent = `CIViC assertions — AMP/ACMG (${assertions.length})`;
                        assertDet.appendChild(assertSum);
                        const aTable = document.createElement('table');
                        aTable.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.79rem;margin-top:4px;';
                        const aThead = aTable.createTHead();
                        const aHrow = aThead.insertRow();
                        ['AMP Level', 'Significance', 'Disease', 'Therapies'].forEach((h) => {
                            const th = document.createElement('th');
                            th.textContent = h;
                            th.style.cssText = 'text-align:left;padding:2px 5px;background:#f3f4f6;font-size:0.77rem;';
                            aHrow.appendChild(th);
                        });
                        const aTbody = aTable.createTBody();
                        assertions.slice(0, 8).forEach((a) => {
                            const tr = aTbody.insertRow();
                            const amp = a.ampLevel || 'N/A';
                            const sig = a.clinicalSignificance || a.significance || 'N/A';
                            const disease = a.disease?.name || 'N/A';
                            const therapies = Array.isArray(a.therapies?.nodes)
                                ? a.therapies.nodes.map((t) => t.name).join(', ')
                                : (Array.isArray(a.therapies) ? a.therapies.join(', ') : 'N/A');
                            [amp, sig, disease, therapies].forEach((val) => {
                                const td = tr.insertCell();
                                td.textContent = val;
                                td.style.cssText = 'padding:2px 5px;border-bottom:1px solid #f0f0f0;vertical-align:top;';
                            });
                        });
                        assertDet.appendChild(aTable);
                        civicApiDiv.appendChild(assertDet);
                    } else {
                        const noAssert = document.createElement('div');
                        noAssert.style.cssText = 'font-size:0.82rem;color:#6b7280;margin-top:2px;';
                        noAssert.textContent = 'No accepted CIViC assertions for this gene.';
                        civicApiDiv.appendChild(noAssert);
                    }

                    const variantCount = allVariants.length;
                    if (variantCount > 0) {
                        const vcEl = document.createElement('div');
                        vcEl.style.cssText = 'font-size:0.8rem;color:#6b7280;margin-top:3px;';
                        vcEl.textContent = `${variantCount} variant${variantCount !== 1 ? 's' : ''} catalogued in CIViC for ${apiGene.name}.`;
                        civicApiDiv.appendChild(vcEl);
                    }
                }).catch(() => {
                    civicApiDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">CIViC API unavailable.</div>';
                });
            }

            // ── Card: OncoKB ───────────────────────────────────────────────
            {
                const oncoCard = document.createElement('div');
                oncoCard.className = 'card';
                const oncoTitle = document.createElement('h3');
                oncoTitle.textContent = 'OncoKB';
                applyCardTheme(oncoCard, 'OncoKB');
                oncoCard.appendChild(oncoTitle);
                const oncoContent = document.createElement('div');
                oncoContent.className = 'card-content';

                const geneLink = document.createElement('a');
                geneLink.href = `https://www.oncokb.org/gene/${encodeURIComponent(gene)}`;
                geneLink.target = '_blank';
                geneLink.rel = 'noopener noreferrer';
                geneLink.textContent = `View ${gene} on OncoKB ↗`;
                oncoContent.appendChild(geneLink);

                if (altType) {
                    const altNote = document.createElement('div');
                    altNote.style.cssText = 'font-size:0.85rem;color:#6b7280;margin-top:0.3rem;';
                    altNote.textContent = `Gene-level view. For specific alteration "${altType}", use the link above to search OncoKB.`;
                    oncoContent.appendChild(altNote);
                }

                oncoCard.appendChild(oncoContent);
                if (cardsContainer) cardsContainer.appendChild(oncoCard);
            }

            // ── Card: PubMed ───────────────────────────────────────────────
            {
                const PUBMED_LIMIT = 5;
                const pmCard = document.createElement('div');
                pmCard.className = 'card';
                const pmTitle = document.createElement('h3');
                pmTitle.textContent = 'PubMed';
                applyCardTheme(pmCard, 'PubMed');
                pmCard.appendChild(pmTitle);
                const pmContent = document.createElement('div');
                pmContent.className = 'card-content';

                const pmSearchTerm = altType ? `${gene} ${altType}` : gene;
                const pmTumorSearchTerm = tumorType ? `${gene} ${tumorType}` : '';
                const pmVariantTumorSearchTerm = (tumorType && altType) ? `${gene} ${altType} ${tumorType}` : '';
                const hasPmTumorTab = !!(pmTumorSearchTerm && pmSearchTerm !== pmTumorSearchTerm);
                const hasPmVariantTumorTab = !!(pmVariantTumorSearchTerm
                    && pmVariantTumorSearchTerm !== pmSearchTerm
                    && pmVariantTumorSearchTerm !== pmTumorSearchTerm);

                const buildPmResultsPanel = (container, searchTerm, extraKey) => {
                    const queryUrl = searchTerm
                        ? `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(searchTerm)}&sort=relevance`
                        : 'https://pubmed.ncbi.nlm.nih.gov/';
                    const linkEl = document.createElement('a');
                    linkEl.href = queryUrl;
                    linkEl.target = '_blank';
                    linkEl.rel = 'noopener noreferrer';
                    linkEl.textContent = 'Search PubMed ↗';
                    container.appendChild(linkEl);
                    if (searchTerm) {
                        const queryLabel = document.createElement('div');
                        queryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
                        queryLabel.textContent = `Query: "${searchTerm}"`;
                        container.appendChild(queryLabel);
                    }
                    const resultsDiv = document.createElement('div');
                    if (searchTerm) {
                        const spinner = document.createElement('div');
                        spinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                        spinner.textContent = 'Loading PubMed results…';
                        resultsDiv.appendChild(spinner);
                    }
                    container.appendChild(resultsDiv);
                    if (searchTerm) {
                        fetchPubmedArticles(searchTerm, PUBMED_LIMIT).then(({ total, articles }) => {
                            geneOnlyAiExtras[extraKey] = { query: searchTerm, total, articles };
                            resultsDiv.innerHTML = '';
                            if (total === 0 || articles.length === 0) {
                                resultsDiv.innerHTML = '<div style="font-size:0.85rem;color:#6b7280;">No PubMed results found.</div>';
                                return;
                            }
                            const countEl = document.createElement('div');
                            countEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:6px;';
                            countEl.textContent = total > PUBMED_LIMIT
                                ? `Showing ${articles.length} of ${total.toLocaleString()} results (most relevant)`
                                : `${articles.length} result${articles.length !== 1 ? 's' : ''}`;
                            resultsDiv.appendChild(countEl);
                            articles.forEach((art) => {
                                const artEl = document.createElement('div');
                                artEl.style.cssText = 'margin-bottom:8px;padding-bottom:8px;border-bottom:1px solid #f0f0f0;font-size:0.82rem;';
                                const titleLink = document.createElement('a');
                                titleLink.href = `https://pubmed.ncbi.nlm.nih.gov/${art.pmid}/`;
                                titleLink.target = '_blank';
                                titleLink.rel = 'noopener noreferrer';
                                titleLink.style.fontWeight = '600';
                                titleLink.textContent = art.title;
                                artEl.appendChild(titleLink);
                                const meta = document.createElement('div');
                                meta.style.cssText = 'color:#6b7280;margin-top:2px;';
                                const parts = [art.authors, art.journal, art.year].filter(Boolean);
                                meta.textContent = parts.join(' · ') + (art.pmid ? ` · PMID ${art.pmid}` : '');
                                artEl.appendChild(meta);
                                if (art.abstract) {
                                    const abstractDetails = document.createElement('details');
                                    abstractDetails.style.cssText = 'margin-top:4px;';
                                    const abstractSummaryEl = document.createElement('summary');
                                    abstractSummaryEl.style.cssText = 'font-size:0.80rem;color:#4b5563;cursor:pointer;padding:2px 0;list-style:revert;';
                                    abstractSummaryEl.textContent = 'Abstract';
                                    abstractDetails.appendChild(abstractSummaryEl);
                                    const abstractText = document.createElement('div');
                                    abstractText.style.cssText = 'font-size:0.80rem;color:#374151;margin-top:4px;line-height:1.5;padding:6px;background:#f9fafb;border-radius:4px;';
                                    abstractText.textContent = art.abstract;
                                    abstractDetails.appendChild(abstractText);
                                    artEl.appendChild(abstractDetails);
                                }
                                resultsDiv.appendChild(artEl);
                            });
                            if (total > PUBMED_LIMIT) {
                                const seeAll = document.createElement('a');
                                seeAll.href = queryUrl;
                                seeAll.target = '_blank';
                                seeAll.rel = 'noopener noreferrer';
                                seeAll.style.fontSize = '0.82rem';
                                seeAll.textContent = `See all ${total.toLocaleString()} results on PubMed ↗`;
                                resultsDiv.appendChild(seeAll);
                            }
                        }).catch(() => {
                            resultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">PubMed unavailable.</div>';
                        });
                    }
                };

                if (hasPmTumorTab) {
                    const tabBar = document.createElement('div');
                    tabBar.className = 'card-tabs';
                    const variantBtn = document.createElement('button');
                    variantBtn.type = 'button';
                    variantBtn.className = 'card-tab-btn active';
                    variantBtn.textContent = altType ? 'Gene + Alteration' : 'Gene';
                    const tumorBtn = document.createElement('button');
                    tumorBtn.type = 'button';
                    tumorBtn.className = 'card-tab-btn';
                    tumorBtn.textContent = 'Gene + Tumor Type';
                    tabBar.appendChild(variantBtn);
                    tabBar.appendChild(tumorBtn);

                    let variantTumorBtn = null;
                    let variantTumorPanel = null;
                    if (hasPmVariantTumorTab) {
                        variantTumorBtn = document.createElement('button');
                        variantTumorBtn.type = 'button';
                        variantTumorBtn.className = 'card-tab-btn';
                        variantTumorBtn.textContent = 'Gene + Alteration + Tumor Type';
                        tabBar.appendChild(variantTumorBtn);
                    }
                    pmContent.appendChild(tabBar);

                    const variantPanel = document.createElement('div');
                    variantPanel.className = 'card-tab-panel active';
                    const tumorPanel = document.createElement('div');
                    tumorPanel.className = 'card-tab-panel';
                    if (hasPmVariantTumorTab) {
                        variantTumorPanel = document.createElement('div');
                        variantTumorPanel.className = 'card-tab-panel';
                    }

                    const tabBtns = [variantBtn, tumorBtn, variantTumorBtn].filter(Boolean);
                    const tabPanels = [variantPanel, tumorPanel, variantTumorPanel].filter(Boolean);
                    const activateTab = (idx) => {
                        tabBtns.forEach((b, i) => b.classList.toggle('active', i === idx));
                        tabPanels.forEach((p, i) => p.classList.toggle('active', i === idx));
                    };
                    variantBtn.addEventListener('click', () => activateTab(0));
                    tumorBtn.addEventListener('click', () => activateTab(1));
                    if (variantTumorBtn) variantTumorBtn.addEventListener('click', () => activateTab(2));

                    buildPmResultsPanel(variantPanel, pmSearchTerm, 'pubmed');
                    buildPmResultsPanel(tumorPanel, pmTumorSearchTerm, 'pubmed_tumor_type');
                    pmContent.appendChild(variantPanel);
                    pmContent.appendChild(tumorPanel);
                    if (variantTumorPanel) {
                        buildPmResultsPanel(variantTumorPanel, pmVariantTumorSearchTerm, 'pubmed_variant_tumor_type');
                        pmContent.appendChild(variantTumorPanel);
                    }
                } else {
                    buildPmResultsPanel(pmContent, pmSearchTerm, 'pubmed');
                }

                pmCard.appendChild(pmContent);
                if (cardsContainer) cardsContainer.appendChild(pmCard);
            }

            // ── Card: FDA-Approved Drugs (by gene) ─────────────────────────
            {
                const fdaCard = document.createElement('div');
                fdaCard.className = 'card';
                const fdaTitle = document.createElement('h3');
                fdaTitle.textContent = 'FDA-Approved Drugs (by gene)';
                applyCardTheme(fdaCard, 'FDA-Approved Drugs (by gene)');
                fdaCard.appendChild(fdaTitle);
                const fdaContent = document.createElement('div');
                fdaContent.className = 'card-content';

                // Tab bar
                const fdaTabBar = document.createElement('div');
                fdaTabBar.className = 'card-tabs';
                const compDxBtn = document.createElement('button');
                compDxBtn.type = 'button';
                compDxBtn.className = 'card-tab-btn active';
                compDxBtn.textContent = 'Companion Dx';
                const openFdaBtn = document.createElement('button');
                openFdaBtn.type = 'button';
                openFdaBtn.className = 'card-tab-btn';
                openFdaBtn.textContent = 'openFDA Labels';
                fdaTabBar.appendChild(compDxBtn);
                fdaTabBar.appendChild(openFdaBtn);
                fdaContent.appendChild(fdaTabBar);

                const compDxPanel = document.createElement('div');
                compDxPanel.className = 'card-tab-panel active';
                const openFdaPanel = document.createElement('div');
                openFdaPanel.className = 'card-tab-panel';

                compDxBtn.addEventListener('click', () => {
                    compDxBtn.classList.add('active'); openFdaBtn.classList.remove('active');
                    compDxPanel.classList.add('active'); openFdaPanel.classList.remove('active');
                });
                openFdaBtn.addEventListener('click', () => {
                    openFdaBtn.classList.add('active'); compDxBtn.classList.remove('active');
                    openFdaPanel.classList.add('active'); compDxPanel.classList.remove('active');
                });

                fdaContent.appendChild(compDxPanel);
                fdaContent.appendChild(openFdaPanel);
                fdaCard.appendChild(fdaContent);
                if (cardsContainer) cardsContainer.appendChild(fdaCard);

                // Companion Dx panel
                const fdaCompDxUrl = 'https://www.fda.gov/medical-devices/in-vitro-diagnostics/list-cleared-or-approved-companion-diagnostic-devices-in-vitro-and-imaging-tools';
                const fdaLinkEl = document.createElement('a');
                fdaLinkEl.href = fdaCompDxUrl;
                fdaLinkEl.target = '_blank';
                fdaLinkEl.rel = 'noopener noreferrer';
                fdaLinkEl.textContent = 'FDA companion diagnostics list ↗';
                compDxPanel.appendChild(fdaLinkEl);

                const fdaQueryLabel = document.createElement('div');
                fdaQueryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 2px;';
                fdaQueryLabel.textContent = `Gene: ${gene}`;
                compDxPanel.appendChild(fdaQueryLabel);

                const fdaDisclaimer = document.createElement('div');
                fdaDisclaimer.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin:0 0 6px;font-style:italic;';
                fdaDisclaimer.textContent = 'Note: FDA list does not always note resistance mutations.';
                compDxPanel.appendChild(fdaDisclaimer);

                const fdaResultsDiv = document.createElement('div');
                const fdaSpinner = document.createElement('div');
                fdaSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                fdaSpinner.textContent = 'Loading FDA drug data…';
                fdaResultsDiv.appendChild(fdaSpinner);
                compDxPanel.appendChild(fdaResultsDiv);

                fetchFdaCompanionDiagnostics(gene).then((records) => {
                    fdaResultsDiv.innerHTML = '';
                    if (!records || records.length === 0) {
                        const noResults = document.createElement('div');
                        noResults.style.cssText = 'font-size:0.85rem;color:#6b7280;';
                        noResults.textContent = `No FDA companion diagnostic records found for ${gene}.`;
                        fdaResultsDiv.appendChild(noResults);
                        return;
                    }
                    const countEl = document.createElement('div');
                    countEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:6px;';
                    countEl.textContent = `${records.length} record${records.length !== 1 ? 's' : ''}`;
                    fdaResultsDiv.appendChild(countEl);

                    const table = document.createElement('table');
                    table.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.82rem;';
                    const thead = document.createElement('thead');
                    const headerRow = document.createElement('tr');
                    ['Drug', 'Disease', 'Biomarker detail'].forEach((col) => {
                        const th = document.createElement('th');
                        th.style.cssText = 'text-align:left;padding:4px 6px;border-bottom:2px solid #fca5a5;background:#fff7f7;color:#7f1d1d;font-weight:600;';
                        th.textContent = col;
                        headerRow.appendChild(th);
                    });
                    thead.appendChild(headerRow);
                    table.appendChild(thead);

                    const FDA_PREVIEW_ROWS = 5;
                    const buildRow = (rec, i) => {
                        const tr = document.createElement('tr');
                        tr.style.background = i % 2 === 0 ? '#fff' : '#fff7f7';
                        const cellStyle = 'padding:4px 6px;vertical-align:top;border-bottom:1px solid #fee2e2;';
                        const tdDrug = document.createElement('td');
                        tdDrug.style.cssText = cellStyle + 'font-weight:600;';
                        const drugsArr = Array.isArray(rec.therapy?.drugs) ? rec.therapy.drugs : [];
                        tdDrug.textContent = drugsArr.length
                            ? drugsArr.map(d => d.trade_name ? `${d.trade_name} (${d.generic_name})` : d.generic_name || d.raw || '').filter(Boolean).join(', ')
                            : (rec.therapy?.raw || '—');
                        const tdDisease = document.createElement('td');
                        tdDisease.style.cssText = cellStyle + 'color:#374151;';
                        const rawDisease = rec.indication?.raw || rec.indication?.disease || '—';
                        tdDisease.textContent = rawDisease.replace(/\s*[-–]\s*(Tissue|Plasma|Blood|Serum|Urine|FFPE|Fresh Frozen|Whole Blood|ctDNA)\s*$/i, '').trim() || rawDisease;
                        const tdDetail = document.createElement('td');
                        tdDetail.style.cssText = cellStyle + 'color:#374151;';
                        tdDetail.textContent = rec.biomarker?.details || rec.biomarker?.name || '—';
                        tr.appendChild(tdDrug);
                        tr.appendChild(tdDisease);
                        tr.appendChild(tdDetail);
                        return tr;
                    };

                    const tbody = document.createElement('tbody');
                    records.slice(0, FDA_PREVIEW_ROWS).forEach((rec, i) => tbody.appendChild(buildRow(rec, i)));
                    table.appendChild(tbody);
                    fdaResultsDiv.appendChild(table);

                    if (records.length > FDA_PREVIEW_ROWS) {
                        const details = document.createElement('details');
                        details.style.cssText = 'margin-top:2px;';
                        const summary = document.createElement('summary');
                        summary.style.cssText = 'font-size:0.82rem;color:#7f1d1d;cursor:pointer;padding:4px 2px;list-style:revert;';
                        summary.textContent = `Show ${records.length - FDA_PREVIEW_ROWS} more…`;
                        details.appendChild(summary);
                        const extraTable = document.createElement('table');
                        extraTable.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.82rem;';
                        const extraTbody = document.createElement('tbody');
                        records.slice(FDA_PREVIEW_ROWS).forEach((rec, i) => extraTbody.appendChild(buildRow(rec, FDA_PREVIEW_ROWS + i)));
                        extraTable.appendChild(extraTbody);
                        details.appendChild(extraTable);
                        fdaResultsDiv.appendChild(details);
                    }

                    const note = document.createElement('div');
                    note.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                    note.textContent = 'Source: FDA companion diagnostics list. Verify against current FDA labeling before clinical use.';
                    fdaResultsDiv.appendChild(note);
                }).catch(() => {
                    fdaResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">FDA drug data unavailable.</div>';
                });

                // openFDA Labels panel
                const ofLinkEl = document.createElement('a');
                ofLinkEl.href = 'https://open.fda.gov/apis/drug/label/';
                ofLinkEl.target = '_blank';
                ofLinkEl.rel = 'noopener noreferrer';
                ofLinkEl.textContent = 'openFDA drug label database ↗';
                openFdaPanel.appendChild(ofLinkEl);

                const ofQueryLabel = document.createElement('div');
                ofQueryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
                ofQueryLabel.textContent = `Searching indications_and_usage for: "${gene}"`;
                openFdaPanel.appendChild(ofQueryLabel);

                const ofResultsDiv = document.createElement('div');
                const ofSpinner = document.createElement('div');
                ofSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                ofSpinner.textContent = 'Loading openFDA results…';
                ofResultsDiv.appendChild(ofSpinner);
                openFdaPanel.appendChild(ofResultsDiv);

                fetchOpenFdaDrugLabels(gene).then(({ total, fetched, excluded, excludedCase, excludedBoundary, excludedNegation, results: ofResults }) => {
                    geneOnlyAiExtras.openfda = { gene, total, fetched, excluded, excludedCase, excludedBoundary, excludedNegation, results: ofResults };
                    ofResultsDiv.innerHTML = '';
                    if (!ofResults || ofResults.length === 0) {
                        ofResultsDiv.innerHTML = `<div style="font-size:0.85rem;color:#6b7280;">No openFDA drug label results found for ${gene}.</div>`;
                        return;
                    }
                    const OF_PREVIEW = 7;
                    const ofCountEl = document.createElement('div');
                    ofCountEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:8px;';
                    const shownOf = `${Math.min(ofResults.length, OF_PREVIEW)} of ${ofResults.length} result${ofResults.length !== 1 ? 's' : ''}`;
                    const excludedReasons = [];
                    if (excludedCase) excludedReasons.push(`${excludedCase} case mismatch`);
                    if (excludedBoundary) excludedReasons.push(`${excludedBoundary} substring of larger word`);
                    if (excludedNegation) excludedReasons.push(`${excludedNegation} negated mention`);
                    const excludedNote = excludedReasons.length ? ` (excluded: ${excludedReasons.join(', ')})` : '';
                    ofCountEl.textContent = shownOf + excludedNote;
                    ofResultsDiv.appendChild(ofCountEl);

                    const buildDrugEl = (item) => {
                        const drugEl = document.createElement('div');
                        drugEl.style.cssText = 'margin-bottom:8px;padding:8px;background:#fff7f7;border:1px solid #fee2e2;border-radius:4px;font-size:0.82rem;';
                        const nameEl = document.createElement('div');
                        nameEl.style.cssText = 'font-weight:700;color:#7f1d1d;margin-bottom:2px;';
                        const brandPart = item.brand_name ? item.brand_name : '';
                        const genericPart = item.generic_name ? (brandPart ? `(${item.generic_name})` : item.generic_name) : '';
                        nameEl.textContent = [brandPart, genericPart].filter(Boolean).join(' ') || 'Unknown drug';
                        drugEl.appendChild(nameEl);
                        const metaParts = [item.manufacturer, item.route].filter(Boolean);
                        if (metaParts.length) {
                            const metaEl = document.createElement('div');
                            metaEl.style.cssText = 'color:#6b7280;margin-bottom:4px;';
                            metaEl.textContent = metaParts.join(' · ');
                            drugEl.appendChild(metaEl);
                        }
                        if (item.indications_and_usage) {
                            const indDetails = document.createElement('details');
                            indDetails.style.cssText = 'margin-top:2px;';
                            const indSummary = document.createElement('summary');
                            indSummary.style.cssText = 'cursor:pointer;color:#991b1b;font-weight:600;font-size:0.8rem;list-style:revert;padding:2px 0;';
                            indSummary.textContent = 'Indications & Usage';
                            indDetails.appendChild(indSummary);
                            const indText = document.createElement('div');
                            indText.style.cssText = 'margin-top:4px;color:#374151;line-height:1.5;white-space:pre-wrap;font-size:0.79rem;max-height:300px;overflow-y:auto;';
                            indText.textContent = item.indications_and_usage;
                            indDetails.appendChild(indText);
                            drugEl.appendChild(indDetails);
                        }
                        return drugEl;
                    };

                    ofResults.slice(0, OF_PREVIEW).forEach((item) => ofResultsDiv.appendChild(buildDrugEl(item)));
                    if (ofResults.length > OF_PREVIEW) {
                        const ofMoreDetails = document.createElement('details');
                        ofMoreDetails.style.cssText = 'margin-top:2px;';
                        const ofMoreSummary = document.createElement('summary');
                        ofMoreSummary.style.cssText = 'font-size:0.82rem;color:#7f1d1d;cursor:pointer;padding:4px 2px;list-style:revert;';
                        ofMoreSummary.textContent = `Show ${ofResults.length - OF_PREVIEW} more…`;
                        ofMoreDetails.appendChild(ofMoreSummary);
                        ofResults.slice(OF_PREVIEW).forEach((item) => ofMoreDetails.appendChild(buildDrugEl(item)));
                        ofResultsDiv.appendChild(ofMoreDetails);
                    }
                    const ofNote = document.createElement('div');
                    ofNote.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                    ofNote.textContent = 'Source: openFDA drug label API. Results show labels mentioning this gene in indications. Verify against current FDA labeling.';
                    ofResultsDiv.appendChild(ofNote);
                }).catch(() => {
                    ofResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">openFDA data unavailable.</div>';
                });
            }

            // ── Card: Clinical Trials ──────────────────────────────────────
            {
                const CT_PREVIEW = 5;
                const ctCard = document.createElement('div');
                ctCard.className = 'card';
                const ctTitle = document.createElement('h3');
                ctTitle.textContent = 'Clinical Trials';
                applyCardTheme(ctCard, 'Clinical Trials');
                ctCard.appendChild(ctTitle);
                const ctContent = document.createElement('div');
                ctContent.className = 'card-content';

                const ctBaseUrl = `https://clinicaltrials.gov/search?query=${encodeURIComponent([gene, tumorType].filter(Boolean).join(' '))}&recrs=b&type=Intr`;
                const ctLinkEl = document.createElement('a');
                ctLinkEl.href = ctBaseUrl;
                ctLinkEl.target = '_blank';
                ctLinkEl.rel = 'noopener noreferrer';
                ctLinkEl.textContent = 'Search ClinicalTrials.gov ↗';
                ctContent.appendChild(ctLinkEl);

                const ctQueryLabel = document.createElement('div');
                ctQueryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
                const queryParts = [gene, tumorType].filter(Boolean);
                ctQueryLabel.textContent = `Query: "${queryParts.join(' + ')}" · Interventional · Recruiting · Phase 2+ · US`;
                ctContent.appendChild(ctQueryLabel);

                const ctResultsDiv = document.createElement('div');
                const ctSpinner = document.createElement('div');
                ctSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                ctSpinner.textContent = 'Loading clinical trials…';
                ctResultsDiv.appendChild(ctSpinner);
                ctContent.appendChild(ctResultsDiv);
                ctCard.appendChild(ctContent);
                if (cardsContainer) cardsContainer.appendChild(ctCard);

                fetchClinicalTrials(gene, tumorType).then(({ total, studies }) => {
                    ctResultsDiv.innerHTML = '';
                    if (total === 0 || studies.length === 0) {
                        ctResultsDiv.innerHTML = '<div style="font-size:0.85rem;color:#6b7280;">No recruiting Phase 2+ interventional trials found in the US.</div>';
                        return;
                    }
                    const countEl = document.createElement('div');
                    countEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:8px;';
                    countEl.textContent = `${total} recruiting trial${total !== 1 ? 's' : ''} found`;
                    ctResultsDiv.appendChild(countEl);

                    const buildTrialEl = (trial) => {
                        const el = document.createElement('div');
                        el.style.cssText = 'margin-bottom:10px;padding-bottom:10px;border-bottom:1px solid #ccfbf1;font-size:0.82rem;';
                        const titleLine = document.createElement('div');
                        if (trial.url) {
                            const link = document.createElement('a');
                            link.href = trial.url;
                            link.target = '_blank';
                            link.rel = 'noopener noreferrer';
                            link.style.fontWeight = '600';
                            link.textContent = trial.title || trial.nctId;
                            titleLine.appendChild(link);
                        } else {
                            titleLine.style.fontWeight = '600';
                            titleLine.textContent = trial.title || trial.nctId;
                        }
                        el.appendChild(titleLine);
                        const meta = document.createElement('div');
                        meta.style.cssText = 'color:#6b7280;margin-top:3px;';
                        const phaseParts = Array.isArray(trial.phases) && trial.phases.length
                            ? trial.phases.map(p => String(p).replace('PHASE', 'Phase ')).join('/')
                            : 'Phase N/A';
                        const drugNames = (trial.interventions || [])
                            .filter(i => i.type === 'DRUG' || i.type === 'BIOLOGICAL' || i.type === 'COMBINATION_PRODUCT')
                            .map(i => i.name)
                            .filter(Boolean)
                            .slice(0, 4);
                        const metaParts = [
                            trial.nctId,
                            phaseParts,
                            drugNames.length ? drugNames.join(', ') : null,
                            trial.usLocationCount ? `${trial.usLocationCount} US site${trial.usLocationCount !== 1 ? 's' : ''}` : null
                        ].filter(Boolean);
                        meta.textContent = metaParts.join(' · ');
                        el.appendChild(meta);
                        const hasExpandable = trial.briefSummary || (trial.conditions && trial.conditions.length > 0) || trial.inclusionCriteria;
                        if (hasExpandable) {
                            const trialDetails = document.createElement('details');
                            trialDetails.style.cssText = 'margin-top:5px;';
                            const trialDetailsSummary = document.createElement('summary');
                            trialDetailsSummary.style.cssText = 'font-size:0.80rem;color:#0f766e;cursor:pointer;padding:2px 0;list-style:revert;';
                            trialDetailsSummary.textContent = 'Summary, conditions & eligibility';
                            trialDetails.appendChild(trialDetailsSummary);
                            const expandContent = document.createElement('div');
                            expandContent.style.cssText = 'font-size:0.80rem;color:#374151;margin-top:4px;line-height:1.5;padding:6px;background:#f0fdfa;border-radius:4px;';
                            if (trial.conditions && trial.conditions.length > 0) {
                                const condLabel = document.createElement('div');
                                condLabel.style.cssText = 'font-weight:600;margin-bottom:2px;';
                                condLabel.textContent = 'Conditions:';
                                expandContent.appendChild(condLabel);
                                const condVal = document.createElement('div');
                                condVal.style.cssText = 'margin-bottom:6px;';
                                condVal.textContent = trial.conditions.join(', ');
                                expandContent.appendChild(condVal);
                            }
                            if (trial.briefSummary) {
                                const summLabel = document.createElement('div');
                                summLabel.style.cssText = 'font-weight:600;margin-bottom:2px;';
                                summLabel.textContent = 'Summary:';
                                expandContent.appendChild(summLabel);
                                const summVal = document.createElement('div');
                                summVal.style.cssText = 'margin-bottom:6px;';
                                summVal.textContent = trial.briefSummary;
                                expandContent.appendChild(summVal);
                            }
                            if (trial.inclusionCriteria) {
                                const inclLabel = document.createElement('div');
                                inclLabel.style.cssText = 'font-weight:600;margin-bottom:2px;';
                                inclLabel.textContent = 'Inclusion Criteria:';
                                expandContent.appendChild(inclLabel);
                                const inclVal = document.createElement('pre');
                                inclVal.style.cssText = 'white-space:pre-wrap;font-family:inherit;margin:0;';
                                inclVal.textContent = trial.inclusionCriteria;
                                expandContent.appendChild(inclVal);
                            }
                            trialDetails.appendChild(expandContent);
                            el.appendChild(trialDetails);
                        }
                        return el;
                    };

                    studies.slice(0, CT_PREVIEW).forEach(trial => ctResultsDiv.appendChild(buildTrialEl(trial)));
                    if (studies.length > CT_PREVIEW) {
                        const moreDetails = document.createElement('details');
                        moreDetails.style.cssText = 'margin-top:4px;';
                        const moreSummary = document.createElement('summary');
                        moreSummary.style.cssText = 'font-size:0.82rem;color:#0f766e;cursor:pointer;padding:4px 2px;list-style:revert;';
                        moreSummary.textContent = `Show ${studies.length - CT_PREVIEW} more…`;
                        moreDetails.appendChild(moreSummary);
                        studies.slice(CT_PREVIEW).forEach(trial => moreDetails.appendChild(buildTrialEl(trial)));
                        ctResultsDiv.appendChild(moreDetails);
                    }
                    const ctNote = document.createElement('div');
                    ctNote.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                    ctNote.textContent = 'Source: ClinicalTrials.gov. Eligibility criteria not evaluated — verify before clinical use.';
                    ctResultsDiv.appendChild(ctNote);
                }).catch(() => {
                    ctResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">Clinical trials data unavailable.</div>';
                });
            }

            // ── Card: Guidelines ───────────────────────────────────────────
            {
                const GUIDELINES_BASE = 'https://drdoubleb.com/guidelines';
                const guidelinesCard = document.createElement('div');
                guidelinesCard.className = 'card';
                const guidelinesTitle = document.createElement('h3');
                guidelinesTitle.textContent = 'Guidelines';
                applyCardTheme(guidelinesCard, 'Guidelines');
                guidelinesCard.appendChild(guidelinesTitle);

                const guidelinesContent = document.createElement('div');
                guidelinesContent.className = 'card-content';

                const guidelinesIntro = document.createElement('p');
                guidelinesIntro.style.cssText = 'font-size:0.85rem;color:#6b7280;margin-bottom:8px;';
                guidelinesIntro.textContent = 'Select a cancer type to retrieve guideline recommendations for this gene.';
                guidelinesContent.appendChild(guidelinesIntro);

                const dropdownRow = document.createElement('div');
                dropdownRow.style.cssText = 'display:flex;align-items:center;gap:8px;margin-bottom:10px;';

                const cancerSelect = document.createElement('select');
                cancerSelect.style.cssText = 'font-size:0.85rem;padding:4px 8px;border:1px solid #d1d5db;border-radius:4px;flex:1;max-width:320px;';

                const defaultOpt = document.createElement('option');
                defaultOpt.value = '';
                defaultOpt.textContent = 'Loading cancer types…';
                defaultOpt.disabled = true;
                defaultOpt.selected = true;
                cancerSelect.appendChild(defaultOpt);
                dropdownRow.appendChild(cancerSelect);
                guidelinesContent.appendChild(dropdownRow);

                const guidelinesResults = document.createElement('div');
                guidelinesContent.appendChild(guidelinesResults);

                guidelinesCard.appendChild(guidelinesContent);
                if (cardsContainer) cardsContainer.appendChild(guidelinesCard);

                (async () => {
                    try {
                        const resp = await fetch(`${GUIDELINES_BASE}/api/cancer-types.php`);
                        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
                        const data = await resp.json();
                        cancerSelect.innerHTML = '';
                        const placeholder = document.createElement('option');
                        placeholder.value = '';
                        placeholder.textContent = `— Select cancer type (${data.count} available) —`;
                        placeholder.disabled = true;
                        placeholder.selected = true;
                        cancerSelect.appendChild(placeholder);
                        (data.cancer_types || []).forEach(ct => {
                            const opt = document.createElement('option');
                            opt.value = ct.name;
                            opt.textContent = ct.name + (ct.record_count ? ` (${ct.record_count} records)` : '');
                            cancerSelect.appendChild(opt);
                        });
                        if (tumorType) {
                            const tumorLower = tumorType.toLowerCase().trim();
                            const matched = (data.cancer_types || []).find(ct =>
                                (ct.aliases || []).some(a => a.toLowerCase() === tumorLower) ||
                                ct.name.toLowerCase() === tumorLower ||
                                ct.name.toLowerCase().includes(tumorLower) ||
                                tumorLower.includes(ct.name.toLowerCase())
                            );
                            if (matched) {
                                cancerSelect.value = matched.name;
                                cancerSelect.dispatchEvent(new Event('change'));
                            }
                        }
                    } catch (err) {
                        cancerSelect.innerHTML = '';
                        const errOpt = document.createElement('option');
                        errOpt.textContent = 'Failed to load cancer types';
                        errOpt.disabled = true;
                        errOpt.selected = true;
                        cancerSelect.appendChild(errOpt);
                    }
                })();

                cancerSelect.addEventListener('change', async () => {
                    const selectedCancer = cancerSelect.value;
                    if (!selectedCancer) return;
                    if (!gene) {
                        guidelinesResults.innerHTML = '<div style="font-size:0.85rem;color:#9ca3af;">No gene available for guideline lookup.</div>';
                        return;
                    }
                    guidelinesResults.innerHTML = '<div style="font-size:0.85rem;color:#6b7280;padding:4px 0;">Loading guidelines…</div>';
                    try {
                        const params = new URLSearchParams({ cancer: selectedCancer, gene });
                        const resp = await fetch(`${GUIDELINES_BASE}/api/search.php?${params}`);
                        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
                        const data = await resp.json();

                        geneOnlyAiExtras.guidelines = {
                            cancer_type: selectedCancer,
                            gene,
                            query: data.query,
                            count: data.count,
                            results: data.results
                        };

                        guidelinesResults.innerHTML = '';

                        if (!data.results || data.results.length === 0) {
                            const noResults = document.createElement('div');
                            noResults.style.cssText = 'font-size:0.85rem;color:#9ca3af;';
                            noResults.textContent = `No guideline records found for ${gene} in ${selectedCancer}.`;
                            guidelinesResults.appendChild(noResults);
                            return;
                        }

                        const countEl = document.createElement('div');
                        countEl.style.cssText = 'font-size:0.85rem;font-weight:600;color:#3f6212;margin-bottom:8px;';
                        countEl.textContent = `${data.count} guideline record${data.count !== 1 ? 's' : ''} for ${gene} in ${selectedCancer}`;
                        guidelinesResults.appendChild(countEl);

                        const renderKV = (obj, container, depth) => {
                            Object.entries(obj).forEach(([key, val]) => {
                                if (val === null || val === undefined || val === '' || (Array.isArray(val) && val.length === 0) || (typeof val === 'object' && !Array.isArray(val) && Object.keys(val).length === 0)) return;
                                const fmtKey = key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
                                if (typeof val === 'object' && !Array.isArray(val)) {
                                    const subTitle = document.createElement('div');
                                    subTitle.style.cssText = `padding-left:${depth * 12}px;font-weight:600;color:#374151;margin-top:3px;`;
                                    subTitle.textContent = fmtKey + ':';
                                    container.appendChild(subTitle);
                                    renderKV(val, container, depth + 1);
                                } else if (Array.isArray(val) && val.some(v => v !== null && typeof v === 'object')) {
                                    const subTitle = document.createElement('div');
                                    subTitle.style.cssText = `padding-left:${depth * 12}px;font-weight:600;color:#374151;margin-top:3px;`;
                                    subTitle.textContent = fmtKey + ':';
                                    container.appendChild(subTitle);
                                    val.forEach(item => {
                                        if (item !== null && typeof item === 'object') {
                                            renderKV(item, container, depth + 1);
                                        } else if (item !== null && item !== '') {
                                            const row = document.createElement('div');
                                            row.style.cssText = `padding-left:${(depth + 1) * 12}px;margin:2px 0;color:#1f2937;`;
                                            row.textContent = String(item);
                                            container.appendChild(row);
                                        }
                                    });
                                } else {
                                    const row = document.createElement('div');
                                    row.style.cssText = `padding-left:${depth * 12}px;margin:2px 0;`;
                                    const labelEl = document.createElement('span');
                                    labelEl.style.cssText = 'font-weight:600;color:#374151;';
                                    labelEl.textContent = fmtKey + ': ';
                                    const valEl = document.createElement('span');
                                    valEl.style.cssText = 'color:#1f2937;';
                                    valEl.textContent = Array.isArray(val) ? val.filter(v => v !== null && v !== '').join(', ') : String(val);
                                    row.appendChild(labelEl);
                                    row.appendChild(valEl);
                                    container.appendChild(row);
                                }
                            });
                        };

                        const renderSection = (title, obj, container) => {
                            if (!obj || typeof obj !== 'object' || Object.keys(obj).length === 0) return;
                            const sectionTitle = document.createElement('div');
                            sectionTitle.style.cssText = 'font-weight:700;color:#166534;margin:6px 0 2px;font-size:0.79rem;text-transform:uppercase;letter-spacing:0.04em;border-top:1px solid #bbf7d0;padding-top:5px;';
                            sectionTitle.textContent = title;
                            container.appendChild(sectionTitle);
                            renderKV(obj, container, 0);
                        };

                        data.results.forEach((record, idx) => {
                            const recordEl = document.createElement('div');
                            recordEl.style.cssText = 'margin-bottom:8px;padding:8px;background:#f0fdf4;border:1px solid #bbf7d0;border-radius:4px;font-size:0.82rem;';
                            const headerEl = document.createElement('div');
                            headerEl.style.cssText = 'font-weight:700;color:#166534;margin-bottom:4px;';
                            const markerDisplay = (record.marker && (record.marker.display_name || record.marker.symbol)) || `Record ${idx + 1}`;
                            const tumorDisplay = (record.tumor && (record.tumor.name || record.tumor.subtype)) || selectedCancer;
                            headerEl.textContent = `${markerDisplay} — ${tumorDisplay}`;
                            recordEl.appendChild(headerEl);

                            const buildKVRow = (label, value) => {
                                if (!value) return;
                                const row = document.createElement('div');
                                row.style.cssText = 'margin-bottom:3px;';
                                const lbl = document.createElement('span');
                                lbl.style.cssText = 'font-weight:600;';
                                lbl.textContent = label + ': ';
                                row.appendChild(lbl);
                                row.appendChild(document.createTextNode(value));
                                recordEl.appendChild(row);
                            };

                            if (record.therapy) {
                                const t = record.therapy;
                                if (t.recommended_therapies && t.recommended_therapies.length > 0) {
                                    buildKVRow('Recommended therapy', t.recommended_therapies.join(', '));
                                }
                                if (t.category) buildKVRow('Category', t.category);
                                if (t.context) buildKVRow('Therapy context', t.context);
                            }
                            if (record.testing) {
                                const te = record.testing;
                                if (te.recommended_test) buildKVRow('Recommended test', te.recommended_test);
                                if (te.context) buildKVRow('Testing context', te.context);
                                if (te.modality) buildKVRow('Modality', te.modality);
                            }

                            const detailsEl = document.createElement('details');
                            detailsEl.style.cssText = 'margin-top:6px;';
                            const detailsSummary = document.createElement('summary');
                            detailsSummary.style.cssText = 'font-size:0.80rem;color:#15803d;cursor:pointer;padding:2px 0;list-style:revert;';
                            detailsSummary.textContent = 'Full record details';
                            detailsEl.appendChild(detailsSummary);
                            const detailsBody = document.createElement('div');
                            detailsBody.style.cssText = 'font-size:0.79rem;padding:4px 2px;margin-top:4px;max-height:350px;overflow-y:auto;line-height:1.6;';
                            if (record.record_id) {
                                const idRow = document.createElement('div');
                                idRow.style.cssText = 'font-size:0.74rem;color:#9ca3af;margin-bottom:2px;';
                                idRow.textContent = `Record ID: ${record.record_id}`;
                                detailsBody.appendChild(idRow);
                            }
                            const sectionDefs = [['Tumor', record.tumor], ['Marker / Biomarker', record.marker], ['Testing', record.testing], ['Therapy', record.therapy], ['Guideline Source', record.guideline_metadata]];
                            sectionDefs.forEach(([title, obj]) => renderSection(title, obj, detailsBody));
                            const skipKeys = new Set(['record_id', 'tumor', 'marker', 'testing', 'therapy', 'guideline_metadata', 'dataset_file', 'dataset_record_index', 'dataset_name']);
                            Object.keys(record).filter(k => !skipKeys.has(k)).forEach(k => {
                                const val = record[k];
                                if (val && typeof val === 'object') {
                                    renderSection(k.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase()), val, detailsBody);
                                } else if (val !== null && val !== undefined && val !== '') {
                                    const row = document.createElement('div');
                                    row.style.cssText = 'margin:2px 0;';
                                    const labelEl = document.createElement('span');
                                    labelEl.style.cssText = 'font-weight:600;color:#374151;';
                                    labelEl.textContent = k.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase()) + ': ';
                                    row.appendChild(labelEl);
                                    row.appendChild(document.createTextNode(String(val)));
                                    detailsBody.appendChild(row);
                                }
                            });
                            detailsEl.appendChild(detailsBody);
                            recordEl.appendChild(detailsEl);
                            guidelinesResults.appendChild(recordEl);
                        });

                        const sourceNote = document.createElement('div');
                        sourceNote.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                        sourceNote.textContent = 'Source: drdoubleb.com/guidelines. For reference only — verify against current published guidelines.';
                        guidelinesResults.appendChild(sourceNote);
                    } catch (err) {
                        guidelinesResults.innerHTML = `<div style="font-size:0.85rem;color:#9ca3af;">Guidelines data unavailable: ${err.message}</div>`;
                        geneOnlyAiExtras.guidelines = { error: err.message, cancer_type: selectedCancer, gene };
                    }
                });
            }

            // ── Card: Optional AI Review ───────────────────────────────────
            {
                const aiCard = document.createElement('div');
                aiCard.className = 'card ai-review-card';
                const aiTitle = document.createElement('h3');
                aiTitle.textContent = 'Optional AI Review';
                applyCardTheme(aiCard, 'Optional AI Review');
                aiCard.appendChild(aiTitle);

                const aiContent = document.createElement('div');
                aiContent.className = 'card-content ai-review-content';

                const aiIntro = document.createElement('p');
                aiIntro.className = 'ai-review-intro';
                aiIntro.textContent = 'Send the retrieved gene-level data to OpenRouter for a structured draft interpretation. No request is sent until you click Run AI review.';
                aiContent.appendChild(aiIntro);

                const controls = document.createElement('div');
                controls.className = 'ai-review-controls';

                const modelLabel = document.createElement('label');
                modelLabel.textContent = 'Model';
                modelLabel.setAttribute('for', 'aiReviewModelSelectGene');
                controls.appendChild(modelLabel);

                const modelSelect = document.createElement('select');
                modelSelect.id = 'aiReviewModelSelectGene';
                OPENROUTER_MODEL_OPTIONS.forEach((modelName) => {
                    const opt = document.createElement('option');
                    opt.value = modelName;
                    opt.textContent = modelName;
                    modelSelect.appendChild(opt);
                });
                controls.appendChild(modelSelect);

                const runButton = document.createElement('button');
                runButton.type = 'button';
                runButton.textContent = 'Run AI review';
                controls.appendChild(runButton);
                aiContent.appendChild(controls);

                const notesWrap = document.createElement('div');
                notesWrap.className = 'ai-review-notes';
                const notesLabel = document.createElement('label');
                notesLabel.setAttribute('for', 'aiReviewNotesGene');
                notesLabel.textContent = 'Any additional notes for AI review';
                notesWrap.appendChild(notesLabel);
                const notesInput = document.createElement('textarea');
                notesInput.id = 'aiReviewNotesGene';
                notesInput.rows = 3;
                notesInput.placeholder = 'Optional — extra clinical context, prior therapies, specific questions, etc.';
                notesWrap.appendChild(notesInput);
                aiContent.appendChild(notesWrap);

                const aiContextInspector = document.createElement('details');
                aiContextInspector.style.cssText = 'margin:6px 0 2px;font-size:0.80rem;';
                const aiContextInspectorSummary = document.createElement('summary');
                aiContextInspectorSummary.style.cssText = 'cursor:pointer;color:#9ca3af;padding:2px 0;list-style:revert;font-size:0.79rem;';
                aiContextInspectorSummary.textContent = 'Context sent to AI (populated after run)';
                aiContextInspector.appendChild(aiContextInspectorSummary);
                const aiContextPre = document.createElement('pre');
                aiContextPre.style.cssText = 'font-size:0.73rem;white-space:pre-wrap;word-break:break-word;background:#f8fafc;border:1px solid #e5e7eb;padding:8px;border-radius:4px;margin-top:4px;max-height:400px;overflow-y:auto;color:#374151;';
                aiContextPre.textContent = 'Run AI review to populate this section.';
                aiContextInspector.appendChild(aiContextPre);
                aiContent.appendChild(aiContextInspector);

                const aiOutput = document.createElement('div');
                aiOutput.className = 'ai-review-output';
                aiContent.appendChild(aiOutput);
                aiCard.appendChild(aiContent);
                if (cardsContainer) cardsContainer.appendChild(aiCard);

                runButton.addEventListener('click', async () => {
                    runButton.disabled = true;
                    const previousText = runButton.textContent;
                    runButton.textContent = 'Running…';
                    aiOutput.innerHTML = '<div class="ai-review-loading">Gathering gene-level context for AI review…</div>';
                    try {
                        // Supplemental context for gene-only mode: skip coordinate-dependent calls
                        const [fdaRecords, clinicalTrialData] = await Promise.all([
                            fetchFdaCompanionDiagnostics(gene).catch(() => []),
                            fetchClinicalTrials(gene, tumorType).catch(() => ({ total: 0, studies: [] }))
                        ]);
                        const pubmedTerm = altType ? `${gene} ${altType}` : gene;
                        // Prefer the PubMed card's already-fetched data over re-fetching.
                        // Re-fetching alongside the card increases NCBI concurrency and can
                        // drop abstracts when efetch hits the rate limit.
                        const cachedPubmed = geneOnlyAiExtras.pubmed;
                        const pubmedPromise = (cachedPubmed && Array.isArray(cachedPubmed.articles) && cachedPubmed.articles.length > 0)
                            ? Promise.resolve({ total: cachedPubmed.total ?? cachedPubmed.articles.length, articles: cachedPubmed.articles })
                            : fetchPubmedArticles(pubmedTerm, 5).catch(() => ({ total: 0, articles: [] }));
                        const [pubmedData, civicData, openFdaData] = await Promise.all([
                            pubmedPromise,
                            fetchCivicApiData(gene, '').catch(() => null),
                            fetchOpenFdaDrugLabels(gene).catch(() => null)
                        ]);
                        const supplementalContext = {
                            ...geneOnlyAiExtras,
                            civic_api: civicData,
                            pubmed: pubmedData,
                            openfda_drug_labels: openFdaData
                        };
                        const userNotes = (notesInput.value || '').trim();
                        const aiContext = {
                            submitted_query: rawInput,
                            gene,
                            alteration_type: altType,
                            tumor_type: tumorType,
                            user_notes: userNotes || undefined,
                            fda_companion_diagnostics_records: fdaRecords,
                            clinical_trials: clinicalTrialData,
                            supplemental_card_data: supplementalContext
                        };
                        aiContextPre.textContent = JSON.stringify(aiContext, null, 2);
                        aiContextInspectorSummary.textContent = 'Context sent to AI';
                        const data = await fetchAiReview(aiContext, modelSelect.value);
                        renderAiReview(data.review, aiOutput);
                    } catch (err) {
                        aiOutput.innerHTML = '';
                        const errorEl = document.createElement('div');
                        errorEl.className = 'ai-review-error';
                        errorEl.textContent = `AI review unavailable: ${err.message}`;
                        aiOutput.appendChild(errorEl);
                    } finally {
                        runButton.disabled = false;
                        runButton.textContent = previousText;
                    }
                });
            }

            return; // Skip the variant annotation pipeline
        }
        // ── End gene-only mode ───────────────────────────────────────────────

        // Parse a triple‑coded protein change from the input query (e.g. p.Val600Glu or VAL600GLU).
        targetProtGlobal = null;
        {
            // Look for p. notation with three‑letter amino acid codes
            const mTriple = query.match(/p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})/);
            if (mTriple) {
                // Combine uppercase triple letters and position
                targetProtGlobal = (mTriple[1] + mTriple[2] + mTriple[3]).toUpperCase();
            }
        }
        // Parse a cDNA change from the input query (e.g. c.178_186del) to prioritise variant candidates.
        targetCdnaGlobal = null;
        {
            const cdMatch = query.match(/c\.[^\s]+/i);
            if (cdMatch) {
                targetCdnaGlobal = cdMatch[0].trim().toLowerCase();
            }
        }
        // Record whether the user's input was genomic (g.) prior to any conversion. We will
        // use this flag later to decide whether to display cDNA/protein from the Ensembl
        // variant recoder or fall back to MyVariant.info annotations.  This avoids a bug
        // where we accidentally treat a non-genomic input (e.g. BRAF V600E) as genomic
        // after conversion to g. notation and therefore ignore recoder-derived transcripts.
        const originalIsGenomic = isGenomicVariant(query);

        // Fetch transcript list via variant_recoder in parallel to other operations.  Always
        // attempt this call regardless of whether the input is genomic, cDNA or protein.
        // The variant_recoder can still return useful transcript annotations for complex
        // genomic variants that MyVariant.info does not index (e.g. delins).  If this
        // request fails, transcriptsFromRecoder will remain an empty array.
        try {
            transcriptsFromRecoder = await getTranscriptsList(query);
        } catch {
            transcriptsFromRecoder = [];
        }
        if (!query) return;
        statusEl.textContent = 'Processing...';
        resultSection.classList.add('hidden');
        try {
            let gVariant = query;
            let recoderData = null;
            let annotation = null;
            // Global candidate variant list. When the Ensembl variant recoder produces multiple genomic
            // representations (e.g. hgvsg or spdi), they will be converted to hg19 and stored here. This
            // array is defined outside of the recoder loop so that it can be referenced later, e.g. when
            // free‑text MyVariant searches fail and we want to try alternate g. notations such as delins.
            let candidateVariants = [];

            // If query is already in genomic HGVS format (chrN:g.posRef>Alt), try fetching directly.
            if (isGenomicVariant(query)) {
                console.log('[DEBUG] Detected genomic HGVS input, attempting direct fetch:', query);
                try {
                    annotation = await fetchMyVariant(query);
                    gVariant = query;
                    console.log('[DEBUG] Direct genomic fetch succeeded for', query);
                } catch (errDirectGenomic) {
                    console.log('[DEBUG] Direct genomic fetch failed, attempting liftover and retry:', errDirectGenomic);
                    try {
                        const lifted = await liftoverHg38ToHg19(query);
                        if (lifted !== query) {
                            console.log('[DEBUG] Liftover of genomic variant:', query, '->', lifted);
                            gVariant = lifted;
                            annotation = await fetchMyVariant(lifted);
                            console.log('[DEBUG] Fetch after liftover succeeded for', lifted);
                        }
                    } catch (liftoverErr) {
                        console.log('[DEBUG] Liftover and fetch after liftover failed:', liftoverErr);
                    }
                }
            }

            // If this appears to be a protein variant (e.g. contains p.) and not already genomic
            // or rsID, attempt a direct MyVariant search first. This often resolves to the
            // correct genomic variant more reliably than using variant_recoder.
            const isProteinVariant = /p\./i.test(query);
            const looksLikeNonGenomic = !isGenomicVariant(query) && !/^rs\d+/i.test(query);
            // We previously attempted a direct MyVariant search for protein variants here, but this rarely
            // returned useful results. The search has been removed in favour of always going through the
            // Ensembl Variant Recoder for non-genomic variants.

            // Attempt a direct MyVariant search for protein or cDNA queries before using
            // the Ensembl Variant Recoder. This improves support for well‑known variants
            // such as BRAF V600E and EGFR L858R. Try the original query and a version
            // with colons replaced by spaces. Only perform this search when the input
            // is non‑genomic and no annotation has been found yet.
            if (!annotation && looksLikeNonGenomic && (isProteinVariant || /c\./i.test(query))) {
                try {
                    // Try direct search using the query as provided (e.g. "BRAF:p.Val600Glu").
                    let hit = await queryMyVariantById(query);
                    if (!hit && query.includes(':')) {
                        // Replace colons with spaces for free‑text lookup (e.g. "BRAF p.Val600Glu").
                        const altQuery = query.replace(/:/g, ' ');
                        hit = await queryMyVariantById(altQuery);
                    }
                    if (hit && hit._id) {
                        console.log('[DEBUG] Direct MyVariant search hit', hit._id);
                        gVariant = hit._id;
                        annotation = await fetchMyVariant(hit._id);
                    }
                } catch (directErr) {
                    console.log('[DEBUG] Direct MyVariant search error:', directErr);
                }
            }

            /*
             * If we still don't have an annotation, fall back to using the Ensembl variant_recoder
             * to resolve to a genomic variant and/or retrieve transcript‑level annotations.
             * Even for genomic HGVS inputs we attempt the recoder because complex delins variants
             * may not be indexed by MyVariant.info and the recoder can still provide useful
             * transcript information. When an annotation is eventually found this block is
             * skipped.
             */
            if (!annotation) {
                if (!isGenomicVariant(query)) {
                    statusEl.textContent = 'Converting variant...';
                } else {
                    statusEl.textContent = 'Recoder fallback…';
                }
                try {
                    recoderData = await fetchVariantRecoder(query);
                } catch (errRecoder) {
                    // Warn in console but don't immediately fail; we will attempt a free-text MyVariant search below.
                    console.warn('Variant recoder failed; will attempt free-text search', errRecoder);
                }
            } else {
                // Annotation already exists; continue to fetch annotation details.
                statusEl.textContent = 'Fetching annotation...';
            }
            if (recoderData && recoderData.length > 0) {
                console.log('[DEBUG] Recoder returned data:', recoderData);
                // Use recoderData to convert to a list of genomic variants (g. notation).
                // We intentionally ignore any rsIDs returned by the recoder to avoid
                // selecting alternate alleles based on SNP identifiers. Instead we build a list of
                // candidate genomic variants from hgvsg and spdi entries, then test them
                // against the expected protein change.
                // candidateVariants is defined in a higher scope (line ~1253) to allow reuse
                // in later fallbacks (e.g. free‑text search). Do not redeclare it here.
                for (const item of recoderData) {
                    for (const key in item) {
                        const sub = item[key];
                        if (!sub) continue;
                        // hgvsg entries may be a string or an array
                        if (sub.hgvsg) {
                            const hgvsgList = Array.isArray(sub.hgvsg) ? sub.hgvsg : [sub.hgvsg];
                            for (const hgvsg of hgvsgList) {
                                try {
                                    const mv = convertHgvsgToMyVariant(hgvsg);
                                    const lifted = await liftoverHg38ToHg19(mv);
                                    console.log('[DEBUG] Converted hgvsg to MV and liftover:', hgvsg, '->', mv, '->', lifted);
                                    candidateVariants.push(lifted);
                                } catch {
                                    // skip conversion errors
                                }
                            }
                        }
                        // spdi entries may be string or array
                        if (sub.spdi) {
                            const spdiList = Array.isArray(sub.spdi) ? sub.spdi : [sub.spdi];
                            for (const spdi of spdiList) {
                                try {
                                    const mv = convertSpdiToMyVariant(spdi);
                                    const lifted = await liftoverHg38ToHg19(mv);
                                    console.log('[DEBUG] Converted SPDI to MV and liftover:', spdi, '->', mv, '->', lifted);
                                    candidateVariants.push(lifted);
                                } catch {
                                    // skip conversion errors
                                }
                            }
                        }
                    }
                }
                // Remove duplicates and ensure at least one candidate exists
                const uniqueCandidates = Array.from(new Set(candidateVariants));
                console.log('[DEBUG] Unique candidate variants (hg19 coordinates):', uniqueCandidates);
                // Expose uniqueCandidates globally for later fallback use (e.g. in free-text search)
                candidateVariants = uniqueCandidates;
                if (uniqueCandidates.length > 0) {
                    // Use targetProtGlobal (if set) to match candidate variants. This value is
                    // parsed from the user's query earlier (e.g. VAL600GLU).
                    let selectedAnn = null;
                    let selectedVar = null;
                    let firstAnn = null;
                    for (const cand of uniqueCandidates) {
                        let ann = null;
                        try {
                            ann = await fetchMyVariant(cand);
                            console.log('[DEBUG] Fetched annotation for candidate', cand);
                        } catch {
                            console.log('[DEBUG] Annotation fetch error for candidate', cand);
                            continue;
                        }
                        if (!ann) continue;
                        if (!firstAnn) {
                            firstAnn = ann;
                            selectedVar = cand;
                        }
                        // Determine whether this annotation matches the desired cDNA or protein change.
                        let cdnaMatch = false;
                        if (typeof targetCdnaGlobal !== 'undefined' && targetCdnaGlobal) {
                            // Check dbnsfp hgvsc list
                            if (ann.dbnsfp && ann.dbnsfp.hgvsc) {
                                const scList = Array.isArray(ann.dbnsfp.hgvsc) ? ann.dbnsfp.hgvsc : [ann.dbnsfp.hgvsc];
                                for (const h of scList) {
                                    const vPart = String(h).split(':').slice(1).join(':').trim().toLowerCase();
                                    if (vPart === targetCdnaGlobal) {
                                        cdnaMatch = true;
                                        break;
                                    }
                                }
                            }
                            // Check snpEff annotations (hgvs_c)
                            if (!cdnaMatch && ann.snpeff && ann.snpeff.ann) {
                                const annList = Array.isArray(ann.snpeff.ann) ? ann.snpeff.ann : [ann.snpeff.ann];
                                for (const a of annList) {
                                    if (a.hgvs_c && String(a.hgvs_c).trim().toLowerCase() === targetCdnaGlobal) {
                                        cdnaMatch = true;
                                        break;
                                    }
                                }
                            }
                            // Check dbsnp gene RNA HGVS strings
                            if (!cdnaMatch && ann.dbsnp && ann.dbsnp.gene) {
                                const geneList = Array.isArray(ann.dbsnp.gene) ? ann.dbsnp.gene : [ann.dbsnp.gene];
                                for (const g of geneList) {
                                    if (g.rnas) {
                                        const rnas = Array.isArray(g.rnas) ? g.rnas : [g.rnas];
                                        for (const r of rnas) {
                                            if (r.hgvs) {
                                                const vPart = String(r.hgvs).split(':').slice(1).join(':').replace(/=/g,'').trim().toLowerCase();
                                                if (vPart === targetCdnaGlobal) {
                                                    cdnaMatch = true;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if (cdnaMatch) break;
                                }
                            }
                        }
                        let matchProt = true;
                        if (targetProtGlobal) {
                            const annStr = JSON.stringify(ann).toUpperCase();
                            const tripleSearch = targetProtGlobal;
                            const singleSearch = tripleToSingle(tripleSearch);
                            matchProt = annStr.includes(tripleSearch) || (singleSearch && annStr.includes(singleSearch));
                        }
                        // Prioritise cDNA matches over protein matches; if cdnaMatch true, use this candidate.
                        if (cdnaMatch || matchProt) {
                            selectedAnn = ann;
                            selectedVar = cand;
                            break;
                        }
                    }
                    if (!selectedAnn && firstAnn) {
                        selectedAnn = firstAnn;
                    }
                    if (selectedAnn) {
                        gVariant = selectedVar;
                        annotation = selectedAnn;
                        console.log('[DEBUG] Selected candidate variant after scanning:', gVariant);
                } else {
                    // Even if no annotation was found for any candidate, use the first candidate variant as the
                    // genomic representation. This ensures that g. is populated with a proper genomic HGVS
                    // notation rather than echoing the original non-genomic query.
                    if (uniqueCandidates && uniqueCandidates.length > 0) {
                        gVariant = uniqueCandidates[0];
                    }
                    }
                }
            }
            // If no annotation yet, or recoder failed, attempt a free-text MyVariant search.
            // However, protein-only variants rarely succeed with a MyVariant.info free-text query.
            // If recoder returned some transcript data but annotation remains null, proceed with a fallback display
            if (!annotation) {
                // Attempt a fallback using recoder transcripts if available. If the initial transcript list
                // is empty but the recoder returned data, extract hgvsc/hgvsp entries directly from
                // recoderData. This covers cases where getTranscriptsList failed but recoderData exists.
                if ((!transcriptsFromRecoder || transcriptsFromRecoder.length === 0) && recoderData && recoderData.length > 0) {
                    try {
                        const recEntry = recoderData[0];
                        let objWithTranscripts = null;
                        if (recEntry.A && typeof recEntry.A === 'object' && (recEntry.A.hgvsc || recEntry.A.hgvsp)) {
                            objWithTranscripts = recEntry.A;
                        } else {
                            for (const [subKey, subVal] of Object.entries(recEntry)) {
                                if (subVal && typeof subVal === 'object' && (subVal.hgvsc || subVal.hgvsp)) {
                                    objWithTranscripts = subVal;
                                    break;
                                }
                            }
                        }
                        if (objWithTranscripts) {
                            const hgvscs = Array.isArray(objWithTranscripts.hgvsc) ? objWithTranscripts.hgvsc : (objWithTranscripts.hgvsc ? [objWithTranscripts.hgvsc] : []);
                            const hgvsp = Array.isArray(objWithTranscripts.hgvsp) ? objWithTranscripts.hgvsp : (objWithTranscripts.hgvsp ? [objWithTranscripts.hgvsp] : []);
                            const len = Math.max(hgvscs.length, hgvsp.length);
                            const list = [];
                            const maxEntries = 200;
                            for (let i = 0; i < Math.min(len, maxEntries); i++) {
                                const sc = hgvscs[i];
                                const sp = hgvsp[i];
                                if (sc) {
                                    const parts = String(sc).split(':');
                                    const transcriptId = parts[0];
                                    const cpart = parts.slice(1).join(':');
                                    let prot = '';
                                    if (sp) {
                                        const pparts = String(sp).split(':');
                                        prot = pparts.slice(1).join(':');
                                    }
                                    list.push({ transcript: transcriptId, cDNA: cpart, protein: prot });
                                }
                            }
                            if (list.length > 0) {
                                transcriptsFromRecoder = list;
                            }
                        }
                    } catch {
                        // ignore errors during extraction
                    }
                }
                if (transcriptsFromRecoder && transcriptsFromRecoder.length > 0) {
                    /*
                     * At this point the variant_recoder successfully returned transcript‑level cDNA/protein
                     * annotations but we were unable to retrieve a genomic annotation from MyVariant.info.
                     * MyVariant.info does not index every large delins variant in g. notation, but it does
                     * sometimes provide minimal annotations keyed off of the transcript:cDNA HGVS string. To
                     * maximise our chances of finding an annotation, attempt a query to MyVariant.info using
                     * the canonical transcript:cDNA combination. If a hit is found, use the resulting
                     * _id (which may be a g. representation) to retrieve a full annotation. Otherwise we
                     * fall back to a minimal annotation using the original query.
                     */
                    let foundViaTranscript = false;
                    try {
                        // Determine canonical transcript index. Supply a dummy source so that
                        // selectCanonicalTranscript can apply its scoring correctly. If the function
                        // throws or is unavailable, simply default to the first transcript.
                        let canonicalIdx = 0;
                        try {
                            const candidatesForCanonical = transcriptsFromRecoder.map(t => {
                                return { transcript: t.transcript, cDNA: t.cDNA, protein: t.protein, source: 'root' };
                            });
                            canonicalIdx = selectCanonicalTranscript(candidatesForCanonical, targetProtGlobal);
                            if (typeof canonicalIdx !== 'number' || canonicalIdx < 0 || canonicalIdx >= transcriptsFromRecoder.length) {
                                canonicalIdx = 0;
                            }
                        } catch {
                            canonicalIdx = 0;
                        }
                        const canonical = transcriptsFromRecoder[canonicalIdx];
                        if (canonical && canonical.transcript && canonical.cDNA) {
                            const transcriptQuery = `${canonical.transcript}:${canonical.cDNA}`;
                            let hit = await queryMyVariantById(transcriptQuery);
                            console.log('[DEBUG] Transcript‑level MyVariant search for', transcriptQuery, 'returned', hit);
                            // If no hit and there is a space or colon, try replacing colon with space as fallback
                            if (!hit && transcriptQuery.includes(':')) {
                                const altQuery = transcriptQuery.replace(/:/g, ' ');
                                hit = await queryMyVariantById(altQuery);
                                console.log('[DEBUG] Transcript‑level MyVariant search for', altQuery, 'returned', hit);
                            }
                            if (hit && hit._id) {
                                // Use the returned _id for a direct MyVariant fetch. This may be a g. or cDNA id.
                                gVariant = hit._id;
                                annotation = await fetchMyVariant(gVariant);
                                foundViaTranscript = true;
                            }
                        }
                    } catch (transcriptSearchErr) {
                        console.log('[DEBUG] Transcript‑level MyVariant search error:', transcriptSearchErr);
                    }
                    if (!foundViaTranscript) {
                        // Either no hit was found or an error occurred. Build a minimal annotation object to
                        // display the transcript information returned by the variant_recoder. Set _id to
                        // the normalised query so that users see the identifier they entered. Also attempt
                        // to populate the gene name from the query to support summary display.
                        annotation = { _id: query };
                        // Attempt to determine the gene symbol for the minimal annotation.  First
                        // extract a leading token from the normalised query (before any colon or space).
                        let geneName = null;
                        {
                            const m = query.match(/^(\w+)[:\s]/);
                            if (m) {
                                geneName = m[1].toUpperCase();
                            }
                        }
                        // If the extracted geneName begins with a chromosome prefix (e.g. "CHR7") or
                        // appears to be an accession (NC_), treat it as unreliable.  Prefer the
                        // gene hint captured from the user's original input, if available.
                        const looksLikeChrom = geneName && /^CHR[0-9XYMT]+$/i.test(geneName);
                        const looksLikeAccession = geneName && /^NC_/.test(geneName);
                        if ((looksLikeChrom || looksLikeAccession || !geneName) && geneHintGlobal) {
                            geneName = geneHintGlobal;
                        }
                        if (geneName) {
                            annotation.dbnsfp = { genename: geneName };
                        }
                    }
                } else if (!isProteinVariant) {
                    // Attempt a free-text MyVariant lookup for non-protein variants
                    statusEl.textContent = 'Searching MyVariant.info...';
                    try {
                        // Attempt free-text search by the query itself. If no hit, try replacing ':' with ' ' as many
                        // free-text queries expect a space between gene and HGVS.
                        let hit = await queryMyVariantById(query);
                        console.log('[DEBUG] Free-text MyVariant search for', query, 'returned', hit);
                        if (!hit && query.includes(':')) {
                            const altQuery = query.replace(/:/g, ' ');
                            hit = await queryMyVariantById(altQuery);
                            console.log('[DEBUG] Free-text MyVariant search for', altQuery, 'returned', hit);
                        }
                        if (hit && hit._id) {
                            gVariant = hit._id;
                            annotation = await fetchMyVariant(gVariant);
                        } else {
                            // If free-text search fails, attempt to use alternate g. notations derived from Ensembl
                            // Some complex delins variants are indexed by MyVariant using a delins range rather than a simple ref>alt substitution.
                            let altFound = false;
                            if (typeof candidateVariants !== 'undefined' && candidateVariants.length > 0) {
                                for (const cv of candidateVariants) {
                                    try {
                                        const ann = await fetchMyVariant(cv);
                                        if (ann && ann._id) {
                                            gVariant = ann._id;
                                            annotation = ann;
                                            altFound = true;
                                            console.log('[DEBUG] Found annotation via alternate candidate variant', cv, ann);
                                            break;
                                        }
                                    } catch {}
                                }
                            }
                            if (!altFound) {
                                // As a last resort, attempt to query the Ensembl VEP HGVS endpoint on the GRCh37 server
                                // to retrieve transcript consequences for this genomic variant. This helps cover cases
                                // where MyVariant.info does not index the deletion/indel but Ensembl VEP still recognises it.
                                try {
                                    // Use the currently normalised gVariant if available; otherwise fall back to the original query.
                                    const hgvsForVep = gVariant || query;
                                    const vepRes = await fetchVepHgvsHg19(hgvsForVep);
                                    const vepConsequences = vepRes.consequences || [];
                                    if (vepConsequences.length > 0) {
                                        // Build a minimal annotation using the first gene symbol from the VEP data. The dbnsfp
                                        // genename field is used for summary display. The _id is set to the original query.
                                        annotation = { _id: query };
                                        let geneSym = await resolveGeneSymbolFromVep(vepConsequences, recoderData);
                                        // If upstream sources still do not provide a usable gene symbol,
                                        // then use the optional user hint as the final fallback.
                                        if ((!geneSym || isChromosomeLikeGeneSymbol(geneSym)) && geneHintGlobal) {
                                            geneSym = geneHintGlobal;
                                        }
                                        if (geneSym) {
                                            annotation.dbnsfp = { genename: geneSym };
                                        }
                                        // Populate transcriptsFromRecoder with transcript identifiers from the VEP data so
                                        // that the summary display can include a list of transcripts. Mark these as coming
                                        // from the VEP source.
                                        transcriptsFromRecoder = vepConsequences.map(c => {
                                            return {
                                                transcript: c.transcript_id || '',
                                                cDNA: '',
                                                protein: '',
                                                source: 'vep'
                                            };
                                        });
                                        altFound = true;
                                    }
                                } catch (vepErr) {
                                    // Ignore VEP errors and fall through to throwing an error below
                                }
                                if (!altFound) {
                                    // If no alternative candidate variant matches were found and VEP fallback failed,
                                    // throw a more helpful error message. The original free‑text search error is
                                    // replaced with guidance suggesting that the provided genomic position or
                                    // reference allele may not exist in the chosen genome build.
                                    throw new Error('Variant not found. Please verify the genomic coordinate and reference allele.');
                                }
                            }
                        }
                    } catch (errSearch) {
                        console.log('[DEBUG] Free-text MyVariant search error:', errSearch);
                        throw new Error(errSearch.message || 'Variant not found');
                    }
                } else {
                    // No annotation and no recoder transcripts for a protein‑only variant
                    throw new Error('Variant not found via Ensembl Variant Recoder');
                }
            }
            // Display annotation
            statusEl.textContent = 'Annotation retrieved';
            summaryTable.innerHTML = '';
            const summaryRows = buildSummary(annotation, gVariant);
            summaryRows.forEach(row => {
                const tr = document.createElement('tr');
                const th = document.createElement('th');
                th.textContent = row.name;
                const td = document.createElement('td');
                td.textContent = String(row.value);
                tr.appendChild(th);
                tr.appendChild(td);
                summaryTable.appendChild(tr);
            });
            // Build detailed sections
            let detailsData = buildDetailsData(annotation, rawInput, gVariant);
            const detailsContainer = document.getElementById('detailsContainer');
            const aiReviewExtras = {};

            // Attempt to fetch extended COSMIC data from a custom API if configured.
            const COSMIC_ENDPOINT = window.COSMIC_API_ENDPOINT || null;
            const COSMIC_META_URL = window.COSMIC_META_URL || null;
            if (COSMIC_ENDPOINT) {
                try {
                    const cosmicPromise = fetchWithTimeout(`${COSMIC_ENDPOINT}?hgvsg=${encodeURIComponent(gVariant)}`, {}, API_TIMEOUT_MS.cosmic);
                    const metaPromise = COSMIC_META_URL
                        ? fetchWithTimeout(COSMIC_META_URL, {}, API_TIMEOUT_MS.cosmicMeta)
                        : Promise.resolve(null);
                    const [cosmicResult, metaResult] = await Promise.allSettled([cosmicPromise, metaPromise]);
                    if (cosmicResult.status === 'fulfilled' && cosmicResult.value.ok) {
                        const cosmicData = await cosmicResult.value.json();
                        // Optionally load cosmic meta to compute frequencies
                        let meta = null;
                        if (metaResult.status === 'fulfilled' && metaResult.value && metaResult.value.ok) {
                            meta = await metaResult.value.json();
                        }
                        const cosmicItems = {};
                        // Basic COSMIC metrics
                        cosmicItems['Total Tumors'] = cosmicData.COSMIC_COUNT;
                        if (cosmicData.COSMIC_PROTEIN) cosmicItems['Protein Change'] = cosmicData.COSMIC_PROTEIN;
                        if (cosmicData.COSMIC_EFFECT) cosmicItems['Effect'] = cosmicData.COSMIC_EFFECT;
                        // Compute frequencies if meta is available
                        let geneNameForDisplay = cosmicData.COSMIC_GENE || 'gene';
                        if (meta) {
                            const totalTumors = meta.total_samples_overall || 1;
                            const geneTumors = cosmicData.COSMIC_SAMPLES_WITH_GENE || 1;
                            const freqOverall = ((cosmicData.COSMIC_COUNT / totalTumors) * 100).toFixed(4);
                            const freqGene = ((cosmicData.COSMIC_COUNT / geneTumors) * 100).toFixed(2);
                            cosmicItems['Frequency (overall)'] = `${freqOverall}% of all tumors`;
                            cosmicItems[`Frequency in ${geneNameForDisplay}`] = `${freqGene}% of tumors with ${geneNameForDisplay} mutations`;
                        }
                        // Gene link to COSMIC analysis page if gene symbol is available
                        if (cosmicData.COSMIC_GENE) {
                            const encodedGene = encodeURIComponent(cosmicData.COSMIC_GENE);
                            const geneLink = `https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=${encodedGene}`;
                            cosmicItems['COSMIC Gene Page'] = { html: `<a href="${geneLink}" target="_blank" rel="noopener noreferrer">View analysis for ${cosmicData.COSMIC_GENE}</a>` };
                        }
                        // Site counts with per-type frequencies and gene-specific frequencies
                        if (cosmicData.COSMIC_SITE_COUNTS) {
                            const siteRows = [];
                            for (const [type, info] of Object.entries(cosmicData.COSMIC_SITE_COUNTS)) {
                                const count = info.count || 0;
                                const samplesWithGeneType = info.samples_with_gene_in_type || 1;
                                let rowText = `${type}: ${count} tumor${count === 1 ? '' : 's'}`;
                                if (meta && meta.total_samples_by_cancer_type && meta.total_samples_by_cancer_type[type]) {
                                    const typeTotal = meta.total_samples_by_cancer_type[type] || 1;
                                    const typeFreq = ((count / typeTotal) * 100).toFixed(2);
                                    const geneFreqType = ((count / samplesWithGeneType) * 100).toFixed(2);
                                    rowText += ` (${typeFreq}% of ${type}, ${geneFreqType}% with ${geneNameForDisplay})`;
                                }
                                siteRows.push(`<li>${rowText}</li>`);
                            }
                            cosmicItems['Site Counts'] = { html: `<ul>${siteRows.join('')}</ul>` };
                        }
                        if (meta && meta.last_updated) {
                            cosmicItems['Dataset Last Updated'] = meta.last_updated;
                        }
                        // Pre-computed frequencies and site breakdown are kept in
                        // details["COSMIC (Extended)"]; the raw counts + meta totals
                        // that produced them are omitted from the AI payload to avoid
                        // duplicating COSMIC data three ways.
                        detailsData.push({ title: 'COSMIC (Extended)', items: cosmicItems });
                    }
                } catch (cosmicErr) {
                    console.warn('COSMIC API error', cosmicErr);
                }
            }

            // Render details
            detailsContainer.innerHTML = '';
            detailsData.forEach(cat => {
                const detailsEl = document.createElement('details');
                const summaryEl = document.createElement('summary');
                summaryEl.textContent = cat.title;
                detailsEl.appendChild(summaryEl);
                const table = document.createElement('table');
                Object.keys(cat.items).forEach(key => {
                    const tr = document.createElement('tr');
                    const th = document.createElement('th');
                    th.textContent = key;
                    const td = document.createElement('td');
                    const value = cat.items[key];
                    if (value && typeof value === 'object' && Object.prototype.hasOwnProperty.call(value, 'html')) {
                        td.innerHTML = value.html;
                    } else {
                        td.textContent = String(value);
                    }
                    tr.appendChild(th);
                    tr.appendChild(td);
                    table.appendChild(tr);
                });
                detailsEl.appendChild(table);
                detailsContainer.appendChild(detailsEl);
            });

            // Build and display cards summarising key annotations
            const cardsContainer = document.getElementById('cardsContainer');
            cardsContainer.innerHTML = '';
            // Helper to extract field from summaryRows
            const getSummaryValue = (name) => {
                const r = summaryRows.find(row => row.name === name);
                return r ? String(r.value) : '';
            };
            // Extract gene(s), cDNA and protein
            let geneNames = getSummaryValue('Gene(s)');
            // Deduplicate gene names to avoid repeated display when multiple transcripts contribute the same gene.
            if (geneNames) {
                const uniq = Array.from(new Set(geneNames.split(',').map(g => g.trim()).filter(Boolean)));
                geneNames = uniq.join(', ');
            }
            // Extract cDNA (hgvsc) and protein (hgvsp) from dbNSFP. If arrays are present, select
            // the variant that best matches the protein change from the user's query. This helps
            // highlight the canonical transcript when multiple isoforms exist.
            // Build HTML strings listing all cDNA and protein notations, with the canonical entry in bold.
            let cDNAHTML = '';
            let transcriptsList = [];
            let protein = '';
            let proteinHTML = '';
            // If transcripts have been returned from the variant recoder and this is not a genomic-only query,
            // use them to build the cDNA/protein list. For genomic variants (original user input
            // was genomic) we ignore recoder results and derive cDNA/protein from the
            // MyVariant.info annotation instead. We must test the original user query rather
            // than the converted genomic variant (gVariant) to avoid discarding useful
            // transcripts for protein or cDNA inputs.
            // If the Ensembl variant recoder returned transcripts, build the cDNA/protein list
            // from these transcripts regardless of whether the input was originally genomic.  This
            // ensures that complex delins variants (which often lack MyVariant annotations) still
            // display transcript and protein information.  Previously, transcripts were ignored
            // when originalIsGenomic was true, leading to missing cDNA/protein details for
            // normalised g. variants.
            if (transcriptsFromRecoder && transcriptsFromRecoder.length > 0) {
                // Use transcript data from the variant recoder when available
                // Determine the canonical cDNA. Prefer RefSeq NM_ transcripts with the
                // lowest accession number (e.g. NM_004333 before NM_001354609). This
                // heuristic more consistently selects widely used canonical isoforms
                // across genes such as EGFR (NM_005228) and BRAF (NM_004333). If no
                // NM_ transcript is found, fall back to choosing the cDNA with the
                // largest numeric coordinate as before.
                // If the user query is a cDNA variant (e.g. contains ":c.") attempt to select
                // the transcript whose cDNA exactly matches the query's variant part. This helps
                // return the expected isoform for inputs like "MSH3 c.178_186del" or "PPM1D c.1518del".
                // Determine a canonical candidate cDNA based on user input and transcript accession/coordinate heuristics.
                let canonicalCandidate = null;
                // Track whether the canonical candidate was set directly from the user's cDNA query.
                let userSpecifiedCandidate = false;
                // Track the lowest RefSeq NM accession number encountered and the smallest positive cDNA coordinate
                // for that accession. These will be used to prefer commonly used isoforms when multiple NM_* transcripts
                // exist for the same gene.
                let minNmAcc = Infinity;
                let minCnumForMinNm = Infinity;
                // 1) If the user's query explicitly contains a cDNA (e.g. "c.1518del"), attempt to select the transcript
                // whose cDNA exactly matches that value (case-insensitive). This helps ensure that when a user enters
                // a specific cDNA (such as "MSH3 c.178_186del" or "PPM1D c.1518del"), the resulting summary will reflect
                // the same nomenclature rather than another isoform.
                {
                    const cdnaMatch = query && query.match(/c\.[^\s]+/i);
                    if (cdnaMatch) {
                        const userCdna = cdnaMatch[0].toLowerCase();
                        for (const t of transcriptsFromRecoder) {
                            if (t.cDNA && t.cDNA.toLowerCase() === userCdna) {
                                canonicalCandidate = t.cDNA;
                                userSpecifiedCandidate = true;
                                break;
                            }
                        }
                    }
                }
                // 2) Prefer RefSeq NM_* transcripts with positive cDNA coordinates. Among transcripts sharing the
                // same NM accession number, choose the transcript with the smallest cDNA coordinate (closest to
                // the N-terminus). This adjustment helps select the widely used canonical isoform when multiple
                // deletions or indels map to the same NM accession at different positions (e.g. PPM1D c.1518del vs c.1747del).
                for (const t of transcriptsFromRecoder) {
                    // Only consider NM_* accessions
                    if (!/^NM_/i.test(t.transcript)) continue;
                    const nmMatch = t.transcript.match(/^NM_0*([0-9]+)(?:\.|$)/i);
                    const cMatch = String(t.cDNA).match(/c\.\s*(-?\d+)/);
                    if (!nmMatch || !cMatch) continue;
                    const nmNum = parseInt(nmMatch[1], 10);
                    const cNum = parseInt(cMatch[1], 10);
                    // Only consider transcripts with positive coordinate (exclude upstream/promoter/intronic positions)
                    if (!userSpecifiedCandidate && !isNaN(nmNum) && !isNaN(cNum) && cNum > 0) {
                        // Skip intronic or upstream positions denoted by '+' or '-' following the coordinate (e.g. c.701+707del).
                        // Only consider exonic positions where the coordinate is a simple integer or range (contains '_' but no '+' or '-').
                        const cdnaAfter = String(t.cDNA).replace(/^c\./i, '');
                        const intronic = /[+-]/.test(cdnaAfter);
                        if (intronic) continue;
                        // If this NM accession is smaller than any seen before, or equal but has a smaller cDNA coordinate,
                        // update canonicalCandidate to this cDNA. This favours widely used isoforms and lower coordinate deletions.
                        if (nmNum < minNmAcc || (nmNum === minNmAcc && cNum < minCnumForMinNm)) {
                            minNmAcc = nmNum;
                            minCnumForMinNm = cNum;
                            canonicalCandidate = t.cDNA;
                        }
                    }
                }
                // 3) If no NM transcripts with positive coordinates were found, fall back to choosing the NM transcript
                // with the smallest accession number, regardless of coordinate sign. This maintains previous behaviour
                // for genes like EGFR where canonical isoforms have upstream coordinates but still follow NM accession ranking.
                if (!canonicalCandidate) {
                    for (const t of transcriptsFromRecoder) {
                        if (!/^NM_/i.test(t.transcript)) continue;
                        const nmMatch = t.transcript.match(/^NM_0*([0-9]+)(?:\.|$)/i);
                        if (!nmMatch) continue;
                        const nmNum = parseInt(nmMatch[1], 10);
                        if (!isNaN(nmNum) && (nmNum < minNmAcc)) {
                            minNmAcc = nmNum;
                            canonicalCandidate = t.cDNA;
                        }
                    }
                }
                // 4) As a last resort, choose the cDNA with the largest absolute numeric coordinate when no NM
                // accessions exist (e.g. genes without RefSeq transcripts). This preserves earlier heuristics.
                if (!canonicalCandidate) {
                    let best = transcriptsFromRecoder[0].cDNA;
                    let bestNum = -Infinity;
                    for (const t of transcriptsFromRecoder) {
                        const m = String(t.cDNA).match(/c\.\s*(-?\d+)/);
                        if (m) {
                            const num = parseInt(m[1], 10);
                            if (!isNaN(num) && num > bestNum) {
                                bestNum = num;
                                best = t.cDNA;
                            }
                        }
                    }
                    canonicalCandidate = best;
                }
                // Build transcriptsList with canonical flag
                transcriptsList = transcriptsFromRecoder.map((t) => {
                    return { ...t, canonical: t.cDNA === canonicalCandidate };
                });
                cDNAHTML = transcriptsList
                    .map((t) => (t.canonical ? `<strong>${t.cDNA}</strong>` : t.cDNA))
                    .join(', ');
                const canonicalEntry = transcriptsList.find((t) => t.canonical);
                // Determine canonical protein: if canonical entry has a protein, use it. Otherwise pick first non-empty protein.
                let canonicalProt = canonicalEntry ? canonicalEntry.protein : '';
                // If the canonical cDNA refers to an intronic or upstream position (contains '+' or '-')
                // then do not display a protein change. Intronic variants do not result in an amino acid
                // alteration and should leave the p. field blank.
                if (canonicalEntry && canonicalEntry.cDNA) {
                    const cdnaNoPrefix = canonicalEntry.cDNA.replace(/^c\./i, '');
                    if (/[+-]/.test(cdnaNoPrefix)) {
                        canonicalProt = '';
                    }
                }
                if (!canonicalProt) {
                    const firstProtEntry = transcriptsList.find(t => t.protein);
                    canonicalProt = firstProtEntry ? firstProtEntry.protein : '';
                }
                // Build protein HTML, bolding the canonical protein. If proteins are missing, skip bolding.
                proteinHTML = transcriptsList
                    .map((t) => {
                        if (!t.protein) return '';
                        return t.protein === canonicalProt ? `<strong>${t.protein}</strong>` : t.protein;
                    })
                    .filter(Boolean)
                    .join(', ');
                protein = canonicalProt || '';
            } else {
                // Fallback: derive from dbNSFP if recoder data is unavailable
                if (annotation.dbnsfp && annotation.dbnsfp.hgvsc) {
                    const hgvscList = Array.isArray(annotation.dbnsfp.hgvsc)
                        ? annotation.dbnsfp.hgvsc
                        : [annotation.dbnsfp.hgvsc];
                    let canonicalCandidate = null;
                    /*
                     * Determine the canonical cDNA among multiple dbNSFP hgvsc entries.  The goal
                     * is to pick the transcript corresponding to the most widely used isoform
                     * rather than arbitrarily choosing the first entry.  The selection logic is:
                     *   1. If all entries contain a numeric coordinate (c.<number>), compute
                     *      these coordinates and select the entry whose coordinate is the median
                     *      of the set.  This heuristic tends to favour canonical isoforms (e.g.
                     *      BRAF c.1799T>A) when alternative isoforms with smaller (e.g. 620) or
                     *      larger (e.g. 1919) coordinates are present.
                     *   2. If multiple entries share the median coordinate, prefer those with
                     *      RefSeq NM_ accession prefixes and the lowest accession number.
                     *   3. If numeric parsing fails for any entry, fall back to preferring NM_
                     *      transcripts with the lowest accession number overall.
                     *   4. If no NM_ accession is present, fall back to the entry with the
                     *      largest numeric coordinate (same as before).
                     */
                    // Parse numeric coordinates from each hgvsc entry.  Use an array equal in
                    // length to hgvscList where non-numeric entries remain null.  This
                    // enables robust handling when some transcripts have intronic or
                    // complex notation (e.g. c.1906-7T>A) without failing the median
                    // selection logic entirely.
                    const coords = new Array(hgvscList.length).fill(null);
                    for (let i = 0; i < hgvscList.length; i++) {
                        const h = hgvscList[i];
                        const m = String(h).match(/c\.\s*(-?\d+)/i);
                        if (m) {
                            const num = parseInt(m[1], 10);
                            if (!isNaN(num)) coords[i] = num;
                        }
                    }
                    // Extract only the numeric coordinates (filter out null) for median
                    const validCoords = coords.filter((c) => c !== null);
                    if (validCoords.length >= 2) {
                        // Compute median of numeric coordinates
                        const sorted = [...validCoords].sort((a, b) => a - b);
                        const medianVal = sorted[Math.floor(sorted.length / 2)];
                        // Find indexes with coordinate equal to medianVal (ties allowed)
                        const candidateIndexes = [];
                        for (let i = 0; i < coords.length; i++) {
                            if (coords[i] === medianVal) {
                                candidateIndexes.push(i);
                            }
                        }
                        if (candidateIndexes.length > 0) {
                            // Among candidates, prefer NM_ transcripts with lowest accession
                            let bestIndex = candidateIndexes[0];
                            let bestAcc = Infinity;
                            for (const idx of candidateIndexes) {
                                const h = hgvscList[idx];
                                const txId = String(h).split(':')[0];
                                if (/^NM_/i.test(txId)) {
                                    const mAcc = txId.match(/^NM_0*([0-9]+)(?:\.|$)/i);
                                    if (mAcc) {
                                        const numAcc = parseInt(mAcc[1], 10);
                                        if (!isNaN(numAcc) && numAcc < bestAcc) {
                                            bestAcc = numAcc;
                                            bestIndex = idx;
                                        }
                                    }
                                }
                            }
                            canonicalCandidate = hgvscList[bestIndex];
                        }
                    }
                    if (!canonicalCandidate) {
                        // Fallback: prefer NM_ transcripts with lowest accession overall
                        let minNmAcc = Infinity;
                        for (const h of hgvscList) {
                            const parts = String(h).split(':');
                            const txId = parts[0];
                            if (/^NM_/i.test(txId)) {
                                const m = txId.match(/^NM_0*([0-9]+)(?:\.|$)/i);
                                if (m) {
                                    const num = parseInt(m[1], 10);
                                    if (!isNaN(num) && num < minNmAcc) {
                                        minNmAcc = num;
                                        canonicalCandidate = h;
                                    }
                                }
                            }
                        }
                    }
                    if (!canonicalCandidate) {
                        // If no NM transcript found, choose the entry with the largest numeric coordinate.
                        let best = hgvscList[0];
                        let bestNum = -Infinity;
                        for (const h of hgvscList) {
                            const m = String(h).match(/c\.(-?\d+)/i);
                            if (m) {
                                const num = parseInt(m[1], 10);
                                if (!isNaN(num) && num > bestNum) {
                                    bestNum = num;
                                    best = h;
                                }
                            }
                        }
                        canonicalCandidate = best;
                    }
                    // Build a highlighted cDNA HTML string
                    cDNAHTML = hgvscList
                        .map((h) => (h === canonicalCandidate ? `<strong>${h}</strong>` : h))
                        .join(', ');
                    // Build a list of transcript mappings (transcript ID -> cDNA, protein). We'll align hgvsc and hgvsp by index.
                    const hgvspListLocal = (annotation.dbnsfp && annotation.dbnsfp.hgvsp) ? (Array.isArray(annotation.dbnsfp.hgvsp) ? annotation.dbnsfp.hgvsp : [annotation.dbnsfp.hgvsp]) : [];
                    for (let i = 0; i < hgvscList.length; i++) {
                        const sc = hgvscList[i];
                        if (!sc) continue;
                        const parts = String(sc).split(':');
                        const transcriptId = parts[0];
                        const cpart = parts.slice(1).join(':');
                        let ppart = '';
                        if (hgvspListLocal[i]) {
                            const pparts = String(hgvspListLocal[i]).split(':');
                            ppart = pparts.slice(1).join(':');
                        }
                        transcriptsList.push({ transcript: transcriptId, cDNA: cpart, protein: ppart, canonical: sc === canonicalCandidate });
                    }
                }
                if (annotation.dbnsfp && annotation.dbnsfp.hgvsp) {
                    const hgvspList = Array.isArray(annotation.dbnsfp.hgvsp) ? annotation.dbnsfp.hgvsp : [annotation.dbnsfp.hgvsp];
                    let canonicalProt = hgvspList[0];
                    // Choose canonical protein: match triple-coded targetProtGlobal if possible
                    if (targetProtGlobal) {
                        for (const p of hgvspList) {
                            const triple = p.replace(/^.*p\./i, '').toUpperCase();
                            if (triple === targetProtGlobal) {
                                canonicalProt = p;
                                break;
                            }
                        }
                    }
                    proteinHTML = hgvspList
                        .map((p) => (p === canonicalProt ? `<strong>${p}</strong>` : p))
                        .join(', ');
                    // Also expose the canonical protein string in a variable used later (e.g. for search links)
                    protein = canonicalProt;
                }
            // If no transcripts were gathered from dbNSFP and snpEff annotations are present,
            // attempt to extract cDNA and protein information from snpEff. This helps capture
            // upstream/promoter variants where dbNSFP does not provide hgvsc/hgvsp (e.g. TERT
            // promoter mutations).
            if (transcriptsList.length === 0 && annotation.snpeff && annotation.snpeff.ann) {
                const annList = Array.isArray(annotation.snpeff.ann) ? annotation.snpeff.ann : [annotation.snpeff.ann];
                const snpEffList = [];
                for (const ann of annList) {
                    if (ann && ann.hgvs_c) {
                        // Determine a transcript ID; snpEff feature_id may include version or be absent.
                        const txId = ann.feature_id || '';
                        const cDNAVal = ann.hgvs_c;
                        // snpEff may provide hgvs_p, otherwise leave empty
                        const protVal = ann.hgvs_p || '';
                        snpEffList.push({ transcript: txId, cDNA: cDNAVal, protein: protVal });
                    }
                }
                if (snpEffList.length > 0) {
                    // Determine canonical cDNA among snpEff annotations. Prefer NM accessions with
                    // the lowest accession number; otherwise pick the first entry.
                    let canonicalCandidateEff = null;
                    let minNmEff = Infinity;
                    for (const t of snpEffList) {
                        if (/^NM_/i.test(t.transcript)) {
                            const m = t.transcript.match(/^NM_0*([0-9]+)(?:\.|$)/i);
                            if (m) {
                                const num = parseInt(m[1], 10);
                                if (!isNaN(num) && num < minNmEff) {
                                    minNmEff = num;
                                    canonicalCandidateEff = t.cDNA;
                                }
                            }
                        }
                    }
                    if (!canonicalCandidateEff) {
                        canonicalCandidateEff = snpEffList[0].cDNA;
                    }
                    transcriptsList = snpEffList.map((t) => ({ ...t, canonical: t.cDNA === canonicalCandidateEff }));
                    // Build cDNA HTML representation
                    cDNAHTML = transcriptsList
                        .map((t) => (t.canonical ? `<strong>${t.cDNA}</strong>` : t.cDNA))
                        .join(', ');
                    // Determine canonical protein: use the protein from the canonical entry if available,
                    // otherwise use the first available protein string.
                    let canonicalProtEff = '';
                    const canonEntry = transcriptsList.find(t => t.canonical);
                    if (canonEntry) {
                        // Only use the protein if the canonical cDNA is not intronic/upstream (contains '+' or '-')
                        const cdnaNoPrefixEff = canonEntry.cDNA.replace(/^c\./i, '');
                        const intronicEff = /[+-]/.test(cdnaNoPrefixEff);
                        if (!intronicEff && canonEntry.protein) {
                            canonicalProtEff = canonEntry.protein;
                        }
                    }
                    if (!canonicalProtEff) {
                        const firstProt = transcriptsList.find(t => t.protein);
                        canonicalProtEff = firstProt ? firstProt.protein : '';
                    }
                    // Build protein HTML (bold the canonical protein) if proteins exist
                    proteinHTML = transcriptsList
                        .map((t) => {
                            if (!t.protein) return '';
                            return t.protein === canonicalProtEff ? `<strong>${t.protein}</strong>` : t.protein;
                        })
                        .filter(Boolean)
                        .join(', ');
                    protein = canonicalProtEff;
                }
            }
            }
            // Additional fallback: if still no transcripts and root-level hgvsc/hgvsp fields exist on the annotation
            // object, use them directly. MyVariant.info sometimes places HGVS coding and protein changes at the
            // top level (e.g. annotation.hgvsc and annotation.hgvsp). In such cases we construct a minimal
            // transcript list and choose a canonical entry (the first entry or one matching the user's protein query).
            if (transcriptsList.length === 0 && annotation && (annotation.hgvsc || annotation.hgvsp)) {
                const hgvscList = annotation.hgvsc ? (Array.isArray(annotation.hgvsc) ? annotation.hgvsc : [annotation.hgvsc]) : [];
                const hgvspList = annotation.hgvsp ? (Array.isArray(annotation.hgvsp) ? annotation.hgvsp : [annotation.hgvsp]) : [];
                let canonicalIdx = 0;
                if (targetProtGlobal && hgvspList.length > 0) {
                    // When a target protein (e.g. from the user's query) is available, select
                    // the matching hgvsp entry as the canonical index. Compare both
                    // three‑letter (triple) and single‑letter conversions for robustness.
                    for (let i = 0; i < hgvspList.length; i++) {
                        const p = hgvspList[i];
                        if (!p) continue;
                        const prot = String(p).replace(/^p\./i, '').toUpperCase();
                        const triple = prot.replace(/[^A-Z0-9]/g, '');
                        const single = tripleToSingle(triple);
                        if (triple === targetProtGlobal || (single && triple === single)) {
                            canonicalIdx = i;
                            break;
                        }
                    }
                } else if (hgvscList.length > 1) {
                    /*
                     * For genomic variants where multiple cDNA notations exist but no specific
                     * protein target was provided, choose the cDNA entry whose numeric
                     * coordinate is the median of all available coordinates.  This heuristic
                     * tends to select the widely used isoform (e.g. c.1799T>A for BRAF V600E)
                     * when other isoforms with smaller or larger coordinates are also present
                     * (e.g. c.620T>A, c.1919T>A).  If coordinate parsing fails for any
                     * reason, the canonicalIdx remains 0 (defaulting to the first entry).
                     */
                    const coords = [];
                    for (const h of hgvscList) {
                        const m = String(h).match(/c\.\s*(-?\d+)/);
                        if (m) {
                            const num = parseInt(m[1], 10);
                            if (!isNaN(num)) coords.push(num);
                        }
                    }
                    if (coords.length === hgvscList.length) {
                        // Compute median coordinate
                        const sorted = [...coords].sort((a, b) => a - b);
                        const medianVal = sorted[Math.floor(sorted.length / 2)];
                        // Select the hgvsc entry whose coordinate is closest to the median
                        let bestDiff = Infinity;
                        let bestIdx = 0;
                        for (let i = 0; i < coords.length; i++) {
                            const diff = Math.abs(coords[i] - medianVal);
                            if (diff < bestDiff) {
                                bestDiff = diff;
                                bestIdx = i;
                            }
                        }
                        canonicalIdx = bestIdx;
                    }
                }
                cDNAHTML = hgvscList
                    .map((h, i) => (i === canonicalIdx ? `<strong>${h}</strong>` : h))
                    .join(', ');
                proteinHTML = hgvspList
                    .map((p, i) => {
                        if (!p) return '';
                        return i === canonicalIdx ? `<strong>${p}</strong>` : p;
                    })
                    .filter(Boolean)
                    .join(', ');
                protein = hgvspList[canonicalIdx] || '';
                // Build a minimal transcriptsList without transcript IDs (empty string)
                transcriptsList = hgvscList.map((h, i) => {
                    return { transcript: '', cDNA: h, protein: hgvspList[i] || '', canonical: i === canonicalIdx };
                });
            }
            // Before computing effect/consequence and rendering the variant card, unify
            // transcripts for genomic inputs. When the original user query was a genomic
            // variant (g.), transcriptsFromRecoder is unused and transcriptsList may
            // contain entries derived from dbNSFP, snpEff or root-level fields. To
            // ensure the canonical transcript is chosen using a unified scoring
            // heuristic across all sources, we rebuild the transcript list from
            // annotation here. This step preserves existing transcriptsList for
            // non-genomic queries (where transcriptsFromRecoder is used).
            if (originalIsGenomic) {
                // When the variant_recoder provided transcripts for a genomic input, we keep
                // those transcripts rather than rebuilding from the MyVariant annotation.  This
                // ensures complex delins variants retain the recoder‑derived cDNA/protein
                // information.  Only unify transcripts from the annotation when no recoder
                // transcripts are available.
                if (!transcriptsFromRecoder || transcriptsFromRecoder.length === 0) {
                    const unifyRes = buildCanonicalFromAnnotation(annotation, targetProtGlobal);
                    // Only overwrite when a non-empty list is returned. Otherwise leave
                    // transcriptsList unchanged.
                    if (unifyRes.transcriptsList && unifyRes.transcriptsList.length > 0) {
                        transcriptsList = unifyRes.transcriptsList;
                        cDNAHTML = unifyRes.cDNAHTML;
                        proteinHTML = unifyRes.proteinHTML;
                        protein = unifyRes.canonicalProtein;
                    }
                }
            }
            // Effect / consequence
            let effect = getSummaryValue('Consequence');
            if (!effect) effect = getSummaryValue('Variant Type');
            // Card: Variant info
            {
                const card = document.createElement('div');
                card.className = 'card';
                const title = document.createElement('h3');
                title.textContent = 'Variant';
                applyCardTheme(card, 'Variant');
                card.appendChild(title);
                const content = document.createElement('div');
                content.className = 'card-content';
                const makeLine = (label, value) => {
                    const span = document.createElement('span');
                    span.innerHTML = `<strong>${label}:</strong> ${value || 'N/A'}`;
                    return span;
                };
                content.appendChild(makeLine('g.', gVariant));
                content.appendChild(makeLine('Gene', geneNames));
                // Determine canonical cDNA and protein from transcriptsList. If not available, fall back to first values.
                let canonicalEntryForDisplay = null;
                if (transcriptsList && transcriptsList.length > 0) {
                    canonicalEntryForDisplay = transcriptsList.find(t => t.canonical) || transcriptsList[0];
                }
                /*
                 * Determine the canonical cDNA (c.) and protein (p.) values to display in the
                 * summary card.  In earlier builds, certain genomic inputs resulted in
                 * populated transcripts lists yet the displayed c. and p. fields were
                 * erroneously reported as "N/A".  This was due to a combination of
                 * canonical entries lacking a protein (intronic variants) and fallback
                 * variables (cDNAHTML/protein) not being consulted when transcripts were
                 * available.  The logic below addresses this by:
                 *   1. Extracting the canonical entry from transcriptsList when present.
                 *      If the canonical entry lacks a protein but another transcript has
                 *      one, that protein is used instead.
                 *   2. Falling back to root-level annotation.hgvsc/hgvsp arrays when no
                 *      transcripts are available or the canonical values remain empty.
                 *   3. Leaving the field blank only when no value can be resolved, allowing
                 *      the downstream makeLine helper to insert "N/A".
                 */
                let canonicalCVal = '';
                let canonicalProtVal = '';
                if (transcriptsList && transcriptsList.length > 0) {
                    const canonicalEntryForDisplay = transcriptsList.find(t => t.canonical) || transcriptsList[0];
                    if (canonicalEntryForDisplay) {
                        canonicalCVal = canonicalEntryForDisplay.cDNA || '';
                        canonicalProtVal = canonicalEntryForDisplay.protein || '';
                    }
                    // If canonical protein is empty, but another transcript has a protein, use that
                    if (!canonicalProtVal) {
                        const firstProtEntry = transcriptsList.find(t => t.protein);
                        if (firstProtEntry) canonicalProtVal = firstProtEntry.protein || '';
                    }
                }
                // Fall back to the first hgvsc/hgvsp at the root level of the annotation when still empty
                if (!canonicalCVal && annotation) {
                    if (annotation.hgvsc) {
                        canonicalCVal = Array.isArray(annotation.hgvsc) ? annotation.hgvsc[0] : annotation.hgvsc;
                    }
                }
                if (!canonicalProtVal && annotation) {
                    if (annotation.hgvsp) {
                        canonicalProtVal = Array.isArray(annotation.hgvsp) ? annotation.hgvsp[0] : annotation.hgvsp;
                    }
                }
                // As a further fallback, use values derived from cDNAHTML/protein variables when
                // canonical values are still empty.  These variables capture the first entry from
                // any formatted lists built earlier (e.g. root‑level hgvsc/hgvsp arrays) and
                // provide a sensible default when annotation properties are unavailable.
                if (!canonicalCVal && cDNAHTML) {
                    const cdClean = cDNAHTML.replace(/<[^>]+>/g, '').split(',')[0].trim();
                    if (cdClean) canonicalCVal = cdClean;
                }
                if (!canonicalProtVal && protein) {
                    canonicalProtVal = protein;
                }
                content.appendChild(makeLine('c.', canonicalCVal));
                content.appendChild(makeLine('p.', formatProteinDisplayWithSingleLetter(canonicalProtVal)));
                content.appendChild(makeLine('Effect', effect));
                const ucscUrl = buildUcscHg19Url(rawInput, gVariant, annotation);
                if (ucscUrl) {
                    content.appendChild(makeLine('UCSC (hg19)', `<a href="${ucscUrl}" target="_blank" rel="noopener noreferrer">Zoom to region</a>`));
                }
                // Append list of transcripts showing cDNA and protein for each transcript in a collapsible details element
                if (transcriptsList && transcriptsList.length > 1) {
                    const detailsEl = document.createElement('details');
                    const summaryEl = document.createElement('summary');
                    summaryEl.textContent = 'Other transcripts';
                    detailsEl.appendChild(summaryEl);
                    const ul = document.createElement('ul');
                    ul.style.marginTop = '0.25em';
                    ul.style.paddingLeft = '1.2em';
                    transcriptsList.forEach((t) => {
                        const li = document.createElement('li');
                        let inner = `${t.transcript}: ${t.cDNA}`;
                        if (t.protein) inner += `, ${t.protein}`;
                        if (t.canonical) {
                            li.innerHTML = `<strong>${inner}</strong>`;
                        } else {
                            li.textContent = inner;
                        }
                        ul.appendChild(li);
                    });
                    detailsEl.appendChild(ul);
                    content.appendChild(detailsEl);
                }
                card.appendChild(content);
                cardsContainer.appendChild(card);
            }
            // Card: ClinVar
            {
                // Build a more comprehensive ClinVar card using data from annotation.clinvar
                const clin = annotation.clinvar;
                // Determine the default ClinVar link: use variant_id when available, otherwise build a search link.
                let clinVarLink = '';
                let variantId = '';
                if (clin && clin.variant_id) {
                    variantId = String(clin.variant_id);
                    // numeric variant_id corresponds to a ClinVar Variation page
                    if (/^\d+$/.test(variantId)) {
                        clinVarLink = `https://www.ncbi.nlm.nih.gov/clinvar/variation/${variantId}/`;
                    } else {
                        clinVarLink = `https://www.ncbi.nlm.nih.gov/clinvar/${variantId}`;
                    }
                } else {
                    // Build a search term using gene and protein when variant_id is unavailable
                    let searchTerm = '';
                    const firstGene = geneNames ? geneNames.split(',')[0].trim() : '';
                    let protSingle = '';
                    if (targetProtGlobal) {
                        const tmp = tripleToSingle(targetProtGlobal);
                        if (tmp) protSingle = tmp;
                    }
                    if (!protSingle && protein) {
                        const tripleMatch = String(protein).match(/([A-Za-z]{3})(\d+)([A-Za-z]{3})/);
                        if (tripleMatch) {
                            const triple = (tripleMatch[1] + tripleMatch[2] + tripleMatch[3]).toUpperCase();
                            const tmp = tripleToSingle(triple);
                            if (tmp) protSingle = tmp;
                        } else {
                            const m2 = String(protein).match(/([A-Za-z])(\d+)([A-Za-z])/);
                            if (m2) protSingle = m2[1].toUpperCase() + m2[2] + m2[3].toUpperCase();
                        }
                    }
                    if (firstGene) {
                        searchTerm = firstGene;
                        if (protSingle) searchTerm += ` ${protSingle}`;
                    }
                    if (!searchTerm) searchTerm = gVariant;
                    clinVarLink = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(searchTerm)}`;
                }
                // Summarize significance categories and conditions from RCV entries
                let sigSummary = [];
                let conditionsList = [];
                let rcvDetails = [];
                if (clin && clin.rcv) {
                    const rcvArr = Array.isArray(clin.rcv) ? clin.rcv : [clin.rcv];
                    const sigCount = {};
                    const condSet = new Set();
                    rcvArr.forEach(rc => {
                        let cs = rc.clinical_significance;
                        let sigStr = null;
                        if (cs) {
                            if (typeof cs === 'object' && cs.description) sigStr = cs.description;
                            else sigStr = cs;
                            if (sigStr) {
                                sigCount[sigStr] = (sigCount[sigStr] || 0) + 1;
                            }
                        }
                        if (rc.conditions && rc.conditions.name) {
                            condSet.add(rc.conditions.name);
                        }
                        // Collect details per accession for optional details view
                        const acc = rc.accession || '';
                        const condName = rc.conditions?.name || '';
                        const review = rc.review_status || '';
                        rcvDetails.push({ accession: acc, significance: sigStr || 'N/A', condition: condName, review });
                    });
                    sigSummary = Object.entries(sigCount).map(([k, v]) => `${k} (${v})`);
                    conditionsList = Array.from(condSet);
                }
                // Build ClinVar card
                const card = document.createElement('div');
                card.className = 'card';
                const title = document.createElement('h3');
                title.textContent = 'ClinVar';
                applyCardTheme(card, 'ClinVar');
                card.appendChild(title);
                const content = document.createElement('div');
                content.className = 'card-content';
                // Variation ID line if present
                if (variantId) {
                    const spanVar = document.createElement('div');
                    spanVar.innerHTML = `<strong>Variation ID:</strong> ${variantId}`;
                    content.appendChild(spanVar);
                }
                // Significance summary
                if (sigSummary.length > 0) {
                    const spanSig = document.createElement('div');
                    spanSig.innerHTML = `<strong>Clinical significance:</strong> ${sigSummary.join('; ')}`;
                    content.appendChild(spanSig);
                } else {
                    const spanSig = document.createElement('div');
                    spanSig.innerHTML = `<strong>Clinical significance:</strong> N/A`;
                    content.appendChild(spanSig);
                }
                // Conditions summary (show up to 3, rest collapsed)
                if (conditionsList.length > 0) {
                    const spanCond = document.createElement('div');
                    const displayConds = conditionsList.slice(0, 3).join(', ');
                    spanCond.innerHTML = `<strong>Conditions:</strong> ${displayConds}${conditionsList.length > 3 ? '…' : ''}`;
                    content.appendChild(spanCond);
                }
                // Somatic / oncogenicity data fetched directly from ClinVar VCV record.
                if (variantId && /^\d+$/.test(variantId)) {
                    try {
                        const cvData = await fetchClinvarVariant(variantId);
                        aiReviewExtras.clinvar_variant_record = cvData;
                        if (cvData) {
                            if (cvData.somatic && cvData.somatic.description) {
                                const somaticDiv = document.createElement('div');
                                somaticDiv.style.marginTop = '0.25rem';
                                somaticDiv.innerHTML = `<strong>Somatic clinical impact:</strong> ${cvData.somatic.description}`;
                                content.appendChild(somaticDiv);
                            }
                            if (cvData.oncogenicity && cvData.oncogenicity.description) {
                                const oncDiv = document.createElement('div');
                                oncDiv.innerHTML = `<strong>Oncogenicity:</strong> ${cvData.oncogenicity.description}`;
                                content.appendChild(oncDiv);
                            }
                            if (cvData.somaticConditions && cvData.somaticConditions.length > 0) {
                                const scDet = document.createElement('details');
                                const scSum = document.createElement('summary');
                                scSum.textContent = `Somatic conditions (${cvData.somaticConditions.length})`;
                                scDet.appendChild(scSum);
                                const scUl = document.createElement('ul');
                                scUl.style.cssText = 'margin-top:0.4rem;font-size:0.82rem;padding-left:1.2rem;line-height:1.6;';
                                cvData.somaticConditions.forEach((sc) => {
                                    const li = document.createElement('li');
                                    const parts = [];
                                    if (sc.condition) parts.push(`<strong>${sc.condition}</strong>`);
                                    // Combine tier + assertion type + clinical significance into one readable string.
                                    const impactParts = [sc.tier, sc.assertionType, sc.clinSig].filter(Boolean);
                                    if (impactParts.length) parts.push(impactParts.join(' — '));
                                    li.innerHTML = parts.join(': ');
                                    scUl.appendChild(li);
                                });
                                scDet.appendChild(scUl);
                                content.appendChild(scDet);
                            }
                        }
                    } catch (e) {
                        console.warn('ClinVar variant fetch failed', e);
                    }
                }
                // Link to ClinVar
                const linkEl = document.createElement('a');
                linkEl.href = clinVarLink;
                linkEl.target = '_blank';
                linkEl.rel = 'noopener noreferrer';
                linkEl.textContent = 'View on ClinVar';
                content.appendChild(linkEl);
                // Nearby ClinVar variants (±30 bp, GRCh37).
                const tuple = buildSpliceAiLookupTuple(rawInput, gVariant);
                if (tuple && tuple.chrom && tuple.pos) {
                    const chr = tuple.chrom.replace(/^chr/i, '');
                    const posNum = Number(tuple.pos);
                    const regionLink = document.createElement('a');
                    regionLink.href = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(`${chr}[Chromosome] AND ${Math.max(1, posNum - 30)}:${posNum + 30}[Base Position for Assembly GRCh37]`)}`;
                    regionLink.target = '_blank';
                    regionLink.rel = 'noopener noreferrer';
                    regionLink.style.display = 'block';
                    regionLink.textContent = 'Search region in ClinVar';
                    content.appendChild(regionLink);
                    try {
                        const { variants: nearby, total: nearbyTotal } = await fetchClinvarRegionVariants(chr, posNum, 30);
                        aiReviewExtras.nearby_clinvar_variants = nearby;

                        // Detect gene strand from snpEff annotation to orient plot 5'→3'.
                        const snpEffAnns = annotation?.snpeff?.ann
                            ? (Array.isArray(annotation.snpeff.ann) ? annotation.snpeff.ann : [annotation.snpeff.ann])
                            : [];
                        const plotMinusStrand = snpEffAnns.some(a => a.strand === '-' || a.strand === -1);

                        // Extract canonical c. position (e.g. c.454 → 454) from the snpEff annotation
                        // that has a protein change; used to compute transcript-agnostic protein positions.
                        const canonicalAnn = snpEffAnns.find(a => a.hgvs_p);
                        const queryHgvsP = canonicalAnn?.hgvs_p || '';
                        const queryProteinPos = parseProteinPos(queryHgvsP);
                        const queryC = parseCdnaCoordinate(canonicalAnn?.hgvs_c || '');

                        if (nearby.length > 0) {
                            const nearbyTitle = document.createElement('div');
                            nearbyTitle.style.cssText = 'font-size:0.86rem;font-weight:600;margin-top:0.5rem;';
                            const truncated = nearbyTotal > nearby.length ? ` of ${nearbyTotal} total` : '';
                            nearbyTitle.textContent = `Nearby ClinVar variants (±30 bp): ${nearby.length}${truncated}`;
                            content.appendChild(nearbyTitle);

                            // Genomic (g.) lollipop plot — display ±10 bp even though data covers ±30
                            const plotWrap = document.createElement('div');
                            plotWrap.style.cssText = 'margin:6px 0 4px;';
                            plotWrap.appendChild(buildLollipopPlot(nearby, posNum, plotMinusStrand, { range: 10 }));
                            content.appendChild(plotWrap);

                            // Protein (p.) lollipop plot — uses genomic-offset positions for alignment
                            if (queryProteinPos !== null && queryC !== null) {
                                const protPlot = buildProteinLollipopPlot(
                                    nearby, queryProteinPos,
                                    { queryGenomicPos: posNum, queryC, minusStrand: plotMinusStrand }
                                );
                                if (protPlot) {
                                    const protPlotWrap = document.createElement('div');
                                    protPlotWrap.style.cssText = 'margin:10px 0 4px;';
                                    const protPlotTitle = document.createElement('div');
                                    protPlotTitle.style.cssText = 'font-size:0.82rem;font-weight:600;margin-bottom:2px;';
                                    protPlotTitle.textContent = 'Protein position (p.)';
                                    protPlotWrap.appendChild(protPlotTitle);
                                    protPlotWrap.appendChild(protPlot);
                                    content.appendChild(protPlotWrap);
                                }
                            }

                            // Condensed variant list (collapsed by default)
                            const listDet = document.createElement('details');
                            listDet.style.marginTop = '4px';
                            const listSum = document.createElement('summary');
                            listSum.style.cssText = 'font-size:0.82rem;font-weight:normal;';
                            listSum.textContent = 'Show variant list';
                            listDet.appendChild(listSum);
                            const ul = document.createElement('ul');
                            ul.style.cssText = 'margin-top:0.4rem;font-size:0.8rem;padding-left:1.2rem;line-height:1.5;';
                            nearby.forEach((v) => {
                                const li = document.createElement('li');
                                const color = getPathogenicityColor(v.germline, v);
                                const sigSpan = `<span style="color:${color};font-weight:600">${v.germline || 'Unknown'}</span>`;
                                const posInfo = v.pos ? ` · pos ${v.pos}` : '';
                                const cvLink = `<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.id}/" target="_blank" rel="noopener noreferrer">${v.id}</a>`;
                                li.innerHTML = `${sigSpan}${posInfo} — ${v.title || cvLink}`;
                                ul.appendChild(li);
                            });
                            listDet.appendChild(ul);
                            content.appendChild(listDet);
                        }
                    } catch (e) {
                        console.warn('ClinVar regional query failed', e);
                        const errNote = document.createElement('div');
                        errNote.style.cssText = 'font-size:0.82rem;color:#9ca3af;margin-top:4px;';
                        errNote.textContent = 'Region data unavailable.';
                        content.appendChild(errNote);
                    }
                }
                // Details: list RCV entries if more than one
                if (rcvDetails.length > 0) {
                    const detailsEl = document.createElement('details');
                    const summaryEl = document.createElement('summary');
                    summaryEl.textContent = 'Show RCV details';
                    detailsEl.appendChild(summaryEl);
                    const ul = document.createElement('ul');
                    ul.style.marginTop = '0.5rem';
                    rcvDetails.forEach((rc) => {
                        const li = document.createElement('li');
                        const parts = [];
                        if (rc.accession) parts.push(rc.accession);
                        if (rc.significance) parts.push(rc.significance);
                        if (rc.condition) parts.push(rc.condition);
                        if (rc.review) parts.push(rc.review);
                        li.textContent = parts.join(' | ');
                        ul.appendChild(li);
                    });
                    detailsEl.appendChild(ul);
                    content.appendChild(detailsEl);
                }
                card.appendChild(content);
                cardsContainer.appendChild(card);
            }
            // Card: CIViC (always shown)
            {
                const entries = Array.isArray(annotation.cgi) ? annotation.cgi : (Array.isArray(annotation.civic) ? annotation.civic : []);
                const legacy = (annotation.civic && typeof annotation.civic === 'object' && !Array.isArray(annotation.civic)) ? annotation.civic : null;
                const hasMyVariantData = entries.length > 0 || legacy;

                const card = document.createElement('div');
                card.className = 'card';
                const civicTitle = document.createElement('h3');
                civicTitle.textContent = 'CIViC';
                applyCardTheme(card, 'CIViC');
                card.appendChild(civicTitle);
                const content = document.createElement('div');
                content.className = 'card-content';

                const addLine = (label, value) => {
                    if (value === null || value === undefined || value === '') return;
                    const div = document.createElement('div');
                    div.style.marginBottom = '0.2rem';
                    div.innerHTML = `<strong>${label}:</strong> ${value}`;
                    content.appendChild(div);
                };

                // Determine gene name and protein change for links and API call
                const civicGene = legacy?.gene?.name || (geneNames ? geneNames.split(',')[0].trim() : '');
                const rawProtein = String(protein || '').replace(/<[^>]+>/g, '');
                const firstProtein = rawProtein.split(',')[0].trim();
                const annotationProtein = firstProtein.includes(':') ? firstProtein.split(':').slice(1).join(':').trim() : firstProtein;
                const bestCivicEntry = findBestCivicEntryForProtein(entries, annotationProtein);
                const civicProtein = legacy?.variant?.name || bestCivicEntry?.protein_change || annotationProtein;

                // Variant link: use a known variant ID immediately if available. Otherwise, keep
                // it hidden until the CIViC API callback resolves an exact variant page.
                const encodedCivicGene = encodeURIComponent(civicGene || '');
                const civicVariantSearchTerm = civicProtein || '';
                const encodedCivicVariantQuery = civicVariantSearchTerm
                    ? `${encodedCivicGene}+${encodeURIComponent(civicVariantSearchTerm)}`
                    : encodedCivicGene;
                const linksDiv = document.createElement('div');
                linksDiv.style.marginBottom = '0.4rem';
                let variantLinkEl = null;
                let variantSepNode = null;
                let geneLinkEl = null;
                const legacyVariantId = legacy?.variant_id;
                if (civicGene) {
                    variantLinkEl = document.createElement('a');
                    if (legacyVariantId) {
                        variantLinkEl.href = `https://civicdb.org/variants/${legacyVariantId}/summary`;
                    } else {
                        variantLinkEl.href = `https://civicdb.org/search?query=${encodedCivicVariantQuery}`;
                        if (!civicVariantSearchTerm) variantLinkEl.style.display = 'none';
                    }
                    variantLinkEl.target = '_blank';
                    variantLinkEl.rel = 'noopener noreferrer';
                    variantLinkEl.textContent = 'View variant on CIViC';
                    linksDiv.appendChild(variantLinkEl);
                    variantSepNode = document.createTextNode((legacyVariantId || civicVariantSearchTerm) ? ' | ' : '');
                    linksDiv.appendChild(variantSepNode);
                    geneLinkEl = document.createElement('a');
                    geneLinkEl.href = `https://civicdb.org/search?query=${encodedCivicGene}`;
                    geneLinkEl.target = '_blank';
                    geneLinkEl.rel = 'noopener noreferrer';
                    geneLinkEl.textContent = 'View gene on CIViC';
                    linksDiv.appendChild(geneLinkEl);
                }
                content.appendChild(linksDiv);

                // Display MyVariant.info CIViC data (clinically prioritised)
                if (hasMyVariantData) {
                    let evidenceItemsForDetails = [];
                    if (legacy) {
                        const mp = legacy?.molecularProfiles;
                        const legacyItems = Array.isArray(mp)
                            ? mp.flatMap((p) => Array.isArray(p?.evidenceItems) ? p.evidenceItems : [])
                            : (Array.isArray(mp?.evidenceItems) ? mp.evidenceItems : []);
                        evidenceItemsForDetails = legacyItems;

                        const diseases = new Set(), drugs = new Set();
                        const byType = {};
                        legacyItems.forEach((item) => {
                            if (item?.disease?.name) diseases.add(String(item.disease.name));
                            if (Array.isArray(item?.therapies)) item.therapies.forEach((t) => { if (t?.name) drugs.add(String(t.name)); });
                            const t = item?.evidenceType || 'OTHER';
                            byType[t] = (byType[t] || 0) + 1;
                        });

                        addLine('Gene', legacy?.gene?.name);
                        addLine('Variant', legacy?.variant?.name || legacy?.name);
                        if (legacyItems.length > 0) {
                            const typeSummary = Object.entries(byType).sort(([,a],[,b]) => b - a)
                                .map(([t, n]) => `${t} (${n})`).join(', ');
                            addLine('Evidence', `${legacyItems.length} items — ${typeSummary}`);
                        }
                        if (diseases.size > 0) addLine('Diseases', Array.from(diseases).slice(0, 4).join(', '));
                        if (drugs.size > 0) addLine('Therapies', Array.from(drugs).slice(0, 4).join(', '));

                        // Highlight top predictive/prognostic evidence
                        const topItems = legacyItems
                            .filter((i) => ['PREDICTIVE', 'PROGNOSTIC', 'DIAGNOSTIC'].includes(String(i?.evidenceType || '').toUpperCase()))
                            .sort((a, b) => {
                                const order = { A: 0, B: 1, C: 2, D: 3, E: 4 };
                                return (order[String(a?.evidenceLevel || '').toUpperCase()] ?? 9) - (order[String(b?.evidenceLevel || '').toUpperCase()] ?? 9);
                            });
                        if (topItems.length > 0) {
                            const ti = topItems[0];
                            const parts = [
                                ti?.disease?.name,
                                Array.isArray(ti?.therapies) ? ti.therapies.map((x) => x.name).join(', ') : null,
                                ti?.significance ? `Sig: ${ti.significance}` : null,
                                ti?.evidenceLevel ? `Level ${ti.evidenceLevel}` : null
                            ].filter(Boolean);
                            addLine('Top clinical evidence', parts.join(' — '));
                        }
                    } else {
                        // Array-format entries (annotation.cgi)
                        evidenceItemsForDetails = entries;
                        const diseases = new Set(), drugs = new Set();
                        const byType = {};
                        entries.forEach((e) => {
                            if (e?.primary_disease) diseases.add(String(e.primary_disease));
                            if (Array.isArray(e?.drugs)) e.drugs.forEach((d) => drugs.add(String(d)));
                            const t = e?.evidence_type || 'OTHER';
                            byType[t] = (byType[t] || 0) + 1;
                        });
                        const typeSummary = Object.entries(byType).sort(([,a],[,b]) => b - a)
                            .map(([t, n]) => `${t} (${n})`).join(', ');
                        if (entries.length > 0) addLine('Evidence', `${entries.length} items — ${typeSummary}`);
                        if (diseases.size > 0) addLine('Diseases', Array.from(diseases).slice(0, 4).join(', '));
                        if (drugs.size > 0) addLine('Therapies', Array.from(drugs).slice(0, 4).join(', '));
                    }

                    // Collapsible evidence table (sorted: high-level predictive first)
                    if (evidenceItemsForDetails.length > 0) {
                        const levelOrder = { A: 0, B: 1, C: 2, D: 3, E: 4 };
                        const typeOrder = { PREDICTIVE: 0, PROGNOSTIC: 1, DIAGNOSTIC: 2, PREDISPOSING: 3 };
                        const sorted = [...evidenceItemsForDetails].sort((a, b) => {
                            const ta = String(a?.evidenceType || a?.evidence_type || '').toUpperCase();
                            const tb = String(b?.evidenceType || b?.evidence_type || '').toUpperCase();
                            const typeDiff = (typeOrder[ta] ?? 9) - (typeOrder[tb] ?? 9);
                            if (typeDiff !== 0) return typeDiff;
                            const la = String(a?.evidenceLevel || a?.evidence_level || '').toUpperCase();
                            const lb = String(b?.evidenceLevel || b?.evidence_level || '').toUpperCase();
                            return (levelOrder[la] ?? 9) - (levelOrder[lb] ?? 9);
                        });

                        const evDetails = document.createElement('details');
                        const evSum = document.createElement('summary');
                        evSum.textContent = `MyVariant evidence details (${evidenceItemsForDetails.length})`;
                        evDetails.appendChild(evSum);
                        const evTable = document.createElement('table');
                        evTable.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.8rem;margin-top:4px;';
                        const thead = evTable.createTHead();
                        const hrow = thead.insertRow();
                        ['Type', 'Lvl', 'Disease', 'Therapies', 'Significance'].forEach((h) => {
                            const th = document.createElement('th');
                            th.textContent = h;
                            th.style.cssText = 'text-align:left;padding:2px 5px;background:#f3f4f6;font-size:0.78rem;';
                            hrow.appendChild(th);
                        });
                        const tbody = evTable.createTBody();
                        sorted.slice(0, 25).forEach((item) => {
                            const tr = tbody.insertRow();
                            const type = item?.evidenceType || item?.evidence_type || 'N/A';
                            const lvl = item?.evidenceLevel || item?.evidence_level || '';
                            const disease = item?.disease?.name || item?.primary_disease || 'N/A';
                            const therapies = Array.isArray(item?.therapies) ? item.therapies.map((t) => t.name).join(', ')
                                : (Array.isArray(item?.drugs) ? item.drugs.join(', ') : 'N/A');
                            const sig = item?.significance || item?.clinicalSignificance || 'N/A';
                            [type, lvl ? `L${lvl}` : 'N/A', disease, therapies, sig].forEach((val) => {
                                const td = tr.insertCell();
                                td.textContent = val;
                                td.style.cssText = 'padding:2px 5px;border-bottom:1px solid #f0f0f0;vertical-align:top;';
                            });
                        });
                        evDetails.appendChild(evTable);
                        if (evidenceItemsForDetails.length > 25) {
                            const more = document.createElement('div');
                            more.style.cssText = 'font-size:0.78rem;color:#666;padding:3px 5px;';
                            more.textContent = `Showing 25 of ${evidenceItemsForDetails.length} items (sorted by evidence type and level).`;
                            evDetails.appendChild(more);
                        }
                        content.appendChild(evDetails);
                    }
                } else {
                    const noData = document.createElement('div');
                    noData.style.cssText = 'color:#6b7280;font-size:0.9rem;margin:0.3rem 0;';
                    noData.textContent = 'No CIViC data in MyVariant annotation for this variant.';
                    content.appendChild(noData);
                }

                // CivicDB API section (async, non-blocking)
                const civicApiDiv = document.createElement('div');
                civicApiDiv.style.marginTop = '0.5rem';
                const civicApiSpinner = document.createElement('div');
                civicApiSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                civicApiSpinner.textContent = 'Loading CIViC API data…';
                civicApiDiv.appendChild(civicApiSpinner);
                content.appendChild(civicApiDiv);
                card.appendChild(content);
                cardsContainer.appendChild(card);

                if (civicGene) {
                    fetchCivicApiData(civicGene, civicProtein).then((civicApiData) => {
                        aiReviewExtras.civic_api = civicApiData;
                        civicApiDiv.innerHTML = '';
                        if (!civicApiData) {
                            civicApiDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">CIViC API not accessible from browser — use the links above to search CIViC directly.</div>';
                            return;
                        }
                        const { gene: apiGene, matchedVariant, assertions } = civicApiData;
                        if (!apiGene) {
                            civicApiDiv.innerHTML = `<div style="font-size:0.82rem;color:#9ca3af;">Gene "${civicGene}" not found in CIViC.</div>`;
                            return;
                        }

                        // Upgrade gene link to the specific feature page now that we have the numeric ID
                        if (apiGene.id && geneLinkEl) {
                            geneLinkEl.href = `https://civicdb.org/features/${apiGene.id}/summary`;
                        }

                        // Upgrade variant link to specific variant page if we matched one, and reveal it if it was hidden
                        if (matchedVariant?.id && variantLinkEl) {
                            variantLinkEl.href = `https://civicdb.org/variants/${matchedVariant.id}/summary`;
                            variantLinkEl.style.display = '';
                            if (variantSepNode) variantSepNode.nodeValue = ' | ';
                            variantLinkEl.textContent = 'View variant on CIViC';
                        }

                        const geneLevelHeading = document.createElement('div');
                        geneLevelHeading.style.cssText = 'font-size:0.9rem;font-weight:600;margin:0.45rem 0 0.25rem;';
                        geneLevelHeading.textContent = 'Gene-level data from CIViC:';
                        civicApiDiv.appendChild(geneLevelHeading);

                        // Show top AMP assertion level prominently
                        if (assertions && assertions.length > 0) {
                            const topAssert = assertions[0];
                            const ampLevel = topAssert.ampLevel || topAssert.significance || '';
                            if (ampLevel) {
                                const ampEl = document.createElement('div');
                                ampEl.style.cssText = 'font-size:0.9rem;font-weight:600;margin-bottom:4px;';
                                ampEl.innerHTML = `<strong>AMP/ACMG tier (CIViC):</strong> ${ampLevel}`;
                                civicApiDiv.appendChild(ampEl);
                            }
                        }

                        // Matched variant info + variant types
                        if (matchedVariant) {
                            const vTypes = Array.isArray(matchedVariant.variantTypes?.nodes)
                                ? matchedVariant.variantTypes.nodes.map((vt) => vt.name).filter(Boolean).join(', ')
                                : '';
                            const varEl = document.createElement('div');
                            varEl.style.fontSize = '0.88rem';
                            varEl.innerHTML = `<strong>CIViC variant:</strong> ${matchedVariant.name}${vTypes ? ` <span style="color:#6b7280">(${vTypes})</span>` : ''} — <a href="https://civicdb.org/variants/${matchedVariant.id}/summary" target="_blank" rel="noopener noreferrer">View ↗</a>`;
                            civicApiDiv.appendChild(varEl);
                        }

                        // Gene description (collapsible)
                        if (apiGene.description) {
                            const descDet = document.createElement('details');
                            const descSum = document.createElement('summary');
                            descSum.style.fontSize = '0.85rem';
                            descSum.textContent = 'CIViC gene description';
                            descDet.appendChild(descSum);
                            const descText = document.createElement('div');
                            descText.style.cssText = 'font-size:0.82rem;padding:4px 0;line-height:1.45;color:#374151;';
                            const desc = String(apiGene.description);
                            descText.textContent = desc.length > 600 ? desc.slice(0, 600) + '…' : desc;
                            descDet.appendChild(descText);
                            civicApiDiv.appendChild(descDet);
                        }

                        // Assertions table (AMP tiering — most clinically important)
                        if (assertions && assertions.length > 0) {
                            const assertDet = document.createElement('details');
                            const assertSum = document.createElement('summary');
                            assertSum.style.fontSize = '0.85rem';
                            assertSum.textContent = `CIViC assertions — AMP/ACMG (${assertions.length})`;
                            assertDet.appendChild(assertSum);
                            const aTable = document.createElement('table');
                            aTable.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.79rem;margin-top:4px;';
                            const aThead = aTable.createTHead();
                            const aHrow = aThead.insertRow();
                            ['AMP Level', 'Significance', 'Disease', 'Therapies'].forEach((h) => {
                                const th = document.createElement('th');
                                th.textContent = h;
                                th.style.cssText = 'text-align:left;padding:2px 5px;background:#f3f4f6;font-size:0.77rem;';
                                aHrow.appendChild(th);
                            });
                            const aTbody = aTable.createTBody();
                            assertions.slice(0, 8).forEach((a) => {
                                const tr = aTbody.insertRow();
                                const amp = a.ampLevel || 'N/A';
                                const sig = a.clinicalSignificance || a.significance || 'N/A';
                                const disease = a.disease?.name || 'N/A';
                                const therapies = Array.isArray(a.therapies?.nodes)
                                    ? a.therapies.nodes.map((t) => t.name).join(', ')
                                    : (Array.isArray(a.therapies) ? a.therapies.join(', ') : 'N/A');
                                [amp, sig, disease, therapies].forEach((val) => {
                                    const td = tr.insertCell();
                                    td.textContent = val;
                                    td.style.cssText = 'padding:2px 5px;border-bottom:1px solid #f0f0f0;vertical-align:top;';
                                });
                            });
                            assertDet.appendChild(aTable);
                            civicApiDiv.appendChild(assertDet);
                        } else {
                            const noAssert = document.createElement('div');
                            noAssert.style.cssText = 'font-size:0.82rem;color:#6b7280;margin-top:2px;';
                            noAssert.textContent = 'No accepted CIViC assertions for this gene.';
                            civicApiDiv.appendChild(noAssert);
                        }

                        // Show total variant count for gene
                        const variantCount = Array.isArray(apiGene.variants?.nodes) ? apiGene.variants.nodes.length : 0;
                        if (variantCount > 0) {
                            const vcEl = document.createElement('div');
                            vcEl.style.cssText = 'font-size:0.8rem;color:#6b7280;margin-top:3px;';
                            vcEl.textContent = `${variantCount} variant${variantCount !== 1 ? 's' : ''} catalogued in CIViC for ${apiGene.name}.`;
                            civicApiDiv.appendChild(vcEl);
                        }
                    }).catch(() => {
                        civicApiDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">CIViC API unavailable.</div>';
                    });
                } else {
                    civicApiDiv.innerHTML = '';
                }
            }
            // Card: gnomAD
            {
                // Build a richer gnomAD card showing summary and optional population details.
                const getGnomadStats = (obj) => {
                    if (!obj) return null;
                    const stats = {};
                    // Overall AF, AC, AN
                    const afVal = obj.af;
                    const acVal = obj.ac;
                    const anVal = obj.an;
                    // Extract numeric values: these fields may be objects keyed by allele.
                    const extractNumeric = (v) => {
                        if (v === null || v === undefined) return null;
                        if (typeof v === 'number') return v;
                        if (typeof v === 'object') {
                            const keys = Object.keys(v);
                            if (keys.length > 0) return v[keys[0]];
                        }
                        return null;
                    };
                    stats.af = extractNumeric(afVal);
                    stats.ac = extractNumeric(acVal);
                    stats.an = extractNumeric(anVal);
                    // Population data: gather af_*, ac_*, an_* for major populations
                    const pops = ['afr','amr','eas','fin','nfe','oth','sas','asj'];
                    const popData = [];
                    pops.forEach((p) => {
                        const afKey = `af_${p}`;
                        const acKey = `ac_${p}`;
                        const anKey = `an_${p}`;
                        let af = obj.af ? obj.af[afKey] : undefined;
                        let ac = obj.ac ? obj.ac[acKey] : undefined;
                        let an = obj.an ? obj.an[anKey] : undefined;
                        // Only include populations with non-null values
                        if (af !== undefined || ac !== undefined || an !== undefined) {
                            // Convert undefined to null
                            af = (af !== undefined ? af : null);
                            ac = (ac !== undefined ? ac : null);
                            an = (an !== undefined ? an : null);
                            popData.push({ pop: p.toUpperCase(), af, ac, an });
                        }
                    });
                    stats.popData = popData;
                    return stats;
                };
                const gnomadExome = annotation.gnomad_exome || annotation.gnomad_exomes;
                const gnomadGenome = annotation.gnomad_genome || annotation.gnomad_genomes;
                const exomeStats = getGnomadStats(gnomadExome);
                const genomeStats = getGnomadStats(gnomadGenome);
                // Build gnomAD link. Prefer raw token input so indels can
                // still open a meaningful gnomAD variant/region page.
                const gnomadLink = buildGnomadVariantUrl(rawInput, gVariant, annotation && annotation._id, annotation && annotation.vcf);
                const card = document.createElement('div');
                card.className = 'card';
                const title = document.createElement('h3');
                title.textContent = 'gnomAD';
                applyCardTheme(card, 'gnomAD');
                card.appendChild(title);
                const content = document.createElement('div');
                content.className = 'card-content';

                // ── v2.1.1 section (GRCh37, via myvariant.info) ──────────────────────
                const v2Header = document.createElement('div');
                v2Header.style.cssText = 'font-size:0.78rem;font-weight:600;color:#6b7280;margin-bottom:0.2rem;';
                v2Header.textContent = 'v2.1.1 · GRCh37 (via myvariant.info)';
                content.appendChild(v2Header);

                // Add overall AF and counts lines
                if (exomeStats) {
                    const exDiv = document.createElement('div');
                    const af = (exomeStats.af != null && !isNaN(exomeStats.af)) ? `${(exomeStats.af * 100).toFixed(4)}%` : 'N/A';
                    const ac = (exomeStats.ac != null ? exomeStats.ac : '—');
                    const an = (exomeStats.an != null ? exomeStats.an : '—');
                    exDiv.innerHTML = `<strong>Exome AF:</strong> ${af} <strong>AC/AN:</strong> ${ac}/${an}`;
                    content.appendChild(exDiv);
                }
                if (genomeStats) {
                    const gnDiv = document.createElement('div');
                    const afg = (genomeStats.af != null && !isNaN(genomeStats.af)) ? `${(genomeStats.af * 100).toFixed(4)}%` : 'N/A';
                    const acg = (genomeStats.ac != null ? genomeStats.ac : '—');
                    const ang = (genomeStats.an != null ? genomeStats.an : '—');
                    gnDiv.innerHTML = `<strong>Genome AF:</strong> ${afg} <strong>AC/AN:</strong> ${acg}/${ang}`;
                    content.appendChild(gnDiv);
                }
                if (!exomeStats && !genomeStats) {
                    const noV2 = document.createElement('div');
                    noV2.style.cssText = 'font-size:0.82rem;color:#9ca3af;';
                    noV2.textContent = 'No v2.1.1 data.';
                    content.appendChild(noV2);
                }
                // Link to gnomAD v2
                if (gnomadLink) {
                    const linkEl = document.createElement('a');
                    linkEl.href = gnomadLink;
                    linkEl.target = '_blank';
                    linkEl.rel = 'noopener noreferrer';
                    linkEl.textContent = 'View on gnomAD (v2.1.1)';
                    content.appendChild(linkEl);
                }
                // Add details for v2 population data if available
                if ((exomeStats && exomeStats.popData && exomeStats.popData.length > 0) || (genomeStats && genomeStats.popData && genomeStats.popData.length > 0)) {
                    const detailsEl = document.createElement('details');
                    const summaryEl = document.createElement('summary');
                    summaryEl.textContent = 'v2.1.1 population details';
                    detailsEl.appendChild(summaryEl);
                    const table = document.createElement('table');
                    table.style.width = '100%';
                    table.style.borderCollapse = 'collapse';
                    const thead = document.createElement('thead');
                    const hdrRow = document.createElement('tr');
                    ['Dataset','Population','AF','AC','AN'].forEach((text) => {
                        const th = document.createElement('th');
                        th.textContent = text;
                        th.style.textAlign = 'left';
                        th.style.padding = '0.25rem 0.5rem';
                        th.style.borderBottom = '1px solid #e0e0e0';
                        hdrRow.appendChild(th);
                    });
                    thead.appendChild(hdrRow);
                    table.appendChild(thead);
                    const tbody = document.createElement('tbody');
                    if (exomeStats && exomeStats.popData) {
                        exomeStats.popData.forEach((pd) => {
                            const tr = document.createElement('tr');
                            ['Exome', pd.pop,
                                (pd.af != null && !isNaN(pd.af)) ? `${(pd.af * 100).toFixed(4)}%` : '—',
                                pd.ac != null ? pd.ac : '—',
                                pd.an != null ? pd.an : '—'
                            ].forEach((val) => {
                                const td = document.createElement('td');
                                td.textContent = val;
                                td.style.padding = '0.25rem 0.5rem';
                                tr.appendChild(td);
                            });
                            tbody.appendChild(tr);
                        });
                    }
                    if (genomeStats && genomeStats.popData) {
                        genomeStats.popData.forEach((pd) => {
                            const tr = document.createElement('tr');
                            ['Genome', pd.pop,
                                (pd.af != null && !isNaN(pd.af)) ? `${(pd.af * 100).toFixed(4)}%` : '—',
                                pd.ac != null ? pd.ac : '—',
                                pd.an != null ? pd.an : '—'
                            ].forEach((val) => {
                                const td = document.createElement('td');
                                td.textContent = val;
                                td.style.padding = '0.25rem 0.5rem';
                                tr.appendChild(td);
                            });
                            tbody.appendChild(tr);
                        });
                    }
                    table.appendChild(tbody);
                    detailsEl.appendChild(table);
                    content.appendChild(detailsEl);
                }

                // ── v4.1 section (GRCh38, direct gnomAD API) ─────────────────────────
                const divider = document.createElement('hr');
                divider.style.cssText = 'margin:0.5rem 0;border:none;border-top:1px solid #e5e7eb;';
                content.appendChild(divider);

                const v4Header = document.createElement('div');
                v4Header.style.cssText = 'font-size:0.78rem;font-weight:600;color:#6b7280;margin-bottom:0.2rem;';
                v4Header.textContent = 'v4.1 · GRCh38 (gnomAD API)';
                content.appendChild(v4Header);

                const v4Section = document.createElement('div');
                const v4Loading = document.createElement('div');
                v4Loading.style.cssText = 'font-size:0.82rem;color:#9ca3af;';
                v4Loading.textContent = 'Loading gnomAD v4.1…';
                v4Section.appendChild(v4Loading);
                content.appendChild(v4Section);

                card.appendChild(content);
                cardsContainer.appendChild(card);

                // Async: query gnomAD v4 API (liftover done server-side), populate v4Section
                (() => {
                    const vcfData = annotation && annotation.vcf;
                    const tuple = buildSpliceAiLookupTuple(rawInput, gVariant);
                    let chromCoord = null, pos37Coord = null, refCoord = null, altCoord = null;
                    if (tuple) {
                        chromCoord = tuple.chrom.replace(/^chr/i, '');
                        pos37Coord = tuple.pos;
                        refCoord = tuple.ref || (vcfData && String(vcfData.ref || '').toUpperCase()) || null;
                        altCoord = tuple.alt || (vcfData && String(vcfData.alt || '').toUpperCase()) || null;
                    }
                    if (!chromCoord || !pos37Coord || !refCoord || !altCoord) {
                        v4Loading.textContent = 'gnomAD v4.1: insufficient coordinate data.';
                        return;
                    }
                    const showV4Msg = (msg) => {
                        v4Section.innerHTML = `<div style="font-size:0.82rem;color:#9ca3af;">${msg}</div>`;
                    };
                    fetchGnomadV4(chromCoord, pos37Coord, refCoord, altCoord).then((result) => {
                        // Drop sex-stratified (_XX/_XY) and 1000 Genomes (1KG:*)
                        // subpopulations from the populations arrays — matches the UI
                        // filter and avoids dumping ~100 redundant rows into the AI payload.
                        const condenseV4Populations = (src) => {
                            if (!src || !Array.isArray(src.populations)) return src;
                            const populations = src.populations.filter((p) => {
                                const id = String(p?.id || '').toUpperCase();
                                if (!id) return false;
                                if (/_XX$|_XY$/.test(id)) return false;
                                if (id.startsWith('1KG:')) return false;
                                return true;
                            });
                            return { ...src, populations };
                        };
                        const condensedV4 = (result && result.data)
                            ? {
                                ...result,
                                data: {
                                    ...result.data,
                                    genome: condenseV4Populations(result.data.genome),
                                    exome: condenseV4Populations(result.data.exome)
                                }
                            }
                            : result;
                        aiReviewExtras.gnomad_v4 = condensedV4;
                        if (!result) { showV4Msg('gnomAD v4.1 unavailable.'); return; }
                        const { status, data: v4data, grch38Id, message, detail } = result;
                        if (status === 'liftover_failed') {
                            showV4Msg(`GRCh37→GRCh38 liftover failed for ${chromCoord}:${pos37Coord}.`);
                            return;
                        }
                        if (status === 'not_found') {
                            showV4Msg(`Not found in gnomAD v4.1${grch38Id ? ` (queried: ${grch38Id})` : ''}.`);
                            return;
                        }
                        if (status === 'api_error' || status === 'error') {
                            const detailSuffix = detail ? ` — ${detail}` : '';
                            showV4Msg(`gnomAD v4.1 error${grch38Id ? ` (${grch38Id})` : ''}: ${message || 'unknown'}${detailSuffix}`);
                            return;
                        }
                        // status === 'found'
                        v4Section.innerHTML = '';
                        // Summary rows for genome and exome
                        const renderV4Summary = (label, src) => {
                            if (!src || (src.af == null && src.ac == null)) return;
                            const row = document.createElement('div');
                            const af = (src.af != null && !isNaN(src.af)) ? `${(src.af * 100).toFixed(4)}%` : 'N/A';
                            const ac = src.ac != null ? src.ac : '—';
                            const an = src.an != null ? src.an : '—';
                            row.innerHTML = `<strong>${label} AF:</strong> ${af} <strong>AC/AN:</strong> ${ac}/${an}`;
                            v4Section.appendChild(row);
                        };
                        renderV4Summary('Exome', v4data.exome);
                        renderV4Summary('Genome', v4data.genome);
                        if (!v4data.exome && !v4data.genome) {
                            const noFreq = document.createElement('div');
                            noFreq.style.cssText = 'font-size:0.82rem;color:#9ca3af;';
                            noFreq.textContent = 'Variant found in gnomAD v4.1 but no frequency data.';
                            v4Section.appendChild(noFreq);
                        }
                        // v4 link
                        if (grch38Id) {
                            const v4Link = document.createElement('a');
                            v4Link.href = `https://gnomad.broadinstitute.org/variant/${grch38Id}?dataset=gnomad_r4`;
                            v4Link.target = '_blank';
                            v4Link.rel = 'noopener noreferrer';
                            v4Link.textContent = 'View on gnomAD (v4.1)';
                            v4Section.appendChild(v4Link);
                        }
                        // Population details
                        const v4Pops = [];
                        const collectV4Pops = (label, src) => {
                            if (!src || !Array.isArray(src.populations)) return;
                            src.populations.forEach((p) => {
                                const id = (p.id || '').toUpperCase();
                                // Skip sex-stratified subgroups (_XX/_XY) and 1000 Genomes subpopulations (1KG:*)
                                if (/_XX$|_XY$/.test(id) || id.startsWith('1KG:')) return;
                                if (p.af != null || p.ac != null) {
                                    // VariantPopulation has no af field — compute from ac/an
                                    const popAf = (p.ac != null && p.an > 0) ? p.ac / p.an : null;
                                    v4Pops.push({ dataset: label, id, af: popAf, ac: p.ac, an: p.an });
                                }
                            });
                        };
                        collectV4Pops('Exome', v4data.exome);
                        collectV4Pops('Genome', v4data.genome);
                        if (v4Pops.length > 0) {
                            const v4Details = document.createElement('details');
                            const v4Summary = document.createElement('summary');
                            v4Summary.textContent = 'v4.1 population details';
                            v4Details.appendChild(v4Summary);
                            const v4Table = document.createElement('table');
                            v4Table.style.width = '100%';
                            v4Table.style.borderCollapse = 'collapse';
                            const v4Thead = document.createElement('thead');
                            const v4HdrRow = document.createElement('tr');
                            ['Dataset','Population','AF','AC','AN'].forEach((text) => {
                                const th = document.createElement('th');
                                th.textContent = text;
                                th.style.textAlign = 'left';
                                th.style.padding = '0.25rem 0.5rem';
                                th.style.borderBottom = '1px solid #e0e0e0';
                                v4HdrRow.appendChild(th);
                            });
                            v4Thead.appendChild(v4HdrRow);
                            v4Table.appendChild(v4Thead);
                            const v4Tbody = document.createElement('tbody');
                            v4Pops.forEach((pd) => {
                                const tr = document.createElement('tr');
                                [pd.dataset, pd.id,
                                    (pd.af != null && !isNaN(pd.af)) ? `${(pd.af * 100).toFixed(4)}%` : '—',
                                    pd.ac != null ? pd.ac : '—',
                                    pd.an != null ? pd.an : '—'
                                ].forEach((val) => {
                                    const td = document.createElement('td');
                                    td.textContent = val;
                                    td.style.padding = '0.25rem 0.5rem';
                                    tr.appendChild(td);
                                });
                                v4Tbody.appendChild(tr);
                            });
                            v4Table.appendChild(v4Tbody);
                            v4Details.appendChild(v4Table);
                            v4Section.appendChild(v4Details);
                        }
                    }).catch(() => {
                        showV4Msg('gnomAD v4.1 unavailable.');
                    });
                })();
            }
            // Card: Predictors
            {
                // Build a comprehensive summary of all functional predictions using emojis.
                const funcCat = detailsData.find(cat => cat.title === 'Functional Predictions');
                const items = funcCat ? funcCat.items : {};
                // Helper to classify a prediction string into an emoji
                const classify = (val, name) => {
                    if (!val) return '❔';
                    const s = String(val).toLowerCase();
                    // If value is numeric only, treat as unknown
                    if (/^\d+(\.\d+)?$/.test(s.trim())) return '❔';
                    // Special handling for MutationAssessor categories: H(high), M(medium), L(low), N(neutral)
                    if (name && name.toLowerCase() === 'mutationassessor') {
                        const m = s.match(/[a-z]/);
                        if (m) {
                            const c = m[0];
                            // High and medium scores are considered damaging/pathogenic
                            if (c === 'h' || c === 'm') return '💥';
                            // Low, benign and neutral are treated as benign
                            if (c === 'l' || c === 'b' || c === 'n') return '😇';
                            // Any other code is unknown
                            return '❔';
                        }
                    }
                    // Special handling for MutationTaster: A/D = disease causing (pathogenic), P/N = polymorphism/benign
                    if (name && name.toLowerCase().includes('mutationtaster')) {
                        // Extract the first letter code (a, d, p, n)
                        const m = s.match(/[adpn]/);
                        if (m) {
                            const c = m[0];
                            if (c === 'a' || c === 'd') return '💥';
                            if (c === 'p' || c === 'n') return '😇';
                        }
                        // If no letter code found, fall back to generic heuristics
                    }
                    // Special handling for AlphaMissense predictions.  The categories can be letters such as
                    // A (pathogenic), B (uncertain), C (benign) or the older P/B/U codes.  We interpret them as:
                    // A/P = pathogenic (💥), B = uncertain (⚠️), C = benign (😇), U or other codes = unknown (❔).
                    if (name && name.toLowerCase().includes('alphamissense')) {
                        const vals = s.split(/[,;\s]+/).filter(Boolean).map(v => v.toLowerCase());
                        // Pathogenic if any token begins with 'a' or 'p'
                        if (vals.some(v => v.startsWith('a') || v.startsWith('p'))) return '💥';
                        // Uncertain if any token begins with 'b'
                        if (vals.some(v => v.startsWith('b'))) return '⚠️';
                        // Benign if any token begins with 'c'
                        if (vals.some(v => v.startsWith('c'))) return '😇';
                        // Unknown or other codes
                        if (vals.some(v => v.startsWith('u'))) return '❔';
                        return '❔';
                    }
                    // Generic heuristics for other prediction descriptions
                    // Pathogenic: deleterious, damaging, harmful, or letter D (whole word)
                    if (/(deleterious|damaging|harmful|\bd\b)/.test(s)) return '💥';
                    // Intermediate: possibly/probably damaging or includes letter P
                    if (/(possibly|probably|intermediate|\bp\b)/.test(s)) return '⚠️';
                    // Benign or tolerated: benign, tolerated, low, neutral, or letters B/T/N
                    if (/(benign|tolerated|low|neutral|\bb\b|\bt\b|\bn\b)/.test(s)) return '😇';
                    return '❔';
                };
                const summaryParts = [];
                Object.entries(items).forEach(([name, val]) => {
                    // If the predictor value contains multiple comma-separated predictions (e.g. multiple transcripts),
                    // classify each component individually and join the resulting emojis.  Otherwise, classify the
                    // single value directly.  Deduplicate repeated emojis to avoid redundant icons.
                    let emojis = [];
                    const valStr = String(val);
                    // Split on commas or semicolons to detect multiple predictions. Trim whitespace around tokens.
                    const tokens = valStr.split(/[,;]+/).map(t => t.trim()).filter(Boolean);
                    if (tokens.length > 1) {
                        tokens.forEach(tok => {
                            // Ignore tokens that start with a non-letter (e.g. numeric scores)
                            if (!/^[A-Za-z]/.test(tok)) return;
                            // Extract the prediction code (first word) in case of extra data like scores in parentheses.
                            const code = tok.split(/\s+/)[0];
                            const emoji = classify(code, name);
                            emojis.push(emoji);
                        });
                    } else {
                        emojis.push(classify(val, name));
                    }
                    summaryParts.push(`<strong>${name}</strong>: ${emojis.join('/')}`);
                });
                const card = document.createElement('div');
                card.className = 'card';
                const title = document.createElement('h3');
                title.textContent = 'Predictors';
                applyCardTheme(card, 'Predictors');
                card.appendChild(title);
                const content = document.createElement('div');
                content.className = 'card-content';
                const summaryDiv = document.createElement('div');
                summaryDiv.className = 'predictor-summary';
                summaryDiv.innerHTML = summaryParts.join(' | ');
                content.appendChild(summaryDiv);
                // Add collapsible details if there are items
                if (Object.keys(items).length > 0) {
                    const detailsEl = document.createElement('details');
                    const summaryEl = document.createElement('summary');
                    summaryEl.textContent = 'Show details';
                    detailsEl.appendChild(summaryEl);
                    const list = document.createElement('ul');
                    Object.entries(items).forEach(([n, v]) => {
                        const li = document.createElement('li');
                        li.innerHTML = `<strong>${n}</strong>: ${v}`;
                        list.appendChild(li);
                    });
                    detailsEl.appendChild(list);
                    content.appendChild(detailsEl);
                }
                card.appendChild(content);
                cardsContainer.appendChild(card);
            }
            // Card: OncoKB
            {
                // Provide Oncogenicity classification from MyVariant if available
                const oncogenic = annotation.oncogenic || (annotation.oncokb && annotation.oncokb.oncogenic) || '';
                // Determine primary gene (first from geneNames)
                const gene = (geneNames || '').split(',')[0].trim();
                let variantLink = '';
                // Use genomic variant to construct hgvsg link if available
                if (gVariant) {
                    const m = String(gVariant).match(/^chr(\w+):g\.(.+)/);
                    if (m) {
                        const chrom = m[1];
                        const rest = m[2];
                        variantLink = `https://www.oncokb.org/hgvsg/${chrom}:g.${rest}`;
                    }
                }
                // Fallback: use gene and protein to construct gene/protein URL
                if (!variantLink && gene && protein) {
                    let prot = protein;
                    if (prot.includes(',')) prot = prot.split(',')[0];
                    prot = prot.replace(/\[|\]|'/g, '').trim();
                    prot = prot.replace(/^p\.?/i, '');
                    variantLink = `https://www.oncokb.org/gene/${gene}/${prot}`;
                }
                const geneLink = gene ? `https://www.oncokb.org/gene/${gene}` : '';
                const card = document.createElement('div');
                card.className = 'card';
                const title = document.createElement('h3');
                title.textContent = 'OncoKB';
                applyCardTheme(card, 'OncoKB');
                card.appendChild(title);
                const content = document.createElement('div');
                content.className = 'card-content';
                const sigSpan = document.createElement('span');
                sigSpan.innerHTML = `<strong>Oncogenicity:</strong> ${oncogenic || 'N/A'}`;
                content.appendChild(sigSpan);
                // Append links
                if (variantLink) {
                    const linkEl = document.createElement('a');
                    linkEl.href = variantLink;
                    linkEl.target = '_blank';
                    linkEl.rel = 'noopener noreferrer';
                    linkEl.textContent = 'Variant Page';
                    content.appendChild(linkEl);
                }
                if (geneLink) {
                    const geneEl = document.createElement('a');
                    geneEl.href = geneLink;
                    geneEl.target = '_blank';
                    geneEl.rel = 'noopener noreferrer';
                    geneEl.textContent = 'Gene Page';
                    if (variantLink) content.appendChild(document.createTextNode(' '));
                    content.appendChild(geneEl);
                }
                card.appendChild(content);
                cardsContainer.appendChild(card);
            }
            // Card: COSMIC
            {
                // Build a COSMIC card showing detailed site counts and frequencies when available.
                const cosmicExt = detailsData.find(cat => cat.title === 'COSMIC (Extended)');
                const cosmicBase = detailsData.find(cat => cat.title === 'COSMIC');
                const card = document.createElement('div');
                card.className = 'card';
                const title = document.createElement('h3');
                title.textContent = 'COSMIC';
                applyCardTheme(card, 'COSMIC');
                card.appendChild(title);
                const content = document.createElement('div');
                content.className = 'card-content';
                if (cosmicExt) {
                    // Show all available metrics from the extended COSMIC annotation.
                    const items = cosmicExt.items;
                    // Total tumors summary
                    if (items['Total Tumors'] !== undefined) {
                        const span = document.createElement('span');
                        span.innerHTML = `<strong>Found in:</strong> ${items['Total Tumors']} tumor${items['Total Tumors'] === 1 ? '' : 's'}`;
                        content.appendChild(span);
                    }
                    // Frequency overall
                    if (items['Frequency (overall)']) {
                        const p = document.createElement('p');
                        p.innerHTML = `<strong>Frequency (overall):</strong> ${items['Frequency (overall)']}`;
                        content.appendChild(p);
                    }
                    // Find the key that starts with "Frequency in" (gene-specific frequency)
                    Object.keys(items).forEach(key => {
                        if (key.startsWith('Frequency in')) {
                            const p = document.createElement('p');
                            p.innerHTML = `<strong>${key}:</strong> ${items[key]}`;
                            content.appendChild(p);
                        }
                    });
                    // Site counts: display within a collapsible details element if HTML is provided
                    if (items['Site Counts'] && typeof items['Site Counts'] === 'object' && items['Site Counts'].html) {
                        const detailsEl = document.createElement('details');
                        const summaryEl = document.createElement('summary');
                        summaryEl.textContent = 'Site counts';
                        detailsEl.appendChild(summaryEl);
                        const div = document.createElement('div');
                        div.innerHTML = items['Site Counts'].html;
                        detailsEl.appendChild(div);
                        content.appendChild(detailsEl);
                    }
                    // Omit the COSMIC gene page link, as external access now requires login.
                } else if (cosmicBase) {
                    // Fallback: use base COSMIC info (mutation frequency or count)
                    const items = cosmicBase.items;
                    if (items['Mutation Frequency'] !== undefined) {
                        const span = document.createElement('span');
                        span.innerHTML = `<strong>Mutation Frequency:</strong> ${items['Mutation Frequency']}`;
                        content.appendChild(span);
                    }
                    // Add any other COSMIC base items except those with html
                    Object.entries(items).forEach(([k,v]) => {
                        if (k === 'Mutation Frequency') return;
                        if (v && typeof v === 'object' && v.html) return;
                        const p = document.createElement('p');
                        p.innerHTML = `<strong>${k}:</strong> ${v}`;
                        content.appendChild(p);
                    });
                } else {
                    // No COSMIC annotation
                    const span = document.createElement('span');
                    span.textContent = 'No COSMIC data available.';
                    content.appendChild(span);
                }
                card.appendChild(content);
                cardsContainer.appendChild(card);
            }
            // Card: TP53 Mutation Database (shown only for TP53 variants)
            if (isTp53Gene(geneNames)) {
                const tp53Card = document.createElement('div');
                tp53Card.className = 'card';
                const tp53Title = document.createElement('h3');
                tp53Title.textContent = 'TP53 Database';
                applyCardTheme(tp53Card, 'TP53 Database');
                tp53Card.appendChild(tp53Title);
                const tp53Content = document.createElement('div');
                tp53Content.className = 'card-content';

                const normalizePlainVariant = (value) => {
                    if (!value) return '';
                    const raw = String(value).trim();
                    const noHtml = raw.replace(/<[^>]+>/g, '');
                    const first = noHtml.split(',')[0].trim();
                    return first.includes(':') ? first.split(':').slice(1).join(':').trim() : first;
                };

                const tp53Cdna = normalizePlainVariant(cDNAHTML);
                const tp53Protein = normalizePlainVariant(protein);
                const tp53Genomic = String(gVariant || '').trim();
                const summary = document.createElement('span');
                summary.innerHTML = `<strong>Detected gene:</strong> TP53`;
                tp53Content.appendChild(summary);

                const variantSummary = document.createElement('span');
                variantSummary.innerHTML = `<strong>Variant:</strong> ${tp53Protein || tp53Cdna || tp53Genomic || 'N/A'}`;
                tp53Content.appendChild(variantSummary);

                const dbHomeUrl = 'https://tp53.cancer.gov/';
                const searchByVariantUrl = 'https://tp53.cancer.gov/search_gene_by_mut';
                const googleTp53Query = encodeURIComponent(`site:tp53.cancer.gov TP53 ${tp53Protein || tp53Cdna || tp53Genomic}`.trim());
                const googleTp53Url = `https://www.google.com/search?q=${googleTp53Query}`;
                const linksLine = document.createElement('span');
                linksLine.innerHTML = `<a href="${dbHomeUrl}" target="_blank" rel="noopener noreferrer">TP53 DB Home</a> | <a href="${searchByVariantUrl}" target="_blank" rel="noopener noreferrer">Search by Gene Variants</a> | <a href="${googleTp53Url}" target="_blank" rel="noopener noreferrer">Variant lookup 🔍</a>`;
                tp53Content.appendChild(linksLine);

                try {
                    const tp53Data = await fetchTp53MutationDatabase({
                        gene: 'TP53',
                        protein: tp53Protein || '',
                        cdna: tp53Cdna || '',
                        genomic: tp53Genomic || '',
                        debug: true
                    });
                    // Strip the `debug` field (dataset fetch attempts, column
                    // listings, etc.) from the AI payload — it's only useful for the
                    // local "Debug info" pane below.
                    if (tp53Data && typeof tp53Data === 'object') {
                        const { debug: _tp53Debug, ...tp53DataForAi } = tp53Data;
                        aiReviewExtras.tp53_mutation_database = tp53DataForAi;
                    } else {
                        aiReviewExtras.tp53_mutation_database = tp53Data;
                    }
                    if (tp53Data && typeof tp53Data === 'object') {
                        const best = tp53Data.matches && tp53Data.matches[0];
                        const path = tp53Data.pathogenicity || (best ? best : null);

                        // Match status banner
                        const statusDiv = document.createElement('div');
                        statusDiv.style.cssText = 'margin:6px 0 10px; padding:6px 10px; border-radius:4px; font-size:0.88rem;';
                        if (tp53Data.match_count > 0 && best) {
                            statusDiv.style.background = '#e8f5e9';
                            statusDiv.style.borderLeft = '3px solid #4caf50';
                            const countLabel = tp53Data.match_count === 1 ? '1 record' : `${tp53Data.match_count} records`;
                            const s = document.createElement('strong');
                            s.textContent = 'Database match: ';
                            statusDiv.appendChild(s);
                            statusDiv.appendChild(document.createTextNode(`found ${countLabel} in TP53 database`));
                        } else {
                            statusDiv.style.background = '#fff8e1';
                            statusDiv.style.borderLeft = '3px solid #ffa726';
                            const s = document.createElement('strong');
                            s.textContent = 'No match: ';
                            statusDiv.appendChild(s);
                            statusDiv.appendChild(document.createTextNode('this exact variant was not found in the TP53 MutationView dataset'));
                        }
                        tp53Content.appendChild(statusDiv);

                        if (best && path) {
                            // Shared helpers
                            const makeSectionHead = (table, label) => {
                                const tr = document.createElement('tr');
                                const td = document.createElement('td');
                                td.colSpan = 2;
                                td.style.cssText = 'padding:6px 0 2px; font-weight:700; font-size:0.82rem; color:#444; border-top:1px solid #e5e5e5;';
                                td.textContent = label;
                                tr.appendChild(td);
                                table.appendChild(tr);
                            };
                            const addRow = (table, label, value) => {
                                if (!value && value !== 0) return;
                                const tr = document.createElement('tr');
                                const td1 = document.createElement('td');
                                td1.style.cssText = 'padding:3px 12px 3px 0; font-weight:600; white-space:nowrap; vertical-align:top; color:#555; font-size:0.82rem;';
                                td1.textContent = label;
                                const td2 = document.createElement('td');
                                td2.style.cssText = 'padding:3px 0; color:#222; font-size:0.84rem;';
                                td2.textContent = String(value);
                                tr.appendChild(td1);
                                tr.appendChild(td2);
                                table.appendChild(tr);
                            };

                            const pathHeading = document.createElement('div');
                            pathHeading.style.cssText = 'font-weight:700; font-size:0.85rem; color:#333; margin-bottom:4px; border-bottom:1px solid #ddd; padding-bottom:2px;';
                            pathHeading.textContent = 'Pathogenicity evidence (best match)';
                            tp53Content.appendChild(pathHeading);

                            const table = document.createElement('table');
                            table.style.cssText = 'width:100%; border-collapse:collapse; margin-bottom:8px;';

                            // — Functional classification —
                            const effect = path.effect || best.effect;
                            const effectGroup = path.effect_group || best.effect_group;
                            const taClass = path.ta_class || best.ta_class;
                            const lofClass = path.lof_class || best.lof_class;
                            const dneClass = path.dne_class || best.dne_class;
                            const structClass = path.structural_class || best.structural_class;
                            const residueFunc = path.residue_function || best.residue_function;
                            const domainFunc = path.domain_function || best.domain_function;
                            const structMotif = path.structural_motif || best.structural_motif;
                            const hasFunctional = effect || effectGroup || taClass || lofClass || dneClass || structClass || residueFunc || domainFunc;
                            if (hasFunctional) {
                                makeSectionHead(table, 'Functional classification');
                                addRow(table, 'Effect', effect);
                                addRow(table, 'Effect group', effectGroup);
                                addRow(table, 'Transactivation class', taClass);
                                addRow(table, 'DN + LOF class', lofClass);
                                addRow(table, 'Dominant-negative class', dneClass);
                                addRow(table, 'Structure/function class', structClass);
                                addRow(table, 'Residue function', residueFunc);
                                addRow(table, 'Domain function', domainFunc);
                                addRow(table, 'Structural motif', structMotif);
                                addRow(table, 'Hot spot', path.hotspot || best.hotspot);
                            }

                            // — Computational predictors —
                            const agvgd = path.agvgd_class || best.agvgd_class;
                            const sift = path.sift_class || best.sift_class;
                            const polyphen2 = path.polyphen2 || best.polyphen2;
                            const bayesDel = path.bayes_del || best.bayes_del;
                            const revel = path.revel || best.revel;
                            let compDetails = null;
                            if (agvgd || sift || polyphen2 || bayesDel || revel) {
                                compDetails = document.createElement('details');
                                const compSummary = document.createElement('summary');
                                compSummary.style.cssText = 'cursor:pointer; font-size:0.82rem; color:#4a5f73; font-weight:600; margin:4px 0;';
                                compSummary.textContent = 'Computational predictions';
                                compDetails.appendChild(compSummary);
                                const compTable = document.createElement('table');
                                compTable.style.cssText = 'border-collapse:collapse; font-size:0.82rem;';
                                addRow(compTable, 'AGVGD class', agvgd);
                                addRow(compTable, 'SIFT class', sift);
                                addRow(compTable, 'PolyPhen-2', polyphen2);
                                addRow(compTable, 'BayesDel', bayesDel);
                                addRow(compTable, 'REVEL', revel);
                                compDetails.appendChild(compTable);
                            }

                            // — Variant characterisation —
                            const variantDetails = document.createElement('details');
                            const variantSummary = document.createElement('summary');
                            variantSummary.style.cssText = 'cursor:pointer; font-size:0.82rem; color:#4a5f73; font-weight:600; margin:4px 0;';
                            variantSummary.textContent = 'Variant';
                            variantDetails.appendChild(variantSummary);
                            const variantTable = document.createElement('table');
                            variantTable.style.cssText = 'border-collapse:collapse; font-size:0.82rem;';
                            addRow(variantTable, 'Mutation type', best.mutation_type);
                            addRow(variantTable, 'Codon / location', [best.codon_number, best.exon_intron].filter(Boolean).join(' — '));
                            addRow(variantTable, 'CpG site', path.cpg_site || best.cpg_site);
                            addRow(variantTable, 'Splice site', path.splice_site || best.splice_site);
                            addRow(variantTable, 'Protein (DB)', best.protein);
                            addRow(variantTable, 'cDNA/Genomic (DB)', best.cdna_or_genomic);
                            variantDetails.appendChild(variantTable);

                            // — Epidemiological evidence —
                            const somatic = best.somatic_count;
                            const germline = best.germline_count;
                            const tcga = best.tcga_icgc_genie_count;
                            if (somatic || germline || tcga) {
                                makeSectionHead(table, 'Epidemiological evidence');
                                addRow(table, 'Somatic occurrences', somatic);
                                addRow(table, 'Germline occurrences', germline);
                                addRow(table, 'TCGA/ICGC/GENIE count', tcga);
                            }

                            // Keep functional/epidemiology visible but collapse computational + variant rows by default.
                            const collapsibleRows = [];
                            Array.from(table.querySelectorAll('tr')).forEach((tr) => {
                                const header = tr.firstChild && tr.firstChild.textContent ? tr.firstChild.textContent : '';
                                if (/Computational predictions \\(collapsed below\\)|Variant \\(collapsed below\\)/.test(header)) {
                                    collapsibleRows.push(tr);
                                    let n = tr.nextSibling;
                                    while (n && !(n.firstChild && n.firstChild.colSpan === 2)) {
                                        collapsibleRows.push(n);
                                        n = n.nextSibling;
                                    }
                                }
                            });
                            collapsibleRows.forEach((r) => { r.style.display = 'none'; });
                            tp53Content.appendChild(table);
                            if (compDetails) tp53Content.appendChild(compDetails);
                            tp53Content.appendChild(variantDetails);

                            // — Collapsible: p53 transactivation target activity —
                            const taTargets = path.ta_targets || best.ta_targets || {};
                            const taEntries = Object.entries(taTargets).filter(([, v]) => v !== '' && v !== null && v !== undefined);
                            if (taEntries.length > 0) {
                                const taDetails = document.createElement('details');
                                taDetails.style.marginTop = '4px';
                                const taSummaryEl = document.createElement('summary');
                                taSummaryEl.style.cssText = 'cursor:pointer; font-size:0.82rem; color:#4a5f73; font-weight:600;';
                                taSummaryEl.textContent = 'p53 transactivation target activity (% of wild-type)';
                                taDetails.appendChild(taSummaryEl);
                                const taNote = document.createElement('div');
                                taNote.style.cssText = 'font-size:0.78rem; color:#666; margin:3px 0 5px;';
                                taNote.textContent = 'Values show mutant p53 activity relative to wild-type (100% = fully active). Lower values indicate loss of transcriptional function.';
                                taDetails.appendChild(taNote);
                                const taTable = document.createElement('table');
                                taTable.style.cssText = 'border-collapse:collapse; font-size:0.82rem;';
                                for (const [gene, val] of taEntries) {
                                    const tr = document.createElement('tr');
                                    const td1 = document.createElement('td');
                                    td1.style.cssText = 'padding:2px 14px 2px 0; font-weight:600; color:#555;';
                                    td1.textContent = gene;
                                    const td2 = document.createElement('td');
                                    const numVal = parseFloat(val);
                                    const pct = Number.isFinite(numVal) ? `${numVal}%` : String(val);
                                    td2.style.cssText = 'padding:2px 0; color:#222;';
                                    if (Number.isFinite(numVal)) {
                                        td2.style.color = numVal >= 50 ? '#2e7d32' : numVal >= 20 ? '#e65100' : '#b71c1c';
                                    }
                                    td2.textContent = pct;
                                    tr.appendChild(td1);
                                    tr.appendChild(td2);
                                    taTable.appendChild(tr);
                                }
                                taDetails.appendChild(taTable);
                                tp53Content.appendChild(taDetails);
                            }

                            // — Collapsible: SpliceAI scores —
                            const spliceAi = path.splice_ai || best.splice_ai || {};
                            const spliceEntries = Object.entries(spliceAi).filter(([, v]) => v !== '' && v !== null && v !== undefined && parseFloat(v) > 0);
                            if (spliceEntries.length > 0) {
                                const saDetails = document.createElement('details');
                                saDetails.style.marginTop = '4px';
                                const saSummaryEl = document.createElement('summary');
                                saSummaryEl.style.cssText = 'cursor:pointer; font-size:0.82rem; color:#4a5f73; font-weight:600;';
                                saSummaryEl.textContent = 'SpliceAI delta scores';
                                saDetails.appendChild(saSummaryEl);
                                const saNote = document.createElement('div');
                                saNote.style.cssText = 'font-size:0.78rem; color:#666; margin:3px 0 5px;';
                                saNote.textContent = 'DS_AL/AG = acceptor loss/gain; DS_DL/DG = donor loss/gain. Scores ≥0.2 suggest splicing impact; ≥0.5 is high confidence.';
                                saDetails.appendChild(saNote);
                                const saTable = document.createElement('table');
                                saTable.style.cssText = 'border-collapse:collapse; font-size:0.82rem;';
                                for (const [key, val] of spliceEntries) {
                                    const tr = document.createElement('tr');
                                    const td1 = document.createElement('td');
                                    td1.style.cssText = 'padding:2px 14px 2px 0; font-weight:600; color:#555;';
                                    td1.textContent = key;
                                    const td2 = document.createElement('td');
                                    td2.style.cssText = 'padding:2px 0; color:#222;';
                                    td2.textContent = String(val);
                                    tr.appendChild(td1);
                                    tr.appendChild(td2);
                                    saTable.appendChild(tr);
                                }
                                saDetails.appendChild(saTable);
                                tp53Content.appendChild(saDetails);
                            }

                            // — Collapsible: all matches —
                            if (tp53Data.matches.length > 1) {
                                const matchDetails = document.createElement('details');
                                matchDetails.style.marginTop = '4px';
                                const matchSummaryEl = document.createElement('summary');
                                matchSummaryEl.style.cssText = 'cursor:pointer; font-size:0.8rem; color:#4a5f73;';
                                matchSummaryEl.textContent = `All matches (${tp53Data.matches.length})`;
                                matchDetails.appendChild(matchSummaryEl);
                                for (const m of tp53Data.matches) {
                                    const mDiv = document.createElement('div');
                                    mDiv.style.cssText = 'padding:4px 0 4px 8px; border-left:2px solid #ddd; margin:4px 0; font-size:0.79rem; color:#333;';
                                    const parts = [
                                        m.protein && `Protein: ${m.protein}`,
                                        m.mutation_type && `Type: ${m.mutation_type}`,
                                        m.exon_intron && `Location: ${m.exon_intron}`,
                                        m.effect && `Effect: ${m.effect}`,
                                        m.effect_group && `Group: ${m.effect_group}`,
                                        m.ta_class && `TA class: ${m.ta_class}`,
                                        m.lof_class && `LOF class: ${m.lof_class}`,
                                        m.hotspot && `Hotspot: ${m.hotspot}`,
                                        m.somatic_count && `Somatic: ${m.somatic_count}`,
                                        m.germline_count && `Germline: ${m.germline_count}`,
                                        m.tcga_icgc_genie_count && `TCGA/ICGC/GENIE: ${m.tcga_icgc_genie_count}`,
                                        m.mut_id && `ID: ${m.mut_id}`,
                                    ].filter(Boolean);
                                    mDiv.textContent = parts.join(' | ') || '(no detail fields)';
                                    matchDetails.appendChild(mDiv);
                                }
                                tp53Content.appendChild(matchDetails);
                            }
                        }

                        // Debug (collapsed, small)
                        if (tp53Data.debug) {
                            const dbgDetails = document.createElement('details');
                            dbgDetails.style.marginTop = '6px';
                            const dbgSummary = document.createElement('summary');
                            dbgSummary.style.cssText = 'cursor:pointer; font-size:0.76rem; color:#aaa;';
                            dbgSummary.textContent = 'Debug info';
                            dbgDetails.appendChild(dbgSummary);
                            const pre = document.createElement('pre');
                            pre.style.cssText = 'white-space:pre-wrap; font-size:0.72rem; color:#666; margin-top:4px;';
                            pre.textContent = JSON.stringify(tp53Data.debug, null, 2);
                            dbgDetails.appendChild(pre);
                            tp53Content.appendChild(dbgDetails);
                        }
                    }
                } catch (tp53Err) {
                    const hint = document.createElement('span');
                    hint.style.fontSize = '0.85rem';
                    hint.style.color = '#4a5f73';
                    hint.textContent = 'Live TP53 API query unavailable (ensure /api/tp53 is deployed and TP53_MUTATION_DATASET_URL is configured if needed).';
                    tp53Content.appendChild(hint);
                    console.warn('TP53 API error', tp53Err);
                }

                tp53Card.appendChild(tp53Content);
                cardsContainer.appendChild(tp53Card);
            }
            let aiReviewGene = '';
            let aiReviewSearchVariantTerm = '';
            let aiReviewCdna = '';
            let aiReviewProtein = '';

            const getAiReviewVariantCoordinates = () => {
                const tuple = buildSpliceAiLookupTuple(rawInput, gVariant);
                const vcfData = annotation && annotation.vcf;
                if (tuple) {
                    return {
                        chrom: tuple.chrom.replace(/^chr/i, ''),
                        pos37: tuple.pos,
                        ref: tuple.ref || (vcfData && String(vcfData.ref || '').toUpperCase()) || null,
                        alt: tuple.alt || (vcfData && String(vcfData.alt || '').toUpperCase()) || null
                    };
                }
                const hg19 = annotation?.hg19 || annotation?.dbsnp?.hg19;
                const chrom = annotation?.chrom || annotation?.cadd?.chrom || annotation?.dbsnp?.chrom;
                const ref = (vcfData && String(vcfData.ref || '').toUpperCase()) || null;
                const alt = (vcfData && String(vcfData.alt || '').toUpperCase()) || null;
                if (hg19?.start !== undefined && chrom) {
                    return { chrom: String(chrom).replace(/^chr/i, ''), pos37: String(hg19.start), ref, alt };
                }
                return { chrom: null, pos37: null, ref, alt };
            };

            const normaliseAiReviewVariantText = (value) => {
                if (!value) return '';
                const raw = String(value).trim().replace(/<[^>]+>/g, '');
                const first = raw.split(',')[0].trim();
                return first.includes(':') ? first.split(':').slice(1).join(':').trim() : first;
            };

            const fetchAiReviewSupplementalContext = async () => {
                const coords = getAiReviewVariantCoordinates();
                const clinvarVariantId = annotation?.clinvar?.variant_id && /^\d+$/.test(String(annotation.clinvar.variant_id))
                    ? String(annotation.clinvar.variant_id)
                    : '';
                const pubmedTerm = [aiReviewGene, aiReviewSearchVariantTerm].filter(Boolean).join(' ');
                const pubmedTumorTerm = tumorType ? [aiReviewGene, tumorType].filter(Boolean).join(' ') : '';
                const pubmedVariantTumorTerm = (tumorType && aiReviewSearchVariantTerm)
                    ? [aiReviewGene, aiReviewSearchVariantTerm, tumorType].filter(Boolean).join(' ')
                    : '';
                // Prefer card-cached PubMed results when available. The cards fetch the same
                // queries on render; re-fetching here adds parallel load on NCBI's eutils, and
                // when efetch trips the rate limit the proxy returns articles with empty
                // abstracts — making the AI payload miss abstracts the cards already have.
                const pubmedFromCacheOrFetch = (key, term) => {
                    const cached = aiReviewExtras[key];
                    if (cached && Array.isArray(cached.articles) && cached.articles.length > 0) {
                        return Promise.resolve({ total: cached.total ?? cached.articles.length, articles: cached.articles });
                    }
                    return term ? fetchPubmedArticles(term, 5) : Promise.resolve({ total: 0, articles: [] });
                };
                const spliceApiVariant = buildSpliceAiApiVariant(rawInput, gVariant, annotation);
                const supplemental = { ...aiReviewExtras };
                const tasks = [
                    ['clinvar_variant_record', clinvarVariantId ? fetchClinvarVariant(clinvarVariantId) : Promise.resolve(null)],
                    ['nearby_clinvar_variants', coords.chrom && coords.pos37 ? fetchClinvarRegionVariants(coords.chrom, coords.pos37, 5).then(r => ({
                        note: 'ClinVar variants within ±5bp of the queried position. The queried variant itself may appear here if it has its own ClinVar entry. Do not use these neighboring variants\' classifications as the classification for the queried variant.',
                        window_bp: 5,
                        variants: r.variants
                    })) : Promise.resolve({ note: 'ClinVar variants within ±5bp of the queried position.', window_bp: 5, variants: [] })],
                    ['civic_api', aiReviewGene ? fetchCivicApiData(aiReviewGene, aiReviewProtein) : Promise.resolve(null)],
                    ['gnomad_v4', coords.chrom && coords.pos37 && coords.ref && coords.alt ? fetchGnomadV4(coords.chrom, coords.pos37, coords.ref, coords.alt) : Promise.resolve(null)],
                    ['spliceai_lookup', spliceApiVariant ? fetchSpliceAiPrediction(spliceApiVariant, { hg: '37', distance: 500, mask: 0, bc: 'basic' }) : Promise.resolve(null)],
                    ['pubmed', pubmedTerm ? pubmedFromCacheOrFetch('pubmed', pubmedTerm) : Promise.resolve({ total: 0, articles: [] })],
                    ['pubmed_tumor_type', pubmedTumorTerm && pubmedTumorTerm !== pubmedTerm ? pubmedFromCacheOrFetch('pubmed_tumor_type', pubmedTumorTerm) : Promise.resolve(null)],
                    ['pubmed_variant_tumor_type', pubmedVariantTumorTerm
                        && pubmedVariantTumorTerm !== pubmedTerm
                        && pubmedVariantTumorTerm !== pubmedTumorTerm
                        ? pubmedFromCacheOrFetch('pubmed_variant_tumor_type', pubmedVariantTumorTerm)
                        : Promise.resolve(null)],
                    ['openfda_drug_labels', aiReviewGene ? fetchOpenFdaDrugLabels(aiReviewGene).catch(() => null) : Promise.resolve(null)],
                    ['tp53_mutation_database', isTp53Gene(geneNames) ? fetchTp53MutationDatabase({
                        gene: 'TP53',
                        protein: normaliseAiReviewVariantText(aiReviewProtein),
                        cdna: normaliseAiReviewVariantText(aiReviewCdna),
                        genomic: String(gVariant || '').trim(),
                        debug: true
                    }) : Promise.resolve(null)]
                ];
                const settled = await Promise.allSettled(tasks.map(([, promise]) => promise));
                settled.forEach((result, idx) => {
                    const key = tasks[idx][0];
                    supplemental[key] = result.status === 'fulfilled'
                        ? result.value
                        : { error: result.reason?.message || String(result.reason || 'Unavailable') };
                });
                // Replace raw SpliceAI payload with compact summary (top 5 transcripts by delta score)
                if (supplemental.spliceai_lookup && !supplemental.spliceai_lookup.error) {
                    const spliceSummary = getSpliceAiScoreSummary(supplemental.spliceai_lookup);
                    const topTranscripts = (spliceSummary.transcripts || [])
                        .filter(t => t.best !== null)
                        .sort((a, b) => (b.best?.value ?? 0) - (a.best?.value ?? 0))
                        .slice(0, 5);
                    supplemental.spliceai_lookup = { best: spliceSummary.best, top_transcripts: topTranscripts };
                }
                // Drop sex-stratified (_XX/_XY) and 1000 Genomes (1KG:*) subpopulations
                // from gnomAD v4 populations — matches the UI filter and avoids dumping
                // ~100 redundant rows into the AI payload.
                if (supplemental.gnomad_v4?.data) {
                    const condenseV4Populations = (src) => {
                        if (!src || !Array.isArray(src.populations)) return src;
                        const populations = src.populations.filter((p) => {
                            const id = String(p?.id || '').toUpperCase();
                            if (!id) return false;
                            if (/_XX$|_XY$/.test(id)) return false;
                            if (id.startsWith('1KG:')) return false;
                            return true;
                        });
                        return { ...src, populations };
                    };
                    supplemental.gnomad_v4 = {
                        ...supplemental.gnomad_v4,
                        data: {
                            ...supplemental.gnomad_v4.data,
                            genome: condenseV4Populations(supplemental.gnomad_v4.data.genome),
                            exome: condenseV4Populations(supplemental.gnomad_v4.data.exome)
                        }
                    };
                }
                // Strip the `debug` field (dataset fetch attempts, sample column listings, etc.)
                // from the TP53 mutation database response — only useful for the local Debug pane.
                if (supplemental.tp53_mutation_database && typeof supplemental.tp53_mutation_database === 'object' && supplemental.tp53_mutation_database.debug) {
                    const { debug: _tp53Debug, ...tp53Rest } = supplemental.tp53_mutation_database;
                    supplemental.tp53_mutation_database = tp53Rest;
                }
                // Strip gene.variants.nodes from CiViC response — all variant names, not relevant for the matched variant
                if (supplemental.civic_api?.gene?.variants?.nodes) {
                    supplemental.civic_api = {
                        ...supplemental.civic_api,
                        gene: {
                            ...supplemental.civic_api.gene,
                            variants: { totalCount: supplemental.civic_api.gene.variants.totalCount }
                        }
                    };
                }
                supplemental.lookup_coordinates = coords;
                supplemental.spliceai_variant = spliceApiVariant;
                supplemental.pubmed_query = pubmedTerm;
                supplemental.clinvar_variant_id = clinvarVariantId;
                return supplemental;
            };

            // Card: Search
            {
                // Derive a single-letter protein code for search queries. Prefer the protein change
                // extracted from the user's query (targetProtGlobal), falling back to the canonical
                // protein string if available. targetProtGlobal contains uppercase triple-coded letters
                // plus position (e.g. VAL600GLU).
                let protSingle = '';
                if (targetProtGlobal) {
                    const tmp = tripleToSingle(targetProtGlobal);
                    if (tmp) protSingle = tmp;
                }
                if (!protSingle && protein) {
                    // Attempt to parse from the canonical protein string
                    const m = String(protein).match(/([A-Za-z]{3})(\d+)([A-Za-z]{3})/);
                    if (m) {
                        const triple = (m[1] + m[2] + m[3]).toUpperCase();
                        const tmp = tripleToSingle(triple);
                        if (tmp) protSingle = tmp;
                    } else {
                        const m2 = String(protein).match(/([A-Za-z])(\d+)([A-Za-z])/);
                        if (m2) protSingle = m2[1].toUpperCase() + m2[2] + m2[3].toUpperCase();
                    }
                }
                // If protein notation is unavailable, fall back to canonical cDNA for search queries.
                const normalizeCdnaSearchTerm = (value) => {
                    if (!value) return '';
                    const raw = String(value).trim();
                    const noHtml = raw.replace(/<[^>]+>/g, '');
                    const first = noHtml.split(',')[0].trim();
                    const fromColon = first.includes(':') ? first.split(':').slice(1).join(':').trim() : first;
                    const m = fromColon.match(/c\.[^\s,;]+/i);
                    return m ? m[0] : (fromColon.startsWith('c.') ? fromColon : '');
                };
                const normalizeProteinSearchTerm = (value) => {
                    if (!value) return '';
                    const raw = String(value).trim();
                    const noHtml = raw.replace(/<[^>]+>/g, '');
                    const first = noHtml.split(',')[0].trim();
                    const fromColon = first.includes(':') ? first.split(':').slice(1).join(':').trim() : first;
                    const m = fromColon.match(/p\.[^\s,;]+/i);
                    if (m) return m[0];
                    return /^p\./i.test(fromColon) ? fromColon : '';
                };
                let protSearch = normalizeProteinSearchTerm(protein);
                if (!protSearch && transcriptsList && transcriptsList.length > 0) {
                    const canonicalTx = transcriptsList.find(t => t.canonical) || transcriptsList[0];
                    protSearch = normalizeProteinSearchTerm(canonicalTx?.protein || '');
                    if (!protSearch) {
                        const firstTxWithProt = transcriptsList.find(t => t.protein);
                        protSearch = normalizeProteinSearchTerm(firstTxWithProt?.protein || '');
                    }
                }

                let cdnaSearch = '';
                if (transcriptsList && transcriptsList.length > 0) {
                    const canonicalTx = transcriptsList.find(t => t.canonical) || transcriptsList[0];
                    cdnaSearch = normalizeCdnaSearchTerm(canonicalTx?.cDNA || '');
                }
                if (!cdnaSearch && annotation?.hgvsc) {
                    const h = Array.isArray(annotation.hgvsc) ? annotation.hgvsc[0] : annotation.hgvsc;
                    cdnaSearch = normalizeCdnaSearchTerm(h);
                }
                if (!cdnaSearch && cDNAHTML) {
                    cdnaSearch = normalizeCdnaSearchTerm(cDNAHTML);
                }
                const searchVariantTerm = protSingle || protSearch || cdnaSearch;

                const genes = geneNames ? geneNames.split(',').map(g => g.trim()).filter(Boolean) : [];
                let firstGene = genes.find(g => !isChromosomeLikeGeneSymbol(g)) || genes[0] || '';
                if ((!firstGene || isChromosomeLikeGeneSymbol(firstGene)) && geneHintGlobal) {
                    firstGene = geneHintGlobal;
                }
                aiReviewGene = firstGene;
                aiReviewSearchVariantTerm = searchVariantTerm;
                aiReviewCdna = cdnaSearch;
                aiReviewProtein = protSearch || protSingle || protein;
                const pathQuery = encodeURIComponent(`pathogenicity of ${firstGene} ${searchVariantTerm}`.trim());
                const clinicalQuery = encodeURIComponent(`clinical significance of ${firstGene} ${searchVariantTerm}`.trim());
                const pathUrl = `https://www.google.com/search?q=${pathQuery}`;
                const clinicalUrl = `https://www.google.com/search?q=${clinicalQuery}`;
                const spliceTuple = buildSpliceAiLookupTuple(rawInput, gVariant);
                // Recover REF from annotation.vcf for MNVs where delins notation drops the reference.
                const spliceRef = (spliceTuple && spliceTuple.ref) || (annotation && annotation.vcf && String(annotation.vcf.ref || '').toUpperCase()) || null;
                const spliceAlt = (spliceTuple && spliceTuple.alt) || (annotation && annotation.vcf && String(annotation.vcf.alt || '').toUpperCase()) || null;
                const spliceVariantText = (spliceTuple && spliceRef && spliceAlt) ? `${spliceTuple.chrom} ${spliceTuple.pos} ${spliceRef} ${spliceAlt}` : '';
                // SpliceAI lookup defaults to hg38 when hg is omitted. Most MyVariant coordinates
                // we surface in this app are hg19/GRCh37, so explicitly request hg=37.
                const spliceAiUrl = spliceVariantText
                    ? `https://spliceailookup.broadinstitute.org/#variant=${encodeURIComponent(spliceVariantText)}&hg=37`
                    : 'https://spliceailookup.broadinstitute.org/#hg=37';

                const card = document.createElement('div');
                card.className = 'card';
                const titleEl = document.createElement('h3');
                titleEl.textContent = 'Search';
                applyCardTheme(card, 'Search');
                card.appendChild(titleEl);
                const content = document.createElement('div');
                content.className = 'card-content';
                const span = document.createElement('span');
                span.innerHTML = `<a href="${pathUrl}" target="_blank" rel="noopener noreferrer">Pathogenicity 🔍</a> | <a href="${clinicalUrl}" target="_blank" rel="noopener noreferrer">Clinical 🔍</a>`;
                content.appendChild(span);
                card.appendChild(content);
                cardsContainer.appendChild(card);

                const spliceCard = document.createElement('div');
                spliceCard.className = 'card';
                const spliceTitle = document.createElement('h3');
                spliceTitle.textContent = 'SpliceAI';
                applyCardTheme(spliceCard, 'SpliceAI');
                spliceCard.appendChild(spliceTitle);
                const spliceContent = document.createElement('div');
                spliceContent.className = 'card-content';
                const spliceLinkLine = document.createElement('span');
                spliceLinkLine.innerHTML = `<a href="${spliceAiUrl}" target="_blank" rel="noopener noreferrer">Open SpliceAI lookup 🔍</a>`;
                spliceContent.appendChild(spliceLinkLine);
                const spliceApiVariant = buildSpliceAiApiVariant(rawInput, gVariant, annotation);
                if (spliceVariantText) {
                    const spliceHint = document.createElement('span');
                    spliceHint.style.fontSize = '0.85rem';
                    spliceHint.style.color = '#4a5f73';
                    spliceHint.textContent = `SpliceAI query: ${spliceVariantText}`;
                    spliceContent.appendChild(spliceHint);
                } else {
                    const spliceHint = document.createElement('span');
                    spliceHint.style.fontSize = '0.85rem';
                    spliceHint.style.color = '#4a5f73';
                    spliceHint.textContent = 'No explicit chr/pos/ref/alt tuple detected; opening SpliceAI home page.';
                    spliceContent.appendChild(spliceHint);
                }
                const spliceResultsDiv = document.createElement('div');
                spliceResultsDiv.style.cssText = 'margin-top:0.45rem;';
                if (spliceApiVariant) {
                    const loading = document.createElement('div');
                    loading.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                    loading.textContent = 'Loading SpliceAI scores…';
                    spliceResultsDiv.appendChild(loading);
                }
                spliceContent.appendChild(spliceResultsDiv);
                spliceCard.appendChild(spliceContent);
                cardsContainer.appendChild(spliceCard);

                if (spliceApiVariant) {
                    fetchSpliceAiPrediction(spliceApiVariant, { hg: '37', distance: 500, mask: 0, bc: 'basic' }).then((spliceData) => {
                        aiReviewExtras.spliceai_lookup = spliceData;
                        spliceResultsDiv.innerHTML = '';
                        const data = spliceData?.data || {};
                        if (data.error) {
                            const errEl = document.createElement('div');
                            errEl.style.cssText = 'font-size:0.82rem;color:#9ca3af;';
                            errEl.textContent = `SpliceAI: ${data.error}`;
                            spliceResultsDiv.appendChild(errEl);
                            return;
                        }
                        const summary = getSpliceAiScoreSummary(spliceData);
                        if (!summary.transcripts.length) {
                            spliceResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">No SpliceAI scores returned.</div>';
                            return;
                        }
                        if (summary.best) {
                            const bestEl = document.createElement('div');
                            bestEl.style.cssText = 'font-size:0.86rem;margin-bottom:4px;';
                            const pct = (summary.best.value * 100).toFixed(1);
                            bestEl.innerHTML = `<strong>Max delta score:</strong> ${summary.best.value.toFixed(3)} (${pct}%) ${summary.best.label}${summary.best.position !== null ? ` at ${summary.best.position}` : ''}${summary.best.transcript ? ` · ${summary.best.transcript}` : ''}`;
                            spliceResultsDiv.appendChild(bestEl);
                        }
                        const tableWrapper = document.createElement('div');
                        tableWrapper.style.cssText = 'overflow-x:auto;margin-top:4px;';
                        const table = document.createElement('table');
                        table.style.cssText = 'min-width:100%;border-collapse:collapse;font-size:0.8rem;white-space:nowrap;';
                        const thead = document.createElement('thead');
                        const hrow = document.createElement('tr');
                        ['Transcript', 'Gene', 'AG', 'AL', 'DG', 'DL'].forEach((h) => {
                            const th = document.createElement('th');
                            th.textContent = h;
                            th.style.cssText = 'text-align:left;padding:2px 5px;border-bottom:1px solid #e5e7eb;';
                            hrow.appendChild(th);
                        });
                        thead.appendChild(hrow);
                        table.appendChild(thead);
                        const tbody = document.createElement('tbody');
                        summary.transcripts.slice(0, 5).forEach((row) => {
                            const tr = document.createElement('tr');
                            const values = [
                                row.transcript || '—',
                                row.gene || '—',
                                ...['AG', 'AL', 'DG', 'DL'].map((label) => {
                                    const d = row.deltas.find(item => item.label === label);
                                    return d && d.value !== null ? d.value.toFixed(3) : '—';
                                })
                            ];
                            values.forEach((val) => {
                                const td = document.createElement('td');
                                td.textContent = val;
                                td.style.cssText = 'padding:2px 5px;border-bottom:1px solid #f3f4f6;vertical-align:top;';
                                tr.appendChild(td);
                            });
                            tbody.appendChild(tr);
                        });
                        table.appendChild(tbody);
                        tableWrapper.appendChild(table);
                        spliceResultsDiv.appendChild(tableWrapper);
                        const note = document.createElement('div');
                        note.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:4px;';
                        note.textContent = 'Source: Broad SpliceAI Lookup API (hg19/GRCh37, distance 500, raw scores). Scores ≥0.2 may suggest splicing impact; verify before clinical use.';
                        spliceResultsDiv.appendChild(note);
                    }).catch((err) => {
                        aiReviewExtras.spliceai_lookup = { error: err.message || 'SpliceAI unavailable', variant: spliceApiVariant };
                        spliceResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">SpliceAI scores unavailable.</div>';
                    });
                }

                // Card: PubMed
                {
                    const PUBMED_LIMIT = 5;
                    const pmCard = document.createElement('div');
                    pmCard.className = 'card';
                    const pmTitle = document.createElement('h3');
                    pmTitle.textContent = 'PubMed';
                    applyCardTheme(pmCard, 'PubMed');
                    pmCard.appendChild(pmTitle);
                    const pmContent = document.createElement('div');
                    pmContent.className = 'card-content';

                    const pmSearchTerm = [firstGene, searchVariantTerm].filter(Boolean).join(' ');
                    const pmTumorSearchTerm = tumorType ? [firstGene, tumorType].filter(Boolean).join(' ') : '';
                    const pmVariantTumorSearchTerm = (tumorType && searchVariantTerm)
                        ? [firstGene, searchVariantTerm, tumorType].filter(Boolean).join(' ')
                        : '';
                    const hasPmTumorTab = !!(pmTumorSearchTerm && pmSearchTerm !== pmTumorSearchTerm);
                    const hasPmVariantTumorTab = !!(pmVariantTumorSearchTerm
                        && pmVariantTumorSearchTerm !== pmSearchTerm
                        && pmVariantTumorSearchTerm !== pmTumorSearchTerm);

                    const buildPmResultsPanel = (container, searchTerm, extraKey) => {
                        const queryUrl = searchTerm
                            ? `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(searchTerm)}&sort=relevance`
                            : 'https://pubmed.ncbi.nlm.nih.gov/';
                        const linkEl = document.createElement('a');
                        linkEl.href = queryUrl;
                        linkEl.target = '_blank';
                        linkEl.rel = 'noopener noreferrer';
                        linkEl.textContent = 'Search PubMed ↗';
                        container.appendChild(linkEl);
                        if (searchTerm) {
                            const queryLabel = document.createElement('div');
                            queryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
                            queryLabel.textContent = `Query: "${searchTerm}"`;
                            container.appendChild(queryLabel);
                        }
                        const resultsDiv = document.createElement('div');
                        if (searchTerm) {
                            const spinner = document.createElement('div');
                            spinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                            spinner.textContent = 'Loading PubMed results…';
                            resultsDiv.appendChild(spinner);
                        }
                        container.appendChild(resultsDiv);
                        if (searchTerm) {
                            fetchPubmedArticles(searchTerm, PUBMED_LIMIT).then(({ total, articles }) => {
                                aiReviewExtras[extraKey] = { query: searchTerm, total, articles };
                                resultsDiv.innerHTML = '';
                                if (total === 0 || articles.length === 0) {
                                    resultsDiv.innerHTML = '<div style="font-size:0.85rem;color:#6b7280;">No PubMed results found.</div>';
                                    return;
                                }
                                const countEl = document.createElement('div');
                                countEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:6px;';
                                countEl.textContent = total > PUBMED_LIMIT
                                    ? `Showing ${articles.length} of ${total.toLocaleString()} results (most relevant)`
                                    : `${articles.length} result${articles.length !== 1 ? 's' : ''}`;
                                resultsDiv.appendChild(countEl);
                                articles.forEach((art) => {
                                    const artEl = document.createElement('div');
                                    artEl.style.cssText = 'margin-bottom:8px;padding-bottom:8px;border-bottom:1px solid #f0f0f0;font-size:0.82rem;';
                                    const titleLink = document.createElement('a');
                                    titleLink.href = `https://pubmed.ncbi.nlm.nih.gov/${art.pmid}/`;
                                    titleLink.target = '_blank';
                                    titleLink.rel = 'noopener noreferrer';
                                    titleLink.style.fontWeight = '600';
                                    titleLink.textContent = art.title;
                                    artEl.appendChild(titleLink);
                                    const meta = document.createElement('div');
                                    meta.style.cssText = 'color:#6b7280;margin-top:2px;';
                                    const parts = [art.authors, art.journal, art.year].filter(Boolean);
                                    meta.textContent = parts.join(' · ') + (art.pmid ? ` · PMID ${art.pmid}` : '');
                                    artEl.appendChild(meta);
                                    if (art.abstract) {
                                        const abstractDetails = document.createElement('details');
                                        abstractDetails.style.cssText = 'margin-top:4px;';
                                        const abstractSummaryEl = document.createElement('summary');
                                        abstractSummaryEl.style.cssText = 'font-size:0.80rem;color:#4b5563;cursor:pointer;padding:2px 0;list-style:revert;';
                                        abstractSummaryEl.textContent = 'Abstract';
                                        abstractDetails.appendChild(abstractSummaryEl);
                                        const abstractText = document.createElement('div');
                                        abstractText.style.cssText = 'font-size:0.80rem;color:#374151;margin-top:4px;line-height:1.5;padding:6px;background:#f9fafb;border-radius:4px;';
                                        abstractText.textContent = art.abstract;
                                        abstractDetails.appendChild(abstractText);
                                        artEl.appendChild(abstractDetails);
                                    }
                                    resultsDiv.appendChild(artEl);
                                });
                                if (total > PUBMED_LIMIT) {
                                    const seeAll = document.createElement('a');
                                    seeAll.href = queryUrl;
                                    seeAll.target = '_blank';
                                    seeAll.rel = 'noopener noreferrer';
                                    seeAll.style.fontSize = '0.82rem';
                                    seeAll.textContent = `See all ${total.toLocaleString()} results on PubMed ↗`;
                                    resultsDiv.appendChild(seeAll);
                                }
                            }).catch(() => {
                                resultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">PubMed unavailable.</div>';
                            });
                        }
                    };

                    if (hasPmTumorTab) {
                        const tabBar = document.createElement('div');
                        tabBar.className = 'card-tabs';
                        const variantBtn = document.createElement('button');
                        variantBtn.type = 'button';
                        variantBtn.className = 'card-tab-btn active';
                        variantBtn.textContent = 'Gene + Variant';
                        const tumorBtn = document.createElement('button');
                        tumorBtn.type = 'button';
                        tumorBtn.className = 'card-tab-btn';
                        tumorBtn.textContent = 'Gene + Tumor Type';
                        tabBar.appendChild(variantBtn);
                        tabBar.appendChild(tumorBtn);

                        let variantTumorBtn = null;
                        let variantTumorPanel = null;
                        if (hasPmVariantTumorTab) {
                            variantTumorBtn = document.createElement('button');
                            variantTumorBtn.type = 'button';
                            variantTumorBtn.className = 'card-tab-btn';
                            variantTumorBtn.textContent = 'Gene + Variant + Tumor Type';
                            tabBar.appendChild(variantTumorBtn);
                        }
                        pmContent.appendChild(tabBar);

                        const variantPanel = document.createElement('div');
                        variantPanel.className = 'card-tab-panel active';
                        const tumorPanel = document.createElement('div');
                        tumorPanel.className = 'card-tab-panel';
                        if (hasPmVariantTumorTab) {
                            variantTumorPanel = document.createElement('div');
                            variantTumorPanel.className = 'card-tab-panel';
                        }

                        const tabBtns = [variantBtn, tumorBtn, variantTumorBtn].filter(Boolean);
                        const tabPanels = [variantPanel, tumorPanel, variantTumorPanel].filter(Boolean);
                        const activateTab = (idx) => {
                            tabBtns.forEach((b, i) => b.classList.toggle('active', i === idx));
                            tabPanels.forEach((p, i) => p.classList.toggle('active', i === idx));
                        };
                        variantBtn.addEventListener('click', () => activateTab(0));
                        tumorBtn.addEventListener('click', () => activateTab(1));
                        if (variantTumorBtn) variantTumorBtn.addEventListener('click', () => activateTab(2));

                        buildPmResultsPanel(variantPanel, pmSearchTerm, 'pubmed');
                        buildPmResultsPanel(tumorPanel, pmTumorSearchTerm, 'pubmed_tumor_type');
                        pmContent.appendChild(variantPanel);
                        pmContent.appendChild(tumorPanel);
                        if (variantTumorPanel) {
                            buildPmResultsPanel(variantTumorPanel, pmVariantTumorSearchTerm, 'pubmed_variant_tumor_type');
                            pmContent.appendChild(variantTumorPanel);
                        }
                    } else {
                        buildPmResultsPanel(pmContent, pmSearchTerm, 'pubmed');
                    }

                    pmCard.appendChild(pmContent);
                    cardsContainer.appendChild(pmCard);
                }

                // Card: FDA-Approved Drugs (by gene)
                {
                    const fdaCard = document.createElement('div');
                    fdaCard.className = 'card';
                    const fdaTitle = document.createElement('h3');
                    fdaTitle.textContent = 'FDA-Approved Drugs (by gene)';
                    applyCardTheme(fdaCard, 'FDA-Approved Drugs (by gene)');
                    fdaCard.appendChild(fdaTitle);
                    const fdaContent = document.createElement('div');
                    fdaContent.className = 'card-content';

                    if (!firstGene) {
                        const noGeneEl = document.createElement('div');
                        noGeneEl.style.cssText = 'font-size:0.85rem;color:#6b7280;';
                        noGeneEl.textContent = 'No gene identified for FDA drug lookup.';
                        fdaContent.appendChild(noGeneEl);
                        fdaCard.appendChild(fdaContent);
                        cardsContainer.appendChild(fdaCard);
                    } else {
                        // Tab bar
                        const fdaTabBar = document.createElement('div');
                        fdaTabBar.className = 'card-tabs';
                        const compDxBtn = document.createElement('button');
                        compDxBtn.type = 'button';
                        compDxBtn.className = 'card-tab-btn active';
                        compDxBtn.textContent = 'Companion Dx';
                        const openFdaBtn = document.createElement('button');
                        openFdaBtn.type = 'button';
                        openFdaBtn.className = 'card-tab-btn';
                        openFdaBtn.textContent = 'openFDA Labels';
                        fdaTabBar.appendChild(compDxBtn);
                        fdaTabBar.appendChild(openFdaBtn);
                        fdaContent.appendChild(fdaTabBar);

                        const compDxPanel = document.createElement('div');
                        compDxPanel.className = 'card-tab-panel active';
                        const openFdaPanel = document.createElement('div');
                        openFdaPanel.className = 'card-tab-panel';

                        compDxBtn.addEventListener('click', () => {
                            compDxBtn.classList.add('active'); openFdaBtn.classList.remove('active');
                            compDxPanel.classList.add('active'); openFdaPanel.classList.remove('active');
                        });
                        openFdaBtn.addEventListener('click', () => {
                            openFdaBtn.classList.add('active'); compDxBtn.classList.remove('active');
                            openFdaPanel.classList.add('active'); compDxPanel.classList.remove('active');
                        });

                        fdaContent.appendChild(compDxPanel);
                        fdaContent.appendChild(openFdaPanel);
                        fdaCard.appendChild(fdaContent);
                        cardsContainer.appendChild(fdaCard);

                        // --- Companion Dx panel (existing logic) ---
                        const fdaCompDxUrl = 'https://www.fda.gov/medical-devices/in-vitro-diagnostics/list-cleared-or-approved-companion-diagnostic-devices-in-vitro-and-imaging-tools';
                        const fdaLinkEl = document.createElement('a');
                        fdaLinkEl.href = fdaCompDxUrl;
                        fdaLinkEl.target = '_blank';
                        fdaLinkEl.rel = 'noopener noreferrer';
                        fdaLinkEl.textContent = 'FDA companion diagnostics list ↗';
                        compDxPanel.appendChild(fdaLinkEl);

                        const fdaQueryLabel = document.createElement('div');
                        fdaQueryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 2px;';
                        fdaQueryLabel.textContent = `Gene: ${firstGene}`;
                        compDxPanel.appendChild(fdaQueryLabel);

                        const fdaDisclaimer = document.createElement('div');
                        fdaDisclaimer.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin:0 0 6px;font-style:italic;';
                        fdaDisclaimer.textContent = 'Note: FDA list does not always note resistance mutations.';
                        compDxPanel.appendChild(fdaDisclaimer);

                        const fdaResultsDiv = document.createElement('div');
                        const fdaSpinner = document.createElement('div');
                        fdaSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                        fdaSpinner.textContent = 'Loading FDA drug data…';
                        fdaResultsDiv.appendChild(fdaSpinner);
                        compDxPanel.appendChild(fdaResultsDiv);

                        fetchFdaCompanionDiagnostics(firstGene).then((records) => {
                            fdaResultsDiv.innerHTML = '';
                            if (!records || records.length === 0) {
                                const noResults = document.createElement('div');
                                noResults.style.cssText = 'font-size:0.85rem;color:#6b7280;';
                                noResults.textContent = `No FDA companion diagnostic records found for ${firstGene}.`;
                                fdaResultsDiv.appendChild(noResults);
                                return;
                            }
                            const countEl = document.createElement('div');
                            countEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:6px;';
                            countEl.textContent = `${records.length} record${records.length !== 1 ? 's' : ''}`;
                            fdaResultsDiv.appendChild(countEl);

                            const table = document.createElement('table');
                            table.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.82rem;';
                            const thead = document.createElement('thead');
                            const headerRow = document.createElement('tr');
                            ['Drug', 'Disease', 'Biomarker detail'].forEach((col) => {
                                const th = document.createElement('th');
                                th.style.cssText = 'text-align:left;padding:4px 6px;border-bottom:2px solid #fca5a5;background:#fff7f7;color:#7f1d1d;font-weight:600;';
                                th.textContent = col;
                                headerRow.appendChild(th);
                            });
                            thead.appendChild(headerRow);
                            table.appendChild(thead);

                            const FDA_PREVIEW_ROWS = 5;
                            const buildRow = (rec, i) => {
                                const tr = document.createElement('tr');
                                tr.style.background = i % 2 === 0 ? '#fff' : '#fff7f7';
                                const cellStyle = 'padding:4px 6px;vertical-align:top;border-bottom:1px solid #fee2e2;';

                                const tdDrug = document.createElement('td');
                                tdDrug.style.cssText = cellStyle + 'font-weight:600;';
                                const drugsArr = Array.isArray(rec.therapy?.drugs) ? rec.therapy.drugs : [];
                                tdDrug.textContent = drugsArr.length
                                    ? drugsArr.map(d => d.trade_name ? `${d.trade_name} (${d.generic_name})` : d.generic_name || d.raw || '').filter(Boolean).join(', ')
                                    : (rec.therapy?.raw || '—');

                                const tdDisease = document.createElement('td');
                                tdDisease.style.cssText = cellStyle + 'color:#374151;';
                                const rawDisease = rec.indication?.raw || rec.indication?.disease || '—';
                                tdDisease.textContent = rawDisease.replace(/\s*[-–]\s*(Tissue|Plasma|Blood|Serum|Urine|FFPE|Fresh Frozen|Whole Blood|ctDNA)\s*$/i, '').trim() || rawDisease;

                                const tdDetail = document.createElement('td');
                                tdDetail.style.cssText = cellStyle + 'color:#374151;';
                                tdDetail.textContent = rec.biomarker?.details || rec.biomarker?.name || '—';

                                tr.appendChild(tdDrug);
                                tr.appendChild(tdDisease);
                                tr.appendChild(tdDetail);
                                return tr;
                            };

                            const tbody = document.createElement('tbody');
                            records.slice(0, FDA_PREVIEW_ROWS).forEach((rec, i) => tbody.appendChild(buildRow(rec, i)));
                            table.appendChild(tbody);
                            fdaResultsDiv.appendChild(table);

                            if (records.length > FDA_PREVIEW_ROWS) {
                                const details = document.createElement('details');
                                details.style.cssText = 'margin-top:2px;';
                                const summary = document.createElement('summary');
                                summary.style.cssText = 'font-size:0.82rem;color:#7f1d1d;cursor:pointer;padding:4px 2px;list-style:revert;';
                                summary.textContent = `Show ${records.length - FDA_PREVIEW_ROWS} more…`;
                                details.appendChild(summary);

                                const extraTable = document.createElement('table');
                                extraTable.style.cssText = 'width:100%;border-collapse:collapse;font-size:0.82rem;';
                                const extraTbody = document.createElement('tbody');
                                records.slice(FDA_PREVIEW_ROWS).forEach((rec, i) => extraTbody.appendChild(buildRow(rec, FDA_PREVIEW_ROWS + i)));
                                extraTable.appendChild(extraTbody);
                                details.appendChild(extraTable);
                                fdaResultsDiv.appendChild(details);
                            }

                            const note = document.createElement('div');
                            note.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                            note.textContent = 'Source: FDA companion diagnostics list. Verify against current FDA labeling before clinical use.';
                            fdaResultsDiv.appendChild(note);
                        }).catch(() => {
                            fdaResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">FDA drug data unavailable.</div>';
                        });

                        // --- openFDA Labels panel ---
                        const ofLinkEl = document.createElement('a');
                        ofLinkEl.href = `https://open.fda.gov/apis/drug/label/`;
                        ofLinkEl.target = '_blank';
                        ofLinkEl.rel = 'noopener noreferrer';
                        ofLinkEl.textContent = 'openFDA drug label database ↗';
                        openFdaPanel.appendChild(ofLinkEl);

                        const ofQueryLabel = document.createElement('div');
                        ofQueryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
                        ofQueryLabel.textContent = `Searching indications_and_usage for: "${firstGene}"`;
                        openFdaPanel.appendChild(ofQueryLabel);

                        const ofResultsDiv = document.createElement('div');
                        const ofSpinner = document.createElement('div');
                        ofSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                        ofSpinner.textContent = 'Loading openFDA results…';
                        ofResultsDiv.appendChild(ofSpinner);
                        openFdaPanel.appendChild(ofResultsDiv);

                        fetchOpenFdaDrugLabels(firstGene).then(({ total, fetched, excluded, excludedCase, excludedBoundary, excludedNegation, results: ofResults }) => {
                            aiReviewExtras.openfda = { gene: firstGene, total, fetched, excluded, excludedCase, excludedBoundary, excludedNegation, results: ofResults };
                            ofResultsDiv.innerHTML = '';
                            if (!ofResults || ofResults.length === 0) {
                                ofResultsDiv.innerHTML = `<div style="font-size:0.85rem;color:#6b7280;">No openFDA drug label results found for ${firstGene}.</div>`;
                                return;
                            }
                            const OF_PREVIEW = 7;
                            const ofCountEl = document.createElement('div');
                            ofCountEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:8px;';
                            const shownOf = `${Math.min(ofResults.length, OF_PREVIEW)} of ${ofResults.length} result${ofResults.length !== 1 ? 's' : ''}`;
                            const excludedReasons = [];
                            if (excludedCase) excludedReasons.push(`${excludedCase} case mismatch`);
                            if (excludedBoundary) excludedReasons.push(`${excludedBoundary} substring of larger word`);
                            if (excludedNegation) excludedReasons.push(`${excludedNegation} negated mention`);
                            const excludedNote = excludedReasons.length ? ` (excluded: ${excludedReasons.join(', ')})` : '';
                            ofCountEl.textContent = shownOf + excludedNote;
                            ofResultsDiv.appendChild(ofCountEl);

                            const buildDrugEl = (item) => {
                                const drugEl = document.createElement('div');
                                drugEl.style.cssText = 'margin-bottom:8px;padding:8px;background:#fff7f7;border:1px solid #fee2e2;border-radius:4px;font-size:0.82rem;';

                                const nameEl = document.createElement('div');
                                nameEl.style.cssText = 'font-weight:700;color:#7f1d1d;margin-bottom:2px;';
                                const brandPart = item.brand_name ? item.brand_name : '';
                                const genericPart = item.generic_name ? (brandPart ? `(${item.generic_name})` : item.generic_name) : '';
                                nameEl.textContent = [brandPart, genericPart].filter(Boolean).join(' ') || 'Unknown drug';
                                drugEl.appendChild(nameEl);

                                const metaParts = [item.manufacturer, item.route].filter(Boolean);
                                if (metaParts.length) {
                                    const metaEl = document.createElement('div');
                                    metaEl.style.cssText = 'color:#6b7280;margin-bottom:4px;';
                                    metaEl.textContent = metaParts.join(' · ');
                                    drugEl.appendChild(metaEl);
                                }

                                if (item.indications_and_usage) {
                                    const indDetails = document.createElement('details');
                                    indDetails.style.cssText = 'margin-top:2px;';
                                    const indSummary = document.createElement('summary');
                                    indSummary.style.cssText = 'cursor:pointer;color:#991b1b;font-weight:600;font-size:0.8rem;list-style:revert;padding:2px 0;';
                                    indSummary.textContent = 'Indications & Usage';
                                    indDetails.appendChild(indSummary);
                                    const indText = document.createElement('div');
                                    indText.style.cssText = 'margin-top:4px;color:#374151;line-height:1.5;white-space:pre-wrap;font-size:0.79rem;max-height:300px;overflow-y:auto;';
                                    indText.textContent = item.indications_and_usage;
                                    indDetails.appendChild(indText);
                                    drugEl.appendChild(indDetails);
                                }

                                return drugEl;
                            };

                            ofResults.slice(0, OF_PREVIEW).forEach((item) => ofResultsDiv.appendChild(buildDrugEl(item)));

                            if (ofResults.length > OF_PREVIEW) {
                                const ofMoreDetails = document.createElement('details');
                                ofMoreDetails.style.cssText = 'margin-top:2px;';
                                const ofMoreSummary = document.createElement('summary');
                                ofMoreSummary.style.cssText = 'font-size:0.82rem;color:#7f1d1d;cursor:pointer;padding:4px 2px;list-style:revert;';
                                ofMoreSummary.textContent = `Show ${ofResults.length - OF_PREVIEW} more…`;
                                ofMoreDetails.appendChild(ofMoreSummary);
                                ofResults.slice(OF_PREVIEW).forEach((item) => ofMoreDetails.appendChild(buildDrugEl(item)));
                                ofResultsDiv.appendChild(ofMoreDetails);
                            }

                            const ofNote = document.createElement('div');
                            ofNote.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                            ofNote.textContent = 'Source: openFDA drug label API. Results show labels mentioning this gene in indications. Verify against current FDA labeling.';
                            ofResultsDiv.appendChild(ofNote);
                        }).catch(() => {
                            ofResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">openFDA data unavailable.</div>';
                        });
                    }
                }

                // Card: Clinical Trials
                {
                    const CT_PREVIEW = 5;
                    const ctCard = document.createElement('div');
                    ctCard.className = 'card';
                    const ctTitle = document.createElement('h3');
                    ctTitle.textContent = 'Clinical Trials';
                    applyCardTheme(ctCard, 'Clinical Trials');
                    ctCard.appendChild(ctTitle);
                    const ctContent = document.createElement('div');
                    ctContent.className = 'card-content';

                    const ctSearchGene = firstGene;
                    const ctBaseUrl = ctSearchGene
                        ? `https://clinicaltrials.gov/search?query=${encodeURIComponent([ctSearchGene, tumorType].filter(Boolean).join(' '))}&recrs=b&type=Intr`
                        : 'https://clinicaltrials.gov/';

                    const ctLinkEl = document.createElement('a');
                    ctLinkEl.href = ctBaseUrl;
                    ctLinkEl.target = '_blank';
                    ctLinkEl.rel = 'noopener noreferrer';
                    ctLinkEl.textContent = 'Search ClinicalTrials.gov ↗';
                    ctContent.appendChild(ctLinkEl);

                    if (ctSearchGene) {
                        const ctQueryLabel = document.createElement('div');
                        ctQueryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
                        const queryParts = [ctSearchGene, tumorType].filter(Boolean);
                        ctQueryLabel.textContent = `Query: "${queryParts.join(' + ')}" · Interventional · Recruiting · Phase 2+ · US`;
                        ctContent.appendChild(ctQueryLabel);
                    }

                    const ctResultsDiv = document.createElement('div');
                    if (ctSearchGene) {
                        const ctSpinner = document.createElement('div');
                        ctSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                        ctSpinner.textContent = 'Loading clinical trials…';
                        ctResultsDiv.appendChild(ctSpinner);
                    }
                    ctContent.appendChild(ctResultsDiv);
                    ctCard.appendChild(ctContent);
                    cardsContainer.appendChild(ctCard);

                    if (ctSearchGene) {
                        fetchClinicalTrials(ctSearchGene, tumorType).then(({ total, studies }) => {
                            ctResultsDiv.innerHTML = '';
                            if (total === 0 || studies.length === 0) {
                                ctResultsDiv.innerHTML = '<div style="font-size:0.85rem;color:#6b7280;">No recruiting Phase 2+ interventional trials found in the US.</div>';
                                return;
                            }
                            const countEl = document.createElement('div');
                            countEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:8px;';
                            countEl.textContent = `${total} recruiting trial${total !== 1 ? 's' : ''} found`;
                            ctResultsDiv.appendChild(countEl);

                            const buildTrialEl = (trial) => {
                                const el = document.createElement('div');
                                el.style.cssText = 'margin-bottom:10px;padding-bottom:10px;border-bottom:1px solid #ccfbf1;font-size:0.82rem;';

                                const titleLine = document.createElement('div');
                                if (trial.url) {
                                    const link = document.createElement('a');
                                    link.href = trial.url;
                                    link.target = '_blank';
                                    link.rel = 'noopener noreferrer';
                                    link.style.fontWeight = '600';
                                    link.textContent = trial.title || trial.nctId;
                                    titleLine.appendChild(link);
                                } else {
                                    titleLine.style.fontWeight = '600';
                                    titleLine.textContent = trial.title || trial.nctId;
                                }
                                el.appendChild(titleLine);

                                const meta = document.createElement('div');
                                meta.style.cssText = 'color:#6b7280;margin-top:3px;';
                                const phaseParts = Array.isArray(trial.phases) && trial.phases.length
                                    ? trial.phases.map(p => String(p).replace('PHASE', 'Phase ')).join('/')
                                    : 'Phase N/A';
                                const drugNames = (trial.interventions || [])
                                    .filter(i => i.type === 'DRUG' || i.type === 'BIOLOGICAL' || i.type === 'COMBINATION_PRODUCT')
                                    .map(i => i.name)
                                    .filter(Boolean)
                                    .slice(0, 4);

                                const metaParts = [
                                    trial.nctId,
                                    phaseParts,
                                    drugNames.length ? drugNames.join(', ') : null,
                                    trial.usLocationCount ? `${trial.usLocationCount} US site${trial.usLocationCount !== 1 ? 's' : ''}` : null
                                ].filter(Boolean);
                                meta.textContent = metaParts.join(' · ');
                                el.appendChild(meta);

                                const hasExpandable = trial.briefSummary || (trial.conditions && trial.conditions.length > 0) || trial.inclusionCriteria;
                                if (hasExpandable) {
                                    const trialDetails = document.createElement('details');
                                    trialDetails.style.cssText = 'margin-top:5px;';
                                    const trialDetailsSummary = document.createElement('summary');
                                    trialDetailsSummary.style.cssText = 'font-size:0.80rem;color:#0f766e;cursor:pointer;padding:2px 0;list-style:revert;';
                                    trialDetailsSummary.textContent = 'Summary, conditions & eligibility';
                                    trialDetails.appendChild(trialDetailsSummary);
                                    const expandContent = document.createElement('div');
                                    expandContent.style.cssText = 'font-size:0.80rem;color:#374151;margin-top:4px;line-height:1.5;padding:6px;background:#f0fdfa;border-radius:4px;';
                                    if (trial.conditions && trial.conditions.length > 0) {
                                        const condLabel = document.createElement('div');
                                        condLabel.style.cssText = 'font-weight:600;margin-bottom:2px;';
                                        condLabel.textContent = 'Conditions:';
                                        expandContent.appendChild(condLabel);
                                        const condVal = document.createElement('div');
                                        condVal.style.cssText = 'margin-bottom:6px;';
                                        condVal.textContent = trial.conditions.join(', ');
                                        expandContent.appendChild(condVal);
                                    }
                                    if (trial.briefSummary) {
                                        const summLabel = document.createElement('div');
                                        summLabel.style.cssText = 'font-weight:600;margin-bottom:2px;';
                                        summLabel.textContent = 'Summary:';
                                        expandContent.appendChild(summLabel);
                                        const summVal = document.createElement('div');
                                        summVal.style.cssText = 'margin-bottom:6px;';
                                        summVal.textContent = trial.briefSummary;
                                        expandContent.appendChild(summVal);
                                    }
                                    if (trial.inclusionCriteria) {
                                        const inclLabel = document.createElement('div');
                                        inclLabel.style.cssText = 'font-weight:600;margin-bottom:2px;';
                                        inclLabel.textContent = 'Inclusion Criteria:';
                                        expandContent.appendChild(inclLabel);
                                        const inclVal = document.createElement('pre');
                                        inclVal.style.cssText = 'white-space:pre-wrap;font-family:inherit;margin:0;';
                                        inclVal.textContent = trial.inclusionCriteria;
                                        expandContent.appendChild(inclVal);
                                    }
                                    trialDetails.appendChild(expandContent);
                                    el.appendChild(trialDetails);
                                }
                                return el;
                            };

                            const previewStudies = studies.slice(0, CT_PREVIEW);
                            previewStudies.forEach(trial => ctResultsDiv.appendChild(buildTrialEl(trial)));

                            if (studies.length > CT_PREVIEW) {
                                const moreDetails = document.createElement('details');
                                moreDetails.style.cssText = 'margin-top:4px;';
                                const moreSummary = document.createElement('summary');
                                moreSummary.style.cssText = 'font-size:0.82rem;color:#0f766e;cursor:pointer;padding:4px 2px;list-style:revert;';
                                moreSummary.textContent = `Show ${studies.length - CT_PREVIEW} more…`;
                                moreDetails.appendChild(moreSummary);
                                studies.slice(CT_PREVIEW).forEach(trial => moreDetails.appendChild(buildTrialEl(trial)));
                                ctResultsDiv.appendChild(moreDetails);
                            }

                            const ctNote = document.createElement('div');
                            ctNote.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                            ctNote.textContent = 'Source: ClinicalTrials.gov. Eligibility criteria not evaluated — verify before clinical use.';
                            ctResultsDiv.appendChild(ctNote);
                        }).catch(() => {
                            ctResultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">Clinical trials data unavailable.</div>';
                        });
                    }
                }
            }
            // Card: Guidelines
            {
                const GUIDELINES_BASE = 'https://drdoubleb.com/guidelines';
                const guidelinesCard = document.createElement('div');
                guidelinesCard.className = 'card';
                const guidelinesTitle = document.createElement('h3');
                guidelinesTitle.textContent = 'Guidelines';
                applyCardTheme(guidelinesCard, 'Guidelines');
                guidelinesCard.appendChild(guidelinesTitle);

                const guidelinesContent = document.createElement('div');
                guidelinesContent.className = 'card-content';

                const guidelinesIntro = document.createElement('p');
                guidelinesIntro.style.cssText = 'font-size:0.85rem;color:#6b7280;margin-bottom:8px;';
                guidelinesIntro.textContent = 'Select a cancer type to retrieve guideline recommendations for this gene.';
                guidelinesContent.appendChild(guidelinesIntro);

                const dropdownRow = document.createElement('div');
                dropdownRow.style.cssText = 'display:flex;align-items:center;gap:8px;margin-bottom:10px;';

                const cancerSelect = document.createElement('select');
                cancerSelect.style.cssText = 'font-size:0.85rem;padding:4px 8px;border:1px solid #d1d5db;border-radius:4px;flex:1;max-width:320px;';

                const defaultOpt = document.createElement('option');
                defaultOpt.value = '';
                defaultOpt.textContent = 'Loading cancer types…';
                defaultOpt.disabled = true;
                defaultOpt.selected = true;
                cancerSelect.appendChild(defaultOpt);
                dropdownRow.appendChild(cancerSelect);
                guidelinesContent.appendChild(dropdownRow);

                const guidelinesResults = document.createElement('div');
                guidelinesContent.appendChild(guidelinesResults);

                guidelinesCard.appendChild(guidelinesContent);
                cardsContainer.appendChild(guidelinesCard);

                (async () => {
                    try {
                        const resp = await fetch(`${GUIDELINES_BASE}/api/cancer-types.php`);
                        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
                        const data = await resp.json();
                        cancerSelect.innerHTML = '';
                        const placeholder = document.createElement('option');
                        placeholder.value = '';
                        placeholder.textContent = `— Select cancer type (${data.count} available) —`;
                        placeholder.disabled = true;
                        placeholder.selected = true;
                        cancerSelect.appendChild(placeholder);
                        (data.cancer_types || []).forEach(ct => {
                            const opt = document.createElement('option');
                            opt.value = ct.name;
                            opt.textContent = ct.name + (ct.record_count ? ` (${ct.record_count} records)` : '');
                            cancerSelect.appendChild(opt);
                        });
                        if (tumorType) {
                            const tumorLower = tumorType.toLowerCase().trim();
                            const matched = (data.cancer_types || []).find(ct =>
                                (ct.aliases || []).some(a => a.toLowerCase() === tumorLower) ||
                                ct.name.toLowerCase() === tumorLower ||
                                ct.name.toLowerCase().includes(tumorLower) ||
                                tumorLower.includes(ct.name.toLowerCase())
                            );
                            if (matched) {
                                cancerSelect.value = matched.name;
                                cancerSelect.dispatchEvent(new Event('change'));
                            }
                        }
                    } catch (err) {
                        cancerSelect.innerHTML = '';
                        const errOpt = document.createElement('option');
                        errOpt.textContent = 'Failed to load cancer types';
                        errOpt.disabled = true;
                        errOpt.selected = true;
                        cancerSelect.appendChild(errOpt);
                    }
                })();

                cancerSelect.addEventListener('change', async () => {
                    const selectedCancer = cancerSelect.value;
                    if (!selectedCancer) return;
                    const gene = aiReviewGene || (geneNames && geneNames[0]) || '';
                    if (!gene) {
                        guidelinesResults.innerHTML = '<div style="font-size:0.85rem;color:#9ca3af;">No gene available for guideline lookup.</div>';
                        return;
                    }
                    guidelinesResults.innerHTML = '<div style="font-size:0.85rem;color:#6b7280;padding:4px 0;">Loading guidelines…</div>';
                    try {
                        const params = new URLSearchParams({ cancer: selectedCancer, gene });
                        const resp = await fetch(`${GUIDELINES_BASE}/api/search.php?${params}`);
                        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
                        const data = await resp.json();

                        aiReviewExtras.guidelines = {
                            cancer_type: selectedCancer,
                            gene,
                            query: data.query,
                            count: data.count,
                            results: data.results
                        };

                        guidelinesResults.innerHTML = '';

                        if (!data.results || data.results.length === 0) {
                            const noResults = document.createElement('div');
                            noResults.style.cssText = 'font-size:0.85rem;color:#9ca3af;';
                            noResults.textContent = `No guideline records found for ${gene} in ${selectedCancer}.`;
                            guidelinesResults.appendChild(noResults);
                            return;
                        }

                        const countEl = document.createElement('div');
                        countEl.style.cssText = 'font-size:0.85rem;font-weight:600;color:#3f6212;margin-bottom:8px;';
                        countEl.textContent = `${data.count} guideline record${data.count !== 1 ? 's' : ''} for ${gene} in ${selectedCancer}`;
                        guidelinesResults.appendChild(countEl);

                        const renderKV = (obj, container, depth) => {
                            Object.entries(obj).forEach(([key, val]) => {
                                if (val === null || val === undefined || val === '' || (Array.isArray(val) && val.length === 0) || (typeof val === 'object' && !Array.isArray(val) && Object.keys(val).length === 0)) return;
                                const fmtKey = key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
                                if (typeof val === 'object' && !Array.isArray(val)) {
                                    const subTitle = document.createElement('div');
                                    subTitle.style.cssText = `padding-left:${depth * 12}px;font-weight:600;color:#374151;margin-top:3px;`;
                                    subTitle.textContent = fmtKey + ':';
                                    container.appendChild(subTitle);
                                    renderKV(val, container, depth + 1);
                                } else if (Array.isArray(val) && val.some(v => v !== null && typeof v === 'object')) {
                                    // Array of objects — render each item as a nested block
                                    const subTitle = document.createElement('div');
                                    subTitle.style.cssText = `padding-left:${depth * 12}px;font-weight:600;color:#374151;margin-top:3px;`;
                                    subTitle.textContent = fmtKey + ':';
                                    container.appendChild(subTitle);
                                    val.forEach(item => {
                                        if (item !== null && typeof item === 'object') {
                                            renderKV(item, container, depth + 1);
                                        } else if (item !== null && item !== '') {
                                            const row = document.createElement('div');
                                            row.style.cssText = `padding-left:${(depth + 1) * 12}px;margin:2px 0;color:#1f2937;`;
                                            row.textContent = String(item);
                                            container.appendChild(row);
                                        }
                                    });
                                } else {
                                    const row = document.createElement('div');
                                    row.style.cssText = `padding-left:${depth * 12}px;margin:2px 0;`;
                                    const labelEl = document.createElement('span');
                                    labelEl.style.cssText = 'font-weight:600;color:#374151;';
                                    labelEl.textContent = fmtKey + ': ';
                                    const valEl = document.createElement('span');
                                    valEl.style.cssText = 'color:#1f2937;';
                                    valEl.textContent = Array.isArray(val) ? val.filter(v => v !== null && v !== '').join(', ') : String(val);
                                    row.appendChild(labelEl);
                                    row.appendChild(valEl);
                                    container.appendChild(row);
                                }
                            });
                        };

                        const renderSection = (title, obj, container) => {
                            if (!obj || typeof obj !== 'object' || Object.keys(obj).length === 0) return;
                            const sectionTitle = document.createElement('div');
                            sectionTitle.style.cssText = 'font-weight:700;color:#166534;margin:6px 0 2px;font-size:0.79rem;text-transform:uppercase;letter-spacing:0.04em;border-top:1px solid #bbf7d0;padding-top:5px;';
                            sectionTitle.textContent = title;
                            container.appendChild(sectionTitle);
                            renderKV(obj, container, 0);
                        };

                        data.results.forEach((record, idx) => {
                            const recordEl = document.createElement('div');
                            recordEl.style.cssText = 'margin-bottom:8px;padding:8px;background:#f0fdf4;border:1px solid #bbf7d0;border-radius:4px;font-size:0.82rem;';

                            const headerEl = document.createElement('div');
                            headerEl.style.cssText = 'font-weight:700;color:#166534;margin-bottom:4px;';
                            const markerDisplay = (record.marker && (record.marker.display_name || record.marker.symbol)) || `Record ${idx + 1}`;
                            const tumorDisplay = (record.tumor && (record.tumor.name || record.tumor.subtype)) || selectedCancer;
                            headerEl.textContent = `${markerDisplay} — ${tumorDisplay}`;
                            recordEl.appendChild(headerEl);

                            const buildKVRow = (label, value) => {
                                if (!value) return;
                                const row = document.createElement('div');
                                row.style.cssText = 'margin-bottom:3px;';
                                const lbl = document.createElement('span');
                                lbl.style.cssText = 'font-weight:600;';
                                lbl.textContent = label + ': ';
                                row.appendChild(lbl);
                                row.appendChild(document.createTextNode(value));
                                recordEl.appendChild(row);
                            };

                            if (record.therapy) {
                                const t = record.therapy;
                                if (t.recommended_therapies && t.recommended_therapies.length > 0) {
                                    buildKVRow('Recommended therapy', t.recommended_therapies.join(', '));
                                }
                                if (t.category) buildKVRow('Category', t.category);
                                if (t.context) buildKVRow('Therapy context', t.context);
                            }
                            if (record.testing) {
                                const te = record.testing;
                                if (te.recommended_test) buildKVRow('Recommended test', te.recommended_test);
                                if (te.context) buildKVRow('Testing context', te.context);
                                if (te.modality) buildKVRow('Modality', te.modality);
                            }

                            const detailsEl = document.createElement('details');
                            detailsEl.style.cssText = 'margin-top:6px;';
                            const detailsSummary = document.createElement('summary');
                            detailsSummary.style.cssText = 'font-size:0.80rem;color:#15803d;cursor:pointer;padding:2px 0;list-style:revert;';
                            detailsSummary.textContent = 'Full record details';
                            detailsEl.appendChild(detailsSummary);

                            const detailsBody = document.createElement('div');
                            detailsBody.style.cssText = 'font-size:0.79rem;padding:4px 2px;margin-top:4px;max-height:350px;overflow-y:auto;line-height:1.6;';
                            if (record.record_id) {
                                const idRow = document.createElement('div');
                                idRow.style.cssText = 'font-size:0.74rem;color:#9ca3af;margin-bottom:2px;';
                                idRow.textContent = `Record ID: ${record.record_id}`;
                                detailsBody.appendChild(idRow);
                            }
                            const sectionDefs = [['Tumor', record.tumor], ['Marker / Biomarker', record.marker], ['Testing', record.testing], ['Therapy', record.therapy], ['Guideline Source', record.guideline_metadata]];
                            sectionDefs.forEach(([title, obj]) => renderSection(title, obj, detailsBody));
                            const skipKeys = new Set(['record_id', 'tumor', 'marker', 'testing', 'therapy', 'guideline_metadata', 'dataset_file', 'dataset_record_index', 'dataset_name']);
                            Object.keys(record).filter(k => !skipKeys.has(k)).forEach(k => {
                                const val = record[k];
                                if (val && typeof val === 'object') {
                                    renderSection(k.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase()), val, detailsBody);
                                } else if (val !== null && val !== undefined && val !== '') {
                                    const row = document.createElement('div');
                                    row.style.cssText = 'margin:2px 0;';
                                    const labelEl = document.createElement('span');
                                    labelEl.style.cssText = 'font-weight:600;color:#374151;';
                                    labelEl.textContent = k.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase()) + ': ';
                                    row.appendChild(labelEl);
                                    row.appendChild(document.createTextNode(String(val)));
                                    detailsBody.appendChild(row);
                                }
                            });
                            detailsEl.appendChild(detailsBody);
                            recordEl.appendChild(detailsEl);

                            guidelinesResults.appendChild(recordEl);
                        });

                        const sourceNote = document.createElement('div');
                        sourceNote.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
                        sourceNote.textContent = 'Source: drdoubleb.com/guidelines. For reference only — verify against current published guidelines.';
                        guidelinesResults.appendChild(sourceNote);

                    } catch (err) {
                        guidelinesResults.innerHTML = `<div style="font-size:0.85rem;color:#9ca3af;">Guidelines data unavailable: ${err.message}</div>`;
                        aiReviewExtras.guidelines = { error: err.message, cancer_type: selectedCancer, gene };
                    }
                });
            }

            // Optional AI review card (manual trigger; sends current annotation context to OpenRouter via backend proxy).
            {
                const aiCard = document.createElement('div');
                aiCard.className = 'card ai-review-card';
                const aiTitle = document.createElement('h3');
                aiTitle.textContent = 'Optional AI Review';
                applyCardTheme(aiCard, 'Optional AI Review');
                aiCard.appendChild(aiTitle);

                const aiContent = document.createElement('div');
                aiContent.className = 'card-content ai-review-content';

                const aiIntro = document.createElement('p');
                aiIntro.className = 'ai-review-intro';
                aiIntro.textContent = 'Send the retrieved variant data to OpenRouter for a structured draft interpretation. No request is sent until you click Run AI review.';
                aiContent.appendChild(aiIntro);

                const controls = document.createElement('div');
                controls.className = 'ai-review-controls';

                const modelLabel = document.createElement('label');
                modelLabel.textContent = 'Model';
                modelLabel.setAttribute('for', 'aiReviewModelSelect');
                controls.appendChild(modelLabel);

                const modelSelect = document.createElement('select');
                modelSelect.id = 'aiReviewModelSelect';
                OPENROUTER_MODEL_OPTIONS.forEach((modelName) => {
                    const opt = document.createElement('option');
                    opt.value = modelName;
                    opt.textContent = modelName;
                    modelSelect.appendChild(opt);
                });
                controls.appendChild(modelSelect);

                const runButton = document.createElement('button');
                runButton.type = 'button';
                runButton.textContent = 'Run AI review';
                controls.appendChild(runButton);

                aiContent.appendChild(controls);

                const notesWrap = document.createElement('div');
                notesWrap.className = 'ai-review-notes';
                const notesLabel = document.createElement('label');
                notesLabel.setAttribute('for', 'aiReviewNotes');
                notesLabel.textContent = 'Any additional notes for AI review';
                notesWrap.appendChild(notesLabel);
                const notesInput = document.createElement('textarea');
                notesInput.id = 'aiReviewNotes';
                notesInput.rows = 3;
                notesInput.placeholder = 'Optional — extra clinical context, prior therapies, specific questions, etc.';
                notesWrap.appendChild(notesInput);
                aiContent.appendChild(notesWrap);

                const aiContextInspector = document.createElement('details');
                aiContextInspector.style.cssText = 'margin:6px 0 2px;font-size:0.80rem;';
                const aiContextInspectorSummary = document.createElement('summary');
                aiContextInspectorSummary.style.cssText = 'cursor:pointer;color:#9ca3af;padding:2px 0;list-style:revert;font-size:0.79rem;';
                aiContextInspectorSummary.textContent = 'Context sent to AI (populated after run)';
                aiContextInspector.appendChild(aiContextInspectorSummary);
                const aiContextPre = document.createElement('pre');
                aiContextPre.style.cssText = 'font-size:0.73rem;white-space:pre-wrap;word-break:break-word;background:#f8fafc;border:1px solid #e5e7eb;padding:8px;border-radius:4px;margin-top:4px;max-height:400px;overflow-y:auto;color:#374151;';
                aiContextPre.textContent = 'Run AI review to populate this section.';
                aiContextInspector.appendChild(aiContextPre);
                aiContent.appendChild(aiContextInspector);

                const aiOutput = document.createElement('div');
                aiOutput.className = 'ai-review-output';
                aiContent.appendChild(aiOutput);
                aiCard.appendChild(aiContent);
                cardsContainer.appendChild(aiCard);

                runButton.addEventListener('click', async () => {
                    runButton.disabled = true;
                    const previousText = runButton.textContent;
                    runButton.textContent = 'Running…';
                    aiOutput.innerHTML = '<div class="ai-review-loading">Gathering FDA, trial, and annotation context for AI review…</div>';
                    try {
                        const [fdaRecords, clinicalTrialData, supplementalContext] = await Promise.all([
                            aiReviewGene ? fetchFdaCompanionDiagnostics(aiReviewGene).catch(() => []) : Promise.resolve([]),
                            aiReviewGene ? fetchClinicalTrials(aiReviewGene, tumorType).catch(() => ({ total: 0, studies: [] })) : Promise.resolve({ total: 0, studies: [] }),
                            fetchAiReviewSupplementalContext().catch((err) => ({ error: err.message || 'Supplemental context unavailable' }))
                        ]);
                        const userNotes = (notesInput.value || '').trim();
                        const aiContext = {
                            submitted_variant: rawInput,
                            normalized_genomic_variant: gVariant,
                            tumor_type: tumorType,
                            user_notes: userNotes || undefined,
                            gene: aiReviewGene,
                            genes: geneNames,
                            selected_variant_term: aiReviewSearchVariantTerm,
                            cdna: aiReviewCdna,
                            protein: aiReviewProtein,
                            summary_rows: summaryRows,
                            details: detailsData,
                            transcripts: transcriptsList,
                            fda_companion_diagnostics_records: fdaRecords,
                            clinical_trials: clinicalTrialData,
                            supplemental_card_data: supplementalContext,
                            myvariant_annotation: (() => {
                                if (!annotation) return null;
                                const { clinvar, civic, gnomad_exome, dbsnp } = annotation;
                                // Use dedicated API results when available; fall back to myvariant.info
                                // only when the dedicated call failed or returned no data.
                                const clinvarOk = supplementalContext?.clinvar_variant_record && !supplementalContext.clinvar_variant_record.error;
                                const civicOk = supplementalContext?.civic_api?.gene != null;
                                const gnomadOk = supplementalContext?.gnomad_v4?.status === 'found';
                                const trimmed = {};
                                if (!clinvarOk && clinvar) trimmed.clinvar = clinvar;
                                if (!civicOk && civic) trimmed.civic = civic;
                                if (!gnomadOk && gnomad_exome) trimmed.gnomad_exome = gnomad_exome;
                                if (dbsnp) {
                                    // Drop noisy fields that don't help variant interpretation:
                                    // - gene.rnas: ~15 transcript entries that all report the reference
                                    //   SPDI (e.g. c.818=), not the variant itself.
                                    // - citations: bare PMID list with no titles/abstracts.
                                    const { citations: _dbsnpCitations, gene: dbsnpGene, ...dbsnpRest } = dbsnp;
                                    let trimmedGene = dbsnpGene;
                                    if (dbsnpGene && typeof dbsnpGene === 'object') {
                                        const { rnas: _dbsnpRnas, ...geneRest } = dbsnpGene;
                                        trimmedGene = geneRest;
                                    }
                                    trimmed.dbsnp = trimmedGene !== undefined
                                        ? { ...dbsnpRest, gene: trimmedGene }
                                        : dbsnpRest;
                                }
                                return Object.keys(trimmed).length > 0 ? trimmed : null;
                            })(),
                            ensembl_recoder: typeof recoderData !== 'undefined' ? recoderData : null
                        };
                        aiContextPre.textContent = JSON.stringify(aiContext, null, 2);
                        aiContextInspectorSummary.textContent = 'Context sent to AI';
                        const data = await fetchAiReview(aiContext, modelSelect.value);
                        renderAiReview(data.review, aiOutput);
                    } catch (err) {
                        aiOutput.innerHTML = '';
                        const errorEl = document.createElement('div');
                        errorEl.className = 'ai-review-error';
                        errorEl.textContent = `AI review unavailable: ${err.message}`;
                        aiOutput.appendChild(errorEl);
                    } finally {
                        runButton.disabled = false;
                        runButton.textContent = previousText;
                    }
                });
            }

            // Show cards and hide legacy tables for a cleaner view
            cardsContainer.classList.remove('hidden');
            summaryTable.classList.add('hidden');
            detailsContainer.classList.add('hidden');
            // Raw output and result section visibility
            rawOutput.textContent = JSON.stringify(annotation, null, 2);
            // Populate Ensembl output with the raw variant recoder response. If none exists (e.g. recoder failed), clear it.
            if (typeof recoderData !== 'undefined' && recoderData) {
                try {
                    ensemblOutput.textContent = JSON.stringify(recoderData, null, 2);
                } catch {
                    ensemblOutput.textContent = '';
                }
            } else {
                ensemblOutput.textContent = '';
            }
            resultSection.classList.remove('hidden');
        } catch (err) {
            statusEl.textContent = 'Error: ' + err.message;
            console.error(err);
        }
    });

    // Auto-run lookup when the page is opened with ?variant=<value> and/or ?tumorType=<value>.
    // This allows direct linking from external tools while preserving the existing normalization/search logic.
    const initialVariant = getVariantFromUrl();
    const initialTumorType = getTumorTypeFromUrl();
    if (initialTumorType && initialTumorType.trim() && tumorTypeInput) {
        tumorTypeInput.value = initialTumorType.trim();
    }
    if (initialVariant && initialVariant.trim()) {
        input.value = initialVariant.trim();
        if (typeof form.requestSubmit === 'function') {
            form.requestSubmit();
        } else {
            form.dispatchEvent(new Event('submit', { cancelable: true, bubbles: true }));
        }
    }
});
