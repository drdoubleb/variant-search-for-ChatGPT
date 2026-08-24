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


// RefSeq chromosome accession VERSIONS distinguish the assembly: GRCh38 bumped
// every chromosome's version by one (chr7 is NC_000007.13 in GRCh37 and .14 in
// GRCh38). This lets recoder/SPDI output declare its own assembly, so the
// pipeline lifts GRCh38 coordinates and — critically — leaves GRCh37 ones
// alone. Blanket-lifting everything silently "mapped" already-hg19 positions
// hundreds of kb away (Ensembl's /map does not error when handed a coordinate
// from the wrong assembly; it just returns a wrong answer).
const GRCH37_NC_VERSIONS = {
    1: 10, 2: 11, 3: 11, 4: 11, 5: 9, 6: 11, 7: 13, 8: 10, 9: 11, 10: 10,
    11: 9, 12: 11, 13: 10, 14: 8, 15: 9, 16: 9, 17: 10, 18: 9, 19: 9, 20: 10,
    21: 8, 22: 10, 23: 10, 24: 9
};

// Returns 'GRCh37' | 'GRCh38' | null for an NC_-prefixed string (a full hgvsg
// or SPDI is fine — only the leading accession is read). Null when the
// accession is unversioned, mitochondrial, or not recognisably either build.
function assemblyFromNcAccession(value) {
    const m = String(value || '').match(/^NC_(\d{6})\.(\d+)/);
    if (!m) return null;
    const num = parseInt(m[1], 10);
    const version = parseInt(m[2], 10);
    const v37 = GRCH37_NC_VERSIONS[num];
    if (!v37) return null;
    if (version === v37) return 'GRCh37';
    if (version === v37 + 1) return 'GRCh38';
    return null;
}

const DEFAULT_BACKEND_API_BASE_URL = 'https://variant-search-for-chat-gpt.vercel.app';
// The canonical production Vercel host. Any *other* *.vercel.app host is a
// preview/branch deployment.
const PRODUCTION_VERCEL_HOST = 'variant-search-for-chat-gpt.vercel.app';

function trimTrailingSlash(value) {
    return String(value || '').trim().replace(/\/+$/, '');
}

// True when the page is served from a Vercel preview/branch deployment
// (a *.vercel.app host other than the pinned production host).
function isVercelPreviewHost() {
    const host = (typeof window !== 'undefined' && window.location && window.location.hostname) || '';
    return /\.vercel\.app$/i.test(host) && host !== PRODUCTION_VERCEL_HOST;
}

function getBackendApiBaseUrl() {
    // An explicit config wins — including an empty string, which means "call the
    // same origin". index.html normally sets this global.
    if (typeof window.BACKEND_API_BASE_URL === 'string') {
        return trimTrailingSlash(window.BACKEND_API_BASE_URL);
    }
    // No explicit config: on a Vercel preview/branch deployment prefer the
    // same-origin API so the preview frontend calls its own serverless functions
    // (which match the loaded build) instead of the pinned production backend.
    // Everywhere else (Hostinger, custom domains, the production apex) use the
    // pinned production backend so no local /api deploy is required.
    if (isVercelPreviewHost()) return '';
    return DEFAULT_BACKEND_API_BASE_URL;
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
// Set when a protein-level query could only be resolved as far as its codon
// (frameshifts and multi-residue delins do not determine a nucleotide change).
// Carries the residue range so the Variant card can say what was and was not
// resolved. Reset at the start of every lookup.
let codonOnlyResolutionGlobal = null;

// Keep third-party API latency from blocking the UI for long periods.
const API_TIMEOUT_MS = {
    // The Ensembl GRCh37 mirror answers a warm query in well under a second but
    // has been measured taking 40 s+ on a cold VEP query and returning 5xx in
    // bursts. These deadlines are per attempt; fetchWithRetry adds the retries.
    myvariant: 10000,
    recoder: 12000,
    vep: 20000,
    lookup: 4000,
    liftover: 4000,
    cosmic: 2500,
    cosmicMeta: 1500,
    tp53: 15000,
    clinvar: 15000,
    civic: 8000,
    // The /api/pubmed function is capped at maxDuration 15 in vercel.json, so
    // waiting longer than that client-side only delays the inevitable error.
    pubmed: 16000,
    fda: 8000,
    bbkb: 8000,
    gnomadV4: 10000,
    clinicalTrials: 20000,
    spliceai: 25000,
    aiReview: 60000,
    cbioportal: 15000,
    openFda: 12000,
    ensemblSequence: 8000,
    ensemblLookup: 10000
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

// Turn a thrown fetch/HTTP error into a short phrase fit for the progress panel
// and for user-facing error text: "HTTP 503", "timed out", "network error".
function describeUpstreamFailure(err) {
    if (!err) return 'failed';
    if (err.status) return `HTTP ${err.status}`;
    const msg = String(err.message || err);
    if (/timed out after/i.test(msg)) return 'timed out';
    if (/Failed to fetch|NetworkError|network error|load failed/i.test(msg)) return 'network error';
    const m = msg.match(/\((\d{3})\)/);
    if (m) return `HTTP ${m[1]}`;
    return msg.length > 60 ? `${msg.slice(0, 57)}…` : msg;
}

/*
 * Build the error shown when a lookup resolves to nothing.
 *
 * The old message was a single fixed string — "Variant not found. Please verify
 * the genomic coordinate and reference allele." — printed no matter the cause.
 * A GRCh37 Ensembl outage, a MyVariant timeout and a genuinely nonexistent
 * variant all rendered as "your coordinate is wrong", which sends users to
 * re-check input that was never the problem. Distinguish the two cases using
 * the source health the progress panel recorded.
 */
function buildLookupFailureError() {
    const failed = LookupProgress.failedSources();
    if (LookupProgress.sawUpstreamOutage()) {
        const err = new Error(
            `Lookup could not complete: ${failed.join(', ')} did not respond. `
            + 'This is an upstream outage, not necessarily a problem with your variant — please retry in a moment.'
        );
        err.upstreamOutage = true;
        return err;
    }
    if (failed.length) {
        return new Error(
            `Lookup could not complete: ${failed.join(', ')} failed. Please retry in a moment.`
        );
    }
    return new Error(
        'No annotation found for this variant. The coordinate resolved, but none of the '
        + 'annotation sources hold a record for it — check the reference allele and that the '
        + 'position is GRCh37/hg19, which is the build this tool uses.'
    );
}

// Collapse a comma-separated symbol list to its unique members, preserving order.
function dedupeSymbolList(value) {
    const seen = new Set();
    return String(value || '')
        .split(',')
        .map((t) => t.trim())
        .filter((t) => t && !seen.has(t.toUpperCase()) && seen.add(t.toUpperCase()))
        .join(', ');
}

// selectCanonicalTranscript scores candidates and can throw on odd shapes; the
// progress panel only needs a best-effort index, never an exception.
function selectCanonicalTranscriptSafe(transcripts) {
    try {
        const idx = selectCanonicalTranscript(
            transcripts.map((t) => ({ transcript: t.transcript, cDNA: t.cDNA, protein: t.protein, source: 'root' })),
            typeof targetProtGlobal !== 'undefined' ? targetProtGlobal : null
        );
        return typeof idx === 'number' && idx >= 0 && idx < transcripts.length ? idx : 0;
    } catch {
        return 0;
    }
}

// HTTP statuses that mean "the server is unwell", not "your request is wrong".
// 429 is included because Ensembl rate-limits bursts and expects a retry.
const RETRYABLE_HTTP_STATUS = new Set([408, 425, 429, 500, 502, 503, 504]);

// True when a thrown fetch error is a transport-level failure (DNS, reset,
// TLS, our own AbortController timeout) rather than a rejected response.
function isRetryableFetchError(err) {
    const msg = String((err && err.message) || err || '');
    return /timed out after|Failed to fetch|NetworkError|network error|load failed|ECONNRESET|terminated/i.test(msg);
}

/*
 * fetchWithTimeout plus bounded retries with exponential backoff and jitter.
 *
 * The public annotation APIs this app depends on — grch37.rest.ensembl.org
 * above all — intermittently return 500/503 and occasionally take tens of
 * seconds to answer a cold query. A single attempt behind a 6-7 s deadline
 * turns those blips into "variant not found", which is both wrong and
 * unactionable. Retrying transport errors and 5xx costs a few seconds in the
 * bad case and rescues the lookup outright in the common one.
 *
 * `onAttempt` reports each attempt to the progress panel so a slow lookup
 * shows what it is waiting on instead of sitting on a frozen status line.
 * The final failure is re-thrown with `.retryable` set, which lets the caller
 * tell "the service is down" apart from "this variant does not exist".
 */
async function fetchWithRetry(url, options = {}, timeoutMs = 6000, retryOpts = {}) {
    const attempts = Math.max(1, retryOpts.attempts ?? 3);
    const baseDelayMs = retryOpts.baseDelayMs ?? 400;
    const onAttempt = typeof retryOpts.onAttempt === 'function' ? retryOpts.onAttempt : null;
    let lastErr = null;
    for (let attempt = 1; attempt <= attempts; attempt++) {
        if (onAttempt) onAttempt({ attempt, attempts });
        try {
            const res = await fetchWithTimeout(url, options, timeoutMs);
            if (RETRYABLE_HTTP_STATUS.has(res.status) && attempt < attempts) {
                lastErr = new Error(`Upstream returned ${res.status}`);
                lastErr.retryable = true;
                lastErr.status = res.status;
            } else {
                return res;
            }
        } catch (err) {
            if (!isRetryableFetchError(err) || attempt === attempts) {
                if (isRetryableFetchError(err)) err.retryable = true;
                throw err;
            }
            lastErr = err;
            lastErr.retryable = true;
        }
        // Exponential backoff with jitter so simultaneous cards do not retry in lockstep.
        const delay = baseDelayMs * Math.pow(2, attempt - 1) * (0.7 + Math.random() * 0.6);
        await new Promise((resolve) => setTimeout(resolve, delay));
    }
    throw lastErr || new Error('Request failed');
}

/*
 * ── Lookup progress panel ───────────────────────────────────────────────────
 *
 * The lookup used to communicate through a single mutable line of text
 * ("Processing…", "Converting variant…"), which hid two things that matter:
 * what the tool actually resolved the input to, and which upstream source
 * failed when nothing came back. A variant that is missing from MyVariant.info
 * looked exactly like an Ensembl outage, and both looked like a typo.
 *
 * LookupProgress renders a compact pipeline instead: the resolved nomenclature
 * up top, then one row per stage with its state, timing and retry count. It
 * also records source health so the final error message can say "Ensembl is
 * returning 503" rather than blaming the user's coordinate.
 */
const LookupProgress = (() => {
    const STEP_LABELS = {
        parse: 'Parse input',
        recoder: 'Resolve nomenclature · Ensembl',
        myvariant: 'Annotations · MyVariant.info',
        vep: 'Consequences · Ensembl VEP',
        evidence: 'Clinical evidence sources'
    };
    const STATE_ICON = { pending: '○', running: '◐', ok: '●', empty: '◍', warn: '▲', fail: '✕' };

    let root = null;
    let stepsEl = null;
    let headerEl = null;
    let bodyEl = null;
    let resolvedEl = null;
    const steps = new Map();
    let startedAt = 0;
    // Frozen at finish() so re-rendering the header (e.g. on collapse) reports
    // how long the lookup took, not how long ago it happened.
    let finishedAt = 0;
    let active = false;
    let collapsed = false;
    // Kept so the collapsed header can still name what was resolved.
    let lastResolved = null;
    let lastHeader = { state: 'running', text: 'Looking up variant…' };

    function ensureDom() {
        if (root && document.body.contains(root)) return;
        const status = document.getElementById('status');
        if (!status || !status.parentNode) return;
        root = document.createElement('section');
        root.id = 'lookupProgress';
        root.className = 'lookup-progress hidden';
        root.setAttribute('aria-live', 'polite');

        // A button, not a div: the header is the collapse control, so it has to
        // be reachable by keyboard and announce its expanded state.
        headerEl = document.createElement('button');
        headerEl.type = 'button';
        headerEl.className = 'lp-header';
        headerEl.setAttribute('aria-controls', 'lookupProgressBody');
        headerEl.addEventListener('click', () => setCollapsed(!collapsed));
        root.appendChild(headerEl);

        bodyEl = document.createElement('div');
        bodyEl.className = 'lp-body';
        bodyEl.id = 'lookupProgressBody';
        root.appendChild(bodyEl);

        resolvedEl = document.createElement('div');
        resolvedEl.className = 'lp-resolved hidden';
        bodyEl.appendChild(resolvedEl);

        stepsEl = document.createElement('ol');
        stepsEl.className = 'lp-steps';
        bodyEl.appendChild(stepsEl);

        status.parentNode.insertBefore(root, status);
    }

    function setCollapsed(next) {
        collapsed = Boolean(next);
        if (!root) return;
        root.classList.toggle('lp-collapsed', collapsed);
        if (bodyEl) bodyEl.classList.toggle('hidden', collapsed);
        renderHeader(lastHeader.state, lastHeader.text);
    }

    // One-line stand-in for the panel body while collapsed: what was resolved,
    // and how many sources answered.
    function collapsedSummary() {
        const parts = [];
        if (lastResolved) {
            const primary = lastResolved.genomic || lastResolved.cdna || lastResolved.protein;
            if (primary) parts.push(primary);
            else if (lastResolved.gene) parts.push(lastResolved.gene);
        }
        const checked = steps.size;
        if (checked) parts.push(`${checked} source${checked === 1 ? '' : 's'} checked`);
        return parts.join(' · ');
    }

    function renderHeader(state, text) {
        if (!headerEl) return;
        lastHeader = { state, text };
        headerEl.className = `lp-header lp-${state}`;
        headerEl.setAttribute('aria-expanded', String(!collapsed));
        headerEl.textContent = '';

        const dot = document.createElement('span');
        dot.className = 'lp-spinner';
        dot.setAttribute('aria-hidden', 'true');
        const label = document.createElement('span');
        label.className = 'lp-header-text';
        label.textContent = text;
        headerEl.appendChild(dot);
        headerEl.appendChild(label);

        // Collapsed, the header is the whole panel, so it carries the summary.
        if (collapsed) {
            const summary = collapsedSummary();
            if (summary) {
                const sum = document.createElement('span');
                sum.className = 'lp-summary';
                sum.textContent = summary;
                headerEl.appendChild(sum);
            }
        }

        if (startedAt) {
            const t = document.createElement('span');
            t.className = 'lp-elapsed';
            t.textContent = `${(((finishedAt || Date.now()) - startedAt) / 1000).toFixed(1)}s`;
            headerEl.appendChild(t);
        }

        const chevron = document.createElement('span');
        chevron.className = 'lp-chevron';
        chevron.setAttribute('aria-hidden', 'true');
        chevron.textContent = collapsed ? '▸' : '▾';
        headerEl.appendChild(chevron);
    }

    function renderStep(key) {
        const step = steps.get(key);
        if (!step || !stepsEl) return;
        let li = stepsEl.querySelector(`[data-step="${key}"]`);
        if (!li) {
            li = document.createElement('li');
            li.className = 'lp-step';
            li.dataset.step = key;
            stepsEl.appendChild(li);
        }
        li.className = `lp-step lp-${step.state}`;
        li.textContent = '';

        const icon = document.createElement('span');
        icon.className = 'lp-icon';
        icon.setAttribute('aria-hidden', 'true');
        icon.textContent = STATE_ICON[step.state] || STATE_ICON.pending;

        const name = document.createElement('span');
        name.className = 'lp-name';
        name.textContent = STEP_LABELS[key] || key;

        const detail = document.createElement('span');
        detail.className = 'lp-detail';
        let detailText = step.detail || '';
        if (step.state === 'running' && step.attempt > 1) {
            detailText = `retry ${step.attempt} of ${step.attempts}…`;
        }
        detail.textContent = detailText;

        const timing = document.createElement('span');
        timing.className = 'lp-timing';
        if (step.startedAt) {
            const endedAt = step.endedAt || Date.now();
            timing.textContent = `${((endedAt - step.startedAt) / 1000).toFixed(1)}s`;
        }

        li.appendChild(icon);
        li.appendChild(name);
        li.appendChild(detail);
        li.appendChild(timing);
    }

    return {
        start() {
            ensureDom();
            steps.clear();
            startedAt = Date.now();
            finishedAt = 0;
            active = true;
            lastResolved = null;
            if (!root) return;
            if (stepsEl) stepsEl.textContent = '';
            if (resolvedEl) { resolvedEl.textContent = ''; resolvedEl.classList.add('hidden'); }
            root.classList.remove('hidden');
            // A new lookup always starts expanded, whatever the last one left behind.
            setCollapsed(false);
            renderHeader('running', 'Looking up variant…');
        },
        // state: running | ok | empty | warn | fail
        step(key, state, detail) {
            if (!active) return;
            ensureDom();
            const prev = steps.get(key) || {};
            steps.set(key, {
                state,
                detail: detail || '',
                attempt: prev.attempt || 1,
                attempts: prev.attempts || 1,
                // Keep the original start so the timing shown spans every retry.
                startedAt: prev.startedAt || Date.now(),
                endedAt: state === 'running' ? null : Date.now(),
                upstream: prev.upstream || false
            });
            renderStep(key);
        },
        // Do not overwrite a recorded failure with a softer state: a helper that
        // swallows its own error and returns an empty list must not be reported
        // as "no data" when the real answer is "the service is down".
        stepUnlessFailed(key, state, detail) {
            const prev = steps.get(key);
            if (prev && prev.state === 'fail') return;
            this.step(key, state, detail);
        },
        stateOf(key) {
            const step = steps.get(key);
            return step ? step.state : null;
        },
        attempt(key, attempt, attempts) {
            if (!active) return;
            const prev = steps.get(key);
            if (!prev) return;
            prev.attempt = attempt;
            prev.attempts = attempts;
            steps.set(key, prev);
            renderStep(key);
        },
        // Show what the tool actually resolved the raw input to.
        resolved({ input, genomic, cdna, protein, gene, consequence, assembly }) {
            if (!active) return;
            ensureDom();
            lastResolved = { genomic, cdna, protein, gene };
            if (!resolvedEl) return;
            resolvedEl.textContent = '';
            const rows = [
                ['You entered', input],
                ['Gene', gene],
                [`Genomic (${assembly || 'GRCh37'})`, genomic],
                ['Coding', cdna],
                ['Protein', protein],
                ['Effect', consequence]
            ].filter(([, v]) => v);
            if (!rows.length) return;
            for (const [label, value] of rows) {
                const row = document.createElement('div');
                row.className = 'lp-resolved-row';
                const k = document.createElement('span');
                k.className = 'lp-key';
                k.textContent = label;
                const v = document.createElement('span');
                v.className = 'lp-value';
                v.textContent = value;
                row.appendChild(k);
                row.appendChild(v);
                resolvedEl.appendChild(row);
            }
            resolvedEl.classList.remove('hidden');
        },
        finish(state, text) {
            if (!active) return;
            finishedAt = Date.now();
            for (const [key, step] of steps.entries()) {
                if (step.state === 'running') {
                    step.state = state === 'fail' ? 'fail' : 'warn';
                    step.endedAt = Date.now();
                    steps.set(key, step);
                    renderStep(key);
                }
            }
            // On success the variant summary below is what the user came for, so
            // fold the panel down to one line and let them reopen it. A failed or
            // degraded lookup stays open: the per-source detail *is* the answer.
            setCollapsed(state === 'ok');
            renderHeader(state, text);
            active = false;
        },
        hide() {
            active = false;
            if (root) root.classList.add('hidden');
        },
        // Which upstream sources failed outright (as opposed to returning no data).
        failedSources() {
            const out = [];
            for (const [key, step] of steps.entries()) {
                if (step.state === 'fail') out.push(STEP_LABELS[key] || key);
            }
            return out;
        },
        sawUpstreamOutage() {
            for (const step of steps.values()) {
                if (step.state === 'fail' && step.upstream) return true;
            }
            return false;
        },
        markUpstreamOutage(key) {
            const step = steps.get(key);
            if (step) { step.upstream = true; steps.set(key, step); }
        }
    };
})();

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


// Parse a genomic HGVS string into a structural descriptor.
//
// Covers every g. form the app can produce or receive, not just substitutions:
//   chr7:g.140453136A>T            → { type: 'sub',    start: 140453136, end: 140453136 }
//   chr16:g.2122948A>-             → { type: 'del',    start: 2122948,   end: 2122948   }
//   chr16:g.2122948_2122950del     → { type: 'del',    start: 2122948,   end: 2122950   }
//   chr16:g.2122948_2122950delAAT  → { type: 'del', refSeq: 'AAT' }
//   chr7:g.55242465_55242479del    → { type: 'del' }
//   chr2:g.29443695_29443696dup    → { type: 'dup' }
//   chr4:g.55593603_55593604insTT  → { type: 'ins' }
//   chr7:g.140453136_140453137delinsAA → { type: 'delins', altSeq: 'AA' }
//   chr3:g.10183694_10183700inv    → { type: 'inv' }
//
// Returns null when the string is not a recognisable genomic HGVS descriptor.
// `refSeq`/`altSeq` are only populated when the notation spells the bases out;
// callers that need the actual alleles fill the gaps from a reference sequence
// (see resolveVcfAlleles) rather than assuming they are present.
function parseGenomicHgvs(gVariant) {
    if (!gVariant) return null;
    const m = String(gVariant).trim().match(/^chr([0-9XYMT]+):g\.(\d+)(?:_(\d+))?(.*)$/i);
    if (!m) return null;
    const chrom = `chr${m[1].toUpperCase()}`;
    const start = Number(m[2]);
    const explicitEnd = m[3] ? Number(m[3]) : null;
    const suffix = m[4].trim();
    const mk = (type, end, refSeq, altSeq) => ({
        chrom,
        start,
        end: Number.isFinite(end) ? end : start,
        type,
        refSeq: refSeq || null,
        altSeq: altSeq || null
    });

    // A bare span with no event ("chr16:g.2136834_2136836") is a region, not a
    // variant. Protein-only queries that cannot be resolved to a nucleotide change
    // are represented this way so position-based cards still work while the
    // allele-dependent ones correctly stay empty.
    if (suffix === '') return mk('region', explicitEnd ?? start, null, null);

    // Substitution — including the "REF>-" / "->ALT" forms that the SPDI converter
    // emits for single-base deletions and insertions.
    let s = suffix.match(/^([A-Za-z-]+)>([A-Za-z-]+)$/);
    if (s) {
        const ref = s[1].toUpperCase();
        const alt = s[2].toUpperCase();
        if (ref === '-' && alt !== '-') return mk('ins', start, null, alt);
        if (alt === '-' && ref !== '-') return mk('del', start + ref.length - 1, ref, null);
        if (ref === '-' || alt === '-') return null;
        return mk(ref.length === 1 && alt.length === 1 ? 'sub' : 'delins', explicitEnd ?? (start + ref.length - 1), ref, alt);
    }
    // delins — must be tested before plain "del" so "delinsAA" is not read as a deletion.
    s = suffix.match(/^del([A-Za-z]*)ins([A-Za-z]+)$/i);
    if (s) return mk('delins', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, s[2].toUpperCase());
    s = suffix.match(/^del([A-Za-z]*)$/i);
    if (s) return mk('del', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, null);
    s = suffix.match(/^dup([A-Za-z]*)$/i);
    if (s) return mk('dup', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, null);
    s = suffix.match(/^ins([A-Za-z]+)$/i);
    if (s) return mk('ins', explicitEnd ?? start, null, s[1].toUpperCase());
    s = suffix.match(/^inv([A-Za-z]*)$/i);
    if (s) return mk('inv', explicitEnd ?? start, s[1] ? s[1].toUpperCase() : null, null);
    return null;
}

/*
 * ── Coordinate-row parsing ──────────────────────────────────────────────────
 *
 * Users paste variants straight out of spreadsheets, VCFs, MAFs and pipeline
 * reports, so a "row" arrives in many shapes:
 *
 *   X  20148674  EIF1AX  T  TG          tab- or space-separated, gene in the middle
 *   EIF1AX  X  20148674  T  TG          MAF column order (Hugo_Symbol first)
 *   X,20148674,EIF1AX,T,TG              CSV paste
 *   X  20148674  EIF1AX  T  TG  hg19    trailing build//extra columns
 *   X  20148674  -  TG                  MAF-style insertion ("-" reference)
 *
 * Only the first shape used to parse. The rest fell through to the generic
 * "GENE p.Change" branch at the end of normalizeVariantInput, which turned
 * "EIF1AX X 20148674 T TG" into the nonsense query "EIF1AX:p.X" and then
 * reported it as a variant that could not be found. Parsing the row properly —
 * and refusing to guess when it is genuinely ambiguous — is the difference
 * between a working lookup and a misleading error.
 */

// Tokens that carry no allele information and only ever confuse the parser.
const COORDINATE_ROW_NOISE = /^(hg18|hg19|hg38|grch3[678]|b3[678]|chr|snp|snv|ins|del|indel|mnp|somatic|germline|het|hom|\+|\.)$/i;

const isDnaAllele = (t) => /^(?:[ACGTN]+|-)$/i.test(String(t));
const isChromToken = (t) => /^(?:chr)?(?:[0-9]{1,2}|X|Y|M|MT)$/i.test(String(t));

// Split a pasted row on any plausible column separator: whitespace (spaces,
// tabs, newlines), commas, semicolons or pipes. Surrounding quotes are stripped
// because spreadsheet exports like to add them.
function splitCoordinateRow(raw) {
    return String(raw || '')
        .trim()
        // Strip thousands separators inside a coordinate ("20,148,674") before
        // commas are treated as column separators. The grouping pattern is
        // required so a genuine CSV row ("1,12345,A,T") is not welded together.
        .replace(/\b\d{1,3}(?:,\d{3})+(?!\d)/g, (m) => m.replace(/,/g, ''))
        .split(/[\s,;|]+/)
        .map((t) => t.replace(/^["']+|["']+$/g, ''))
        .filter(Boolean);
}

// True when the input looks like a pasted coordinate row rather than HGVS.
// Used to refuse the "GENE p.Change" fallback for rows we could not parse,
// so an unparseable row reports what is wrong instead of querying garbage.
function looksLikeCoordinateRow(raw) {
    const toks = splitCoordinateRow(raw).filter((t) => !COORDINATE_ROW_NOISE.test(t));
    if (toks.length < 4) return false;
    if (/:[gcp]\./i.test(String(raw))) return false;
    for (let i = 0; i < toks.length - 1; i++) {
        if (isChromToken(toks[i]) && /^\d[\d,]*$/.test(toks[i + 1])) return true;
    }
    return false;
}

/*
 * Parse a pasted coordinate row into { chrom, pos, ref, alt, gene }.
 *
 * Strategy: anchor on the chromosome/position pair wherever it sits in the row,
 * then read the alleles out of the tokens that follow. When more than two
 * DNA-like tokens follow the position the last two win, because a gene symbol
 * that happens to spell DNA (TTN, CAT, TAT) is always written before the
 * alleles in every layout above. Returns null when nothing matches, so callers
 * can fall through rather than act on a guess.
 */
function parseCoordinateRow(raw) {
    const all = splitCoordinateRow(raw);
    if (all.length < 4) return null;
    const toks = all.filter((t) => !COORDINATE_ROW_NOISE.test(t));
    if (toks.length < 4) return null;

    // Anchor: a chromosome token immediately followed by a numeric position.
    let chromIdx = -1;
    for (let i = 0; i < toks.length - 1; i++) {
        if (isChromToken(toks[i]) && /^\d[\d,]*$/.test(toks[i + 1])) { chromIdx = i; break; }
    }
    if (chromIdx === -1) return null;

    const chrom = toks[chromIdx].replace(/^chr/i, '').toUpperCase();
    const pos = toks[chromIdx + 1].replace(/,/g, '');
    if (!/^\d+$/.test(pos)) return null;

    const after = toks.slice(chromIdx + 2);
    const dnaIdx = [];
    after.forEach((t, i) => { if (isDnaAllele(t)) dnaIdx.push(i); });
    if (dnaIdx.length < 2) return null;
    // Last two DNA-like tokens are the alleles (see note above).
    const refIdx = dnaIdx[dnaIdx.length - 2];
    const altIdx = dnaIdx[dnaIdx.length - 1];
    // Alleles are adjacent columns in every supported layout; anything else is
    // too ambiguous to guess at.
    if (altIdx !== refIdx + 1) return null;

    // "-" and "." are the MAF/VCF spellings of an absent allele.
    const clean = (t) => (t === '-' || t === '.' ? '' : t.toUpperCase());
    const ref = clean(after[refIdx]);
    const alt = clean(after[altIdx]);
    if (!ref && !alt) return null;

    // Gene: any non-allele alphabetic token, wherever it sits in the row.
    const geneCandidates = toks
        .filter((t, i) => i !== chromIdx && i !== chromIdx + 1)
        .filter((t) => /^[A-Za-z][A-Za-z0-9-]*$/.test(t) && !isChromToken(t));
    const consumed = new Set([after[refIdx], after[altIdx]]);
    const gene = geneCandidates.find((t) => !consumed.has(t) || !isDnaAllele(t)) || null;

    return { chrom, pos, ref, alt, gene: gene ? gene.toUpperCase() : null };
}

// Render a parsed coordinate row as normalised GRCh37 genomic HGVS.
// Trims shared prefix/suffix so VCF-anchored indels ("T" -> "TG") become the
// minimal HGVS event ("g.20148674_20148675insG") the annotation APIs expect.
function coordinateRowToGenomicHgvs(row) {
    if (!row) return null;
    const { chrom, ref, alt } = row;
    const pos = parseInt(row.pos, 10);
    if (!Number.isFinite(pos)) return null;

    if (ref.length === 1 && alt.length === 1) return `chr${chrom}:g.${pos}${ref}>${alt}`;
    // Pure insertion (MAF "-" reference): the new bases sit between pos and pos+1.
    if (!ref) return `chr${chrom}:g.${pos}_${pos + 1}ins${alt}`;
    // Pure deletion (MAF "-" alternate): ref spells out the deleted bases.
    if (!alt) {
        return ref.length === 1
            ? `chr${chrom}:g.${pos}del`
            : `chr${chrom}:g.${pos}_${pos + ref.length - 1}del`;
    }
    // Multi-base ref/alt: trim the shared prefix and suffix to the minimal event.
    // This handles VCF-anchored indels AND same-length pairs carrying context
    // bases ("CT">"CA" is really a T>A substitution) — the annotation APIs index
    // the minimal form, so an untrimmed delins would come back "not found".
    let r = ref;
    let a = alt;
    let start = pos;
    let end = pos + ref.length - 1;
    while (r.length && a.length && r[0] === a[0]) { r = r.slice(1); a = a.slice(1); start += 1; }
    while (r.length && a.length && r[r.length - 1] === a[a.length - 1]) {
        r = r.slice(0, -1); a = a.slice(0, -1); end -= 1;
    }
    // Identical ref and alt describe no variant at all — refuse rather than
    // sending a nonsense delins to the annotation APIs.
    if (!r && !a) return null;
    if (!a) {
        return start === end ? `chr${chrom}:g.${start}del` : `chr${chrom}:g.${start}_${end}del`;
    }
    if (!r) {
        // Insertion between the flanking bases left after trimming.
        return `chr${chrom}:g.${end}_${end + 1}ins${a}`;
    }
    if (r.length === 1 && a.length === 1) return `chr${chrom}:g.${start}${r}>${a}`;
    return start === end
        ? `chr${chrom}:g.${start}delins${a}`
        : `chr${chrom}:g.${start}_${end}delins${a}`;
}

// Build a genomic coordinate tuple from either raw user input (preferred) or a
// normalised g. variant. Returns { chrom, pos, ref, alt, start, end, type } when
// enough information is available, else null.
//
// `ref`/`alt` are null for notations that do not spell the alleles out (plain
// del/dup/inv). Position data is still returned in that case, so coordinate-only
// consumers — the UCSC link, the ClinVar region pull, the nearby-variant plot —
// keep working for indels; allele-dependent consumers resolve the bases against
// the reference sequence via resolveVcfAlleles().
function buildVariantCoordinateTuple(rawInput, gVariant) {
    // Reuse the row parser the query normaliser uses, so every layout that can be
    // looked up also produces working gnomAD/UCSC/SpliceAI coordinates rather
    // than falling back to the less precise g.-string parse.
    const parseTokenInput = (raw) => {
        const row = parseCoordinateRow(raw);
        if (!row) return null;
        let startNum = Number(row.pos);
        if (!Number.isFinite(startNum)) return null;
        let ref = row.ref;
        let alt = row.alt;
        // Same-length pairs may carry shared context bases ("CT">"CA" is really a
        // T>A one base along). VCF SNVs/MNVs need no anchor base, so trim to the
        // minimal representation — the one gnomAD and SpliceAI index. Indels keep
        // their anchored VCF form untouched.
        if (ref.length === alt.length && ref.length > 1) {
            while (ref.length > 1 && ref[0] === alt[0]) {
                ref = ref.slice(1); alt = alt.slice(1); startNum += 1;
            }
            while (ref.length > 1 && ref[ref.length - 1] === alt[alt.length - 1]) {
                ref = ref.slice(0, -1); alt = alt.slice(0, -1);
            }
        }
        return {
            chrom: `chr${row.chrom}`, pos: String(startNum), ref, alt,
            start: startNum, end: startNum + Math.max(ref.length, 1) - 1,
            type: ref.length === 1 && alt.length === 1 ? 'sub' : 'delins'
        };
    };
    // Parse any genomic HGVS form (sub / del / dup / ins / inv / delins).
    const parseHgvs = (gv) => {
        const desc = parseGenomicHgvs(gv);
        if (!desc) return null;
        return {
            chrom: desc.chrom,
            pos: String(desc.start),
            ref: desc.type === 'sub' ? desc.refSeq : null,
            alt: desc.type === 'sub' ? desc.altSeq : null,
            start: desc.start,
            end: desc.end,
            type: desc.type
        };
    };
    return parseTokenInput(rawInput) || parseHgvs(gVariant) || null;
}

// Build a GRCh37 genomic HGVS string from a VEP result object.
//
// VEP resolves gene-level cDNA HGVS ("TSC2:c.4952delA") that the variant recoder
// can miss, and its response already carries the genomic location — but the app
// used to read only the consequence off it, leaving the g. field showing the
// user's own query and every coordinate-derived card dark.
//
// VEP's coordinate conventions: a deletion's start..end span the deleted bases; an
// insertion is reported with start = end + 1 and a "-" reference.
function buildGenomicHgvsFromVep(vepResult) {
    if (!vepResult) return '';
    const chrom = String(vepResult.seq_region_name || '').replace(/^chr/i, '').toUpperCase();
    const start = Number(vepResult.start);
    const end = Number(vepResult.end);
    const alleles = String(vepResult.allele_string || '').split('/');
    if (!/^[0-9XYMT]+$/.test(chrom) || !Number.isFinite(start) || !Number.isFinite(end)) return '';
    const ref = String(alleles[0] || '').toUpperCase();
    const alt = String(alleles[1] || '').toUpperCase();
    if (!ref || !alt) return '';
    if (ref === '-') {
        // Insertion between `end` and `start` (VEP reports start = end + 1).
        if (!/^[ACGTN]+$/.test(alt)) return '';
        return `chr${chrom}:g.${end}_${start}ins${alt}`;
    }
    if (alt === '-') {
        return start === end
            ? `chr${chrom}:g.${start}del`
            : `chr${chrom}:g.${start}_${end}del`;
    }
    if (!/^[ACGTN]+$/.test(ref) || !/^[ACGTN]+$/.test(alt)) return '';
    if (ref.length === 1 && alt.length === 1) return `chr${chrom}:g.${start}${ref}>${alt}`;
    return start === end
        ? `chr${chrom}:g.${start}delins${alt}`
        : `chr${chrom}:g.${start}_${end}delins${alt}`;
}

function reverseComplementSeq(seq) {
    const comp = { A: 'T', T: 'A', C: 'G', G: 'C', N: 'N' };
    return String(seq || '').toUpperCase().split('').reverse().map((b) => comp[b] || 'N').join('');
}

// Fetch GRCh37 reference bases for an inclusive 1-based interval.
async function fetchGrch37ReferenceSequence(chrom, start, end) {
    const c = String(chrom || '').replace(/^chr/i, '');
    if (!c || !Number.isFinite(start) || !Number.isFinite(end) || end < start) return null;
    const url = `https://grch37.rest.ensembl.org/sequence/region/human/${c}:${start}..${end}?content-type=text/plain`;
    const res = await fetchWithTimeout(url, { headers: { 'Accept': 'text/plain' } }, API_TIMEOUT_MS.ensemblSequence);
    if (!res.ok) throw new Error(`Ensembl sequence request failed (${res.status})`);
    const seq = (await res.text()).trim().toUpperCase();
    if (!/^[ACGTN]+$/.test(seq) || seq.length !== end - start + 1) {
        throw new Error('Unexpected Ensembl sequence payload');
    }
    return seq;
}

// Turn an HGVS g. descriptor into anchored VCF-style alleles using a reference
// window. `windowStart` is the 1-based coordinate of seq[0].
//
// HGVS describes indels by the affected span and shifts them 3'; VCF (and therefore
// gnomAD, SpliceAI and every VCF-derived resource) anchors them on the preceding
// base. Without this conversion an indel has no ref/alt at all, which is why the
// allele-dependent cards go blank for anything more complex than a substitution.
function vcfAllelesFromGenomicHgvs(desc, windowStart, seq) {
    if (!desc || !seq) return null;
    const baseAt = (p) => seq[p - windowStart] || '';
    const sliceAt = (a, b) => seq.slice(a - windowStart, b - windowStart + 1);
    const { start, end, type } = desc;
    switch (type) {
        case 'sub': {
            const ref = desc.refSeq || baseAt(start);
            if (!ref || !desc.altSeq) return null;
            return { pos: start, ref, alt: desc.altSeq };
        }
        case 'delins': {
            const ref = sliceAt(start, end) || desc.refSeq;
            if (!ref || !desc.altSeq) return null;
            return { pos: start, ref, alt: desc.altSeq };
        }
        case 'del': {
            const deleted = sliceAt(start, end);
            const anchor = baseAt(start - 1);
            if (!deleted || !anchor) return null;
            return { pos: start - 1, ref: anchor + deleted, alt: anchor };
        }
        case 'dup': {
            // The duplicated copy is inserted immediately after `end`.
            const duplicated = sliceAt(start, end);
            const anchor = baseAt(end);
            if (!duplicated || !anchor) return null;
            return { pos: end, ref: anchor, alt: anchor + duplicated };
        }
        case 'ins': {
            // g.START_ENDins uses START as the base preceding the insertion point.
            const anchor = baseAt(start);
            if (!anchor || !desc.altSeq) return null;
            return { pos: start, ref: anchor, alt: anchor + desc.altSeq };
        }
        case 'inv': {
            const ref = sliceAt(start, end);
            if (!ref) return null;
            return { pos: start, ref, alt: reverseComplementSeq(ref) };
        }
        default:
            return null;
    }
}

// Left-align VCF alleles against the reference (the bcftools/vt normalisation).
// gnomAD stores the 5'-most representation while HGVS uses the 3'-most, so an
// indel in a repeat looks like a different variant unless it is normalised —
// e.g. CFTR F508del is g.117199646_117199648del in HGVS but 117199644 ATCT>A in
// gnomAD. Returns the input unchanged when the reference window runs out.
function leftNormalizeVcfAlleles(pos, ref, alt, windowStart, seq) {
    let p = Number(pos);
    let r = String(ref || '').toUpperCase();
    let a = String(alt || '').toUpperCase();
    if (!Number.isFinite(p) || !r || !a) return { pos, ref, alt };
    // Bounded so a malformed window can never spin here.
    for (let i = 0; i < 500; i++) {
        if (r.length > 0 && a.length > 0 && r[r.length - 1] === a[a.length - 1]) {
            r = r.slice(0, -1);
            a = a.slice(0, -1);
            continue;
        }
        if (r.length === 0 || a.length === 0) {
            const prevBase = seq && seq[p - 1 - windowStart];
            if (!prevBase) break;
            r = prevBase + r;
            a = prevBase + a;
            p -= 1;
            continue;
        }
        break;
    }
    // A fully consumed allele means we walked off the window; keep the input.
    if (!r || !a) return { pos, ref, alt };
    return { pos: p, ref: r, alt: a };
}

// Widest indel we will try to express as VCF alleles. Larger events are structural
// and are not represented as ref/alt strings by gnomAD or SpliceAI anyway.
const MAX_VCF_ALLELE_SPAN = 1000;
// Padding fetched on each side of the variant so left-alignment has room to shift.
const VCF_NORMALIZATION_PAD = 60;

// Resolve GRCh37 VCF-style alleles for the queried variant.
//
// Order of preference:
//   1. alleles already present in the tuple (raw token input, or a substitution)
//   2. MyVariant.info's `vcf` block, when the variant is indexed there
//   3. the reference sequence, for indels that carry no alleles in their notation
//
// Returns { chrom, pos, ref, alt, source } or null. Never throws.
async function resolveVcfAlleles(rawInput, gVariant, annotation) {
    const tuple = buildVariantCoordinateTuple(rawInput, gVariant);
    if (!tuple || !tuple.chrom) return null;
    const chrom = tuple.chrom.replace(/^chr/i, '');
    if (tuple.ref && tuple.alt && tuple.ref !== '-' && tuple.alt !== '-') {
        return { chrom, pos: String(tuple.pos), ref: tuple.ref, alt: tuple.alt, source: 'notation' };
    }
    const vcfData = annotation && annotation.vcf;
    const vcfRef = vcfData && String(vcfData.ref || '').toUpperCase();
    const vcfAlt = vcfData && String(vcfData.alt || '').toUpperCase();
    if (vcfRef && vcfAlt && /^[ACGTN]+$/.test(vcfRef) && /^[ACGTN]+$/.test(vcfAlt)) {
        const vcfPos = vcfData.position ?? vcfData.pos ?? tuple.pos;
        return { chrom, pos: String(vcfPos), ref: vcfRef, alt: vcfAlt, source: 'myvariant' };
    }
    const desc = parseGenomicHgvs(gVariant);
    if (!desc) return null;
    if (desc.end - desc.start > MAX_VCF_ALLELE_SPAN) return null;
    const windowStart = Math.max(1, desc.start - VCF_NORMALIZATION_PAD);
    const windowEnd = desc.end + VCF_NORMALIZATION_PAD;
    try {
        const seq = await fetchGrch37ReferenceSequence(chrom, windowStart, windowEnd);
        if (!seq) return null;
        const anchored = vcfAllelesFromGenomicHgvs(desc, windowStart, seq);
        if (!anchored) return null;
        const normalized = leftNormalizeVcfAlleles(anchored.pos, anchored.ref, anchored.alt, windowStart, seq);
        return { chrom, pos: String(normalized.pos), ref: normalized.ref, alt: normalized.alt, source: 'reference' };
    } catch (err) {
        console.warn('VCF allele resolution failed', err);
        return null;
    }
}

// Build a gnomAD variant-browser URL from either raw token input (preferred) or
// a simple genomic substitution HGVS string.
// vcfData: optional annotation.vcf object {ref, alt, position} from myvariant.info —
// used to recover REF for MNVs where the delins gVariant notation drops the reference.
// resolvedAlleles: optional { chrom, pos, ref, alt } from resolveVcfAlleles(), which
// recovers VCF alleles for indels whose HGVS notation does not spell them out.
function buildGnomadVariantUrl(rawInput, gVariant, annotationId, vcfData, resolvedAlleles) {
    if (resolvedAlleles && resolvedAlleles.chrom && resolvedAlleles.pos && resolvedAlleles.ref && resolvedAlleles.alt) {
        const { chrom, pos, ref, alt } = resolvedAlleles;
        return `https://gnomad.broadinstitute.org/variant/${chrom}-${pos}-${ref}-${alt}?dataset=gnomad_r2_1`;
    }
    // Raw token input (e.g. "10 89692913 PTEN G GG") is preferred because it
    // preserves REF/ALT for indels even when gVariant is normalised to ins/delins.
    const tuple = buildVariantCoordinateTuple(rawInput, gVariant);
    if (tuple) {
        const chrom = tuple.chrom.replace(/^chr/i, '');
        const ref = tuple.ref || (vcfData && String(vcfData.ref || '').toUpperCase()) || null;
        const alt = tuple.alt || (vcfData && String(vcfData.alt || '').toUpperCase()) || null;
        if (ref && alt) {
            return `https://gnomad.broadinstitute.org/variant/${chrom}-${tuple.pos}-${ref}-${alt}?dataset=gnomad_r2_1`;
        }
        // No alleles available — fall back to a region view spanning the variant so
        // the button still appears (and, for an indel, actually covers the event).
        const startNum = Number(tuple.start ?? tuple.pos);
        const endNum = Number(tuple.end ?? tuple.pos);
        return `https://gnomad.broadinstitute.org/region/${chrom}-${startNum}-${Math.max(startNum, endNum)}?dataset=gnomad_r2_1`;
    }
    // Fallback: use annotation/gVariant only for simple substitutions.
    const source = annotationId || gVariant || '';
    const m = String(source).match(/^chr([0-9XYMT]+):g\.(\d+)([A-Za-z-]+)>([A-Za-z-]+)$/i);
    if (!m) return '';
    return `https://gnomad.broadinstitute.org/variant/${m[1]}-${m[2]}-${m[3].toUpperCase()}-${m[4].toUpperCase()}?dataset=gnomad_r2_1`;
}

// Build a UCSC Genome Browser link (hg19/GRCh37) centered on the variant, with a
// minimum 20-nt window (10 upstream, 9 downstream) that widens to keep a
// multi-base event and 10 nt of flank in view.
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
        const lo = Math.min(s, e);
        const hi = Math.max(s, e);
        const center = Math.round((lo + hi) / 2);
        let winStart = Math.max(1, center - 10);
        let winEnd = winStart + 19;
        // Widen for multi-base events so the whole span plus flanks stays visible.
        if (lo - 10 < winStart) winStart = Math.max(1, lo - 10);
        if (hi + 10 > winEnd) winEnd = hi + 10;
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

    const tuple = buildVariantCoordinateTuple(rawInput, gVariant);
    if (tuple && tuple.chrom && tuple.pos) {
        const ucscChrom = toUcscChrom(tuple.chrom);
        const spanStart = Number(tuple.start ?? tuple.pos);
        const spanEnd = Number(tuple.end ?? tuple.pos);
        const win = toTwentyNtWindow(spanStart, spanEnd);
        if (ucscChrom && Number.isFinite(spanStart) && Number.isFinite(spanEnd) && win) {
            const region = `${ucscChrom}:${win.start}-${win.end}`;
            const highlightRegion = `${ucscChrom}:${Math.min(spanStart, spanEnd)}-${Math.max(spanStart, spanEnd)}`;
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

// Convert a single-letter amino acid code to its three-letter form (Title case,
// e.g. 'G' -> 'Gly', '*' -> 'Ter'). Returns null for unknown codes.
function aaSingleToThree(aa) {
    const map = {
        A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly',
        H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser',
        T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter'
    };
    return map[String(aa || '').toUpperCase()] || null;
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

    // In-frame deletions, duplications, insertions and delins (single residue or range):
    //   Lys18del -> K18del, Lys18_Arg20del -> K18_R20del,
    //   Lys18_Arg20delinsGlyPro -> K18_R20delinsGP, Met1? -> unchanged.
    // Replace each three-letter amino-acid token with its single-letter code while leaving
    // the del/dup/ins/delins keywords and positions intact. Only apply when the body actually
    // contains such a keyword and an amino-acid-plus-position token, to avoid touching other forms.
    if (/(?:del|dup|ins)/i.test(body) && /[A-Za-z]{3}\d/.test(body)) {
        const converted = body.replace(
            /Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter|Stop/gi,
            (token) => aaThreeToSingle(token) || token
        );
        if (converted !== body) return converted;
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

// Parse the affected residue range out of a protein HGVS body.
//
// Handles single- and three-letter codes across substitutions, nonsense,
// frameshifts, in-frame del/dup/ins/delins and extensions:
//   p.N1651Mfs*21             -> { start: 1651, end: 1651 }
//   p.Asn1651MetfsTer21       -> { start: 1651, end: 1651 }
//   p.L773_I774delinsF        -> { start: 773,  end: 774  }
//   p.Leu773_Ile774delinsPhe  -> { start: 773,  end: 774  }
//   p.Lys18del                -> { start: 18,   end: 18   }
// Returns null when no residue position can be read.
function parseProteinResidueRange(proteinChange) {
    const body = String(proteinChange || '').trim().replace(/^p\./i, '').replace(/[()]/g, '');
    if (!body) return null;
    // A range is written "<aa><pos>_<aa><pos>"; take the first and last positions.
    const positions = body.match(/(?:[A-Za-z]{3}|[A-Za-z])(\d+)/g);
    if (!positions || positions.length === 0) return null;
    const nums = positions
        .map((tok) => Number((tok.match(/(\d+)$/) || [])[1]))
        .filter((n) => Number.isFinite(n) && n > 0);
    if (nums.length === 0) return null;
    // "fs*21" / "extTer5" trailing counts are lengths, not residue positions, and the
    // regex above cannot match them (they have no preceding amino-acid token), so the
    // remaining numbers are all genuine residue coordinates.
    return { start: Math.min(...nums), end: Math.max(...nums) };
}

// Map a residue range on a gene's canonical protein to GRCh37 genomic coordinates.
//
// Protein notation does not determine a nucleotide change — a frameshift such as
// p.N1651Mfs*21 can arise from many different indels, and Ensembl refuses these
// outright ("Frameshifts are not supported for HGVS protein input"). What *is*
// determined is the codon, so resolve that much and let the position-only cards
// (UCSC, ClinVar region, nearby variants) work instead of failing the whole lookup.
//
// Returns { chrom, start, end, transcript, protein } or null. Never throws.
async function resolveProteinCodonRegion(gene, startAa, endAa) {
    const safeGene = String(gene || '').trim();
    if (!/^[A-Za-z0-9][A-Za-z0-9._-]*$/.test(safeGene)) return null;
    if (!Number.isFinite(startAa) || !Number.isFinite(endAa) || startAa < 1 || endAa < startAa) return null;
    try {
        const lookupUrl = `https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/${encodeURIComponent(safeGene)}?expand=1&content-type=application/json`;
        const lookupRes = await fetchWithTimeout(lookupUrl, { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.ensemblLookup);
        if (!lookupRes.ok) return null;
        const geneData = await lookupRes.json();
        const transcripts = Array.isArray(geneData?.Transcript) ? geneData.Transcript : [];
        const canonicalId = String(geneData?.canonical_transcript || '').split('.')[0];
        const canonical = transcripts.find((t) => t?.id === canonicalId && t?.Translation)
            || transcripts.find((t) => t?.is_canonical && t?.Translation)
            || transcripts.find((t) => t?.Translation);
        const translationId = canonical?.Translation?.id;
        if (!translationId) return null;
        // A residue past the end of the protein means the query does not belong to
        // this transcript; mapping it would silently point at the wrong locus.
        const proteinLength = Number(canonical.Translation.length);
        if (Number.isFinite(proteinLength) && endAa > proteinLength) return null;

        const mapUrl = `https://grch37.rest.ensembl.org/map/translation/${encodeURIComponent(translationId)}/${startAa}..${endAa}?content-type=application/json`;
        const mapRes = await fetchWithTimeout(mapUrl, { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.ensemblLookup);
        if (!mapRes.ok) return null;
        const mapData = await mapRes.json();
        const mappings = Array.isArray(mapData?.mappings) ? mapData.mappings : [];
        if (mappings.length === 0) return null;
        // A residue range spanning an intron maps to several blocks; take the full extent.
        const chrom = String(mappings[0].seq_region_name || '').replace(/^chr/i, '').toUpperCase();
        const starts = mappings.map((m) => Number(m.start)).filter(Number.isFinite);
        const ends = mappings.map((m) => Number(m.end)).filter(Number.isFinite);
        if (!/^[0-9XYMT]+$/.test(chrom) || starts.length === 0 || ends.length === 0) return null;
        return {
            chrom: `chr${chrom}`,
            start: Math.min(...starts),
            end: Math.max(...ends),
            transcript: canonical.id,
            protein: translationId
        };
    } catch (err) {
        console.warn('Protein codon-region resolution failed', err);
        return null;
    }
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
    // Gather candidates from dbNSFP hgvsc/hgvsp arrays. These entries are BARE
    // c./p. strings with no accession ("c.1799T>A"), so an accession-less entry
    // is the cDNA itself — splitting on ':' used to shove it into the transcript
    // slot and leave cDNA empty. hgvsp is index-paired only when the arrays
    // actually align; dbNSFP's hgvsc and hgvsp can have different lengths.
    if (annotation.dbnsfp && annotation.dbnsfp.hgvsc) {
        const hgvsc = Array.isArray(annotation.dbnsfp.hgvsc) ? annotation.dbnsfp.hgvsc : [annotation.dbnsfp.hgvsc];
        const hgvsp = annotation.dbnsfp.hgvsp ? (Array.isArray(annotation.dbnsfp.hgvsp) ? annotation.dbnsfp.hgvsp : [annotation.dbnsfp.hgvsp]) : [];
        const aligned = hgvsp.length === hgvsc.length;
        for (let i = 0; i < hgvsc.length; i++) {
            const sc = hgvsc[i];
            if (!sc) continue;
            const scStr = String(sc);
            const scColon = scStr.indexOf(':');
            const txId = scColon !== -1 ? scStr.slice(0, scColon) : '';
            const cpart = scColon !== -1 ? scStr.slice(scColon + 1) : scStr;
            let ppart = '';
            if (aligned && hgvsp[i]) {
                const spStr = String(hgvsp[i]);
                const spColon = spStr.indexOf(':');
                ppart = spColon !== -1 ? spStr.slice(spColon + 1) : spStr;
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
    // Both ref and alt empty describes no change at all — refuse rather than
    // emit malformed HGVS ("g.123>") that downstream parsers would choke on.
    throw new Error(`SPDI describes no sequence change: ${spdi}`);
}

// The GRCh37 mirror is preferred: complex clinical variants are catalogued
// relative to hg19/GRCh37, and its hgvsg/spdi come back in GRCh37 directly.
// The main host is a degraded-mode fallback for when the mirror is down (it
// 500s in bursts) — its coordinates are GRCh38, which the candidate conversion
// detects from the NC_ accession version and lifts back to hg19.
const RECODER_HOSTS = ['https://grch37.rest.ensembl.org', 'https://rest.ensembl.org'];

async function fetchVariantRecoder(query) {
    const encoded = encodeURIComponent(query);
    const attemptsPerHost = 2;
    let lastErr = null;
    for (let hostIdx = 0; hostIdx < RECODER_HOSTS.length; hostIdx++) {
        const url = `${RECODER_HOSTS[hostIdx]}/variant_recoder/human/${encoded}?content-type=application/json`;
        let response;
        try {
            response = await fetchWithRetry(url, {
                headers: {
                    'Accept': 'application/json'
                }
            }, API_TIMEOUT_MS.recoder, {
                attempts: attemptsPerHost,
                onAttempt: ({ attempt }) => LookupProgress.attempt(
                    'recoder', hostIdx * attemptsPerHost + attempt, RECODER_HOSTS.length * attemptsPerHost)
            });
        } catch (err) {
            // Transport-level failure after this host's retries — try the next host.
            lastErr = err;
            if (err && err.retryable && hostIdx < RECODER_HOSTS.length - 1) continue;
            throw err;
        }
        if (!response.ok) {
            const text = await response.text();
            const err = new Error(`Variant recoder request failed (${response.status}): ${text}`);
            err.status = response.status;
            // A 4xx here means Ensembl parsed the request and rejected the variant —
            // a real answer, no point asking another host. A 5xx means the service
            // is unwell and the variant may be perfectly fine, so fail over.
            err.retryable = RETRYABLE_HTTP_STATUS.has(response.status);
            lastErr = err;
            if (err.retryable && hostIdx < RECODER_HOSTS.length - 1) continue;
            throw err;
        }
        return response.json();
    }
    throw lastErr || new Error('Variant recoder request failed');
}

// Ensembl REST hosts that serve the /map assembly-conversion endpoint. The
// GRCh37 mirror is a separate deployment answering the same paths in both
// directions, and it has stayed up through outages that took the main host's
// /map and /vep endpoints down entirely — so it is worth trying before giving up.
const ENSEMBL_MAP_HOSTS = ['https://rest.ensembl.org', 'https://grch37.rest.ensembl.org'];

// Liftover a genomic variant from hg38 to hg19 using Ensembl REST API. If the
// conversion fails or the variant is already hg19, returns the original
// variant. Expects variant in form 'chr7:g.140753336A>T'.
//
// Returning the input unchanged on failure means callers proceed with an hg38
// coordinate labelled hg19, so an outage on one host quietly degrades lookup
// accuracy rather than erroring — hence the failover across hosts.
async function liftoverHg38ToHg19(variant) {
    const m = variant.match(/^chr([\w]+):g\.(\d+)([A-Za-z-]+)>([A-Za-z-]+)/);
    if (!m) return variant;
    const chrom = m[1];
    const pos = parseInt(m[2], 10);
    const ref = m[3];
    const alt = m[4];
    // Build URL for liftover using Ensembl map endpoint
    const path = `/map/human/GRCh38/${chrom}:${pos}..${pos}:1/GRCh37?content-type=application/json`;
    for (const host of ENSEMBL_MAP_HOSTS) {
        try {
            const res = await fetchWithTimeout(`${host}${path}`, { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.liftover);
            if (!res.ok) continue; // transient or host-specific — try the next host
            const data = await res.json();
            const mappings = (data && data.mappings) || [];
            // Only accept a mapping that landed on the same chromosome in GRCh37;
            // a segment can map to a patch or scaffold, and taking mappings[0]
            // blind turns that into a confident wrong coordinate.
            const mapped = mappings
                .map((entry) => entry && entry.mapped)
                .find((mp) => mp
                    && String(mp.assembly || '').toUpperCase() === 'GRCH37'
                    && String(mp.seq_region_name || '').replace(/^chr/i, '') === String(chrom).replace(/^chr/i, '')
                    && Number.isInteger(mp.start) && mp.start > 0);
            if (mapped) return `chr${chrom}:g.${mapped.start}${ref}>${alt}`;
            // The host answered cleanly with nothing to offer; another host will
            // not know better.
            return variant;
        } catch (err) {
            // Network error or timeout — fall through to the next host.
        }
    }
    return variant;
}



// Read a known GRCh38 start coordinate out of a myvariant.info annotation.
// dbNSFP and dbSNP both carry one for variants they index, which lets the gnomAD
// v4 proxy skip its Ensembl liftover — faster, and it keeps the card working when
// Ensembl's assembly-mapping API is down.
function extractHg38Start(annotation) {
    if (!annotation) return null;
    const candidates = [
        annotation.hg38 && annotation.hg38.start,
        annotation.dbnsfp && annotation.dbnsfp.hg38 && annotation.dbnsfp.hg38.start,
        annotation.dbsnp && annotation.dbsnp.hg38 && annotation.dbsnp.hg38.start,
    ];
    for (const value of candidates) {
        const n = Number.parseInt(value, 10);
        if (Number.isInteger(n) && n > 0) return n;
    }
    return null;
}

// Fetch gnomAD v4.1 (GRCh38) data for a variant.
// Liftover GRCh37→GRCh38 is done server-side by the proxy; pass pos38 when the
// annotation already supplies the GRCh38 start so the proxy can skip it.
// Returns { status, data?, grch38Id?, message? } or null on hard failure.
async function fetchGnomadV4(chrom, pos37, ref, alt, pos38) {
    const c = String(chrom).replace(/^chr/i, '');
    if (!c || !pos37 || !ref || !alt) return null;
    const endpoint = getConfiguredApiEndpoint('GNOMAD_V4_API_ENDPOINT', '/api/gnomad-v4');
    const params = new URLSearchParams({ chrom: c, pos37: String(pos37), ref: ref.toUpperCase(), alt: alt.toUpperCase() });
    const hinted = Number.parseInt(pos38, 10);
    if (Number.isInteger(hinted) && hinted > 0) params.set('pos38', String(hinted));
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

async function fetchClinvarProteinVariants(gene, changeForms) {
    const safeGene = String(gene || '').replace(/[^A-Za-z0-9\-_.]/g, '');
    const change = (Array.isArray(changeForms) ? changeForms : [changeForms])
        .map((s) => String(s || '').replace(/[^A-Za-z0-9*]/g, ''))
        .filter(Boolean)
        .join(',');
    if (!safeGene || !change) return { variants: [], total: 0 };
    const endpoint = getConfiguredApiEndpoint('CLINVAR_PROTEIN_API_ENDPOINT', '/api/clinvar-protein');
    const params = new URLSearchParams({ gene: safeGene, change });
    const res = await fetchWithTimeout(`${endpoint}?${params}`, {}, API_TIMEOUT_MS.clinvar);
    if (!res.ok) throw new Error(`ClinVar protein proxy error: ${res.status}`);
    const data = await res.json();
    if (data?.error) throw new Error(data.error);
    return { variants: data.variants || [], total: data.total ?? (data.variants || []).length };
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

// Return the Watson-Crick complement of a single DNA base, or '' if not A/C/G/T.
function complementBase(b) {
    return { A: 'T', T: 'A', C: 'G', G: 'C' }[String(b || '').toUpperCase()] || '';
}

// Normalise a cDNA HGVS string for comparison across sources. Strips any accession
// prefix and the optionally spelled-out reference bases, so the user's
// "c.2319_2321delAAT" and ClinVar's "c.2319_2321del" compare equal. Inserted bases
// are kept, because insAA and insTT are genuinely different variants.
function normalizeCdnaForMatch(cdna) {
    let s = String(cdna || '').trim().toLowerCase().replace(/\s+/g, '');
    const colon = s.lastIndexOf(':');
    if (colon !== -1) s = s.slice(colon + 1);
    if (!/^[cn]\./.test(s)) return '';
    s = s.replace(/^[cn]\./, 'c.');
    s = s.replace(/del[acgt]+(?=ins)/, 'del');
    s = s.replace(/(del|dup|inv)[acgt]+$/, '$1');
    return s;
}

// Locate the ClinVar region-pull entry that corresponds to the exact queried
// variant (same genomic position and substitution). This is used as a fallback
// when MyVariant.info carries no `clinvar` block for the variant even though
// ClinVar itself has an entry — the per-variant lookup otherwise shows "N/A".
//
// variants: array of region-pull records ({ id, title, variationName, pos, germline }).
// tuple:    queried variant ({ pos, ref, alt }) from buildVariantCoordinateTuple.
//
// Matching is intentionally conservative: only single-nucleotide substitutions
// are recovered, and the substitution bases in the ClinVar title must match the
// queried alleles (allowing for transcript-strand complementation on minus-strand
// genes). When two SNVs share a position (e.g. G>C and G>A) the allele check
// disambiguates; if no entry's alleles match, no fallback is applied.
//
// cdnaForms: optional list of the queried variant's cDNA HGVS strings. Indels have
// no ref/alt to compare, and their genomic coordinates differ between HGVS (3'-shifted)
// and ClinVar's own representation, so they are matched on the c. notation in the
// record title instead — which is exactly how a user would recognise the record.
function findExactClinvarRegionMatch(variants, tuple, cdnaForms) {
    if (!Array.isArray(variants) || !tuple || !tuple.pos) return null;
    const pos = Number(tuple.pos);
    const ref = String(tuple.ref || '').toUpperCase();
    const alt = String(tuple.alt || '').toUpperCase();
    if (!/^[ACGT]$/.test(ref) || !/^[ACGT]$/.test(alt)) {
        // Non-SNV: fall back to cDNA-notation matching when forms were supplied.
        const wanted = new Set(
            (Array.isArray(cdnaForms) ? cdnaForms : [cdnaForms])
                .map(normalizeCdnaForMatch)
                .filter(Boolean)
        );
        if (wanted.size === 0) return null;
        for (const v of variants) {
            const text = [v.title, v.variationName].filter(Boolean).join(' ');
            const tokens = text.match(/c\.[A-Za-z0-9_>+*\-]+/gi) || [];
            if (tokens.some((t) => wanted.has(normalizeCdnaForMatch(t)))) return v;
        }
        return null;
    }

    const refC = complementBase(ref);
    const altC = complementBase(alt);
    const candidates = variants.filter((v) => Number(v.pos) === pos);

    for (const v of candidates) {
        const text = [v.title, v.variationName].filter(Boolean).join(' ');
        const m = text.match(/([ACGT])\s*>\s*([ACGT])/i);
        if (!m) continue;
        const tRef = m[1].toUpperCase();
        const tAlt = m[2].toUpperCase();
        if ((tRef === ref && tAlt === alt) || (tRef === refC && tAlt === altC)) {
            return v;
        }
    }
    return null;
}

function parseProteinPos(text) {
    // Handle both p.Arg273His and p.(Arg273His) (predicted-effect parenthesised form)
    const m = String(text || '').match(/\bp\.\(?(?:[A-Za-z*?]{1,5})?(\d+)/i);
    return m ? Number(m[1]) : null;
}

// Parse a missense protein change from arbitrary text (a ClinVar title, an HGVS
// p. string, or a bare "G34R") into a normalised single-letter key {ref, pos, alt}.
// Accepts three-letter (Gly34Arg / p.(Gly34Arg)) and single-letter (G34R) forms.
// Returns null for anything that is not a simple amino-acid substitution
// (frameshifts, indels, synonymous, etc.) — those are not "same protein change".
function parseProteinChange(text) {
    const s = String(text || '');
    // Three-letter form, optionally wrapped in p.( ) — e.g. "p.(Gly34Arg)".
    let m = s.match(/p\.\(?([A-Za-z]{3})(\d+)([A-Za-z]{3})\)?/);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]);
        if (ref && alt) return { ref, pos: Number(m[2]), alt };
    }
    // Bare three-letter form without the "p." prefix — e.g. "Gly34Arg".
    m = s.match(/\b([A-Za-z]{3})(\d+)([A-Za-z]{3})\b/);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]);
        if (ref && alt) return { ref, pos: Number(m[2]), alt };
    }
    // Single-letter form — e.g. "G34R" or "p.G34R".
    m = s.match(/\bp?\.?([A-Z])(\d+)([A-Z*])\b/);
    if (m) return { ref: m[1].toUpperCase(), pos: Number(m[2]), alt: m[3].toUpperCase() };
    return null;
}

// True when two parsed protein changes describe the same amino-acid substitution.
function sameProteinChange(a, b) {
    return !!a && !!b && a.ref === b.ref && a.pos === b.pos && a.alt === b.alt;
}

// Apply a stable thematic class to each card so CSS can render distinct colors per section.
function applyCardTheme(cardEl, cardTitle) {
    if (!cardEl || !cardTitle) return;
    const key = String(cardTitle).toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
    if (key) cardEl.setAttribute('data-card', key);
}

// Escape a string for safe interpolation into innerHTML.
function escapeHtml(text) {
    return String(text ?? '').replace(/[&<>"']/g, (ch) => (
        { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch]
    ));
}

// Return a URL that is safe to place in an href, or '' when it isn't. Only http(s)
// and mailto schemes are allowed, so untrusted values (e.g. a CIViC source URL or a
// model-supplied link) can't smuggle in `javascript:`/`data:` schemes. The result is
// NOT HTML-escaped — escape it separately when interpolating into an attribute.
function safeUrl(url) {
    const raw = String(url ?? '').trim();
    if (!raw) return '';
    try {
        const parsed = new URL(raw, window.location.href);
        return /^(https?:|mailto:)$/i.test(parsed.protocol) ? raw : '';
    } catch {
        return '';
    }
}

// Render a ClinVar review status as its 0–4 gold-star confidence rating
// (matching the stars on the ClinVar website), or '' when no status is given.
// e.g. "criteria provided, single submitter" → ★☆☆☆.
function clinvarReviewStars(status) {
    const s = String(status || '').toLowerCase();
    if (!s) return '';
    let filled = 0;
    // Check the zero-star "no assertion / no classification" statuses first: the
    // string "no assertion criteria provided" contains "criteria provided" and
    // would otherwise be miscounted as one star.
    if (s.includes('no assertion') || s.includes('no classification')) filled = 0;
    else if (s.includes('practice guideline')) filled = 4;
    else if (s.includes('expert panel')) filled = 3;
    else if (s.includes('multiple submitters')) filled = 2; // "…, no conflicts"
    else if (s.includes('criteria provided')) filled = 1; // single submitter or conflicting
    else filled = 0;
    const stars = '★'.repeat(filled) + '☆'.repeat(4 - filled);
    // display:inline overrides the card-content `span { display:block }` rule so
    // the stars sit inline next to the classification text.
    return `<span style="color:#f59e0b;display:inline" title="${escapeHtml(status)}">${stars}</span>`;
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
 * Query the Biomarker–therapy (BBKB) API for FDA-approved oncology
 * biomarker/therapy combinations associated with the given gene.
 */
async function fetchBbkbBiomarkerTherapies(gene) {
    if (!gene) return { total_matched: 0, returned: 0, results: [] };
    const params = new URLSearchParams({ gene });
    const url = `https://www.drdoubleb.com/BBKB/api/biomarker.php?${params}`;
    let resp;
    try {
        resp = await fetchWithTimeout(url, {}, API_TIMEOUT_MS.bbkb);
    } catch (e) {
        console.warn('BBKB biomarker-therapy fetch failed', e);
        return { total_matched: 0, returned: 0, results: [] };
    }
    if (!resp.ok) throw new Error(`BBKB biomarker-therapy API error: ${resp.status}`);
    const data = await resp.json();
    return {
        query: data.query || { gene },
        matched_genes: data.matched_genes || [],
        total_matched: data.total_matched || 0,
        returned: data.returned || 0,
        offset: data.offset || 0,
        limit: data.limit || 100,
        generated: data.generated || '',
        results: Array.isArray(data.results) ? data.results : []
    };
}

function renderBbkbBiomarkerTherapies(panel, gene, aiExtras) {
    const linkEl = document.createElement('a');
    linkEl.href = 'https://www.drdoubleb.com/BBKB';
    linkEl.target = '_blank';
    linkEl.rel = 'noopener noreferrer';
    linkEl.textContent = 'BBKB biomarker–therapy database ↗';
    panel.appendChild(linkEl);

    const queryLabel = document.createElement('div');
    queryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
    queryLabel.textContent = `Gene: ${gene}`;
    panel.appendChild(queryLabel);

    const resultsDiv = document.createElement('div');
    const spinner = document.createElement('div');
    spinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
    spinner.textContent = 'Loading BBKB therapy results…';
    resultsDiv.appendChild(spinner);
    panel.appendChild(resultsDiv);

    fetchBbkbBiomarkerTherapies(gene).then((data) => {
        if (aiExtras) aiExtras.bbkb_biomarker_therapies = { gene, ...data };
        const records = data.results || [];
        resultsDiv.innerHTML = '';
        if (!records.length) {
            resultsDiv.innerHTML = `<div style="font-size:0.85rem;color:#6b7280;">No BBKB biomarker–therapy records found for ${escapeHtml(gene)}.</div>`;
            return;
        }
        const countEl = document.createElement('div');
        countEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:8px;';
        countEl.textContent = `${records.length} record${records.length !== 1 ? 's' : ''}${data.total_matched && data.total_matched !== records.length ? ` (${data.total_matched} total matched)` : ''}`;
        resultsDiv.appendChild(countEl);

        const BBKB_PREVIEW = 7;
        const buildRecordEl = (rec) => {
            const wrap = document.createElement('div');
            wrap.style.cssText = 'margin-bottom:8px;padding:8px;background:#fff7f7;border:1px solid #fee2e2;border-radius:4px;font-size:0.82rem;';
            const therapies = Array.isArray(rec.therapy) ? rec.therapy.join(', ') : (rec.therapy || 'Unknown therapy');
            const generics = Array.isArray(rec.generic_components) && rec.generic_components.length ? ` (${rec.generic_components.join(', ')})` : '';
            const title = document.createElement('div');
            title.style.cssText = 'font-weight:700;color:#7f1d1d;margin-bottom:2px;';
            title.textContent = `${therapies}${generics}`;
            wrap.appendChild(title);

            const meta = [rec.tumor_type, rec.alteration, rec.therapeutic_context].filter(Boolean);
            if (meta.length) {
                const metaEl = document.createElement('div');
                metaEl.style.cssText = 'color:#374151;margin-bottom:4px;';
                metaEl.textContent = meta.join(' · ');
                wrap.appendChild(metaEl);
            }
            const extras = [rec.setting, rec.approval_status, Array.isArray(rec.fda_application_numbers) ? rec.fda_application_numbers.join(', ') : ''].filter(Boolean);
            if (extras.length) {
                const extraEl = document.createElement('div');
                extraEl.style.cssText = 'color:#6b7280;font-size:0.78rem;';
                extraEl.textContent = extras.join(' · ');
                wrap.appendChild(extraEl);
            }
            return wrap;
        };
        records.slice(0, BBKB_PREVIEW).forEach((rec) => resultsDiv.appendChild(buildRecordEl(rec)));
        if (records.length > BBKB_PREVIEW) {
            const details = document.createElement('details');
            details.style.cssText = 'margin-top:2px;';
            const summary = document.createElement('summary');
            summary.style.cssText = 'font-size:0.82rem;color:#7f1d1d;cursor:pointer;padding:4px 2px;list-style:revert;';
            summary.textContent = `Show ${records.length - BBKB_PREVIEW} more…`;
            details.appendChild(summary);
            records.slice(BBKB_PREVIEW).forEach((rec) => details.appendChild(buildRecordEl(rec)));
            resultsDiv.appendChild(details);
        }
        const note = document.createElement('div');
        note.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
        note.textContent = 'Source: BBKB biomarker–therapy API. Verify against current FDA prescribing information before clinical use.';
        resultsDiv.appendChild(note);
    }).catch(() => {
        resultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">BBKB data unavailable.</div>';
    });
}

/**
 * Query the FDA companion diagnostics API for structured drug/disease/biomarker rows
 * for the given gene. Returns the raw records array from the API response.
 */
async function fetchFdaCompanionDiagnostics(gene) {
    if (!gene) return [];
    const params = new URLSearchParams({ gene });
    const url = `https://www.drdoubleb.com/companion/gene_api.php?${params}`;
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

// ── openFDA gene synonyms & false-positive exceptions ────────────────────
// By default openFDA matching is an "all-caps" filter: the queried gene is
// upper-cased (e.g. "KIT") and a record is kept only when that exact upper-case
// symbol appears as a whole word in indications_and_usage. That is deliberate —
// it stops lowercase English words ("kit", "met") and case-variant non-gene
// uses from matching.
//
// A few genes, though, are written in FDA labels under aliases that are NOT the
// all-caps symbol, so the all-caps filter wrongly drops the label we want. KIT
// is the motivating case: the Gleevec/imatinib label never says "KIT"; it says
// "Kit (CD117)-positive" and "c-Kit mutation". OPENFDA_GENE_SYNONYMS lets us add
// per-gene aliases that are matched case-SENSITIVELY *exactly as written* (so
// "c-Kit" matches "c-Kit", not "C-KIT"). This bypasses the all-caps rule for
// that one symbol while leaving the general filter intact for every other gene.
// (A bare "Kit" alias was tried but produced false positives and was removed —
// "Kit (CD117)" is still matched via the CD117 alias.)
//
// OPENFDA_FALSE_POSITIVE_PHRASES lists case-INSENSITIVE phrases that contain the
// gene's letters but are not the gene — e.g. "BOWEL PREP KIT" (a colonoscopy
// prep, not the KIT receptor). A record is dropped when every gene/synonym hit
// falls inside one of these phrases.
//
// KEEP IN SYNC with tests/openfda-filter.test.js.
const OPENFDA_GENE_SYNONYMS = {
    KIT: ['c-Kit', 'CD117'],
};
const OPENFDA_FALSE_POSITIVE_PHRASES = {
    KIT: ['bowel prep kit'],
};

function openFdaSynonymsFor(gene) {
    return OPENFDA_GENE_SYNONYMS[gene] || [];
}

// All whole-word, case-sensitive match spans of `term` in `text`.
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

// Case-insensitive match spans of each false-positive phrase.
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

// Decide whether an FDA indications_and_usage string matches the queried gene.
// Returns '' to keep, or an exclusion reason: 'case' (neither the all-caps
// symbol nor any verbatim synonym is present), 'boundary' (present only inside a
// larger word), 'falsePositive' (present only inside a known false-positive
// phrase), or 'negation' (every symbol occurrence sits in a negation span).
// For genes with no synonyms and no false-positive phrases this reduces exactly
// to the original three-step filter. KEEP IN SYNC with tests/openfda-filter.test.js.
function openFdaRecordExclusionReason(ind, gene) {
    if (!ind || !gene) return 'case';
    const synonyms = openFdaSynonymsFor(gene);
    const fpSpans = openFdaFalsePositiveSpans(ind, OPENFDA_FALSE_POSITIVE_PHRASES[gene] || []);
    const outsideFp = (s) => !fpSpans.some(f => s.start >= f.start && s.end <= f.end);

    // (1) case-sensitive presence of the all-caps symbol OR a verbatim synonym.
    if (!ind.includes(gene) && !synonyms.some(t => ind.includes(t))) return 'case';

    // Whole-word hits, dropping any that sit wholly inside a false-positive phrase.
    const geneSpans = openFdaWordSpans(ind, gene).filter(outsideFp);
    const synSpans = synonyms.flatMap(t => openFdaWordSpans(ind, t)).filter(outsideFp);

    // (2) no surviving whole-word hit → substring of a larger word, or only
    //     inside a false-positive phrase ("BOWEL PREP KIT").
    if (geneSpans.length === 0 && synSpans.length === 0) {
        const hadRawWord = openFdaGeneAppearsAsWord(ind, gene)
            || synonyms.some(t => openFdaGeneAppearsAsWord(ind, t));
        return (hadRawWord && fpSpans.length) ? 'falsePositive' : 'boundary';
    }

    // (3) negation — drop only when the all-caps symbol is the sole evidence and
    //     every one of its occurrences is negated. A positive synonym hit (e.g.
    //     "c-Kit mutation") rescues the record; FDA negation boilerplate is
    //     written with all-caps symbols, so synonyms aren't negation-checked.
    if (synSpans.length === 0 && openFdaGeneOnlyInNegativeContext(ind, gene)) return 'negation';
    return '';
}

async function fetchOpenFdaDrugLabels(gene) {
    if (!gene) return { total: 0, fetched: 0, excluded: 0, excludedCase: 0, excludedBoundary: 0, excludedFalsePositive: 0, excludedNegation: 0, results: [] };
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
        if (!firstData) return { total: 0, fetched: 0, excluded: 0, excludedCase: 0, excludedBoundary: 0, excludedFalsePositive: 0, excludedNegation: 0, results: [] };

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
        let excludedFalsePositive = 0;
        let excludedNegation = 0;
        const results = [];
        for (const item of raw) {
            const ind = item.indications_and_usage?.[0] || '';
            const reason = openFdaRecordExclusionReason(ind, gene);
            if (reason === 'case') { excludedCase++; continue; }
            if (reason === 'boundary') { excludedBoundary++; continue; }
            if (reason === 'falsePositive') { excludedFalsePositive++; continue; }
            if (reason === 'negation') { excludedNegation++; continue; }
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
            excluded: excludedCase + excludedBoundary + excludedFalsePositive + excludedNegation,
            excludedCase,
            excludedBoundary,
            excludedFalsePositive,
            excludedNegation,
            results,
        };
    } catch (e) {
        console.warn('openFDA fetch failed', e);
        throw e;
    }
}

// openFDA can return many label records for drug-rich genes — HER2/ERBB2 assembles
// ~957 KB, almost entirely openFDA — which risks Vercel's request ceiling and model
// context limits. Cap the number of records for the AI payload, but keep each record's
// full text: `indications_and_usage` is the crucial biomarker-approval section and must
// not be truncated. When records are dropped, add an explicit note so the model (and
// anyone reading the context inspector) knows the list was capped. Returns a NEW object
// so the on-screen openFDA card's data is left untouched.
function condenseOpenFdaForAi(data, { maxRecords = 40 } = {}) {
    if (!data || !Array.isArray(data.results)) return data;
    if (data.results.length <= maxRecords) return data;
    const results = data.results.slice(0, maxRecords);
    const omitted = data.results.length - results.length;
    return {
        ...data,
        results,
        openfda_records_truncated: {
            shown: results.length,
            total: data.results.length,
            omitted,
            note: `openFDA returned ${data.results.length} label records; only the first ${results.length} are included here to bound payload size. Full indications_and_usage text is retained for the included records.`
        }
    };
}

// `variantQuery` ({ gene, variant, extra? }) switches the proxy to its
// LitVar2/PubTator3 path: variant-normalized literature (papers writing the
// same change as p.Val600Glu, c.1799T>A or an rsID all match), relevance-ranked.
// The proxy falls back to the plain `searchTerm` search when no variant entity
// matches, so passing both is always safe.
async function fetchPubmedArticles(searchTerm, limit = 5, variantQuery = null) {
    const hasVariant = variantQuery && String(variantQuery.variant || '').trim();
    if (!searchTerm && !hasVariant) return { total: 0, articles: [] };
    const params = new URLSearchParams({ limit: String(limit) });
    if (searchTerm) params.set('term', searchTerm);
    if (hasVariant) {
        if (variantQuery.gene) params.set('gene', variantQuery.gene);
        params.set('variant', variantQuery.variant);
        if (variantQuery.extra) params.set('extra', variantQuery.extra);
    }
    const endpoint = getConfiguredApiEndpoint('PUBMED_API_ENDPOINT', '/api/pubmed');
    const res = await fetchWithTimeout(appendQueryParams(endpoint, params), {}, API_TIMEOUT_MS.pubmed);
    if (!res.ok) throw new Error(`PubMed proxy error: ${res.status}`);
    const data = await res.json();
    if (data?.error) throw new Error(data.error);
    return { total: data.total || 0, articles: data.articles || [], litvar: data.litvar || null, backend: data.backend || 'term' };
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

/*
 * ── Shared card builders ────────────────────────────────────────────────────
 *
 * The gene-only mode and the full variant mode render the same PubMed, FDA
 * drugs, Clinical Trials and Guidelines cards. These used to be two verbatim
 * copies (~800 lines) that drifted apart one forgotten edit at a time; each
 * card is now built once here and parameterised on the few things that differ
 * (tab labels, search terms, which extras/cache object to populate).
 */

// PubMed card. `tabs` is 1-3 entries of { label, term, extraKey, variantQuery? };
// with a single entry no tab bar is rendered. `variantQuery` ({ gene, variant,
// extra? }) routes that tab through the LitVar2 variant-normalized search, with
// the term search as fallback. Results are stashed under extras[extraKey] for
// the AI-review payload.
function createPubmedCard({ container, tabs, extras }) {
    const PUBMED_LIMIT = 5;
    const pmCard = document.createElement('div');
    pmCard.className = 'card';
    const pmTitle = document.createElement('h3');
    pmTitle.textContent = 'PubMed';
    applyCardTheme(pmCard, 'PubMed');
    pmCard.appendChild(pmTitle);
    const pmContent = document.createElement('div');
    pmContent.className = 'card-content';

    const buildPmResultsPanel = (panel, searchTerm, extraKey, variantQuery) => {
        const queryUrl = searchTerm
            ? `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(searchTerm)}&sort=relevance`
            : 'https://pubmed.ncbi.nlm.nih.gov/';
        const linkEl = document.createElement('a');
        linkEl.href = queryUrl;
        linkEl.target = '_blank';
        linkEl.rel = 'noopener noreferrer';
        linkEl.textContent = 'Search PubMed ↗';
        panel.appendChild(linkEl);
        if (searchTerm) {
            const queryLabel = document.createElement('div');
            queryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
            queryLabel.textContent = `Query: "${searchTerm}"`;
            panel.appendChild(queryLabel);
        }
        const resultsDiv = document.createElement('div');
        if (searchTerm || variantQuery) {
            const spinner = document.createElement('div');
            spinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
            spinner.textContent = 'Loading PubMed results…';
            resultsDiv.appendChild(spinner);
        }
        panel.appendChild(resultsDiv);
        if (searchTerm || variantQuery) {
            fetchPubmedArticles(searchTerm, PUBMED_LIMIT, variantQuery).then(({ total, articles, litvar, backend }) => {
                extras[extraKey] = { query: searchTerm, total, articles, backend, litvar };
                resultsDiv.innerHTML = '';
                // Provenance line: which spelling-normalized entity matched. A
                // LitVar2 hit means every nomenclature form of the variant was
                // searched, not just the literal query text above.
                if (litvar && backend === 'litvar') {
                    const lvNote = document.createElement('div');
                    lvNote.style.cssText = 'font-size:0.78rem;color:#0369a1;margin-bottom:6px;';
                    lvNote.textContent = `Variant-matched via LitVar2: ${litvar.name}${litvar.rsid ? ` (${litvar.rsid})` : ''} — all nomenclature spellings, ranked by relevance`;
                    resultsDiv.appendChild(lvNote);
                }
                if (total === 0 || articles.length === 0) {
                    const noRes = document.createElement('div');
                    noRes.style.cssText = 'font-size:0.85rem;color:#6b7280;';
                    noRes.textContent = 'No PubMed results found.';
                    resultsDiv.appendChild(noRes);
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
                    seeAll.href = (litvar && litvar.url) ? litvar.url : queryUrl;
                    seeAll.target = '_blank';
                    seeAll.rel = 'noopener noreferrer';
                    seeAll.style.fontSize = '0.82rem';
                    seeAll.textContent = (litvar && litvar.url)
                        ? `See all ${total.toLocaleString()} variant-matched articles on LitVar2 ↗`
                        : `See all ${total.toLocaleString()} results on PubMed ↗`;
                    resultsDiv.appendChild(seeAll);
                }
            }).catch(() => {
                resultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">PubMed unavailable.</div>';
            });
        }
    };

    if (tabs.length > 1) {
        const tabBar = document.createElement('div');
        tabBar.className = 'card-tabs';
        const tabBtns = [];
        const tabPanels = [];
        const activateTab = (idx) => {
            tabBtns.forEach((b, i) => b.classList.toggle('active', i === idx));
            tabPanels.forEach((p, i) => p.classList.toggle('active', i === idx));
        };
        tabs.forEach((tab, i) => {
            const btn = document.createElement('button');
            btn.type = 'button';
            btn.className = i === 0 ? 'card-tab-btn active' : 'card-tab-btn';
            btn.textContent = tab.label;
            btn.addEventListener('click', () => activateTab(i));
            tabBar.appendChild(btn);
            tabBtns.push(btn);
            const panel = document.createElement('div');
            panel.className = i === 0 ? 'card-tab-panel active' : 'card-tab-panel';
            tabPanels.push(panel);
        });
        pmContent.appendChild(tabBar);
        tabs.forEach((tab, i) => {
            buildPmResultsPanel(tabPanels[i], tab.term, tab.extraKey, tab.variantQuery || null);
            pmContent.appendChild(tabPanels[i]);
        });
    } else {
        buildPmResultsPanel(pmContent, tabs[0].term, tabs[0].extraKey, tabs[0].variantQuery || null);
    }

    pmCard.appendChild(pmContent);
    if (container) container.appendChild(pmCard);
    return pmCard;
}

// FDA-Approved Drugs card: Companion Dx / openFDA Labels / BBKB tabs.
// Companion-Dx records land in cardCache.fda_companion_diagnostics_records and
// the openFDA payload in extras.openfda so the AI-review buildContext can reuse
// them instead of re-fetching.
function createFdaDrugsCard({ container, gene, extras, cardCache }) {
    const fdaCard = document.createElement('div');
    fdaCard.className = 'card';
    const fdaTitle = document.createElement('h3');
    fdaTitle.textContent = 'FDA-Approved Drugs (by gene)';
    applyCardTheme(fdaCard, 'FDA-Approved Drugs (by gene)');
    fdaCard.appendChild(fdaTitle);
    const fdaContent = document.createElement('div');
    fdaContent.className = 'card-content';

    if (!gene) {
        const noGeneEl = document.createElement('div');
        noGeneEl.style.cssText = 'font-size:0.85rem;color:#6b7280;';
        noGeneEl.textContent = 'No gene identified for FDA drug lookup.';
        fdaContent.appendChild(noGeneEl);
        fdaCard.appendChild(fdaContent);
        if (container) container.appendChild(fdaCard);
        return fdaCard;
    }

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
    const bbkbBtn = document.createElement('button');
    bbkbBtn.type = 'button';
    bbkbBtn.className = 'card-tab-btn';
    bbkbBtn.textContent = 'BBKB';
    fdaTabBar.appendChild(compDxBtn);
    fdaTabBar.appendChild(openFdaBtn);
    fdaTabBar.appendChild(bbkbBtn);
    fdaContent.appendChild(fdaTabBar);

    const compDxPanel = document.createElement('div');
    compDxPanel.className = 'card-tab-panel active';
    const openFdaPanel = document.createElement('div');
    openFdaPanel.className = 'card-tab-panel';
    const bbkbPanel = document.createElement('div');
    bbkbPanel.className = 'card-tab-panel';

    compDxBtn.addEventListener('click', () => {
        compDxBtn.classList.add('active'); openFdaBtn.classList.remove('active'); bbkbBtn.classList.remove('active');
        compDxPanel.classList.add('active'); openFdaPanel.classList.remove('active'); bbkbPanel.classList.remove('active');
    });
    openFdaBtn.addEventListener('click', () => {
        openFdaBtn.classList.add('active'); compDxBtn.classList.remove('active'); bbkbBtn.classList.remove('active');
        openFdaPanel.classList.add('active'); compDxPanel.classList.remove('active'); bbkbPanel.classList.remove('active');
    });
    bbkbBtn.addEventListener('click', () => {
        bbkbBtn.classList.add('active'); compDxBtn.classList.remove('active'); openFdaBtn.classList.remove('active');
        bbkbPanel.classList.add('active'); compDxPanel.classList.remove('active'); openFdaPanel.classList.remove('active');
    });

    fdaContent.appendChild(compDxPanel);
    fdaContent.appendChild(openFdaPanel);
    fdaContent.appendChild(bbkbPanel);
    fdaCard.appendChild(fdaContent);
    if (container) container.appendChild(fdaCard);

    // --- Companion Dx panel ---
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
        if (cardCache) cardCache.fda_companion_diagnostics_records = records;
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

    // --- BBKB panel ---
    renderBbkbBiomarkerTherapies(bbkbPanel, gene, extras);

    // --- openFDA Labels panel ---
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

    fetchOpenFdaDrugLabels(gene).then(({ total, fetched, excluded, excludedCase, excludedBoundary, excludedFalsePositive, excludedNegation, results: ofResults }) => {
        extras.openfda = { gene, total, fetched, excluded, excludedCase, excludedBoundary, excludedFalsePositive, excludedNegation, results: ofResults };
        ofResultsDiv.innerHTML = '';
        if (!ofResults || ofResults.length === 0) {
            ofResultsDiv.innerHTML = `<div style="font-size:0.85rem;color:#6b7280;">No openFDA drug label results found for ${escapeHtml(gene)}.</div>`;
            return;
        }
        const OF_PREVIEW = 7;
        const ofCountEl = document.createElement('div');
        ofCountEl.style.cssText = 'font-size:0.85rem;font-weight:600;margin-bottom:8px;';
        const shownOf = `${Math.min(ofResults.length, OF_PREVIEW)} of ${ofResults.length} result${ofResults.length !== 1 ? 's' : ''}`;
        const excludedReasons = [];
        if (excludedCase) excludedReasons.push(`${excludedCase} case mismatch`);
        if (excludedBoundary) excludedReasons.push(`${excludedBoundary} substring of larger word`);
        if (excludedFalsePositive) excludedReasons.push(`${excludedFalsePositive} false-positive phrase`);
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

    return fdaCard;
}

// Clinical Trials card. When `gene` is empty only the search link renders
// (matching the variant mode's no-gene behaviour). Fetched studies are stashed
// in cardCache.clinical_trials for the AI-review buildContext.
function createClinicalTrialsCard({ container, gene, tumorType, cardCache }) {
    const CT_PREVIEW = 5;
    const ctCard = document.createElement('div');
    ctCard.className = 'card';
    const ctTitle = document.createElement('h3');
    ctTitle.textContent = 'Clinical Trials';
    applyCardTheme(ctCard, 'Clinical Trials');
    ctCard.appendChild(ctTitle);
    const ctContent = document.createElement('div');
    ctContent.className = 'card-content';

    const ctBaseUrl = gene
        ? `https://clinicaltrials.gov/search?query=${encodeURIComponent([gene, tumorType].filter(Boolean).join(' '))}&recrs=b&type=Intr`
        : 'https://clinicaltrials.gov/';

    const ctLinkEl = document.createElement('a');
    ctLinkEl.href = ctBaseUrl;
    ctLinkEl.target = '_blank';
    ctLinkEl.rel = 'noopener noreferrer';
    ctLinkEl.textContent = 'Search ClinicalTrials.gov ↗';
    ctContent.appendChild(ctLinkEl);

    if (gene) {
        const ctQueryLabel = document.createElement('div');
        ctQueryLabel.style.cssText = 'font-size:0.8rem;color:#6b7280;margin:2px 0 6px;';
        const queryParts = [gene, tumorType].filter(Boolean);
        ctQueryLabel.textContent = `Query: "${queryParts.join(' + ')}" · Interventional · Recruiting · Phase 2+ · US`;
        ctContent.appendChild(ctQueryLabel);
    }

    const ctResultsDiv = document.createElement('div');
    if (gene) {
        const ctSpinner = document.createElement('div');
        ctSpinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
        ctSpinner.textContent = 'Loading clinical trials…';
        ctResultsDiv.appendChild(ctSpinner);
    }
    ctContent.appendChild(ctResultsDiv);
    ctCard.appendChild(ctContent);
    if (container) container.appendChild(ctCard);

    if (gene) {
        fetchClinicalTrials(gene, tumorType).then(({ total, studies }) => {
            if (cardCache) cardCache.clinical_trials = { total, studies };
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

    return ctCard;
}

/*
 * ── Cancer Prevalence card (cBioPortal) ─────────────────────────────────────
 *
 * "How often is this variant seen, and in which tumor types?" answered from
 * cBioPortal's public API (CORS-open — called directly from the browser, no
 * serverless function spent). v1 queries one large pan-cancer cohort,
 * MSK-IMPACT 2017 (10,945 prospectively sequenced tumors), which is frozen —
 * hence the hard-coded denominator — and clearly labels the source.
 */
const CBIOPORTAL_API = 'https://www.cbioportal.org/api';
const CBIOPORTAL_STUDY = {
    id: 'msk_impact_2017',
    name: 'MSK-IMPACT 2017',
    sampleListId: 'msk_impact_2017_all',
    profileId: 'msk_impact_2017_mutations',
    totalSamples: 10945, // allSampleCount of the frozen study
    url: 'https://www.cbioportal.org/study/summary?id=msk_impact_2017'
};

async function fetchCbioportalGeneMutations(gene) {
    const geneRes = await fetchWithTimeout(
        `${CBIOPORTAL_API}/genes/${encodeURIComponent(gene)}`,
        { headers: { 'Accept': 'application/json' } }, API_TIMEOUT_MS.cbioportal);
    if (!geneRes.ok) return null; // unknown symbol → no card data
    const geneData = await geneRes.json();
    const entrez = Number(geneData && geneData.entrezGeneId);
    if (!Number.isInteger(entrez) || entrez <= 0) return null;
    const mutRes = await fetchWithTimeout(
        `${CBIOPORTAL_API}/molecular-profiles/${CBIOPORTAL_STUDY.profileId}/mutations/fetch?projection=SUMMARY`,
        {
            method: 'POST',
            headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
            body: JSON.stringify({ sampleListId: CBIOPORTAL_STUDY.sampleListId, entrezGeneIds: [entrez] })
        }, API_TIMEOUT_MS.cbioportal);
    if (!mutRes.ok) throw new Error(`cBioPortal mutations fetch failed (${mutRes.status})`);
    const mutations = await mutRes.json();
    return { entrez, mutations: Array.isArray(mutations) ? mutations : [] };
}

// Resolve CANCER_TYPE for a set of sample IDs. Capped so a very common gene
// (TP53) cannot post a giant identifier list; callers note when capped.
async function fetchCbioportalCancerTypes(sampleIds, cap = 1000) {
    const ids = sampleIds.slice(0, cap);
    if (ids.length === 0) return {};
    const res = await fetchWithTimeout(
        `${CBIOPORTAL_API}/clinical-data/fetch?clinicalDataType=SAMPLE&projection=SUMMARY`,
        {
            method: 'POST',
            headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
            body: JSON.stringify({
                attributeIds: ['CANCER_TYPE'],
                identifiers: ids.map((id) => ({ entityId: id, studyId: CBIOPORTAL_STUDY.id }))
            })
        }, API_TIMEOUT_MS.cbioportal);
    if (!res.ok) throw new Error(`cBioPortal clinical-data fetch failed (${res.status})`);
    const rows = await res.json();
    const bySample = {};
    (Array.isArray(rows) ? rows : []).forEach((r) => {
        if (r && r.sampleId && r.value) bySample[r.sampleId] = r.value;
    });
    return bySample;
}

// Cancer Prevalence card. `proteinChange` (optional) is matched against the
// cohort's protein changes after single-letter normalization; without it the
// card shows gene-level prevalence only. Results land in
// extras.cbioportal_prevalence for the AI-review payload.
function createCbioportalCard({ container, gene, proteinChange, extras }) {
    if (!gene) return null;
    const card = document.createElement('div');
    card.className = 'card';
    const title = document.createElement('h3');
    title.textContent = 'Cancer Prevalence';
    applyCardTheme(card, 'Cancer Prevalence');
    card.appendChild(title);
    const content = document.createElement('div');
    content.className = 'card-content';

    const sourceLine = document.createElement('div');
    sourceLine.style.cssText = 'font-size:0.78rem;color:#6b7280;margin-bottom:4px;';
    sourceLine.textContent = `${CBIOPORTAL_STUDY.name} · ${CBIOPORTAL_STUDY.totalSamples.toLocaleString()} sequenced tumors (cBioPortal)`;
    content.appendChild(sourceLine);

    const linkEl = document.createElement('a');
    linkEl.href = `https://www.cbioportal.org/results?cancer_study_list=${encodeURIComponent(CBIOPORTAL_STUDY.id)}&gene_list=${encodeURIComponent(gene)}`;
    linkEl.target = '_blank';
    linkEl.rel = 'noopener noreferrer';
    linkEl.textContent = `View ${gene} on cBioPortal ↗`;
    content.appendChild(linkEl);

    const resultsDiv = document.createElement('div');
    resultsDiv.style.cssText = 'margin-top:0.4rem;';
    const spinner = document.createElement('div');
    spinner.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
    spinner.textContent = 'Loading cohort prevalence…';
    resultsDiv.appendChild(spinner);
    content.appendChild(resultsDiv);
    card.appendChild(content);
    if (container) container.appendChild(card);

    const pct = (n, d) => (d > 0 ? `${((n / d) * 100).toFixed(n / d < 0.01 ? 2 : 1)}%` : 'n/a');
    const line = (html) => {
        const div = document.createElement('div');
        div.style.cssText = 'font-size:0.86rem;margin-bottom:2px;';
        div.innerHTML = html;
        return div;
    };

    (async () => {
        const data = await fetchCbioportalGeneMutations(gene);
        if (!data) {
            resultsDiv.innerHTML = `<div style="font-size:0.82rem;color:#9ca3af;">Gene "${escapeHtml(gene)}" not found in cBioPortal.</div>`;
            return;
        }
        const { mutations } = data;
        const total = CBIOPORTAL_STUDY.totalSamples;
        const geneSamples = new Set(mutations.map((m) => m.sampleId).filter(Boolean));

        const targetNorm = normaliseProteinForCivicMatch(proteinChange || '');
        const matchedSampleSet = new Set();
        const changeCounts = {};
        mutations.forEach((m) => {
            const change = String(m.proteinChange || '').trim();
            if (change) changeCounts[change] = (changeCounts[change] || 0) + 1;
            if (targetNorm && normaliseProteinForCivicMatch(change) === targetNorm && m.sampleId) {
                matchedSampleSet.add(m.sampleId);
            }
        });

        resultsDiv.innerHTML = '';
        resultsDiv.appendChild(line(`<strong>${escapeHtml(gene)} mutated:</strong> ${geneSamples.size.toLocaleString()} of ${total.toLocaleString()} tumors (${pct(geneSamples.size, total)})`));
        if (targetNorm) {
            if (matchedSampleSet.size > 0) {
                resultsDiv.appendChild(line(`<strong>${escapeHtml(targetNorm)} specifically:</strong> ${matchedSampleSet.size.toLocaleString()} tumors (${pct(matchedSampleSet.size, total)} of cohort, ${pct(matchedSampleSet.size, geneSamples.size)} of ${escapeHtml(gene)}-mutant)`));
            } else {
                resultsDiv.appendChild(line(`<strong>${escapeHtml(targetNorm)}:</strong> not seen in this cohort`));
            }
        }

        // Top protein changes in the gene, for context.
        const topChanges = Object.entries(changeCounts).sort((a, b) => b[1] - a[1]).slice(0, 5);
        if (topChanges.length > 0) {
            const topDet = document.createElement('details');
            const topSum = document.createElement('summary');
            topSum.style.cssText = 'font-size:0.82rem;cursor:pointer;';
            topSum.textContent = `Most frequent ${gene} changes in this cohort`;
            topDet.appendChild(topSum);
            const ul = document.createElement('ul');
            ul.style.cssText = 'margin-top:0.3rem;font-size:0.82rem;padding-left:1.2rem;line-height:1.5;';
            topChanges.forEach(([change, n]) => {
                const li = document.createElement('li');
                li.textContent = `${change} — ${n.toLocaleString()}`;
                ul.appendChild(li);
            });
            topDet.appendChild(ul);
            resultsDiv.appendChild(topDet);
        }

        // Tumor-type breakdown — for the exact variant when matched, else the gene.
        const breakdownIds = matchedSampleSet.size > 0 ? [...matchedSampleSet] : [...geneSamples];
        const breakdownLabel = matchedSampleSet.size > 0 ? targetNorm : `${gene}-mutant`;
        let breakdown = [];
        let breakdownCapped = false;
        try {
            const CAP = 1000;
            breakdownCapped = breakdownIds.length > CAP;
            const bySample = await fetchCbioportalCancerTypes(breakdownIds, CAP);
            // Count only the samples we asked about, whatever the response holds.
            const wanted = new Set(breakdownIds.slice(0, CAP));
            const counts = {};
            Object.entries(bySample).forEach(([sid, ct]) => {
                if (wanted.has(sid)) counts[ct] = (counts[ct] || 0) + 1;
            });
            breakdown = Object.entries(counts).sort((a, b) => b[1] - a[1]);
        } catch (e) {
            console.warn('cBioPortal cancer-type breakdown failed', e);
        }
        if (breakdown.length > 0) {
            const bkTitle = document.createElement('div');
            bkTitle.style.cssText = 'font-size:0.84rem;font-weight:600;margin-top:0.4rem;';
            bkTitle.textContent = `Tumor types with ${breakdownLabel}${breakdownCapped ? ' (first 1,000 samples)' : ''}:`;
            resultsDiv.appendChild(bkTitle);
            const shown = breakdown.slice(0, 8);
            const bkDiv = document.createElement('div');
            bkDiv.style.cssText = 'font-size:0.82rem;color:#374151;line-height:1.5;';
            bkDiv.textContent = shown.map(([ct, n]) => `${ct} (${n.toLocaleString()})`).join(' · ')
                + (breakdown.length > shown.length ? ` · +${breakdown.length - shown.length} more` : '');
            resultsDiv.appendChild(bkDiv);
        }

        const note = document.createElement('div');
        note.style.cssText = 'font-size:0.75rem;color:#9ca3af;margin-top:6px;';
        note.textContent = 'Source: cBioPortal (MSK-IMPACT 2017, targeted panel of a single institution\'s advanced-cancer patients — not population prevalence). Research use only.';
        resultsDiv.appendChild(note);

        extras.cbioportal_prevalence = {
            study: { id: CBIOPORTAL_STUDY.id, name: CBIOPORTAL_STUDY.name, total_samples: total },
            gene,
            gene_mutant_samples: geneSamples.size,
            variant: targetNorm || null,
            variant_samples: targetNorm ? matchedSampleSet.size : null,
            top_protein_changes: topChanges.map(([change, n]) => ({ change, count: n })),
            cancer_type_breakdown: breakdown.map(([cancer_type, n]) => ({ cancer_type, count: n })),
            cancer_type_breakdown_for: breakdownLabel,
            cancer_type_breakdown_capped: breakdownCapped,
            note: 'Single-institution advanced-cancer sequencing cohort; frequencies are cohort prevalence, not population or incidence figures.'
        };
    })().catch((err) => {
        console.warn('cBioPortal prevalence card failed', err);
        resultsDiv.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">cBioPortal unavailable.</div>';
    });

    return card;
}

// Guidelines card. `getGene` is resolved at selection time (the variant mode's
// gene may arrive after render); results land in extras.guidelines.
function createGuidelinesCard({ container, tumorType, getGene, extras }) {
    const GUIDELINES_BASE = 'https://www.drdoubleb.com/guidelines';
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
    if (container) container.appendChild(guidelinesCard);

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
        const gene = getGene();
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

            extras.guidelines = {
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

            const renderKV = (obj, kvContainer, depth) => {
                Object.entries(obj).forEach(([key, val]) => {
                    if (val === null || val === undefined || val === '' || (Array.isArray(val) && val.length === 0) || (typeof val === 'object' && !Array.isArray(val) && Object.keys(val).length === 0)) return;
                    const fmtKey = key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
                    if (typeof val === 'object' && !Array.isArray(val)) {
                        const subTitle = document.createElement('div');
                        subTitle.style.cssText = `padding-left:${depth * 12}px;font-weight:600;color:#374151;margin-top:3px;`;
                        subTitle.textContent = fmtKey + ':';
                        kvContainer.appendChild(subTitle);
                        renderKV(val, kvContainer, depth + 1);
                    } else if (Array.isArray(val) && val.some(v => v !== null && typeof v === 'object')) {
                        // Array of objects — render each item as a nested block
                        const subTitle = document.createElement('div');
                        subTitle.style.cssText = `padding-left:${depth * 12}px;font-weight:600;color:#374151;margin-top:3px;`;
                        subTitle.textContent = fmtKey + ':';
                        kvContainer.appendChild(subTitle);
                        val.forEach(item => {
                            if (item !== null && typeof item === 'object') {
                                renderKV(item, kvContainer, depth + 1);
                            } else if (item !== null && item !== '') {
                                const row = document.createElement('div');
                                row.style.cssText = `padding-left:${(depth + 1) * 12}px;margin:2px 0;color:#1f2937;`;
                                row.textContent = String(item);
                                kvContainer.appendChild(row);
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
                        kvContainer.appendChild(row);
                    }
                });
            };

            const renderSection = (title, obj, sectionContainer) => {
                if (!obj || typeof obj !== 'object' || Object.keys(obj).length === 0) return;
                const sectionTitle = document.createElement('div');
                sectionTitle.style.cssText = 'font-weight:700;color:#166534;margin:6px 0 2px;font-size:0.79rem;text-transform:uppercase;letter-spacing:0.04em;border-top:1px solid #bbf7d0;padding-top:5px;';
                sectionTitle.textContent = title;
                sectionContainer.appendChild(sectionTitle);
                renderKV(obj, sectionContainer, 0);
            };

            data.results.forEach((record, idx) => {
                const recordEl = document.createElement('div');
                recordEl.style.cssText = 'margin-bottom:8px;padding:8px;background:#f0fdf4;border:1px solid #bbf7d0;border-radius:4px;font-size:0.82rem;';
                const headerEl = document.createElement('div');
                headerEl.style.cssText = 'font-weight:700;color:#166534;margin-bottom:4px;';
                const markerDisplay = (record.biomarker && (record.biomarker.display_name || record.biomarker.symbol)) || `Record ${idx + 1}`;
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

                if (Array.isArray(record.therapy) && record.therapy.length > 0) {
                    buildKVRow('Recommended therapy', record.therapy.map(t => t.name).filter(Boolean).join(', '));
                    const cat = record.therapy.map(t => t.evidence_or_category).filter(Boolean).join(', ');
                    if (cat) buildKVRow('Category', cat);
                    const setting = record.therapy.map(t => t.setting).filter(Boolean).join(', ');
                    if (setting) buildKVRow('Therapy setting', setting);
                }
                if (record.testing) {
                    const te = record.testing;
                    if (te.recommended) buildKVRow('Recommended test', te.recommended);
                    if (Array.isArray(te.methods) && te.methods.length > 0) {
                        buildKVRow('Method', te.methods.map(m => m.method).filter(Boolean).join(', '));
                    }
                }
                if (record.clinical_significance && record.clinical_significance.summary) {
                    buildKVRow('Clinical significance', record.clinical_significance.summary);
                }

                const detailsEl = document.createElement('details');
                detailsEl.style.cssText = 'margin-top:6px;';
                const detailsSummary = document.createElement('summary');
                detailsSummary.style.cssText = 'font-size:0.80rem;color:#15803d;cursor:pointer;padding:2px 0;list-style:revert;';
                detailsSummary.textContent = 'Full record details';
                detailsEl.appendChild(detailsSummary);
                const detailsBody = document.createElement('div');
                detailsBody.style.cssText = 'font-size:0.79rem;padding:4px 2px;margin-top:4px;max-height:350px;overflow-y:auto;line-height:1.6;';
                if (record.id) {
                    const idRow = document.createElement('div');
                    idRow.style.cssText = 'font-size:0.74rem;color:#9ca3af;margin-bottom:2px;';
                    idRow.textContent = `Record ID: ${record.id}`;
                    detailsBody.appendChild(idRow);
                }
                const sectionDefs = [['Tumor', record.tumor], ['Biomarker', record.biomarker], ['Clinical Significance', record.clinical_significance], ['Testing', record.testing], ['Therapy', record.therapy], ['Eligibility Context', record.eligibility_context], ['Interpretive States', record.interpretive_states], ['Guideline Source', record.guideline_metadata]];
                sectionDefs.forEach(([title, obj]) => renderSection(title, obj, detailsBody));
                const skipKeys = new Set(['id', 'tumor', 'biomarker', 'clinical_significance', 'testing', 'therapy', 'eligibility_context', 'interpretive_states', 'guideline_metadata', 'dataset_file', 'dataset_record_index', 'dataset_name']);
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
            guidelinesResults.innerHTML = `<div style="font-size:0.85rem;color:#9ca3af;">Guidelines data unavailable: ${escapeHtml(err.message)}</div>`;
            extras.guidelines = { error: err.message, cancer_type: selectedCancer, gene };
        }
    });

    return guidelinesCard;
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

async function fetchAiReview(context, model, extra = {}) {
    const endpoint = getConfiguredApiEndpoint('AI_REVIEW_API_ENDPOINT', '/api/ai-review');
    const payload = { model, context };
    // Optional bring-your-own OpenRouter key and Cloudflare Turnstile token. Both are
    // omitted when empty so the backend uses the owner key / skips verification.
    if (extra.apiKey) payload.openrouter_api_key = extra.apiKey;
    if (extra.turnstileToken) payload.turnstile_token = extra.turnstileToken;
    const res = await fetchWithTimeout(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
    }, API_TIMEOUT_MS.aiReview);
    const data = await res.json().catch(() => ({}));
    if (!res.ok || data?.error) {
        const detail = data?.detail ? `: ${data.detail}` : '';
        throw new Error(`${data?.error || `AI review request failed (${res.status})`}${detail}`);
    }
    return data;
}

// Fetch the fully assembled prompt from the backend without calling any model.
// Powers the "Copy prompt" button so users can paste it into their own LLM.
async function fetchAiReviewPrompt(context, model) {
    const endpoint = getConfiguredApiEndpoint('AI_REVIEW_API_ENDPOINT', '/api/ai-review');
    const res = await fetchWithTimeout(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ model, context, mode: 'prompt' })
    }, API_TIMEOUT_MS.aiReview);
    const data = await res.json().catch(() => ({}));
    if (!res.ok || data?.error) {
        throw new Error(data?.error || `Prompt request failed (${res.status})`);
    }
    return data.prompt || '';
}

// sessionStorage helpers for the optional BYO OpenRouter key. Stored only for the
// current tab (cleared when it closes) so returning users don't re-paste each lookup.
const AI_REVIEW_KEY_STORAGE = 'aiReviewOpenRouterKey';
function getStoredUserKey() {
    try { return sessionStorage.getItem(AI_REVIEW_KEY_STORAGE) || ''; } catch { return ''; }
}
function setStoredUserKey(value) {
    try {
        if (value) sessionStorage.setItem(AI_REVIEW_KEY_STORAGE, value);
        else sessionStorage.removeItem(AI_REVIEW_KEY_STORAGE);
    } catch { /* storage unavailable — ignore */ }
}

// Lazily load the Cloudflare Turnstile script, only when a site key is configured.
let _turnstileScriptPromise = null;
function loadTurnstileScript() {
    if (window.turnstile) return Promise.resolve();
    if (_turnstileScriptPromise) return _turnstileScriptPromise;
    _turnstileScriptPromise = new Promise((resolve, reject) => {
        const s = document.createElement('script');
        s.src = 'https://challenges.cloudflare.com/turnstile/v0/api.js?render=explicit';
        s.async = true;
        s.defer = true;
        s.onload = () => resolve();
        s.onerror = () => reject(new Error('Failed to load Turnstile'));
        document.head.appendChild(s);
    });
    return _turnstileScriptPromise;
}

// Build the shared "Optional AI Review" card. Both the gene-only and full-variant
// flows use this so the model picker, notes, BYO-key field, copy-prompt button,
// Turnstile widget, context inspector, and run handler are defined once.
// `buildContext` is an async callback returning the context object to send.
function createAiReviewCard({ container, idSuffix, introText, loadingText, buildContext }) {
    const aiCard = document.createElement('div');
    aiCard.className = 'card ai-review-card';
    aiCard.setAttribute('data-card', 'optional-ai-review'); // mirrors applyCardTheme()
    const aiTitle = document.createElement('h3');
    aiTitle.textContent = 'Optional AI Review';
    aiCard.appendChild(aiTitle);

    const aiContent = document.createElement('div');
    aiContent.className = 'card-content ai-review-content';

    const aiIntro = document.createElement('p');
    aiIntro.className = 'ai-review-intro';
    aiIntro.textContent = introText;
    aiContent.appendChild(aiIntro);

    const controls = document.createElement('div');
    controls.className = 'ai-review-controls';

    const modelLabel = document.createElement('label');
    modelLabel.textContent = 'Model';
    modelLabel.setAttribute('for', `aiReviewModelSelect${idSuffix}`);
    controls.appendChild(modelLabel);

    const modelSelect = document.createElement('select');
    modelSelect.id = `aiReviewModelSelect${idSuffix}`;
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

    const copyPromptButton = document.createElement('button');
    copyPromptButton.type = 'button';
    copyPromptButton.className = 'ai-review-secondary-btn';
    copyPromptButton.textContent = 'Copy prompt';
    copyPromptButton.title = 'Copy the exact prompt so you can paste it into your own LLM';
    controls.appendChild(copyPromptButton);

    aiContent.appendChild(controls);

    const notesWrap = document.createElement('div');
    notesWrap.className = 'ai-review-notes';
    const notesLabel = document.createElement('label');
    notesLabel.setAttribute('for', `aiReviewNotes${idSuffix}`);
    notesLabel.textContent = 'Any additional notes for AI review';
    notesWrap.appendChild(notesLabel);
    const notesInput = document.createElement('textarea');
    notesInput.id = `aiReviewNotes${idSuffix}`;
    notesInput.rows = 3;
    notesInput.placeholder = 'Optional — extra clinical context, prior therapies, specific questions, etc.';
    notesWrap.appendChild(notesInput);
    aiContent.appendChild(notesWrap);

    // Optional bring-your-own OpenRouter key (collapsed by default).
    const keyDetails = document.createElement('details');
    keyDetails.className = 'ai-review-key';
    const keySummary = document.createElement('summary');
    keySummary.textContent = 'Use your own OpenRouter API key (optional)';
    keyDetails.appendChild(keySummary);
    const keyRow = document.createElement('div');
    keyRow.className = 'ai-review-key-row';
    const keyInput = document.createElement('input');
    keyInput.type = 'password';
    keyInput.id = `aiReviewKey${idSuffix}`;
    keyInput.autocomplete = 'off';
    keyInput.placeholder = 'sk-or-...';
    keyInput.value = getStoredUserKey();
    keyRow.appendChild(keyInput);
    const rememberWrap = document.createElement('label');
    rememberWrap.className = 'ai-review-key-remember';
    const rememberChk = document.createElement('input');
    rememberChk.type = 'checkbox';
    rememberChk.checked = Boolean(getStoredUserKey());
    rememberWrap.appendChild(rememberChk);
    rememberWrap.appendChild(document.createTextNode(' Remember in this tab'));
    keyRow.appendChild(rememberWrap);
    keyDetails.appendChild(keyRow);
    const keyNote = document.createElement('div');
    keyNote.className = 'ai-review-key-note';
    keyNote.textContent = 'Your key is sent only with your own review request and, if remembered, stored only in this browser tab (cleared when it closes). Requests using your own key are not rate-limited by this site.';
    keyDetails.appendChild(keyNote);
    if (getStoredUserKey()) keyDetails.open = true;
    aiContent.appendChild(keyDetails);

    // Turnstile widget mount point — only used when a site key is configured.
    const turnstileHost = document.createElement('div');
    turnstileHost.className = 'ai-review-turnstile';
    aiContent.appendChild(turnstileHost);
    let currentTurnstileToken = '';
    let turnstileWidgetId = null;
    if (window.TURNSTILE_SITE_KEY) {
        loadTurnstileScript().then(() => {
            if (!window.turnstile) return;
            try {
                turnstileWidgetId = window.turnstile.render(turnstileHost, {
                    sitekey: window.TURNSTILE_SITE_KEY,
                    callback: (token) => { currentTurnstileToken = token; },
                    'expired-callback': () => { currentTurnstileToken = ''; },
                    'error-callback': () => { currentTurnstileToken = ''; }
                });
            } catch { /* rendering failed — treat as no challenge */ }
        }).catch(() => { /* script blocked — backend skips when no secret set */ });
    }

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

    // Collapsible prompt output — filled by "Copy prompt" and shown for manual copy
    // when the clipboard API is unavailable/blocked.
    const promptDetails = document.createElement('details');
    promptDetails.className = 'ai-review-prompt-view';
    promptDetails.style.cssText = 'margin:2px 0;font-size:0.80rem;';
    const promptSummary = document.createElement('summary');
    promptSummary.style.cssText = 'cursor:pointer;color:#9ca3af;padding:2px 0;font-size:0.79rem;';
    promptSummary.textContent = 'Full prompt (for your own LLM)';
    promptDetails.appendChild(promptSummary);
    const promptTextarea = document.createElement('textarea');
    promptTextarea.readOnly = true;
    promptTextarea.rows = 8;
    promptTextarea.style.cssText = 'width:100%;box-sizing:border-box;font:inherit;font-size:0.73rem;line-height:1.4;background:#f8fafc;border:1px solid #e5e7eb;border-radius:4px;padding:8px;margin-top:4px;color:#374151;';
    promptTextarea.placeholder = 'Click "Copy prompt" to populate.';
    promptDetails.appendChild(promptTextarea);
    aiContent.appendChild(promptDetails);

    const aiOutput = document.createElement('div');
    aiOutput.className = 'ai-review-output';
    aiContent.appendChild(aiOutput);
    aiCard.appendChild(aiContent);
    if (container) container.appendChild(aiCard);

    const persistKey = () => {
        const key = (keyInput.value || '').trim();
        if (rememberChk.checked && key) setStoredUserKey(key);
        else setStoredUserKey('');
    };

    copyPromptButton.addEventListener('click', async () => {
        copyPromptButton.disabled = true;
        const prev = copyPromptButton.textContent;
        copyPromptButton.textContent = 'Building…';
        try {
            const aiContext = await buildContext((notesInput.value || '').trim());
            aiContextPre.textContent = JSON.stringify(aiContext, null, 2);
            aiContextInspectorSummary.textContent = 'Context sent to AI';
            const prompt = await fetchAiReviewPrompt(aiContext, modelSelect.value);
            promptTextarea.value = prompt;
            let copied = false;
            try {
                await navigator.clipboard.writeText(prompt);
                copied = true;
            } catch { /* clipboard blocked — fall back to manual copy below */ }
            if (copied) {
                copyPromptButton.textContent = 'Copied!';
            } else {
                promptDetails.open = true;
                promptTextarea.focus();
                promptTextarea.select();
                copyPromptButton.textContent = 'Copy from box below';
            }
        } catch (err) {
            copyPromptButton.textContent = 'Prompt failed';
            promptDetails.open = true;
            promptTextarea.value = `Unable to build prompt: ${err.message}`;
        } finally {
            setTimeout(() => { copyPromptButton.textContent = prev; }, 2000);
            copyPromptButton.disabled = false;
        }
    });

    runButton.addEventListener('click', async () => {
        runButton.disabled = true;
        const previousText = runButton.textContent;
        runButton.textContent = 'Running…';
        aiOutput.innerHTML = `<div class="ai-review-loading">${loadingText}</div>`;
        try {
            const userKey = (keyInput.value || '').trim();
            persistKey();
            // Require a Turnstile token only when a site key is configured and the user
            // is spending the owner key (BYO-key requests bypass the challenge).
            let turnstileToken = '';
            if (window.TURNSTILE_SITE_KEY && !userKey) {
                turnstileToken = currentTurnstileToken;
                if (!turnstileToken) throw new Error('Please complete the verification challenge above, then run again.');
            }
            const aiContext = await buildContext((notesInput.value || '').trim());
            aiContextPre.textContent = JSON.stringify(aiContext, null, 2);
            aiContextInspectorSummary.textContent = 'Context sent to AI';
            const data = await fetchAiReview(aiContext, modelSelect.value, {
                apiKey: userKey || undefined,
                turnstileToken: turnstileToken || undefined
            });
            renderAiReview(data.review, aiOutput);
            // Turnstile tokens are single-use; reset so a follow-up run gets a fresh one.
            if (window.turnstile && turnstileWidgetId != null) {
                try { window.turnstile.reset(turnstileWidgetId); } catch { /* ignore */ }
                currentTurnstileToken = '';
            }
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

    return aiCard;
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
            const href = safeUrl(value.url);
            if (href) {
                const a = document.createElement('a');
                a.href = href;
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


// resolvedAlleles: optional { chrom, pos, ref, alt } from resolveVcfAlleles().
function buildSpliceAiApiVariant(rawInput, gVariant, annotation, resolvedAlleles) {
    if (resolvedAlleles && resolvedAlleles.chrom && resolvedAlleles.pos && resolvedAlleles.ref && resolvedAlleles.alt) {
        return `chr${String(resolvedAlleles.chrom).replace(/^chr/i, '')}-${resolvedAlleles.pos}-${resolvedAlleles.ref}-${resolvedAlleles.alt}`;
    }
    const tuple = buildVariantCoordinateTuple(rawInput, gVariant);
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
    const response = await fetchWithRetry(url, {}, API_TIMEOUT_MS.myvariant, {
        attempts: 2,
        onAttempt: ({ attempt, attempts }) => LookupProgress.attempt('myvariant', attempt, attempts)
    });
    if (!response.ok) {
        const text = await response.text();
        const err = new Error(`MyVariant.info request failed (${response.status}): ${text}`);
        err.status = response.status;
        // 404 is the normal answer for an indel MyVariant does not index — it is
        // "no record here", not a failure, and must not be reported as an outage.
        err.notFound = response.status === 404;
        err.retryable = RETRYABLE_HTTP_STATUS.has(response.status);
        throw err;
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

    // Rank consequences so that a gene actually altered by the variant is preferred over a
    // neighbouring gene that VEP only reports as an upstream/downstream/intergenic variant.
    // Without this, a nearby gene (e.g. RRP9 listed as an upstream_gene_variant) can mask the
    // gene that truly carries the change (e.g. PARP3 with an inframe deletion).
    const isNeighbourOnly = (c) => {
        const terms = Array.isArray(c?.consequence_terms) ? c.consequence_terms : [];
        if (terms.length === 0) return false;
        return terms.every(t => /^(upstream_gene_variant|downstream_gene_variant|intergenic_variant)$/i.test(String(t)));
    };
    const hasCodingEffect = (c) => Boolean(c?.hgvsp || c?.amino_acids || c?.protein_start != null || c?.hgvsc);
    const scoreConsequence = (c) => (hasCodingEffect(c) ? 2 : 0) + (isNeighbourOnly(c) ? 0 : 1);
    const rankedConsequences = [...consequences].sort((a, b) => scoreConsequence(b) - scoreConsequence(a));

    // 1) Prefer any direct gene_symbol that is not chromosome-like, from the highest-ranked
    // (variant-overlapping) consequences first.
    for (const c of rankedConsequences) {
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
    // Request HGVS notations (hgvs=1) so hgvsc/hgvsp strings are returned per transcript, and
    // prefer the RefSeq transcript set (refseq=1). The rest of the app selects and displays the
    // canonical isoform by RefSeq NM_ accession, so RefSeq notations (e.g. NM_005485.6:c.52_54del,
    // NP_005476.4:p.Lys18del) match the nomenclature users expect. Fall back to the default
    // Ensembl transcript set when RefSeq returns no transcript consequences.
    const requestVep = async (useRefseq) => {
        const params = useRefseq ? 'hgvs=1&refseq=1' : 'hgvs=1';
        const url = `https://grch37.rest.ensembl.org/vep/human/hgvs/${encodeURIComponent(hgvs)}?content-type=application/json&${params}`;
        const response = await fetchWithRetry(url, {
            headers: {
                'Accept': 'application/json'
            }
        }, API_TIMEOUT_MS.vep, {
            attempts: 2,
            onAttempt: ({ attempt, attempts }) => LookupProgress.attempt('vep', attempt, attempts)
        });
        if (!response.ok) {
            const text = await response.text();
            const err = new Error(`VEP HGVS request failed (${response.status}): ${text}`);
            err.status = response.status;
            err.retryable = RETRYABLE_HTTP_STATUS.has(response.status);
            throw err;
        }
        const data = await response.json();
        if (!Array.isArray(data) || data.length === 0) {
            throw new Error('No VEP HGVS data found');
        }
        return data;
    };
    let data;
    try {
        data = await requestVep(true);
        const cons = (data[0] && data[0].transcript_consequences) || [];
        if (cons.length === 0) {
            data = await requestVep(false);
        }
    } catch (refseqErr) {
        // If the RefSeq request failed outright, fall back to the default Ensembl transcript set.
        data = await requestVep(false);
    }
    // Extract transcript consequences list from the first result
    const first = data[0];
    const consequences = first.transcript_consequences || [];
    return {
        vepData: data,
        consequences,
        mostSevere: first.most_severe_consequence || '',
        gVariant: buildGenomicHgvsFromVep(first)
    };
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
    const response = await fetchWithRetry(url, {}, API_TIMEOUT_MS.myvariant, {
        attempts: 2,
        onAttempt: ({ attempt, attempts }) => LookupProgress.attempt('myvariant', attempt, attempts)
    });
    if (!response.ok) {
        const text = await response.text();
        const err = new Error(`MyVariant.info query request failed (${response.status}): ${text}`);
        err.status = response.status;
        err.retryable = RETRYABLE_HTTP_STATUS.has(response.status);
        throw err;
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
            if (c.source && c.source.url) {
                const srcUrl = safeUrl(c.source.url);
                if (srcUrl) civic['CIViC Source'] = { html: `<a href="${escapeHtml(srcUrl)}" target="_blank" rel="noopener noreferrer">${escapeHtml(c.source.url)}</a>` };
            }
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
    // Raw variant_recoder response from the most recent call, keyed by the query
    // it answered. getTranscriptsList() runs at the start of every lookup, and
    // the pipeline used to call fetchVariantRecoder AGAIN for the same query a
    // moment later — a duplicated 12s-class (40s+ cold) Ensembl round-trip on
    // most non-genomic lookups. Cache both success and failure so the second
    // consumer reuses the answer instead of re-asking (or re-retrying) Ensembl.
    let lastRecoderResult = { query: null, data: null, error: null };


    /**
     * Fetch transcripts for a given variant using the Ensembl variant recoder. Returns an array of
     * objects with transcript, cDNA and protein fields. If no transcripts can be obtained, returns an empty array.
     * @param {string} q Variant query (e.g. protein HGVS, cDNA or genomic variant).
     */
    async function getTranscriptsList(q) {
        try {
            const recResults = await fetchVariantRecoder(q);
            lastRecoderResult = { query: q, data: recResults, error: null };
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
            // Record it though: an empty transcript list because Ensembl is down
            // must not be presented the same way as one because the variant is
            // genuinely unknown.
            lastRecoderResult = { query: q, data: null, error: err };
            LookupProgress.step('recoder', 'fail', describeUpstreamFailure(err));
            if (err && err.retryable) LookupProgress.markUpstreamOutage('recoder');
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

    // Monotonic search counter, shared across submits, for the overlapping-search guard.
    let searchSeq = 0;
    form.addEventListener('submit', async (e) => {
        e.preventDefault();
        // Guard against overlapping searches. Each submit takes the next sequence
        // number; any render/append that runs after a newer submit started is stale and
        // bails via isCurrentSearch(). Without this, a slow earlier search can resolve
        // after a newer one and overwrite the newer results with stale cards.
        const mySeq = ++searchSeq;
        const isCurrentSearch = () => mySeq === searchSeq;
        // Reset any previous gene hint.  Without resetting, a prior search that
        // included a gene symbol could incorrectly influence the next query.
        geneHintGlobal = null;
        alterationTypeGlobal = null;
        isGeneOnlyMode = false;
        codonOnlyResolutionGlobal = null;

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
            // Handle pasted coordinate rows: "chr7 140453136 A T", "7 140453136 BRAF A T",
            // MAF column order, CSV separators and trailing build columns all land here.
            {
                const row = parseCoordinateRow(s);
                if (row) {
                    if (row.gene && !isChromosomeLikeGeneSymbol(row.gene)) {
                        geneHintGlobal = row.gene;
                    }
                    const hgvs = coordinateRowToGenomicHgvs(row);
                    if (hgvs) return hgvs;
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
            // A row that looks like coordinates but did not parse must not fall through
            // to the "GENE p.Change" branch below — that is what turned
            // "EIF1AX X 20148674 T TG" into the query "EIF1AX:p.X" and then reported
            // it as a variant that does not exist. Hand the raw text back so the
            // caller can say what is actually wrong with the row.
            if (looksLikeCoordinateRow(s)) {
                return s;
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

        LookupProgress.start();
        // A coordinate row that survived normalisation unchanged is one the parser
        // could not read. Say exactly that instead of sending the raw text to the
        // annotation APIs and reporting the miss as "variant not found".
        if (looksLikeCoordinateRow(query)) {
            LookupProgress.step('parse', 'fail', 'unrecognised column layout');
            LookupProgress.finish('fail', 'Could not read that row');
            statusEl.textContent = 'Error: could not read that as a variant row. Expected chromosome, '
                + 'position, then reference and alternate alleles — for example "X 20148674 EIF1AX T TG". '
                + 'HGVS also works: chrX:g.20148675dup or NM_001412.4:c.388dup.';
            resultSection.classList.add('hidden');
            return;
        }
        LookupProgress.step('parse', 'ok', query !== rawInput ? `read as ${query}` : 'HGVS input');
        if (query !== rawInput) {
            LookupProgress.resolved({ input: rawInput, genomic: isGenomicVariant(query) ? query : '', gene: geneHintGlobal });
        }

        // ── Gene-only / gene+descriptor mode ────────────────────────────────
        // When the input is just a gene symbol or a gene + alteration descriptor
        // (e.g. "BRAF" or "BRAF amplification"), skip the variant annotation
        // pipeline entirely and render a focused set of gene-level cards instead.
        if (isGeneOnlyMode && geneHintGlobal) {
            const gene = geneHintGlobal;
            const altType = alterationTypeGlobal; // e.g. "amplification" or null

            // Gene-level mode renders its own cards and runs no variant pipeline,
            // so the per-source panel has nothing to report.
            LookupProgress.hide();
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
            // Card results reused by the AI-review buildContext but not spread into
            // supplemental_card_data (they appear top-level in the payload already).
            const geneOnlyCardCache = {};

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
                        civicApiDiv.innerHTML = `<div style="font-size:0.82rem;color:#9ca3af;">Gene "${escapeHtml(gene)}" not found in CIViC.</div>`;
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
                            ampEl.innerHTML = `<strong>AMP/ACMG tier (CIViC):</strong> ${escapeHtml(ampLevel)}`;
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
                const pmSearchTerm = altType ? `${gene} ${altType}` : gene;
                const pmTumorSearchTerm = tumorType ? `${gene} ${tumorType}` : '';
                const pmVariantTumorSearchTerm = (tumorType && altType) ? `${gene} ${altType} ${tumorType}` : '';
                const pmTabs = [{ label: altType ? 'Gene + Alteration' : 'Gene', term: pmSearchTerm, extraKey: 'pubmed' }];
                if (pmTumorSearchTerm && pmSearchTerm !== pmTumorSearchTerm) {
                    pmTabs.push({ label: 'Gene + Tumor Type', term: pmTumorSearchTerm, extraKey: 'pubmed_tumor_type' });
                    if (pmVariantTumorSearchTerm
                        && pmVariantTumorSearchTerm !== pmSearchTerm
                        && pmVariantTumorSearchTerm !== pmTumorSearchTerm) {
                        pmTabs.push({ label: 'Gene + Alteration + Tumor Type', term: pmVariantTumorSearchTerm, extraKey: 'pubmed_variant_tumor_type' });
                    }
                }
                createPubmedCard({ container: cardsContainer, tabs: pmTabs, extras: geneOnlyAiExtras });
            }

            // ── Card: FDA-Approved Drugs (by gene) ─────────────────────────
            createFdaDrugsCard({ container: cardsContainer, gene, extras: geneOnlyAiExtras, cardCache: geneOnlyCardCache });

            // ── Card: Clinical Trials ──────────────────────────────────────
            createClinicalTrialsCard({ container: cardsContainer, gene, tumorType, cardCache: geneOnlyCardCache });

            // ── Card: Cancer Prevalence (cBioPortal cohort frequency) ──────
            createCbioportalCard({ container: cardsContainer, gene, proteinChange: '', extras: geneOnlyAiExtras });

            // ── Card: Guidelines ───────────────────────────────────────────
            createGuidelinesCard({ container: cardsContainer, tumorType, getGene: () => gene, extras: geneOnlyAiExtras });

            // ── Card: Optional AI Review (gene-only mode) ──────────────────
            createAiReviewCard({
                container: cardsContainer,
                idSuffix: 'Gene',
                introText: 'Send the retrieved gene-level data to OpenRouter for a structured draft interpretation. No request is sent until you click Run AI review.',
                loadingText: 'Gathering gene-level context for AI review…',
                buildContext: async (userNotes) => {
                    // Supplemental context for gene-only mode: skip coordinate-dependent
                    // calls, and reuse what the cards already fetched — a fresh fetch is
                    // the fallback, not the default.
                    const fdaPromise = Array.isArray(geneOnlyCardCache.fda_companion_diagnostics_records)
                        ? Promise.resolve(geneOnlyCardCache.fda_companion_diagnostics_records)
                        : fetchFdaCompanionDiagnostics(gene).catch(() => []);
                    const bbkbPromise = (geneOnlyAiExtras.bbkb_biomarker_therapies && Array.isArray(geneOnlyAiExtras.bbkb_biomarker_therapies.results))
                        ? Promise.resolve(geneOnlyAiExtras.bbkb_biomarker_therapies)
                        : fetchBbkbBiomarkerTherapies(gene).catch(() => ({ total_matched: 0, returned: 0, results: [] }));
                    const trialsPromise = (geneOnlyCardCache.clinical_trials && Array.isArray(geneOnlyCardCache.clinical_trials.studies))
                        ? Promise.resolve(geneOnlyCardCache.clinical_trials)
                        : fetchClinicalTrials(gene, tumorType).catch(() => ({ total: 0, studies: [] }));
                    const [fdaRecords, bbkbTherapies, clinicalTrialData] = await Promise.all([
                        fdaPromise,
                        bbkbPromise,
                        trialsPromise
                    ]);
                    const pubmedTerm = altType ? `${gene} ${altType}` : gene;
                    // Prefer the PubMed card's already-fetched data over re-fetching.
                    // Re-fetching alongside the card increases NCBI concurrency and can
                    // drop abstracts when efetch hits the rate limit.
                    const cachedPubmed = geneOnlyAiExtras.pubmed;
                    const pubmedPromise = (cachedPubmed && Array.isArray(cachedPubmed.articles) && cachedPubmed.articles.length > 0)
                        ? Promise.resolve({ total: cachedPubmed.total ?? cachedPubmed.articles.length, articles: cachedPubmed.articles })
                        : fetchPubmedArticles(pubmedTerm, 5).catch(() => ({ total: 0, articles: [] }));
                    const civicPromise = geneOnlyAiExtras.civic_api
                        ? Promise.resolve(geneOnlyAiExtras.civic_api)
                        : fetchCivicApiData(gene, '').catch(() => null);
                    const openFdaPromise = (geneOnlyAiExtras.openfda && Array.isArray(geneOnlyAiExtras.openfda.results))
                        ? Promise.resolve(geneOnlyAiExtras.openfda)
                        : fetchOpenFdaDrugLabels(gene).catch(() => null);
                    const [pubmedData, civicData, openFdaData] = await Promise.all([
                        pubmedPromise,
                        civicPromise,
                        openFdaPromise
                    ]);
                    const supplementalContext = {
                        ...geneOnlyAiExtras,
                        civic_api: civicData,
                        pubmed: pubmedData,
                        openfda_drug_labels: condenseOpenFdaForAi(openFdaData)
                    };
                    // The card's copy under `openfda` would duplicate `openfda_drug_labels`
                    // in the payload — for drug-rich genes that is the dominant context,
                    // doubled. Keep the one canonical key.
                    if (supplementalContext.openfda_drug_labels && supplementalContext.openfda) {
                        delete supplementalContext.openfda;
                    } else if (supplementalContext.openfda) {
                        supplementalContext.openfda = condenseOpenFdaForAi(supplementalContext.openfda);
                    }
                    return {
                        submitted_query: rawInput,
                        gene,
                        alteration_type: altType,
                        tumor_type: tumorType,
                        user_notes: userNotes || undefined,
                        fda_companion_diagnostics_records: fdaRecords,
                        bbkb_biomarker_therapies: bbkbTherapies,
                        clinical_trials: clinicalTrialData,
                        supplemental_card_data: supplementalContext
                    };
                }
            });

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
        LookupProgress.step('recoder', 'running');
        // Start the recoder WITHOUT awaiting it. The whole pipeline used to block
        // on this call before even trying the direct MyVariant search — during an
        // Ensembl outage that meant ~37s of retries/timeouts in front of a result
        // that the direct search then produced in ~2s. The promise is awaited only
        // where its answer is actually needed: in full before the recoder-candidate
        // fallback, and with a short grace period before rendering a result that
        // was resolved without it (so a healthy Ensembl still contributes RefSeq
        // transcript nomenclature, but a sick one cannot hold the page hostage).
        transcriptsFromRecoder = [];
        const transcriptsPromise = getTranscriptsList(query).then((list) => {
            // Guard the shared state: a stale search's late-resolving recoder must
            // not overwrite the transcripts or progress panel of a newer search.
            if (!isCurrentSearch()) return list;
            transcriptsFromRecoder = list;
            LookupProgress.stepUnlessFailed('recoder', list.length ? 'ok' : 'empty',
                list.length
                    ? `${list.length} transcript${list.length === 1 ? '' : 's'}`
                    : 'no transcript match');
            return list;
        }).catch((recoderErr) => {
            // getTranscriptsList swallows its own errors, so this is a safety net.
            if (!isCurrentSearch()) return [];
            transcriptsFromRecoder = [];
            LookupProgress.step('recoder', 'fail', describeUpstreamFailure(recoderErr));
            if (recoderErr && recoderErr.retryable) LookupProgress.markUpstreamOutage('recoder');
            return [];
        });
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
                LookupProgress.step('myvariant', 'running');
                try {
                    annotation = await fetchMyVariant(query);
                    gVariant = query;
                    LookupProgress.step('myvariant', 'ok', 'annotation found');
                    console.log('[DEBUG] Direct genomic fetch succeeded for', query);
                } catch (errDirectGenomic) {
                    // 404 is the routine answer for an indel MyVariant does not index.
                    LookupProgress.step('myvariant',
                        errDirectGenomic && errDirectGenomic.notFound ? 'empty' : 'fail',
                        errDirectGenomic && errDirectGenomic.notFound
                            ? 'no record for this variant'
                            : describeUpstreamFailure(errDirectGenomic));
                    if (errDirectGenomic && errDirectGenomic.retryable && !errDirectGenomic.notFound) {
                        LookupProgress.markUpstreamOutage('myvariant');
                    }
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
                // Nothing resolved without the recoder, so now its answer is genuinely
                // needed — wait for the in-flight call (this also guarantees
                // lastRecoderResult below reflects this query's attempt).
                await transcriptsPromise;
                // getTranscriptsList already asked the recoder for this exact query at
                // the start of the lookup (and fetchWithRetry exhausted its retries if
                // Ensembl was unwell) — reuse that answer instead of a second 12s-class
                // round trip. A cached failure stays a failure: recoderData remains
                // null and the free-text MyVariant fallback below takes over.
                if (lastRecoderResult.query === query) {
                    recoderData = lastRecoderResult.data;
                    if (!recoderData && lastRecoderResult.error) {
                        console.warn('Variant recoder failed earlier in this lookup; skipping retry', lastRecoderResult.error);
                    }
                } else {
                    try {
                        recoderData = await fetchVariantRecoder(query);
                    } catch (errRecoder) {
                        // Warn in console but don't immediately fail; we will attempt a free-text MyVariant search below.
                        console.warn('Variant recoder failed; will attempt free-text search', errRecoder);
                    }
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
                // Convert each hgvsg/SPDI entry to MyVariant notation, lifting to hg19
                // ONLY when its accession says GRCh38. The GRCh37 recoder mirror emits
                // GRCh37 accessions (NC_000007.13), and the previous unconditional
                // liftoverHg38ToHg19() "mapped" those already-hg19 positions to a locus
                // hundreds of kb away — Ensembl's /map does not error on a coordinate
                // from the wrong assembly, it just answers wrongly. Every substitution
                // candidate was corrupted this way (indels escaped only because the
                // liftover helper's regex matches substitutions alone). Unknown
                // accessions conservatively keep the old lift-it behaviour.
                const convertCandidate = async (raw, converter, label) => {
                    try {
                        const mv = converter(raw);
                        const assembly = assemblyFromNcAccession(raw);
                        if (assembly === 'GRCh37') {
                            console.log(`[DEBUG] Converted ${label} to MV:`, raw, '->', mv, '(already GRCh37)');
                            return mv;
                        }
                        const lifted = await liftoverHg38ToHg19(mv);
                        // A known-GRCh38 coordinate that came back unchanged was NOT
                        // converted (the helper only lifts substitutions, and returns
                        // its input on failure) — drop it rather than carry an hg38
                        // position mislabelled as hg19 into the cards.
                        if (assembly === 'GRCh38' && lifted === mv) {
                            console.log(`[DEBUG] Dropping unliftable GRCh38 ${label} candidate:`, raw);
                            return null;
                        }
                        console.log(`[DEBUG] Converted ${label} to MV and liftover:`, raw, '->', mv, '->', lifted);
                        return lifted;
                    } catch {
                        return null; // skip conversion errors
                    }
                };
                {
                    // Collect conversions first, then run them concurrently — the lifted
                    // (GRCh38) entries each cost an Ensembl /map round-trip and used to
                    // run strictly in series.
                    const conversions = [];
                    for (const item of recoderData) {
                        for (const key in item) {
                            const sub = item[key];
                            if (!sub) continue;
                            // hgvsg / spdi entries may be a string or an array
                            if (sub.hgvsg) {
                                const hgvsgList = Array.isArray(sub.hgvsg) ? sub.hgvsg : [sub.hgvsg];
                                for (const hgvsg of hgvsgList) {
                                    conversions.push(convertCandidate(hgvsg, convertHgvsgToMyVariant, 'hgvsg'));
                                }
                            }
                            if (sub.spdi) {
                                const spdiList = Array.isArray(sub.spdi) ? sub.spdi : [sub.spdi];
                                for (const spdi of spdiList) {
                                    conversions.push(convertCandidate(spdi, convertSpdiToMyVariant, 'SPDI'));
                                }
                            }
                        }
                    }
                    candidateVariants.push(...(await Promise.all(conversions)).filter(Boolean));
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
                    // Fetch every candidate's annotation concurrently (typically 1-5);
                    // the selection walk below still runs in the recoder's priority
                    // order, so the outcome matches the old sequential scan minus the
                    // serial MyVariant round-trips.
                    const candidateAnns = await Promise.all(uniqueCandidates.map(async (cand) => {
                        try {
                            const fetched = await fetchMyVariant(cand);
                            console.log('[DEBUG] Fetched annotation for candidate', cand);
                            return fetched;
                        } catch {
                            console.log('[DEBUG] Annotation fetch error for candidate', cand);
                            return null;
                        }
                    }));
                    for (let candIdx = 0; candIdx < uniqueCandidates.length; candIdx++) {
                        const cand = uniqueCandidates[candIdx];
                        const ann = candidateAnns[candIdx];
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
                        // MyVariant.info has nothing for this variant, so the summary would
                        // show "Effect: N/A" — common for indels it does not index. VEP knows
                        // the consequence for any variant the recoder could resolve, so ask it
                        // rather than leaving the field blank. It also reports the genomic
                        // location, which recovers g. when the recoder gave us no candidate.
                        try {
                            const vepRes = await fetchVepHgvsHg19(isGenomicVariant(gVariant) ? gVariant : query);
                            if (vepRes && vepRes.mostSevere) {
                                annotation.cadd = Object.assign({}, annotation.cadd, { consequence: vepRes.mostSevere });
                            }
                            if (!isGenomicVariant(gVariant) && vepRes && vepRes.gVariant) {
                                gVariant = vepRes.gVariant;
                            }
                        } catch (vepEffectErr) {
                            console.log('[DEBUG] VEP consequence lookup for minimal annotation failed:', vepEffectErr);
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
                                LookupProgress.step('vep', 'running');
                                try {
                                    // Use the currently normalised gVariant if available; otherwise fall back to the original query.
                                    const hgvsForVep = gVariant || query;
                                    const vepRes = await fetchVepHgvsHg19(hgvsForVep);
                                    const vepConsequences = vepRes.consequences || [];
                                    LookupProgress.step('vep', vepConsequences.length ? 'ok' : 'empty',
                                        vepConsequences.length ? vepRes.mostSevere || 'resolved' : 'no consequences');
                                    if (vepConsequences.length > 0) {
                                        // Build a minimal annotation using the first gene symbol from the VEP data. The dbnsfp
                                        // genename field is used for summary display. The _id is set to the original query.
                                        annotation = { _id: query };
                                        // VEP resolved the variant, so it knows where it is. Without this
                                        // the g. field keeps echoing the user's cDNA query and every
                                        // coordinate-derived card (UCSC, ClinVar region, gnomAD, SpliceAI)
                                        // stays dark even though the variant was successfully resolved.
                                        if (!isGenomicVariant(gVariant) && vepRes.gVariant) {
                                            gVariant = vepRes.gVariant;
                                        }
                                        let geneSym = await resolveGeneSymbolFromVep(vepConsequences, recoderData);
                                        // If upstream sources still do not provide a usable gene symbol,
                                        // then use the optional user hint as the final fallback.
                                        if ((!geneSym || isChromosomeLikeGeneSymbol(geneSym)) && geneHintGlobal) {
                                            geneSym = geneHintGlobal;
                                        }
                                        // Last resort: a gene-prefixed query ("TSC2:c.4952delA") names the
                                        // gene outright, which beats leaving the card blank.
                                        if (!geneSym || isChromosomeLikeGeneSymbol(geneSym)) {
                                            const queryGene = (query.match(/^([A-Za-z][A-Za-z0-9-]*)\s*:/) || [])[1];
                                            if (queryGene && !isChromosomeLikeGeneSymbol(queryGene)) geneSym = queryGene.toUpperCase();
                                        }
                                        if (geneSym) {
                                            annotation.dbnsfp = { genename: geneSym };
                                        }
                                        // Expose the most severe consequence so the summary can display an Effect
                                        // (e.g. "inframe_deletion") instead of N/A for VEP-only variants.
                                        if (vepRes.mostSevere) {
                                            annotation.cadd = Object.assign({}, annotation.cadd, { consequence: vepRes.mostSevere });
                                        }
                                        // Populate transcriptsFromRecoder with the cDNA (hgvsc) and protein (hgvsp)
                                        // notations returned by VEP so the summary can display c./p. and select a
                                        // canonical isoform. Split the "accession:change" strings into a transcript
                                        // accession plus the c./p. portion, matching the shape used elsewhere. Skip
                                        // consequences without an hgvsc (e.g. neighbouring upstream/downstream genes).
                                        transcriptsFromRecoder = vepConsequences
                                            .filter(c => c && c.hgvsc)
                                            .map(c => {
                                                const cParts = String(c.hgvsc).split(':');
                                                const transcriptId = cParts[0] || (c.transcript_id || '');
                                                const cDNA = cParts.slice(1).join(':');
                                                let protein = '';
                                                if (c.hgvsp) {
                                                    const pParts = String(c.hgvsp).split(':');
                                                    protein = pParts.slice(1).join(':');
                                                }
                                                return {
                                                    transcript: transcriptId,
                                                    cDNA: cDNA,
                                                    protein: protein,
                                                    source: 'vep'
                                                };
                                            });
                                        altFound = true;
                                    }
                                } catch (vepErr) {
                                    // Record the reason before falling through: the error thrown
                                    // below depends on whether this was an outage or a real miss.
                                    LookupProgress.step('vep', 'fail', describeUpstreamFailure(vepErr));
                                    if (vepErr && vepErr.retryable) LookupProgress.markUpstreamOutage('vep');
                                }
                                if (!altFound) {
                                    // Nothing resolved. Which message is correct depends entirely on
                                    // *why* nothing resolved, so ask the progress panel rather than
                                    // assuming the user mistyped a coordinate: an Ensembl 503 used to
                                    // be reported as an invalid reference allele.
                                    throw buildLookupFailureError();
                                }
                            }
                        }
                    } catch (errSearch) {
                        console.log('[DEBUG] Free-text MyVariant search error:', errSearch);
                        // Re-wrapping used to drop the outage flag, so an Ensembl 503 was
                        // reported to the user as a plain lookup failure. Keep the original
                        // error when it already carries a diagnosis.
                        if (errSearch && (errSearch.upstreamOutage || errSearch.retryable)) throw errSearch;
                        throw new Error(errSearch.message || 'Variant not found');
                    }
                } else {
                    // Protein-only query the recoder could not resolve. This is not a
                    // transient failure for a whole class of input: Ensembl rejects
                    // frameshifts outright ("Frameshifts are not supported for HGVS
                    // protein input") and cannot derive a nucleotide change for a
                    // multi-residue delins — because none is determined. p.N1651Mfs*21
                    // can arise from many different indels.
                    //
                    // The codon *is* determined, so resolve that much rather than failing
                    // the lookup outright: the user still gets the right locus, the UCSC
                    // link, ClinVar's variants in that region, and every gene-level card.
                    // The allele-dependent cards stay empty, which is the honest answer.
                    const proteinBody = (query.match(/p\.[^\s:]+/i) || [])[0] || '';
                    const residues = parseProteinResidueRange(proteinBody);
                    const geneForCodon = geneHintGlobal
                        || (query.match(/^([A-Za-z][A-Za-z0-9-]*)\s*:/) || [])[1]
                        || '';
                    let codonRegion = null;
                    if (residues && geneForCodon) {
                        statusEl.textContent = 'Resolving codon position…';
                        codonRegion = await resolveProteinCodonRegion(geneForCodon, residues.start, residues.end);
                    }
                    if (!codonRegion) {
                        throw new Error(
                            `Could not resolve ${query} to a genomic position. Protein notation does not `
                            + 'determine a nucleotide change — Ensembl does not support frameshift protein '
                            + 'input, and a delins spanning several residues has no unique coding change. '
                            + 'Try the cDNA (c.) or genomic (g.) notation for this variant.'
                        );
                    }
                    // A bare span, not a variant: parseGenomicHgvs reads this as type
                    // 'region', so position-based cards light up and allele-based ones
                    // correctly report nothing.
                    gVariant = `chr${codonRegion.chrom.replace(/^chr/i, '')}:g.${codonRegion.start}_${codonRegion.end}`;
                    codonOnlyResolutionGlobal = {
                        proteinChange: proteinBody,
                        residueStart: residues.start,
                        residueEnd: residues.end,
                        transcript: codonRegion.transcript,
                        protein: codonRegion.protein
                    };
                    annotation = { _id: query };
                    if (geneForCodon) annotation.dbnsfp = { genename: geneForCodon.toUpperCase() };
                    transcriptsFromRecoder = [{
                        transcript: codonRegion.transcript,
                        cDNA: '',
                        protein: proteinBody
                    }];
                }
            }
            // The recoder may still be in flight when the annotation was resolved
            // without it (direct MyVariant hit). Give it a short grace so a healthy
            // Ensembl still contributes RefSeq transcript nomenclature to the
            // Variant card, but never let a sick one hold the finished result page
            // hostage — the dbNSFP fallback renders the correct canonical without it.
            await Promise.race([
                transcriptsPromise,
                new Promise((resolve) => setTimeout(resolve, 6000))
            ]);
            // A newer search may have started while this one's annotation was in
            // flight; from here on every write hits shared UI (progress panel,
            // status line, cards), so a stale run must stop before it overwrites
            // the newer lookup's output.
            if (!isCurrentSearch()) return;
            // Display annotation
            statusEl.textContent = 'Annotation retrieved';
            // Show what the input actually resolved to. Users paste a coordinate row
            // and get back a page of cards; without this they cannot tell which
            // transcript, which HGVS form, or even which build the tool settled on.
            {
                const canonical = Array.isArray(transcriptsFromRecoder) && transcriptsFromRecoder.length
                    ? transcriptsFromRecoder[selectCanonicalTranscriptSafe(transcriptsFromRecoder)]
                    : null;
                const summaryForPanel = buildSummary(annotation, gVariant);
                const rowValue = (name) => {
                    const row = summaryForPanel.find((r) => r.name === name);
                    return row ? String(row.value) : '';
                };
                LookupProgress.resolved({
                    input: rawInput,
                    genomic: isGenomicVariant(gVariant) ? gVariant : '',
                    cdna: canonical && canonical.cDNA ? `${canonical.transcript}:${canonical.cDNA}` : '',
                    protein: canonical && canonical.protein ? canonical.protein : '',
                    // The summary's gene row can repeat a symbol once per transcript
                    // ("BRAF, BRAF,BRAF"); the panel wants it named once.
                    gene: dedupeSymbolList(rowValue('Gene(s)')) || geneHintGlobal || '',
                    consequence: rowValue('Consequence'),
                    assembly: 'GRCh37/hg19'
                });
                LookupProgress.finish('ok', 'Variant resolved');
                // The panel header now carries the outcome, so the legacy status
                // line below it would just repeat "Annotation retrieved".
                statusEl.textContent = '';
            }
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
            // Card results reused by the AI-review buildContext but NOT spread into
            // supplemental_card_data (they already appear top-level in the payload,
            // so keeping them out of aiReviewExtras avoids sending two copies).
            const cardDataCache = {};

            // GRCh37 VCF-style alleles for this lookup, resolved once and shared by
            // every allele-dependent consumer (gnomAD v2 link, gnomAD v4 API, SpliceAI,
            // AI review context). For indels whose HGVS notation carries no alleles
            // this reads the reference sequence, so the promise is memoised here rather
            // than re-derived per card.
            let vcfAllelesPromise = null;
            const getResolvedVcfAlleles = () => {
                if (!vcfAllelesPromise) {
                    vcfAllelesPromise = resolveVcfAlleles(rawInput, gVariant, annotation).catch(() => null);
                }
                return vcfAllelesPromise;
            };

            // Attempt to fetch extended COSMIC data from a custom API if configured.
            const COSMIC_ENDPOINT = window.COSMIC_API_ENDPOINT || null;
            const COSMIC_META_URL = window.COSMIC_META_URL || null;
            // A codon-region pseudo-variant is not an hgvsg COSMIC can match.
            if (COSMIC_ENDPOINT && !codonOnlyResolutionGlobal) {
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
                            cosmicItems['COSMIC Gene Page'] = { html: `<a href="${geneLink}" target="_blank" rel="noopener noreferrer">View analysis for ${escapeHtml(cosmicData.COSMIC_GENE)}</a>` };
                        }
                        // Site counts with per-type frequencies and gene-specific frequencies
                        if (cosmicData.COSMIC_SITE_COUNTS) {
                            const siteRows = [];
                            for (const [type, info] of Object.entries(cosmicData.COSMIC_SITE_COUNTS)) {
                                const count = info.count || 0;
                                const samplesWithGeneType = info.samples_with_gene_in_type || 1;
                                let rowText = `${escapeHtml(type)}: ${count} tumor${count === 1 ? '' : 's'}`;
                                if (meta && meta.total_samples_by_cancer_type && meta.total_samples_by_cancer_type[type]) {
                                    const typeTotal = meta.total_samples_by_cancer_type[type] || 1;
                                    const typeFreq = ((count / typeTotal) * 100).toFixed(2);
                                    const geneFreqType = ((count / samplesWithGeneType) * 100).toFixed(2);
                                    rowText += ` (${typeFreq}% of ${escapeHtml(type)}, ${geneFreqType}% with ${escapeHtml(geneNameForDisplay)})`;
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

            // Build and display cards summarising key annotations. Bail if a newer
            // search started while this one's annotation was being fetched — otherwise
            // this stale run would clear and rebuild over the newer results.
            if (!isCurrentSearch()) return;
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
                    // Build a list of transcript mappings (transcript ID -> cDNA, protein).
                    //
                    // MyVariant's dbnsfp.hgvsc entries are BARE c. strings with no
                    // accession ("c.1799T>A"), so splitting on ':' used to shove the
                    // whole string into the transcript slot and leave cDNA empty — the
                    // Variant card's canonical entry then had no cDNA, the display fell
                    // back to the FIRST list entry, and BRAF V600E rendered as the
                    // alternate-isoform c.620T>A whenever the recoder was down. Treat an
                    // accession-less entry as the cDNA itself. hgvsp is index-paired
                    // only when the arrays actually align — dbNSFP's hgvsc and hgvsp
                    // arrays can have different lengths (V600E: 3 vs 4), and pairing
                    // misaligned arrays attaches the wrong protein to a transcript.
                    const hgvspListLocal = (annotation.dbnsfp && annotation.dbnsfp.hgvsp) ? (Array.isArray(annotation.dbnsfp.hgvsp) ? annotation.dbnsfp.hgvsp : [annotation.dbnsfp.hgvsp]) : [];
                    const hgvspAligned = hgvspListLocal.length === hgvscList.length;
                    for (let i = 0; i < hgvscList.length; i++) {
                        const sc = hgvscList[i];
                        if (!sc) continue;
                        const scStr = String(sc);
                        const scColon = scStr.indexOf(':');
                        const transcriptId = scColon !== -1 ? scStr.slice(0, scColon) : '';
                        const cpart = scColon !== -1 ? scStr.slice(scColon + 1) : scStr;
                        let ppart = '';
                        if (hgvspAligned && hgvspListLocal[i]) {
                            const spStr = String(hgvspListLocal[i]);
                            const spColon = spStr.indexOf(':');
                            ppart = spColon !== -1 ? spStr.slice(spColon + 1) : spStr;
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
                    span.innerHTML = `<strong>${escapeHtml(label)}:</strong> ${value ? escapeHtml(value) : 'N/A'}`;
                    return span;
                };
                if (codonOnlyResolutionGlobal) {
                    // Not a variant coordinate — the span of the affected codon(s).
                    content.appendChild(makeLine('Codon region (hg19)', gVariant));
                    const codonNote = document.createElement('span');
                    codonNote.style.cssText = 'font-size:0.78rem;color:#92400e;background:#fef3c7;padding:3px 6px;border-radius:4px;';
                    codonNote.textContent = `Resolved to residue ${codonOnlyResolutionGlobal.residueStart}`
                        + (codonOnlyResolutionGlobal.residueEnd !== codonOnlyResolutionGlobal.residueStart
                            ? `–${codonOnlyResolutionGlobal.residueEnd}` : '')
                        + ` on ${codonOnlyResolutionGlobal.transcript} only. Protein notation does not determine the`
                        + ' nucleotide change, so allele-based results (gnomAD, SpliceAI, an exact ClinVar record)'
                        + ' are unavailable — enter the c. or g. notation for those.';
                    content.appendChild(codonNote);
                } else {
                    content.appendChild(makeLine('g.', gVariant));
                }
                // VCF-style CHROM-POS-REF-ALT (GRCh37) — the form gnomAD, SpliceAI and
                // most VCF-derived tools index by, and the one to paste into them.
                // Resolved asynchronously because an indel's alleles are read off the
                // reference; rendered as a placeholder so the line keeps its position.
                const vcfLine = document.createElement('span');
                vcfLine.innerHTML = '<strong>VCF (hg19):</strong> <span style="color:#9ca3af;">resolving…</span>';
                content.appendChild(vcfLine);
                getResolvedVcfAlleles().then((alleles) => {
                    if (!alleles || !alleles.ref || !alleles.alt) {
                        vcfLine.innerHTML = '<strong>VCF (hg19):</strong> N/A';
                        return;
                    }
                    const vcfId = `${alleles.chrom}-${alleles.pos}-${alleles.ref}-${alleles.alt}`;
                    vcfLine.innerHTML = `<strong>VCF (hg19):</strong> ${escapeHtml(vcfId)}`;
                    // Only flag the shift when the position actually moved — otherwise the
                    // caveat reads as a warning about a variant that has nothing to warn about.
                    const gStart = parseGenomicHgvs(gVariant)?.start;
                    if (alleles.source === 'reference' && gStart !== undefined && Number(alleles.pos) !== gStart) {
                        // Explain the coordinate mismatch a reader would otherwise take
                        // for a bug: HGVS shifts indels 3′, VCF left-aligns them.
                        vcfLine.title = 'Left-aligned against the GRCh37 reference. HGVS places an indel at its 3′-most position while VCF (and gnomAD) use the 5′-most one, so this position can differ from the g. notation above.';
                        const note = document.createElement('span');
                        note.style.cssText = 'font-size:0.75rem;color:#6b7280;';
                        note.textContent = '(left-aligned)';
                        vcfLine.appendChild(document.createTextNode(' '));
                        vcfLine.appendChild(note);
                    }
                });
                content.appendChild(makeLine('Gene', geneNames));
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
                    // The canonical entry is the <strong>-wrapped one — prefer it over
                    // whatever happens to be first in the list (dbNSFP's array order is
                    // arbitrary; taking [0] once showed BRAF V600E as c.620T>A).
                    const strongMatch = cDNAHTML.match(/<strong>([^<]*)<\/strong>/i);
                    const cdClean = (strongMatch ? strongMatch[1] : cDNAHTML.replace(/<[^>]+>/g, '').split(',')[0]).trim();
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
                    // Build with DOM nodes rather than makeLine (which escapes its value) so the
                    // link renders as a real anchor.
                    const ucscSpan = document.createElement('span');
                    const ucscLabel = document.createElement('strong');
                    ucscLabel.textContent = 'UCSC (hg19):';
                    ucscSpan.appendChild(ucscLabel);
                    ucscSpan.appendChild(document.createTextNode(' '));
                    const ucscLink = document.createElement('a');
                    ucscLink.href = ucscUrl;
                    ucscLink.target = '_blank';
                    ucscLink.rel = 'noopener noreferrer';
                    ucscLink.textContent = 'Zoom to region';
                    ucscSpan.appendChild(ucscLink);
                    content.appendChild(ucscSpan);
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
                        // dbNSFP entries have no transcript accession — render just the
                        // cDNA rather than a dangling ": c.1799T>A".
                        let inner = [t.transcript, t.cDNA].filter(Boolean).join(': ');
                        if (t.protein) inner += `, ${t.protein}`;
                        if (t.canonical) {
                            li.innerHTML = `<strong>${escapeHtml(inner)}</strong>`;
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
                // Fallback recovery from the live ClinVar region query.
                //
                // MyVariant.info sometimes carries no `clinvar` block for the queried
                // variant even though ClinVar itself has an entry (its mirror lags live
                // ClinVar, particularly for newer somatic records). When that happens
                // there is no variant_id, so significance renders as "N/A" and the
                // per-variant VCV lookup below is skipped — yet the same variant shows
                // up in the region pull. Recover the ClinVar variation ID by matching the
                // region results to the queried variant on exact position + alleles, so
                // the significance line, somatic/oncogenicity data and variation link all
                // populate as they would for a MyVariant-indexed variant.
                let prefetchedNearby = null;       // reused by the nearby-variants plot and the same-protein section
                let recoveredSignificance = '';    // germline classification recovered from the region pull
                let recoveredFromRegion = false;
                {
                    const recoveryTuple = buildVariantCoordinateTuple(rawInput, gVariant);
                    const haveNumericVariantId = variantId && /^\d+$/.test(variantId);
                    // cDNA forms used to recognise indel records, which carry no
                    // comparable ref/alt: the user's own c. notation plus every
                    // transcript's c. notation we resolved for this variant.
                    const recoveryCdnaForms = [
                        typeof targetCdnaGlobal !== 'undefined' ? targetCdnaGlobal : '',
                        ...(Array.isArray(transcriptsList) ? transcriptsList.map((t) => t && t.cDNA) : [])
                    ].filter(Boolean);
                    if (!haveNumericVariantId && sigSummary.length === 0
                        && recoveryTuple && recoveryTuple.chrom && recoveryTuple.pos) {
                        try {
                            const chrR = recoveryTuple.chrom.replace(/^chr/i, '');
                            const posR = Number(recoveryTuple.pos);
                            prefetchedNearby = await fetchClinvarRegionVariants(chrR, posR, 30);
                            const match = findExactClinvarRegionMatch(prefetchedNearby.variants, recoveryTuple, recoveryCdnaForms);
                            if (match && match.id && /^\d+$/.test(String(match.id))) {
                                variantId = String(match.id);
                                clinVarLink = `https://www.ncbi.nlm.nih.gov/clinvar/variation/${variantId}/`;
                                recoveredSignificance = match.germline || '';
                                recoveredFromRegion = true;
                            }
                        } catch (e) {
                            console.warn('ClinVar exact-match recovery failed', e);
                        }
                    }
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
                    spanVar.innerHTML = `<strong>Variation ID:</strong> ${escapeHtml(variantId)}`;
                    content.appendChild(spanVar);
                }
                // Significance summary. Kept in a reference so the per-variant VCV
                // lookup below can fill it in when the significance was only recoverable
                // from the live ClinVar record (and not from MyVariant.info).
                const spanSig = document.createElement('div');
                let haveSignificance = false;
                if (sigSummary.length > 0) {
                    spanSig.innerHTML = `<strong>Clinical significance:</strong> ${sigSummary.map(escapeHtml).join('; ')}`;
                    haveSignificance = true;
                } else if (recoveredSignificance) {
                    spanSig.innerHTML = `<strong>Clinical significance:</strong> ${escapeHtml(recoveredSignificance)}`;
                    haveSignificance = true;
                } else {
                    spanSig.innerHTML = `<strong>Clinical significance:</strong> N/A`;
                }
                content.appendChild(spanSig);
                // Note when the variant was matched via the live ClinVar region query
                // rather than the MyVariant.info annotation, for transparency.
                if (recoveredFromRegion) {
                    const recNote = document.createElement('div');
                    recNote.style.cssText = 'font-size:0.78rem;color:#6b7280;margin-top:1px;';
                    recNote.textContent = "Matched via live ClinVar query — MyVariant.info's ClinVar mirror has no record for this allele (expected for somatic-only entries).";
                    content.appendChild(recNote);
                }
                // Conditions summary (show up to 3, rest collapsed)
                if (conditionsList.length > 0) {
                    const spanCond = document.createElement('div');
                    const displayConds = conditionsList.slice(0, 3).map(escapeHtml).join(', ');
                    spanCond.innerHTML = `<strong>Conditions:</strong> ${displayConds}${conditionsList.length > 3 ? '…' : ''}`;
                    content.appendChild(spanCond);
                }
                // Somatic / oncogenicity data fetched directly from ClinVar VCV record.
                if (variantId && /^\d+$/.test(variantId)) {
                    try {
                        const cvData = await fetchClinvarVariant(variantId);
                        if (!isCurrentSearch()) return; // a newer search superseded this one
                        aiReviewExtras.clinvar_variant_record = cvData;
                        if (cvData) {
                            // If significance is still unresolved (no MyVariant rcv and the
                            // region pull carried no germline classification), fall back to
                            // the germline classification from the full VCV record.
                            if (!haveSignificance && cvData.germline && cvData.germline.description) {
                                spanSig.innerHTML = `<strong>Clinical significance:</strong> ${escapeHtml(cvData.germline.description)}`;
                                haveSignificance = true;
                            }
                            if (cvData.somatic && cvData.somatic.description) {
                                const somaticDiv = document.createElement('div');
                                somaticDiv.style.marginTop = '0.25rem';
                                somaticDiv.innerHTML = `<strong>Somatic clinical impact:</strong> ${escapeHtml(cvData.somatic.description)}`;
                                content.appendChild(somaticDiv);
                            }
                            if (cvData.oncogenicity && cvData.oncogenicity.description) {
                                const oncDiv = document.createElement('div');
                                oncDiv.innerHTML = `<strong>Oncogenicity:</strong> ${escapeHtml(cvData.oncogenicity.description)}`;
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
                                    if (sc.condition) parts.push(`<strong>${escapeHtml(sc.condition)}</strong>`);
                                    // Combine tier + assertion type + clinical significance into one readable string.
                                    const impactParts = [sc.tier, sc.assertionType, sc.clinSig].filter(Boolean);
                                    if (impactParts.length) parts.push(escapeHtml(impactParts.join(' — ')));
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
                // Same protein change in ClinVar (any underlying nucleotide change).
                //
                // Shown up here — above the nearby-region section — because same
                // amino-acid evidence is usually more relevant than nearby variants.
                // ClinVar groups by protein consequence ("CTNNB1 G34R" returns both
                // c.100G>C and c.100G>A), and the tumor/condition associations for a
                // Tier I somatic change are often split across those alternate
                // nucleotide changes, so we list each record's conditions here.
                {
                    const spAnns = annotation?.snpeff?.ann
                        ? (Array.isArray(annotation.snpeff.ann) ? annotation.snpeff.ann : [annotation.snpeff.ann])
                        : [];
                    const hgvsP = (spAnns.find((a) => a.hgvs_p) || {}).hgvs_p || '';
                    const queryPC = parseProteinChange(hgvsP)
                        || parseProteinChange(targetProtGlobal)
                        || parseProteinChange(protein);
                    const firstGene = geneNames ? geneNames.split(',')[0].trim() : '';
                    if (queryPC && firstGene) {
                        const threeRef = aaSingleToThree(queryPC.ref);
                        const threeAlt = aaSingleToThree(queryPC.alt);
                        const changeForms = [];
                        if (threeRef && threeAlt) changeForms.push(`${threeRef}${queryPC.pos}${threeAlt}`);
                        changeForms.push(`${queryPC.ref}${queryPC.pos}${queryPC.alt}`);
                        const pcKey = `p.${queryPC.ref}${queryPC.pos}${queryPC.alt}`;
                        try {
                            // Primary source: the gene+change proxy. It can fail transiently
                            // (NCBI eUtils rate limits), so degrade to the region pull, which
                            // for same-codon alternate alleles already holds the match. Reuse a
                            // prefetched region pull when present; otherwise only reach for one
                            // when the protein search returned nothing.
                            let protMatches = [];
                            try {
                                ({ variants: protMatches } = await fetchClinvarProteinVariants(firstGene, changeForms));
                            } catch (protErr) {
                                console.warn('ClinVar protein-change query failed; using region records only', protErr);
                            }
                            let fallbackVariants = (prefetchedNearby && prefetchedNearby.variants) || [];
                            if (protMatches.length === 0 && fallbackVariants.length === 0) {
                                try {
                                    const tup = buildVariantCoordinateTuple(rawInput, gVariant);
                                    if (tup && tup.chrom && tup.pos) {
                                        prefetchedNearby = await fetchClinvarRegionVariants(tup.chrom.replace(/^chr/i, ''), Number(tup.pos), 30);
                                        fallbackVariants = prefetchedNearby.variants || [];
                                    }
                                } catch (regErr) {
                                    console.warn('ClinVar region fallback for protein match failed', regErr);
                                }
                            }
                            // Keep only exact protein-change matches (the free-text search is
                            // recall-oriented); region records catch same-codon alternate
                            // alleles the protein search can under-recall. protMatches are
                            // richer (they carry conditions), so they win de-duplication.
                            const seen = new Set();
                            const sameProt = [];
                            [...protMatches, ...fallbackVariants].forEach((v) => {
                                const pc = parseProteinChange(v.title) || parseProteinChange(v.variationName);
                                if (sameProteinChange(pc, queryPC) && !seen.has(String(v.id))) {
                                    seen.add(String(v.id));
                                    sameProt.push(v);
                                }
                            });
                            if (sameProt.length > 0) {
                                // Queried allele first, then the rest.
                                sameProt.sort((a, b) => {
                                    const aq = variantId && String(a.id) === String(variantId) ? 0 : 1;
                                    const bq = variantId && String(b.id) === String(variantId) ? 0 : 1;
                                    return aq - bq;
                                });
                                const block = document.createElement('div');
                                block.style.cssText = 'margin-top:8px;';
                                const header = document.createElement('div');
                                header.style.cssText = 'font-size:0.9rem;font-weight:600;';
                                header.textContent = `Same protein change (${pcKey}) in ClinVar: ${sameProt.length}`;
                                block.appendChild(header);
                                const note = document.createElement('div');
                                note.style.cssText = 'font-size:0.75rem;color:#6b7280;margin:1px 0 4px;';
                                // The "belongs to a different nucleotide change" caveat only makes
                                // sense when the list actually contains an alternate allele. For a
                                // protein-level query (e.g. "BRAF V600E") that resolves to a single
                                // record — the queried variant itself — it would read as a false
                                // mismatch warning, so drop it in that case.
                                const queriedInList = variantId ? sameProt.filter((v) => String(v.id) === String(variantId)).length : 0;
                                const alternateNtCount = sameProt.length - queriedInList;
                                let noteText = 'Records with the same amino-acid change; conditions/tumors are listed per record.';
                                if (alternateNtCount > 0) {
                                    noteText += ' An alternate nucleotide change’s classification is that record’s own and may differ from the queried variant.';
                                }
                                note.textContent = noteText;
                                block.appendChild(note);
                                sameProt.forEach((v) => {
                                    const isQueried = variantId && String(v.id) === String(variantId);
                                    const cLabel = (String(v.title || '').match(/c\.[A-Za-z0-9>_+*\-]+/) || [])[0]
                                        || v.title || `Variation ${v.id}`;
                                    const entry = document.createElement('div');
                                    entry.style.cssText = `margin:0 0 6px;padding-left:0.55rem;border-left:3px solid ${isQueried ? '#6366f1' : '#e5e7eb'};`;
                                    const head = document.createElement('div');
                                    head.style.cssText = 'font-size:0.83rem;font-weight:600;';
                                    head.innerHTML = `<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.id}/" target="_blank" rel="noopener noreferrer">${escapeHtml(cLabel)}</a>`
                                        + (isQueried ? ' <em style="font-weight:normal">(queried variant)</em>' : '');
                                    entry.appendChild(head);
                                    const clsBits = [];
                                    if (v.germline) {
                                        const gcolor = getPathogenicityColor(v.germline, v);
                                        const stars = clinvarReviewStars(v.review);
                                        clsBits.push(`Germline: <span style="color:${gcolor};font-weight:600;display:inline">${escapeHtml(v.germline)}</span>${stars ? ' ' + stars : ''}`);
                                    }
                                    if (v.somatic) clsBits.push(`Somatic: <strong>${escapeHtml(v.somatic)}</strong>`);
                                    if (v.oncogenicity) clsBits.push(`Oncogenicity: <strong>${escapeHtml(v.oncogenicity)}</strong>`);
                                    if (clsBits.length === 0) clsBits.push('<span style="color:#6b7280;display:inline">See ClinVar</span>');
                                    // One classification per line so the narrow card wraps cleanly.
                                    clsBits.forEach((bit) => {
                                        const clsDiv = document.createElement('div');
                                        clsDiv.style.cssText = 'font-size:0.8rem;';
                                        clsDiv.innerHTML = bit;
                                        entry.appendChild(clsDiv);
                                    });
                                    if (Array.isArray(v.conditions) && v.conditions.length > 0) {
                                        const condDiv = document.createElement('div');
                                        condDiv.style.cssText = 'font-size:0.78rem;color:#4b5563;';
                                        condDiv.innerHTML = `<strong>Conditions:</strong> ${v.conditions.map(escapeHtml).join('; ')}`;
                                        entry.appendChild(condDiv);
                                    }
                                    block.appendChild(entry);
                                });
                                content.appendChild(block);
                                // Expose to the AI review context.
                                aiReviewExtras.same_protein_change_clinvar_variants = sameProt.map((v) => ({
                                    id: v.id, title: v.title, classification: v.germline || null,
                                    somatic_clinical_impact: v.somatic || null, oncogenicity: v.oncogenicity || null,
                                    review_status: v.review || null, conditions: v.conditions || [], pos: v.pos ?? null
                                }));
                            }
                        } catch (e) {
                            console.warn('ClinVar protein-change section failed', e);
                        }
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
                const tuple = buildVariantCoordinateTuple(rawInput, gVariant);
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
                        // Reuse the region pull already performed during fallback recovery
                        // (same chrom/pos/window) to avoid a duplicate eutils call.
                        const { variants: nearby, total: nearbyTotal } = prefetchedNearby
                            || await fetchClinvarRegionVariants(chr, posNum, 30);
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
                                const sigSpan = `<span style="color:${color};font-weight:600">${escapeHtml(v.germline || 'Unknown')}</span>`;
                                const posInfo = v.pos ? ` · pos ${escapeHtml(String(v.pos))}` : '';
                                const cvLink = `<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${encodeURIComponent(v.id)}/" target="_blank" rel="noopener noreferrer">${escapeHtml(String(v.id))}</a>`;
                                li.innerHTML = `${sigSpan}${posInfo} — ${v.title ? escapeHtml(v.title) : cvLink}`;
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
                    div.innerHTML = `<strong>${escapeHtml(label)}:</strong> ${escapeHtml(value)}`;
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
                            civicApiDiv.innerHTML = `<div style="font-size:0.82rem;color:#9ca3af;">Gene "${escapeHtml(civicGene)}" not found in CIViC.</div>`;
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
                                ampEl.innerHTML = `<strong>AMP/ACMG tier (CIViC):</strong> ${escapeHtml(ampLevel)}`;
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
                            varEl.innerHTML = `<strong>CIViC variant:</strong> ${escapeHtml(matchedVariant.name)}${vTypes ? ` <span style="color:#6b7280">(${escapeHtml(vTypes)})</span>` : ''} — <a href="https://civicdb.org/variants/${encodeURIComponent(matchedVariant.id)}/summary" target="_blank" rel="noopener noreferrer">View ↗</a>`;
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
                // Link to gnomAD v2. Rendered synchronously from whatever the notation
                // gives us (variant page, or a region view for an indel), then upgraded
                // to the exact variant page once the reference-derived alleles land.
                const gnomadLinkEl = document.createElement('a');
                gnomadLinkEl.target = '_blank';
                gnomadLinkEl.rel = 'noopener noreferrer';
                gnomadLinkEl.textContent = 'View on gnomAD (v2.1.1)';
                gnomadLinkEl.href = gnomadLink || '';
                if (!gnomadLink) gnomadLinkEl.style.display = 'none';
                content.appendChild(gnomadLinkEl);
                getResolvedVcfAlleles().then((alleles) => {
                    if (!alleles) return;
                    const upgraded = buildGnomadVariantUrl(rawInput, gVariant, annotation && annotation._id, annotation && annotation.vcf, alleles);
                    if (!upgraded) return;
                    gnomadLinkEl.href = upgraded;
                    gnomadLinkEl.style.display = '';
                });
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
                (async () => {
                    const showV4Msg = (msg) => {
                        v4Section.innerHTML = `<div style="font-size:0.82rem;color:#9ca3af;">${escapeHtml(msg)}</div>`;
                    };
                    // Indels arrive as "g.START_ENDdel"-style notation with no alleles;
                    // resolveVcfAlleles() reads them off the reference and left-aligns
                    // them, which is the representation gnomAD indexes by.
                    const alleles = await getResolvedVcfAlleles();
                    const chromCoord = alleles && alleles.chrom;
                    const pos37Coord = alleles && alleles.pos;
                    const refCoord = alleles && alleles.ref;
                    const altCoord = alleles && alleles.alt;
                    if (!chromCoord || !pos37Coord || !refCoord || !altCoord) {
                        // A codon-only resolution has a position but no alleles by
                        // definition; say so rather than implying a lookup problem.
                        v4Loading.textContent = codonOnlyResolutionGlobal
                            ? 'gnomAD v4.1: needs a specific allele — protein notation resolves only to the codon.'
                            : 'gnomAD v4.1: insufficient coordinate data.';
                        return;
                    }
                    // The annotation usually already carries the GRCh38 start; handing
                    // it to the proxy avoids a liftover round-trip that can fail.
                    const knownPos38 = extractHg38Start(annotation);
                    return fetchGnomadV4(chromCoord, pos37Coord, refCoord, altCoord, knownPos38).then((result) => {
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
                        const grch37Id = `${chromCoord}-${pos37Coord}-${refCoord}-${altCoord}`;
                        if (status === 'liftover_unavailable') {
                            // Ensembl is down, not the variant's fault — say which,
                            // so "try again later" reads as real advice.
                            showV4Msg(`gnomAD v4.1 unavailable: the Ensembl liftover service is not responding, so ${chromCoord}:${pos37Coord} could not be converted to GRCh38. This is usually temporary — try again shortly.`);
                            return;
                        }
                        if (status === 'liftover_failed') {
                            showV4Msg(`GRCh37→GRCh38 liftover failed: ${chromCoord}:${pos37Coord} has no unambiguous GRCh38 equivalent.`);
                            return;
                        }
                        if (status === 'not_found') {
                            showV4Msg(`Not found in gnomAD v4.1 (queried: ${grch38Id || grch37Id}${grch38Id ? `, from GRCh37 ${grch37Id}` : ''}).`);
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
                })().catch(() => {
                    v4Section.innerHTML = '<div style="font-size:0.82rem;color:#9ca3af;">gnomAD v4.1 unavailable.</div>';
                });
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
                    summaryParts.push(`<strong>${escapeHtml(name)}</strong>: ${emojis.join('/')}`);
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
                        li.innerHTML = `<strong>${escapeHtml(n)}</strong>: ${escapeHtml(v)}`;
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
                if (gVariant && !codonOnlyResolutionGlobal) {
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
                sigSpan.innerHTML = `<strong>Oncogenicity:</strong> ${escapeHtml(oncogenic || 'N/A')}`;
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
                        span.innerHTML = `<strong>Found in:</strong> ${escapeHtml(items['Total Tumors'])} tumor${items['Total Tumors'] === 1 ? '' : 's'}`;
                        content.appendChild(span);
                    }
                    // Frequency overall
                    if (items['Frequency (overall)']) {
                        const p = document.createElement('p');
                        p.innerHTML = `<strong>Frequency (overall):</strong> ${escapeHtml(items['Frequency (overall)'])}`;
                        content.appendChild(p);
                    }
                    // Find the key that starts with "Frequency in" (gene-specific frequency)
                    Object.keys(items).forEach(key => {
                        if (key.startsWith('Frequency in')) {
                            const p = document.createElement('p');
                            p.innerHTML = `<strong>${escapeHtml(key)}:</strong> ${escapeHtml(items[key])}`;
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
                        span.innerHTML = `<strong>Mutation Frequency:</strong> ${escapeHtml(items['Mutation Frequency'])}`;
                        content.appendChild(span);
                    }
                    // Add any other COSMIC base items except those with html
                    Object.entries(items).forEach(([k,v]) => {
                        if (k === 'Mutation Frequency') return;
                        if (v && typeof v === 'object' && v.html) return;
                        const p = document.createElement('p');
                        p.innerHTML = `<strong>${escapeHtml(k)}:</strong> ${escapeHtml(v)}`;
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
            // Card: Cancer Prevalence (cBioPortal cohort frequency)
            {
                const genesForCbio = geneNames ? geneNames.split(',').map((g) => g.trim()).filter(Boolean) : [];
                let cbioGene = genesForCbio.find((g) => !isChromosomeLikeGeneSymbol(g)) || '';
                if (!cbioGene && geneHintGlobal && !isChromosomeLikeGeneSymbol(geneHintGlobal)) cbioGene = geneHintGlobal;
                createCbioportalCard({
                    container: cardsContainer,
                    gene: cbioGene,
                    proteinChange: protein || targetProtGlobal || '',
                    extras: aiReviewExtras
                });
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
                variantSummary.innerHTML = `<strong>Variant:</strong> ${escapeHtml(tp53Protein || tp53Cdna || tp53Genomic || 'N/A')}`;
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
                    if (!isCurrentSearch()) return; // a newer search superseded this one
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

            const getAiReviewVariantCoordinates = async () => {
                const tuple = buildVariantCoordinateTuple(rawInput, gVariant);
                const vcfData = annotation && annotation.vcf;
                if (tuple) {
                    // Prefer reference-resolved alleles so indels carry a usable
                    // ref/alt into the gnomAD/SpliceAI parts of the AI payload.
                    const alleles = await getResolvedVcfAlleles();
                    if (alleles && alleles.ref && alleles.alt) {
                        return { chrom: alleles.chrom, pos37: alleles.pos, ref: alleles.ref, alt: alleles.alt };
                    }
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
                const coords = await getAiReviewVariantCoordinates();
                const clinvarVariantId = annotation?.clinvar?.variant_id && /^\d+$/.test(String(annotation.clinvar.variant_id))
                    ? String(annotation.clinvar.variant_id)
                    : '';
                const pubmedTerm = [aiReviewGene, aiReviewSearchVariantTerm].filter(Boolean).join(' ');
                const pubmedTumorTerm = tumorType ? [aiReviewGene, tumorType].filter(Boolean).join(' ') : '';
                const pubmedVariantTumorTerm = (tumorType && aiReviewSearchVariantTerm)
                    ? [aiReviewGene, aiReviewSearchVariantTerm, tumorType].filter(Boolean).join(' ')
                    : '';
                // Prefer card-cached results when available. The cards fetch the same
                // data on render; re-fetching everything on every "Run AI review" click
                // re-hammers rate-limited upstreams (NCBI eutils above all — when efetch
                // trips the limit the proxy returns articles with EMPTY abstracts, so the
                // AI payload used to miss abstracts the cards already had) and adds
                // seconds of latency per run. `cachedOr` reuses what a card stored in
                // aiReviewExtras and only fetches when the card's data is absent/errored.
                const cachedOr = (key, isValid, fetcher) => {
                    const cached = aiReviewExtras[key];
                    if (cached != null && isValid(cached)) return Promise.resolve(cached);
                    return fetcher();
                };
                const pubmedFromCacheOrFetch = (key, term, variantQuery = null) => cachedOr(
                    key,
                    (c) => Array.isArray(c.articles) && c.articles.length > 0,
                    () => (term || variantQuery ? fetchPubmedArticles(term, 5, variantQuery) : Promise.resolve({ total: 0, articles: [] }))
                );
                const aiPubmedVariantQuery = aiReviewSearchVariantTerm
                    ? { gene: aiReviewGene, variant: aiReviewSearchVariantTerm }
                    : null;
                const spliceApiVariant = buildSpliceAiApiVariant(rawInput, gVariant, annotation, await getResolvedVcfAlleles());
                const supplemental = { ...aiReviewExtras };
                const tasks = [
                    // The ClinVar card may have recovered a variation ID from the region
                    // pull that the MyVariant annotation lacks — reusing its record (not
                    // just saving a fetch) also stops that recovery from being clobbered
                    // with null here.
                    ['clinvar_variant_record', cachedOr('clinvar_variant_record',
                        (c) => typeof c === 'object' && !c.error,
                        () => (clinvarVariantId ? fetchClinvarVariant(clinvarVariantId) : Promise.resolve(null)))],
                    ['nearby_clinvar_variants', coords.chrom && coords.pos37 ? fetchClinvarRegionVariants(coords.chrom, coords.pos37, 5).then(r => ({
                        note: 'ClinVar variants within ±5bp of the queried position. The queried variant itself may appear here if it has its own ClinVar entry. Do not use these neighboring variants\' classifications as the classification for the queried variant.',
                        window_bp: 5,
                        variants: r.variants
                    })) : Promise.resolve({ note: 'ClinVar variants within ±5bp of the queried position.', window_bp: 5, variants: [] })],
                    ['civic_api', cachedOr('civic_api',
                        (c) => typeof c === 'object' && !c.error,
                        () => (aiReviewGene ? fetchCivicApiData(aiReviewGene, aiReviewProtein) : Promise.resolve(null)))],
                    ['gnomad_v4', cachedOr('gnomad_v4',
                        (c) => typeof c === 'object' && Boolean(c.status),
                        () => (coords.chrom && coords.pos37 && coords.ref && coords.alt ? fetchGnomadV4(coords.chrom, coords.pos37, coords.ref, coords.alt, extractHg38Start(annotation)) : Promise.resolve(null)))],
                    ['spliceai_lookup', cachedOr('spliceai_lookup',
                        (c) => typeof c === 'object' && !c.error && Boolean(c.data),
                        () => (spliceApiVariant ? fetchSpliceAiPrediction(spliceApiVariant, { hg: '37', distance: 500, mask: 0, bc: 'basic' }) : Promise.resolve(null)))],
                    ['pubmed', pubmedTerm ? pubmedFromCacheOrFetch('pubmed', pubmedTerm, aiPubmedVariantQuery) : Promise.resolve({ total: 0, articles: [] })],
                    ['pubmed_tumor_type', pubmedTumorTerm && pubmedTumorTerm !== pubmedTerm ? pubmedFromCacheOrFetch('pubmed_tumor_type', pubmedTumorTerm) : Promise.resolve(null)],
                    ['pubmed_variant_tumor_type', pubmedVariantTumorTerm
                        && pubmedVariantTumorTerm !== pubmedTerm
                        && pubmedVariantTumorTerm !== pubmedTumorTerm
                        ? pubmedFromCacheOrFetch('pubmed_variant_tumor_type', pubmedVariantTumorTerm,
                            aiPubmedVariantQuery ? { ...aiPubmedVariantQuery, extra: tumorType } : null)
                        : Promise.resolve(null)],
                    // The openFDA card stores the same payload shape under `openfda`;
                    // reuse it rather than re-paging up to 10 requests of label text.
                    ['openfda_drug_labels', cachedOr('openfda',
                        (c) => Array.isArray(c.results),
                        () => (aiReviewGene ? fetchOpenFdaDrugLabels(aiReviewGene).catch(() => null) : Promise.resolve(null)))],
                    ['bbkb_biomarker_therapies', cachedOr('bbkb_biomarker_therapies',
                        (c) => Array.isArray(c.results),
                        () => (aiReviewGene ? fetchBbkbBiomarkerTherapies(aiReviewGene).catch(() => ({ total_matched: 0, returned: 0, results: [] })) : Promise.resolve(null)))],
                    ['tp53_mutation_database', cachedOr('tp53_mutation_database',
                        (c) => typeof c === 'object' && !c.error,
                        () => (isTp53Gene(geneNames) ? fetchTp53MutationDatabase({
                            gene: 'TP53',
                            protein: normaliseAiReviewVariantText(aiReviewProtein),
                            cdna: normaliseAiReviewVariantText(aiReviewCdna),
                            genomic: String(gVariant || '').trim(),
                            debug: true
                        }) : Promise.resolve(null)))]
                ];
                const settled = await Promise.allSettled(tasks.map(([, promise]) => promise));
                settled.forEach((result, idx) => {
                    const key = tasks[idx][0];
                    supplemental[key] = result.status === 'fulfilled'
                        ? result.value
                        : { error: result.reason?.message || String(result.reason || 'Unavailable') };
                });
                // The openFDA card's copy (spread in as `openfda`) would duplicate
                // `openfda_drug_labels` byte-for-byte — for drug-rich genes that is the
                // dominant part of the payload, doubled. Keep the one canonical key.
                if (supplemental.openfda_drug_labels && supplemental.openfda) delete supplemental.openfda;
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
                const spliceTuple = buildVariantCoordinateTuple(rawInput, gVariant);
                // Recover REF from annotation.vcf for MNVs where delins notation drops the reference.
                const spliceRef = (spliceTuple && spliceTuple.ref) || (annotation && annotation.vcf && String(annotation.vcf.ref || '').toUpperCase()) || null;
                const spliceAlt = (spliceTuple && spliceTuple.alt) || (annotation && annotation.vcf && String(annotation.vcf.alt || '').toUpperCase()) || null;
                const spliceVariantText = (spliceTuple && spliceRef && spliceAlt) ? `${spliceTuple.chrom} ${spliceTuple.pos} ${spliceRef} ${spliceAlt}` : '';
                // SpliceAI lookup defaults to hg38 when hg is omitted. Most MyVariant coordinates
                // we surface in this app are hg19/GRCh37, so explicitly request hg=37.
                const buildSpliceAiPageUrl = (variantText) => (variantText
                    ? `https://spliceailookup.broadinstitute.org/#variant=${encodeURIComponent(variantText)}&hg=37`
                    : 'https://spliceailookup.broadinstitute.org/#hg=37');
                const spliceAiUrl = buildSpliceAiPageUrl(spliceVariantText);

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
                let spliceApiVariant = buildSpliceAiApiVariant(rawInput, gVariant, annotation);
                const spliceHint = document.createElement('span');
                spliceHint.style.fontSize = '0.85rem';
                spliceHint.style.color = '#4a5f73';
                spliceHint.textContent = spliceVariantText
                    ? `SpliceAI query: ${spliceVariantText}`
                    : 'Resolving genomic alleles for the SpliceAI query…';
                spliceContent.appendChild(spliceHint);
                const spliceResultsDiv = document.createElement('div');
                spliceResultsDiv.style.cssText = 'margin-top:0.45rem;';
                const showSpliceLoading = () => {
                    spliceResultsDiv.innerHTML = '';
                    const loading = document.createElement('div');
                    loading.style.cssText = 'font-size:0.82rem;color:#6b7280;font-style:italic;';
                    loading.textContent = 'Loading SpliceAI scores…';
                    spliceResultsDiv.appendChild(loading);
                };
                if (spliceApiVariant) showSpliceLoading();
                spliceContent.appendChild(spliceResultsDiv);
                spliceCard.appendChild(spliceContent);
                cardsContainer.appendChild(spliceCard);

                // Indels reach here with no ref/alt in their notation. Resolve them off
                // the reference so the SpliceAI query — and the deep link — work for
                // anything more complex than a substitution.
                const spliceReady = getResolvedVcfAlleles().then((alleles) => {
                    if (spliceApiVariant || !alleles) {
                        if (!spliceApiVariant && !spliceVariantText) {
                            spliceHint.textContent = 'No chr/pos/ref/alt tuple could be resolved; opening SpliceAI home page.';
                        }
                        return;
                    }
                    const resolvedText = `chr${alleles.chrom} ${alleles.pos} ${alleles.ref} ${alleles.alt}`;
                    spliceApiVariant = buildSpliceAiApiVariant(rawInput, gVariant, annotation, alleles);
                    spliceHint.textContent = `SpliceAI query: ${resolvedText}`;
                    spliceLinkLine.innerHTML = `<a href="${buildSpliceAiPageUrl(resolvedText)}" target="_blank" rel="noopener noreferrer">Open SpliceAI lookup 🔍</a>`;
                    if (spliceApiVariant) showSpliceLoading();
                });

                spliceReady.then(() => {
                    if (!spliceApiVariant) return;
                    return fetchSpliceAiPrediction(spliceApiVariant, { hg: '37', distance: 500, mask: 0, bc: 'basic' }).then((spliceData) => {
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
                            bestEl.innerHTML = `<strong>Max delta score:</strong> ${summary.best.value.toFixed(3)} (${pct}%) ${escapeHtml(summary.best.label)}${summary.best.position !== null ? ` at ${escapeHtml(String(summary.best.position))}` : ''}${summary.best.transcript ? ` · ${escapeHtml(summary.best.transcript)}` : ''}`;
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
                });

                // Card: PubMed
                {
                    const pmSearchTerm = [firstGene, searchVariantTerm].filter(Boolean).join(' ');
                    const pmTumorSearchTerm = tumorType ? [firstGene, tumorType].filter(Boolean).join(' ') : '';
                    const pmVariantTumorSearchTerm = (tumorType && searchVariantTerm)
                        ? [firstGene, searchVariantTerm, tumorType].filter(Boolean).join(' ')
                        : '';
                    // Variant tabs go through LitVar2 (all nomenclature spellings)
                    // with the term search as automatic fallback; the tumor-only
                    // tab stays a plain term search.
                    const pmVariantQuery = searchVariantTerm
                        ? { gene: firstGene, variant: searchVariantTerm }
                        : null;
                    const pmTabs = [{ label: 'Gene + Variant', term: pmSearchTerm, extraKey: 'pubmed', variantQuery: pmVariantQuery }];
                    if (pmTumorSearchTerm && pmSearchTerm !== pmTumorSearchTerm) {
                        pmTabs.push({ label: 'Gene + Tumor Type', term: pmTumorSearchTerm, extraKey: 'pubmed_tumor_type' });
                        if (pmVariantTumorSearchTerm
                            && pmVariantTumorSearchTerm !== pmSearchTerm
                            && pmVariantTumorSearchTerm !== pmTumorSearchTerm) {
                            pmTabs.push({
                                label: 'Gene + Variant + Tumor Type',
                                term: pmVariantTumorSearchTerm,
                                extraKey: 'pubmed_variant_tumor_type',
                                variantQuery: pmVariantQuery ? { ...pmVariantQuery, extra: tumorType } : null
                            });
                        }
                    }
                    createPubmedCard({ container: cardsContainer, tabs: pmTabs, extras: aiReviewExtras });
                }

                // Card: FDA-Approved Drugs (by gene)
                createFdaDrugsCard({ container: cardsContainer, gene: firstGene, extras: aiReviewExtras, cardCache: cardDataCache });

                // Card: Clinical Trials
                createClinicalTrialsCard({ container: cardsContainer, gene: firstGene, tumorType, cardCache: cardDataCache });
            }
            // Card: Guidelines
            createGuidelinesCard({
                container: cardsContainer,
                tumorType,
                getGene: () => aiReviewGene || (geneNames ? geneNames.split(',')[0].trim() : ''),
                extras: aiReviewExtras
            });

            // Optional AI review card (manual trigger; sends current annotation context to OpenRouter via backend proxy).
            createAiReviewCard({
                container: cardsContainer,
                idSuffix: '',
                introText: 'Send the retrieved variant data to OpenRouter for a structured draft interpretation. No request is sent until you click Run AI review.',
                loadingText: 'Gathering FDA, trial, and annotation context for AI review…',
                buildContext: async (userNotes) => {
                    // Reuse what the cards already fetched; fall back to a fresh fetch
                    // only when a card's data never arrived (or errored).
                    const fdaPromise = Array.isArray(cardDataCache.fda_companion_diagnostics_records)
                        ? Promise.resolve(cardDataCache.fda_companion_diagnostics_records)
                        : (aiReviewGene ? fetchFdaCompanionDiagnostics(aiReviewGene).catch(() => []) : Promise.resolve([]));
                    const bbkbPromise = (aiReviewExtras.bbkb_biomarker_therapies && Array.isArray(aiReviewExtras.bbkb_biomarker_therapies.results))
                        ? Promise.resolve(aiReviewExtras.bbkb_biomarker_therapies)
                        : (aiReviewGene ? fetchBbkbBiomarkerTherapies(aiReviewGene).catch(() => ({ total_matched: 0, returned: 0, results: [] })) : Promise.resolve(null));
                    const trialsPromise = (cardDataCache.clinical_trials && Array.isArray(cardDataCache.clinical_trials.studies))
                        ? Promise.resolve(cardDataCache.clinical_trials)
                        : (aiReviewGene ? fetchClinicalTrials(aiReviewGene, tumorType).catch(() => ({ total: 0, studies: [] })) : Promise.resolve({ total: 0, studies: [] }));
                    const [fdaRecords, bbkbTherapies, clinicalTrialData, supplementalContext] = await Promise.all([
                        fdaPromise,
                        bbkbPromise,
                        trialsPromise,
                        fetchAiReviewSupplementalContext().catch((err) => ({ error: err.message || 'Supplemental context unavailable' }))
                    ]);
                    // Cap the openFDA record count (the dominant payload for drug-rich genes)
                    // on the shallow supplemental copy, keeping each record's full text and
                    // leaving the on-screen card data intact. supplementalContext is fresh per run.
                    if (supplementalContext && typeof supplementalContext === 'object') {
                        if (supplementalContext.openfda_drug_labels) supplementalContext.openfda_drug_labels = condenseOpenFdaForAi(supplementalContext.openfda_drug_labels);
                        if (supplementalContext.openfda) supplementalContext.openfda = condenseOpenFdaForAi(supplementalContext.openfda);
                    }
                    return {
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
                        bbkb_biomarker_therapies: bbkbTherapies,
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
                }
            });

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
            console.error(err);
            // A stale run's failure must not overwrite the status/progress panel
            // of a newer search that is already underway.
            if (!isCurrentSearch()) return;
            statusEl.textContent = 'Error: ' + err.message;
            LookupProgress.finish('fail', err.upstreamOutage ? 'Upstream service unavailable' : 'Lookup failed');
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
