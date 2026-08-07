// Shared GRCh37 ↔ GRCh38 coordinate liftover, backed by the Ensembl REST
// assembly-mapping API.
//
// The `_` prefix keeps Vercel from treating this file as a routable serverless
// function — it is imported by the sibling handlers only.
//
// Why this exists: the gnomAD v4 card needs a GRCh38 coordinate, and the only
// liftover the app had was a single un-retried call to rest.ensembl.org/map.
// That host is a single point of failure and it does fall over — during one
// outage every /map and /vep request answered HTTP 500 (or simply hung) for
// hours while /info/ping stayed green, so the app reported
// "GRCh37→GRCh38 liftover failed for 4:143094904" on a coordinate that maps
// perfectly well, to 4:142173751. The message blamed the variant for an Ensembl
// outage, and no amount of retrying by the user would help.
//
// Three things make that failure mode much rarer, and honest when it happens:
//
//   1. Fail over to the GRCh37 mirror (grch37.rest.ensembl.org), a separate
//      deployment that answers the same /map paths in both directions and
//      stayed healthy throughout the outage above.
//   2. Hedge rather than wait. A sick host does not politely return an error —
//      it hangs, and a strictly sequential "try A, then try B" spends the whole
//      budget waiting on A. So the mirror is started once the primary has been
//      quiet for HEDGE_DELAY_MS (or immediately, if the primary already failed),
//      and the first usable answer wins. A healthy primary answers well inside
//      the hedge delay, so the mirror is normally never called.
//   3. Report "the provider is unavailable" separately from "this coordinate
//      genuinely has no equivalent in the target assembly". Those are different
//      answers and the UI should not present the first as the second.
//
// Time budget: api/gnomad-v4.js runs with maxDuration 15 and still has a gnomAD
// GraphQL query to make afterwards, so the whole ladder works to an overall
// deadline (default 6s) rather than letting attempts accumulate. That is
// deliberately below the 8s single-shot timeout this replaced.

const MAPPING_HOSTS = [
    'https://rest.ensembl.org',
    'https://grch37.rest.ensembl.org',
];

const RETRYABLE_STATUS = new Set([429, 500, 502, 503, 504]);
const MAX_ROUNDS = 2;
const DEFAULT_ATTEMPT_TIMEOUT_MS = 3000;
const DEFAULT_DEADLINE_MS = 6000;
const DEFAULT_HEDGE_DELAY_MS = 1200;

function sleep(ms) {
    return new Promise((resolve) => setTimeout(resolve, ms));
}

// Normalise "chr7"/"7"/7 to the bare name Ensembl expects.
export function normalizeChrom(chrom) {
    return String(chrom == null ? '' : chrom).trim().replace(/^chr/i, '');
}

// Pick the mapping that actually answers the question we asked. Ensembl returns
// an array and the old code took [0] blind; a segment can map to a patch or
// scaffold or to another chromosome, and taking that silently produces a
// plausible-looking but wrong coordinate. Require the target assembly and the
// same chromosome, and treat several disagreeing candidates as no answer rather
// than guessing between them.
function selectMapping(mappings, chrom, targetAssembly) {
    if (!Array.isArray(mappings)) return null;
    const wanted = mappings.filter((m) => {
        const mapped = m && m.mapped;
        if (!mapped) return false;
        if (String(mapped.assembly || '').toUpperCase() !== targetAssembly.toUpperCase()) return false;
        if (normalizeChrom(mapped.seq_region_name) !== chrom) return false;
        return Number.isInteger(mapped.start) && mapped.start > 0;
    });
    if (wanted.length === 0) return null;
    const starts = new Set(wanted.map((m) => m.mapped.start));
    if (starts.size > 1) return null; // ambiguous — refuse to pick a winner
    return wanted[0].mapped;
}

// One request to one host.
// Resolves { ok: true, pos, host }
//        | { ok: false, definitive: true, ... }  — host answered; the answer is no
//        | { ok: false, definitive: false, ... } — host failed; someone else may know
async function queryHost(host, path, chrom, targetAssembly, timeoutMs) {
    let res;
    try {
        res = await fetch(`${host}${path}`, {
            headers: { Accept: 'application/json' },
            signal: AbortSignal.timeout(timeoutMs),
        });
    } catch (err) {
        // Network error or timeout — always worth another host/round.
        return { ok: false, definitive: false, detail: `${host}: ${err && err.message ? err.message : String(err)}` };
    }

    if (!res.ok) {
        return {
            ok: false,
            definitive: false,
            detail: `${host}: HTTP ${res.status}`,
            retryAfter: RETRYABLE_STATUS.has(res.status) ? Number(res.headers.get('retry-after')) : NaN,
        };
    }

    let data;
    try {
        data = await res.json();
    } catch {
        // A 200 carrying non-JSON means the host is serving an error page, not
        // that the coordinate is unmappable.
        return { ok: false, definitive: false, detail: `${host}: malformed JSON response` };
    }

    const mapped = selectMapping(data && data.mappings, chrom, targetAssembly);
    if (mapped) return { ok: true, pos: mapped.start, host };

    // The provider answered cleanly and had nothing for us. That is a real
    // "unmappable" — stop here rather than shopping it around.
    return { ok: false, definitive: true, detail: `${host}: no unambiguous ${targetAssembly} mapping` };
}

// Run one hedged round across the hosts: start the first, start the next once
// the previous has failed or gone quiet for hedgeDelayMs, and resolve with the
// first usable answer. Falls back to the last failure if every host fails.
function hedgedRound(hosts, run, hedgeDelayMs) {
    return new Promise((resolve) => {
        const started = new Set();
        const timers = [];
        let completed = 0;
        let settled = false;
        let lastFailure = null;

        const finish = (result) => {
            if (settled) return;
            settled = true;
            timers.forEach(clearTimeout);
            resolve(result);
        };

        const launch = (i) => {
            if (settled || i >= hosts.length || started.has(i)) return;
            started.add(i);

            run(hosts[i]).then((result) => {
                completed++;
                if (settled) return;
                if (result.ok || result.definitive) { finish(result); return; }
                lastFailure = result;
                launch(i + 1); // failed fast — bring the next host forward
                if (completed === started.size && started.size === hosts.length) finish(result);
            });

            if (i + 1 < hosts.length) {
                timers.push(setTimeout(() => launch(i + 1), hedgeDelayMs));
            }
        };

        launch(0);
        if (hosts.length === 0) finish({ ok: false, definitive: false, detail: 'no hosts configured' });
        else if (settled === false && started.size === 0) finish(lastFailure || { ok: false, definitive: false, detail: 'no hosts tried' });
    });
}

// Map a single coordinate between assemblies.
//
// Returns one of:
//   { ok: true,  pos, host }               — mapped; `pos` is the target-assembly start
//   { ok: false, reason: 'unmapped' }      — provider answered; the coordinate has no
//                                            unambiguous equivalent in the target
//   { ok: false, reason: 'unavailable', detail } — every host/attempt failed
//
// Never throws: callers turn the outcome into a user-facing status.
export async function liftoverPosition(chrom, pos, {
    from = 'GRCh37',
    to = 'GRCh38',
    attemptTimeoutMs = DEFAULT_ATTEMPT_TIMEOUT_MS,
    deadlineMs = DEFAULT_DEADLINE_MS,
    hedgeDelayMs = DEFAULT_HEDGE_DELAY_MS,
    hosts = MAPPING_HOSTS,
    now = () => Date.now(),
} = {}) {
    const c = normalizeChrom(chrom);
    const p = Number.parseInt(pos, 10);
    if (!c || !Number.isInteger(p) || p <= 0) {
        return { ok: false, reason: 'unmapped', detail: `Invalid coordinate ${chrom}:${pos}` };
    }

    const path = `/map/human/${from}/${c}:${p}..${p}:1/${to}?content-type=application/json`;
    const deadline = now() + deadlineMs;
    let lastDetail = null;

    for (let round = 0; round < MAX_ROUNDS; round++) {
        const remaining = deadline - now();
        if (remaining <= 0) break;

        const result = await hedgedRound(
            hosts,
            (host) => queryHost(host, path, c, to, Math.min(attemptTimeoutMs, remaining)),
            Math.min(hedgeDelayMs, remaining),
        );

        if (result.ok) return { ok: true, pos: result.pos, host: result.host };
        if (result.definitive) {
            return {
                ok: false,
                reason: 'unmapped',
                detail: `${from} ${c}:${p} has no unambiguous ${to} equivalent`,
            };
        }

        lastDetail = result.detail || lastDetail;
        if (round === MAX_ROUNDS - 1) break;

        // Exponential backoff with jitter, never below a Retry-After hint.
        let waitMs = 250 * 2 ** round + Math.floor(Math.random() * 150);
        if (Number.isFinite(result.retryAfter) && result.retryAfter > 0) {
            waitMs = Math.max(waitMs, Math.min(result.retryAfter * 1000, 2000));
        }
        if (now() + waitMs >= deadline) break;
        await sleep(waitMs);
    }

    return {
        ok: false,
        reason: 'unavailable',
        detail: lastDetail || 'Ensembl assembly-mapping API unreachable',
    };
}
