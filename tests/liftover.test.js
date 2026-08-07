/*
 * Tests for the shared assembly-liftover helper (api/_liftover.js) and the
 * gnomAD v4 proxy's use of it.
 *
 * Context: the gnomAD v4 card reported
 *   "GRCh37→GRCh38 liftover failed for 4:143094904"
 * for a coordinate that maps fine (to 4:142173751). The cause was not the
 * variant — rest.ensembl.org was answering HTTP 500 to every /map request — but
 * the proxy made one un-retried call to one host and collapsed every failure
 * into a single "liftover_failed" status, so an outage was indistinguishable
 * from an unmappable coordinate.
 *
 * These tests pin the three behaviours that fix it: retry/failover, the
 * outage-vs-unmappable distinction, and the caller-supplied pos38 fast path.
 *
 * global.fetch is stubbed so no network is touched. Run with:
 *   node tests/liftover.test.js
 */

import { liftoverPosition, normalizeChrom } from '../api/_liftover.js';
import handler from '../api/gnomad-v4.js';

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) { passed++; }
    else { failed++; console.error(`  ✗ ${name}`); }
}

function mapResponse(chrom, start, assembly = 'GRCh38') {
    return {
        ok: true,
        status: 200,
        headers: { get: () => null },
        json: async () => ({
            mappings: [{
                original: { assembly: 'GRCh37', seq_region_name: chrom, start: 1, end: 1, strand: 1 },
                mapped: { assembly, seq_region_name: chrom, start, end: start, strand: 1 },
            }],
        }),
        text: async () => '',
    };
}

function errorResponse(status, headers = {}) {
    const h = new Map(Object.entries(headers).map(([k, v]) => [k.toLowerCase(), String(v)]));
    return {
        ok: false,
        status,
        headers: { get: (k) => (h.has(k.toLowerCase()) ? h.get(k.toLowerCase()) : null) },
        json: async () => { throw new Error('not json'); },
        text: async () => 'error page',
    };
}

// Speed the tests up: neutralise the backoff sleeps.
const realSetTimeout = global.setTimeout;
global.setTimeout = (fn) => realSetTimeout(fn, 0);

// A monotonic fake clock so deadline logic is deterministic.
function fakeClock() {
    let t = 0;
    return { now: () => t, advance: (ms) => { t += ms; } };
}

async function run() {
    // 1. Happy path: first host answers, coordinate comes back mapped.
    {
        const calls = [];
        global.fetch = async (url) => { calls.push(url); return mapResponse('4', 142173751); };
        const out = await liftoverPosition('chr4', 143094904);
        check('maps 4:143094904 (GRCh37) to 142173751 (GRCh38)', out.ok && out.pos === 142173751);
        check('strips the chr prefix before calling Ensembl', calls[0].includes('/map/human/GRCh37/4:143094904..143094904:1/GRCh38'));
        check('needs only one request when the first host answers', calls.length === 1);
    }

    // 2. The regression: the primary host 500s on every request (the real outage).
    //    The GRCh37 mirror answers, so the lookup must still succeed.
    {
        const calls = [];
        global.fetch = async (url) => {
            calls.push(url);
            if (url.startsWith('https://rest.ensembl.org')) return errorResponse(500);
            return mapResponse('4', 142173751);
        };
        const out = await liftoverPosition('4', 143094904);
        check('fails over to the GRCh37 mirror when the primary host 500s', out.ok && out.pos === 142173751);
        check('failover cost one extra request, not a full retry ladder', calls.length === 2);
        check('reports which host answered', out.host === 'https://grch37.rest.ensembl.org');
    }

    // 3. Every host down => "unavailable", explicitly NOT "unmapped".
    {
        let n = 0;
        global.fetch = async () => { n++; return errorResponse(503); };
        const out = await liftoverPosition('4', 143094904, { now: fakeClock().now });
        check('all hosts failing yields reason "unavailable"', !out.ok && out.reason === 'unavailable');
        check('an outage is never reported as an unmappable coordinate', out.reason !== 'unmapped');
        check('retries across both hosts before giving up', n === 4);
        check('carries a detail explaining what failed', /HTTP 503/.test(out.detail || ''));
    }

    // 4. Transient network errors are retried too, and a later attempt can win.
    {
        let n = 0;
        global.fetch = async () => {
            n++;
            if (n < 3) throw new Error('network timeout');
            return mapResponse('7', 140753336);
        };
        const out = await liftoverPosition('7', 140453136);
        check('recovers after transient network errors', out.ok && out.pos === 140753336);
    }

    // 5. A clean 200 with no usable mapping is a genuine "unmapped" — and must
    //    stop immediately rather than burning retries on a definitive answer.
    {
        let n = 0;
        global.fetch = async () => {
            n++;
            return { ok: true, status: 200, headers: { get: () => null }, json: async () => ({ mappings: [] }), text: async () => '' };
        };
        const out = await liftoverPosition('4', 143094904);
        check('empty mappings yields reason "unmapped"', !out.ok && out.reason === 'unmapped');
        check('does not retry a definitive negative answer', n === 1);
    }

    // 6. Mappings that answer a different question are rejected: wrong assembly,
    //    wrong chromosome, or two candidates that disagree.
    {
        global.fetch = async () => mapResponse('4', 142173751, 'GRCh37');
        const wrongAssembly = await liftoverPosition('4', 143094904);
        check('rejects a mapping in the wrong assembly', !wrongAssembly.ok && wrongAssembly.reason === 'unmapped');

        global.fetch = async () => mapResponse('11', 142173751);
        const wrongChrom = await liftoverPosition('4', 143094904);
        check('rejects a mapping that landed on another chromosome', !wrongChrom.ok && wrongChrom.reason === 'unmapped');

        global.fetch = async () => ({
            ok: true,
            status: 200,
            headers: { get: () => null },
            json: async () => ({
                mappings: [
                    { mapped: { assembly: 'GRCh38', seq_region_name: '4', start: 142173751, end: 142173751, strand: 1 } },
                    { mapped: { assembly: 'GRCh38', seq_region_name: '4', start: 999999, end: 999999, strand: 1 } },
                ],
            }),
            text: async () => '',
        });
        const ambiguous = await liftoverPosition('4', 143094904);
        check('refuses to pick between disagreeing mappings', !ambiguous.ok && ambiguous.reason === 'unmapped');
    }

    // 7. A 200 carrying an HTML error page must not be read as "unmappable".
    {
        let n = 0;
        global.fetch = async () => {
            n++;
            if (n === 1) return { ok: true, status: 200, headers: { get: () => null }, json: async () => { throw new Error('bad json'); }, text: async () => '<html>' };
            return mapResponse('4', 142173751);
        };
        const out = await liftoverPosition('4', 143094904);
        check('a malformed 200 body is treated as a host failure, not an answer', out.ok && out.pos === 142173751);
    }

    // 8. Bad input never reaches the network.
    {
        let touched = false;
        global.fetch = async () => { touched = true; return mapResponse('4', 1); };
        const out = await liftoverPosition('', 'not-a-number');
        check('invalid coordinates short-circuit without a request', !out.ok && !touched);
    }

    // 9. The deadline is honoured — a hung provider cannot outlive the budget.
    {
        const clock = fakeClock();
        global.fetch = async () => { clock.advance(3000); return errorResponse(500); };
        let n = 0;
        const counting = global.fetch;
        global.fetch = async (...a) => { n++; return counting(...a); };
        const out = await liftoverPosition('4', 143094904, { deadlineMs: 6000, now: clock.now });
        check('stops once the deadline is spent', !out.ok && out.reason === 'unavailable' && n <= 2);
    }

    // 9b. The failure mode that motivated hedging: a sick host does not return an
    //     error, it hangs. Sequential failover would spend the whole budget waiting
    //     on it, so the mirror must be started while the primary is still pending.
    {
        global.setTimeout = realSetTimeout; // hedge delay must actually elapse
        global.fetch = async (url) => {
            if (String(url).startsWith('https://rest.ensembl.org')) {
                await new Promise((r) => realSetTimeout(r, 5000)); // hangs past the deadline
                return errorResponse(500);
            }
            await new Promise((r) => realSetTimeout(r, 50));
            return mapResponse('4', 142173751);
        };
        const t0 = Date.now();
        const out = await liftoverPosition('4', 143094904, { hedgeDelayMs: 100, attemptTimeoutMs: 3000, deadlineMs: 6000 });
        const elapsed = Date.now() - t0;
        check('a hanging primary does not block the mirror', out.ok && out.pos === 142173751);
        check('answers via the hedge without waiting out the primary', elapsed < 1000);
        global.setTimeout = (fn) => realSetTimeout(fn, 0);
    }

    // 9c. A healthy primary answers before the hedge fires, so the mirror is
    //     never called — failover must not double the load in the normal case.
    {
        global.setTimeout = realSetTimeout;
        const hostsHit = new Set();
        global.fetch = async (url) => {
            hostsHit.add(new URL(String(url)).host);
            return mapResponse('4', 142173751);
        };
        const out = await liftoverPosition('4', 143094904, { hedgeDelayMs: 1200 });
        check('healthy primary answers', out.ok && out.host === 'https://rest.ensembl.org');
        check('healthy primary means the mirror is never contacted', hostsHit.size === 1);
        global.setTimeout = (fn) => realSetTimeout(fn, 0);
    }

    check('normalizeChrom strips chr and trims', normalizeChrom(' chr17 ') === '17' && normalizeChrom(7) === '7');

    // --- gnomAD v4 handler integration ---

    function makeRes() {
        const captured = {};
        return {
            captured,
            setHeader() {},
            status(code) { captured.code = code; return this; },
            json(body) { captured.body = body; return this; },
            end() { return this; },
        };
    }

    // 10. A caller-supplied pos38 skips Ensembl entirely.
    {
        const urls = [];
        global.fetch = async (url, opts) => {
            urls.push(String(url));
            if (String(url).includes('gnomad')) {
                return { ok: true, status: 200, text: async () => JSON.stringify({ data: { variant: { variant_id: '4-142173751-G-A', chrom: '4', pos: 142173751, ref: 'G', alt: 'A', exome: { ac: 1, an: 250752, af: 0.000004 }, genome: null } } }) };
            }
            return errorResponse(500); // Ensembl is down; must not matter
        };
        const res = makeRes();
        await handler({ method: 'GET', query: { chrom: '4', pos37: '143094904', ref: 'G', alt: 'A', pos38: '142173751' } }, res);
        check('pos38 hint yields a found result while Ensembl is down', res.captured.body.status === 'found');
        check('pos38 hint makes no Ensembl request at all', urls.every((u) => !u.includes('ensembl')));
        check('pos38 hint builds the right GRCh38 variant id', res.captured.body.grch38Id === '4-142173751-G-A');
    }

    // 11. Ensembl down with no hint => liftover_unavailable (the honest status).
    {
        global.fetch = async (url) => {
            if (String(url).includes('ensembl')) return errorResponse(500);
            return { ok: true, status: 200, text: async () => JSON.stringify({ data: { variant: null } }) };
        };
        const res = makeRes();
        await handler({ method: 'GET', query: { chrom: '4', pos37: '143094904', ref: 'G', alt: 'A' } }, res);
        check('Ensembl outage surfaces as liftover_unavailable', res.captured.body.status === 'liftover_unavailable');
        check('liftover_unavailable still returns HTTP 200 so the card degrades', res.captured.code === 200);
        check('liftover_unavailable message names the provider, not the variant', /Ensembl/.test(res.captured.body.message));
    }

    // 12. Ensembl healthy but the coordinate has no equivalent => liftover_failed.
    {
        global.fetch = async (url) => {
            if (String(url).includes('ensembl')) {
                return { ok: true, status: 200, headers: { get: () => null }, json: async () => ({ mappings: [] }), text: async () => '' };
            }
            return { ok: true, status: 200, text: async () => JSON.stringify({ data: { variant: null } }) };
        };
        const res = makeRes();
        await handler({ method: 'GET', query: { chrom: '4', pos37: '143094904', ref: 'G', alt: 'A' } }, res);
        check('an unmappable coordinate still reports liftover_failed', res.captured.body.status === 'liftover_failed');
    }

    // 13. No hint, Ensembl healthy: the lifted position drives the gnomAD query.
    {
        let gnomadBody = null;
        global.fetch = async (url, opts) => {
            if (String(url).includes('ensembl')) return mapResponse('4', 142173751);
            gnomadBody = JSON.parse(opts.body);
            return { ok: true, status: 200, text: async () => JSON.stringify({ data: { variant: { variant_id: '4-142173751-G-A', chrom: '4', pos: 142173751, ref: 'G', alt: 'A', exome: null, genome: null } } }) };
        };
        const res = makeRes();
        await handler({ method: 'GET', query: { chrom: 'chr4', pos37: '143094904', ref: 'g', alt: 'a' } }, res);
        check('lifted position feeds the gnomAD variant id', gnomadBody.variables.variantId === '4-142173751-G-A');
        check('found status returned after a real liftover', res.captured.body.status === 'found');
    }

    // 14. gnomAD signals an absent variant as a GraphQL error, not a null. That
    //     must read as "not found", not as an API error.
    {
        global.fetch = async (url) => {
            if (String(url).includes('ensembl')) return mapResponse('4', 142173751);
            return { ok: true, status: 200, text: async () => JSON.stringify({ errors: [{ message: 'Variant not found' }], data: { variant: null } }) };
        };
        const res = makeRes();
        await handler({ method: 'GET', query: { chrom: '4', pos37: '143094904', ref: 'G', alt: 'A' } }, res);
        check('"Variant not found" reports as not_found, not api_error', res.captured.body.status === 'not_found');
        check('not_found still reports the GRCh38 id that was queried', res.captured.body.grch38Id === '4-142173751-G-A');
    }

    // 15. And because that check no longer returns early, the indel region scan is
    //     reachable again — the fallback that covers a liftover landing an indel a
    //     few bases from where gnomAD indexes it.
    {
        const posted = [];
        global.fetch = async (url, opts) => {
            if (String(url).includes('ensembl')) return mapResponse('4', 142173751);
            const body = JSON.parse(opts.body);
            posted.push(body.operationName);
            if (body.operationName === 'GnomadVariant') {
                return { ok: true, status: 200, text: async () => JSON.stringify({ errors: [{ message: 'Variant not found' }], data: { variant: null } }) };
            }
            return { ok: true, status: 200, text: async () => JSON.stringify({ data: { region: { variants: [
                { variant_id: '4-142173753-GAT-G', pos: 142173753, ref: 'GAT', alt: 'G', exome: { ac: 3, an: 100000, af: 0.00003 }, genome: null },
            ] } } }) };
        };
        const res = makeRes();
        await handler({ method: 'GET', query: { chrom: '4', pos37: '143094904', ref: 'GAT', alt: 'G' } }, res);
        check('an indel miss now falls through to the region scan', posted.includes('GnomadRegion'));
        check('region scan recovers the shifted indel', res.captured.body.status === 'found' && res.captured.body.grch38Id === '4-142173753-GAT-G');
        check('region-scan recovery is disclosed via matchedVia', /region scan/.test(res.captured.body.matchedVia || ''));
        check('indel results still carry the inexact-liftover caveat', /may be inexact/.test(res.captured.body.caveat || ''));
    }

    // 16. A genuine API error is still an API error.
    {
        global.fetch = async (url) => {
            if (String(url).includes('ensembl')) return mapResponse('4', 142173751);
            return { ok: true, status: 200, text: async () => JSON.stringify({ errors: [{ message: 'Cannot query field "bogus"' }] }) };
        };
        const res = makeRes();
        await handler({ method: 'GET', query: { chrom: '4', pos37: '143094904', ref: 'G', alt: 'A' } }, res);
        check('real GraphQL errors still surface as api_error', res.captured.body.status === 'api_error');
    }

    console.log(`\nLiftover tests: ${passed} passed, ${failed} failed`);
    if (failed > 0) process.exit(1);
}

run().catch((err) => { console.error(err); process.exit(1); });
