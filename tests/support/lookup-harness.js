/*
 * Replay harness for end-to-end lookup tests.
 *
 * Runs the REAL form handler from script.js (LinkeDOM + vm, no browser) against
 * a "cassette": the answers the live services gave when the scenario was
 * recorded. Nothing here reaches the network.
 *
 * A cassette is
 *   {
 *     input:     the text typed into the form,
 *     responses: { "<host><path><query>": { status, body } },   core upstreams
 *     reference: [ { chrom, start, seq } ],                     GRCh37 windows
 *     api:       { "<name>?<query>": body }      answers of the site's own /api
 *   }                                            proxies, where one decides the
 *                                                lookup (ClinVar's other names)
 * Keys are the request URL with the scheme dropped and percent-encoding undone,
 * exactly as they appear in the returned `requests` list.
 *
 * What a request gets:
 *   1. `route(url, key, respond, options, { document })` if the test supplies
 *      one and it answers (the page's document is there for a route that has
 *      to answer at a particular point in the lookup);
 *   2. the cassette's recorded response for that exact URL;
 *   3. a slice of a recorded reference window, for Ensembl sequence requests
 *      — and for UCSC's, the page's second source of those bases, where the
 *      test says UCSC is up (`ucsc: true`). By default it is not: a scenario
 *      that takes Ensembl's reference away is one without the reference;
 *      (MyVariant's batch POST is answered from the same per-id recordings);
 *   4. `backend[name](url)` for the site's own /api/<name> proxies;
 *   5. otherwise the failure that source gives when it has nothing: 404 for a
 *      MyVariant id, an empty hit list for a MyVariant search, 503 for Ensembl
 *      and 404 for every evidence source (those fail soft and have no say in
 *      which variant this is).
 * A recorded status of 0 is a request that never completed (timeout).
 */
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { parseHTML } from 'linkedom';

const script = readFileSync(new URL('../../script.js', import.meta.url), 'utf8');
const html = readFileSync(new URL('../../index.html', import.meta.url), 'utf8');

// In a query string "+" is a space (URLSearchParams writes it that way).
export const requestKey = (url) => `${url.hostname}${decodeURIComponent(url.pathname)}${decodeURIComponent(url.search.replace(/\+/g, ' '))}`;

export function referenceSlice(windows, chrom, start, end) {
    const win = [].concat(windows || []).find(w => String(w.chrom) === String(chrom)
        && w.start <= start && w.start + w.seq.length - 1 >= end);
    return win ? win.seq.slice(start - win.start, end - win.start + 1) : null;
}

// LinkeDOM implements no table-section helpers; two evidence cards build their
// tables with them.
function addTableHelpers(document) {
    const helpers = [['createTHead', 'thead'], ['createTBody', 'tbody'], ['insertRow', 'tr'], ['insertCell', 'td']];
    for (const tag of ['table', 'thead', 'tbody', 'tr']) {
        const proto = Object.getPrototypeOf(document.createElement(tag));
        for (const [method, child] of helpers) {
            if (!proto[method]) {
                proto[method] = function () {
                    const el = this.ownerDocument.createElement(child);
                    this.appendChild(el);
                    return el;
                };
            }
        }
    }
}

export async function runLookup(cassette, {
    input = cassette.input, tumorType = '', grch38 = false,
    route = () => undefined, backend = {}, without = [], settleMs = 25,
    // Every wait the page sets is cut to this many milliseconds. A test that
    // needs something to land DURING one of the page's waits widens it.
    timerCap = 5,
    // A second lookup on the same page, submitted once the first has finished:
    // { input }. What is returned then describes the page after the second.
    followUp = null,
    // UCSC's sequence service answers from the cassette's reference windows.
    ucsc = false
} = {}) {
    const { document } = parseHTML(html);
    addTableHelpers(document);
    const form = document.getElementById('variantForm');
    const listen = form.addEventListener.bind(form);
    let submit;
    form.addEventListener = (event, fn, ...args) => {
        if (event === 'submit') submit = fn;
        return listen(event, fn, ...args);
    };
    const requests = [], errors = [], timers = [];
    const respond = (data, status = 200) => ({
        ok: status >= 200 && status < 300, status,
        json: async () => (typeof data === 'string' ? JSON.parse(data) : structuredClone(data)),
        text: async () => (typeof data === 'string' ? data : JSON.stringify(data))
    });
    // `without` drops recorded answers (a substring or RegExp per entry) so a
    // test can take one service away from an otherwise complete recording.
    const withheld = key => without.some(w => (w instanceof RegExp ? w.test(key) : key.includes(w)));
    const fetch = async (address, options = {}) => {
        const url = new URL(String(address), 'https://drdoubleb.com/variant/');
        const key = requestKey(url);
        // (Its place in the log: another request may be logged while `route`
        // is awaited, and the batch POST below rewrites its own entry.)
        const logged = requests.push(key) - 1;
        const routed = await route(url, key, respond, options, { document });
        if (routed) return routed;
        const path = decodeURIComponent(url.pathname);
        // MyVariant's batch endpoint: one POST for several ids, answered from the
        // per-id recordings. Logged with the ids it asked for.
        if (url.hostname === 'myvariant.info' && path === '/v1/variant' && String(options.method).toUpperCase() === 'POST') {
            const ids = new URLSearchParams(String(options.body || '')).get('ids').split(',');
            requests[logged] = `${key} [POST ids=${ids.join(',')}]`;
            if (withheld(key)) return respond({ error: 'Service Unavailable' }, 503);
            return respond(ids.map((id) => {
                const recorded = cassette.responses?.[`myvariant.info/v1/variant/${id}`];
                return recorded && recorded.status === 200 ? { query: id, ...recorded.body } : { query: id, notfound: true };
            }));
        }
        const api = path.match(/\/api\/([a-z0-9-]+)$/);
        if (api && /vercel\.app$/.test(url.hostname)) {
            if (backend[api[1]]) {
                const answer = await backend[api[1]](url, options);
                return answer ? respond(answer) : respond({ error: 'not recorded' }, 404);
            }
            // A proxy answer recorded with the scenario (cassette.api), keyed
            // "name?query" the way `page.api(name)` reports requests.
            const recordedApi = withheld(key) ? null : cassette.api?.[`${api[1]}${decodeURIComponent(url.search.replace(/\+/g, ' '))}`];
            return recordedApi ? respond(recordedApi) : respond({ error: 'not recorded' }, 404);
        }
        const recorded = withheld(key) ? null : cassette.responses?.[key];
        if (recorded) {
            if (!recorded.status) throw Object.assign(new Error('The operation was aborted'), { name: 'AbortError' });
            return respond(recorded.body, recorded.status);
        }
        const seq = path.match(/\/sequence\/region\/human\/([^:]+):(\d+)\.\.(\d+)$/);
        if (seq && !withheld(key)) {
            const bases = referenceSlice(cassette.reference, seq[1], Number(seq[2]), Number(seq[3]));
            if (bases !== null) return respond(bases);
        }
        // UCSC's copy of the assembly, in the shape its API answers in: the
        // start counted from 0, chromosomes named "chr…", and (as it writes
        // repeats) lower case — here all of it, so that a page reading the
        // bases as given would stumble.
        if (ucsc && url.hostname === 'api.genome.ucsc.edu' && path === '/getData/sequence' && !withheld(key)) {
            const asked = url.searchParams;
            const start = Number(asked.get('start'));
            const end = Number(asked.get('end'));
            const bases = referenceSlice(cassette.reference, String(asked.get('chrom')).replace(/^chr/, ''), start + 1, end);
            if (bases !== null) {
                return respond({ downloadTime: '2026:11:18T14:40:07Z', downloadTimeStamp: 1795012807,
                    genome: asked.get('genome'), chrom: asked.get('chrom'), start, end, dna: bases.toLowerCase() });
            }
        }
        if (url.hostname === 'myvariant.info') {
            return path.startsWith('/v1/variant/') ? respond({ success: false, error: 'Not Found.' }, 404) : respond({ hits: [] });
        }
        if (/ensembl\.org$/.test(url.hostname)) return respond({ error: 'Service Unavailable' }, 503);
        return respond({}, 404);
    };
    const window = { location: new URL('https://drdoubleb.com/variant/'),
        history: { replaceState() {} }, matchMedia: () => ({ matches: false }) };
    const context = vm.createContext({ document, window, fetch,
        console: { log() {}, warn() {}, error(e) { errors.push(String(e)); } },
        URL, URLSearchParams, AbortController, performance, structuredClone,
        // Timeouts and retry back-offs collapse to a few milliseconds.
        setTimeout(fn, ms, ...args) { const timer = setTimeout(fn, Math.min(ms, timerCap), ...args); timers.push(timer); return timer; },
        clearTimeout
    });
    try {
        vm.runInContext(script, context);
        vm.runInContext('initVariantSearchUi()', context);
        document.getElementById('variantInput').value = input;
        if (tumorType) document.getElementById('tumorTypeInput').value = tumorType;
        if (grch38) document.getElementById('grch38Toggle').checked = true;
        await submit({ preventDefault() {} });
        if (followUp) {
            document.getElementById('variantInput').value = followUp.input;
            await submit({ preventDefault() {} });
        }
        await new Promise(resolve => setTimeout(resolve, settleMs));
        const text = el => (el ? el.textContent.replace(/\s+/g, ' ').trim() : '');
        const cards = Object.fromEntries([...document.querySelectorAll('.card')]
            .map(el => [el.querySelector('h3')?.textContent || '', text(el)]));
        const links = [...document.querySelectorAll('a')].map(a => ({ text: text(a), href: a.getAttribute('href') || '' }));
        return {
            status: document.getElementById('status').textContent,
            hidden: document.getElementById('resultSection').classList.contains('hidden'),
            panel: text(document.getElementById('lookupProgress')),
            variant: cards.Variant || '', clinvar: cards.ClinVar || '',
            cards, links, requests, errors, document,
            // The search half of a link: what a click would ask that site for.
            searches: site => links.filter(l => l.href.includes(site))
                .map(l => decodeURIComponent((l.href.split('?')[1] || '').replace(/\+/g, ' '))),
            // Requests to one of the site's own /api proxies, as "name?query".
            api: name => requests.filter(r => r.includes(`/api/${name}?`)).map(r => r.slice(r.indexOf(`/api/${name}?`) + 5))
        };
    } finally {
        timers.forEach(clearTimeout);
    }
}

export function assertResolved(assert, page, label = 'lookup') {
    assert.equal(page.status, '', `${label}: ${page.errors.join('\n')}`);
    assert.equal(page.hidden, false, `${label}: result rendered`);
}
