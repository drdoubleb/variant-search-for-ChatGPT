// Public API regression tests: reject malformed/oversized inputs without
// reaching external data services or consuming OpenRouter credits.
import assert from 'node:assert/strict';
import { boundedQueryText, positiveQueryInteger, validatedChromosome, validatedGene, validatedAllele } from '../api/_input.js';
import gnomadHandler, { parseGnomadQuery } from '../api/gnomad-v4.js';
import clinvarGeneHandler from '../api/clinvar-gene.js';
import clinvarProteinHandler from '../api/clinvar-protein.js';
import clinvarRegionHandler from '../api/clinvar-region.js';
import pubmedHandler from '../api/pubmed.js';
import aiReviewHandler from '../api/ai-review.js';

function fakeResponse() {
    return {
        statusCode: 200,
        payload: null,
        headers: {},
        setHeader(name, value) { this.headers[name] = value; return this; },
        status(code) { this.statusCode = code; return this; },
        json(value) { this.payload = value; return this; },
        end() { return this; }
    };
}

async function expectStatus(handler, query, status, options = {}) {
    const req = {
        method: 'GET',
        headers: { origin: 'https://www.drdoubleb.com' },
        query,
        ...options
    };
    const res = fakeResponse();
    await handler(req, res);
    assert.equal(res.statusCode, status, JSON.stringify({ query, payload: res.payload }));
}

// Helpers don't silently coerce arrays/objects from repeated query params.
assert.equal(boundedQueryText(['BRAF'], 64), null);
assert.equal(boundedQueryText(' BRAF ', 64), 'BRAF');
assert.equal(boundedQueryText('A'.repeat(65), 64), null);
assert.equal(positiveQueryInteger('140453136'), 140453136);
assert.equal(positiveQueryInteger('140453136oops'), null);
assert.equal(positiveQueryInteger('-1'), null);
assert.equal(positiveQueryInteger('0'), null);
assert.equal(validatedChromosome('chr7'), '7');
assert.equal(validatedChromosome('chrX'), 'X');
assert.equal(validatedChromosome('chr23'), null);
assert.equal(validatedGene('MT-ND1'), 'MT-ND1');
assert.equal(validatedGene('../BRAF'), null);
assert.equal(validatedAllele('atg'), 'ATG');

// Standard variant lookups still parse, including optional GRCh38 position.
const good = { chrom: 'chr7', pos37: '140453136', ref: 'a', alt: 't', pos38: '140753336' };
assert.deepEqual(parseGnomadQuery(good), {
    ok: true, chrom: '7', pos37: '140453136', ref: 'A', alt: 'T', pos38: '140753336'
});
assert.equal(parseGnomadQuery({ ...good, pos38: undefined }).ok, true);
assert.equal(parseGnomadQuery({ ...good, chrom: '../etc' }).ok, false);
assert.equal(parseGnomadQuery({ ...good, pos37: '3.14' }).ok, false);
assert.equal(parseGnomadQuery({ ...good, pos38: '123bad' }).ok, false);
assert.equal(parseGnomadQuery({ ...good, ref: 'A'.repeat(1001) }).ok, false);
assert.equal(parseGnomadQuery({ ...good, chrom: ['7', 'X'] }).ok, false);

await expectStatus(gnomadHandler, { ...good, pos37: 'abc' }, 400);
await expectStatus(gnomadHandler, good, 405, { method: 'POST' });
await expectStatus(clinvarGeneHandler, { gene: 'BRAF'.repeat(100) }, 400);
await expectStatus(clinvarGeneHandler, { gene: ['BRAF', 'KRAS'] }, 400);
await expectStatus(clinvarProteinHandler, { gene: 'BRAF', change: 'A'.repeat(129) }, 400);
await expectStatus(clinvarProteinHandler, { gene: 'BRAF', change: 'G34R', retmax: '200', extra: 'ignored' }, 400, { method: 'PUT' });
await expectStatus(clinvarRegionHandler, { chrom: 'chr23', pos: '140453136' }, 400);
await expectStatus(clinvarRegionHandler, { chrom: 'chr7', pos: '140453136bogus' }, 400);
await expectStatus(clinvarRegionHandler, { chrom: '7', pos: '140453136', window: '1000000000000' }, 400);
await expectStatus(pubmedHandler, { term: 'A'.repeat(513) }, 400);
await expectStatus(pubmedHandler, { term: 'BRAF', limit: '3oops' }, 400);
await expectStatus(pubmedHandler, { gene: 'TP53', variant: ['R175H', 'R248Q'] }, 400);

// Content-Length may be absent/spoofed; the parsed JSON size must also be checked.
await expectStatus(aiReviewHandler, {}, 413, {
    method: 'POST', body: { mode: 'prompt', context: { note: 'X'.repeat(4 * 1024 * 1024) } },
    headers: { origin: 'https://www.drdoubleb.com', 'content-length': '0' }
});
await expectStatus(aiReviewHandler, {}, 400, {
    method: 'POST', body: [],
    headers: { origin: 'https://www.drdoubleb.com' }
});
console.log('Security input-validation tests passed');
