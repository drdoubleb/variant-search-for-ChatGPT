/*
 * Names ClinVar knows a variant by.
 *
 * "BRCA1 5382insC", "BRCA1 185delAG", "BRCA2 6174delT", "MTHFR C677T": what
 * these variants were called before HGVS, and still what a report says. None
 * of them is HGVS and no nomenclature service resolves them:
 *
 *   BRCA1 5382insC             sent on as "BRCA1:p.5382insC"; lookup failed
 *   BRCA1 185delAG             read as a protein change, "p.185delAG"; failed
 *   BRCA1 p.Gln1756ProfsTer74  Ensembl takes no frameshift protein input; the
 *                              page fell back to "residues 74-1756", a 47 kb
 *                              stretch of BRCA1
 *   MTHFR C677T                read as p.Cys677Thr; failed
 *
 * ClinVar lists them as other names of the variant. A name that will not
 * resolve is looked up there; when exactly one variant of the gene goes by it,
 * the lookup runs again for that variant's HGVS name and says so.
 *
 * The scenarios replay answers recorded from the live services, including the
 * site's own ClinVar proxy.
 *
 * Run with: node tests/other-names.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';
import clinvarProteinHandler from '../api/clinvar-protein.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const ncbi = JSON.parse(readFileSync(new URL('./fixtures/clinvar-names.json', import.meta.url)));

// ── One spelling per name ──────────────────────────────────────────────────
for (const [a, b] of [
    ['5382insC', 'c.5382insC'], ['5382insC', '5382insc'], ['185delAG', 'NM_007294.4(BRCA1):c.185delAG'],
    ['p.Gln1756ProfsTer74', 'NP_009225.1:p.Gln1756ProfsTer74'], ['p.Gln1756ProfsTer74', 'p.(Gln1756ProfsTer74)'],
    ['p.Gln1756ProfsTer74', 'p.Gln1756ProfsX74'], ['Q1756Pfs*74', 'p.Gln1756ProfsTer74'], ['Gln1756Profs*74', 'Q1756PfsX74'],
    ['p.Glu23fs*17', 'Glu23fsTer17'], ['R213X', 'p.Arg213Ter'], ['R213*', 'Arg213*'], ['V600E', 'p.Val600Glu']
]) {
    assert.equal(h.variantNameKey(a), h.variantNameKey(b), `${a} is ${b}`);
}
for (const [a, b] of [['5382insC', '5385insC'], ['185delAG', '187delAG'], ['Q1756Pfs*74', 'Q1756fs'], ['V600E', 'V600K'], ['c.68_69del', '68_69delAG']]) {
    assert.notEqual(h.variantNameKey(a), h.variantNameKey(b), `${a} is not ${b}`);
}
// The old way of writing a nucleotide substitution, and ClinVar's OMIM-style lists.
const shares = (a, b) => Array.from(h.variantNameKeys(a)).some(key => h.variantNameKeys(b).has(key));
assert(shares('C677T', '677C>T'));
assert(shares('C677T', '677C-T'));
assert(shares('C677T', 'MTHFR, 677C-T, ALA222VAL (rs1801133)'), 'as ClinVar lists it');
assert(shares('A222V', 'MTHFR, 677C-T, ALA222VAL (rs1801133)'));
assert(!shares('C677T', 'C667T'), 'ClinVar\'s own typo for this variant is not it');
assert(!shares('C677T', 'A1298C'));
assert(!shares('C677T', 'MTHFR'));

// ── What ClinVar's search is asked ─────────────────────────────────────────
assert.equal(h.clinvarSearchTextForName('5382insC'), '5382insC');
assert.equal(h.clinvarSearchTextForName('c.5382insC'), '5382insC');
assert.equal(h.clinvarSearchTextForName('p.Gln1756ProfsTer74'), 'Gln1756ProfsTer74');
assert.equal(h.clinvarSearchTextForName('Q1756Pfs*74'), 'Gln1756ProfsTer74', 'written out: "*" is a wildcard to ClinVar');
assert.equal(h.clinvarSearchTextForName('Q1756fs'), 'Gln1756fs');
assert.equal(h.clinvarSearchTextForName('C677T'), 'C677T');
assert.equal(h.clinvarSearchTextForName(''), '');
assert.deepEqual(h.parseClinvarTitle('NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs)'), { accession: 'NM_007294.4', gene: 'BRCA1', change: 'c.5266dup' });
assert.deepEqual(h.parseClinvarTitle('NM_002875.5(RAD51):c.-98G>C'), { accession: 'NM_002875.5', gene: 'RAD51', change: 'c.-98G>C' });
assert.equal(h.parseClinvarTitle('NC_012920.1(MT-ATP6):m.8993T>G'), null, 'no transcript-level name to look up');
assert.equal(h.parseClinvarTitle(''), null);

// ── The proxy: ClinVar's search, with each record's other names ────────────
async function proxy(query) {
    const realFetch = globalThis.fetch;
    const asked = [];
    globalThis.fetch = async (address) => {
        const url = String(address);
        asked.push(url);
        if (url.includes('esearch.fcgi')) return Response.json(ncbi.esearch);
        if (url.includes('esummary.fcgi')) return Response.json(ncbi.esummary);
        throw new Error(`unexpected upstream request ${url}`);
    };
    try {
        const res = { statusCode: 200, payload: null, setHeader() { return this; }, status(c) { this.statusCode = c; return this; }, json(v) { this.payload = v; return this; }, end() { return this; } };
        await clinvarProteinHandler({ method: 'GET', headers: { origin: 'https://www.drdoubleb.com' }, query }, res);
        return { status: res.statusCode, payload: res.payload, asked };
    } finally {
        globalThis.fetch = realFetch;
    }
}
const searched = await proxy({ gene: 'BRCA1', change: '5382insC', retmax: '25' });
assert.equal(searched.status, 200);
assert.match(decodeURIComponent(searched.asked[0]), /term=BRCA1 5382insC/);
assert.match(searched.asked[0], /retmax=25/);
// As recorded: five records, one of them the variant, one in another gene.
assert.deepEqual(searched.payload.variants.map(v => v.id), ['37560', '17677', '17662', '13128', '9325']);
const founder = searched.payload.variants.find(v => v.id === '17677');
assert.equal(founder.title, 'NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs)');
assert(founder.aliases.includes('5382insC') && founder.aliases.includes('5385insC'));
assert.equal(founder.variationCount, 1);

// ── Which record goes by the name ──────────────────────────────────────────
{
    const variants = searched.payload.variants;
    const expected = { id: '17677', title: 'NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs)', hgvs: 'NM_007294.4:c.5266dup' };
    for (const name of ['5382insC', 'c.5382insC', '5382insc', '5385insC', 'p.Gln1756ProfsTer74', 'Q1756Pfs*74', 'p.(Gln1756ProfsTer74)']) {
        assert.deepEqual(h.clinvarVariantNamed(variants, 'BRCA1', name), expected, name);
    }
    assert.deepEqual(h.clinvarVariantNamed(variants, 'brca1', '185delAG').id, '17662');
    assert.equal(h.clinvarVariantNamed(variants, 'BRCA1', '4154delA').hgvs, 'NM_007294.4:c.4035del');
    // Found by the search but not called that; not in this gene; not a name.
    assert.equal(h.clinvarVariantNamed(variants, 'BRCA1', '5382insG'), null);
    assert.equal(h.clinvarVariantNamed(variants, 'BRCA1', '6174delT'), null, 'that is a BRCA2 variant');
    assert.equal(h.clinvarVariantNamed(variants, 'BRCA2', '6174delT').id, '9325');
    assert.equal(h.clinvarVariantNamed(variants, 'BRCA1', '135G-C'), null, 'RAD51');
    assert.equal(h.clinvarVariantNamed(variants, 'BRCA1', ''), null);
    assert.equal(h.clinvarVariantNamed([], 'BRCA1', '5382insC'), null);
    // Two variants under one name: not for this code to choose between.
    const twin = { ...founder, id: '999999' };
    assert.equal(h.clinvarVariantNamed([...variants, twin], 'BRCA1', '5382insC'), null);
    // A haplotype that contains the variant is not the variant.
    const alone = variants.filter(v => v.id !== '17677');
    assert.equal(h.clinvarVariantNamed([...alone, { ...founder, variationCount: 2 }], 'BRCA1', '5382insC'), null);
    // A backend that does not send the names yet (deployed before this change): nothing matches.
    const { aliases, variationCount, ...older } = founder;
    assert.equal(h.clinvarVariantNamed([older], 'BRCA1', '5382insC'), null);
}
{
    const calls = [];
    const fetcher = async (gene, forms, retmax) => { calls.push([gene, forms, retmax]); return { variants: searched.payload.variants }; };
    assert.equal((await h.resolveClinvarOtherName('BRCA1', 'c.5382insC', fetcher)).hgvs, 'NM_007294.4:c.5266dup');
    assert.deepEqual(calls, [['BRCA1', ['5382insC'], 25]], 'asked for the few records that match, not the gene\'s first 200');
    assert.equal(await h.resolveClinvarOtherName('BRCA1', '9999delT', fetcher), null, 'not listed');
    const warn = console.warn;
    console.warn = () => {};
    assert.equal(await h.resolveClinvarOtherName('BRCA1', '5382insC', async () => { throw new Error('502'); }), undefined, 'could not ask');
    console.warn = warn;
    assert.equal(await h.resolveClinvarOtherName('', '5382insC', fetcher), null);
}

// ── A frameshift's stop distance is not a residue ──────────────────────────
assert.deepEqual(h.parseProteinResidueRange('p.Gln1756ProfsTer74'), { start: 1756, end: 1756 }, 'was 74-1756');
assert.deepEqual(h.parseProteinResidueRange('p.N1651Mfs*21'), { start: 1651, end: 1651 });
assert.deepEqual(h.parseProteinResidueRange('p.Asn1651MetfsTer21'), { start: 1651, end: 1651 });
assert.deepEqual(h.parseProteinResidueRange('p.Pro153fs'), { start: 153, end: 153 });
assert.deepEqual(h.parseProteinResidueRange('p.Ter394GlnextTer12'), { start: 394, end: 394 });
assert.deepEqual(h.parseProteinResidueRange('p.Leu773_Ile774delinsPhe'), { start: 773, end: 774 });
assert.deepEqual(h.parseProteinResidueRange('p.Lys18del'), { start: 18, end: 18 });
// ── NCBI alongside a slow Ensembl ──────────────────────────────────────────
// (Deadlines shortened for the test; the real ones are 2.5 s and 12 s.)
{
    const realFetch = globalThis.fetch;
    const real = { head: h.API_TIMEOUT_MS.recoderHeadStart, recoder: h.API_TIMEOUT_MS.recoder };
    h.API_TIMEOUT_MS.recoderHeadStart = 40;
    h.API_TIMEOUT_MS.recoder = 150;
    const placed = { data: { spdis: [{ seq_id: 'NC_000017.10', position: 41276044, deleted_sequence: 'CTCT', inserted_sequence: 'CT' }] } };
    const ensemblAnswer = [{ '-': { input: 'NM_007294.4:c.68_69del', hgvsg: ['NC_000017.10:g.41276047_41276048del'], hgvsc: ['NM_007294.4:c.68_69del'] } }];
    const run = async (ensembl) => {
        const asked = [];
        globalThis.fetch = (address, options = {}) => {
            const url = String(address);
            asked.push(url.includes('api.ncbi.nlm.nih.gov') ? 'ncbi' : 'ensembl');
            if (url.includes('api.ncbi.nlm.nih.gov')) return Promise.resolve(Response.json(placed));
            return ensembl(options);
        };
        const started = Date.now();
        try {
            return { answer: await h.fetchVariantRecoder('NM_007294.4:c.68_69del'), asked, ms: Date.now() - started };
        } catch (error) {
            return { error, asked, ms: Date.now() - started };
        }
    };
    // Takes longer than the head start. (It then answers 400, which ends
    // Ensembl's retries, so nothing of this case outlives the test.)
    const slowly = () => new Promise(resolve => setTimeout(() => resolve(new Response('{"error":"late"}', { status: 400 })), 110));
    try {
        // Ensembl answers within its head start: NCBI is not asked.
        const prompt = await run(() => Promise.resolve(Response.json(ensemblAnswer)));
        assert.deepEqual(prompt.answer, ensemblAnswer);
        assert.deepEqual(prompt.asked, ['ensembl']);
        // Ensembl still working after its head start: NCBI's answer is taken,
        // long before Ensembl's four attempts would have run out.
        const slow = await run(slowly);
        assert.equal(slow.answer.__vsFallback, true);
        assert.deepEqual(Object.values(slow.answer[0])[0].spdi.length > 0, true);
        assert(slow.asked.includes('ncbi'));
        assert(slow.ms < 105, `answered after the head start, without waiting for Ensembl (${slow.ms} ms)`);
        await new Promise(resolve => setTimeout(resolve, 90));
        // Ensembl never answers. Once NCBI's answer has been taken nothing reads
        // Ensembl's, so its remaining attempts are called off: asked once, not
        // the four times (two per host) an unanswered lookup is allowed.
        const hangs = (options) => new Promise((_, reject) => {
            options.signal.addEventListener('abort', () => reject(Object.assign(new Error('aborted'), { name: 'AbortError' })));
        });
        const hung = await run(hangs);
        assert.equal(hung.answer.__vsFallback, true);
        // Past the first deadline (150 ms) and the longest backoff before a retry (520 ms).
        await new Promise(resolve => setTimeout(resolve, 800));
        assert.deepEqual(hung.asked.filter(source => source === 'ensembl'), ['ensembl'], 'Ensembl was not asked again');
        // Ensembl rejects the name at once: NCBI is asked straight away, and
        // when it has nothing either, it is Ensembl's reason that is reported.
        globalThis.fetch = (address) => Promise.resolve(String(address).includes('api.ncbi.nlm.nih.gov')
            ? new Response('{}', { status: 400 })
            : new Response('{"error":"Unable to parse HGVS notation"}', { status: 400 }));
        await assert.rejects(h.fetchVariantRecoder('NM_007294.4:c.5382insC'), err => err.status === 400 && /Unable to parse/.test(err.message));
    } finally {
        globalThis.fetch = realFetch;
        h.API_TIMEOUT_MS.recoderHeadStart = real.head;
        h.API_TIMEOUT_MS.recoder = real.recoder;
    }
}
console.log('Other-name unit checks passed');

const nameSearches = page => page.api('clinvar-protein').filter(r => /retmax=25/.test(r));
const everywhere = page => [page.status, page.panel, page.variant, ...page.requests].join('\n');

// 1. The Ashkenazi founder mutation by the name it is known by, four ways.
for (const typed of ['BRCA1 5382insC', 'BRCA1 c.5382insC', 'BRCA1 p.Gln1756ProfsTer74', 'BRCA1 Q1756Pfs*74']) {
    const page = await runLookup(cassettes.brca1OtherNames, { input: typed });
    assertResolved(assert, page, typed);
    const name = typed.split(' ')[1];
    assert.match(page.variant, /g\.: chr17:g\.41209082dup/, typed);
    assert.match(page.variant, /Transcript: NM_007294\.4 \(MANE Select\)c\.: c\.5266dupp\.: p\.Gln1756ProfsTer74/, typed);
    assert.match(page.clinvar, /Variation ID: 17677/, typed);
    // It says what the name was taken to mean, in the panel and on the card.
    assert(page.panel.includes(`Other names · ClinVar"${name}" is another name for NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs) (Variation 17677)`), typed);
    assert(page.document.querySelector('.other-name-note').textContent.includes(
        `"${name}" is another name ClinVar lists for NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs) (Variation 17677)`), typed);
    assert.match(page.panel, /You entered\s*BRCA1 /, `${typed}: what was typed is still shown`);
    assert.equal(nameSearches(page).length, 1, `${typed}: ClinVar's names are consulted once`);
    assert.doesNotMatch(everywhere(page), /p\.5382insC|Codon region|residue 74/, typed);
}
console.log('ok other names: 5382insC, c.5382insC and the protein-level frameshift are BRCA1 c.5266dup');

// 2. 185delAG. Was: read as "p.185delAG".
{
    const page = await runLookup(cassettes.brca1Legacy185delAG);
    assertResolved(assert, page, 'BRCA1 185delAG');
    assert.match(page.variant, /g\.: chr17:g\.41276047_41276048del/);
    assert.match(page.variant, /c\.: c\.68_69delp\.: p\.Glu23ValfsTer17/);
    assert.match(page.panel, /"185delAG" is another name for NM_007294\.4\(BRCA1\):c\.68_69del \(p\.Glu23fs\) \(Variation 17662\)/);
    assert.deepEqual(nameSearches(page), ['clinvar-protein?gene=BRCA1&change=185delAG&retmax=25']);
    assert.doesNotMatch(everywhere(page), /p\.185delAG/);
    // Nothing is asked of Ensembl or MyVariant about a name that is not HGVS.
    assert(!page.requests.some(r => /185delAG/.test(r) && !/clinvar-protein/.test(r)));
    console.log('ok other names: 185delAG is BRCA1 c.68_69del');
}

// 3. A nucleotide name that looks like a protein change. Read as p.Cys677Thr
//    first, as before; when that is no variant, ClinVar's names are asked.
{
    const page = await runLookup(cassettes.mthfrC677T);
    assertResolved(assert, page, 'MTHFR C677T');
    assert.match(page.variant, /g\.: chr1:g\.11856378G>A/);
    assert.match(page.variant, /Transcript: NM_005957\.5 \(MANE Select\)c\.: c\.665C>Tp\.: p\.Ala222Val/);
    assert.match(page.panel, /"C677T" is another name for NM_005957\.5\(MTHFR\):c\.665C>T \(p\.Ala222Val\) \(Variation 3520\)/);
    console.log('ok other names: MTHFR C677T is c.665C>T');
}

// 4. A name ClinVar does not list.
{
    const page = await runLookup(cassettes.brca1UnknownName);
    assert.equal(page.status, 'Error: "BRCA1 9999delT" is not a name this tool can read. A variant is written in HGVS — '
        + 'BRCA1 c.…, NM_…:c.…, chr7:g.… — or as a coordinate row. ClinVar lists no BRCA1 variant under that name either.');
    assert.equal(page.hidden, true);
    assert.match(page.panel, /✕Parse input\s*not HGVS/);
    assert.match(page.panel, /◍Other names · ClinVar\s*no single BRCA1 variant listed under "9999delT"/);
    assert.deepEqual(page.requests.filter(r => /myvariant|ensembl/.test(r)), [], 'no other service is asked about it');
    console.log('ok other names: an unknown name is refused, and says what would be read');
}

// 5. ClinVar's names cannot be consulted: no claim is made about what it lists.
{
    const page = await runLookup(cassettes.brca1OtherNames, { backend: { 'clinvar-protein': () => null } });
    assert.equal(page.status, 'Error: "BRCA1 5382insC" is not a name this tool can read. A variant is written in HGVS — '
        + 'BRCA1 c.…, NM_…:c.…, chr7:g.… — or as a coordinate row.');
    assert.match(page.panel, /▲Other names · ClinVar\s*could not be consulted/);
    console.log('ok other names: proxy down is "could not be consulted", not "not listed"');
}

// 6. The backend as deployed before this change sends no names: nothing is
//    resolved, and nothing is wrongly resolved either.
{
    const older = structuredClone(cassettes.brca1OtherNames.api['clinvar-protein?gene=BRCA1&change=5382insC&retmax=25']);
    for (const variant of older.variants) { delete variant.aliases; delete variant.variationCount; }
    const page = await runLookup(cassettes.brca1OtherNames, { backend: { 'clinvar-protein': () => older } });
    assert.match(page.status, /^Error: "BRCA1 5382insC" is not a name this tool can read/);
    assert.equal(page.hidden, true);
    console.log('ok other names: an older backend degrades to the plain error');
}

// 7. A name that resolves the ordinary way never consults ClinVar's names.
{
    const page = await runLookup(cassettes.kras35);
    assertResolved(assert, page, 'KRAS c.35G>A');
    assert.deepEqual(nameSearches(page), []);
    assert.equal(page.document.querySelector('.other-name-note'), null);
    console.log('ok other names: not consulted for a name that resolves');
}

console.log('Other-name tests passed');
