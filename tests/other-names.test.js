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
 * Found in review, and fixed here:
 *
 *   CFTR 621+1G>T   searched as "6211GT", a word ClinVar has never seen. It
 *                   answered for CFTR alone - 25 of 9,631 records - and on the
 *                   strength of that the page said "ClinVar lists no CFTR
 *                   variant under that name either". It lists one.
 *   CFTR F508del    Ensembl's allele for p.Phe508del is c.1522_1524del; the
 *                   variant known as F508del is c.1521_1523del. ClinVar says
 *                   which, in the 96th of 96 records (newest first).
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
// The marks a name is written with are what ClinVar matches it by.
assert.equal(h.clinvarSearchTextForName('621+1G>T'), '621+1G>T', 'was "6211GT"');
assert.equal(h.clinvarSearchTextForName('c.1585-1G>A'), '1585-1G>A');
assert.equal(h.clinvarSearchTextForName('3849+10kbC>T'), '3849+10kbC>T');
assert.equal(h.clinvarSearchTextForName('IVS8-5T'), 'IVS8-5T');
assert.equal(h.clinvarSearchTextForName('c.68_69del'), '68_69del');
assert.equal(h.clinvarSearchTextForName('R213*'), 'R213', 'never the wildcard');
assert.equal(h.clinvarSearchTextForName('"621+1G>T" OR x[gene]'), '621+1G>TORxgene', 'nothing that is search syntax');
// A protein-level indel over a range, in the codes ClinVar's records carry.
assert.equal(h.clinvarSearchTextForName('E746_A750del'), 'Glu746_Ala750del');
assert.equal(h.clinvarSearchTextForName('p.E746_A750del'), 'Glu746_Ala750del');
assert.equal(h.clinvarSearchTextForName('L747_P753delinsS'), 'Leu747_Pro753delinsSer');
assert.equal(h.clinvarSearchTextForName('A763_Y764insFQEA'), 'Ala763_Tyr764insPheGlnGluAla');
assert.equal(h.clinvarSearchTextForName('Y772_A775dup'), 'Tyr772_Ala775dup');
assert.equal(h.clinvarSearchTextForName('p.Glu746_Ala750del'), 'Glu746_Ala750del');
assert.equal(h.clinvarSearchTextForName('F508del'), 'F508del', 'a name in its own right: as typed');
assert.equal(h.clinvarSearchTextForName('B746_A750del'), 'B746_A750del', 'not residues: left alone');
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
const searched = await proxy({ gene: 'BRCA1', change: '5382insC', retmax: '200' });
assert.equal(searched.status, 200);
assert.match(decodeURIComponent(searched.asked[0]), /term=BRCA1 5382insC/);
assert.match(searched.asked[0], /retmax=200/);
assert.deepEqual(searched.payload.notFound, [], 'ClinVar knew every word of the search');
assert.equal(searched.payload.total, 5);
{
    // A name is searched as it is written. "+" and ">" used to be stripped at
    // this end too, and ClinVar has never seen "6211GT".
    const marked = await proxy({ gene: 'CFTR', change: '621+1G>T', retmax: '200' });
    assert.match(decodeURIComponent(marked.asked[0]), /term=CFTR 621\+1G>T(&|$)/);
    assert(marked.asked[0].includes('621%2B1G%3ET'), 'and encoded, so "+" is not read as a space');
    const ranged = await proxy({ gene: 'BRCA1', change: '68_69delAG,IVS5-11T>G' });
    assert.match(decodeURIComponent(ranged.asked[0]), /term=BRCA1 68_69delAG(&|$)/);
    // Nothing that is search syntax gets through, and a form of marks alone is no form.
    const hostile = await proxy({ gene: 'BRCA1', change: '5382insC"[gene] OR (TP53)' });
    assert.match(decodeURIComponent(hostile.asked[0]), /term=BRCA1 5382insCgeneORTP53(&|$)/);
    assert.equal((await proxy({ gene: 'BRCA1', change: '+->' })).status, 400);
}
{
    // A word ClinVar has never seen is dropped from the search, which then
    // answers for the gene alone: 85,659 records for "BRCA1 9999delT". Those
    // are not variants with the change asked for. The proxy says which word
    // was dropped and returns nothing - without fetching 200 summaries.
    const realFetch = globalThis.fetch;
    const asked = [];
    globalThis.fetch = async (address) => {
        asked.push(String(address));
        if (String(address).includes('esearch.fcgi')) {
            return Response.json({ esearchresult: { count: '85659', idlist: ['4690000', '4689999', '4689998'],
                errorlist: { phrasesnotfound: ['9999delT'], fieldsnotfound: [] }, querytranslation: 'BRCA1[All Fields]' } });
        }
        throw new Error('no summaries are fetched for records that were not asked for');
    };
    try {
        const res = { statusCode: 200, payload: null, setHeader() { return this; }, status(c) { this.statusCode = c; return this; }, json(v) { this.payload = v; return this; }, end() { return this; } };
        await clinvarProteinHandler({ method: 'GET', headers: { origin: 'https://www.drdoubleb.com' }, query: { gene: 'BRCA1', change: '9999delT', retmax: '200' } }, res);
        assert.equal(res.statusCode, 200);
        assert.deepEqual(res.payload, { variants: [], total: 0, notFound: ['9999delT'] });
        assert.equal(asked.length, 1);
        // The same for the protein-change search this proxy was written for.
        const same = { statusCode: 200, payload: null, setHeader() { return this; }, status(c) { this.statusCode = c; return this; }, json(v) { this.payload = v; return this; }, end() { return this; } };
        globalThis.fetch = async () => Response.json({ esearchresult: { count: '756', idlist: ['1', '2'], errorlist: { phrasesnotfound: ['G12K'] } } });
        await clinvarProteinHandler({ method: 'GET', headers: {}, query: { gene: 'KRAS', change: 'Gly12Lys,G12K' } }, same);
        assert.deepEqual(same.payload, { variants: [], total: 0, notFound: ['G12K'] });
    } finally {
        globalThis.fetch = realFetch;
    }
}
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
    const found = await h.resolveClinvarOtherName('BRCA1', 'c.5382insC', fetcher);
    assert.equal(found.hgvs, 'NM_007294.4:c.5266dup');
    assert.equal(found.basis, 'name');
    // ClinVar returns the newest record first, so a long-standing name's
    // variant comes last: everything is asked for, not the first 25.
    assert.equal(h.CLINVAR_NAME_SEARCH_SIZE, 200);
    assert.deepEqual(calls, [['BRCA1', ['5382insC'], 200]]);
    assert.equal(await h.resolveClinvarOtherName('BRCA1', '9999delT', fetcher), null, 'not listed');
    const warn = console.warn;
    console.warn = () => {};
    assert.equal(await h.resolveClinvarOtherName('BRCA1', '5382insC', async () => { throw new Error('502'); }), undefined, 'could not ask');
    console.warn = warn;
    assert.equal(await h.resolveClinvarOtherName('', '5382insC', fetcher), null);
}

// ── "ClinVar lists none" only on ClinVar's word ────────────────────────────
{
    const variants = searched.payload.variants;
    const ask = async (name, answer, gene = 'BRCA1') => {
        const told = {};
        const warn = console.warn;
        console.warn = () => {};
        const result = await h.resolveClinvarOtherName(gene, name, async () => (typeof answer === 'function' ? answer() : answer), told);
        console.warn = warn;
        return { result, told };
    };
    const others = variants.filter(v => v.id !== '17677');

    // Everything the search found is in hand, names and all, and none is it.
    let { result, told } = await ask('5382insG', { variants, total: 5 });
    assert.equal(result, null);
    assert.equal(told.outcome, 'none');
    assert.equal(h.otherNameOutcomeText(told, 'BRCA1', '5382insG'), 'no BRCA1 variant listed under "5382insG"');

    // ClinVar says it has never seen the text. What it returned instead - the
    // gene's newest records, 25 of 85,659 - does not come into it.
    ({ result, told } = await ask('9999delT', { variants: others, total: 85659, notFound: ['9999delT'] }));
    assert.equal(result, null);
    assert.equal(told.outcome, 'none');
    ({ result, told } = await ask('9999DELT', { variants: [], total: 0, notFound: ['9999delT'] }));
    assert.equal(told.outcome, 'none', 'whatever the case');

    // The list was cut short and the name is not in the part that came: that
    // settles nothing. (This is what "CFTR 6211GT" was: 25 of 9,631.)
    ({ result, told } = await ask('5382insC', { variants: others, total: 96 }));
    assert.equal(result, undefined, 'not "ClinVar lists none"');
    assert.equal(told.outcome, 'partial');
    assert.equal(h.otherNameOutcomeText(told, 'BRCA1', '5382insC'),
        'ClinVar returned 4 of 96 records for "5382insC" and none of those goes by it — not settled');
    // ...whereas a name found in a list cut short is found.
    ({ result, told } = await ask('5382insC', { variants, total: 96 }));
    assert.equal(result.id, '17677');
    assert.equal(told.outcome, 'found');

    // A backend deployed before it sent the names: nothing can be read off it.
    const bare = variants.map(({ aliases, variationCount, ...rest }) => rest);
    ({ result, told } = await ask('5382insC', { variants: bare, total: 5 }));
    assert.equal(result, undefined);
    assert.equal(told.outcome, 'no-names');
    assert.match(h.otherNameOutcomeText(told, 'BRCA1', '5382insC'), /does not pass on ClinVar's other names — could not be consulted/);

    // The proxy cannot be reached.
    ({ result, told } = await ask('5382insC', () => { throw new Error('502'); }));
    assert.equal(result, undefined);
    assert.equal(told.outcome, 'failed');
    assert.equal(h.otherNameOutcomeText(told, 'BRCA1', '5382insC'), 'could not be consulted');

    // Two variants go by the name: neither is chosen, both are named.
    const twin = { ...variants.find(v => v.id === '17677'), id: '999999', title: 'NM_007294.4(BRCA1):c.5267dup (p.Gln1756fs)' };
    ({ result, told } = await ask('5382insC', { variants: [...variants, twin], total: 6 }));
    assert.equal(result, null);
    assert.equal(told.outcome, 'several');
    assert.deepEqual(told.candidates.map(v => v.hgvs), ['NM_007294.4:c.5266dup', 'NM_007294.4:c.5267dup']);
    assert.equal(h.otherNameOutcomeText(told, 'BRCA1', '5382insC'),
        '2 BRCA1 variants are listed under "5382insC" — not for this page to choose between');
}
{
    // A protein change is a description several variants can share; a name is
    // one variant's. As ClinVar holds CFTR (trimmed to what is read):
    const cftr = [
        { id: '4072070', title: 'NM_000492.4(CFTR):c.1522_1524del (p.Phe508del)', aliases: [], variationCount: 1 },
        { id: '634837', title: 'NM_000492.3(CFTR):c.[1521_1523delCTT;3080T>C]', aliases: ['F508del', 'deltaF508'], variationCount: 2 },
        { id: '7105', title: 'NM_000492.3(CFTR):c.1521_1523del (p.Phe508del)', aliases: ['F508del', 'deltaF508', 'Phe508del', '[delta]F508'], variationCount: 1 }
    ];
    const ask = async (gene, name, variants, total = variants.length) => {
        const told = {};
        return { result: await h.resolveClinvarOtherName(gene, name, async () => ({ variants, total }), told), told };
    };
    // Two variants have the protein change; one of them is CALLED F508del.
    for (const name of ['F508del', 'deltaF508', 'p.Phe508del', 'Phe508del']) {
        const { result, told } = await ask('CFTR', name, cftr);
        assert.deepEqual(result, { id: '7105', title: 'NM_000492.3(CFTR):c.1521_1523del (p.Phe508del)', hgvs: 'NM_000492.3:c.1521_1523del', basis: 'name' }, name);
        assert.equal(told.outcome, 'found', name);
    }
    assert.deepEqual(h.clinvarVariantsDescribed(cftr, 'CFTR', 'F508del').map(v => v.id), ['4072070', '7105']);
    // The haplotype that contains it is not it, as before.
    assert.deepEqual(h.clinvarVariantsNamed(cftr, 'CFTR', 'F508del').map(v => v.id), ['7105']);

    // No record lists the name; three have it as their protein change.
    const egfr = [
        { id: '177620', title: 'NM_005228.5(EGFR):c.2236_2250del (p.Glu746_Ala750del)', aliases: [], variationCount: 1 },
        { id: '163343', title: 'NM_005228.5(EGFR):c.2235_2249del (p.Glu746_Ala750del)', aliases: [], variationCount: 1 },
        { id: '45233', title: 'NM_005228.5(EGFR):c.2232_2249delinsAAA (p.Glu746_Ala750del)', aliases: [], variationCount: 1 }
    ];
    for (const name of ['p.Glu746_Ala750del', 'E746_A750del']) {
        const { result, told } = await ask('EGFR', name, egfr);
        assert.equal(result, null, name);
        assert.equal(told.outcome, 'described', name);
        assert.deepEqual(told.candidates.map(v => v.hgvs), ['NM_005228.5:c.2236_2250del', 'NM_005228.5:c.2235_2249del', 'NM_005228.5:c.2232_2249delinsAAA']);
        assert.equal(h.otherNameOutcomeText(told, 'EGFR', name),
            'not a name ClinVar lists; 3 EGFR variants have this protein change — it does not say which');
    }
    // Exactly one variant with the protein change, and it is nobody's name: that one.
    let { result, told } = await ask('EGFR', 'p.Glu746_Ala750del', egfr.slice(0, 1));
    assert.deepEqual(result, { id: '177620', title: egfr[0].title, hgvs: 'NM_005228.5:c.2236_2250del', basis: 'protein' });
    assert.equal(told.outcome, 'found');
    // ...but not out of a list cut short, where there may be another.
    ({ result, told } = await ask('EGFR', 'p.Glu746_Ala750del', egfr.slice(0, 1), 40));
    assert.equal(result, undefined);
    assert.equal(told.outcome, 'partial');
    // A different protein change in the title is not this one.
    ({ result, told } = await ask('EGFR', 'p.Glu746_Ala750del', [{ ...egfr[0], title: 'NM_005228.5(EGFR):c.2237_2251del (p.Glu746_Thr751delinsAla)' }]));
    assert.equal(result, null);
    assert.equal(told.outcome, 'none');

    // How each is announced.
    assert.equal(h.otherNameNote({ name: 'F508del', title: 'T', id: '7105', basis: 'name' }), '"F508del" is another name for T (Variation 7105)');
    assert.equal(h.otherNameNote({ name: 'F508del', title: 'T', id: '7105' }, true),
        '"F508del" is another name ClinVar lists for T (Variation 7105), which is what is shown here.');
    assert.equal(h.otherNameNote({ name: 'p.X', title: 'T', id: '1', basis: 'protein' }, true),
        '"p.X" is the protein change of one variant in ClinVar, T (Variation 1), which is what is shown here.');
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

const nameSearches = page => page.api('clinvar-protein').filter(r => /retmax=200/.test(r));
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
    assert.deepEqual(nameSearches(page), ['clinvar-protein?gene=BRCA1&change=185delAG&retmax=200']);
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

// 4. A name ClinVar does not list — on ClinVar's word: it reports the text as
//    one it has never seen (as recorded).
{
    assert.deepEqual(cassettes.brca1UnknownName.api['clinvar-protein?gene=BRCA1&change=9999delT&retmax=200'],
        { variants: [], total: 0, notFound: ['9999delT'] });
    const page = await runLookup(cassettes.brca1UnknownName);
    assert.equal(page.status, 'Error: "BRCA1 9999delT" is not a name this tool can read. A variant is written in HGVS — '
        + 'BRCA1 c.…, NM_…:c.…, chr7:g.… — or as a coordinate row. ClinVar lists no BRCA1 variant under that name either.');
    assert.equal(page.hidden, true);
    assert.match(page.panel, /✕Parse input\s*not HGVS/);
    assert.match(page.panel, /◍Other names · ClinVar\s*no BRCA1 variant listed under "9999delT"/);
    assert.deepEqual(page.requests.filter(r => /myvariant|ensembl/.test(r)), [], 'no other service is asked about it');
    console.log('ok other names: an unknown name is refused, and says what would be read');
}

// 4b. The same sentence used to be said on no better ground than a list that
//     did not contain the name: the gene's first 25 records, of 85,659, as the
//     proxy returned them before it said which word ClinVar had dropped. A
//     list cut short says nothing about what ClinVar lists.
{
    const cutShort = { total: 85659, variants: cassettes.brca1OtherNames.api['clinvar-protein?gene=BRCA1&change=5382insC&retmax=200']
        .variants.filter(v => v.id !== '17677') };
    const page = await runLookup(cassettes.brca1UnknownName, { backend: { 'clinvar-protein': () => cutShort } });
    assert.equal(page.status, 'Error: "BRCA1 9999delT" is not a name this tool can read. A variant is written in HGVS — '
        + 'BRCA1 c.…, NM_…:c.…, chr7:g.… — or as a coordinate row.', 'no claim about ClinVar');
    assert.match(page.panel, /▲Other names · ClinVar\s*ClinVar returned 4 of 85659 records for "9999delT" and none of those goes by it — not settled/);
    console.log('ok other names: a list cut short is not "ClinVar lists none"');
}

// 4c. Two variants under one name: said, with both, and neither is chosen.
{
    const answer = structuredClone(cassettes.brca1OtherNames.api['clinvar-protein?gene=BRCA1&change=5382insC&retmax=200']);
    const founder = answer.variants.find(v => v.id === '17677');
    answer.variants.push({ ...founder, id: '999999', title: 'NM_007294.4(BRCA1):c.5267dup (p.Gln1756fs)' });
    answer.total += 1;
    const page = await runLookup(cassettes.brca1OtherNames, { backend: { 'clinvar-protein': () => answer } });
    assert.equal(page.status, 'Error: "BRCA1 5382insC" is not a name this tool can read. A variant is written in HGVS — '
        + 'BRCA1 c.…, NM_…:c.…, chr7:g.… — or as a coordinate row. ClinVar lists 2 BRCA1 variants under that name, so it does not say '
        + 'which is meant: NM_007294.4:c.5266dup, NM_007294.4:c.5267dup.');
    assert.match(page.panel, /◍Other names · ClinVar\s*2 BRCA1 variants are listed under "5382insC" — not for this page to choose between/);
    console.log('ok other names: two variants under one name are both named, neither chosen');
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
//    resolved, nothing is wrongly resolved — and nothing is said about what
//    ClinVar lists, which that answer cannot show. (It used to end in
//    "ClinVar lists no BRCA1 variant under that name either": for 5382insC.)
{
    const older = structuredClone(cassettes.brca1OtherNames.api['clinvar-protein?gene=BRCA1&change=5382insC&retmax=200']);
    for (const variant of older.variants) { delete variant.aliases; delete variant.variationCount; }
    delete older.notFound;
    const page = await runLookup(cassettes.brca1OtherNames, { backend: { 'clinvar-protein': () => older } });
    assert.equal(page.status, 'Error: "BRCA1 5382insC" is not a name this tool can read. A variant is written in HGVS — '
        + 'BRCA1 c.…, NM_…:c.…, chr7:g.… — or as a coordinate row.');
    assert.equal(page.hidden, true);
    assert.match(page.panel, /▲Other names · ClinVar\s*this backend does not pass on ClinVar's other names — could not be consulted/);
    console.log('ok other names: an older backend degrades to the plain error, with no claim about ClinVar');
}

// 7. A name that resolves the ordinary way never consults ClinVar's names.
{
    const page = await runLookup(cassettes.kras35);
    assertResolved(assert, page, 'KRAS c.35G>A');
    assert.deepEqual(nameSearches(page), []);
    assert.equal(page.document.querySelector('.other-name-note'), null);
    console.log('ok other names: not consulted for a name that resolves');
}

// 8. The lookup runs itself again for the variant ClinVar names — and must run
//    for what was SUBMITTED. It used to read the form a second time: with the
//    next variant already typed into the box while ClinVar was being asked, the
//    page showed BRCA1 c.5266dup and ClinVar 17677 beside
//    "VCF (hg19): 7-140453136-A-T", and sent that position to gnomAD and SpliceAI.
{
    const page = await runLookup(cassettes.brca1OtherNames, {
        route: (url, key, respond, options, { document }) => {
            // The user types on while the names are being looked up.
            if (/\/api\/clinvar-protein/.test(key) && /retmax=200/.test(key)) {
                document.getElementById('variantInput').value = '7 140453136 A T';
                document.getElementById('tumorTypeInput').value = 'melanoma';
            }
            return undefined;
        }
    });
    assertResolved(assert, page, 'BRCA1 5382insC, box edited meanwhile');
    assert.match(page.variant, /g\.: chr17:g\.41209082dup/);
    assert.match(page.variant, /VCF \(hg19\): 17-41209079-T-TG/);
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.match(page.panel, /You entered\s*BRCA1 5382insC/);
    assert.doesNotMatch(everywhere(page), /140453136/, 'nothing about the variant that was only typed');
    assert.doesNotMatch(page.requests.join('\n'), /melanoma/i, 'nor the tumour type that was only typed');
    assert.equal(nameSearches(page).length, 1);
    console.log('ok other names: the second pass is for what was submitted, not for what is in the box by then');
}

// 9. The second pass failing ends the matter: ClinVar's names are not asked
//    again, and the error is the second pass's own.
{
    const page = await runLookup(cassettes.brca1OtherNames, { without: [/ensembl\.org/, /myvariant\.info/, /api\.ncbi\.nlm\.nih\.gov\/variation/] });
    assert.match(page.status, /^Error: Lookup could not complete: .*did not respond/);
    assert.equal(nameSearches(page).length, 1, 'asked once, on the first pass');
    assert.match(page.panel, /Other names · ClinVar"5382insC" is another name for NM_007294\.4\(BRCA1\):c\.5266dup/);
    // The same when the second pass fails for a reason of its own rather than
    // an outage — the case in which the first pass DOES go to ClinVar's names.
    let asks = 0;
    const rejected = await runLookup(cassettes.brca1OtherNames, {
        route: (url, key, respond) => {
            // (Were the second pass to ask again, each answer would start a third:
            // cut that short so this fails rather than hangs.)
            if (/\/api\/clinvar-protein/.test(key) && /retmax=200/.test(key) && ++asks > 3) return respond({ error: 'asked again' }, 500);
            // Every service answers — that it has nothing, or will not read the name.
            if (/variant_recoder|\/vep\/|api\.ncbi\.nlm\.nih\.gov\/variation/.test(key)) return respond({ error: 'Unable to parse HGVS notation' }, 400);
            if (/myvariant\.info\/v1\/query/.test(key)) return respond({ total: 0, hits: [] });
            if (/myvariant\.info\/v1\/variant/.test(key)) return respond({ success: false, error: 'Not Found.' }, 404);
            return undefined;
        }
    });
    assert.match(rejected.status, /^Error: /);
    assert.doesNotMatch(rejected.status, /did not respond/, 'not an outage this time');
    assert.equal(nameSearches(rejected).length, 1, 'asked once: the second pass does not ask again');
    console.log('ok other names: a failed second pass is reported, not retried');
}

// 10. A name written with marks. "621+1G>T" went to ClinVar as "6211GT", which
//     it has never seen; the page then said ClinVar listed no such CFTR variant.
{
    const page = await runLookup(cassettes.cftr621);
    assertResolved(assert, page, 'CFTR 621+1G>T');
    assert.deepEqual(nameSearches(page), ['clinvar-protein?gene=CFTR&change=621+1G>T&retmax=200'], 'searched as written');
    assert.match(page.variant, /g\.: chr7:g\.117171169G>T/);
    assert.match(page.variant, /c\.: c\.489\+1G>T/);
    assert.match(page.clinvar, /Variation ID: 38799/);
    assert(page.document.querySelector('.other-name-note').textContent.includes(
        '"621+1G>T" is another name ClinVar lists for NM_000492.4(CFTR):c.489+1G>T (Variation 38799)'));
    assert(!page.requests.some(r => /6211GT/.test(r)), 'nothing is asked about the name with its marks removed');
    assert(!page.requests.some(r => /621\+1G>T/.test(r) && !/clinvar-protein/.test(r)), 'and Ensembl and MyVariant are not asked about a name that is not HGVS');
    console.log('ok other names: CFTR 621+1G>T is c.489+1G>T');
}

// 11. The old way of naming an intronic change by its intron. "IVS8-5T" was
//     read as a protein change — residue 8 — and that codon was shown.
{
    const nothing = { variants: [], total: 0, notFound: ['IVS8-5T'] };
    const page = await runLookup(cassettes.cftr621, { input: 'CFTR IVS8-5T', backend: { 'clinvar-protein': () => nothing } });
    assert.equal(page.status, 'Error: "CFTR IVS8-5T" is not a name this tool can read. A variant is written in HGVS — '
        + 'CFTR c.…, NM_…:c.…, chr7:g.… — or as a coordinate row. ClinVar lists no CFTR variant under that name either.');
    assert.equal(page.hidden, true);
    assert.deepEqual(nameSearches(page), ['clinvar-protein?gene=CFTR&change=IVS8-5T&retmax=200']);
    assert.deepEqual(page.requests.filter(r => /myvariant|ensembl/.test(r)), [], 'not sent anywhere as "p.IVS8-5T"');
    // ...and one ClinVar does list resolves like any other name.
    const listed = { total: 1, notFound: [], variants: [{ id: '38799', title: 'NM_000492.4(CFTR):c.489+1G>T', aliases: ['IVS4+1G>T', '621+1G>T'], variationCount: 1 }] };
    const known = await runLookup(cassettes.cftr621, { input: 'CFTR IVS4+1G>T', backend: { 'clinvar-protein': () => listed } });
    assertResolved(assert, known, 'CFTR IVS4+1G>T');
    assert.match(known.variant, /c\.: c\.489\+1G>T/);
    console.log('ok other names: an IVS name is a name, not a protein change');
}

// 12. "CFTR F508del". Ensembl's allele for p.Phe508del is c.1522_1524del (TTT
//     removed; rs121909001, which MyVariant holds — so it was shown). The
//     variant everyone means is c.1521_1523del (CTT): another allele, the same
//     protein. ClinVar lists "F508del" as a name of that one only — in the
//     last of the 96 records its search returns, newest first.
{
    const answer = cassettes.cftrF508delName.api['clinvar-protein?gene=CFTR&change=F508del&retmax=200'];
    assert.equal(answer.variants.length, 96);
    assert(answer.variants.findIndex(v => v.id === '7105') > 25, 'well past the 25 that used to be asked for');
    const page = await runLookup(cassettes.cftrF508delName);
    assertResolved(assert, page, 'CFTR F508del');
    assert.match(page.variant, /c\.: c\.1521_1523del/);
    assert.match(page.variant, /p\.: p\.Phe508del/);
    assert.match(page.variant, /VCF \(hg19\): 7-117199644-ATCT-A/);
    assert.match(page.clinvar, /Variation ID: 7105/);
    assert.doesNotMatch(page.variant, /1522_1524|117199647_117199649/);
    assert(page.document.querySelector('.other-name-note').textContent.includes(
        '"F508del" is another name ClinVar lists for NM_000492.3(CFTR):c.1521_1523del (p.Phe508del) (Variation 7105)'));
    assert.equal(nameSearches(page).length, 1);
    console.log('ok other names: CFTR F508del is c.1521_1523del, not the other allele with that protein change');
}

// 13. Where the name is nobody's and several variants have the protein change,
//     Ensembl's allele is still what is shown — and the card says there are
//     others, with each a click away. (Here: ClinVar's CFTR answer with the
//     names taken off the founder variant, which leaves two that are described
//     as p.Phe508del and none that is called F508del.)
{
    const answer = structuredClone(cassettes.cftrF508delName.api['clinvar-protein?gene=CFTR&change=F508del&retmax=200']);
    for (const variant of answer.variants) if (variant.variationCount === 1) variant.aliases = [];
    const page = await runLookup(cassettes.cftrF508delName, { backend: { 'clinvar-protein': () => answer } });
    assertResolved(assert, page, 'CFTR F508del, a name of no variant');
    assert.match(page.variant, /c\.: c\.1522_1524del/, 'the allele Ensembl gives');
    assert.match(page.panel, /Other names · ClinVar\s*not a name ClinVar lists; 2 CFTR variants have this protein change — it does not say which/);
    const note = page.document.querySelector('.protein-change-shared-note');
    assert.equal(note.textContent, '2 CFTR variants in ClinVar have the protein change "F508del", none of them the allele Ensembl gives for it, '
        + 'which is what is shown: NM_000492.4:c.1522_1524del, NM_000492.3:c.1521_1523del.');
    assert.deepEqual([...note.querySelectorAll('a')].map(a => decodeURIComponent(a.getAttribute('href'))),
        ['?variant=NM_000492.4:c.1522_1524del', '?variant=NM_000492.3:c.1521_1523del']);
    assert.equal(page.document.querySelector('.other-name-note'), null);
    console.log('ok other names: a protein change several variants share is said to be shared');
}

// 13b. The allele in hand can be the very one ClinVar knows by the name. Then
//      there is nothing to look up a second time, and nothing to announce as
//      a different variant. (ClinVar's answer with the name moved to the
//      c.1522_1524del record, and MyVariant's record of that allele carrying
//      that ClinVar id.)
{
    const cassette = cassettes.cftrF508delName;
    const answer = structuredClone(cassette.api['clinvar-protein?gene=CFTR&change=F508del&retmax=200']);
    for (const variant of answer.variants) {
        if (variant.variationCount === 1) variant.aliases = variant.id === '4072070' ? ['F508del'] : [];
    }
    const twinKey = 'myvariant.info/v1/variant/chr7:g.117199647_117199649del';
    const twin = { ...cassette.responses[twinKey].body, clinvar: { variant_id: 4072070, gene: { symbol: 'CFTR' } } };
    const page = await runLookup(cassette, {
        backend: { 'clinvar-protein': () => answer },
        route: (url, key, respond) => (key === twinKey ? respond(twin) : undefined)
    });
    assertResolved(assert, page, 'CFTR F508del, the allele in hand is the named one');
    assert.match(page.variant, /c\.: c\.1522_1524del/);
    assert.match(page.panel, /●Other names · ClinVar\s*ClinVar knows this variant as "F508del"/);
    assert.equal(page.document.querySelector('.other-name-note'), null, 'not announced as another variant');
    assert(!page.requests.some(r => /variant_recoder\/human\/NM_000492/.test(r)), 'and not looked up a second time');
    assert.equal(nameSearches(page).length, 1);
    console.log('ok other names: the variant in hand being the named one, nothing is run again');
}

// 14. A substitution MyVariant holds is what a substitution's name means:
//     ClinVar's names are not consulted for it.
{
    for (const cassette of [cassettes.brafV600E, cassettes.tp53Stop, cassettes.krasQ61X]) {
        const page = await runLookup(cassette);
        assertResolved(assert, page, cassette.input);
        assert.deepEqual(nameSearches(page), [], cassette.input);
    }
    console.log('ok other names: not consulted for a substitution on record');
}

console.log('Other-name tests passed');
