/*
 * A frameshift is not the missense at its first residue.
 *
 * p.Gln1756ProfsTer74 opens the way a substitution does, and the code that
 * shortened protein changes for searching read only that far. BRCA1 c.5266dup
 * (5382insC) was looked up in ClinVar's "same protein change" search, in
 * LitVar/PubMed and in the Google links as "BRCA1 Q1756P" — a different variant
 * with its own ClinVar records and its own paper. A frameshift MyVariant writes
 * without the new residue ("p.Pro153fs") came out worse: "TP53 O153F".
 *
 * The scenarios replay answers recorded from the live services
 * (tests/fixtures/lookup-cassettes.json) through the real form handler.
 *
 * Run with: node tests/protein-search-term.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';
import pubmedHandler, { variantPhraseDropped } from '../api/pubmed.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));

// ── Only a substitution is a substitution ──────────────────────────────────
for (const text of [
    'p.Gln1756ProfsTer74', 'p.Gln1756Profs*74', 'p.Q1756Pfs*74', 'p.Pro153fs', 'p.Q1756fs',
    'NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs)',
    'p.Lys18del', 'p.K18del', 'p.Glu746_Ala750del', 'p.Gly12delinsAlaSer', 'p.Ter394GlnextTer12', 'p.Val600=', 'p.Met1?'
]) {
    assert.equal(h.parseProteinChange(text), null, `${text} is not an amino-acid substitution`);
}
assert.deepEqual(h.parseProteinChange('p.Gly34Arg'), { ref: 'G', pos: 34, alt: 'R' });
assert.deepEqual(h.parseProteinChange('p.(Gly34Arg)'), { ref: 'G', pos: 34, alt: 'R' });
assert.deepEqual(h.parseProteinChange('NM_001904.4(CTNNB1):c.100G>C (p.Gly34Arg)'), { ref: 'G', pos: 34, alt: 'R' });
assert.deepEqual(h.parseProteinChange('GLY34ARG'), { ref: 'G', pos: 34, alt: 'R' });
assert.deepEqual(h.parseProteinChange('G34R'), { ref: 'G', pos: 34, alt: 'R' });
// A stop, in each of the ways it is written.
for (const text of ['p.Arg213Ter', 'p.Arg213*', 'p.R213*', 'R213*', 'NM_000546.6(TP53):c.637C>T (p.Arg213Ter)']) {
    assert.deepEqual(h.parseProteinChange(text), { ref: 'R', pos: 213, alt: '*' }, text);
}

// ── What a protein change is searched as ───────────────────────────────────
assert.equal(h.proteinSearchTerm('p.Val600Glu'), 'V600E');
assert.equal(h.proteinSearchTerm('VAL600GLU'), 'V600E', 'the form the query target is kept in');
assert.equal(h.proteinSearchTerm('NP_004324.2:p.Val600Glu'), 'V600E');
for (const [text, term] of [
    ['p.Gln1756ProfsTer74', 'Q1756fs'], ['NP_009225.1:p.Gln1756ProfsTer74', 'Q1756fs'], ['p.Gln1756Profs*74', 'Q1756fs'],
    ['p.Q1756Pfs*74', 'Q1756fs'], ['p.Gln1756fs', 'Q1756fs'], ['p.Gln1756ProfsTer?', 'Q1756fs'],
    ['p.Pro153fs', 'P153fs'], ['p.Pro153AlafsTer28', 'P153fs'], ['p.Asp362GlnfsTer6', 'D362fs']
]) {
    assert.equal(h.proteinSearchTerm(text), term, `${text} is searched as ${term}`);
}
// A stop goes by the one-letter name all three searches know it under.
for (const text of ['p.Arg213Ter', 'p.Arg213*', 'p.R213*', 'ARG213TER', 'NP_000537.3:p.Arg213Ter']) {
    assert.equal(h.proteinSearchTerm(text), 'R213X', `${text} is searched as R213X`);
}
// Not shortened: these are searched by their full p. string, as before.
for (const text of ['p.Glu746_Ala750del', 'p.Lys18del', 'p.K18del', 'p.Val600=', 'p.Met1?', 'p.Ter394Gln', '', null]) {
    assert.equal(h.proteinSearchTerm(text), '', `${text} has no short search name`);
}

// ── A substitution typed by hand ───────────────────────────────────────────
for (const typed of ['V600E', 'v600e', 'val600glu', 'VAL600GLU', 'Val600Glu']) {
    assert.equal(h.normalizeProteinSubstitution(typed), 'Val600Glu', typed);
}
for (const typed of ['R213*', 'R213X', 'r213x', 'Arg213*', 'R213Ter', 'ARG213TER', 'Arg213Ter']) {
    assert.equal(h.normalizeProteinSubstitution(typed), 'Arg213Ter', `${typed} is a stop at 213`);
}
for (const typed of ['Gln1756ProfsTer74', 'Q1756Pfs*74', 'Glu746_Ala750del', 'K18del', '185delAG', 'V600=', 'M1?', 'B12Z']) {
    assert.equal(h.normalizeProteinSubstitution(typed), typed, `${typed} is left as typed`);
}
console.log('Protein search-term unit checks passed');

// ── PubMed drops what it does not know ─────────────────────────────────────
assert.equal(variantPhraseDropped('D362fs', ['D362fs']), true);
assert.equal(variantPhraseDropped('p.Arg213Ter', ['p.Arg213Ter']), true);
assert.equal(variantPhraseDropped('Q1756Pfs*74', ['Q1756Pfs*74']), true);
assert.equal(variantPhraseDropped('c.455dup', ['c.455dup']), true);
assert.equal(variantPhraseDropped('p.Glu746_Ala750del', ['p.Glu746']), false, 'the rest of the variant was still searched');
assert.equal(variantPhraseDropped('V600E', []), false);
assert.equal(variantPhraseDropped('V600E', undefined), false);
assert.equal(variantPhraseDropped('V600E', ['hgsoc']), false, 'some other word was dropped, not the variant');
assert.equal(variantPhraseDropped('', ['D362fs']), false, 'no variant was asked about');
// What says which variant it is, is the position. With that dropped, PubMed
// goes on to search for the base that is left: PTEN AND "G" (1,491 papers),
// CFTR AND "T" (8,572); "A" it ignores altogether (CFTR alone: 15,803).
assert.equal(variantPhraseDropped('c.1026+59T>G', ['c.1026+59T']), true);
assert.equal(variantPhraseDropped('c.489+1G>T', ['c.489+1G']), true);
assert.equal(variantPhraseDropped('c.1585-1G>A', ['c.1585-1G']), true);
assert.equal(variantPhraseDropped('c.5266_5277del', ['c.5266_5277']), true, 'a bare "del" is not the variant either');
assert.equal(variantPhraseDropped('c.1083_1084insCA', ['1084insCA']), false, 'a position is still being searched');
assert.equal(variantPhraseDropped('c.35G>A', []), false, 'PubMed knew it as a phrase');
// A variant that never had a position is judged as before.
assert.equal(variantPhraseDropped('exon 14 skipping', ['skipping']), false);
assert.equal(variantPhraseDropped('amplification', ['amplification']), true);
assert.equal(variantPhraseDropped('MET amplification', ['metamp']), false);

function fakeResponse() {
    return {
        statusCode: 200, payload: null, headers: {},
        setHeader(name, value) { this.headers[name] = value; return this; },
        status(code) { this.statusCode = code; return this; },
        json(value) { this.payload = value; return this; },
        end() { return this; }
    };
}
async function pubmed(query, upstream) {
    const realFetch = globalThis.fetch;
    const asked = [];
    globalThis.fetch = async (address) => {
        const url = String(address);
        asked.push(url);
        const body = upstream(url);
        return new Response(typeof body === 'string' ? body : JSON.stringify(body), { status: 200 });
    };
    try {
        const res = fakeResponse();
        await pubmedHandler({ method: 'GET', headers: { origin: 'https://www.drdoubleb.com' }, query }, res);
        return { status: res.statusCode, payload: res.payload, asked };
    } finally {
        globalThis.fetch = realFetch;
    }
}
{
    // As recorded: LitVar has no entity for this rare frameshift, and PubMed,
    // not knowing "D362fs", answers with all 17,385 BRCA2 papers.
    const out = await pubmed({ term: 'BRCA2 D362fs', gene: 'BRCA2', variant: 'D362fs', limit: '5' }, (url) => {
        if (url.includes('/entity/autocomplete/')) return [];
        if (url.includes('esearch.fcgi')) {
            return { esearchresult: { count: '17385', idlist: ['39999999', '39999998', '39999997'],
                errorlist: { phrasesnotfound: ['D362fs'], fieldsnotfound: [] } } };
        }
        throw new Error(`unexpected upstream request ${url}`);
    });
    assert.equal(out.status, 200);
    assert.deepEqual(out.payload, { total: 0, articles: [], backend: 'term', notFound: ['D362fs'] });
    assert(!out.asked.some(u => /esummary|efetch/.test(u)), 'the gene-only papers are not fetched');
}
{
    // As recorded: PubMed does not know "c.1026+59T", and answers for PTEN and
    // the letter G.
    const out = await pubmed({ term: 'PTEN c.1026+59T>G', gene: 'PTEN', variant: 'c.1026+59T>G', limit: '5' }, (url) => {
        if (url.includes('/entity/autocomplete/')) return [];
        if (url.includes('esearch.fcgi')) {
            return { esearchresult: { count: '1491', idlist: ['39000001', '39000002'], querytranslation: '"PTEN"[All Fields] AND "G"[All Fields]',
                errorlist: { phrasesnotfound: ['c.1026+59T'], fieldsnotfound: [] } } };
        }
        throw new Error(`unexpected upstream request ${url}`);
    });
    assert.deepEqual(out.payload, { total: 0, articles: [], backend: 'term', notFound: ['c.1026+59T'] });
    assert(!out.asked.some(u => /esummary|efetch/.test(u)), 'the papers on PTEN and "G" are not fetched');
}
{
    // A variant PubMed does know is answered as before.
    const out = await pubmed({ term: 'KRAS G12C', gene: 'KRAS', variant: 'G12C', limit: '5' }, (url) => {
        if (url.includes('/entity/autocomplete/')) return [];
        if (url.includes('esearch.fcgi')) return { esearchresult: { count: '2109', idlist: ['31189530'] } };
        if (url.includes('esummary.fcgi')) {
            return { result: { 31189530: { title: 'The KRAS(G12C) inhibitor', source: 'Nature', pubdate: '2019 Nov', authors: [{ name: 'Canon J' }] } } };
        }
        if (url.includes('efetch.fcgi')) return '<PubmedArticle><PMID Version="1">31189530</PMID><AbstractText>Abstract.</AbstractText></PubmedArticle>';
        throw new Error(`unexpected upstream request ${url}`);
    });
    assert.equal(out.payload.total, 2109);
    assert.equal(out.payload.articles.length, 1);
    assert.equal(out.payload.articles[0].pmid, '31189530');
    assert.equal(out.payload.notFound, undefined);
}
{
    // Part of the variant surviving is still a search for it.
    const out = await pubmed({ term: 'EGFR p.Glu746_Ala750del', gene: 'EGFR', variant: 'p.Glu746_Ala750del', limit: '5' }, (url) => {
        if (url.includes('/entity/autocomplete/')) return [];
        if (url.includes('esearch.fcgi')) {
            return { esearchresult: { count: '7', idlist: ['34824841'], errorlist: { phrasesnotfound: ['p.Glu746'], fieldsnotfound: [] } } };
        }
        if (url.includes('esummary.fcgi')) return { result: { 34824841: { title: 'EGFR exon 19 deletions', source: 'J', pubdate: '2021' } } };
        if (url.includes('efetch.fcgi')) return '';
        throw new Error(`unexpected upstream request ${url}`);
    });
    assert.equal(out.payload.total, 7);
    assert.equal(out.payload.articles.length, 1);
}
{
    // A term-only search (the tumor-type tab, gene-level mode) names no variant
    // and is left exactly as it was.
    const out = await pubmed({ term: 'BRAF hgsoc', limit: '5' }, (url) => {
        if (url.includes('esearch.fcgi')) {
            return { esearchresult: { count: '900', idlist: ['1'], errorlist: { phrasesnotfound: ['hgsoc'], fieldsnotfound: [] } } };
        }
        if (url.includes('esummary.fcgi')) return { result: { 1: { title: 'A BRAF paper', source: 'J', pubdate: '2020' } } };
        if (url.includes('efetch.fcgi')) return '';
        throw new Error(`unexpected upstream request ${url}`);
    });
    assert.equal(out.payload.total, 900);
}
console.log('ok pubmed proxy: a variant PubMed has never indexed is not answered with the gene\'s papers');

// ── Full lookups ───────────────────────────────────────────────────────────
// The ClinVar card's main link: the record when there is one, else a search.
const clinvarLink = page => decodeURIComponent(page.links.find(l => l.text === 'View on ClinVar')?.href || '');
// Every search the page starts or links to for the variant, in one string.
const everySearch = page => [
    ...page.api('pubmed'), ...page.api('clinvar-protein'),
    ...page.searches('google.com/search'), ...page.searches('ncbi.nlm.nih.gov/clinvar/?term='),
    ...page.searches('pubmed.ncbi.nlm.nih.gov')
].join('\n');

// 1. BRCA1 c.5266dup (5382insC), p.Gln1756ProfsTer74. Was searched as Q1756P.
{
    const page = await runLookup(cassettes.brca1Dup);
    assertResolved(assert, page, 'BRCA1 c.5266dup');
    assert.match(page.variant, /p\.: p\.Gln1756ProfsTer74 \(p\.Q1756Pfs\*74\)/);
    assert.match(page.variant, /frameshift_variant/);
    assert.deepEqual(page.api('clinvar-protein'), [], 'no "same protein change" search: a frameshift is not a substitution');
    assert.doesNotMatch(page.clinvar, /Same protein change/);
    const pubmedAsked = page.api('pubmed');
    assert(pubmedAsked.length > 0, 'the PubMed card searched');
    for (const asked of pubmedAsked) {
        assert.match(asked, /term=BRCA1 Q1756fs(&|$)/);
        assert.match(asked, /gene=BRCA1&variant=Q1756fs(&|$)/);
    }
    assert.deepEqual(page.searches('google.com/search'),
        ['q=pathogenicity of BRCA1 Q1756fs', 'q=clinical significance of BRCA1 Q1756fs']);
    assert.doesNotMatch(everySearch(page), /Q1756P/, 'the missense Q1756P is never searched for');
    console.log('ok frameshift search: BRCA1 c.5266dup is searched as Q1756fs');
}

// 2. A frameshift ClinVar has no record of: the ClinVar link is a search, and
//    it used to be for the missense too.
{
    const page = await runLookup(cassettes.brca2InsCA);
    assertResolved(assert, page, 'BRCA2 insCA');
    assert.match(page.variant, /p\.: p\.Asp362GlnfsTer6/);
    assert.equal(clinvarLink(page), 'https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA2 D362fs');
    assert.deepEqual(page.api('clinvar-protein'), []);
    for (const asked of page.api('pubmed')) assert.match(asked, /gene=BRCA2&variant=D362fs(&|$)/);
    assert.doesNotMatch(everySearch(page), /D362[A-Z](?!fs)/, 'no missense at 362 is searched for');
    console.log('ok frameshift search: BRCA2 c.1083_1084insCA is searched as D362fs');
}

// 3. The same, when PubMed has nothing: the card says so instead of listing
//    the gene's most relevant papers under "Gene + Variant".
{
    const page = await runLookup(cassettes.brca2InsCA, {
        backend: { pubmed: () => ({ total: 0, articles: [], backend: 'term', notFound: ['D362fs'] }) }
    });
    assertResolved(assert, page, 'BRCA2 insCA, PubMed empty');
    assert.match(page.cards.PubMed, /No PubMed results found — no article in PubMed uses "D362fs"\./);
    console.log('ok frameshift search: an unknown term is named in the PubMed card');
}

// 4. MyVariant's own frameshift notation carries no new residue. With Ensembl
//    not answering that is all there is, and it became "TP53 O153F".
{
    const page = await runLookup(cassettes.tp53Dup, { without: [/ensembl\.org/] });
    assertResolved(assert, page, 'TP53 c.455dup, MyVariant only');
    assert.match(page.variant, /p\.: p\.Pro153fs/);
    assert.equal(clinvarLink(page), 'https://www.ncbi.nlm.nih.gov/clinvar/?term=TP53 P153fs');
    for (const asked of page.api('pubmed')) assert.match(asked, /gene=TP53&variant=P153fs(&|$)/);
    assert.doesNotMatch(everySearch(page), /O153F/);
    assert.doesNotMatch(everySearch(page), /P153[A-Z](?!fs)/);
    console.log('ok frameshift search: p.Pro153fs is searched as P153fs');
}

// 5. ...and with Ensembl answering, the fuller name gives the same search.
{
    const page = await runLookup(cassettes.tp53Dup);
    assertResolved(assert, page, 'TP53 c.455dup');
    assert.match(page.variant, /p\.: p\.Pro153AlafsTer28/);
    for (const asked of page.api('pubmed')) assert.match(asked, /gene=TP53&variant=P153fs(&|$)/);
    assert.doesNotMatch(everySearch(page), /P153A(?!fs)/);
    console.log('ok frameshift search: p.Pro153AlafsTer28 is searched as P153fs');
}

// 6. A stop typed the short way. "TP53 R213*" was sent on as "p.R213*", which
//    nothing indexes; the page came back with no coordinates at all.
for (const typed of ['TP53 R213*', 'TP53 R213X', 'tp53 r213x', 'TP53:p.R213*', 'TP53 p.Arg213*']) {
    const page = await runLookup(cassettes.tp53Stop, { input: typed });
    assertResolved(assert, page, typed);
    assert.match(page.panel, /read as TP53:p\.Arg213Ter/, typed);
    assert.match(page.variant, /g\.: chr17:g\.7578212G>A/, typed);
    assert.match(page.variant, /c\.: c\.637C>T/, typed);
    assert.match(page.variant, /p\.: p\.Arg213\* \(p\.R213\*\)/, typed);
    assert.match(page.variant, /stop_gained/, typed);
    // A stop IS a substitution: the same-protein search still runs for it.
    assert.deepEqual(page.api('clinvar-protein'), ['clinvar-protein?gene=TP53&change=Arg213Ter,R213*'], typed);
    // snpEff writes it "p.Arg213*", which LitVar cannot read and ClinVar's
    // search takes for a wildcard; it is searched under the name they know.
    for (const asked of page.api('pubmed')) assert.match(asked, /term=TP53 R213X&gene=TP53&variant=R213X(&|$)/, typed);
    assert.deepEqual(page.searches('google.com/search').slice(-2),
        ['q=pathogenicity of TP53 R213X', 'q=clinical significance of TP53 R213X'], typed);
}
console.log('ok nonsense input: R213*, R213X and Arg213* all resolve as p.Arg213Ter');

// 7. A missense is searched exactly as it always was.
{
    const page = await runLookup(cassettes.brafV600E);
    assertResolved(assert, page, 'BRAF V600E');
    assert.match(page.variant, /p\.: p\.Val600Glu \(p\.V600E\)/);
    assert.deepEqual(page.api('clinvar-protein'), ['clinvar-protein?gene=BRAF&change=Val600Glu,V600E']);
    for (const asked of page.api('pubmed')) {
        assert.match(asked, /term=BRAF V600E(&|$)/);
        assert.match(asked, /gene=BRAF&variant=V600E(&|$)/);
    }
    assert.deepEqual(page.searches('google.com/search'),
        ['q=pathogenicity of BRAF V600E', 'q=clinical significance of BRAF V600E']);
    console.log('ok missense search: BRAF V600E unchanged');
}

console.log('Protein search-term tests passed');
