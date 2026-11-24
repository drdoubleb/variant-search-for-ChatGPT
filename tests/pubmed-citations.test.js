import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import handler, { parsePubmedRecords } from '../api/pubmed.js';
import { captureResponse } from './support/api-response.js';

const xml = readFileSync(new URL('./fixtures/pubmed-tp53-functional.xml', import.meta.url), 'utf8');
const records = parsePubmedRecords(xml);
assert.deepEqual([...records.keys()].sort(), ['12826609', '29979965', '30224644']);
assert.equal(records.get('12826609').year, '2003');
assert.match(records.get('12826609').authors, /^Kato S/);
assert.match(records.get('29979965').journal, /Mol Cell/);
assert.match(records.get('30224644').title, /Mutational processes/);
assert.ok([...records.values()].every(article => article.abstract));
const structured = parsePubmedRecords(`<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>1</PMID><Article>
<ArticleTitle>TP53 &amp; <i>p53</i></ArticleTitle><Abstract><AbstractText Label="RESULTS">5 &lt; 7 and <b>valid</b>.</AbstractText></Abstract>
</Article></MedlineCitation></PubmedArticle></PubmedArticleSet>`);
assert.equal(structured.get('1').title, 'TP53 & p53');
assert.equal(structured.get('1').abstract, 'RESULTS: 5 < 7 and valid.');
assert.throws(() => parsePubmedRecords('<!DOCTYPE PubmedArticleSet [<!ENTITY bad SYSTEM "https://example.invalid/ignored-entity">]><PubmedArticleSet/>'), /document type/);

const realFetch = globalThis.fetch;
try {
    const requests = [];
    let response = new Response(xml);
    globalThis.fetch = async url => { requests.push(new URL(url)); return response.clone(); };
    const call = query => handler({ method: 'GET', headers: {}, query }, captureResponse());
    for (const pmids of ['', '0', '-1', '01234', '123,x', '123,', '1&retmax=500', '1000000000', ['1', '2'],
        Array.from({ length: 21 }, (_, i) => String(i + 1)).join(',')]) {
        assert.equal((await call({ pmids })).statusCode, 400);
    }
    assert.equal((await call({ pmids: '1', term: 'TP53' })).statusCode, 400);
    assert.equal(requests.length, 0);
    const result = await call({ pmids: '29979965,12826609,29979965,9999,30224644' });
    assert.equal(result.statusCode, 200);
    assert.equal(requests.length, 1);
    assert.equal(requests[0].searchParams.get('id'), '29979965,12826609,9999,30224644');
    assert.equal(result.body.backend, 'pmid');
    assert.deepEqual(result.body.articles.map(article => article.pmid), ['29979965', '12826609', '30224644']);
    assert.deepEqual(result.body.missingPmids, ['9999']);
    assert.match(result.headers['Cache-Control'], /no-store/);
    response = new Response('<PubmedArticleSet><ERROR>Invalid request</ERROR></PubmedArticleSet>');
    const malformed = await call({ pmids: '12826609' });
    assert.equal(malformed.statusCode, 502);
    assert.match(malformed.headers['Cache-Control'], /no-store/);
    response = new Response('', { status: 404 });
    assert.equal((await call({ pmids: '12826609' })).statusCode, 502);
} finally {
    globalThis.fetch = realFetch;
}
console.log('PubMed explicit PMID tests passed');
