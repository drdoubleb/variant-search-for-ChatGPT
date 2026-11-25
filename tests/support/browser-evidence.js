import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import clinvar from '../../api/clinvar-variant.js';
import pubmed, { parsePubmedRecords } from '../../api/pubmed.js';
import tp53 from '../../api/tp53.js';
import aiReview from '../../api/ai-review.js';
import { captureResponse } from './api-response.js';
import { checkM237lRelease, checkReleaseAssets } from './m237l-release.js';

const fixture = name => readFileSync(new URL(`../fixtures/${name}`, import.meta.url), 'utf8');
const m237lAnnotation = JSON.parse(fixture('myvariant-m237l.json'));
const r175hAnnotation = JSON.parse(fixture('lookup-cassettes.json')).rs28934578.responses['myvariant.info/v1/variant/chr17:g.7578406C>T'].body;
const clinvarXml = fixture('clinvar-2774649.xml');
const markedClinvarXml = clinvarXml.replace('</Comment>', ' &lt;strong&gt;literal&lt;/strong&gt;</Comment>');
const pubmedXml = fixture('pubmed-tp53-functional.xml');
const pubmedRecords = parsePubmedRecords(pubmedXml);
const tp53Csv = fixture('tp53-mutationview-r21.csv');
const tp53Map = JSON.parse(fixture('ensembl-gene-maps.json')).TP53;
const datasetUrl = 'https://tp53.cancer.gov/static/data/MutationView_r21.csv';
const transcript = { seq_id: 'NM_000546.6', position: 850, deleted_sequence: 'A', inserted_sequence: 'T' };
const grch37 = { seq_id: 'NC_000017.10', position: 7577571, deleted_sequence: 'T', inserted_sequence: 'A' };
const grch38 = { seq_id: 'NC_000017.11', position: 7674253, deleted_sequence: 'T', inserted_sequence: 'A' };
const studyPmids = ['12826609', '29979965', '30224644'];

function promptContext(prompt) {
    return JSON.parse(prompt.split('Context JSON:\n\n')[1].split('\n\nEnd of supplied context.')[0]);
}

function pageEvidence(context, { citationFailure = false } = {}) {
    const supplemental = context.supplemental_card_data;
    const [submission] = supplemental.clinvar_variant_record.submissions;
    assert.equal(submission.accession, 'SCV004359989');
    assert.equal(submission.contributesToAggregate, true);
    assert.match(submission.classificationComments[0], /inconsistent results in human cell growth suppression assays/);
    assert.match(submission.classificationComments[0], /functional in human cell proliferation assays/);
    assert.match(submission.classificationComments[0], /partially reduced transactivation activity/);
    assert.deepEqual(submission.citations.filter(c => c.scope === 'variant_evidence').map(c => c.pmid).sort(), studyPmids);
    const papers = [...supplemental.pubmed.articles, ...supplemental.clinvar_cited_publications.articles];
    assert.deepEqual(papers.map(a => a.pmid).sort(), studyPmids);
    assert.ok(papers.every(a => a.clinvar_citations[0].accession === 'SCV004359989.2'));
    if (citationFailure) {
        assert.equal(supplemental.clinvar_cited_publications.pubmedLookupUnavailableCount, 2);
        assert.ok(supplemental.clinvar_cited_publications.articles.every(a => a.abstract_status === 'unavailable'));
    } else {
        assert.ok(papers.every(a => a.abstract));
    }
    assert.equal(supplemental.tp53_mutation_database.selection.match_type, 'exact_genomic');
    assert.equal(supplemental.tp53_mutation_database.pathogenicity.somatic_count, '3');
    assert.equal(supplemental.tp53_mutation_database.pathogenicity.splice_ai.DS_AG, '0.08');
    const studies = supplemental.tp53_mutation_database.assay_comparison.studies;
    assert.equal(studies.find(study => study.study.startsWith('Giacomelli')).classifications[0].value, 'Unclassified');
    assert.match(studies.at(-1).variant, /c\.709A>T.*302465/);
}

export async function runEvidenceSmoke(newPage, base, check) {
    const originalFetch = globalThis.fetch;
    const originalDatasetUrl = process.env.TP53_MUTATION_DATASET_URL;
    process.env.TP53_MUTATION_DATASET_URL = datasetUrl;
    const requestedPmids = [];
    const modelPrompts = [];
    let citationFailure = false;
    globalThis.fetch = async (address, options) => {
        const url = new URL(address);
        if (url.hostname === 'eutils.ncbi.nlm.nih.gov' && url.searchParams.get('db') === 'clinvar') {
            assert.equal(url.searchParams.get('id'), '2774649');
            return new Response(markedClinvarXml);
        }
        if (url.hostname === 'eutils.ncbi.nlm.nih.gov' && url.searchParams.get('db') === 'pubmed') {
            requestedPmids.push(url.searchParams.get('id'));
            return citationFailure ? new Response('', { status: 404 }) : new Response(pubmedXml);
        }
        if (url.href === datasetUrl) return new Response(tp53Csv, { headers: { 'Content-Type': 'text/csv' } });
        if (url.hostname === 'openrouter.ai') {
            modelPrompts.push(JSON.parse(options.body).messages[1].content);
            const modelReview = {
                germline_pathogenicity: { classification: 'Uncertain significance', condition: 'Example condition', rationale: 'Germline fixture' },
                somatic_oncogenicity: { classification: 'Insufficient evidence', rationale: 'Somatic fixture' },
                resistance_or_lack_of_benefit: [{
                    therapy: 'Example therapy', tumor_type: 'Example tumor', biomarker_context: 'Example biomarker',
                    evidence: 'Example source', url: 'https://example.invalid/resistance'
                }]
            };
            return Response.json({ choices: [{ message: { content: JSON.stringify(modelReview) } }] });
        }
        throw new Error(`Unexpected server request to ${url.hostname}${url.pathname}`);
    };
    try {
        for (const scenario of ['complete', 'pubmed-failure', 'superseded', 'r175h']) {
            citationFailure = scenario === 'pubmed-failure';
            const annotation = scenario === 'r175h' ? r175hAnnotation : m237lAnnotation;
            const page = await newPage();
            await page.setViewportSize({ width: 390, height: 844 });
            const requests = { clinvar: 0, search: 0, pubmed: 0, tp53: [], ai: [] };
            const gated = scenario === 'complete' || scenario === 'superseded';
            let release = () => {};
            let settledClinvar = () => {};
            const arrived = new Promise(resolve => { settledClinvar = resolve; });
            const withheld = new Promise(resolve => { release = resolve; });
            const json = (route, body) => route.fulfill({ json: body });
            const serve = async (route, handler, url) => {
                const request = route.request();
                const result = await handler({
                    method: request.method(), headers: { origin: 'https://drdoubleb.com' },
                    query: Object.fromEntries(url.searchParams), body: request.postData() ? request.postDataJSON() : undefined
                }, captureResponse());
                await route.fulfill({ status: result.statusCode, json: result.body });
            };
            await page.route('**/*', async route => {
                const request = route.request();
                const url = new URL(request.url());
                const path = decodeURIComponent(url.pathname);
                if (url.hostname === 'grch37.rest.ensembl.org') {
                    if (scenario === 'complete' && path.startsWith('/lookup/')) {
                        return route.fulfill({ status: 500, contentType: 'text/html', body: 'Ensembl unavailable' });
                    }
                    if (path === '/lookup/symbol/homo_sapiens/TP53') return json(route, tp53Map.gene);
                    if (path === `/lookup/id/${tp53Map.transcript.id}` && url.searchParams.get('expand') === '1') {
                        return json(route, tp53Map.transcript);
                    }
                    if (path === `/overlap/translation/${tp53Map.transcript.Translation.id}`) {
                        if (scenario === 'superseded') await withheld;
                        return json(route, tp53Map.domains);
                    }
                }
                if (scenario === 'complete' && url.hostname === 'grch37-archive.rest.ensembl.org') {
                    if (path === '/lookup/symbol/homo_sapiens/TP53' && url.searchParams.get('expand') === '1') {
                        return json(route, tp53Map.archiveGene);
                    }
                    if (path === `/overlap/translation/${tp53Map.transcript.Translation.id}`) return json(route, tp53Map.domains);
                }
                if (url.hostname === 'myvariant.info' && path === `/v1/variant/${annotation._id}`) return json(route, annotation);
                if (url.hostname === 'api.ncbi.nlm.nih.gov' && path.includes('/hgvs/NM_000546.6:c.709A>T/')) {
                    return json(route, { data: { spdis: [transcript], input_hgvs_validity: 'valid' } });
                }
                if (url.hostname === 'api.ncbi.nlm.nih.gov' && path.includes('/spdi/NM_000546.6:850:A:T/')) {
                    return json(route, { data: { spdis: [transcript, grch38, grch37] } });
                }
                if (path === '/api/clinvar-variant' && url.searchParams.get('id') === '2774649') {
                    requests.clinvar++;
                    if (gated) await withheld;
                    await serve(route, clinvar, url);
                    settledClinvar();
                    return;
                }
                if (path === '/api/pubmed' && url.searchParams.has('pmids')) {
                    requests.pubmed++;
                    return serve(route, pubmed, url);
                }
                if (path === '/api/pubmed' && url.searchParams.get('gene') === 'TP53') {
                    requests.search++;
                    return json(route, { backend: 'litvar', total: 1, articles: [pubmedRecords.get('12826609')] });
                }
                if (path === '/api/tp53') {
                    requests.tp53.push(request.postDataJSON());
                    return serve(route, tp53, url);
                }
                if (path === '/api/ai-review') {
                    requests.ai.push(request.postDataJSON());
                    return serve(route, aiReview, url);
                }
                return route.fallback();
            });
            try {
                const variant = encodeURIComponent(scenario === 'r175h' ? annotation._id : 'NM_000546.6:c.709A>T');
                await page.goto(`${base}/?variant=${variant}`, { waitUntil: 'domcontentloaded' });
                const tp53Card = page.locator('[data-card="tp53-database"]');
                if (scenario === 'r175h') {
                    const comparison = tp53Card.getByText(/^Assay comparison · p\.(?:Arg175His|R175H)$/);
                    await comparison.waitFor({ timeout: 30000 });
                    assert.equal(await comparison.evaluate(summary => summary.parentElement.open), false);
                    assert.ok(await tp53Card.getByText('DNE · Yes | DNE/LOF class · DNE_LOF', { exact: true }).isVisible());
                    await comparison.press('Enter');
                    assert.ok(await tp53Card.getByRole('heading', { name: 'Kato et al. (2003)' }).isVisible());
                    assert.ok(await tp53Card.getByRole('heading', { name: 'Kakudo et al. (2005)' }).isVisible());
                    assert.ok(await tp53Card.getByRole('cell', { name: '-2.2 % of wild type', exact: true }).isVisible());
                    assert.ok(await tp53Card.getByRole('cell', { name: '0 % of wild type', exact: true }).isVisible());
                    assert.ok(await tp53Card.getByRole('cell', { name: /DNE_LOF.*DNE and LOF/ }).isVisible());
                    assert.ok(await tp53Card.locator('.tp53-assay-label[data-tone="red"]').filter({ hasText: 'Non-functional' }).isVisible());
                    assert.ok(await tp53Card.getByRole('cell', { name: /1\.791 Z-score.*Meets DNE cutoff/ }).isVisible());
                    assert.ok(await tp53Card.getByText(/c\.524G>A.*construct 200956/).isVisible());
                    assert.ok(await tp53Card.getByRole('cell', { name: /-0\.145 log₂ scale.*LOF in this assay/ }).isVisible());
                    check('evidence: R175H displays its own Giacomelli Z-scores and exact Kotler construct', true);
                    continue;
                }
                const ai = page.locator('[data-card="optional-ai-review"]');
                const copy = ai.getByRole('button', { name: 'Copy prompt', exact: true });
                const run = ai.getByRole('button', { name: 'Run AI review', exact: true });
                const submissions = page.locator('[data-card="clinvar"] summary').filter({ hasText: 'Submission evidence (1)' });
                await copy.waitFor({ timeout: 30000 });
                assert.equal(requests.clinvar, 1);
                if (gated) assert.ok(await page.getByText('Loading ClinVar submission evidence…').isVisible());
                await ai.locator('details.ai-review-key > summary').click();
                await page.locator('#aiReviewKey').fill('sk-browser-smoke-fixture-1234');
                await run.click();
                await copy.click();
                if (gated) assert.equal(requests.ai.length, 0, 'both buttons wait for the pending ClinVar request');

                if (scenario === 'superseded') {
                    const map = page.locator('#geneMap');
                    await map.getByText('Loading Ensembl protein domains…', { exact: true }).waitFor();
                    await page.locator('#variantInput').fill('BRAF V600E');
                    await page.locator('#variantForm button[type="submit"]').click();
                    await page.locator('[data-card="variant"]').getByText(/p\.Val600Glu/).first().waitFor({ timeout: 30000 });
                    release();
                    await arrived;
                    await page.waitForFunction(() => !document.querySelector('[data-card="tp53-database"]'));
                    // The superseded requests have time to resume after the response.
                    await page.waitForTimeout(250);
                    assert.equal(requests.ai.length, 0);
                    assert.equal(requests.pubmed, 0);
                    assert.equal(await submissions.count(), 0);
                    assert.equal(await map.getByText('TP53 gene map', { exact: true }).count(), 0);
                    assert.equal(await map.getByText(/p53, DNA-binding domain/).count(), 0);
                    check('evidence: superseded lookup neither rendered nor sent stale evidence', true);
                    continue;
                }
                release();
                await arrived;
                await submissions.waitFor({ timeout: 10000 });
                const copiedPrompt = ai.locator('.ai-review-prompt-view textarea');
                await page.waitForFunction(() => document.querySelector('.ai-review-prompt-view textarea').value.includes('End of supplied context.'));
                await ai.getByText('Insufficient disease context', { exact: true }).waitFor();
                await ai.getByText('Germline pathogenicity', { exact: true }).waitFor();
                await ai.getByText('Somatic oncogenicity', { exact: true }).waitFor();
                assert.equal(await ai.getByRole('link', { name: 'Example therapy' }).getAttribute('href'), 'https://example.invalid/resistance');
                assert.ok(await ai.getByText(/Example tumor.*Example biomarker.*Example source/).isVisible());
                const prompts = [await copiedPrompt.inputValue(), modelPrompts.at(-1)];
                for (const prompt of prompts) {
                    pageEvidence(promptContext(prompt), { citationFailure });
                    assert.match(prompt, /untrusted evidence/);
                    assert.match(prompt, /assay system, expression level, and controls/);
                }
                assert.deepEqual(requests.ai.map(request => request.mode || 'run').sort(), ['prompt', 'run']);
                assert.equal(requests.clinvar, 1);
                assert.equal(requests.search, 1);
                assert.equal(requests.pubmed, 1);
                assert.equal(requestedPmids.at(-1), '29979965,30224644');
                assert.equal(requests.tp53.length, 1);
                assert.match(requests.tp53[0].cdna, /^NM_000546\.[56]:c\.709A>T$/);
                const allele = requests.tp53[0].genomic_allele;
                assert.deepEqual({ ...allele, pos: Number(allele.pos) }, { assembly: 'GRCh37', chrom: '17', pos: 7577572, ref: 'T', alt: 'A' });
                assert.equal(await submissions.evaluate(summary => summary.parentElement.open), false);
                const citations = page.getByText('ClinVar-cited publications (3)', { exact: true });
                assert.equal(await citations.evaluate(summary => summary.parentElement.open), false);
                if (citationFailure) {
                    await citations.click();
                    assert.ok(await page.getByText(/could not be retrieved for 2 cited PMIDs/).isVisible());
                    check('evidence: failed PubMed metadata preserved the rationale and citation links in both prompts', true);
                    continue;
                }
                await checkReleaseAssets(page);
                await checkM237lRelease(page);
                const map = page.getByRole('region', { name: 'TP53 gene map' });
                const exonPlot = map.getByRole('img', { name: /11 exons.*chr17:7,577,572/ });
                const proteinPlot = map.getByRole('img', { name: /393 amino acids.*residue 237/ });
                await proteinPlot.waitFor();
                assert.ok(await map.getByText(/Canonical transcript · Ensembl GRCh37 archive/).isVisible());
                assert.ok(await map.getByRole('heading', { name: /Exons · GRCh37/ }).isVisible());
                const recordLink = map.getByRole('link', { name: /^ENST00000269305\.4/ });
                assert.equal(new URL(await recordLink.getAttribute('href')).hostname, 'grch37-archive.rest.ensembl.org');
                check('gene map: the GRCh37 archive supplies the transcript and domains after a primary lookup fails', true);
                assert.ok(await exonPlot.locator('.gene-map-pointer').isVisible());
                assert.ok(await proteinPlot.locator('.gene-map-pointer').evaluate(el => {
                    const track = el.parentElement.getBoundingClientRect();
                    const pointer = el.getBoundingClientRect();
                    const fraction = (pointer.left + pointer.width / 2 - track.left) / track.width;
                    return fraction > 0.60 && fraction < 0.61;
                }));
                assert.ok(await map.getByRole('link', { name: /p53, DNA-binding domain.*95–288/ }).isVisible());
                const exonLabels = exonPlot.locator('.gene-map-inline-label:visible');
                const domainLabels = proteinPlot.locator('.gene-map-inline-label:visible');
                assert.deepEqual(await domainLabels.allTextContents(), [
                    'p53 transactivation domain', 'p53, DNA-binding domain', 'p53, tetramerisation domain'
                ]);
                assert.ok(await exonLabels.filter({ hasText: /^7$/ }).isVisible());
                // The same in-track labels must fit after resizing the result.
                await page.setViewportSize({ width: 320, height: 844 });
                await page.waitForFunction(() => {
                    const plots = [...document.querySelectorAll('#geneMap .gene-map-diagram')];
                    return plots.length === 2 && plots.every(plot => {
                        const frame = plot.getBoundingClientRect();
                        const boxes = [...plot.querySelectorAll('.gene-map-inline-label.gene-map-label-placed')]
                            .map(label => label.getBoundingClientRect());
                        return boxes.length > 0 && boxes.every((box, index) => box.left >= frame.left - 1 && box.right <= frame.right + 1
                            && boxes.slice(index + 1).every(other => box.right <= other.left + 1 || box.left >= other.right - 1
                                || box.bottom <= other.top + 1 || box.top >= other.bottom - 1));
                    });
                });
                await page.setViewportSize({ width: 390, height: 844 });
                check('gene map: exon and domain labels fit in the tracks while the linked legend stays visible', true);
                assert.ok(await map.evaluate(el => el.scrollWidth <= el.clientWidth
                    && el.getBoundingClientRect().right <= innerWidth
                    && Boolean(el.compareDocumentPosition(document.querySelector('#cardsContainer')) & Node.DOCUMENT_POSITION_FOLLOWING)));
                check('gene map: M237L is marked on the reverse-strand exon track and in the DNA-binding domain at 390px', true);
                assert.ok(await tp53Card.getByText('DNE · Unavailable | DNE/LOF class · Unclassified', { exact: true }).isVisible());
                check('release smoke: expected assets, exact M237L and live-evidence card contract', true);
                const secondCopy = ai.getByRole('button', { name: 'Copy prompt', exact: true });
                await Promise.all([
                    page.waitForResponse(response => response.url().includes('/api/ai-review') && response.request().postDataJSON().mode === 'prompt'),
                    secondCopy.click()
                ]);
                assert.equal(requests.clinvar, 1);
                assert.equal(requests.pubmed, 1);
                assert.equal(requests.search, 1);

                await submissions.press('Enter');
                const scv = page.locator('[data-card="clinvar"] summary').filter({ hasText: /Color Health/ });
                assert.equal(await scv.evaluate(summary => summary.parentElement.open), false);
                await scv.press('Enter');
                assert.ok(await page.getByRole('link', { name: 'PMID 12826609 ↗' }).isVisible());
                assert.ok(await page.getByRole('link', { name: 'ClinVar submission SCV004359989.2 ↗' }).isVisible());
                const comment = page.locator('[data-card="clinvar"] .evidence-comment').filter({ hasText: /partially reduced transactivation/ });
                assert.ok(await comment.isVisible());
                assert.ok((await comment.textContent()).includes('<strong>literal</strong>'));
                assert.equal(await comment.locator('strong').count(), 0);
                const comparison = tp53Card.getByText('Assay comparison · p.Met237Leu', { exact: true });
                assert.equal(await comparison.evaluate(summary => summary.parentElement.open), false);
                await comparison.press('Enter');
                assert.ok(await tp53Card.getByRole('heading', { name: 'Kato et al. (2003)' }).isVisible());
                assert.ok(await tp53Card.getByRole('cell', { name: '78.2 % of wild type', exact: true }).isVisible());
                assert.ok(await tp53Card.getByRole('cell', { name: 'Unclassified', exact: true }).isVisible());
                assert.ok(await tp53Card.getByText(/c\.709A>T.*construct 302465/).isVisible());
                assert.ok(await tp53Card.getByRole('cell', { name: '-0.156 Z-score', exact: true }).isVisible());
                for (const [measurement, label, tone, background] of [
                    ['Median across eight promoters', 'Partially functional', 'orange', 'rgb(255, 237, 213)'],
                    ['p53 WT · nutlin-3', 'Meets DNE cutoff', 'orange', 'rgb(255, 237, 213)'],
                    ['p53 null · etoposide', 'No LOF in this assay', 'green', 'rgb(220, 252, 231)'],
                    ['H1299 relative fitness score', 'No LOF in this assay', 'green', 'rgb(220, 252, 231)']
                ]) {
                    const row = tp53Card.getByRole('row').filter({ has: page.getByRole('rowheader', { name: measurement, exact: true }) });
                    const badge = row.getByText(label, { exact: true });
                    assert.ok(await badge.isVisible());
                    assert.equal(await badge.getAttribute('data-tone'), tone);
                    assert.equal(await badge.evaluate(el => getComputedStyle(el).backgroundColor), background);
                }
                assert.ok(await tp53Card.getByRole('cell', { name: /-1\.666 log₂ scale.*No LOF in this assay/ }).isVisible());
                assert.ok(await tp53Card.getByText(/NCI WT nutlin-3 Z-score.*≥ 0\.61/).isVisible());
                assert.ok(await tp53Card.getByText(/NCI null etoposide Z-score.*≤ −0\.21/).isVisible());
                assert.ok(await tp53Card.getByText(/ClinGen H1299 RFS.*< −1 .*≥ −1/).isVisible());
                assert.ok(await tp53Card.getByRole('link', { name: 'ClinGen TP53 rules v2.4 ↗' }).isVisible());
                const predictions = page.getByText('Computational predictions', { exact: true }).last();
                await predictions.click();
                assert.ok(await page.getByText('Structure-based computational prediction', { exact: true }).isVisible());
                for (const card of ['clinvar', 'tp53-database', 'pubmed']) {
                    const fits = await page.locator(`[data-card="${card}"]`).evaluate(el => el.scrollWidth <= el.clientWidth && el.getBoundingClientRect().right <= innerWidth);
                    assert.ok(fits, `${card} fits a 390px viewport`);
                }
                await submissions.press('Enter');
                assert.equal(await submissions.evaluate(summary => summary.parentElement.open), false);
                check('evidence: both AI paths used the unopened sections and reused the pending and cached requests', true);
                check('evidence: safe rendering, keyboard toggles and narrow layout at 390px', true);
            } finally {
                release();
                await page.close();
            }
        }
    } finally {
        globalThis.fetch = originalFetch;
        if (originalDatasetUrl === undefined) delete process.env.TP53_MUTATION_DATASET_URL;
        else process.env.TP53_MUTATION_DATASET_URL = originalDatasetUrl;
    }
}
