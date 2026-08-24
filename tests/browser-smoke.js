// Hermetic browser smoke test for the card UI (gene-only and variant modes).
//
// Not part of `npm test` (it needs a browser): run it manually around DOM
// refactors. All external hosts are mocked with Playwright routes: MyVariant
// answers with a canned (real) BRAF V600E payload, Ensembl is DOWN (aborted) —
// which also re-verifies the recoder-outage canonical-transcript fix — and the
// backend proxies return empty-but-valid shapes, so no network access is
// needed and results are deterministic.
//
// Run:
//   npm i --no-save playwright-core        (once; or have playwright installed)
//   python3 -m http.server 8123 &          (serve the repo root)
//   CHROMIUM_PATH=/path/to/chromium node tests/browser-smoke.js
// CHROMIUM_PATH defaults to /opt/pw-browsers/chromium (the Claude Code web
// sandbox's preinstalled build); locally, `npx playwright install chromium`
// and point CHROMIUM_PATH at it, or omit executablePath by editing below.
const { chromium } = require('playwright-core');
const CHROMIUM_PATH = process.env.CHROMIUM_PATH || '/opt/pw-browsers/chromium';

const BASE = 'http://127.0.0.1:8123';
const results = [];
let failures = 0;
function check(name, cond, detail) {
    results.push(`${cond ? 'ok  ' : 'FAIL'} ${name}${!cond && detail ? ` — ${detail}` : ''}`);
    if (!cond) failures++;
}

// Real dbNSFP/snpEff shape fetched live from MyVariant for chr7:g.140453136A>T.
const BRAF_ANNOTATION = {
    _id: 'chr7:g.140453136A>T',
    chrom: '7',
    hg19: { start: 140453136, end: 140453136 },
    vcf: { ref: 'A', alt: 'T', position: '140453136' },
    cadd: { consequence: 'NON_SYNONYMOUS', chrom: '7', phred: 32 },
    clinvar: { variant_id: 13961, gene: { symbol: 'BRAF' } },
    dbnsfp: {
        genename: ['BRAF', 'BRAF', 'BRAF', 'BRAF'],
        hgvsc: ['c.620T>A', 'c.1919T>A', 'c.1799T>A'],
        hgvsp: ['p.Val640Glu', 'p.Val600Glu', 'p.Val207Glu', 'p.V600E'],
        ensembl: { transcriptid: ['ENST00000496384', 'ENST00000644969', 'ENST00000646891', 'ENST00000288602'] }
    },
    snpeff: { ann: { feature_id: 'NM_004333.4', hgvs_c: 'c.1799T>A', hgvs_p: 'p.Val600Glu', genename: 'BRAF' } }
};

async function installRoutes(page) {
    await page.route('**/*', (route) => {
        const url = route.request().url();
        if (url.startsWith(BASE) || url.includes('127.0.0.1')) return route.continue();
        const json = (body) => route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify(body) });
        // The assembly-map endpoint answers (the GRCh38 toggle depends on it);
        // the rest of Ensembl stays down. V600E GRCh38→GRCh37.
        if (/ensembl\.org\/map\/human\/GRCh38\/7:140753336\.\.140753336/.test(url)) return json({
            mappings: [{ mapped: { assembly: 'GRCh37', seq_region_name: '7', start: 140453136, end: 140453136, strand: 1 } }]
        });
        if (/ensembl\.org/.test(url)) return route.abort(); // Ensembl outage
        if (/myvariant\.info\/v1\/query/.test(url)) return json({ hits: [BRAF_ANNOTATION] });
        if (/myvariant\.info\/v1\/variant/.test(url)) {
            // Only the true hg19 id resolves — an unlifted GRCh38 coordinate must
            // NOT find an annotation, or the toggle test would pass vacuously.
            return /140453136/.test(url)
                ? json(BRAF_ANNOTATION)
                : route.fulfill({ status: 404, contentType: 'application/json', body: JSON.stringify({ success: false, error: 'not found' }) });
        }
        if (/cancer-types\.php/.test(url)) return json({ count: 1, cancer_types: [{ name: 'Melanoma', record_count: 3, aliases: ['melanoma'] }] });
        if (/guidelines\/api\/search\.php/.test(url)) return json({ count: 0, results: [], query: {} });
        if (/cbioportal\.org\/api\/genes\//.test(url)) return json({ entrezGeneId: 673, hugoGeneSymbol: 'BRAF' });
        if (/cbioportal\.org\/api\/molecular-profiles\/.*\/mutations\/fetch/.test(url)) return json([
            { sampleId: 'S1', proteinChange: 'V600E' },
            { sampleId: 'S2', proteinChange: 'V600E' },
            { sampleId: 'S3', proteinChange: 'G469A' }
        ]);
        if (/cbioportal\.org\/api\/clinical-data\/fetch/.test(url)) return json([
            { sampleId: 'S1', clinicalAttributeId: 'CANCER_TYPE', value: 'Melanoma' },
            { sampleId: 'S2', clinicalAttributeId: 'CANCER_TYPE', value: 'Thyroid Cancer' },
            { sampleId: 'S3', clinicalAttributeId: 'CANCER_TYPE', value: 'Melanoma' }
        ]);
        if (/api\/pubmed/.test(url)) {
            // Variant queries get the LitVar2-shaped answer; term queries stay empty.
            return url.includes('variant=')
                ? json({
                    total: 12,
                    articles: [{ pmid: '1', title: 'Mock LitVar paper', authors: 'A, B et al.', journal: 'J', year: '2024', abstract: '' }],
                    litvar: { id: '@VARIANT_p.V600E_BRAF_human', name: 'p.V600E', rsid: 'rs113488022', url: 'https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?variant=litvar%40rs113488022%23%23' },
                    backend: 'litvar'
                })
                : json({ total: 0, articles: [], backend: 'term' });
        }
        if (/api\/clinvar/.test(url)) return json({ variants: [], total: 0 });
        if (/api\/clinicaltrials/.test(url)) return json({ total: 0, studies: [] });
        if (/api\.fda\.gov/.test(url)) return json({ meta: { results: { total: 0 } }, results: [] });
        if (/erepo\.clinicalgenome\.org/.test(url)) return json({
            variantInterpretations: [{
                gene: { label: 'BRAF' },
                '@id': 'https://erepo.genome.network/evrepo/api/interpretation/CA123/MONDO:1/001',
                caid: 'CAR:CA123',
                publishedDate: '2024-01-01',
                guidelines: [{ outcome: { label: 'Pathogenic' }, agents: [{ affiliation: 'Mock VCEP', outcome: { label: 'Pathogenic' } }] }]
            }]
        });
        if (/gnomad\.broadinstitute\.org\/api/.test(url)) return json({
            data: { gene: { gnomad_constraint: { pli: 0.9999, oe_lof: 0.1, oe_lof_upper: 0.209, mis_z: 3.72, oe_mis: 0.487, syn_z: -0.22 } } }
        });
        if (/dgidb\.org\/api\/graphql/.test(url)) return json({
            data: { genes: { nodes: [{ name: 'BRAF', interactions: [
                { drug: { name: 'DABRAFENIB', approved: true }, interactionScore: 0.5, interactionTypes: [{ type: 'inhibitor' }], sources: [{ sourceDbName: 'CKB' }] },
                { drug: { name: 'MOCKDRUG', approved: false }, interactionScore: 0.1, interactionTypes: [], sources: [] }
            ] }] } }
        });
        return json({});
    });
}

(async () => {
    const browser = await chromium.launch({
        executablePath: CHROMIUM_PATH,
        headless: true,
        args: ['--no-sandbox']
    });
    const ctx = await browser.newContext();

    const pageErrors = [];
    const newPage = async () => {
        const page = await ctx.newPage();
        page.on('pageerror', (err) => pageErrors.push(String(err)));
        await installRoutes(page);
        return page;
    };

    // ── Gene-only mode ──────────────────────────────────────────────────────
    {
        const page = await newPage();
        await page.goto(`${BASE}/?variant=BRAF`, { waitUntil: 'domcontentloaded' });
        await page.waitForSelector('#cardsContainer .card', { timeout: 30000 });
        await page.waitForTimeout(3000);
        const titles = await page.$$eval('#cardsContainer .card > h3', els => els.map(e => e.textContent));
        check('gene-only: card set', JSON.stringify(titles) === JSON.stringify(
            ['CIViC', 'OncoKB', 'PubMed', 'FDA-Approved Drugs (by gene)', 'Clinical Trials', 'Cancer Prevalence', 'Guidelines', 'Optional AI Review']),
            JSON.stringify(titles));
        const cbioGeneText = await page.$eval('[data-card="cancer-prevalence"]', el => el.textContent);
        check('gene-only: prevalence gene stats', /BRAF mutated:\s*3 of 10,945/.test(cbioGeneText.replace(/\s+/g, ' ')), cbioGeneText.slice(0, 200));
        const pmText = await page.$eval('[data-card="pubmed"]', el => el.textContent);
        check('gene-only: PubMed panel', /Search PubMed/.test(pmText) && /Query: "BRAF"/.test(pmText), pmText.slice(0, 120));
        const fdaTabs = await page.$$eval('[data-card="fda-approved-drugs-by-gene"] .card-tab-btn', els => els.map(e => e.textContent));
        check('gene-only: FDA tabs', JSON.stringify(fdaTabs) === JSON.stringify(['Companion Dx', 'openFDA Labels', 'BBKB', 'DGIdb']), JSON.stringify(fdaTabs));
        await page.click('[data-card="fda-approved-drugs-by-gene"] .card-tab-btn:nth-child(2)');
        const activeTab = await page.$eval('[data-card="fda-approved-drugs-by-gene"] .card-tab-btn.active', el => el.textContent);
        check('gene-only: FDA tab switch', activeTab === 'openFDA Labels', activeTab);
        const dgidbText = (await page.$eval('[data-card="fda-approved-drugs-by-gene"]', el => el.textContent)).replace(/\s+/g, ' ');
        check('gene-only: DGIdb panel populated', /DABRAFENIB.*approved · inhibitor · score 0\.50/.test(dgidbText), dgidbText.slice(0, 300));
        const ctText = await page.$eval('[data-card="clinical-trials"]', el => el.textContent);
        check('gene-only: Trials panel', /Search ClinicalTrials.gov/.test(ctText) && /Query: "BRAF"/.test(ctText), ctText.slice(0, 120));
        const glText = await page.$eval('[data-card="guidelines"]', el => el.textContent);
        check('gene-only: Guidelines rendered', /Select a cancer type/.test(glText), glText.slice(0, 120));
        await page.close();
    }

    // ── Variant mode with tumor type; Ensembl down (recoder-outage path) ────
    {
        const page = await newPage();
        await page.goto(`${BASE}/?variant=braf%20v600e&tumorType=melanoma`, { waitUntil: 'domcontentloaded' });
        await page.waitForSelector('[data-card="variant"]', { timeout: 60000 });
        await page.waitForSelector('[data-card="guidelines"]', { timeout: 30000 });
        await page.waitForTimeout(3000);
        const titles = await page.$$eval('#cardsContainer .card > h3', els => els.map(e => e.textContent));
        check('variant: extracted cards present',
            ['PubMed', 'FDA-Approved Drugs (by gene)', 'Clinical Trials', 'Cancer Prevalence', 'Guidelines', 'Optional AI Review']
                .every(t => titles.includes(t)), JSON.stringify(titles));
        const cbioText = (await page.$eval('[data-card="cancer-prevalence"]', el => el.textContent)).replace(/\s+/g, ' ');
        check('variant: prevalence variant stats', /V600E specifically:\s*2 tumors/.test(cbioText), cbioText.slice(0, 250));
        check('variant: prevalence tumor-type breakdown (only requested samples counted)',
            /Tumor types with V600E:.*Melanoma \(1\) · Thyroid Cancer \(1\)/.test(cbioText), cbioText.slice(0, 320));
        const pmLitvarText = (await page.$eval('[data-card="pubmed"]', el => el.textContent)).replace(/\s+/g, ' ');
        check('variant: LitVar2 provenance line', /Variant-matched via LitVar2: p\.V600E \(rs113488022\)/.test(pmLitvarText), pmLitvarText.slice(0, 250));
        const clinvarText = (await page.$eval('[data-card="clinvar"]', el => el.textContent)).replace(/\s+/g, ' ');
        check('variant: ClinGen expert-panel line', /ClinGen expert panel: Pathogenic — Mock VCEP \(2024-01-01\)/.test(clinvarText), clinvarText.slice(0, 300));
        const gnomadText = (await page.$eval('[data-card="gnomad"]', el => el.textContent)).replace(/\s+/g, ' ');
        check('variant: gnomAD gene constraint', /Gene constraint · BRAF.*pLI 1\.00 · LOEUF 0\.21 · missense Z 3\.72/.test(gnomadText), gnomadText.slice(0, 300));
        const variantText = await page.$eval('[data-card="variant"]', el => el.textContent);
        check('variant: canonical cDNA is c.1799T>A with recoder down',
            /c\.:\s*c\.1799T>A/.test(variantText.replace(/\s+/g, ' ')), variantText.slice(0, 250));
        check('variant: canonical protein p.Val600Glu', /p\.Val600Glu/.test(variantText), variantText.slice(0, 250));
        const pmTabs = await page.$$eval('[data-card="pubmed"] .card-tab-btn', els => els.map(e => e.textContent));
        check('variant: PubMed 3 tabs with tumor type', JSON.stringify(pmTabs) === JSON.stringify(
            ['Gene + Variant', 'Gene + Tumor Type', 'Gene + Variant + Tumor Type']), JSON.stringify(pmTabs));
        const pmText = await page.$eval('[data-card="pubmed"]', el => el.textContent);
        check('variant: PubMed default query', /Query: "BRAF V600E"/.test(pmText), pmText.slice(0, 160));
        const fdaText = await page.$eval('[data-card="fda-approved-drugs-by-gene"]', el => el.textContent);
        check('variant: FDA gene label', /Gene: BRAF/.test(fdaText), fdaText.slice(0, 120));
        const ctText = await page.$eval('[data-card="clinical-trials"]', el => el.textContent);
        check('variant: Trials query includes tumor type', /Query: "BRAF \+ melanoma"/.test(ctText), ctText.slice(0, 160));
        const glText = await page.$eval('[data-card="guidelines"]', el => el.textContent);
        check('variant: Guidelines auto-selected tumor type ran a search',
            /No guideline records found for BRAF in Melanoma/.test(glText), glText.slice(0, 200));
        await page.close();
    }

    // ── GRCh38 input toggle: hg38 coordinate lifted before the pipeline ────
    {
        const page = await newPage();
        await page.goto(`${BASE}/?variant=chr7%3Ag.140753336A%3ET&assembly=grch38`, { waitUntil: 'domcontentloaded' });
        check('grch38: checkbox pre-checked from ?assembly=grch38',
            await page.$eval('#grch38Toggle', el => el.checked));
        await page.waitForSelector('[data-card="variant"]', { timeout: 60000 });
        await page.waitForTimeout(2000);
        const lpText = (await page.$eval('#lookupProgress', el => el.textContent)).replace(/\s+/g, ' ');
        check('grch38: liftover step shows the mapping',
            /140753336.*GRCh38.*→.*140453136/.test(lpText), lpText.slice(0, 300));
        const variantText = (await page.$eval('[data-card="variant"]', el => el.textContent)).replace(/\s+/g, ' ');
        check('grch38: pipeline ran on the GRCh37 locus (canonical c.1799T>A)',
            /c\.1799T>A/.test(variantText), variantText.slice(0, 250));
        const gUrl = page.url();
        check('grch38: assembly param kept in the shareable URL', /assembly=grch38/.test(gUrl), gUrl);
        await page.close();
    }

    await browser.close();
    check('no uncaught page errors', pageErrors.length === 0, pageErrors.join(' | ').slice(0, 400));
    console.log(results.join('\n'));
    console.log(failures ? `\n${failures} FAILED` : '\nAll smoke checks passed');
    process.exit(failures ? 1 : 0);
})().catch((err) => { console.error('SMOKE ERROR:', err); process.exit(2); });
