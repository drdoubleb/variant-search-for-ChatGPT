/*
 * A gene named beside coordinates.
 *
 * A pasted row gives the allele by its coordinates and, usually, a gene column:
 * "12 25398285 KRAS2 C CACC". So does "ATP6 m.8993T>G". The coordinates place
 * the allele; the gene is a label that came with them. It was used to choose
 * among the transcripts Ensembl's VEP gives for the allele — which is right
 * where two genes overlap (EGFR and EGFR-AS1) — and where it matched none of
 * them, all of them were dropped:
 *
 *   12 25398285 KRAS2 C CACC   Gene: KRAS2 · c.36_38dup · p.: N/A · Effect: N/A
 *   12 25398284 TP53 C T       Gene: KRAS  · c.35G>A    · p.: N/A · Effect: N/A
 *   ATP6 m.8993T>G             Gene: MT-ATP6 · c.: N/A · p.: N/A
 *
 * KRAS2 is an older name for KRAS, ATP6 is how MT-ATP6 is usually written, and
 * TP53 is simply not where chr12:25398284 is. In all three the allele is known
 * exactly, and so is what it does. Where the gene given is not one the allele
 * lies in, the allele's own gene takes its place — for the transcripts, for the
 * card and for every gene-level source — and the card says which of the two it
 * was: another name for the gene (NCBI's gene records say so), or a different
 * gene, which is worth a look at the row.
 *
 * The scenarios replay answers recorded from the live services.
 *
 * Run with: node tests/gene-beside-coordinates.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const note = page => [...page.document.querySelectorAll('.gene-given-note')].map(el => el.textContent);
const wait = ms => new Promise(resolve => setTimeout(resolve, ms));

// ── The pieces ─────────────────────────────────────────────────────────────
{
    const given = h.geneGivenBesideCoordinates;
    // Another name for the gene: NCBI lists KRAS for "KRAS2".
    assert.deepEqual(given('KRAS2', ['KRAS', 'KRAS'], 'KRAS', ['KRAS']), { given: 'KRAS2', gene: 'KRAS', alias: true });
    assert.deepEqual(given(' kras2 ', ['KRAS'], 'KRAS', ['kras']), { given: 'kras2', gene: 'KRAS', alias: true }, 'whatever the case');
    // (A name two genes have had: MLL2 was KMT2D's, and is listed for KMT2B too.)
    assert.deepEqual(given('MLL2', ['KMT2D'], 'KMT2D', ['KMT2D', 'KMT2B']), { given: 'MLL2', gene: 'KMT2D', alias: true });
    // A mitochondrial gene written without its prefix.
    assert.deepEqual(given('ATP6', ['MT-ATP6'], 'MT-ATP6', []), { given: 'ATP6', gene: 'MT-ATP6', alias: true });
    assert.deepEqual(given('ATP6', ['MT-ATP6'], 'MT-ATP6', null), { given: 'ATP6', gene: 'MT-ATP6', alias: true });
    // A different gene — or a name nobody has said is this gene's.
    assert.deepEqual(given('TP53', ['KRAS'], 'KRAS', []), { given: 'TP53', gene: 'KRAS', alias: false });
    assert.deepEqual(given('KRAS2', ['KRAS'], 'KRAS', null), { given: 'KRAS2', gene: 'KRAS', alias: false }, 'NCBI could not say');
    assert.deepEqual(given('KRAS2', ['KRAS'], 'KRAS', ['NRAS']), { given: 'KRAS2', gene: 'KRAS', alias: false }, 'another gene\'s name');
    // The gene given is one the allele lies in: it is what was asked for, and
    // chooses among overlapping genes as it always did.
    assert.equal(given('EGFR-AS1', ['EGFR', 'EGFR', 'EGFR-AS1'], 'EGFR', []), null);
    assert.equal(given('kras', ['KRAS'], 'KRAS', []), null);
    // Nothing says otherwise: VEP placed the allele in no transcript (a
    // promoter variant keeps the gene it came with), or named no gene.
    assert.equal(given('TERT', [], 'TERT', []), null);
    assert.equal(given('TERT', [], '', []), null);
    assert.equal(given('TERT', [], 'CLPTM1L', ['CLPTM1L']), null, 'no transcript row names a gene: none of them says the name is wrong');
    assert.equal(given('TERT', ['', null], 'CLPTM1L', []), null);
    assert.equal(given('KRAS2', ['KRAS'], '', ['KRAS']), null);
    assert.equal(given('KRAS2', ['CHR12'], 'CHR12', ['KRAS']), null, 'a chromosome is not a gene');
    // The gene that takes the name's place is one whose transcripts the
    // allele is in — there are no rows to show for any other.
    assert.equal(given('TP53', ['KRAS'], 'NRAS', []), null);
    assert.equal(given('KRAS2', ['KRAS'], 'NRAS', ['NRAS']), null);
    assert.equal(given('', ['KRAS'], 'KRAS', ['KRAS']), null);
    assert.equal(given(null, null, null), null);

    // Whose name a name is: where it is no gene's symbol, the genes NCBI
    // answers with are the ones it is another name for.
    const report = (...products) => ({ reports: products.map(([symbol, ...transcripts]) => ({ product: { symbol,
        transcripts: transcripts.map(([accession_version, select_category, ensembl_transcript]) => ({ accession_version, select_category, ensembl_transcript })) } })) });
    const kras = ['KRAS', ['NM_004985.5', 'MANE_SELECT', 'ENST00000311936.8'], ['NM_033360.4', 'MANE_PLUS_CLINICAL', 'ENST00000256078.10']];
    assert.deepEqual(h.clinicalTranscriptsFromReport(report(kras), 'KRAS2'),
        { select: '', selectEnsembl: '', plusClinical: [], refseqSelect: '', ensembl: [], knownAs: ['KRAS'] });
    assert.deepEqual(h.clinicalTranscriptsFromReport(report(['KMT2D'], ['KMT2B'], ['KMT2D']), 'mll2').knownAs, ['KMT2D', 'KMT2B']);
    // A name that IS a gene's symbol is that gene, whatever else comes with it
    // ("TTN" brings TTR first, which has it as another name).
    assert.deepEqual(h.clinicalTranscriptsFromReport(report(['TTR'], ['TTN', ['NM_001267550.2', 'MANE_SELECT', 'ENST00000589042.5']]), 'TTN').knownAs, []);
    assert.deepEqual(h.clinicalTranscriptsFromReport(report(kras), 'KRAS').knownAs, []);
    assert.deepEqual(h.clinicalTranscriptsFromReport({}, 'NOTAGENE').knownAs, []);
    assert.deepEqual(h.clinicalTranscriptsFromReport(report(['<img onerror=x>'], ['KRAS']), 'KRAS2').knownAs, ['KRAS'], 'only what reads as a symbol');
    console.log('gene-beside-coordinates unit checks passed');
}

// 1. A row whose gene column is an older name for the gene. Was: "Gene: KRAS2 ·
//    c.36_38dup · p.: N/A · Effect: N/A", with CIViC, ClinVar and the rest asked
//    about a gene called KRAS2.
{
    const page = await runLookup(cassettes.geneAliasRow, {});
    assertResolved(assert, page, '12 25398285 KRAS2 C CACC');
    assert.match(page.variant, /g\.: chr12:g\.25398284_25398286dupVCF \(hg19\): 12-25398280-G-GCCA \(left-aligned\)Gene: KRASGiven as KRAS2, which is another name for KRAS\./);
    assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.36_38dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/);
    assert.deepEqual(note(page), ['Given as KRAS2, which is another name for KRAS.']);
    assert.match(page.panel, /You entered\s*12 25398285 KRAS2 C CACC\s*Gene\s*KRAS\s*Genomic/, 'the panel names the gene the allele lies in');
    assert.match(page.panel, /Parse input\s*read as chr12:g\.25398285_25398286insACC · gene given as KRAS2, another name for KRAS/);
    // Every gene-level source is asked about KRAS. The name given is asked
    // about twice only: whose it is (NCBI), and — before anything was known —
    // for its MANE transcript.
    assert.deepEqual(page.api('civic'), ['civic?gene=KRAS&protein=Gly13dup']);
    assert.deepEqual(page.requests.filter(r => /KRAS2/.test(r)).sort(), [
        'api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/KRAS2/taxon/9606/product_report',
        'rest.ensembl.org/lookup/symbol/homo_sapiens/KRAS2?content-type=application/json']);
    assert(page.requests.filter(r => /[?&]gene=/.test(r)).every(r => /[?&]gene=KRAS(?:&|$)/.test(r)), 'no source is asked about another gene');
    console.log('ok gene beside coordinates: an older name for the gene (KRAS2) is shown as the gene, with its protein change');
}

// 2. A row whose gene column is a different gene altogether. Was: "c.35G>A ·
//    p.: N/A · Effect: N/A" for KRAS G12D. The allele is shown as what it is,
//    and the card says the row names another gene.
{
    const page = await runLookup(cassettes.geneWrongRow, {});
    assertResolved(assert, page, '12 25398284 TP53 C T');
    assert.match(page.variant, /g\.: chr12:g\.25398284C>TVCF \(hg19\): 12-25398284-C-TGene: KRAS/);
    assert.match(page.variant, /Transcript: NM_004985\.5 \(MANE Select\)c\.: c\.35G>Ap\.: p\.Gly12Asp \(p\.G12D\)Effect: missense_variant/);
    assert.deepEqual(note(page), ['The gene given with these coordinates is TP53, but the allele lies in KRAS (Ensembl), which is what is shown. Check the row.']);
    assert.match(page.panel, /gene given as TP53 — the allele lies in KRAS, which is what is shown/);
    assert.deepEqual(page.api('civic'), ['civic?gene=KRAS&protein=Gly12Asp']);
    assert(page.requests.filter(r => /[?&]gene=/.test(r)).every(r => /[?&]gene=KRAS(?:&|$)/.test(r)), 'nothing gene-level is asked about TP53');
    // (The record's own annotations are left as they are: the only gene label
    // replaced is one made up for an allele that has no record.)
    assert.match(page.cards.Predictors, /AlphaMissense: P/);
    console.log('ok gene beside coordinates: a different gene in the gene column does not cost the allele its protein change');
}

// 3. A mitochondrial gene written without "MT-". Was: no transcript, "c.: N/A ·
//    p.: N/A · Effect: NON_SYNONYMOUS" for m.8993T>G, MT-ATP6 p.Leu156Arg.
{
    const page = await runLookup(cassettes.geneMitoShort, {});
    assertResolved(assert, page, 'ATP6 m.8993T>G');
    assert.match(page.variant, /g\.: chrMT:g\.8993T>GVCF \(hg19\): MT-8993-T-GGene: MT-ATP6Given as ATP6, which is another name for MT-ATP6\./);
    assert.match(page.variant, /Transcript: ENST00000361899\.2c\.: c\.467T>Gp\.: p\.Leu156Arg \(p\.L156R\)Effect: missense_variant/);
    assert.deepEqual(note(page), ['Given as ATP6, which is another name for MT-ATP6.']);
    console.log('ok gene beside coordinates: ATP6 is MT-ATP6');
}

// 4. Whether a name is another name for the gene is NCBI's to say. When it
//    does not answer, the allele is shown all the same — as a gene that differs
//    from the one given, which is all that is known.
{
    const page = await runLookup(cassettes.geneAliasRow, { without: [/datasets\/v2\/gene\/symbol/] });
    assertResolved(assert, page, 'KRAS2 row, NCBI not answering');
    assert.match(page.variant, /Gene: KRASThe gene given with these coordinates is KRAS2, but the allele lies in KRAS \(Ensembl\), which is what is shown\. Check the row\./);
    assert.match(page.variant, /c\.: c\.36_38dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/);
    // ...and nor is the result held for it: an answer that never completes.
    const never = new Promise(() => {});
    const stalled = await Promise.race([
        runLookup(cassettes.geneAliasRow, { route: (url, key) => (/datasets\/v2\/gene\/symbol/.test(key) ? never : undefined) }),
        wait(8000).then(() => null)]);
    assert(stalled, 'the lookup finishes');
    assertResolved(assert, stalled, 'KRAS2 row, NCBI\'s answer never completing');
    assert.match(stalled.variant, /Gene: KRASThe gene given with these coordinates is KRAS2.*p\.: p\.Gly13dup/);
    console.log('ok gene beside coordinates: without NCBI the allele is still shown, as a gene that differs from the one given');
}

// 5. Nothing changes where there is nothing to say: the gene given is the
//    allele's, or none was given — and where VEP does not answer, nothing is
//    concluded about the gene that was.
{
    for (const input of ['12 25398285 KRAS C CACC', '12 25398285 C CACC']) {
        const page = await runLookup(cassettes.geneAliasRow, { input });
        assertResolved(assert, page, input);
        assert.match(page.variant, /Gene: KRAS(?:Transcript|c\.)/, input);
        assert.match(page.variant, /c\.: c\.36_38dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/, input);
        assert.deepEqual(note(page), [], input);
        assert.doesNotMatch(page.panel, /gene given as/, input);
    }
    const silent = await runLookup(cassettes.geneAliasRow, { without: [/\/vep\//] });
    assertResolved(assert, silent, 'KRAS2 row, VEP not answering');
    assert.match(silent.variant, /Gene: KRAS2Transcript: NM_004985\.5c\.: c\.36_38dupp\.: N\/A/);
    assert.deepEqual(note(silent), []);
    console.log('ok gene beside coordinates: nothing is said where the gene given is the allele\'s, or VEP has not placed the allele');
}

// 6. An overtaken lookup leaves the next one alone. Whose name "KRAS2" is,
//    NCBI is still being asked when "BRAF V600E" is looked up; its answer
//    arrives when the second lookup has its transcript rows and is waiting for
//    its MANE label. The first lookup, no longer the one on the page, would
//    go on to put "its" rows where the page keeps them: BRAF V600E came out as
//    "Gene: BRAF · NM_004985.5 · c.36_38dup · p.Gly13dup", and CIViC was asked
//    about BRAF Gly13dup.
{
    const first = cassettes.geneAliasRow;
    const second = cassettes.brafV600E;
    const both = { input: first.input, reference: [...first.reference, ...second.reference],
        api: { ...(first.api || {}), ...(second.api || {}) }, responses: { ...first.responses, ...second.responses } };
    let release;
    const held = new Promise((resolve) => { release = resolve; });
    let lastRequest = Date.now();
    // (No request for a while: the page is waiting for an answer that is held.)
    const quiet = async (ms) => { while (Date.now() - lastRequest < ms) await wait(10); };
    let firstAsked = false;
    let submitted = false;
    let secondHasRows = false;
    let answeredWhileSecondWaited = false;
    const page = await runLookup(both, {
        // (The page's waits for NCBI run at their real length: both lookups
        // are held inside one.)
        settleMs: 1500, timerCap: 4000,
        route: async (url, key, respond, options, { document }) => {
            lastRequest = Date.now();
            // NCBI's answer about KRAS2 is held until the second lookup has its rows.
            if (/datasets\/v2\/gene\/symbol\/KRAS2\//.test(key)) {
                await held;
                return undefined;
            }
            // The first lookup asks VEP about the allele; once it has the rows
            // it asks whose name KRAS2 is, and waits. Then the next is submitted.
            if (!firstAsked && /\/vep\/human\/hgvs\/12:g\.25398284_25398286dup/.test(key)) {
                firstAsked = true;
                quiet(150).then(() => {
                    submitted = true;
                    document.getElementById('variantInput').value = 'BRAF V600E';
                    document.getElementById('variantForm').dispatchEvent(new document.defaultView.Event('submit', { cancelable: true }));
                });
            }
            if (submitted && /\/vep\/human\/hgvs\/7:g\.140453136A>T/.test(key)) secondHasRows = true;
            // The second lookup has its rows and waits for its MANE label —
            // which comes only after NCBI has answered the first one.
            if (submitted && /lookup\/symbol\/homo_sapiens\/BRAF|datasets\/v2\/gene\/symbol\/BRAF\//.test(key)) {
                for (let i = 0; i < 600 && !secondHasRows; i++) await wait(5);
                await quiet(150);
                answeredWhileSecondWaited = secondHasRows;
                release();
                await wait(150);
                return undefined;
            }
            return undefined;
        }
    });
    assert(firstAsked && submitted && answeredWhileSecondWaited,
        'NCBI answered the first lookup while the second waited with its rows');
    assertResolved(assert, page, 'BRAF V600E, submitted while the KRAS2 row waited for NCBI');
    assert.match(page.panel, /You entered\s*BRAF V600E\s*Gene\s*BRAF.*Coding\s*NM_004333\.6:c\.1799T>A\s*Protein\s*p\.Val600Glu/);
    assert.match(page.variant, /Gene: BRAFTranscript: NM_004333\.6 \(MANE Select\)c\.: c\.1799T>Ap\.: p\.Val600Glu \(p\.V600E\)Effect: missense_variant/);
    assert.deepEqual(page.api('civic'), ['civic?gene=BRAF&protein=Val600Glu']);
    assert.deepEqual(note(page), []);
    assert.doesNotMatch(page.panel + page.variant, /KRAS|NM_004985|Gly13dup|gene given as/, 'nothing of the lookup that was overtaken');
    console.log('ok gene beside coordinates: an overtaken lookup leaves the next one\'s gene and rows alone');
}

console.log('Gene beside coordinates tests passed');
