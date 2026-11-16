/*
 * Allele identity for indels: the same variant spelled two ways must resolve,
 * and two different variants spelled alike must not be confused.
 *
 * Three production failures are pinned here, each reproduced against the live
 * services before the fix (fixtures are those recorded responses):
 *
 *   1. "X 20156693 EIF1AX C CAGATTCATTCTCATTCTTACCCCTGCGTCT" failed with
 *      "Variant identity mismatch: annotation belongs to a different genomic
 *      allele". Ensembl resolves it — to chrX:g.20156694_20156723dup, the HGVS
 *      spelling of the same 30-bp insertion — and the identity check compared
 *      the two strings. Every indel that Ensembl re-spells and MyVariant does
 *      not index went the same way, including the README's own example
 *      "X 20148674 EIF1AX T TG".
 *
 *   2. "chr13:g.32906698_32906699insCA" (BRCA2) returned a DIFFERENT variant,
 *      silently: g.32906699_32906700insCA, c.1084_1085insCA, ClinVar 2505518
 *      "Likely pathogenic". The SPDI → HGVS conversion placed every insertion
 *      one base too far 3', and that neighbour happens to be in MyVariant.
 *
 *   3. "M 8993 T G" failed with "Variant identity mismatch: the resolved allele
 *      does not match the genomic input". It resolves to chrMT:g.8993T>G, and
 *      chrM / chrMT were compared as different chromosomes.
 *
 * Run with: node tests/indel-identity.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { parseHTML } from 'linkedom';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const fixture = JSON.parse(readFileSync(new URL('./fixtures/indel-identity.json', import.meta.url)));
const { eif1axDup, eif1axReadme, brca2Ins, cftr, mito } = fixture;

// Recorded GRCh37 blocks; hand any sub-range back the way Ensembl would.
const referenceSlice = (ref, start, end) => {
    assert(start >= ref.start && end < ref.start + ref.seq.length, `fixture does not cover ${start}..${end}`);
    return ref.seq.slice(start - ref.start, end - ref.start + 1);
};
const same = (a, b, ref) => h.sameGenomicAllele(a, b, ref.start, ref.seq);

// ── SPDI → HGVS ────────────────────────────────────────────────────────────
// An SPDI position is a 0-based interbase coordinate, so an insertion at p
// sits between 1-based bases p and p+1. Ensembl echoes a genomic insertion
// back at exactly that position — converting its SPDI must reproduce the input.
for (const fx of [eif1axDup, eif1axReadme, brca2Ins]) {
    const allele = Object.values(fx.recoder[0])[0];
    const spdi = [].concat(allele.spdi).find(s => s.startsWith('NC_'));
    assert.equal(h.convertSpdiToMyVariant(spdi), fx.query, `${spdi} must round-trip to the input`);
}
assert.equal(h.convertSpdiToMyVariant('NC_000013.10:32906698::CA'), 'chr13:g.32906698_32906699insCA');
assert.notEqual(h.convertSpdiToMyVariant('NC_000013.10:32906698::CA'), brca2Ins.neighbourId,
    'the one-base-3\' reading is ClinVar 2505518, a different allele');
// The other SPDI shapes were always right and must stay that way.
assert.equal(h.convertSpdiToMyVariant('NC_000007.13:140453135:A:T'), 'chr7:g.140453136A>T');
assert.equal(h.convertSpdiToMyVariant('NC_000007.13:117199644:TCT:'), 'chr7:g.117199645_117199647del');
assert.equal(h.convertSpdiToMyVariant('NC_000013.10:32914437:T:'), 'chr13:g.32914438T>-');
assert.equal(h.convertSpdiToMyVariant('NC_000007.13:140453135:AC:TT'), 'chr7:g.140453136AC>TT');
assert.throws(() => h.convertSpdiToMyVariant('NC_000007.13:140453135::'), /no sequence change/);
// NCBI's contextual "G→GG" trims to two placements of the same duplicated G;
// both must now be real spellings of chrX:g.20148675dup, not one base beyond.
{
    const forms = h.minimalSpdiForms({ seq_id: 'NC_000023.10', position: 20148674, deleted_sequence: 'G', inserted_sequence: 'GG' })
        .map(h.convertSpdiToMyVariant);
    assert.deepEqual(forms, ['chrX:g.20148675_20148676insG', 'chrX:g.20148674_20148675insG']);
    for (const form of forms) assert.equal(same(form, 'chrX:g.20148675dup', eif1axReadme.reference), true, form);
}

// ── Same allele, different spelling ────────────────────────────────────────
assert.equal(same(eif1axDup.query, eif1axDup.resolved, eif1axDup.reference), true, 'reported variant: ins ≡ dup');
assert.equal(same(eif1axReadme.query, eif1axReadme.resolved, eif1axReadme.reference), true, 'README example: ins ≡ dup');
// Ensembl shifts CFTR F508del 3' out of its VCF-minimal position.
assert.equal(same('chr7:g.117199645_117199647del', 'chr7:g.117199646_117199648del', cftr.reference), true, 'F508del shifts 3\'');
// A delins that restates an unchanged base is the substitution inside it.
assert.equal(same('chr7:g.117199644_117199645delinsAG', 'chr7:g.117199645T>G', cftr.reference), true);
assert.equal(same('chr7:g.117199645T>G', 'chr7:g.117199645T>G', cftr.reference), true);

// ── Different allele, similar spelling ─────────────────────────────────────
assert.equal(same(brca2Ins.query, brca2Ins.neighbourId, brca2Ins.reference), false, 'BRCA2 insCA one base 3\' is another variant');
assert.equal(same(eif1axDup.query, 'chrX:g.20156694_20156695insAGATTCATTCTCATTCTTACCCCTGCGTCT', eif1axDup.reference), false,
    'the off-by-one candidate the old conversion produced for the reported variant');
// One base 5' of F508del is I507del — the textbook pair of neighbouring CF alleles.
assert.equal(same('chr7:g.117199645_117199647del', 'chr7:g.117199644_117199646del', cftr.reference), false, 'F508del ≠ I507del');
assert.equal(same('chr7:g.117199645T>G', 'chr7:g.117199645T>C', cftr.reference), false);
assert.equal(same('chr7:g.117199645T>G', 'chr8:g.117199645T>G', cftr.reference), false, 'different chromosome');

// ── "Cannot tell" is never "different" ─────────────────────────────────────
assert.equal(same('chr7:g.117199645A>G', 'chr7:g.117199645T>G', cftr.reference), null, 'stated REF contradicts the reference');
assert.equal(same('chr7:g.117199645_117199647delAAA', 'chr7:g.117199646_117199648del', cftr.reference), null, 'spelled-out deletion is not what is there');
assert.equal(same('chr7:g.117199645_117199647del', 'chr7:g.117209646_117209648del', cftr.reference), null, 'outside the window');
assert.equal(same('chr7:g.117199645_117199647', 'chr7:g.117199646_117199648del', cftr.reference), null, 'a bare span is a locus, not an allele');
assert.equal(same('CFTR:c.1521_1523del', 'chr7:g.117199646_117199648del', cftr.reference), null, 'not genomic');
assert.equal(h.applyGenomicHgvsToReference(h.parseGenomicHgvs('chr7:g.117199646_117199648del'), cftr.reference.start, cftr.reference.seq),
    cftr.reference.seq.slice(0, 96) + cftr.reference.seq.slice(99), 'deletion removes exactly the named bases');

// ── compareGenomicAlleles: the async wrapper around the reference fetch ────
{
    const calls = [];
    const referenceFor = fx => async (chrom, start, end) => {
        calls.push({ chrom, start, end });
        assert.equal(chrom.replace(/^chr/i, ''), fx.reference.chrom);
        return referenceSlice(fx.reference, start, end);
    };
    const compare = (fx, resolved, fetchReference = referenceFor(fx)) => h.compareGenomicAlleles(fx.query, resolved, fetchReference);

    assert.deepEqual(await compare(eif1axDup, eif1axDup.resolved), { status: 'same', detail: '' });
    assert.equal(calls.length, 1, 'one reference request');
    assert(calls[0].start < 20156693 && calls[0].end > 20156723, 'window covers both spellings');
    assert((await compare(eif1axDup, 'NC_000023.10:g.20156694_20156723dup')).status === 'same', 'NC_ GRCh37 form is compared too');
    assert.equal((await compare(eif1axReadme, eif1axReadme.resolved)).status, 'same');
    assert.equal((await compare(brca2Ins, brca2Ins.neighbourId)).status, 'different');

    // Nothing to check → no request at all.
    calls.length = 0;
    assert.equal((await compare(eif1axDup, eif1axDup.query)).status, 'not-needed');
    assert.equal((await compare(eif1axDup, eif1axDup.query.toUpperCase().replace('CHRX', 'chrX'))).status, 'not-needed');
    assert.equal((await h.compareGenomicAlleles('EIF1AX:c.34_63dup', eif1axDup.resolved, referenceFor(eif1axDup))).status, 'not-needed');
    assert.equal((await compare(eif1axDup, 'NC_000023.11:g.20138575_20138604dup')).status, 'not-needed', 'GRCh38 is never compared as hg19');
    // How the chromosome is NAMED is not a re-spelling of the variant.
    assert.equal((await h.compareGenomicAlleles('7:g.117199645_117199647del', 'chr7:g.117199645_117199647del', referenceFor(cftr))).status, 'not-needed');
    assert.equal((await h.compareGenomicAlleles(mito.query, mito.resolved, referenceFor(cftr))).status, 'not-needed', 'chrM is chrMT');
    // An LRG_/NG_ coordinate gives a chromosomal variant nothing to be held to.
    assert.equal((await h.compareGenomicAlleles('LRG_663:g.139821_139823del', 'chr7:g.117199646_117199648del', referenceFor(cftr))).status, 'not-needed');
    // Verdicts that need no sequence either.
    assert.equal((await compare(eif1axDup, 'chrX:g.20256694_20256723dup')).status, 'different', '100 kb away is another locus');
    assert.equal((await compare(eif1axDup, 'chr7:g.20156694_20156723dup')).status, 'different');
    assert.equal((await h.compareGenomicAlleles('chrX:g.20156693_20196693del', 'chrX:g.20156694_20196694del', referenceFor(eif1axDup))).status,
        'unverified', 'a 40 kb event is not compared base by base');
    assert.equal(calls.length, 0, 'none of the above may fetch the reference');
    // ...but a chr-less input that WAS re-spelled is still checked.
    assert.equal((await h.compareGenomicAlleles('7:g.117199645_117199647del', 'chr7:g.117199646_117199648del', referenceFor(cftr))).status, 'same');
    assert.equal((await h.compareGenomicAlleles('7:g.117199645_117199647del', 'chr7:g.117199644_117199646del', referenceFor(cftr))).status, 'different');

    // An unreachable or unusable reference downgrades to "unverified"; it never throws.
    const warn = console.warn;
    console.warn = () => {}; // the helper logs the outage; keep the test output clean
    try {
        assert.deepEqual(await compare(eif1axDup, eif1axDup.resolved, async () => { throw new Error('Ensembl sequence request failed (503)'); }),
            { status: 'unverified', detail: 'reference sequence unavailable' });
    } finally {
        console.warn = warn;
    }
    assert.equal((await compare(eif1axDup, eif1axDup.resolved, async () => null)).status, 'unverified');
    assert.equal((await compare(eif1axDup, eif1axDup.resolved, async () => 'ACGT')).status, 'unverified', 'short payload');
}

// ── The spelling-level check ───────────────────────────────────────────────
// The placeholder built when MyVariant has no record is filed under the query.
// Before the fix this threw for every re-spelled indel.
h.assertResolvedVariantIdentity(eif1axDup.query, eif1axDup.resolved, { _id: eif1axDup.query });
h.assertResolvedVariantIdentity(eif1axReadme.query, eif1axReadme.resolved, { _id: eif1axReadme.query, dbnsfp: { genename: 'EIF1AX' } });
h.assertResolvedVariantIdentity(eif1axDup.query, eif1axDup.resolved, { _id: eif1axDup.resolved });
h.assertResolvedVariantIdentity(eif1axDup.query, eif1axDup.resolved, null);
// ...and it still refuses an annotation filed under neither.
assert.throws(() => h.assertResolvedVariantIdentity(eif1axDup.query, eif1axDup.resolved, { _id: 'chrX:g.20156695_20156724dup' }),
    /annotation belongs to a different genomic allele/);
assert.throws(() => h.assertResolvedVariantIdentity('EIF1AX:c.34_63dup', eif1axDup.resolved, { _id: eif1axDup.query }),
    /annotation belongs to a different genomic allele/, 'a transcript query vouches for no genomic spelling');
assert.throws(() => h.assertResolvedVariantIdentity('chr7:g.117199645T>G', 'chr7:g.117199645T>C', { _id: 'chr7:g.117199645T>G' }),
    /does not match the genomic input/, 'substitutions stay exact');
// One chromosome, two names: UCSC's chrM is the chrMT of Ensembl and MyVariant.
// Before the fix "M 8993 T G" was refused as a different allele from itself.
h.assertResolvedVariantIdentity(mito.query, mito.resolved, { _id: mito.resolved });
assert.equal(h.selectMyVariantHit([mito.record], mito.query), mito.record);
assert.throws(() => h.assertResolvedVariantIdentity(mito.query, 'chrMT:g.8993T>C', { _id: 'chrMT:g.8993T>C' }),
    /does not match the genomic input/, 'the alleles still have to agree');
assert.throws(() => h.assertResolvedVariantIdentity(mito.query, 'chr1:g.8993T>G', { _id: 'chr1:g.8993T>G' }),
    /does not match the genomic input/, 'only M and MT are the same chromosome');
console.log('Indel identity unit checks passed');

// ── Full lookups through the real form handler ─────────────────────────────
const script = readFileSync(new URL('../script.js', import.meta.url), 'utf8');
const html = readFileSync(new URL('../index.html', import.meta.url), 'utf8');

// `fx` supplies the recorded upstream answers. `myvariant` maps ids MyVariant
// holds to their records (everything else is a 404, as it was live).
async function lookup(input, fx, { myvariant = {}, recoder = fx.recoder, referenceDown = false } = {}) {
    const { document } = parseHTML(html);
    const form = document.getElementById('variantForm');
    const listen = form.addEventListener.bind(form);
    let submit;
    form.addEventListener = (event, fn, ...args) => {
        if (event === 'submit') submit = fn;
        return listen(event, fn, ...args);
    };
    const requests = [], errors = [], timers = [];
    const response = (data, status = 200) => ({ ok: status === 200, status,
        json: async () => structuredClone(data),
        text: async () => (typeof data === 'string' ? data : JSON.stringify(data)) });
    const fetch = async address => {
        const url = new URL(address);
        const path = decodeURIComponent(url.pathname);
        requests.push(`${url.hostname}${path}`);
        if (url.hostname === 'myvariant.info') {
            if (!path.startsWith('/v1/variant/')) return response({ hits: [] });
            const id = path.slice('/v1/variant/'.length);
            return myvariant[id] ? response(myvariant[id]) : response({ success: false, error: 'Not Found.' }, 404);
        }
        if (path.includes('/variant_recoder/')) {
            return path.endsWith(`/${fx.query}`) ? response(recoder) : response({ error: 'unexpected recoder query' }, 400);
        }
        if (path.includes('/vep/human/hgvs/')) {
            return path.endsWith(`/${fx.vep[0].input}`) ? response(fx.vep) : response({ error: 'unexpected VEP query' }, 400);
        }
        if (path.includes('/sequence/region/human/')) {
            if (referenceDown) return response('Service Unavailable', 503);
            const [, chrom, start, end] = path.match(/human\/([^:]+):(\d+)\.\.(\d+)$/);
            assert.equal(chrom, fx.reference.chrom);
            return response(referenceSlice(fx.reference, Number(start), Number(end)));
        }
        if (path.includes('/lookup/symbol/')) return response(fx.mane);
        // Evidence cards fail soft; they have no say in which variant this is.
        return response({}, 404);
    };
    const window = { location: new URL('https://drdoubleb.com/variant/'),
        history: { replaceState() {} }, matchMedia: () => ({ matches: false }) };
    const context = vm.createContext({ document, window, fetch,
        console: { log() {}, warn() {}, error(e) { errors.push(String(e)); } },
        URL, URLSearchParams, AbortController, performance,
        setTimeout(fn, ms) { const timer = setTimeout(fn, Math.min(ms, 5)); timers.push(timer); return timer; },
        clearTimeout
    });
    try {
        vm.runInContext(script, context);
        vm.runInContext('initVariantSearchUi()', context);
        document.getElementById('variantInput').value = input;
        await submit({ preventDefault() {} });
        await new Promise(resolve => setTimeout(resolve, 15));
        const card = name => [...document.querySelectorAll('.card')].find(el => el.querySelector('h3')?.textContent === name);
        return {
            status: document.getElementById('status').textContent,
            hidden: document.getElementById('resultSection').classList.contains('hidden'),
            panel: document.getElementById('lookupProgress').textContent,
            variant: card('Variant')?.textContent || '',
            clinvar: card('ClinVar')?.textContent || '',
            requests, errors
        };
    } finally {
        timers.forEach(clearTimeout);
    }
}
const resolved = (page, label) => {
    assert.equal(page.status, '', `${label}: ${page.errors.join('\n')}`);
    assert.equal(page.hidden, false, `${label}: result rendered`);
};

// 1. The reported row. Was: "Variant identity mismatch".
{
    const page = await lookup(eif1axDup.row, eif1axDup);
    resolved(page, 'reported row');
    assert.match(page.variant, /g\.: chrX:g\.20156694_20156723dup/);
    assert.match(page.variant, /X-20156693-C-CAGATTCATTCTCATTCTTACCCCTGCGTCT/, 'VCF alleles are the ones pasted');
    assert.match(page.variant, /Gene: EIF1AX/);
    assert.match(page.variant, /Transcript: NM_001412\.4 \(MANE Select\)/);
    assert.match(page.variant, /c\.: c\.34_63dup/);
    assert.match(page.variant, /p\.: p\.Arg12_Ser21dup/);
    assert.match(page.variant, /inframe_insertion/);
    assert.match(page.panel, /Confirm allele · GRCh37 reference/);
    assert.match(page.panel, /chrX:g\.20156694_20156723dup is the same allele as the input/);
    assert(!page.requests.some(r => r.includes('20156694_20156695ins')), 'the one-base-off candidate is never formed');
    assert.equal(page.requests.filter(r => r.includes('/sequence/region/')).length, 1, 'one reference request');
    console.log('ok indel lookup: reported EIF1AX 30-bp duplication');
}

// 2. Same variant typed as HGVS, and with the reference endpoint down: it must
//    still resolve (Ensembl's own mapping is trusted) but claim no confirmation.
{
    const page = await lookup(eif1axDup.query, eif1axDup, { referenceDown: true });
    resolved(page, 'reference down');
    assert.match(page.variant, /g\.: chrX:g\.20156694_20156723dup/);
    assert.match(page.variant, /c\.: c\.34_63dup/);
    assert.match(page.panel, /not confirmed — reference sequence unavailable/);
    assert.doesNotMatch(page.panel, /same allele as the input/);
    console.log('ok indel lookup: reference unavailable');
}

// 3. The README's documented example. Was: "Variant identity mismatch".
{
    const page = await lookup(eif1axReadme.row, eif1axReadme);
    resolved(page, 'README row');
    assert.match(page.variant, /g\.: chrX:g\.20148675dup/);
    assert.match(page.variant, /X-20148674-T-TG/);
    assert.match(page.variant, /c\.: c\.388dup/);
    assert.match(page.variant, /p\.: p\.Gln130ProfsTer3/);
    assert.match(page.panel, /chrX:g\.20148675dup is the same allele as the input/);
    console.log('ok indel lookup: README example X 20148674 EIF1AX T TG');
}

// 4. BRCA2 insCA with the neighbouring ClinVar record present in MyVariant,
//    exactly as it is live. Was: that neighbour, presented as the answer.
{
    const page = await lookup(brca2Ins.row, brca2Ins, { myvariant: { [brca2Ins.neighbourId]: brca2Ins.neighbour } });
    resolved(page, 'BRCA2 row');
    assert.match(page.variant, /g\.: chr13:g\.32906698_32906699insCA/);
    assert.match(page.variant, /c\.: c\.1083_1084insCA/);
    assert.match(page.variant, /p\.: p\.Asp362GlnfsTer6/);
    assert.doesNotMatch(page.variant, /32906699_32906700|1084_1085|Asp362Ala/, 'nothing from the neighbouring allele');
    assert.doesNotMatch(page.clinvar, /2505518|Likely pathogenic/, 'the neighbour\'s ClinVar record is not this variant\'s');
    assert(!page.requests.some(r => r.includes(brca2Ins.neighbourId)), 'the neighbouring allele is never even requested');
    assert.doesNotMatch(page.panel, /Confirm allele/, 'spelling unchanged — nothing to confirm');
    console.log('ok indel lookup: BRCA2 insertion keeps its own allele');
}

// 5. Backstop. If a different allele is ever adopted for a genomic input again
//    — here the recoder is made to answer with the neighbour — the lookup has
//    to stop rather than present it. Indels used to be exempt from this check.
{
    const tampered = structuredClone(brca2Ins.recoder);
    const allele = Object.values(tampered[0])[0];
    allele.hgvsg = ['NC_000013.10:g.32906699_32906700insCA'];
    allele.spdi = [];
    const page = await lookup(brca2Ins.row, brca2Ins, {
        recoder: tampered, myvariant: { [brca2Ins.neighbourId]: brca2Ins.neighbour }
    });
    assert.equal(page.status, 'Error: Variant identity mismatch: the resolved allele does not match the genomic input.');
    assert.equal(page.hidden, true, 'a different allele is never rendered');
    assert.match(page.panel, /chr13:g\.32906699_32906700insCA is not the allele entered/);
    console.log('ok indel lookup: a provably different allele is refused');
}

// 6. m.8993T>G (NARP/Leigh, ClinVar 9641) pasted with the hg19 chromosome name.
//    MyVariant and Ensembl call it MT. Was: "Variant identity mismatch".
{
    const page = await lookup(mito.row, mito, { myvariant: { [mito.resolved]: mito.record } });
    resolved(page, 'mitochondrial row');
    assert.match(page.variant, /g\.: chrMT:g\.8993T>G/);
    assert.match(page.clinvar, /Variation ID: 9641/);
    assert.doesNotMatch(page.panel, /Confirm allele/, 'a chromosome alias is not a re-spelling');
    assert(!page.requests.some(r => r.includes('/sequence/region/')), 'so it costs no reference request');
    console.log('ok indel lookup: chrM input resolves to its chrMT record');
}
