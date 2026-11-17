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

// ── Which variant a gene-symbol name means ─────────────────────────────────
// MyVariant's hits for "BRCA1 c.5266dup", in the order it ranks them. Each is
// c.5266dup on SOME RefSeq transcript; only ClinVar 17677 is c.5266dup on the
// transcript ClinVar names it by (5382insC). It ranks fourth, and the first
// hit used to be the answer.
{
    const { brca1Dup } = fixture;
    const hits = brca1Dup.searchHits;
    assert.deepEqual(hits.map(hit => hit.clinvar.variant_id), [266551, 55618, 2583814, 17677, 548237, 55526], 'fixture keeps MyVariant\'s ranking');
    for (const hit of hits) assert(h.hitMatchesVariantQuery(hit, 'BRCA1 c.5266dup'), `${hit._id} matches on some transcript`);
    for (const name of ['BRCA1 c.5266dup', 'BRCA1:c.5266dup', 'BRCA1:c.5266dupC']) {
        assert.equal(h.selectMyVariantHit(hits, name)?.clinvar.variant_id, 17677, name);
    }
    // An explicit accession was always unambiguous and still is.
    assert.equal(h.selectMyVariantHit(hits, 'NM_007294.4:c.5266dup').clinvar.variant_id, 17677);
    assert.equal(h.selectMyVariantHit(hits, 'NM_007294.4:c.5389dup').clinvar.variant_id, 266551);
    assert.equal(h.selectMyVariantHit(hits, 'NM_001407664.1:c.5266dup').clinvar.variant_id, 266551, 'the alternative transcript is honoured when named');
    // Without the variant ClinVar names this way, the name still resolves — a
    // preference, not a filter — on the lowest-numbered NM_ it matched.
    const others = hits.filter(hit => hit.clinvar.variant_id !== 17677);
    const fallback = h.selectMyVariantHit(others, 'BRCA1 c.5266dup');
    assert(fallback && h.variantNameRank(fallback, 'BRCA1 c.5266dup')[0] === 0);
    assert.equal(h.variantNameRank(hits[3], 'BRCA1 c.5266dup')[0], 1);
    assert.equal(h.variantNameRank(hits[3], 'BRCA1:c.5266dupC')[0], 1, 'trailing bases do not change the name');
    assert.deepEqual(h.variantNameRank(hits[0], 'chr17:g.41201154_41201155dup'), [0, 0, 0], 'genomic queries are not ranked by name');
    // The same for protein names: an isoform-shifted p. match must not outrank
    // the variant ClinVar calls by that name.
    const shifted = { _id: 'chr7:g.140453256A>T', dbnsfp: { genename: 'BRAF', hgvsp: ['p.Val600Glu', 'p.Val560Glu'] },
        clinvar: { gene: { symbol: 'BRAF' }, rcv: { preferred_name: 'NM_004333.6(BRAF):c.1679T>A (p.Val560Glu)' } } };
    const v600e = { _id: 'chr7:g.140453136A>T', dbnsfp: { genename: 'BRAF', hgvsp: ['p.Val640Glu', 'p.Val600Glu'] },
        clinvar: { gene: { symbol: 'BRAF' }, rcv: [{ preferred_name: 'NM_004333.6(BRAF):c.1799T>A (p.Val600Glu)' }] } };
    assert.equal(h.selectMyVariantHit([shifted, v600e], 'BRAF:p.V600E'), v600e);
    assert.equal(h.selectMyVariantHit([shifted, v600e], 'BRAF:p.Val560Glu'), shifted);
}

// ── A MyVariant id that contradicts the record's own ClinVar HGVS ──────────
{
    const { brca1Dup, egfrDup } = fixture;
    const compareWith = fx => (a, b) => h.compareGenomicAlleles(a, b, async (chrom, start, end) => referenceSlice(fx.reference, start, end));
    const neverCompared = () => { throw new Error('nothing here needs the reference'); };

    // 5382insC: filed under the insertion flanks + "dup", which as HGVS is a
    // 2-bp duplication. ClinVar's own HGVS in the same record is the 1-bp one.
    assert.equal(h.clinvarGenomicAllele(brca1Dup.record), brca1Dup.resolved);
    assert.equal(same(brca1Dup.recordId, brca1Dup.resolved, brca1Dup.reference), false, 'the id is not the allele ClinVar describes');
    assert.equal(h.recordFoundThroughClinvar(brca1Dup.record, 'BRCA1 c.5266dup'), true);
    {
        const out = await h.reconcileClinvarAllele(brca1Dup.record, true, compareWith(brca1Dup));
        assert.equal(out.status, 'rekeyed');
        assert.equal(out.record._id, 'chr17:g.41209082dup');
        assert.equal(out.record._myvariant_id, brca1Dup.recordId, 'provenance kept');
        assert.equal(out.record.clinvar.variant_id, 17677);
        assert.equal(out.record.dbnsfp.genename, 'BRCA1');
        assert.equal(out.record.dbsnp, undefined, 'the dbSNP block is rs2051500205, the literal 2-bp dup — not this variant');
    }
    {
        // Looked up as the id's own allele, the ClinVar entry is someone else's.
        const out = await h.reconcileClinvarAllele(brca1Dup.record, false, compareWith(brca1Dup));
        assert.equal(out.status, 'clinvar-dropped');
        assert.equal(out.record.clinvar, undefined);
        assert.equal(out.record.dbsnp.rsid, 'rs2051500205', 'what was merged in by id stays');
        assert.equal(out.record._id, brca1Dup.recordId);
    }
    // EGFR exon 20 duplication: ClinVar-only record, and MyVariant's ClinVar
    // block names the antisense gene.
    {
        assert.equal(egfrDup.record.clinvar.gene.symbol, 'EGFR-AS1');
        const out = await h.reconcileClinvarAllele(egfrDup.record, true, compareWith(egfrDup));
        assert.equal(out.record._id, 'chr7:g.55249002_55249010dup');
        assert.equal(out.record.dbnsfp.genename, 'EGFR', 'gene taken from ClinVar\'s preferred name');
        const literal = await h.reconcileClinvarAllele(egfrDup.record, false, compareWith(egfrDup));
        assert.equal(literal.status, 'clinvar-dropped');
        assert.equal(literal.record, null, 'nothing else was filed under that id');
    }
    // BRCA2 insCA: a record with its own vcf block, one base off ClinVar 2505518.
    {
        const record = brca2Ins.neighbour;
        assert.equal(h.clinvarGenomicAllele(record), 'chr13:g.32906700_32906701insCA');
        assert.equal(h.recordFoundThroughClinvar(record, 'NM_000059.4:c.1085_1086insCA'), true);
        assert.equal(h.recordFoundThroughClinvar(record, 'NM_000059.4:c.1084_1085insCA'), false);
        const literal = await h.reconcileClinvarAllele(record, false, compareWith(brca2Ins));
        assert.equal(literal.status, 'clinvar-dropped');
        assert.deepEqual(literal.record.vcf, record.vcf, 'its own annotations stay');
        assert.equal(literal.record.clinvar, undefined);
        const named = await h.reconcileClinvarAllele(record, true, compareWith(brca2Ins));
        assert.equal(named.record._id, 'chr13:g.32906700_32906701insCA');
        assert.equal(named.record.vcf, undefined, 'the vcf block is the other allele\'s');
        assert.equal(named.record.clinvar.variant_id, 2505518);
    }
    // Consistent records are returned untouched.
    {
        const sameAllele = { _id: 'chrX:g.20148674_20148675insG', vcf: { ref: 'T', alt: 'TG' },
            clinvar: { variant_id: 1, hgvs: { genomic: ['NC_000023.11:g.20130557dup', 'NC_000023.10:g.20148675dup'] } } };
        const out = await h.reconcileClinvarAllele(sameAllele, true, compareWith(eif1axReadme));
        assert.equal(out.status, 'consistent');
        assert.equal(out.record, sameAllele);
        assert.equal(h.clinvarGenomicAllele(sameAllele), 'chrX:g.20148675dup', 'GRCh38 entries are skipped');
        // Deletions and substitutions cannot go wrong this way and cost nothing.
        const del = { _id: 'chr7:g.117199645_117199647del', clinvar: { hgvs: { genomic: ['NC_000007.13:g.117199646_117199648del'] } } };
        assert.equal((await h.reconcileClinvarAllele(del, true, neverCompared)).record, del);
        const snv = { _id: 'chr7:g.140453136A>T', clinvar: { hgvs: { genomic: ['NC_000007.13:g.140453136A>T'] } } };
        assert.equal((await h.reconcileClinvarAllele(snv, true, neverCompared)).record, snv);
        assert.equal((await h.reconcileClinvarAllele({ _id: 'chrX:g.20148674_20148675insG' }, true, neverCompared)).status, 'consistent', 'no ClinVar block');
        // No reference, no verdict: left exactly as it was.
        const unverified = await h.reconcileClinvarAllele(brca1Dup.record, true, async () => ({ status: 'unverified', detail: 'reference sequence unavailable' }));
        assert.equal(unverified.status, 'unverified');
        assert.equal(unverified.record, brca1Dup.record);
    }
}
// ── Reference bases the input states ───────────────────────────────────────
{
    const stated = h.statedReferenceBases;
    assert.deepEqual(stated('16 2122946 TSC2 TTAA T', 'chr16:g.2122947_2122949del'), { chrom: 'chr16', start: 2122946, bases: 'TTAA' });
    assert.deepEqual(stated('X\t20148674\tEIF1AX\tT\tTG', 'chrX:g.20148674_20148675insG'), { chrom: 'chrX', start: 20148674, bases: 'T' });
    assert.deepEqual(stated('7 140453136 BRAF AC TT', 'chr7:g.140453136_140453137delinsTT'), { chrom: 'chr7', start: 140453136, bases: 'AC' });
    assert.deepEqual(stated('chr7:g.117199645_117199647delAAA', 'chr7:g.117199645_117199647delAAA'), { chrom: 'chr7', start: 117199645, bases: 'AAA' });
    // Nothing stated, or a substitution (its REF is in the MyVariant id, and
    // Ensembl rejects a wrong one itself).
    assert.equal(stated('chr7:g.117199645_117199647del', 'chr7:g.117199645_117199647del'), null);
    assert.equal(stated('7 140453136 BRAF A T', 'chr7:g.140453136A>T'), null);
    assert.equal(stated('chr7:g.140453136A>T', 'chr7:g.140453136A>T'), null);
    assert.equal(stated('X 20148674 EIF1AX - G', 'chrX:g.20148674_20148675insG'), null, 'MAF "-" states no bases');
    assert.equal(stated('BRCA1 c.5266dup', 'BRCA1:c.5266dup'), null);

    const reference = async (chrom, start, end) => referenceSlice(cftr.reference, start, end);
    assert.deepEqual(await h.verifyStatedReference({ chrom: 'chr7', start: 117199645, bases: 'TCT' }, reference), { status: 'matches', actual: 'TCT' });
    assert.deepEqual(await h.verifyStatedReference({ chrom: 'chr7', start: 117199645, bases: 'AAA' }, reference), { status: 'mismatch', actual: 'TCT' });
    assert.equal((await h.verifyStatedReference({ chrom: 'chr7', start: 117199645, bases: 'TNT' }, reference)).status, 'unverified');
    assert.equal((await h.verifyStatedReference({ chrom: 'chr7', start: 117199645, bases: 'TCT' }, async () => null)).status, 'unverified');
    const warn = console.warn;
    console.warn = () => {};
    try {
        assert.equal((await h.verifyStatedReference({ chrom: 'chr7', start: 117199645, bases: 'TCT' },
            async () => { throw new Error('Ensembl sequence request failed (503)'); })).status, 'unverified', 'an outage is not a mismatch');
    } finally {
        console.warn = warn;
    }
}

// ── Who may answer, and how an answer is read ──────────────────────────────
// A chromosomal coordinate is only ever asked of the GRCh37 mirror: the main
// host would read the same numbers as GRCh38.
assert.deepEqual(h.recoderHostsFor('chr9:g.21971058del'), ['https://grch37.rest.ensembl.org']);
assert.deepEqual(h.recoderHostsFor('chrX:g.20156693_20156694insAGATT'), ['https://grch37.rest.ensembl.org']);
for (const assemblyFree of ['NM_000077.5:c.300del', 'CDKN2A:c.300del', 'BRAF:p.Val600Glu', 'rs80357906', 'chrMT:g.8993T>G', 'chrM:g.8993T>G']) {
    assert.equal(h.recoderHostsFor(assemblyFree).length, 2, `${assemblyFree} still fails over`);
}
// Ensembl writes the mitochondrial genome "m.", not "g.".
assert.equal(h.convertHgvsgToMyVariant('NC_012920.1:m.8993T>G'), 'chrMT:g.8993T>G');
assert.equal(h.convertHgvsgToMyVariant('NC_000023.10:g.20156694_20156723dup'), 'chrX:g.20156694_20156723dup');
assert.throws(() => h.convertHgvsgToMyVariant('LRG_292:g.160921dup'), /Cannot map accession/);
assert.throws(() => h.convertHgvsgToMyVariant('NC_012920.1'), /Cannot read genomic HGVS/);
console.log('Indel identity unit checks passed');

// ── Full lookups through the real form handler ─────────────────────────────
const script = readFileSync(new URL('../script.js', import.meta.url), 'utf8');
const html = readFileSync(new URL('../index.html', import.meta.url), 'utf8');

// `fx` supplies the recorded upstream answers. `myvariant` maps ids MyVariant
// holds to their records (everything else is a 404, as it was live); `search`
// answers its free-text endpoint; a null `recoder` is Ensembl not answering;
// `route` may answer any request first (return undefined to fall through).
async function lookup(input, fx, {
    myvariant = {}, search = () => [], recoder = fx.recoder, referenceDown = false, route = () => undefined, settleMs = 15
} = {}) {
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
        requests.push(`${url.hostname}${path}${decodeURIComponent(url.search)}`);
        const routed = route(url, path, response);
        if (routed) return routed;
        if (url.hostname === 'myvariant.info') {
            if (!path.startsWith('/v1/variant/')) return response({ hits: search(url.searchParams.get('q')) || [] });
            const id = path.slice('/v1/variant/'.length);
            return myvariant[id] ? response(myvariant[id]) : response({ success: false, error: 'Not Found.' }, 404);
        }
        if (path.includes('/variant_recoder/')) {
            if (!recoder) return response({ error: 'Service Unavailable' }, 503);
            // Ensembl is asked about the mitochondrial chromosome as "MT:".
            const asked = fx.query.replace(/^chrMT?:g\./i, 'MT:g.');
            return path.endsWith(`/${asked}`) ? response(recoder) : response({ error: 'unexpected recoder query' }, 400);
        }
        if (path.includes('/vep/human/hgvs/')) {
            return path.endsWith(`/${fx.vep[0].input}`) ? response(fx.vep) : response({ error: 'unexpected VEP query' }, 400);
        }
        if (path.includes('/sequence/region/human/')) {
            if (referenceDown) return response('Service Unavailable', 503);
            const [, chrom, start, end] = path.match(/human\/([^:]+):(\d+)\.\.(\d+)$/);
            assert.equal(chrom, fx.reference.chrom);
            // The wide window read to find an allele's other MyVariant ids is
            // not in these fixtures; that step then reports "not checked".
            if (Number(end) - Number(start) > 250) return response('Service Unavailable', 503);
            return response(referenceSlice(fx.reference, Number(start), Number(end)));
        }
        // The MANE lookup is two requests: the gene (for its canonical
        // transcript), then that transcript with its MANE record. Both are
        // answered from the one gene record these fixtures hold.
        const maneTranscript = (fx.mane.Transcript || []).find(t => [].concat(t.MANE || []).some(m => m.type === 'MANE_Select'));
        if (path.includes('/lookup/symbol/')) {
            return maneTranscript ? response({ canonical_transcript: `${maneTranscript.id}.1` }) : response({}, 404);
        }
        if (path.includes('/lookup/id/')) {
            return maneTranscript && path.endsWith(`/${maneTranscript.id}`) ? response(maneTranscript) : response({}, 404);
        }
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
        await new Promise(resolve => setTimeout(resolve, settleMs));
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
    // Three reference reads: the row's REF base and the stretch both spellings
    // are applied to (both started early), and the wider one that shows which
    // other ids MyVariant may file this allele under.
    assert.deepEqual(page.requests.filter(r => r.includes('/sequence/region/')).map(r => r.match(/human\/([^?]+)/)[1]).sort(),
        ['X:20156544..20156873', 'X:20156691..20156725', 'X:20156693..20156693']);
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
    // HGVS g. input carries no gene. What precedes ":g." is a chromosome, not a
    // symbol: this used to render "Gene: CHRX", lose the protein change and
    // effect (the transcript pairing filtered on "CHRX"), and send "CHRX" to
    // every gene-level source. VEP names the gene.
    assert.match(page.variant, /Gene: EIF1AX/);
    assert.match(page.variant, /Transcript: NM_001412\.4 \(MANE Select\)/);
    assert.match(page.variant, /p\.: p\.Arg12_Ser21dup/);
    assert.match(page.variant, /inframe_insertion/);
    assert.doesNotMatch(page.variant, /Gene: CHR/i);
    assert.match(page.panel, /GeneEIF1AX/, 'and the panel names it too');
    assert(page.requests.some(r => /[?&]gene=EIF1AX\b/.test(r)), 'gene-level cards ask about EIF1AX');
    assert(!page.requests.some(r => /[?&](?:gene|term)=CHR/i.test(r)), 'never about "CHRX"');
    console.log('ok indel lookup: reference unavailable, gene named by VEP');
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

// 7. "BRCA1 c.5266dup" — 5382insC. Was: ClinVar 266551 (c.5389dup on the
//    standard transcript), because it ranks first of six hits that are
//    c.5266dup on some transcript. And once the right record is chosen, its
//    MyVariant id is the insertion flanks + "dup": read as HGVS that gave VCF
//    17-41209078-C-CTG, a 2-bp insertion, for gnomAD and SpliceAI.
//    Ensembl does not answer a gene-symbol BRCA1 query in time; it is down here.
{
    const { brca1Dup } = fixture;
    const records = Object.fromEntries(brca1Dup.searchHits.map(hit => [hit._id, hit]));
    records[brca1Dup.recordId] = brca1Dup.record;
    const page = await lookup(brca1Dup.input, brca1Dup, {
        recoder: null, myvariant: records, settleMs: 80,
        search: q => (/^BRCA1[ :]c\.5266dup$/.test(q) ? brca1Dup.searchHits : [])
    });
    resolved(page, 'BRCA1 c.5266dup');
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.doesNotMatch(page.clinvar + page.variant, /266551|5389dup/);
    assert(!page.requests.some(r => r.includes('41201154')), 'the first-ranked hit is not even fetched');
    assert.match(page.variant, /g\.: chr17:g\.41209082dup/, 'ClinVar\'s own HGVS, not MyVariant\'s id');
    assert.doesNotMatch(page.variant, /41209079_41209080dup/);
    assert.match(page.variant, /VCF \(hg19\): 17-41209079-T-TG/);
    assert.match(page.variant, /Gene: BRCA1/);
    assert.match(page.variant, /Transcript: NM_007294\.4 \(MANE Select\)/);
    assert.match(page.variant, /c\.: c\.5266dup/);
    assert.match(page.variant, /p\.: p\.Gln1756ProfsTer74/);
    assert.match(page.panel, /allele chr17:g\.41209082dup from ClinVar's own HGVS \(MyVariant files it under chr17:g\.41209079_41209080dup\)/);
    const downstream = page.requests.filter(r => /\/api\/(?:gnomad-v4|spliceai)/.test(r));
    assert(downstream.some(r => /gnomad-v4\?chrom=17&pos37=41209079&ref=T&alt=TG\b/.test(r)), downstream.join('\n'));
    assert(downstream.some(r => r.includes('variant=chr17-41209079-T-TG&')), downstream.join('\n'));
    assert(!downstream.some(r => /41209078|CTG/.test(r)), 'never the 2-bp insertion');
    console.log('ok indel lookup: BRCA1 c.5266dup is 5382insC, with ClinVar\'s allele');
}

// 8. EGFR exon 20 duplication by accession. Was: VCF 7-55248997-G-GAT (2 bp)
//    for a 9-bp duplication, and "EGFR-AS1" as the gene for every gene card.
{
    const { egfrDup } = fixture;
    const page = await lookup(egfrDup.input, egfrDup, {
        myvariant: { [egfrDup.recordId]: egfrDup.record }, settleMs: 80,
        search: q => (q === 'NM_005228.5 c.2300_2308dup' ? egfrDup.searchHits : [])
    });
    resolved(page, 'EGFR exon 20 dup');
    assert.match(page.variant, /g\.: chr7:g\.55249002_55249010dup/);
    assert.match(page.variant, /VCF \(hg19\): 7-55248998-A-ATGGCCAGCG/);
    assert.match(page.variant, /Gene: EGFR, EGFR-AS1/, 'the gene ClinVar names the variant in comes first');
    assert.match(page.variant, /Transcript: NM_005228\.5 \(MANE Select\)/);
    assert.match(page.variant, /c\.: c\.2300_2308dup/);
    assert.match(page.variant, /p\.: p\.Ala767_Val769dup/);
    assert.doesNotMatch(page.variant, /Transcript: NR_|c\.: n\./, 'not the antisense transcript');
    assert.match(page.clinvar, /Variation ID: 177678/);
    assert(page.requests.some(r => /\/api\/civic\?gene=EGFR&/.test(r)), 'gene cards ask about EGFR');
    assert(!page.requests.some(r => /[?&](?:gene|term)=EGFR-AS1/.test(r)), 'not about the antisense gene');
    assert(page.requests.some(r => r.includes('variant=chr7-55248998-A-ATGGCCAGCG&')));
    assert(!page.requests.some(r => /55248997|G-GAT/.test(r)), 'never the 2-bp insertion');
    console.log('ok indel lookup: EGFR exon 20 duplication, right allele and right gene');
}

// 9. The other direction. "13 32906699 BRCA2 G GCA" is literally the allele
//    MyVariant files ClinVar 2505518 under — and it is not that variant
//    (2505518 is c.1085_1086insCA; this is c.1084_1085insCA). Was: ClinVar
//    2505518 "Likely pathogenic" shown for it.
{
    const fx = { ...brca2Ins, ...fixture.brca2Literal };
    const page = await lookup(fx.row, fx, { myvariant: { [brca2Ins.neighbourId]: brca2Ins.neighbour } });
    resolved(page, 'literal BRCA2 row');
    assert.match(page.variant, /g\.: chr13:g\.32906699_32906700insCA/);
    assert.match(page.variant, /c\.: c\.1084_1085insCA/);
    assert.match(page.variant, /p\.: p\.Asp362AlafsTer6/);
    assert.doesNotMatch(page.clinvar, /2505518|Likely pathogenic/, 'another variant\'s ClinVar entry');
    assert.match(page.panel, /MyVariant\.infoannotation found/, 'the rest of the record is still used');
    console.log('ok indel lookup: an allele does not inherit a neighbour\'s ClinVar entry');
}

// 10. Bases the input states that GRCh37 does not hold. A pasted row whose REF
//     is wrong for this build (GRCh38 coordinates with the box unticked is the
//     usual way) used to resolve to whatever sits at those numbers in GRCh37.
{
    const wrongRef = eif1axDup.row.replace('\tC\tC', '\tG\tG');
    assert.notEqual(wrongRef, eif1axDup.row);
    const page = await lookup(wrongRef, eif1axDup);
    assert.equal(page.status, 'Error: Reference mismatch: the input gives G at chrX:20156693, but GRCh37/hg19 has C there. '
        + 'If these coordinates are GRCh38/hg38, tick the GRCh38 box and look it up again; otherwise check the position and the reference allele.');
    assert.equal(page.hidden, true, 'nothing is rendered for a variant this build does not have');
    assert.match(page.panel, /input gives G at chrX:20156693; GRCh37 has C/);
    // Said at once: the verdict is in before MyVariant's first answer, so the
    // lookup stops there instead of resolving and annotating the wrong locus
    // first (a real lookup spent nine seconds doing that).
    assert(!page.requests.some(r => r.includes('/vep/') || r.includes('/v1/query')), 'nothing more is fetched for it');
    console.log('ok indel lookup: a row whose REF is not the reference is refused');
}
// ...and bases spelled out in HGVS, which Ensembl ignores: it answers
// g.117199645_117199647delAAA as F508del although the reference there is TCT.
{
    const withBases = bases => {
        const query = `chr7:g.117199645_117199647del${bases}`;
        const recoder = structuredClone(cftr.recoder);
        Object.values(recoder[0])[0].input = query;
        return { ...cftr, query, recoder };
    };
    const wrong = await lookup('chr7:g.117199645_117199647delAAA', withBases('AAA'));
    assert.match(wrong.status, /^Error: Reference mismatch: the input gives AAA at chr7:117199645-117199647, but GRCh37\/hg19 has TCT there\./);
    assert.equal(wrong.hidden, true);
    const right = await lookup('chr7:g.117199645_117199647delTCT', withBases('TCT'));
    resolved(right, 'correctly spelled-out deletion');
    assert.match(right.variant, /g\.: chr7:g\.117199646_117199648del/);
    assert.match(right.variant, /p\.: p\.Phe508del/);
    assert.match(right.panel, /chr7:g\.117199646_117199648del is the same allele as the input/);
    console.log('ok indel lookup: spelled-out bases must be the reference\'s');
}

// 11. GRCh37 mirror down. The recoder used to fail over to rest.ensembl.org,
//     which read the hg19 coordinate as GRCh38; its SPDI was then "lifted back"
//     one base and chr9:g.21971057C>- — the deletion of a different base —
//     rendered as the result. A genomic query now fails over to nobody.
{
    const { cdkn2aFailover } = fixture;
    const fx = { query: cdkn2aFailover.query, vep: [{ input: 'unused' }], mane: {}, reference: { chrom: '9', start: 0, seq: '' } };
    const page = await lookup(cdkn2aFailover.row, fx, {
        route: (url, path, response) => {
            if (url.hostname === 'grch37.rest.ensembl.org') return response({ error: 'Service Unavailable' }, 503);
            if (url.hostname !== 'rest.ensembl.org') return undefined;
            if (path.includes('/variant_recoder/')) return response(cdkn2aFailover.grch38Recoder);
            if (path.includes('/map/human/GRCh38/')) return response(cdkn2aFailover.map);
            return undefined;
        }
    });
    assert(!page.requests.some(r => r.startsWith('rest.ensembl.org/variant_recoder/')), 'the GRCh38 host is not asked about an hg19 coordinate');
    assert.match(page.status, /^Error: Lookup could not complete: .*did not respond\. This is an upstream outage/);
    assert.equal(page.hidden, true);
    assert.doesNotMatch(page.panel + page.variant, /21971057C>-/, 'never the neighbouring base');
    console.log('ok indel lookup: mirror outage is reported, not answered from the wrong assembly');
}

// 12. A mitochondrial variant MyVariant does not hold. Ensembl answers
//     NC_012920.1:m.16300A>C; split on ":g." that became "chrMT:g.undefined",
//     which was looked up and then sent to VEP.
{
    const { mitoNovel } = fixture;
    const page = await lookup(mitoNovel.row, mitoNovel);
    resolved(page, 'novel mitochondrial row');
    assert.match(page.variant, /g\.: chrMT:g\.16300A>C/);
    assert(!page.requests.some(r => r.includes('undefined')), page.requests.filter(r => r.includes('undefined')).join('\n'));
    console.log('ok indel lookup: mitochondrial "m." answer is read as a variant');
}
