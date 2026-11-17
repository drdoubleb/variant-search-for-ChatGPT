/*
 * One allele, several MyVariant ids.
 *
 * MyVariant keys a record by the spelling its source used. An indel in a repeat
 * is spelled differently by each family of sources, so one allele is several
 * records and a lookup showed whichever it landed on:
 *
 *   rs80357906 (BRCA1 5382insC)   no ClinVar entry — ClinVar's is under another id
 *   BRCA1 c.5266dup               no population frequency — gnomAD's is under a third
 *   17 7578474 TP53 C CG          no ClinVar entry, no rsID
 *   CFTR F508del                  no rsID
 *
 * The other ids are derived from the reference, fetched in one request, and
 * each is proven to be the same allele before anything is taken from it.
 *
 * The scenarios replay answers recorded from the live services.
 *
 * Run with: node tests/sibling-records.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
// The widest recorded stretch of reference in a cassette.
const windowOf = cassette => cassette.reference.reduce((a, b) => (b.seq.length > a.seq.length ? b : a));
const recordOf = (cassette, id) => {
    const recorded = cassette.responses[`myvariant.info/v1/variant/${id}`];
    return recorded && recorded.status === 200 ? structuredClone(recorded.body) : null;
};

// ── How far an indel slides, and what that makes its ids ────────────────────
{
    // BRCA1 5382insC: a G added to GGG at 41209080-41209082.
    const { start, seq } = windowOf(cassettes.rs80357906);
    assert.equal(seq.slice(41209079 - start, 41209084 - start), 'TGGGA');
    const expected = { kind: 'ins', chrom: 'chr17', left: { pos: 41209079, seq: 'G' }, right: { pos: 41209082, seq: 'G' }, complete: true };
    // Every way of writing the allele gives the same two ends.
    for (const spelling of ['chr17:g.41209082dup', 'chr17:g.41209079_41209080insG', 'chr17:g.41209081_41209082insG',
        'chr17:g.41209082_41209083insG', 'chr17:g.41209080dup']) {
        assert.deepEqual(h.indelPlacements(spelling, start, seq), expected, spelling);
    }
    assert.deepEqual(h.siblingIdCandidates('chr17:g.41209082dup', start, seq), [
        'chr17:g.41209079_41209080insG',    // left-aligned: gnomAD, ExAC, snpEff
        'chr17:g.41209082_41209083insG',    // 3'-most, as an insertion
        'chr17:g.41209082dup',              // 3'-most, as HGVS writes it: dbSNP
        'chr17:g.41209079_41209080dup'      // MyVariant's own id for the ClinVar record
    ]);
    // A different allele next door does not share them.
    assert.deepEqual(h.indelPlacements('chr17:g.41209079_41209080insT', start, seq).right, { pos: 41209079, seq: 'T' });
    // What is not an indel, or does not fit the reference, has no siblings.
    assert.equal(h.indelPlacements('chr17:g.41209082G>A', start, seq), null);
    assert.equal(h.indelPlacements('chr17:g.41209080_41209081delinsAT', start, seq), null);
    assert.equal(h.indelPlacements('chr17:g.41209080delA', start, seq), null, 'states a base the reference does not hold');
    assert.equal(h.indelPlacements('chr17:g.99999999dup', start, seq), null, 'outside the window');
    assert.deepEqual(h.siblingIdCandidates('chr17:g.41209082G>A', start, seq), []);
    // A slide that runs off the window says so.
    const narrow = seq.slice(41209081 - start, 41209086 - start);
    assert.equal(h.indelPlacements('chr17:g.41209082dup', 41209081, narrow).complete, false);
}
{
    // CFTR F508del: CTT out of ATCTTT. Left-aligned 645-647, 3'-most 646-648.
    const { start, seq } = windowOf(cassettes.cftrF508del);
    for (const spelling of ['chr7:g.117199645_117199647del', 'chr7:g.117199646_117199648del', 'chr7:g.117199646_117199648delCTT']) {
        assert.deepEqual(h.indelPlacements(spelling, start, seq),
            { kind: 'del', chrom: 'chr7', left: { start: 117199645, end: 117199647 }, right: { start: 117199646, end: 117199648 }, complete: true }, spelling);
    }
    assert.deepEqual(h.siblingIdCandidates('chr7:g.117199646_117199648del', start, seq),
        ['chr7:g.117199645_117199647del', 'chr7:g.117199646_117199648del']);
    // Ile507del, one codon back, is another allele and another pair of ids.
    assert.equal(h.sameGenomicAllele('chr7:g.117199644_117199646del', 'chr7:g.117199645_117199647del', start, seq), false);
    assert(!h.siblingIdCandidates('chr7:g.117199644_117199646del', start, seq).includes('chr7:g.117199645_117199647del'));
}
{
    // BRCA2 c.1085_1086insCA: left-aligned it is insAC one base back. MyVariant
    // files the ClinVar record under the LEFT position with the un-shifted bases.
    const { start, seq } = windowOf(cassettes.brca2ClinVar2505518);
    assert.deepEqual(h.siblingIdCandidates('chr13:g.32906700_32906701insCA', start, seq), [
        'chr13:g.32906699_32906700insAC', 'chr13:g.32906700_32906701insCA',
        'chr13:g.32906699_32906700insCA', 'chr13:g.32906699_32906700dup'
    ]);
    // An insertion that cannot slide has one spelling (and the "dup" id).
    assert.deepEqual(h.siblingIdCandidates('chr13:g.32906699_32906700insCA', start, seq),
        ['chr13:g.32906699_32906700insCA', 'chr13:g.32906699_32906700dup']);
}

// ── What a record under another id may contribute ──────────────────────────
{
    const { start, seq } = windowOf(cassettes.brca1Dup);
    const allele = 'chr17:g.41209082dup';
    // The left-aligned spelling is this allele: everything it holds applies.
    const leftAligned = recordOf(cassettes.brca1Dup, 'chr17:g.41209079_41209080insG');
    assert.deepEqual(Object.keys(h.blocksForAllele(leftAligned, allele, start, seq)).sort(),
        ['exac', 'gnomad_exome', 'gnomad_genome', 'snpeff', 'vcf']);
    // MyVariant's "dup" id reads as a 2-bp duplication — another allele. Only
    // its ClinVar block is about this one, and ClinVar's own HGVS says so.
    const clinvarId = recordOf(cassettes.rs80357906, 'chr17:g.41209079_41209080dup');
    assert.equal(clinvarId.dbsnp.rsid, 'rs2051500205', 'the dbSNP block there is the 2-bp duplication\'s');
    assert.equal(h.sameGenomicAllele(clinvarId._id, allele, start, seq), false);
    assert.deepEqual(Object.keys(h.blocksForAllele(clinvarId, allele, start, seq)), ['clinvar']);
    // Asked about the 2-bp duplication itself, the same record gives its dbSNP
    // block and not ClinVar's.
    assert.deepEqual(Object.keys(h.blocksForAllele(clinvarId, 'chr17:g.41209079_41209080dup', start, seq)), ['dbsnp']);
    // A ClinVar block with no GRCh37 HGVS of its own is taken at its id's word.
    const bare = { _id: 'chr17:g.41209079_41209080insG', clinvar: { variant_id: 1 } };
    assert.deepEqual(Object.keys(h.blocksForAllele(bare, allele, start, seq)), ['clinvar']);
    assert.deepEqual(Object.keys(h.blocksForAllele({ ...bare, _id: 'chr17:g.41209079_41209080dup' }, allele, start, seq)), []);
    // Several ClinVar entries under one id: each is judged on its own HGVS.
    const two = { _id: 'chr17:g.41209079_41209080dup', clinvar: [
        { variant_id: 17677, hgvs: { genomic: ['NC_000017.10:g.41209082dup'] } },
        { variant_id: 999, hgvs: { genomic: ['NC_000017.10:g.41209079_41209080dup'] } }
    ] };
    assert.deepEqual(h.blocksForAllele(two, allele, start, seq).clinvar.map(c => c.variant_id), [17677]);
}
{
    // The hybrid id. Looked up as ClinVar's allele it gives ClinVar and nothing
    // else — its snpEff and vcf were computed for the id's literal allele.
    const { start, seq } = windowOf(cassettes.brca2ClinVar2505518);
    const hybrid = recordOf(cassettes.brca2ClinVar2505518, 'chr13:g.32906699_32906700insCA');
    assert.deepEqual(Object.keys(hybrid).filter(k => !k.startsWith('_')).sort(), ['chrom', 'clinvar', 'hg19', 'observed', 'snpeff', 'vcf']);
    assert.deepEqual(Object.keys(h.blocksForAllele(hybrid, 'chr13:g.32906700_32906701insCA', start, seq)), ['clinvar']);
    assert.deepEqual(Object.keys(h.blocksForAllele(hybrid, 'chr13:g.32906699_32906700insCA', start, seq)).sort(), ['snpeff', 'vcf']);
}

// ── Gathering ──────────────────────────────────────────────────────────────
{
    const win = windowOf(cassettes.rs80357906);
    const reference = async (chrom, start, end) => {
        assert.equal(chrom, 'chr17');
        return win.seq.slice(start - win.start, end - win.start + 1);
    };
    const asked = [];
    const records = async (ids) => {
        asked.push(ids);
        return ids.map(id => recordOf(cassettes.rs80357906, id) || recordOf(cassettes.brca1Dup, id)).filter(Boolean);
    };
    const primary = recordOf(cassettes.rs80357906, 'chr17:g.41209082dup');
    assert.deepEqual(Object.keys(primary).filter(k => !k.startsWith('_')).sort(), ['chrom', 'dbsnp', 'observed'], 'a dbSNP-only record');
    const out = await h.gatherAlleleRecords(primary, 'chr17:g.41209082dup', reference, records);
    assert.equal(out.status, 'merged');
    assert.deepEqual(asked, [['chr17:g.41209079_41209080insG', 'chr17:g.41209082_41209083insG', 'chr17:g.41209079_41209080dup']],
        'one request; the id already held is not asked for again');
    assert.equal(out.record._id, 'chr17:g.41209082dup');
    assert.equal(out.record.dbsnp.rsid, 'rs80357906', 'what the primary had is kept');
    assert.equal(out.record.clinvar.variant_id, 17677);
    assert.equal(typeof out.record.gnomad_exome.af.af, 'number');
    assert.deepEqual(out.record._also_filed_under, out.added);
    assert.deepEqual(out.added['chr17:g.41209079_41209080dup'], ['clinvar']);
    assert(out.added['chr17:g.41209079_41209080insG'].includes('gnomad_exome'));
    assert.equal(primary.clinvar, undefined, 'the primary itself is not modified');
    assert.match(h.describeSiblingBlocks(out.added),
        /^chr17:g\.41209079_41209080insG \(.*gnomAD.*snpEff.*\), chr17:g\.41209079_41209080dup \(ClinVar\)$/);
    assert.doesNotMatch(h.describeSiblingBlocks(out.added), /gnomAD, gnomAD/, 'exome and genome are named once');

    // A block the primary already has is never replaced.
    const withOwn = { ...primary, clinvar: { variant_id: 'mine' } };
    const kept = await h.gatherAlleleRecords(withOwn, 'chr17:g.41209082dup', reference, records);
    assert.equal(kept.record.clinvar.variant_id, 'mine');
    assert.equal(kept.added['chr17:g.41209079_41209080dup'], undefined);

    // A record re-keyed to ClinVar's allele was fetched under MyVariant's id
    // for it: the allele's own id is still to be asked, the borrowed one is not.
    asked.length = 0;
    const rekeyed = { _id: 'chr17:g.41209082dup', _myvariant_id: 'chr17:g.41209079_41209080dup', clinvar: { variant_id: 17677 } };
    const fromRekeyed = await h.gatherAlleleRecords(rekeyed, 'chr17:g.41209082dup', reference, records);
    assert.deepEqual(asked, [['chr17:g.41209079_41209080insG', 'chr17:g.41209082_41209083insG', 'chr17:g.41209082dup']]);
    assert.equal(fromRekeyed.record.dbsnp.rsid, 'rs80357906');
    assert.equal(fromRekeyed.record.clinvar.variant_id, 17677);

    // No record under the id asked for: what the others hold is still found.
    const fromNothing = await h.gatherAlleleRecords(null, 'chr17:g.41209082dup', reference, records);
    assert.equal(fromNothing.status, 'merged');
    assert.equal(fromNothing.record._id, 'chr17:g.41209082dup');
    assert.equal(fromNothing.record.clinvar.variant_id, 17677);

    // One spelling, one record: nothing is fetched.
    asked.length = 0;
    for (const single of ['chr17:g.41209082G>A', 'chr17:g.41209080_41209081delinsAT', 'NM_007294.4:c.5266dup', '']) {
        const none = await h.gatherAlleleRecords(primary, single, () => { throw new Error('no reference read'); }, records);
        assert.equal(none.status, 'not-needed', single);
        assert.equal(none.record, primary);
    }
    assert.deepEqual(asked, []);

    // "Could not check" is reported, never guessed around, and never throws.
    // (The page logs why to the console; not wanted in the test output.)
    const warn = console.warn;
    console.warn = () => {};
    const noReference = await h.gatherAlleleRecords(primary, 'chr17:g.41209082dup', async () => { throw new Error('503'); }, records);
    assert.deepEqual([noReference.status, noReference.detail, noReference.record], ['unverified', 'reference sequence unavailable', primary]);
    assert.deepEqual(asked, [], 'no ids are tried without the reference that derives and proves them');
    const noBatch = await h.gatherAlleleRecords(primary, 'chr17:g.41209082dup', reference,
        async () => { throw Object.assign(new Error('MyVariant.info batch request failed (503)'), { status: 503 }); });
    assert.deepEqual([noBatch.status, noBatch.detail, noBatch.record], ['unverified', 'HTTP 503', primary]);
    console.warn = warn;
    const nothingMore = await h.gatherAlleleRecords(primary, 'chr17:g.41209082dup', reference, async () => []);
    assert.deepEqual([nothingMore.status, nothingMore.record], ['none', primary]);
}

// ── A position belongs with the alleles it was written for ─────────────────
{
    // The record's vcf block is anchored at ITS position, not the HGVS start.
    const vcf = { position: '41209079', ref: 'T', alt: 'TG' };
    const hgvs = h.buildVariantCoordinateTuple('', 'chr17:g.41209082dup');
    assert.deepEqual(h.tupleVcfAlleles(hgvs, vcf), { pos: '41209079', ref: 'T', alt: 'TG' });
    assert.match(h.buildGnomadVariantUrl('', 'chr17:g.41209082dup', '', vcf), /\/variant\/17-41209079-T-TG\?/);
    assert.equal(h.buildSpliceAiApiVariant('', 'chr17:g.41209082dup', { vcf }), 'chr17-41209079-T-TG');
    // F508del typed as HGVS went to gnomAD as 7-117199645-ATCT-A.
    const f508 = { position: '117199644', ref: 'ATCT', alt: 'A' };
    assert.match(h.buildGnomadVariantUrl('', 'chr7:g.117199645_117199647del', '', f508), /\/variant\/7-117199644-ATCT-A\?/);
    // A pasted row carries its own alleles and keeps them.
    const row = h.buildVariantCoordinateTuple('17\t41209079\tT\tTG', 'chr17:g.41209082dup');
    assert.deepEqual(h.tupleVcfAlleles(row, { position: '1', ref: 'A', alt: 'C' }), { pos: '41209079', ref: 'T', alt: 'TG' });
    // The case this was first written for: an MNV, where the two positions coincide.
    assert.deepEqual(h.tupleVcfAlleles(h.buildVariantCoordinateTuple('', 'chr7:g.140453136_140453137delinsTT'),
        { position: '140453136', ref: 'AC', alt: 'TT' }), { pos: '140453136', ref: 'AC', alt: 'TT' });
    assert.equal(h.tupleVcfAlleles(hgvs, null), null);
    assert.equal(h.tupleVcfAlleles(hgvs, { position: '1', ref: '-', alt: 'G' }), null);
}

// ── ...and the alleles belong in the spelling the resources file them under ──
{
    // MyVariant writes a record's vcf block for the id it sits under. Under an
    // indel's 3'-most id that is not the left-aligned spelling gnomAD uses.
    const cftr = windowOf(cassettes.cftrF508del);
    const realign = (vcf, genomic, win = cftr) => h.leftAlignedVcfBlock(vcf, genomic, win.start, win.seq);
    const filed = recordOf(cassettes.cftrF508del, 'chr7:g.117199646_117199648del').vcf;
    assert.deepEqual([filed.position, filed.ref, filed.alt], ['117199645', 'TCTT', 'T'], 'as MyVariant has it');
    assert.deepEqual(realign(filed, 'chr7:g.117199646_117199648del'), { ...filed, position: '117199644', ref: 'ATCT', alt: 'A' });
    // Whichever spelling the allele itself is given in.
    assert.deepEqual(realign(filed, 'chr7:g.117199645_117199647del'), { ...filed, position: '117199644', ref: 'ATCT', alt: 'A' });
    // Already left-aligned, not an indel, or another allele's block: left alone.
    assert.equal(realign({ position: '117199644', ref: 'ATCT', alt: 'A' }, 'chr7:g.117199646_117199648del'), null);
    assert.equal(realign({ position: '117199645', ref: 'T', alt: 'C' }, 'chr7:g.117199646_117199648del'), null);
    assert.equal(realign({ position: '117199645', ref: 'TCT', alt: 'T' }, 'chr7:g.117199646_117199648del'), null, 'a 2-bp deletion is not F508del');
    assert.equal(realign({ position: '117199650', ref: 'TGG', alt: 'T' }, 'chr7:g.117199646_117199648del'), null);
    assert.equal(realign(null, 'chr7:g.117199646_117199648del'), null);
    assert.equal(realign(filed, 'chr7:g.117199646_117199648'), null, 'a bare span names no allele');
    // An insertion: 5382insC under its 3'-most id would read 17-41209081-G-GG.
    const brca1 = windowOf(cassettes.rs80357906);
    assert.deepEqual(realign({ position: '41209081', ref: 'G', alt: 'GG' }, 'chr17:g.41209082dup', brca1),
        { position: '41209079', ref: 'T', alt: 'TG' });
    assert.equal(realign({ position: '41209079', ref: 'T', alt: 'TG' }, 'chr17:g.41209082dup', brca1), null);

    // gatherAlleleRecords hands the record back with the block re-anchored —
    // also when the other ids had nothing to add — and says what it was.
    const reference = async (chrom, start, end) => cftr.seq.slice(start - cftr.start, end - cftr.start + 1);
    const primary = recordOf(cassettes.cftrF508del, 'chr7:g.117199646_117199648del');
    const none = await h.gatherAlleleRecords(primary, 'chr7:g.117199646_117199648del', reference, async () => []);
    assert.equal(none.status, 'none');
    assert.equal(none.vcfRespelled, true);
    assert.deepEqual([none.record.vcf.position, none.record.vcf.ref, none.record.vcf.alt], ['117199644', 'ATCT', 'A']);
    assert.deepEqual(none.record._vcf_as_filed, filed, 'what MyVariant has stays in the raw JSON');
    assert.deepEqual(primary.vcf, filed, 'the record passed in is not changed');
    const sibling = recordOf(cassettes.cftrF508del, 'chr7:g.117199645_117199647del');
    const merged = await h.gatherAlleleRecords(primary, 'chr7:g.117199646_117199648del', reference, async () => [sibling]);
    assert.equal(merged.status, 'merged');
    assert.equal(merged.record.vcf.position, '117199644');
    assert.equal(merged.record.clinvar.variant_id, 7105);
    // Nothing to re-anchor: no flag, and the very same record.
    const left = await h.gatherAlleleRecords(sibling, 'chr7:g.117199645_117199647del', reference, async () => []);
    assert.equal(left.vcfRespelled, undefined);
    assert.equal(left.record, sibling);
}
console.log('Sibling-record unit checks passed');

const batchIds = page => page.requests.filter(r => r.includes('[POST ids=')).map(r => r.match(/ids=([^\]]+)/)[1].split(','));
const rawJson = page => JSON.parse(page.document.getElementById('rawOutput').textContent);

// 1. rs80357906 lands on a dbSNP-only record. Was: "Clinical significance: N/A"
//    and no gnomAD card data, for the Ashkenazi founder mutation.
{
    const page = await runLookup(cassettes.rs80357906);
    assertResolved(assert, page, 'rs80357906');
    assert.match(page.variant, /g\.: chr17:g\.41209082dup/);
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.match(page.clinvar, /Clinical significance: Pathogenic \(14\)/);
    assert.match(page.cards.gnomAD, /Exome AF: 0\.0183% AC\/AN: 46\/251494/);
    assert.match(page.variant, /VCF \(hg19\): 17-41209079-T-TG/);
    assert.match(page.panel, /also filed under chr17:g\.41209079_41209080insG \([^)]*gnomAD[^)]*\), chr17:g\.41209079_41209080dup \(ClinVar\)/);
    assert.deepEqual(batchIds(page), [['chr17:g.41209079_41209080insG', 'chr17:g.41209082_41209083insG', 'chr17:g.41209079_41209080dup']]);
    const raw = rawJson(page);
    assert.equal(raw._id, 'chr17:g.41209082dup');
    assert.equal(raw.dbsnp.rsid, 'rs80357906', 'not the 2-bp duplication\'s rs2051500205');
    assert.deepEqual(Object.keys(raw._also_filed_under), ['chr17:g.41209079_41209080insG', 'chr17:g.41209079_41209080dup']);
    console.log('ok sibling records: rs80357906 gets ClinVar 17677 and gnomAD');
}

// 2. The same variant by name lands on the ClinVar record instead. Was: no
//    population frequency, no rsID.
{
    const page = await runLookup(cassettes.brca1Dup);
    assertResolved(assert, page, 'BRCA1 c.5266dup');
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.match(page.cards.gnomAD, /Exome AF: 0\.0183% AC\/AN: 46\/251494/);
    assert.equal(rawJson(page).dbsnp.rsid, 'rs80357906');
    assert.deepEqual(batchIds(page), [['chr17:g.41209079_41209080insG', 'chr17:g.41209082_41209083insG', 'chr17:g.41209082dup']]);
    // The SpliceAI link pairs the record's alleles with the record's position.
    assert.match(page.cards.SpliceAI, /SpliceAI query: chr17 41209079 T TG/);
    assert.doesNotMatch(page.cards.SpliceAI, /41209082 T TG/);
    console.log('ok sibling records: BRCA1 c.5266dup gets gnomAD and rs80357906');
}

// 3. A pasted row lands on the gnomAD record. Was: "Clinical significance: N/A".
{
    const page = await runLookup(cassettes.tp53Dup);
    assertResolved(assert, page, 'TP53 c.455dup');
    assert.match(page.variant, /g\.: chr17:g\.7578474_7578475insG/);
    assert.match(page.clinvar, /Variation ID: 182958/);
    assert.match(page.clinvar, /Clinical significance: Pathogenic \(4\)/);
    assert.equal(rawJson(page).dbsnp.rsid, 'rs730882019');
    assert.match(page.panel, /also filed under chr17:g\.7578479dup \(dbSNP\), chr17:g\.7578474_7578475dup \(ClinVar\)/);
    console.log('ok sibling records: 17 7578474 TP53 C CG gets ClinVar 182958');
}

// 4. F508del's rsID is under the 3'-most spelling.
{
    const page = await runLookup(cassettes.cftrF508del);
    assertResolved(assert, page, 'F508del');
    assert.match(page.clinvar, /Variation ID: 7105/);
    assert.equal(rawJson(page).dbsnp.rsid, 'rs113993960');
    assert.match(page.panel, /also filed under chr7:g\.117199646_117199648del \(dbSNP\)/);
    console.log('ok sibling records: F508del gets rs113993960');
}

// 5. MyVariant holds nothing under either real spelling of BRCA2
//    c.1085_1086insCA — only ClinVar's record, under an id of its own making.
{
    const page = await runLookup(cassettes.brca2ClinVar2505518);
    assertResolved(assert, page, 'BRCA2 c.1085_1086insCA');
    assert.match(page.variant, /c\.: c\.1085_1086insCA/);
    assert.match(page.clinvar, /Variation ID: 2505518/);
    assert.match(page.panel, /no record for this variant · also filed under chr13:g\.32906699_32906700insCA \(ClinVar\)/);
    // Only ClinVar's block: that id's snpEff and vcf describe another allele.
    const raw = rawJson(page);
    assert.equal(raw.snpeff, undefined);
    assert.equal(raw.vcf, undefined);
    console.log('ok sibling records: an allele with no record of its own still finds its ClinVar entry');
}

// 6. ...and the allele that id literally spells does not get that entry.
{
    const page = await runLookup(cassettes.brca2Neighbour);
    assertResolved(assert, page, 'BRCA2 c.1084_1085insCA');
    assert.match(page.variant, /g\.: chr13:g\.32906699_32906700insCA/);
    assert.match(page.variant, /c\.: c\.1084_1085insCA/);
    assert.doesNotMatch(page.clinvar, /2505518/);
    assert.match(page.clinvar, /Clinical significance: N\/A/);
    assert.equal(rawJson(page).clinvar, undefined);
    console.log('ok sibling records: a neighbouring allele does not inherit it');
}

// 7. Without the reference nothing can be derived or proven: the lookup still
//    resolves, with what the one record holds, and says what was not checked.
{
    const page = await runLookup(cassettes.rs80357906, { without: [/\/sequence\/region\//] });
    assertResolved(assert, page, 'rs80357906, no reference');
    assert.match(page.variant, /g\.: chr17:g\.41209082dup/);
    assert.match(page.clinvar, /Clinical significance: N\/A/);
    assert.deepEqual(batchIds(page), []);
    assert.match(page.panel, /its other ids for this allele were not checked — reference sequence unavailable/);
    console.log('ok sibling records: no reference, no merge — and it says so');
}

// 8. MyVariant's batch endpoint failing is reported the same way.
{
    const page = await runLookup(cassettes.rs80357906, {
        route: (url, key, respond, options) => (String(options.method).toUpperCase() === 'POST' && url.hostname === 'myvariant.info'
            ? respond({ error: 'Service Unavailable' }, 503) : undefined)
    });
    assertResolved(assert, page, 'rs80357906, batch down');
    assert.match(page.clinvar, /Clinical significance: N\/A/);
    assert.match(page.panel, /its other ids for this allele were not checked — HTTP 503/);
    console.log('ok sibling records: batch endpoint down is reported, not fatal');
}

// 9. A substitution has one id: no reference read, no batch request.
{
    const page = await runLookup(cassettes.brafV600E);
    assertResolved(assert, page, 'BRAF V600E');
    assert.deepEqual(batchIds(page), []);
    assert.deepEqual(page.requests.filter(r => r.includes('/sequence/region/')), []);
    console.log('ok sibling records: a substitution costs nothing extra');
}

// 10. Ensembl's sequence service down (it was, for part of the live run these
//     scenarios come from). Whether two spellings are one allele is normally
//     read off the reference; without it the answer used to be "cannot tell"
//     for every pair, and scenario 6's neighbour was shown ClinVar 2505518 —
//     "Likely pathogenic" — after all. Some pairs need no reference: insCA one
//     position along would have to read insAC to be the same allele.
{
    const page = await runLookup(cassettes.brca2Neighbour, { without: [/\/sequence\/region\//] });
    assertResolved(assert, page, 'BRCA2 c.1084_1085insCA, no reference');
    assert.match(page.variant, /g\.: chr13:g\.32906699_32906700insCA/);
    assert.doesNotMatch(page.clinvar, /2505518|Likely pathogenic/);
    assert.match(page.clinvar, /Clinical significance: N\/A/);
    assert.equal(rawJson(page).clinvar, undefined);
    console.log('ok sibling records: no reference, and the neighbouring allele still does not inherit it');
}

// 11. The same outage, by name. MyVariant's id for ClinVar 17677 is a 2-bp
//     duplication and ClinVar's own HGVS a 1-bp one: different alleles, whatever
//     the reference holds. Was: "g.: chr17:g.41209079_41209080dup" with a
//     warning that the id could not be checked.
{
    const page = await runLookup(cassettes.brca1Dup, { without: [/\/sequence\/region\//] });
    assertResolved(assert, page, 'BRCA1 c.5266dup, no reference');
    assert.match(page.variant, /g\.: chr17:g\.41209082dup/);
    assert.doesNotMatch(page.variant, /41209079_41209080dup/);
    assert.match(page.clinvar, /Variation ID: 17677/);
    assert.match(page.panel, /ClinVar record found — allele chr17:g\.41209082dup from ClinVar's own HGVS \(MyVariant files it under chr17:g\.41209079_41209080dup\)/);
    assert.doesNotMatch(page.panel, /could not be checked against/);
    // What could not be done without the reference is still said.
    assert.match(page.panel, /its other ids for this allele were not checked — reference sequence unavailable/);
    console.log('ok sibling records: no reference, and ClinVar\'s own allele is still the one shown');
}

// 12. VEP's RefSeq transcript set failing while its Ensembl one answers (the
//     same live run; the rows are Ensembl's, as recorded). VEP is asked here
//     without being needed, so its answer must not displace the RefSeq row the
//     record has: the Variant card came out on ENST00000544455.1.
{
    const row = (transcript_id, hgvsc, hgvsp, extra = {}) => ({ gene_symbol: 'BRCA2', gene_id: 'ENSG00000139618', transcript_id,
        consequence_terms: ['frameshift_variant'], hgvsc, hgvsp, biotype: 'protein_coding', strand: 1, impact: 'HIGH', variant_allele: 'CA', ...extra });
    const ensemblSet = [{ input: '13:g.32906699_32906700insCA', id: '13:g.32906699_32906700insCA', assembly_name: 'GRCh37',
        seq_region_name: '13', start: 32906700, end: 32906699, strand: 1, allele_string: '-/CA', most_severe_consequence: 'frameshift_variant',
        transcript_consequences: [
            row('ENST00000380152', 'ENST00000380152.3:c.1084_1085insCA', 'ENSP00000369497.3:p.Asp362AlafsTer6'),
            row('ENST00000530893', 'ENST00000530893.2:c.715_716insCA', 'ENSP00000435699.2:p.Asp239AlafsTer6'),
            row('ENST00000544455', 'ENST00000544455.1:c.1084_1085insCA', 'ENSP00000439902.1:p.Asp362AlafsTer6', { canonical: 1 })
        ] }];
    const vep = [];
    const page = await runLookup(cassettes.brca2Neighbour, {
        route: (url, key, respond) => {
            if (!/\/vep\/human\/hgvs\//.test(key)) return undefined;
            vep.push(/refseq=1/.test(key) ? 'RefSeq' : 'Ensembl');
            return /refseq=1/.test(key) ? respond({ error: 'Internal Server Error' }, 500) : respond(ensemblSet);
        }
    });
    assertResolved(assert, page, 'BRCA2 c.1084_1085insCA, RefSeq VEP down');
    assert.equal(vep.at(-1), 'Ensembl', 'VEP fell back to its own transcript set');
    assert.match(page.variant, /Transcript: NM_000059\.3c\.: c\.1084_1085insCA/);
    assert.doesNotMatch(page.variant, /Transcript: ENST/);
    assert.match(page.panel, /Consequences · Ensembl VEPRefSeq transcripts unavailable — the record's own are shown/);
    // RefSeq answering that it has no transcript there is not that outage: then
    // Ensembl's rows are what there is to show.
    const noRefSeq = await runLookup(cassettes.brca2Neighbour, {
        route: (url, key, respond) => {
            if (!/\/vep\/human\/hgvs\//.test(key)) return undefined;
            return /refseq=1/.test(key) ? respond([{ ...ensemblSet[0], transcript_consequences: [] }]) : respond(ensemblSet);
        }
    });
    assertResolved(assert, noRefSeq, 'BRCA2 c.1084_1085insCA, no RefSeq transcript');
    assert.match(noRefSeq.variant, /Transcript: ENST00000544455\.1c\.: c\.1084_1085insCA/);
    assert.match(noRefSeq.panel, /●Consequences · Ensembl VEP3 transcripts/);
    // Where the record has no transcript rows of its own, Ensembl's are the
    // answer: there is no RefSeq row for them to displace. (rs80357906 lands on
    // a dbSNP-only record; with the reference withheld nothing is added to it.
    // Ensembl's set answers as the recorded RefSeq one does, on its own
    // transcript.)
    const refseqKey = 'grch37.rest.ensembl.org/vep/human/hgvs/17:g.41209082dup?content-type=application/json&hgvs=1&refseq=1&mane=1&canonical=1';
    const recorded = cassettes.rs80357906.responses[refseqKey].body[0];
    const ensemblOnly = [{ ...recorded, transcript_consequences: [{ ...recorded.transcript_consequences[0], gene_id: 'ENSG00000012048',
        transcript_id: 'ENST00000357654', hgvsc: 'ENST00000357654.3:c.5266dup', hgvsp: 'ENSP00000350283.3:p.Gln1756ProfsTer74', canonical: 1 }] }];
    const sets = [];
    const dbsnpOnly = await runLookup(cassettes.rs80357906, {
        route: (url, key, respond) => {
            if (!/\/vep\/human\/hgvs\//.test(key)) return undefined;
            sets.push(/refseq=1/.test(key) ? 'RefSeq' : 'Ensembl');
            return /refseq=1/.test(key) ? respond({ error: 'Internal Server Error' }, 500) : respond(ensemblOnly);
        },
        without: [/\/sequence\/region\//]
    });
    assertResolved(assert, dbsnpOnly, 'rs80357906, RefSeq VEP down');
    assert.equal(sets.at(-1), 'Ensembl', 'VEP fell back to its own transcript set');
    assert.match(dbsnpOnly.variant, /Transcript: ENST00000357654(?:\.\d+)?c\.: c\.5266dupp\.: p\.Gln1756ProfsTer74/);
    assert.match(dbsnpOnly.panel, /●Consequences · Ensembl VEP1 transcript/);
    assert.doesNotMatch(dbsnpOnly.panel, /the record's own are shown/);
    console.log('ok sibling records: an Ensembl-only VEP answer does not displace the record\'s RefSeq row');
}

// 13. The lookup lands on the record under the 3'-most id (through the HGVS g.
//     name, or rs113993960) and its vcf block is anchored there. gnomAD v4.1
//     holds F508del as 7-117199644-ATCT-A and answered "Not found" when asked
//     for 7-117199645-TCTT-T; the v2 link opened a "variant not found" page.
{
    const page = await runLookup(cassettes.cftrF508del, { input: 'chr7:g.117199646_117199648del' });
    assertResolved(assert, page, 'F508del by its HGVS g. name');
    assert.match(page.variant, /g\.: chr7:g\.117199646_117199648delVCF \(hg19\): 7-117199644-ATCT-A/);
    assert.deepEqual(page.api('gnomad-v4'), ['gnomad-v4?chrom=7&pos37=117199644&ref=ATCT&alt=A']);
    assert.deepEqual(page.api('spliceai').map(r => r.match(/variant=([^&]+)/)[1]), ['chr7-117199644-ATCT-A']);
    assert(page.links.some(l => /gnomad\.broadinstitute\.org\/variant\/7-117199644-ATCT-A\?dataset=gnomad_r2_1$/.test(l.href)));
    assert(!page.links.some(l => /117199645-TCTT-T/.test(l.href)), 'no link under the other spelling');
    assert.deepEqual(rawJson(page)._vcf_as_filed, { alt: 'T', position: '117199645', ref: 'TCTT' });
    // 185delAG, by its legacy name: 17-41276044-ACT-A, not 17-41276046-TCT-T.
    const legacy = await runLookup(cassettes.brca1Legacy185delAG);
    assertResolved(assert, legacy, 'BRCA1 185delAG');
    assert.match(legacy.variant, /VCF \(hg19\): 17-41276044-ACT-A/);
    assert.deepEqual(legacy.api('gnomad-v4'), ['gnomad-v4?chrom=17&pos37=41276044&ref=ACT&alt=A']);
    // A row pasted with its own alleles keeps them.
    const row = await runLookup(cassettes.cftrF508del);
    assert.deepEqual(row.api('gnomad-v4'), ['gnomad-v4?chrom=7&pos37=117199644&ref=ATCT&alt=A']);
    assert.equal(rawJson(row)._vcf_as_filed, undefined);
    console.log('ok sibling records: the VCF alleles are left-aligned whichever record the lookup landed on');
}

console.log('Sibling-record tests passed');
