/*
 * When Ensembl's GRCh37 service is down.
 *
 * Every lookup needs Ensembl's VEP for what a variant does on each transcript,
 * and a variant nobody else holds a record of needs it to be shown at all. It
 * was asked of one host, grch37.rest.ensembl.org, and on the day that host
 * answered HTTP 500 from every endpoint — its "ping" included; a browser calls
 * that a network error, because the error page carries no CORS header —
 *
 *   17 7576910 TP53 G GGT
 *     ✕ Resolve nomenclature · Ensembl   network error
 *     ◍ Annotations · MyVariant.info     no record for this variant
 *     ✕ Consequences · Ensembl VEP       network error
 *     Error: Lookup could not complete: … did not respond. This is an
 *     upstream outage, not necessarily a problem with your variant
 *
 * Ensembl keeps an archive of the same services for the same assembly
 * (grch37-archive.rest.ensembl.org, release 99), on other machines; it went on
 * answering. For a genomic ALLELE — a position means the same in both — VEP is
 * now asked there when the current service fails, or has not answered after a
 * head start, and the result says so: the archive's transcripts are in the
 * versions of 2020, and none of them is flagged MANE Select.
 *
 * The scenarios replay answers recorded from the live services during the
 * outage: the current service's are its HTTP 500s.
 *
 * Run with: node tests/ensembl-archive.test.js
 */
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runLookup, assertResolved } from './support/lookup-harness.js';

await import('../script.js');
const h = globalThis.__variantSearchHelpers;
const cassettes = JSON.parse(readFileSync(new URL('./fixtures/lookup-cassettes.json', import.meta.url)));
const wait = ms => new Promise(resolve => setTimeout(resolve, ms));

// ── One answer from two sources ────────────────────────────────────────────
{
    const after = (ms, value) => new Promise(resolve => setTimeout(() => resolve(value), ms));
    const failing = (ms, message) => new Promise((_, reject) => setTimeout(() => reject(new Error(message)), ms));
    const either = h.answerFromEither;
    let asked = 0;
    const second = (ms, value) => () => { asked += 1; return after(ms, value); };
    let began = Date.now();
    const took = () => Date.now() - began;

    // The first answers within its head start: the second never hears of it.
    assert.deepEqual(await either(after(20, 'first'), second(20, 'second'), 400), { answer: 'first', second: false });
    await wait(450);
    assert.equal(asked, 0);
    // The first has no answer to give: the second is asked at once.
    began = Date.now();
    assert.deepEqual(await either(after(20, null), second(40, 'second'), 5000), { answer: 'second', second: true });
    assert(took() < 2000, 'without waiting out the head start');
    assert.equal(asked, 1);
    // ...and has none either.
    assert.deepEqual(await either(after(20, null), second(40, null), 5000), { answer: null, second: false });
    assert.deepEqual(await either(after(20, undefined), second(40, undefined), 5000), { answer: null, second: false });
    // The first is slow: the second is asked when the head start is over, and
    // whichever of the two answers first is taken.
    began = Date.now();
    asked = 0;
    assert.deepEqual(await either(after(2500, 'first'), second(20, 'second'), 250), { answer: 'second', second: true });
    assert(took() >= 240 && took() < 1500, 'the second was asked after the head start, and the first was not waited for');
    began = Date.now();
    assert.deepEqual(await either(after(600, 'first'), second(2500, 'second'), 250), { answer: 'first', second: false });
    assert(took() < 1800);
    assert.equal(asked, 2, 'asked both times');
    // One of the two has nothing: the other is waited for.
    assert.deepEqual(await either(after(500, null), second(600, 'second'), 250), { answer: 'second', second: true });
    assert.deepEqual(await either(after(900, 'first'), second(50, null), 250), { answer: 'first', second: false });
    assert.deepEqual(await either(after(500, null), second(100, null), 250), { answer: null, second: false });
    // (An answer is anything but nothing.)
    assert.deepEqual(await either(after(10, 0), second(10, 'second'), 400), { answer: 0, second: false });
    assert.deepEqual(await either(after(10, null), second(10, ''), 400), { answer: '', second: true });
    // The first refuses — rejects: that is its answer, and the second is not asked to overrule it.
    asked = 0;
    await assert.rejects(either(failing(20, 'refused'), second(10, 'second'), 400), /^Error: refused$/);
    assert.equal(asked, 0);
    await assert.rejects(either(failing(500, 'refused late'), second(2500, 'second'), 200), /^Error: refused late$/);
    assert.equal(asked, 1, 'it had been asked by then, and is not waited for');
    console.log('ensembl-archive: one answer from two sources — unit checks passed');
}

// ── VEP from the current service, or from the archive ──────────────────────
const CURRENT = 'https://grch37.rest.ensembl.org/vep/human/hgvs/';
const ARCHIVE = 'https://grch37-archive.rest.ensembl.org/vep/human/hgvs/';
const refseqSet = '?content-type=application/json&hgvs=1&refseq=1&mane=1&canonical=1';
const ensemblSet = '?content-type=application/json&hgvs=1&mane=1&canonical=1';
{
    const realFetch = globalThis.fetch;
    const realWarn = console.warn;
    const real = { head: h.API_TIMEOUT_MS.vepArchiveHeadStart };
    console.warn = () => {};
    const allele = 'chr17:g.7576910_7576911insGT';
    const asAsked = encodeURIComponent('17:g.7576910_7576911insGT');
    const vepAnswer = (transcript, protein) => [{ input: '17:g.7576910_7576911insGT', assembly_name: 'GRCh37', seq_region_name: '17',
        start: 7576911, end: 7576910, strand: 1, allele_string: '-/GT', most_severe_consequence: 'frameshift_variant',
        transcript_consequences: [{ transcript_id: transcript, hgvsc: `${transcript}:c.934_935dup`, hgvsp: `${protein}:p.Ser313ProfsTer33`,
            gene_symbol: 'TP53', consequence_terms: ['frameshift_variant'], strand: -1, canonical: 1 }] }];
    const current = vepAnswer('NM_000546.6', 'NP_000537.3');
    const archived = vepAnswer('NM_000546.5', 'NP_000537.3');
    const json = (body, status = 200) => () => new Response(JSON.stringify(body), { status });
    // (As Ensembl was that day: an HTML page, HTTP 500, for every request.)
    const down = () => new Response('<!doctype html><html lang="en"><head><title>Internal Server Error</title></head></html>', { status: 500 });
    const unreachable = () => { throw new TypeError('Failed to fetch'); };
    const refusal = json({ error: 'Unable to parse HGVS notation \'17:g.7576910_7576911insGT\': Reference allele extracted from 17 does not match' }, 400);
    const after = (ms, answer) => (url, options) => new Promise((resolve, reject) => {
        const timer = setTimeout(() => resolve(answer()), ms);
        options.signal?.addEventListener('abort', () => { clearTimeout(timer); reject(Object.assign(new Error('This operation was aborted'), { name: 'AbortError' })); });
    });
    const serve = (hosts) => {
        const asked = [];
        globalThis.fetch = async (address, options = {}) => {
            const url = String(address);
            asked.push(url);
            return (url.startsWith(ARCHIVE) ? hosts.archive : hosts.current)(url, options);
        };
        return asked;
    };
    const ofArchive = asked => asked.filter(url => url.startsWith(ARCHIVE));
    try {
        // The current service answers: the archive is not asked.
        let asked = serve({ current: json(current), archive: json(archived) });
        let vep = await h.fetchVepHgvsHg19(allele);
        assert.equal(vep.archive, false);
        assert.equal(vep.consequences[0].transcript_id, 'NM_000546.6');
        assert.deepEqual([vep.mostSevere, vep.refseqFailed], ['frameshift_variant', false]);
        assert.deepEqual(asked, [`${CURRENT}${asAsked}${refseqSet}`]);

        // The current service is down — HTTP 500 for each of its attempts: the
        // archive is asked the same question, and what comes back is marked as its.
        asked = serve({ current: down, archive: json(archived) });
        vep = await h.fetchVepHgvsHg19(allele);
        assert.equal(vep.archive, true);
        assert.equal(vep.consequences[0].transcript_id, 'NM_000546.5');
        assert.equal(vep.mostSevere, 'frameshift_variant');
        assert.deepEqual(ofArchive(asked), [`${ARCHIVE}${asAsked}${refseqSet}`]);
        assert.equal(asked.length - 1, 4, 'the current service was tried: twice for each transcript set');
        // ...or cannot be reached at all, which is what a browser makes of that page.
        asked = serve({ current: unreachable, archive: json(archived) });
        vep = await h.fetchVepHgvsHg19(allele);
        assert.equal(vep.archive, true);
        assert.equal(ofArchive(asked).length, 1);

        // The rows of such an answer say where they came from; other rows are as they were.
        const rows = h.transcriptsFromVep(vep.consequences, vep);
        assert.deepEqual(rows.map(row => [row.transcript, row.cDNA, row.protein, row.source, row.archive]),
            [['NM_000546.5', 'c.934_935dup', 'p.Ser313ProfsTer33', 'vep', true]]);
        for (const usual of [h.transcriptsFromVep(vep.consequences), h.transcriptsFromVep(vep.consequences, { archive: false }), h.transcriptsFromVep(vep.consequences, null)]) {
            assert.equal('archive' in usual[0], false);
            assert.equal(usual[0].source, 'vep');
        }

        // The current service REFUSES the allele: that is its answer. The archive is not asked.
        asked = serve({ current: refusal, archive: json(archived) });
        await assert.rejects(h.fetchVepHgvsHg19(allele), (err) => err.status === 400 && !err.retryable && /Reference allele extracted/.test(err.message));
        assert.deepEqual(ofArchive(asked), []);
        // Both are down: the failure told is the current service's.
        asked = serve({ current: down, archive: down });
        await assert.rejects(h.fetchVepHgvsHg19(allele), (err) => err.status === 500 && err.retryable === true);
        assert(ofArchive(asked).length >= 2, 'the archive was tried');
        asked = serve({ current: unreachable, archive: unreachable });
        await assert.rejects(h.fetchVepHgvsHg19(allele), (err) => /Failed to fetch/.test(err.message) && err.retryable === true);
        // The archive's refusal is not taken for the current service's: it is still an outage.
        asked = serve({ current: down, archive: refusal });
        await assert.rejects(h.fetchVepHgvsHg19(allele), (err) => err.status === 500 && err.retryable === true);
        // (An answer with nothing in it is no answer either.)
        asked = serve({ current: down, archive: json([]) });
        await assert.rejects(h.fetchVepHgvsHg19(allele), (err) => err.status === 500);

        // A NAME is read on the transcripts a service has, and the archive has
        // none of the newer RefSeq versions: it is not asked about names.
        for (const name of ['NM_000546.6:c.934_935dup', 'TP53:c.934_935dup', 'NP_000537.3:p.Arg175His', 'ENSP00000269305.4:p.Arg175His', 'rs28934578', 'LRG_321t1:c.524G>A']) {
            asked = serve({ current: down, archive: json(archived) });
            await assert.rejects(h.fetchVepHgvsHg19(name), (err) => err.status === 500, name);
            assert.deepEqual(ofArchive(asked), [], name);
        }
        // ...except a name on one of Ensembl's own transcripts: GRCh37's are the
        // same in the archive, version for version.
        for (const name of ['ENST00000269305.4:c.934_935dup', 'ENST00000269305:c.934_935dup', 'enst00000269305.4:c.934_935dup', 'ENST00000269305.4:n.1124_1125dup']) {
            asked = serve({ current: down, archive: json(archived) });
            assert.equal((await h.fetchVepHgvsHg19(name)).archive, true, name);
            assert.deepEqual(ofArchive(asked), [`${ARCHIVE}${encodeURIComponent(name)}${refseqSet}`], name);
        }
        // Every chromosome's alleles are, the mitochondrial genome's on Ensembl's own transcripts.
        for (const [genomic, asUrl, set] of [['chrX:g.20156693_20156694insAGA', 'X:g.20156693_20156694insAGA', refseqSet],
            ['chr22:g.100A>T', '22:g.100A>T', refseqSet], ['chrMT:g.8993T>G', 'MT:g.8993T>G', ensemblSet]]) {
            asked = serve({ current: down, archive: json(archived) });
            assert.equal((await h.fetchVepHgvsHg19(genomic)).archive, true, genomic);
            assert.deepEqual(ofArchive(asked), [`${ARCHIVE}${encodeURIComponent(asUrl)}${set}`], genomic);
        }
        // Rows the archive cannot be relied on for are left out: VEP's own flag
        // that a RefSeq transcript does not match the genome before or within its
        // coding sequence. (PTEN, as the archive answered: NM_000314.6 a base off.)
        const pten = (transcript, c, flags) => ({ transcript_id: transcript, hgvsc: `${transcript}:${c}`, hgvsp: 'NP_000305.3:p.Pro367GlnfsTer49',
            gene_symbol: 'PTEN', consequence_terms: ['frameshift_variant'], strand: 1, ...(flags ? { refseq_match: flags } : {}) });
        const ptenRows = [pten('NM_000314.4', 'c.1100del'), pten('NM_000314.6', 'c.1101del', ['rseq_mrna_nonmatch', 'rseq_5p_mismatch']),
            pten('NM_001304717.2', 'c.1619del', ['rseq_mrna_nonmatch', 'rseq_cds_mismatch']), pten('NM_000314.5', 'c.1100del', ['rseq_mrna_nonmatch', 'rseq_3p_mismatch'])];
        assert.deepEqual(ptenRows.map(h.archiveRowReliable), [true, false, false, true], 'a mismatch after the coding sequence moves no coding position');
        assert.deepEqual([{ refseq_match: ['rseq_mrna_match'] }, { refseq_match: ['rseq_ens_match_cds', 'rseq_mrna_nonmatch'] }, { refseq_match: 'rseq_5p_mismatch' },
            { refseq_match: [] }, {}, null].map(h.archiveRowReliable), [true, true, false, true, true, true]);
        const ptenAnswer = [{ ...archived[0], most_severe_consequence: 'frameshift_variant', transcript_consequences: ptenRows }];
        asked = serve({ current: down, archive: json(ptenAnswer) });
        vep = await h.fetchVepHgvsHg19('chr10:g.89725116del');
        assert.deepEqual(vep.consequences.map(row => row.hgvsc), ['NM_000314.4:c.1100del', 'NM_000314.5:c.1100del']);
        assert.equal(ofArchive(asked).length, 1, 'what is left names the change: nobody else is asked');
        // ...of the archive's answer only: the current service's rows are all its own.
        asked = serve({ current: json(ptenAnswer), archive: json(archived) });
        assert.equal((await h.fetchVepHgvsHg19('chr10:g.89725116del')).consequences.length, 4);
        // None left that names a coding change: Ensembl's own transcripts are asked for.
        const upstream = { transcript_id: 'NM_000999.1', gene_symbol: 'OTHER', consequence_terms: ['upstream_gene_variant'], strand: 1 };
        for (const flaggedOnly of [[ptenRows[1], ptenRows[2]], [ptenRows[1], upstream]]) {
            asked = serve({ current: down, archive: url => (url.includes('refseq=1') ? json([{ ...archived[0], transcript_consequences: flaggedOnly }])()
                : json(vepAnswer('ENST00000371953.3', 'ENSP00000361021.3'))()) });
            vep = await h.fetchVepHgvsHg19('chr10:g.89725116del');
            assert.deepEqual([vep.archive, vep.refseqFailed, vep.consequences.map(row => row.transcript_id)], [true, false, ['ENST00000371953.3']]);
            assert.deepEqual(ofArchive(asked).map(url => url.includes('refseq=1')), [true, false]);
        }
        // (Rows without a coding name left, and none taken out: that is RefSeq's answer, as it always was.)
        asked = serve({ current: down, archive: json([{ ...archived[0], transcript_consequences: [upstream] }]) });
        vep = await h.fetchVepHgvsHg19('chr10:g.89725116del');
        assert.deepEqual(vep.consequences.map(row => row.transcript_id), ['NM_000999.1']);
        assert.equal(ofArchive(asked).length, 1);

        // Where the archive's RefSeq set has no transcript here, its own set is asked, as the current service's is.
        asked = serve({ current: down, archive: url => (url.includes('refseq=1') ? json([{ ...archived[0], transcript_consequences: [] }])() : json(vepAnswer('ENST00000269305.4', 'ENSP00000269305.4'))()) });
        vep = await h.fetchVepHgvsHg19(allele);
        assert.deepEqual([vep.archive, vep.refseqFailed, vep.consequences[0].transcript_id], [true, false, 'ENST00000269305.4']);
        assert.deepEqual(ofArchive(asked).map(url => url.includes('refseq=1')), [true, false]);

        // The current service is given a head start; the archive is asked when it is over.
        h.API_TIMEOUT_MS.vepArchiveHeadStart = 250;
        // In good time: the archive never hears of it.
        asked = serve({ current: after(30, json(current)), archive: json(archived) });
        assert.equal((await h.fetchVepHgvsHg19(allele)).archive, false);
        await wait(400);
        assert.deepEqual(ofArchive(asked), []);
        // Silent: the archive's answer is taken, and the current service is not waited for.
        asked = serve({ current: after(3000, json(current)), archive: after(30, json(archived)) });
        let began = Date.now();
        vep = await h.fetchVepHgvsHg19(allele);
        assert.equal(vep.archive, true);
        assert(Date.now() - began < 1800, 'not waited for');
        // Slow, and still the first of the two: the current service's answer is the one shown.
        asked = serve({ current: after(600, json(current)), archive: after(3000, json(archived)) });
        began = Date.now();
        vep = await h.fetchVepHgvsHg19(allele);
        assert.deepEqual([vep.archive, vep.consequences[0].transcript_id], [false, 'NM_000546.6']);
        assert.equal(ofArchive(asked).length, 1, 'the archive had been asked');
        assert(Date.now() - began < 2200);
        // Slow, and then refusing: the refusal stands, whatever the archive would say.
        asked = serve({ current: after(500, refusal), archive: after(3000, json(archived)) });
        await assert.rejects(h.fetchVepHgvsHg19(allele), (err) => err.status === 400);
    } finally {
        globalThis.fetch = realFetch;
        console.warn = realWarn;
        h.API_TIMEOUT_MS.vepArchiveHeadStart = real.head;
    }
    console.log('ensembl-archive: VEP from the archive — unit checks passed');
}

// ── Who may read a name ────────────────────────────────────────────────────
{
    const hosts = query => h.recoderHostsFor(query).map(host => host.replace('https://', ''));
    const mirror = ['grch37.rest.ensembl.org'];
    const both = ['grch37.rest.ensembl.org', 'rest.ensembl.org'];
    // A name on one of Ensembl's own transcripts or proteins is read on GRCh37's
    // version of it. The main host has another version, and reads the name on
    // that whatever version is asked for.
    for (const name of ['ENST00000288602.6:c.1799T>A', 'ENST00000288602:c.1799T>A', 'enst00000288602.6:c.1799T>A',
        'ENST00000361899.2:n.10A>T', 'ENSP00000288602.6:p.Val600Glu', ' ENST00000288602.6:c.1799T>A']) {
        assert.deepEqual(hosts(name), mirror, name);
    }
    // RefSeq's transcripts and proteins are the same sequence under the same
    // accession and version on both; a gene symbol and an rsID name no version.
    for (const name of ['NM_004333.6:c.1799T>A', 'NP_004324.2:p.Val600Glu', 'BRAF:c.1799T>A', 'BRAF:p.Val600Glu', 'rs113488022', 'LRG_299t1:c.1799T>A']) {
        assert.deepEqual(hosts(name), both, name);
    }
    // (As before: a chromosomal position is GRCh37's alone; the mitochondrion's is the same in both.)
    assert.deepEqual(hosts('chr7:g.140453136A>T'), mirror);
    assert.deepEqual(hosts('chrMT:g.8993T>G'), both);
    assert.deepEqual([h.namesEnsemblTranscript('ENST00000288602.6:c.1799T>A'), h.namesEnsemblTranscript('ENST00000288602:n.5A>T'),
        h.namesEnsemblTranscript('ENSP00000288602.6:p.Val600Glu'), h.namesEnsemblTranscript('NM_004333.6:c.1799T>A'),
        h.namesEnsemblTranscript('ENST00000288602.6'), h.namesEnsemblTranscript('7:g.140453136A>T'), h.namesEnsemblTranscript(null)],
    [true, true, false, false, false, false, false]);
    console.log('ensembl-archive: who may read a name — unit checks passed');
}

// ── The archive's transcripts beside MANE Select ───────────────────────────
{
    const row = (transcript, extra = {}) => ({ transcript, cDNA: 'c.35G>A', protein: 'p.Gly12Asp', gene: 'KRAS', source: 'vep', archive: true, ...extra });
    const mane = ['NM_004985.5', 'ENST00000311936.8'];
    const marked = h.markManeTranscripts([row('NM_004985.3'), row('NM_004985.4'), row('NM_033360.2'), row('NM_004985.5'),
        row('ENST00000311936.3'), row('NM_004985.4', { archive: undefined })], mane);
    // The MANE Select transcript in another version is not MANE Select — but it is that transcript.
    assert.deepEqual(marked.map(t => [t.transcript, t.maneSelect, t.maneSelectVersion]), [
        ['NM_004985.3', false, 'NM_004985.5'],
        ['NM_004985.4', false, 'NM_004985.5'],
        ['NM_033360.2', false, undefined],
        ['NM_004985.5', true, undefined],
        ['ENST00000311936.3', false, 'ENST00000311936.8'],
        // (Only of the archive's rows is this said: a record's own rows are left as they always were.)
        ['NM_004985.4', false, undefined]
    ]);
    assert.equal('maneSelectVersion' in marked[2], false);
    assert.deepEqual(h.markManeTranscripts([row('NM_004985.4')], []).map(t => [t.maneSelect, 'maneSelectVersion' in t]), [[false, false]], 'MANE not known');

    // It is the row shown: before a transcript with a lower accession, and the latest of several versions.
    const pick = rows => rows[h.selectCanonicalTranscript(rows)].transcript;
    const nf1 = h.markManeTranscripts([row('NM_000267.3', { canonicalSource: true }), row('NM_001042492.2'), row('NM_001128147.2')], ['NM_001042492.3']);
    assert.equal(pick(nf1), 'NM_001042492.2', 'NF1: the MANE Select transcript, not the lowest accession');
    assert.equal(pick(h.markManeTranscripts([row('NM_000267.3', { canonicalSource: true }), row('NM_001042492.2')], [])), 'NM_000267.3', 'as before, where MANE is not known');
    assert.equal(pick(marked.slice(0, 3)), 'NM_004985.4', 'the latest version there is');
    assert.equal(pick(h.markManeTranscripts([row('NM_004985.10'), row('NM_004985.9')], ['NM_004985.11'])), 'NM_004985.10', 'versions are numbers');
    assert.equal(pick(marked), 'NM_004985.5', 'MANE Select itself comes first');
    const typed = h.markManeTranscripts([row('NM_004985.4'), row('NM_033360.2', { protein: 'p.Gly12Val' })], mane);
    assert.equal(typed[h.selectCanonicalTranscript(typed, 'Gly12Val')].transcript, 'NM_033360.2', 'the protein change asked for still comes first');

    // What the card says beside such a row.
    assert.equal(h.ensemblArchiveNote(marked[1]), 'Named by Ensembl\'s GRCh37 archive (release 99, 2020): the current Ensembl service did not answer. '
        + 'NM_004985.4 is another version of the MANE Select transcript, NM_004985.5: where the two versions differ, so does the numbering.');
    assert.equal(h.ensemblArchiveNote(marked[2]), 'Named by Ensembl\'s GRCh37 archive (release 99, 2020): the current Ensembl service did not answer. '
        + 'Its transcripts are in the versions of that release, which can be numbered differently from today\'s, and it does not say which is MANE Select.');
    console.log('ensembl-archive: transcripts beside MANE Select — unit checks passed');
}

// ── GRCh38 coordinates are mapped by the archive too ───────────────────────
{
    const realFetch = globalThis.fetch;
    const mapped = { mappings: [{ original: { assembly: 'GRCh38', seq_region_name: '7', start: 140753336, end: 140753336, strand: 1 },
        mapped: { assembly: 'GRCh37', seq_region_name: '7', start: 140453136, end: 140453136, strand: 1, coord_system: 'chromosome' } }] };
    const path = '/map/human/GRCh38/7:140753336..140753336:1/GRCh37?content-type=application/json';
    try {
        const asked = [];
        globalThis.fetch = async (address) => {
            const url = String(address);
            asked.push(url);
            if (url.startsWith('https://rest.ensembl.org')) throw new TypeError('Failed to fetch');
            if (url.startsWith('https://grch37.rest.ensembl.org')) return new Response('<!doctype html>', { status: 500 });
            return new Response(JSON.stringify(mapped));
        };
        assert.equal(await h.mapGrch38GenomicToGrch37('chr7:g.140753336A>T'), 'chr7:g.140453136A>T');
        assert.deepEqual(asked, [`https://rest.ensembl.org${path}`, `https://grch37.rest.ensembl.org${path}`, `https://grch37-archive.rest.ensembl.org${path}`]);
        // (A host that answers is not second-guessed: the archive is for when the others do not.)
        asked.length = 0;
        globalThis.fetch = async (address) => { asked.push(String(address)); return new Response(JSON.stringify(mapped)); };
        assert.equal(await h.mapGrch38GenomicToGrch37('chr7:g.140753336A>T'), 'chr7:g.140453136A>T');
        assert.deepEqual(asked, [`https://rest.ensembl.org${path}`]);
        globalThis.fetch = async () => new Response('<!doctype html>', { status: 500 });
        await assert.rejects(h.mapGrch38GenomicToGrch37('chr7:g.140753336A>T'), /Ensembl assembly map returned 500/);
    } finally {
        globalThis.fetch = realFetch;
    }
    console.log('ensembl-archive: GRCh38 coordinates — unit checks passed');
}

const archiveNote = page => [...page.document.querySelectorAll('.archive-transcript-note')].map(el => el.textContent);
const vepRequests = (page, host) => page.requests.filter(r => r.startsWith(`${host}/vep/human/hgvs/`)).map(r => r.slice(`${host}/vep/human/hgvs/`.length));
const ofArchive = page => vepRequests(page, 'grch37-archive.rest.ensembl.org');
const ofCurrent = page => vepRequests(page, 'grch37.rest.ensembl.org');
const sets = { refseq: '?content-type=application/json&hgvs=1&refseq=1&mane=1&canonical=1', ensembl: '?content-type=application/json&hgvs=1&mane=1&canonical=1' };
const said = 'answered by Ensembl\'s GRCh37 archive (release 99) — the current service did not answer; transcript versions are the archive\'s';
// (The archive down as well: nobody at Ensembl answers.)
const archiveDown = [/grch37-archive\.rest\.ensembl\.org/];
const outage = 'Error: Lookup could not complete: Resolve nomenclature · Ensembl, Consequences · Ensembl VEP did not respond. '
    + 'This is an upstream outage, not necessarily a problem with your variant — please retry in a moment.';

// 1. The row, as recorded: every request to the current service answered HTTP
//    500, MyVariant holding no record under the allele's id. Was: "Lookup could
//    not complete". The archive names it, and the card and the panel say it did.
{
    const page = await runLookup(cassettes.tp53InsArchive, {});
    assertResolved(assert, page, '17 7576910 TP53 G GGT');
    assert.match(page.variant, /g\.: chr17:g\.7576910_7576911insGTVCF \(hg19\): 17-7576910-G-GGTGene: TP53Transcript: NM_000546\.5Named by/);
    assert.match(page.variant, /c\.: c\.934_935dupp\.: p\.Ser313ProfsTer33 \(p\.S313Pfs\*33\)Effect: frameshift_variant/);
    assert.deepEqual(archiveNote(page), ['Named by Ensembl\'s GRCh37 archive (release 99, 2020): the current Ensembl service did not answer. '
        + 'NM_000546.5 is another version of the MANE Select transcript, NM_000546.6: where the two versions differ, so does the numbering.']);
    assert.doesNotMatch(page.variant, /\(MANE Select\)/, 'a version MANE Select is not, is not called it');
    assert.match(page.panel, /^Variant resolved/);
    assert.match(page.panel, /Coding\s*NM_000546\.5:c\.934_935dup\s*Protein\s*p\.Ser313ProfsTer33\s*Effect\s*frameshift_variant/);
    assert.match(page.panel, /✕Resolve nomenclature · Ensembl\s*HTTP 500/, 'what did not answer is still shown as that');
    assert(page.panel.includes(`▲Consequences · Ensembl VEPframeshift_variant · ${said}`), page.panel);
    assert.deepEqual(ofArchive(page), [`17:g.7576910_7576911insGT${sets.refseq}`], 'asked once, of its RefSeq set');
    assert(ofCurrent(page).length >= 1, 'the current service was asked first');
    assert.deepEqual(page.api('civic'), ['civic?gene=TP53&protein=Ser313ProfsTer33']);
    console.log('ok ensembl archive: with Ensembl\'s GRCh37 service down, the TP53 row is named by its archive — and says so');
}

// 2. Nobody at Ensembl answering: the outage it is, in the words it had.
{
    const page = await runLookup(cassettes.tp53InsArchive, { without: archiveDown });
    assert.equal(page.status, outage);
    assert.equal(page.hidden, true);
    assert.match(page.panel, /^Upstream service unavailable/);
    assert.match(page.panel, /✕Consequences · Ensembl VEP\s*HTTP 500/);
    assert.deepEqual(archiveNote(page), []);
    assert(ofArchive(page).length >= 2, 'the archive was tried');
    // ...and the same when what the archive sends is no answer.
    const empty = await runLookup(cassettes.tp53InsArchive, { route: (url, key, respond) => (/grch37-archive.*\/vep\//.test(key) ? respond([]) : undefined) });
    assert.equal(empty.status, outage);
    const refused = await runLookup(cassettes.tp53InsArchive, { route: (url, key, respond) => (/grch37-archive.*\/vep\//.test(key)
        ? respond({ error: 'Unable to parse HGVS notation' }, 400) : undefined) });
    assert.equal(refused.status, outage, 'the archive\'s refusal is not taken for the current service\'s');
    console.log('ok ensembl archive: with the archive down as well, it is the outage it was');
}

// 3. A variant with a record (KRAS G12D): the record's own rows used to be all
//    that was left — on whatever versions it was filed with, protein names
//    without their accession. The archive's rows are VEP's; of two versions of
//    the MANE Select transcript the later is the one shown.
{
    const page = await runLookup(cassettes.krasG12dArchive, {});
    assertResolved(assert, page, '12 25398284 KRAS C T');
    assert.match(page.variant, /Gene: KRASTranscript: NM_004985\.4Named by.*c\.: c\.35G>Ap\.: p\.Gly12Asp \(p\.G12D\)Effect: missense_variant/);
    assert.deepEqual(archiveNote(page), ['Named by Ensembl\'s GRCh37 archive (release 99, 2020): the current Ensembl service did not answer. '
        + 'NM_004985.4 is another version of the MANE Select transcript, NM_004985.5: where the two versions differ, so does the numbering.']);
    assert.match(page.variant, /NM_004985\.3: c\.35G>A, NP_004976\.2:p\.Gly12Asp — missense_variant/);
    assert(page.panel.includes(`▲Consequences · Ensembl VEP5 transcripts · ${said}`), page.panel);
    assert.match(page.panel, /●Annotations · MyVariant\.info\s*annotation found/);
    // Without the archive: the record's rows, as before, and nothing said of an archive.
    const without = await runLookup(cassettes.krasG12dArchive, { without: archiveDown });
    assertResolved(assert, without, 'KRAS G12D row, nobody at Ensembl answering');
    assert.match(without.variant, /Transcript: NM_004985\.4c\.: c\.35G>Ap\.: p\.Gly12Asp/);
    assert.deepEqual(archiveNote(without), []);
    assert.doesNotMatch(without.panel, /archive/);
    console.log('ok ensembl archive: a variant with a record gets VEP\'s rows from the archive, on the latest version it has of MANE Select');
}

// 4. A name. NCBI places "KRAS c.34_36dup" (Ensembl's recoder being down too),
//    and the allele is put to the archive — both ways it is written, as for
//    any in-frame insertion on the reverse strand. Was: "Lookup could not
//    complete: Consequences · Ensembl VEP did not respond".
{
    const page = await runLookup(cassettes.krasDupArchive, {});
    assertResolved(assert, page, 'KRAS c.34_36dup');
    assert.match(page.variant, /g\.: chr12:g\.25398286_25398287insCCAVCF \(hg19\): 12-25398280-G-GCCA \(left-aligned\)Gene: KRASTranscript: NM_004985\.4Named by/);
    assert.match(page.variant, /c\.: c\.36_38dupp\.: p\.Gly13dup \(p\.G13dup\)Effect: inframe_insertion/);
    assert.match(page.panel, /Parse input\s*read as KRAS:c\.34_36dup · on NM_004985\.5, the MANE Select transcript of KRAS/);
    assert.match(page.panel, /Resolve nomenclature · Ensembl\s*variant placed by NCBI Variation Services — Ensembl had not answered/);
    assert.match(page.panel, /reverse-strand transcripts asked as chr12:g\.25398280_25398281insCCA — the same allele, left-aligned/);
    assert(page.panel.includes(said), page.panel);
    assert.deepEqual(ofArchive(page).sort(), [`12:g.25398280_25398281insCCA${sets.refseq}`, `12:g.25398286_25398287insCCA${sets.refseq}`]);
    // (The name itself is not put to the archive, which has no NM_004985.5 to read it on.)
    assert(!page.requests.some(r => /grch37-archive.*(?:NM_|KRAS:)/.test(r)), 'only the allele is asked of it');
    const without = await runLookup(cassettes.krasDupArchive, { without: archiveDown });
    assert.equal(without.status, 'Error: Lookup could not complete: Consequences · Ensembl VEP did not respond. '
        + 'This is an upstream outage, not necessarily a problem with your variant — please retry in a moment.');
    console.log('ok ensembl archive: a name NCBI places is named by the archive, both ways an in-frame insertion is written');
}

// 5. A name read on a transcript that is not MANE Select (NF1's MANE Plus
//    Clinical, NM_000267.4): the row for it is the archive's NM_000267.3, the
//    same transcript — and the note does not call it a version of MANE Select.
{
    const page = await runLookup(cassettes.nf1DelArchive, {});
    assertResolved(assert, page, 'NF1 c.6789_6792del');
    assert.match(page.variant, /NF1 c\.6789_6792del is read on NM_000267\.4, a MANE Plus Clinical transcript of NF1\./);
    assert.match(page.variant, /Gene: NF1Transcript: NM_000267\.3Named by.*c\.: c\.6789_6792delp\.: p\.Tyr2264ThrfsTer5 \(p\.Y2264Tfs\*5\)Effect: frameshift_variant/);
    assert.deepEqual(archiveNote(page), ['Named by Ensembl\'s GRCh37 archive (release 99, 2020): the current Ensembl service did not answer. '
        + 'Its transcripts are in the versions of that release, which can be numbered differently from today\'s, and it does not say which is MANE Select.']);
    assert.match(page.variant, /NM_001042492\.2: c\.6852_6855del, NP_001035957\.1:p\.Tyr2285ThrfsTer5/, 'the MANE Select transcript, in the archive\'s version, is listed');
    console.log('ok ensembl archive: a name read on another transcript keeps that transcript');
}

// 6. The mitochondrial genome: Ensembl's own transcript set, on the archive as on the current service.
{
    const page = await runLookup(cassettes.mtAtp6Archive, {});
    assertResolved(assert, page, 'MT-ATP6 m.8993T>G');
    assert.match(page.variant, /g\.: chrMT:g\.8993T>G.*Gene: MT-ATP6Transcript: ENST00000361899\.2Named by.*c\.: c\.467T>Gp\.: p\.Leu156Arg \(p\.L156R\)Effect: missense_variant/);
    assert.deepEqual(ofArchive(page), [`MT:g.8993T>G${sets.ensembl}`]);
    assert(page.panel.includes(said), page.panel);
    console.log('ok ensembl archive: a mitochondrial variant is named on Ensembl\'s own transcripts');
}

// 7. The current service silent rather than failing — a request that neither
//    answers nor ends. The archive is asked when the head start is over, and
//    the result does not wait for the current service to be given up on.
{
    let held = 0;
    const page = await runLookup(cassettes.tp53InsArchive, { timerCap: 60, settleMs: 600,
        route: (url, key) => {
            if (!key.startsWith('grch37.rest.ensembl.org/vep/')) return undefined;
            held += 1;
            return new Promise(() => {});
        } });
    assertResolved(assert, page, 'TP53 row, the current service silent');
    assert.match(page.variant, /Transcript: NM_000546\.5Named by.*c\.: c\.934_935dupp\.: p\.Ser313ProfsTer33/);
    assert.equal(held, 1, 'its one request is still out');
    assert.deepEqual(ofArchive(page), [`17:g.7576910_7576911insGT${sets.refseq}`]);
    console.log('ok ensembl archive: a current service that does not answer is not waited out');
}

// 8. The current service answering: the archive is not asked, and nothing is said of it.
{
    for (const name of ['eif1axDup30', 'krasDup', 'cftrF508del', 'brafV600E', 'geneAliasRow', 'mtNd4R340H']) {
        const page = await runLookup(cassettes[name], {});
        assertResolved(assert, page, name);
        assert(ofCurrent(page).length > 0, `${name}: VEP is asked`);
        assert.deepEqual(ofArchive(page), [], `${name}: not of the archive`);
        assert.deepEqual(archiveNote(page), [], name);
        assert.doesNotMatch(page.panel + page.variant, /archive/i, name);
    }
    console.log('ok ensembl archive: while the current service answers, the archive is not asked');
}

// 9. The current service REFUSING the allele (HTTP 400: an insertion between
//    positions five apart). That is its answer: the archive is not asked, and
//    the refusal is told as before.
{
    const page = await runLookup(cassettes.refusedNames, {});
    assert.match(page.status, /^Error: chr12:g\.25398280_25398285insCCA is not a readable insertion: /);
    assert.deepEqual(ofArchive(page), []);
    // ...also where only VEP is left to refuse it (the recoder not answering).
    const vepOnly = await runLookup(cassettes.refusedNames, { without: [/variant_recoder/] });
    assert.match(vepOnly.status, /^Error: chr12:g\.25398280_25398285insCCA is not a readable insertion: /);
    assert.deepEqual(ofArchive(vepOnly), []);
    console.log('ok ensembl archive: a refusal of the current service is its answer — the archive is not asked');
}

// 10. A name on one of Ensembl's own transcripts: "ENST00000288602.6:c.1799T>A",
//     BRAF V600E on GRCh37's version of the transcript. With the mirror down
//     the name went to the main, GRCh38, host, which reads it on ITS version
//     (12) of that transcript, where the position is another base — and the
//     lookup said "names a reference base that is not there: ENST00000288602.6
//     has C at that position, not T". (That answer is in the recording, as the
//     main host gave it the same day.) The name is not sent there; the archive
//     has the same versions GRCh37 has, and reads it.
{
    const page = await runLookup(cassettes.enstBrafArchive, {});
    assertResolved(assert, page, 'ENST00000288602.6:c.1799T>A');
    assert.match(page.variant, /g\.: chr7:g\.140453136A>TVCF \(hg19\): 7-140453136-A-TGene: BRAFTranscript: NM_004333\.4Named by/);
    assert.match(page.variant, /c\.: c\.1799T>Ap\.: p\.Val600Glu \(p\.V600E\)Effect: missense_variant/);
    assert(!page.requests.some(r => r.startsWith('rest.ensembl.org/variant_recoder/')), 'the name is not read on the GRCh38 host');
    assert.deepEqual(ofArchive(page).map(r => r.split('?')[0]), ['ENST00000288602.6:c.1799T>A', '7:g.140453136A>T'], 'the name, then the allele it is');
    assert(page.panel.includes(said), page.panel);
    assert.doesNotMatch(page.status + page.panel, /reference base that is not there|has C there/);
    // With the archive down as well it is an outage, and is called one — not a fault in the name.
    const without = await runLookup(cassettes.enstBrafArchive, { without: archiveDown });
    assert.equal(without.status, outage);
    assert(!without.requests.some(r => r.startsWith('rest.ensembl.org/variant_recoder/')));
    console.log('ok ensembl archive: a name on an Ensembl transcript is read on GRCh37\'s version of it, not on the GRCh38 host\'s');
}

// 11. A transcript the archive counts a base off. PTEN's NM_000314.6 differs
//     from GRCh37 by one base in its 5' UTR; the VEP of that release names
//     "chr10:g.89725116del" c.1101del on it — and c.1100del on NM_000314.4, on
//     Ensembl's own transcript, and (NCBI, the current service) on NM_000314.8.
//     VEP flags the row itself, and it is left out: the later version is not
//     the one shown where it cannot be relied on.
{
    const page = await runLookup(cassettes.ptenDelArchive, {});
    assertResolved(assert, page, 'chr10:g.89725116del');
    assert.match(page.variant, /Gene: PTENTranscript: NM_000314\.4Named by.*c\.: c\.1100delp\.: p\.Pro367GlnfsTer49 \(p\.P367Qfs\*49\)Effect: frameshift_variant/);
    assert.deepEqual(archiveNote(page), ['Named by Ensembl\'s GRCh37 archive (release 99, 2020): the current Ensembl service did not answer. '
        + 'NM_000314.4 is another version of the MANE Select transcript, NM_000314.8: where the two versions differ, so does the numbering.']);
    assert.doesNotMatch(page.variant + page.panel, /NM_000314\.6|c\.1101del|NM_001304717|NM_001304718/, 'the rows that cannot be relied on are not shown at all');
    // (They are in what the archive sent.)
    const sent = cassettes.ptenDelArchive.responses[`grch37-archive.rest.ensembl.org/vep/human/hgvs/10:g.89725116del${sets.refseq}`].body[0].transcript_consequences;
    assert.deepEqual(sent.map(row => [row.hgvsc, (row.refseq_match || []).join(',')]), [
        ['NM_000314.4:c.1100del', ''],
        ['NM_000314.6:c.1101del', 'rseq_mrna_nonmatch,rseq_5p_mismatch'],
        ['NM_001304717.2:c.1619del', 'rseq_mrna_nonmatch,rseq_cds_mismatch'],
        ['NM_001304718.1:c.510del', 'rseq_mrna_nonmatch,rseq_5p_mismatch']
    ]);
    console.log('ok ensembl archive: a transcript the archive counts a base off (PTEN NM_000314.6) is left out');
}

console.log('Ensembl archive tests passed');
