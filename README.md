# Variant Search (Beta)

## Variant identity and transcript consistency

Exact genomic alleles are preserved through annotation fallback. A missing
MyVariant record is not evidence of a different genome assembly, and free-text
search hits are accepted only after matching the requested HGVS change and
gene/accession. Transcript fallback hits must also match the resolved genomic
allele. Candidate recoder alleles are filtered against the query, and a final
identity check rejects inconsistent coordinates or annotations — by exact
spelling for substitutions, and by allele for everything else (see
[Indels: one allele, several spellings](#indels-one-allele-several-spellings)).

The Variant card displays the selected transcript accession and keeps its
cDNA, protein and effect together. Ensembl recoder and dbNSFP arrays are not
paired by index. VEP and snpEff provide transcript-specific pairs; unavailable
protein annotations remain blank. Explicit transcript/cDNA queries are honored.
For genomic queries, MANE Select is preferred using a cached, bounded Ensembl
gene lookup. That lookup supplies transcript metadata only, never GRCh38 variant
coordinates. If MANE metadata is unavailable, the RefSeq accession fallback
includes UTR transcripts rather than preferring the longest coding isoform.
Ensembl's GRCh37 canonical flag is not treated as equivalent to MANE Select.

Regression example: both `PTEN c.-8T>G` and `10 89624219 PTEN T G` resolve to
`chr10:g.89624219T>G`, `NM_000314.8:c.-8T>G` (5′ UTR), and ClinVar 4757098.
The alternative `NM_001304717.5:c.512T>G / NP_001291646.4:p.Leu171Arg` is
shown with its own transcript. Missing annotation must never select ClinVar
189414, 1750073, or the same-position T>C allele 631126.

### Indels: one allele, several spellings

A substitution has exactly one spelling. An indel has several: HGVS shifts it
to the 3′-most position it can occupy and calls an insertion that repeats its
neighbour a duplication, VCF left-aligns, and Ensembl answers in whichever form
it was not asked in. A pasted row and Ensembl's answer for it routinely differ
as text while being the same variant:

```
X 20156693 EIF1AX C CAGATTCATTCTCATTCTTACCCCTGCGTCT
  queried as   chrX:g.20156693_20156694insAGATTCATTCTCATTCTTACCCCTGCGTCT
  resolved as  chrX:g.20156694_20156723dup   NM_001412.4:c.34_63dup   p.Arg12_Ser21dup
```

So spellings cannot decide identity for indels in either direction. Comparing
them as strings refused that lookup outright ("Variant identity mismatch:
annotation belongs to a different genomic allele") along with every other indel
that Ensembl re-spells and MyVariant does not index — while two insertions one
base apart, which *are* different alleles, passed.

`compareGenomicAlleles` settles it against the sequence instead: whenever a
genomic input comes back spelled differently, both descriptions are applied to
the same stretch of GRCh37 reference (Ensembl `/sequence/region`, one small
request, started as soon as the recoder answers) and the results compared. No
shifting or alignment heuristics are involved. The outcome is a row in the
lookup progress panel, **Confirm allele · GRCh37 reference**:

| Outcome | Panel | Effect |
| --- | --- | --- |
| same sequence | `●` "… is the same allele as the input" | lookup proceeds |
| different sequence | `✕` "… is not the allele entered" | hard error — never rendered |
| reference unreachable, or event wider than 10 kb | `▲` "not confirmed — …" | lookup proceeds on the GRCh37 mirror's mapping, without claiming a confirmation |

The identity checks read a few dozen bases and normally have them in well under
a second. They sit in front of the result, so they run on a 5-second leash
(`API_TIMEOUT_MS.identityReference`) and report "not checked" rather than hold
the page behind a hung Ensembl.

The recoder's SPDI output is part of the same story. An SPDI position is a
0-based *interbase* coordinate, so an insertion at `p` lies between 1-based
bases `p` and `p+1`. `convertSpdiToMyVariant` used to read it as `p+1`/`p+2`,
one base too far 3′ — a different allele unless the inserted bases repeat the
base they hop over. That candidate was then looked up, and adopted whenever
MyVariant happened to hold it.

Regression examples (live responses recorded in
`tests/fixtures/indel-identity.json`):

- The row above, and the documented `X 20148674 EIF1AX T TG` →
  `chrX:g.20148675dup` / `NM_001412.4:c.388dup`, resolve instead of failing.
- `13 32906698 BRCA2 T TCA` resolves to its own allele,
  `chr13:g.32906698_32906699insCA` / `NM_000059.4:c.1083_1084insCA` /
  `p.Asp362GlnfsTer6`. It must never again be answered with the neighbouring
  `g.32906699_32906700insCA` (`c.1084_1085insCA`, ClinVar 2505518).
- `chr7:g.117199645_117199647del` and `…646_…648del` are both CFTR F508del;
  `…644_…646del`, one base 5′, is I507del.
- `M 8993 T G` resolves to `chrMT:g.8993T>G` (ClinVar 9641). How a chromosome
  is *named* is not part of the allele either: `7:g.` is `chr7:g.`, and UCSC's
  `chrM` is the `chrMT` of Ensembl and MyVariant. Identity comparisons
  normalise both, so the hg19 spelling is no longer refused as a mismatch.
  (`chrM` is read as rCRS, NC_012920 — what MyVariant and Ensembl hold. The
  `chrM` of UCSC's *original* hg19 is the older NC_001807, numbered 1–2 bp
  differently in places.) Ensembl writes mitochondrial answers `m.` rather
  than `g.`; `M 16300 A C`, which MyVariant does not hold, resolves too.

### Reference bases the input states

A pasted row has a REF column and HGVS can spell its bases out (`delAAT`).
Neither used to reach the lookup: the row is reduced to its minimal event, and
Ensembl ignores spelled-out bases on anything but a substitution — it answers
`chr7:g.117199645_117199647delAAA` as F508del although GRCh37 has `TCT` there.
An indel entered in the wrong build therefore resolved to whatever sits at those
numbers in GRCh37: a GRCh38 TP53 deletion pasted without the GRCh38 box
(`17 7675087 TP53 GCGC G`) landed 97 kb away.

`statedReferenceBases` / `verifyStatedReference` read those bases off GRCh37
alongside the lookup. A mismatch is a hard error that names the position, both
sequences and the GRCh38 box:

> Reference mismatch: the input gives GCGC at chr17:7675087-7675090, but
> GRCh37/hg19 has TGCA there. If these coordinates are GRCh38/hg38, tick the
> GRCh38 box and look it up again; otherwise check the position and the
> reference allele.

Substitutions are left out: their REF is part of the MyVariant id, Ensembl
rejects a wrong one itself, and checking them would put an Ensembl request on
the commonest lookup there is. If the reference cannot be read, nothing is
refused on this ground.

### MyVariant's id is not always the allele ClinVar describes

`_id` is MyVariant's own construction, and for the records it builds from
ClinVar it is often not the allele ClinVar describes. Comparing the id with the
GRCh37 genomic HGVS in the *same record's* ClinVar block, over 200 ClinVar
records of each kind (BRCA1, BRCA2, TP53, MLH1, APC):

| ClinVar type | id names a different allele |
| --- | --- |
| insertion | 49 / 200 — position from one normalisation, inserted bases from another |
| duplication | 199 / 200 — written `g.A_Bdup` with the insertion *flanks*, which as HGVS duplicates those two bases instead |
| deletion, delins | 0 / 200 each |

The site read the id as genomic HGVS. BRCA1 c.5266dupC (5382insC, ClinVar 17677)
is filed under `chr17:g.41209079_41209080dup` while its ClinVar block says
`NC_000017.10:g.41209082dup`, so the page showed the right ClinVar record beside
**VCF `17-41209078-C-CTG`** — a 2-bp insertion — and asked gnomAD and SpliceAI
about that allele rather than `17-41209079-T-TG`.

`reconcileClinvarAllele` compares the two against the reference whenever an
insertion or duplication id is spelled differently from its ClinVar HGVS, and
puts the record right for whatever the lookup was for:

- **Found by name through its ClinVar entry** (a c., p. or gene-symbol query):
  re-keyed to ClinVar's allele, keeping only the ClinVar block. Everything else
  in the record — vcf, snpEff, population frequencies, dbSNP — was merged in by
  id and describes the other allele. The progress panel says so: *"allele
  chr17:g.41209082dup from ClinVar's own HGVS (MyVariant files it under
  chr17:g.41209079_41209080dup)"*, and the raw-JSON panel keeps the original id
  as `_myvariant_id`.
- **Fetched as the id's own allele** (a genomic query, or a recoder candidate):
  the ClinVar block belongs to another variant and is removed.
  `13 32906699 BRCA2 G GCA` is literally the allele MyVariant files ClinVar
  2505518 under, and is not that variant.
- **Reference unreachable**: Ensembl's sequence endpoint is the one this
  comparison depends on, and it does go down (it returned HTTP 500 for part of a
  live run of these very cases). "Cannot tell" used to be the answer for every
  pair, so `13 32906699 BRCA2 G GCA` was shown ClinVar 2505518 — *Likely
  pathogenic* — after all, and `BRCA1 c.5266dup` came out as
  `g.41209079_41209080dup`. Two things need no reference
  (`provablyDifferentAlleles`), and between them they cover both ways MyVariant
  mis-files a record:
  - alleles that change the sequence length by different amounts differ — a
    2-bp duplication is not a 1-bp one;
  - an insertion slid *k* positions along a repeat has its bases rotated by *k*,
    so `insCA` one position on is `insAC` if it is anything, never `insCA`.

  Anything else — an insertion against the duplication it may be, two deletions
  of one length — is still left as it was, and the panel says the id could not be
  checked. The proof only ever says "different": on random repetitive sequence it
  is checked never to fire for two spellings the reference calls one allele.

A re-keyed record also takes its gene from ClinVar's preferred name. MyVariant
keeps one of the genes ClinVar lists, and for EGFR exon 20 that is the antisense
gene: `NM_005228.5:c.2300_2308dup` was labelled `EGFR-AS1` and sent to CIViC,
PubMed and the other gene cards as that.

### One allele, several MyVariant records

MyVariant keys a record by the spelling its source used, and its sources spell
an indel in a repeat differently. One allele is therefore several records, each
holding only what one family of sources knew:

| | id | holds |
| --- | --- | --- |
| BRCA1 c.5266dup (5382insC) | `chr17:g.41209082dup` — 3′-most, HGVS | dbSNP: rs80357906 |
| | `chr17:g.41209079_41209080insG` — left-aligned, VCF | gnomAD, ExAC, snpEff, vcf |
| | `chr17:g.41209079_41209080dup` — MyVariant's own | ClinVar: 17677 |
| CFTR F508del | `chr7:g.117199645_117199647del` — left-aligned | ClinVar, gnomAD, ExAC, snpEff |
| | `chr7:g.117199646_117199648del` — 3′-most | dbSNP: rs113993960 |

A lookup showed whichever record it landed on:

| Input | was missing | now |
| --- | --- | --- |
| `rs80357906` | ClinVar ("Clinical significance: N/A"), gnomAD | ClinVar 17677 · Pathogenic (14) · gnomAD exome AF 0.0183% |
| `BRCA1 c.5266dup` | gnomAD, rsID | both |
| `17 7578474 TP53 C CG` | ClinVar, rsID | ClinVar 182958 · rs730882019 |
| `7 117199644 CFTR ATCT A` (F508del) | rsID | rs113993960 |
| `chr13:g.32906700_32906701insCA` | everything — no record under either real spelling | ClinVar 2505518 |

`gatherAlleleRecords` collects the rest. The other ids are **derived, not
guessed**: the reference says how far the indel slides, and the two ends of
that stretch are the left-aligned and the 3′-most placement
(`indelPlacements`). On top of those come the two ids MyVariant makes up for
ClinVar records, both measured — over ClinVar records in BRCA1, BRCA2, TP53,
MLH1 and APC, 80 of 80 duplications are filed as `g.<left>_<left+1>dup` and 49
of 49 mis-filed insertions as `g.<left>_<left+1>ins<3′-most bases>`. All the
candidates go to MyVariant in one batch request.

**Nothing is merged on the strength of a spelling.** Each record that comes back
is applied to the reference and must come out as the allele being looked up
(`blocksForAllele`). The ids MyVariant makes up read, literally, as *other*
alleles, so from those only the ClinVar block is taken, and only when ClinVar's
own HGVS inside it is this allele — the dbSNP block beside ClinVar 17677 is
rs2051500205, the 2-bp duplication its id spells, and is left where it is.
The same rule keeps `13 32906699 BRCA2 G GCA` from inheriting ClinVar 2505518,
which MyVariant files under exactly that id. A block the lookup already has is
never replaced; what was borrowed, and from where, is listed in the progress
panel (`also filed under chr17:g.41209079_41209080insG (gnomAD, ExAC, snpEff,
vcf), chr17:g.41209079_41209080dup (ClinVar)`) and in the raw JSON under
`_also_filed_under`.

Without the reference nothing can be derived or proven, so nothing is merged and
the panel says `its other ids for this allele were not checked`. The step runs
alongside the wait for Ensembl and is given 6 s; a substitution has one id and
costs nothing.

A related slip, fixed with it: links that took REF/ALT from the record's `vcf`
block paired them with the first base of the HGVS event rather than the position
that block is anchored at — F508del typed as HGVS went to gnomAD as
`7-117199645-ATCT-A` (it is at 117199644).

**The VCF alleles are left-aligned whichever record the lookup landed on.**
MyVariant writes a record's `vcf` block for the id it sits under, so the record
under the 3′-most id is anchored there, and that is not the spelling
VCF-derived resources use:

| Reached through | VCF alleles were | gnomAD v4.1 said | now |
| --- | --- | --- | --- |
| `chr7:g.117199646_117199648del`, `rs113993960` (F508del) | `7-117199645-TCTT-T` | Not found | `7-117199644-ATCT-A` — 18,037 / 1,460,114 |
| `BRCA1 185delAG`, `NM_007294.4:c.68_69del` | `17-41276046-TCT-T` | Not found | `17-41276044-ACT-A` — 169 / 1,459,934 |

(The gnomAD v2 link opened a "variant not found" page for the same reason.) The
reference is already in hand when the records are gathered, so the block is
re-anchored there (`leftAlignedVcfBlock`) — only if, left-aligned, it is the
allele being looked up — and every consumer sees the one spelling: the VCF line,
both gnomAD halves, SpliceAI and the AI review. What MyVariant has is kept in the
raw JSON as `_vcf_as_filed`.

See `tests/sibling-records.test.js`.

### Which variant a gene-symbol name means

A gene symbol does not choose a transcript, and c./p. numbering differs between
transcripts. `BRCA1 c.5266dup` is c.5266dup on *some* RefSeq transcript for six
different ClinVar variants. MyVariant's search ranks 5382insC (ClinVar 17677,
`NM_007294.4:c.5266dup`) fourth; the first hit — ClinVar 266551, which on the
standard transcript is `c.5389dup` — used to be presented as the answer.

`selectMyVariantHit` ranks the matching hits by name (`variantNameRank`): a hit
whose ClinVar *preferred name* carries the queried change first (ClinVar writes
that name on the MANE Select transcript), then the lowest-numbered `NM_` the
change matched on, then the existing ClinVar-presence tie-break. It is a
preference rather than a filter, an explicit accession is honoured as before,
and protein names are ranked the same way.

**A coding name on a gene symbol is read on the gene's MANE Select transcript.**
When ClinVar itself names a hit `GENE c.X` that is the end of it, and Ensembl's
recoder is not asked at all (it does not answer `BRCA1:c.5266dup` within a
minute). Otherwise the ranking above was only a fallback — and with no hit at
all the name went to Ensembl, which picks a transcript of its own:

| `PTEN c.512T>G` | |
| --- | --- |
| before | `chr10:g.89720934T>G`, which is **c.1026+59T>G — an intron** — on every RefSeq transcript of PTEN, presented as the answer |
| now | *PTEN c.512T>G is not a change on NM_000314.8, the MANE Select transcript of PTEN: the reference base at that position is A, not T. The name may be counted on another transcript — enter it with that transcript's accession (NM_001304717.5:c.512T>G — the base there is T).* |

So the lookup first settles which transcript is standard
(`standardTranscriptForGene`) and from there on is for `NM_…:c.X`, which every
later step treats exactly; the progress panel says `read as
BRCA2:c.1083_1084insCA · on NM_000059.4, the MANE Select transcript of BRCA2`.
Two sources answer that question:

- **Ensembl's MANE record**, the authority. It is now fetched as two small
  requests — the gene, for its canonical transcript, then that one transcript.
  The single "expanded" request used before is 236 KB for BRCA1 and 818 KB for
  TTN, took 18 s and 13 s when measured, and so missed its 4 s deadline for
  exactly the large, much-studied genes: that is also why "(MANE Select)" was
  so often missing from the Variant card.
- **The transcript ClinVar names the gene's variants on**
  (`fetchClinvarNamingTranscript`, a 0.3 s MyVariant query), used when Ensembl
  has not answered in 3 s. Checked against MANE for 20 genes — BRCA1/2, TP53,
  PTEN, EGFR, KRAS, CDKN2A, MUTYH, TTN and others — it was the MANE Select,
  version included, in all 20.

When a service *rejects* the name (HTTP 400: it read it and said no) that is
reported as what it is. A reference base that is not there used to come out as
"Ensembl failed. Please retry in a moment."; it now names the base that is
there, for a gene-symbol name and for an accession typed with the wrong base
alike, and suggests transcripts on which MyVariant's other records show the
stated base to be right.

If neither source can say which transcript is standard, the old behaviour is all
that is left — and is said: the panel warns that the name was read on whichever
transcript matched first, and a typed substitution that is on none of the
transcripts shown gets a note beside the `c.` line.

See `tests/gene-symbol-name.test.js`.

### A genomic query names no gene

What precedes `:g.` is a chromosome or an accession. Reading it as a symbol
rendered `Gene: CHR16` for `chr16:g.2122948_2122950del`, filtered VEP's
transcript rows on "CHR16" — losing the protein change and the effect — and sent
"CHR16" to every gene-level source. When no gene column was pasted the gene now
comes from VEP.

For a variant that lies between genes VEP lists every gene within 5 kb, in
transcript order, and the first used to be taken: `m.16300A>C`, 277 bp upstream
of MT-TP, was labelled **MT-ND4** — 4 kb away, and a gene it lies *downstream*
of — beside the effect `upstream_gene_variant`. The gene named is now one the
reported effect is about, and of several, the nearest
(`resolveGeneSymbolFromVep`). A gene the variant is actually in still comes
before any neighbour.

### Mitochondrial variants

`m.8993T>G` — MT-ATP6 p.Leu156Arg, NARP/Leigh syndrome — was shown as
`Gene: MT-ATP6, MT-ATP8, RNR2-END · Transcript: RNR2-END · c.: n.8993T>G ·
Effect: intergenic_region`. Each source the page read is wrong about the
mitochondrial genome in its own way, and each needed its own answer:

| Source | What it says | What the page does |
| --- | --- | --- |
| snpEff (in MyVariant) | No gene models beyond the two rRNAs: everything else is `intergenic_region` against a made-up feature, `RNR2-END` (nine of nine well-known variants checked) | The block is set aside for `chrMT` records (`withoutUnusableBlocks`); the reason is left in the raw JSON under `_not_used` |
| ClinVar's gene (in MyVariant) | ClinVar's list is whatever submitters sent — 26 genes for m.8993T>G — and MyVariant keeps one: MT-ATP8 for that variant, MT-ND1 for m.1555A>G in MT-RNR1 | Used for the gene line only when CADD, which works from the position, names none |
| VEP, RefSeq transcripts | RefSeq has no accessions for the 37 mitochondrial genes; VEP invents them — gene `ATP6`, transcript `ATP6.1` | VEP is asked for Ensembl's own set on `MT`: `MT-ATP6`, `ENST00000361899` |
| Ensembl recoder | Answers `NC_012920.1:m.8993T>G`; cannot parse MyVariant's `chrMT:g.…` | Asked as `MT:g.…`; `m.` on NC_012920 is read as the same allele (the sequence is identical in both assemblies, so it is never lifted) |

The result is `Gene: MT-ATP6 · ENST00000361899.2 · c.467T>G · p.Leu156Arg ·
missense_variant`; m.1555A>G is `MT-RNR1 n.908A>G` and m.3243A>G `MT-TL1
n.14A>G`. With Ensembl down nothing has computed a coding change, so none is
shown — and CADD still names the gene.

Every common spelling reads as the same query, `chrMT:g.8993T>G`:

```
M 8993 T G        NC_012920.1:m.8993T>G        m.8993T>G        MT-ATP6 m.8993T>G
chrM:g.8993T>G    chrMT:g.8993T>G              MT:g.8993T>G
```

Only the bare `m.` form and the rows used to work. ClinVar's own name for the
variant, `NC_012920.1:m.8993T>G`, was read as the protein change `p.M.8993T>G`,
and `chrM` — the hg19 name — was asked of MyVariant as such, a guaranteed miss
that was recovered only when Ensembl answered and re-spelled it.

**gnomAD.** gnomAD keeps mitochondrial variants in a call set of their own,
counted by heteroplasmy, which neither half of the gnomAD card read. v2 has
none, and the v4 proxy looks up nuclear variants only: asked about one on `M` it
attempted a liftover, Ensembl rejected the chromosome name, and the card said

> gnomAD v4.1 unavailable: the Ensembl liftover service is not responding, so
> M:8993 could not be converted to GRCh38. This is usually temporary — try again
> shortly.

for every mitochondrial variant, every time (typed as `chrMT:g.…` the name was
accepted and the answer was *"Not found in gnomAD v4.1"*, for m.1555A>G with its
63 homoplasmic carriers) — beside a "View on gnomAD (v2.1.1)" link to a page
reading *"Mitochondrial variants are not available in gnomAD v2"*. The card now
reads the call set from gnomAD's own CORS-open API, the one gene constraint
already comes from (`fetchGnomadMitochondrial`; no liftover — the sequence is the
same in both builds):

```
Mitochondrial call set (gnomAD API)
Homoplasmic AF: 0.1117%  AC/AN: 63/56401
Heteroplasmic AF: 0.0195%  AC/AN: 11/56401
Homoplasmic: in at least 95% of a carrier's mitochondrial copies; heteroplasmic: in 10–95%. Highest level seen: 100%.
View on gnomAD (mitochondrial)
```

gnomAD not answering is reported as *could not be read*, never as *not in
gnomAD*, and the AI review receives the same counts with what they mean.

Two things found here apply to every chromosome:

- **A genomic query's transcripts survive a re-spelling.** The record for
  m.3271del is filed as `g.3271del`; Ensembl answers `g.3273del`, the same
  deletion at the other end of a run of T's. The recoder's transcript rows were
  matched to the allele by spelling, found no match and were discarded. A
  genomic query has exactly one allele in the recoder's answer — its own — so
  the rows are kept; a name (`c.`, `p.`) can still come back as several alleles
  and is still matched by id.
- **snpEff's protein-structure contacts are not transcripts.** TP53 R213*
  carries 762 `structural_interaction_variant` rows
  (`4QO1:B_169-B_213:NM_001126114.2`), each of which was listed under "Other
  transcripts".

See `tests/mitochondrial.test.js`.

### rsID lookups

An rsID used to be handed to Ensembl's variant recoder and waited for. Measured
live, with this page's own timeouts:

| rsID | before | now |
| --- | --- | --- |
| `rs80357906` (BRCA1 5382insC) | 49.5 s, then `c.: N/A · p.: N/A · Effect: delins` | 6.5 s · NM_007294.4 c.5266dup · p.Gln1756ProfsTer74 |
| `rs113488022` (BRAF V600E) | **p.Val600Ala** | p.Val600Glu, and a note naming the other two alleles |
| `rs28934578` (TP53 R175H) | 51.8 s, then `p.: N/A · Effect: N/A` | 8.9 s · p.Arg175His |
| `rs1042522` (TP53 P72R) | 52.1 s, then `p.: N/A · Effect: N/A` | 6.1 s · p.Pro72Arg |

Three separate things were wrong.

**The wait.** Ensembl takes 5 to 50 s over an rsID. The page allows 12 s, tried
four times, and then fell back to the MyVariant search it could have started
with. MyVariant is now asked first: it says at once which alleles it files under
the rsID. The recoder keeps running in the background and is no longer needed.

**Which allele.** An rsID names a site, not an allele, and rs113488022 has
three. The one shown was the first that Ensembl listed — and Ensembl lists them
in a different order from one call to the next (`T, G, C` in one recording of
this very query, `G, T, C` in another), so the rsID everyone uses for V600E came
back as V600A or V600E at random. The alleles MyVariant holds are now ordered by
how well documented each is — ClinVar submissions, then population frequency
(`rankRsidAlleles`) — and because that is a choice made on the user's behalf,
the Variant card says so and links the others:
`rs113488022 names 3 alleles. This is the best documented one; the others are
chr7:g.140453136A>C, chr7:g.140453136A>G.`

**No transcript.** The record an rsID lands on can be one of MyVariant's
dbSNP-only ones (an rsID and a gene, nothing else), and VEP was only asked when
the recoder had supplied transcripts to pair — which, having timed out, it had
not. VEP works from the allele alone and answers in a second or two, so it is
now asked as soon as the allele is known, for every lookup; once it has
answered, the recoder is not waited for. When the recoder contributed nothing,
VEP is given 8 s from when it was asked (`API_TIMEOUT_MS.vepUnprompted`) and the
result is shown without it after that. VEP answers from Ensembl's own transcript
set when its RefSeq one fails; asked unprompted, such an answer does not
displace RefSeq rows the record already has (during an Ensembl blip
`13 32906699 BRCA2 G GCA` came out on `ENST00000544455.1` while the record's
`NM_000059.3` row went unused), and the panel says *"RefSeq transcripts
unavailable — the record's own are shown"*.

Not waiting for the recoder means it answers later — typically a second or two
after VEP, in the middle of what remains of the lookup. Its rows carry no protein
change, and arriving then they replaced VEP's: in a real browser F508del by its
g. name came out `c.1521_1523del · p.: N/A · Effect: N/A` whenever Ensembl's
recoder took 3 s rather than 1 (and a late *failure* emptied the rows
altogether). Once the rows a result is built from have been decided, the
recoder's outcome is only reported: *"answered after the result was ready — not
needed: the variant was resolved without it"*.

Without VEP the record's own rows are used, and there a bare `c.` from ClinVar's
list no longer displaces snpEff's row for the same transcript: ClinVar lists the
current version (NM_000546.6) beside snpEff's older one (.5), the bare row was
selected for being the MANE version, and that is where `p.: N/A, Effect: N/A`
came from.

See `tests/rsid-lookup.test.js`. Its recordings were made with the page's real
timeouts, so a request Ensembl did not answer in time replays as a timeout.

`npm test` includes real-response identity fixtures and full form/card regression
checks using LinkeDOM, with external requests mocked. It covers both PTEN input
paths, explicit transcript selection, VEP/MANE outages, BRAF V600E while the
recoder is unavailable, and every case in the sections above
(`tests/indel-identity.test.js`: thirteen end-to-end scenarios, each of which
fails on the code before its fix). LinkeDOM is a development-only test
dependency.

Deployment for this frontend fix: upload **`script.js` and `index.html` together**
to the Hostinger `/variant/` directory after merging. The HTML contains the
matching script content hash to invalidate caches. No API-handler changes or
stylesheet upload are required.

## Accepted input formats

Coordinates are interpreted as **GRCh37/hg19 by default** — that is the assembly
every resolution endpoint in this app targets (`grch37.rest.ensembl.org`), so an
hg38 position passed through unconverted would either miss or silently resolve
to the wrong locus. GRCh38 input is supported through an explicit toggle — see
[GRCh38 input toggle](#grch38-input-toggle).

HGVS is accepted in genomic, coding and protein form:

```
chrX:g.20148675dup      NM_001412.4:c.388dup      EIF1AX:p.Gln130ProfsTer3
```

So are an rsID (`rs113488022` — see [rsID lookups](#rsid-lookups) for a site
with several alleles), mitochondrial `m.` notation (`m.8993T>G`,
`NC_012920.1:m.8993T>G`), a stop written the short way (`TP53 R213*`,
`R213X`), and the older names ClinVar lists a variant under (`BRCA1 5382insC`,
`BRCA1 185delAG` — see [Names that are not HGVS](#names-that-are-not-hgvs)).

Pasted coordinate rows are also accepted, in the shapes people actually paste
out of VCFs, MAFs, spreadsheets and pipeline reports:

| Input | Notes |
| --- | --- |
| `X  20148674  EIF1AX  T  TG` | chrom, pos, gene, ref, alt — tabs or spaces |
| `X  20148674  T  TG` | gene column optional |
| `EIF1AX  X  20148674  T  TG` | MAF column order (Hugo_Symbol first) |
| `X,20148674,EIF1AX,T,TG` | CSV/semicolon/pipe separated, quoted cells fine |
| `X  20,148,674  EIF1AX  T  TG` | thousands separators in the coordinate |
| `X  20148674  EIF1AX  T  TG  hg19` | trailing build/zygosity columns ignored |
| `X  20148674  EIF1AX  -  G` | MAF-style `-` for an absent allele |

Parsing anchors on the chromosome/position pair wherever it sits in the row and
reads the alleles from the columns that follow. When more than two DNA-like
columns follow the position, the **last two** are the alleles — a gene symbol
that happens to spell DNA (`TTN`, `CAT`, `TAT`) is always written before them.

A row that looks like coordinates but cannot be parsed is rejected with a
message naming the expected layout. It is deliberately **not** passed through to
the `GENE p.Change` interpretation: doing so used to turn
`EIF1AX X 20148674 T TG` into the query `EIF1AX:p.X` and then report it back as
a variant that could not be found.

VCF-anchored indels are reduced to the minimal HGVS event before lookup, so
`X 20148674 EIF1AX T TG` is queried as `chrX:g.20148674_20148675insG` and
resolves to `chrX:g.20148675dup` / `NM_001412.4:c.388dup`.

See `tests/coordinate-row.test.js`.

## GRCh38 input toggle

A checkbox under the search box — **"My coordinates are GRCh38/hg38"** — lets
genomic-coordinate input (chr:g. HGVS, pasted VCF/MAF rows) be entered in
GRCh38. Unchecked is the long-standing GRCh37/hg19 default. The choice is kept
in the shareable URL as `?assembly=grch38` (also accepted: `hg38`, `38`).

How it works:

- When checked, the normalised genomic query is mapped **once, up front**
  through Ensembl's `/map/human/GRCh38/…/GRCh37` endpoint (with the GRCh37
  mirror as a failover host), and the entire pipeline then runs on the GRCh37
  coordinate exactly as if it had been typed. The mapping appears as its own
  step in the lookup progress panel.
- The mapping is **strict**: it refuses coordinates that land on a different
  chromosome (patch/scaffold), change span length, flip strand, or split into
  multiple segments — and a refused or failed mapping is a hard error, not a
  silent fall-through. An hg38 position labelled hg19 would confidently
  annotate the wrong locus, which is worse than no answer.
- An `NC_` accession declares its own assembly through its version
  (`NC_000007.13` = GRCh37, `.14` = GRCh38) and **overrides the checkbox in
  both directions** — the accession is unambiguous, the checkbox can be stale.
  `NC_` genomic HGVS of either build is converted to the chr form and, when
  GRCh38, lifted the same way.
- Transcript- and protein-level input (`NM_…:c.…`, `GENE:p.…`, gene-only) has
  no assembly and is untouched by the checkbox.
- Because the lifted input pins down both builds' positions, the GRCh38 start
  is handed to the gnomAD v4 proxy as its `pos38` hint, skipping the
  server-side liftover round-trip for that variant.

Only positions are remapped; the event description (`A>T`, `del`, `insG`, …)
is carried over byte-for-byte. Helpers: `splitGenomicHgvsPositions`,
`ncGenomicToChr`, `mapGrch38GenomicToGrch37` in `script.js`; unit tests in
`tests/genomic-coordinates.test.js`, browser coverage in
`tests/browser-smoke.js`.

## Upstream reliability and the lookup progress panel

Every nomenclature lookup depends on `grch37.rest.ensembl.org`. That host
intermittently returns 500/503 and has been measured taking **20-40 seconds** to
answer a cold `variant_recoder` or `vep` query.

The variant recoder now **fails over to `rest.ensembl.org`** when the GRCh37
mirror is unwell (transport errors and retryable 5xx only — a 4xx is a real
answer and is not shopped around). The main host's coordinates are GRCh38, so
candidate conversion reads the assembly off the RefSeq chromosome accession
version (`NC_000007.13` = GRCh37, `.14` = GRCh38) and lifts **only** GRCh38
candidates back to hg19. This assembly check also fixed a latent bug: the
GRCh37 mirror's substitution candidates were previously *always* passed through
the hg38→hg19 liftover, and Ensembl's `/map` does not error when handed a
coordinate from the wrong assembly — it silently returns a position hundreds of
kb away. A GRCh38 candidate the client-side helper cannot lift (indels; the
helper handles substitutions only) is dropped rather than carried forward
mislabelled as hg19.

The failover is for **transcript-level queries only** (`recoderHostsFor`). A
chromosomal coordinate means nothing without its assembly, every genomic query
that reaches the recoder is GRCh37, and the main host reads the same numbers as
GRCh38 — another locus — without complaint: asked for `chr9:g.21971058del` it
returns SPDI `NC_000009.12:21971057:C:`, which the candidate conversion would
then "lift back" to `chr9:g.21971057C>-`, the deletion of a different base. So
a genomic query fails over to nobody, and with the mirror down the lookup
reports an outage. The mitochondrion, whose sequence is the same in both
assemblies, still fails over.

When Ensembl has not answered an accessioned name within 2.5 s — or has failed
outright — the recoder also asks **NCBI Variation Services**
(`api.ncbi.nlm.nih.gov/variation/v0`, CORS-open, independent infrastructure), and
the first usable answer is taken. It used to be asked only after both Ensembl
hosts had failed, four 12 s attempts later: `NM_007294.4:c.68_69del` (BRCA1
185delAG) waited 53 s for an answer NCBI gives in one. Ensembl still goes first
and is preferred when it is prompt, because its answer carries transcript rows.
`/hgvs/{query}/contextuals` places the variant as a SPDI, and
`/spdi/{spdi}/all_equivalent_contextual` remaps it to the GRCh37 `NC_`
chromosome placement. The scope is deliberately narrow — accessioned nucleotide
HGVS only (`NM_/NC_/NG_/NR_` with `c./g./n.`); gene-symbol and protein queries
still need Ensembl. NCBI's contextual SPDIs are VCF-anchored ("TAAT→T" for a
deletion), so they are trimmed to the minimal event in **both** trim orders
(`minimalSpdiForms`) — in a repeat region the order decides which flank
survives, and MyVariant indexes only its normalised form — and fed into the
existing candidate-conversion path as a recoder-shaped response with no
transcript nomenclature. The progress panel says "variant placed by NCBI
Variation Services — Ensembl had not answered" rather than pretending Ensembl
did. Once NCBI's answer has been taken, Ensembl's remaining attempts are called
off (`shouldStop` in `fetchWithRetry`): nothing reads their answer, and a
service slow enough to be overtaken should not be asked three more times.

The app previously made a single attempt behind a 6-7 second deadline, so a
slow or unwell upstream produced:

> Variant not found. Please verify the genomic coordinate and reference allele.

— which blames the user's input for an outage, and sends them to re-check a
coordinate that was never wrong.

Two changes address this:

- **`fetchWithRetry`** wraps `fetchWithTimeout` with bounded retries and
  exponential backoff plus jitter. Transport errors and retryable statuses
  (408, 425, 429, 500, 502, 503, 504) are retried; 4xx is not, because a 404
  from MyVariant.info is a real answer — it holds no record for many indels.
  Per-attempt deadlines were raised to match measured latency (recoder 12 s,
  VEP 20 s, MyVariant 10 s). See `tests/fetch-retry.test.js`.

- **The lookup progress panel** (`LookupProgress` in `script.js`, `.lookup-progress`
  in `style.css`) replaces the single mutable status line. It shows what the
  input resolved to — gene, genomic/coding/protein HGVS, effect, assembly — and
  one row per upstream source with its state, timing and retry count.

The panel distinguishes three outcomes that used to look identical:

| Glyph | State | Meaning |
| --- | --- | --- |
| `●` | ok | the source answered and had data |
| `◍` | empty | the source answered and holds no record — routine for indels |
| `✕` | fail | the source did not answer (status, timeout or network error) |

The final error message is derived from that record, so an Ensembl outage now
reads "…did not respond. This is an upstream outage, not necessarily a problem
with your variant — please retry in a moment" rather than blaming the
coordinate.

## TP53 backend (Vercel)

This repo now includes a serverless endpoint at `api/tp53.js` for TP53 mutation database lookups.

### Endpoint

- `POST /api/tp53`
- JSON payload:

```json
{
  "gene": "TP53",
  "protein": "p.R175H",
  "cdna": "c.524G>A",
  "genomic": "chr17:g.7673803C>T"
}
```

### Response (example)

The API returns:
- `match_count`
- `classification`
- `prevalence`
- `matches` (top matched rows)
- `dataset_url` used by the proxy
- optional `debug` object (when `debug: true` is passed in the request payload)

Protein matching normalizes HGVS three-letter amino-acid notation to single-letter
notation (for example `p.Arg273His` → `R273H`) before scoring. A row scores on an
**exact** protein token match, an exact cDNA substring, or an exact genomic position;
same-codon variants (e.g. `R175C` when you queried `R175H`) are surfaced in `matches`
for context (flagged `same_codon_match: true`) but do not count toward `match_count`
(they are reported separately as `related_codon_count`) and never drive the
pathogenicity summary.

### Dataset source links

The backend loads the dataset directly from the NCI TP53 database's CSV files under
`https://tp53.cancer.gov/static/data/` (e.g. `MutationView_r21.csv`). It does **not**
use the `view_data?bq_view_name=*` URLs, which render HTML pages rather than CSV.

To survive release bumps (r21 → r22 → …), the proxy first scrapes the current-release
CSV links from the download page `https://tp53.cancer.gov/get_tp53data`; if that can't
be reached it falls back to the hard-coded r21 URLs. If a refresh fails but a dataset
was previously loaded, the last-good copy keeps being served instead of erroring.

### Optional env var override

Set this in Vercel Project Settings → Environment Variables:

- `TP53_MUTATION_DATASET_URL` = a direct CSV URL that takes priority over discovery and the built-in list.

If not set, the function discovers the current release and falls back to the built-in URLs.

### Troubleshooting mode

Pass `debug: true` in the POST body to receive:
- normalized protein query value used for matching
- dataset fetch/download attempt log (`dataset_attempts`)
- sample protein values parsed from rows

## Optional OpenRouter AI review

This repo includes an optional AI interpretation endpoint at `api/ai-review.js`. The webpage adds a bottom "Optional AI Review" card after a variant lookup. The card does not send data automatically; users choose a model and click **Run AI review** to send the retrieved annotation, transcript, FDA companion-diagnostic, and clinical-trial context to the backend proxy.

### Endpoint

- `POST /api/ai-review`
- Requires the Vercel environment variable `OPENROUTER_API_KEY` (unless the caller supplies their own key — see BYO-key below).
- Available model choices are kept in sync between the frontend selector and backend allowlist: `openai/gpt-6-luna` (default), `openai/gpt-5-mini`, `openai/gpt-5-nano`, `openai/gpt-oss-120b`, `anthropic/claude-haiku-5.5`, `google/gemini-2.5-flash-lite`, `google/gemini-3.1-flash-lite`, `deepseek/deepseek-v4-flash`, `deepseek/deepseek-v4-pro`, `qwen/qwen3.7-flash`, `minimax/minimax-m2.7`, and `nvidia/nemotron-3-super-120b-a12b:free`. The list favors cheap models with a 1M-token window because review contexts are input-heavy (hundreds of thousands of tokens for rich genes); `gpt-oss-120b` has a 131k window and cannot take the largest payloads. The `:free` Nemotron route costs nothing but is served by third-party free-tier providers whose data-use terms differ from paid routes — check OpenRouter's provider policy before relying on it for sensitive context.
- JSON payload:

```json
{
  "model": "openai/gpt-6-luna",
  "context": {
    "submitted_variant": "BRAF V600E",
    "tumor_type": "melanoma",
    "gene": "BRAF"
  }
}
```

The endpoint prompts OpenRouter for strict JSON containing pathogenicity, specific AMP tier/subtier, FDA-approved therapies, clinical trials, a short interpretation summary, and limitations. AI output is intended for research and education only and must be verified against current FDA labeling, guidelines, curated databases, and trial eligibility criteria before clinical use.

### Abuse protection for the owner's key

The proxy protects the owner's `OPENROUTER_API_KEY` from being drained by anonymous callers. Every layer is **graceful** — it stays off until its environment variables are set, so the site keeps working before you provision anything:

- **Function duration:** `vercel.json` sets `maxDuration` 300s for `api/ai-review.js` (the Hobby-plan ceiling with Fluid Compute enabled) and the handler aborts the OpenRouter call at 290s, so slow models on large contexts return this handler's own error rather than a platform timeout.
- **Per-request caps (always on):** the completion is capped at `max_tokens` 8000, the context is clamped to 1,500,000 chars, user notes to 4,000 chars, and request bodies over ~4 MB (just under Vercel's serverless request-body ceiling) are rejected with `413`. These bound the cost of any single call. The context clamp is generous because rich genes assemble very large payloads (e.g. HER2/ERBB2 amplification ≈ 957 KB, almost all openFDA label text). To keep such payloads under the ceiling, the frontend caps the openFDA **record count** (first 40) for the AI query while keeping each record's **full** `indications_and_usage` text; when records are dropped, an `openfda_records_truncated` note is added to the context so the model knows the list was capped. Large contexts may still exceed some models' windows — prefer a large-window model.
- **Origin allowlist:** browser requests are restricted to the live site plus test hosts. Defaults: `https://drdoubleb.com`, `https://www.drdoubleb.com`, `https://variant-search-for-chat-gpt.vercel.app`, and this project's own previews (`variant-search-for-chat-gpt-*.vercel.app`) — a bare `*.vercel.app` was previously allowed, which admitted *anyone's* Vercel deployment. Override with `AI_ALLOWED_ORIGINS` (comma-separated; entries may embed `*` wildcards that match within one DNS label; a bare `*` disables the check). Requests with no `Origin` header (curl, server-to-server) are allowed through to the other layers.
- **Cloudflare Turnstile:** set `TURNSTILE_SECRET_KEY` to require a bot-challenge token (`turnstile_token` in the body) on owner-key requests. Unset → skipped.
- **Rate limiting (Upstash Redis):** set `UPSTASH_REDIS_REST_URL` and `UPSTASH_REDIS_REST_TOKEN` to enable per-IP and a global daily circuit breaker. Defaults: `AI_RL_IP_PER_MIN` 12, `AI_RL_IP_PER_DAY` 120, `AI_RL_GLOBAL_PER_DAY` 1500. On limit, returns `429` (per-IP) or `503` (global daily) with a `Retry-After` header. Unset → skipped.

### Why AI reviews fail, and the token budget

`max_tokens` is shared with the model's **reasoning** tokens on reasoning models (`gpt-*`, `deepseek-v4-*`, `qwen3.7-flash`). A model that thinks for 1,500 tokens has only the remainder left for the answer, so a low cap makes the JSON stop mid-array — or come back empty when reasoning consumes the whole budget. That is why the cap is 8000 rather than 3000, and why the same variant can succeed on one model and fail on another.

The proxy now separates the causes instead of surfacing a raw parser message like `Expected ',' or ']' after array element in JSON at position 4348`:

| Condition | Message |
| --- | --- |
| `finish_reason: 'length'` | ran out of output tokens — try another model or a smaller context |
| empty content | spent its whole budget on internal reasoning |
| no `{` in content | replied with text instead of JSON |
| unbalanced `{...}` | JSON response was incomplete or malformed |

Both prompts also **restate their output format after the context payload**. The payload can be hundreds of thousands of tokens, so a format instruction placed before it is far from the end; restating it last keeps weak models from drifting into prose. `tests/ai-review-response.test.js` covers the failure branches and `tests/ai-review-prompt.test.js` pins the ordering.

### Bring-your-own key (BYO-key)

Callers may supply their own OpenRouter key as `openrouter_api_key` in the POST body (format `sk-...`). BYO-key requests pay their own cost, so they **bypass Turnstile and the rate limiter** (a useful valve when the shared key or a user's own balance runs out). A malformed key returns `400`. The success response includes `used_user_key: true|false`.

### Copy-the-prompt mode

Send `{ "mode": "prompt", "context": { ... } }` to receive the fully assembled prompt (`{ prompt, model }`) without calling any model — this powers the "copy prompt" button so users can paste it into their own LLM. It requires no key and incurs no cost (origin-checked only).

Because a person pasting the prompt into a general chat LLM wants the readable answer (not the strict JSON the site parses into cards), this mode returns a **separate, self-contained prompt** that asks for a Markdown-formatted, human-readable interpretation covering the same fields and tiering rules. It does not mention JSON at all.

#### Prompt file layout

The two prompts share their clinical guidance so they cannot drift apart:

| File | Role |
| --- | --- |
| `api/_ai-review-prompt-core.js` | **Shared clinical guidance** — AMP/ASCO/CAP tier definitions, hard tiering rules, controlled vocabularies, interpretation rules, and the self-check. Format-neutral: it never mentions JSON, schemas, or arrays. **Tiering and evidence edits go here** and take effect in both prompts. |
| `api/_ai-review-prompt.js` | Core + a **JSON** output-format section (the schema, field semantics). Used by "Run AI review". |
| `api/_ai-review-prompt-human.js` | Core + a **Markdown** output-format section (numbered sections for a reader). Used by "Copy prompt". |

All three are underscore-prefixed because Vercel turns every non-underscore `.js` file in `api/` into a Serverless Function, and the Hobby plan caps a deployment at **12**. These are imported modules, not endpoints — same convention as `_ncbi.js`, `_ratelimit.js`, and `_turnstile.js`. `tests/vercel-function-count.test.js` enforces the budget and checks that every counted file really does default-export a handler, so a missing underscore fails locally instead of at deploy time.

Earlier versions used one JSON template plus a trailing "disregard the JSON instructions" override for copy-prompt mode. The resulting self-contradiction confused weaker models, so each prompt is now internally consistent end to end. `tests/ai-review-prompt.test.js` enforces both properties: that the clinical core is byte-identical in the two prompts, and that the human prompt contains no JSON/schema instruction or override phrasing.

### Optional dependencies

Rate limiting uses [`@upstash/ratelimit`](https://github.com/upstash/ratelimit-js) and `@upstash/redis` (declared in `package.json`). They are imported lazily and only when the `UPSTASH_*` env vars are present, so deployments that don't enable rate limiting don't need them at runtime.

The AI review payload also includes a `supplemental_card_data` object populated from live card lookups when available, including direct ClinVar VCV/nearby-variant results, CIViC API assertions, gnomAD v4, SpliceAI Lookup scores, PubMed article previews, COSMIC extended data, and the TP53 mutation database for TP53 variants.


## What a protein change is searched as

The ClinVar search link, the PubMed/LitVar card, the Google links and the AI
review all look the variant up under a short name for its protein change.
`proteinSearchTerm()` is the one place that name is made:

| Protein change | Searched as | |
| --- | --- | --- |
| `p.Val600Glu` | `V600E` | |
| `p.Gln1756ProfsTer74`, `p.Pro153fs` | `Q1756fs`, `P153fs` | ClinVar's own protein-change name for a frameshift, and the form LitVar resolves |
| `p.Arg213Ter`, `p.Arg213*` | `R213X` | one ClinVar record, the LitVar entity that holds the literature, and a term PubMed knows |
| in-frame indels, extensions, synonymous | the full `p.` string | |

A frameshift opens the way a substitution does, and the old shortening read only
that far: `p.Gln1756ProfsTer74` (BRCA1 c.5266dup, 5382insC) was searched as the
missense **Q1756P** — a different variant with its own ClinVar records and its
own paper — and a frameshift written without the new residue, `p.Pro153fs`,
came out as "TP53 **O153F**". The spellings in the table were chosen by asking
the services: `BRCA1 Q1756Pfs*74` finds nothing in LitVar and, the asterisk
being a wildcard, 85,659 records in ClinVar; `TP53 p.Arg213*` (snpEff's
spelling of a stop) is 6,089 ClinVar records and no LitVar entity.

Only a substitution — missense or nonsense — gets the ClinVar card's "Same
protein change" list. Two frameshifts at one residue are different proteins
unless they shift the frame the same way, and ClinVar's titles (`p.Gln1756fs`)
do not say which way.

**PubMed drops what it does not know.** Asked for `BRCA2 D362fs`, PubMed finds
no article using "D362fs", silently discards the phrase and answers with the
17,385 papers on BRCA2 — which were listed under "Gene + Variant" and handed to
the AI review as literature on the variant. `esearch` names what it discarded
(`errorlist.phrasesnotfound`); when that is the whole variant, `api/pubmed.js`
returns no articles and the phrase, and the card says
`no article in PubMed uses "D362fs"`. Part of a variant surviving
(`p.Glu746_Ala750del` loses `p.Glu746` and still matches on `Ala750del`) is left
alone, as is any search that names no variant.

**A stop typed by hand.** `TP53 R213*`, `R213X` and `p.Arg213*` are all read as
`TP53:p.Arg213Ter`; the first two used to be passed on as typed and came back
"could not resolve".

See `tests/protein-search-term.test.js`.

## Variant-normalized literature (LitVar2)

A free-text PubMed search for "BRAF V600E" misses papers that write the same
variant as `p.Val600Glu`, `c.1799T>A`, or `rs113488022`. When the lookup has a
specific variant, the PubMed card's variant tabs go through NCBI's **LitVar2**
(PubTator3 API) instead: the query is resolved to a variant entity, and the
entity search returns relevance-ranked articles across every nomenclature
spelling — optionally refined with the tumor type. The card shows a provenance
line naming the matched entity/rsID and deep-links to LitVar2; when no entity
matches (or the service is down) the plain term search runs unchanged. This
lives inside `api/pubmed.js` (`gene`/`variant`/`extra` params), so no extra
serverless function is spent.

## Cancer Prevalence card (cBioPortal)

Answers "how often is this variant seen, and in which tumor types?" from
cBioPortal's public API, which is CORS-open and therefore called directly from
the browser (no serverless function). v1 queries the **MSK-IMPACT 2017** cohort
(10,945 prospectively sequenced advanced-cancer tumors): gene mutation
frequency, the exact protein change's frequency (single-letter-normalized
match), the gene's most frequent changes, and a tumor-type breakdown for the
matched variant. The card labels the cohort explicitly — these are
single-institution advanced-cancer cohort frequencies, not population
prevalence. Results are included in the AI-review payload under
`supplemental_card_data.cbioportal_prevalence`.

## ClinGen expert-panel classifications (eRepo)

When a variant has a ClinGen VCEP (expert panel) classification, the ClinVar
card shows it — classification, panel, publication date, and a link to the
Evidence Repository record. VCEP assertions apply the ACMG criteria via a
gene-specific expert process, so they are higher-confidence than aggregate
ClinVar where they exist (most variants have none; the line only renders on a
hit). The lookup is a single CORS-open GET against
`erepo.clinicalgenome.org` using the GRCh37 `NC_` genomic HGVS built from the
resolved coordinate, so it matches regardless of which transcript version the
panel used. Included in the AI payload as
`supplemental_card_data.clingen_erepo`.

## gnomAD gene constraint

The gnomAD card now ends with a gene-constraint line — pLI, LOEUF
(`oe_lof_upper`), missense Z, and o/e missense from gnomAD v2 — fetched with
one CORS-open GraphQL query. Constraint says how tolerant the gene is to
loss-of-function and missense variation, which is crucial context for
truncating variants. Included in the AI payload as
`supplemental_card_data.gnomad_constraint`.

## DGIdb drug–gene interactions

The FDA drugs card gains a fourth tab, **DGIdb**: aggregated drug–gene
interaction claims (approved *and* investigational compounds) from DGIdb's
CORS-open GraphQL API, sorted by interaction score with approval status,
interaction type, and source counts. Clearly labelled as interaction claims —
not efficacy or approval evidence. Deliberately **excluded from the AI-review
payload**: DGIdb's aggregated investigational/preclinical claims would dilute
the clinical sources the review is meant to weigh — the tab is for human
browsing.

## Proxy origin allowlist

Every serverless proxy (not just `ai-review`) now checks the browser `Origin`
header against a shared allowlist (`api/_origin.js`): the live site, the
production Vercel host, and this project's own `variant-search-for-chat-gpt-*`
previews. Before this, any third-party website's JavaScript could use the
deployment — and its `NCBI_API_KEY` quota and Vercel invocation budget — as its
backend. Requests with **no** `Origin` header (curl, server-to-server
integrations, GPT actions) are deliberately allowed through; only third-party
browser embedding is blocked. Override with `PROXY_ALLOWED_ORIGINS`
(comma-separated; falls back to `AI_ALLOWED_ORIGINS`; `*` disables the check).

## SpliceAI Lookup proxy

The SpliceAI card uses `api/spliceai.js` as a lightweight proxy to the Broad Institute SpliceAI Lookup API. The proxy accepts variants in `chr-pos-ref-alt` format and forwards interactive-use requests with hg, distance, mask, and Gencode (`bc`) parameters. Returned scores are displayed in the SpliceAI card and included in the AI review payload under `supplemental_card_data.spliceai_lookup`.

## Genomic coordinates and VCF alleles

Several cards need more than a gene and a protein change: the UCSC link and the
ClinVar region pull need a genomic position, while gnomAD, SpliceAI and the
gnomAD variant pages need VCF-style `chr-pos-ref-alt` alleles.

Only substitutions carry their alleles in the notation. A multi-base indel comes
back from the Ensembl variant recoder as a bare span — `TSC2 c.2319_2321delAAT`
normalises to `chr16:g.2122948_2122950del`, with no REF or ALT anywhere in the
string — and MyVariant.info does not index every such variant, so its `vcf` block
is often absent too. Coordinates are therefore resolved in two stages:

1. **`parseGenomicHgvs` / `buildVariantCoordinateTuple`** parse every genomic HGVS
   form (`sub`, `del`, `dup`, `ins`, `inv`, `delins`, and the `REF>-` / `->ALT`
   spellings the SPDI converter emits) into a chromosome plus a start/end span.
   Alleles are reported as `null` rather than guessed. Position-only consumers —
   the UCSC hg19 link, the ClinVar region search and nearby-variant plot, the
   gnomAD region view — work off this alone.
2. **`resolveVcfAlleles`** fills in the alleles for allele-dependent consumers,
   preferring what the notation already carries, then MyVariant's `vcf` block, and
   finally reading the bases off the GRCh37 reference (Ensembl
   `/sequence/region`, ±60 bp) and anchoring them VCF-style on the preceding base.
   The result is then **left-aligned**, because HGVS shifts indels 3′ while VCF —
   and therefore gnomAD — indexes the 5′-most representation. CFTR F508del is
   `g.117199646_117199648del` in HGVS but `7-117199644-ATCT-A` in gnomAD; without
   left-alignment the lookup misses. Events wider than 1 kb are skipped, as they
   are not represented as ref/alt strings by those resources anyway.

The resolution runs once per lookup and is shared by the Variant card's
**VCF (hg19)** line, the gnomAD v2 link, the gnomAD v4.1 query, the SpliceAI card
and the AI review context. Cards render immediately and upgrade in place when the
alleles land.

The Variant card shows the resolved `CHROM-POS-REF-ALT` next to the g. notation,
since that is the form gnomAD, SpliceAI and most VCF-derived tools index by. When
left-alignment moved the position it is marked `(left-aligned)` with a tooltip, so
a coordinate that disagrees with the 3'-shifted g. notation above it does not read
as a bug.

`api/gnomad-v4.js` lifts only the anchor position GRCh37→GRCh38 and reuses the
alleles, which is exact for SNVs but can land an indel a few bases off how gnomAD
indexes it. When the constructed ID misses for a non-SNV, the proxy scans a small
region around the lifted position for a record carrying the same alleles before
reporting `not_found`; a match is reported with a `matchedVia` note.

Note that gnomAD reports an absent variant as a GraphQL *error* (`Variant not
found`) rather than a null result, so that particular message is treated as
absence rather than as an API failure — otherwise the region scan above is never
reached and the card shows an error where "not found" belongs.

### Liftover (`api/_liftover.js`)

The GRCh38 coordinate is resolved in this order:

1. **A caller-supplied `pos38`.** myvariant.info returns a GRCh38 start in
   `dbnsfp.hg38` / `dbsnp.hg38` for most indexed variants, and the frontend
   already holds that annotation, so it passes it to the proxy. This skips the
   network round-trip entirely and is the common path.
2. **Ensembl `/map`**, across two hosts — `rest.ensembl.org` and the
   `grch37.rest.ensembl.org` mirror, which serves the same paths from a separate
   deployment. Requests are *hedged*: the mirror starts once the primary has
   failed or gone quiet for ~1.2s, and the first usable answer wins. A healthy
   primary answers well inside that window, so the mirror is normally never hit.
   The whole ladder works to a 6s deadline, below the 15s function budget it
   shares with the gnomAD query.

This replaced a single un-retried call to `rest.ensembl.org`. That host is a real
single point of failure: during one outage it answered HTTP 500 or hung on every
`/map` and `/vep` request for hours while `/info/ping` stayed green, and the card
reported "liftover failed for 4:143094904" — a coordinate that maps fine, to
4:142173751.

Outcomes are also reported distinctly. `liftover_unavailable` means the provider
could not be reached and retrying later will probably work; `liftover_failed`
means Ensembl answered and the coordinate has no unambiguous GRCh38 equivalent.
Collapsing the two blames the variant for someone else's downtime. Mappings are
validated against the requested assembly and chromosome, and several disagreeing
candidates count as no answer — the previous code took `mappings[0]` blind, which
turns a patch/scaffold hit into a confident wrong coordinate.

The client-side hg38→hg19 helper in `script.js` fails over across the same two
hosts. It returns its input unchanged when conversion fails, so an outage there
silently mislabels an hg38 coordinate as hg19 rather than erroring.

ClinVar recovery from the region pull matches SNVs on position plus alleles (with
strand complementation). Indels have no comparable alleles and their coordinates
differ between HGVS and ClinVar's own representation, so they are matched on the
c. notation in the record title instead, normalised so that the user's
`c.2319_2321delAAT` and ClinVar's `c.2319_2321del` compare equal.

## Resolving c. and p. queries to a genomic position

The Ensembl variant recoder is the primary route from a non-genomic query to a
genomic coordinate, but it does not cover everything, and the fallbacks matter.

**cDNA (c.) queries.** When the recoder is unavailable or returns no usable
candidate, the app falls back to Ensembl VEP, which resolves gene-level cDNA HGVS
(`TSC2:c.4952delA`) that the recoder can miss. VEP's response carries the genomic
location, so `buildGenomicHgvsFromVep` turns it into a g. string and the lookup
keeps a real coordinate. Previously only the consequence was read off that
response, so `g.` kept echoing the user's own query and every coordinate-derived
card stayed dark even though the variant had been resolved.

**Protein (p.) queries.** These are fundamentally harder, and for two shapes they
are impossible:

- Frameshifts are rejected outright — *"Frameshifts are not supported for HGVS
  protein input"*. `p.N1651Mfs*21` can arise from many different indels.
- A delins spanning several residues has no unique coding change either
  (`p.L773_I774delinsF` → *"Could not determine nucleotide change from peptide
  change"*), in single- or three-letter form, against a gene or a protein
  accession.

Simple substitutions (`p.Val600Glu`) do resolve, so those are unaffected.

What *is* determined by protein notation is the codon. Rather than failing the
whole lookup, `resolveProteinCodonRegion` maps the residue range to GRCh37
coordinates via the gene's canonical translation (Ensembl `/lookup/symbol` then
`/map/translation`), and the lookup continues with `gVariant` set to a bare span
such as `chr16:g.2136834_2136836`. `parseGenomicHgvs` reads that as type `region`:
position-based cards (UCSC, ClinVar region search, nearby variants, gene-level
CIViC/OncoKB/FDA/trials/PubMed) work, while allele-based ones (VCF line, gnomAD
v4.1, SpliceAI) correctly report nothing. The Variant card labels the line
**Codon region (hg19)** and states which residues resolved and why the exact
nucleotide change is unavailable. A residue beyond the end of the canonical
protein is rejected rather than silently mapped to the wrong locus.

Queries that cannot even be resolved to a codon fail with a message naming the
reason and pointing at c./g. notation, instead of the previous bare
"Variant not found via Ensembl Variant Recoder".

Before settling for the codon, the name is looked up among ClinVar's other names
(next section): `BRCA1 p.Gln1756ProfsTer74` is one ClinVar variant and nothing
else, so it resolves to that variant, not to a region. (The region it used to
get was wrong as well: "Ter74", the distance to the new stop, was read as a
residue, giving "residues 74–1756" — 47 kb of BRCA1.)

## Names that are not HGVS

"BRCA1 5382insC", "BRCA1 185delAG", "BRCA2 6174delT", "MTHFR C677T": what these
variants were called before HGVS, counted from a different first base, and still
what a clinic letter says. No nomenclature service resolves them.

| Input | before | now |
| --- | --- | --- |
| `BRCA1 5382insC` | sent on as `BRCA1:p.5382insC`; "could not resolve" | NM_007294.4:c.5266dup · ClinVar 17677 |
| `BRCA1 c.5382insC` | "Lookup could not complete … please retry" | the same |
| `BRCA1 185delAG` | read as a protein change, `p.185delAG` | NM_007294.4:c.68_69del · ClinVar 17662 |
| `BRCA2 6174delT` | `p.6174delT` | NM_000059.4:c.5946del · ClinVar 9325 |
| `MTHFR C677T` | read as p.Cys677Thr; "could not resolve" | NM_005957.5:c.665C>T (p.Ala222Val) · ClinVar 3520 |
| `BRCA1 p.Gln1756ProfsTer74`, `BRCA1 Q1756Pfs*74` | a 47 kb "codon region" | NM_007294.4:c.5266dup · ClinVar 17677 |

ClinVar lists them: every record carries the other names its variant goes by
(`aliases` — "5382insC", "5385insC", "p.Gln1756ProfsX74", and OMIM-style strings
such as "MTHFR, 677C-T, ALA222VAL (rs1801133)"). When a `GENE name` input will
not resolve, `resolveClinvarOtherName` asks the existing `/api/clinvar-protein`
proxy for `GENE name` (it is a plain ClinVar search; it now returns each
record's `aliases` and `variationCount`) and accepts a record only if

- it is a single variant (not a haplotype) named in that gene, and
- one of its other names **is** the typed name, compared in one spelling
  (`variantNameKeys`: "c." and the sequence dropped, a protein change in
  single-letter codes with one mark for a stop, and `C677T` = `677C>T` =
  `677C-T`), and
- it is the only such record. Two variants under one name are not for the page
  to choose between.

The lookup then runs again for that record's HGVS name, and says so in the
progress panel and on the Variant card: *"5382insC" is another name ClinVar
lists for NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs) (Variation 17677), which
is what is shown here.* What was typed stays in "You entered" and in the URL.

This is a last resort, asked once per lookup and only after the ordinary route
has failed — a name that resolves never consults it — and not when the failure
was a service being down. A token that starts with a digit is no longer given a
`p.`: `BRCA1 9999delT` is refused at once as *not a name this tool can read …
ClinVar lists no BRCA1 variant under that name either*, without a request to
Ensembl or MyVariant. If the proxy cannot be reached the panel says "could not
be consulted" and no claim is made about what ClinVar lists; a backend deployed
before this change sends no `aliases`, so nothing resolves and nothing is
wrongly resolved.

See `tests/other-names.test.js`.

## Test suite imports the real code

`npm test` runs every `tests/*.test.js` with plain Node — no framework. The
tests used to re-declare the helpers they exercised ("copied verbatim — KEEP
IN SYNC" blocks), which drifted from production by design pressure alone.
They now import the real implementations:

- **`script.js`** stays a single classic `<script>` for the browser (no build
  step, one file to copy to Hostinger — `export` statements would be a
  SyntaxError there), but Node imports it fine as a module. A guarded block at
  the bottom stashes the pure helpers on `globalThis.__variantSearchHelpers`
  when running under Node; tests `await import('../script.js')` and destructure
  from there. The `DOMContentLoaded` registration is skipped when `document`
  does not exist, and `fetchWithRetry` is built through a `makeFetchWithRetry`
  factory so the retry loop can be driven with an injected fetch.
- **`api/*.js`** are ES modules already; the handlers now also carry named
  exports for their pure helpers (`export { buildAssertions, … }`), which the
  Vercel runtime ignores. `api/tp53.js` gained `extractDatasetUrls` — the pure
  parsing core of `discoverDatasetUrls` — as a real seam instead of a test-only
  re-implementation.

**Whole lookups are replayed, not mocked by hand.** `tests/support/lookup-harness.js`
runs the real form handler (LinkeDOM + `vm`, no browser) against a *cassette*:
every answer MyVariant and Ensembl gave when the scenario was run live, keyed by
request URL, plus the stretch of GRCh37 reference it read
(`tests/fixtures/lookup-cassettes.json`). A request the recording does not hold
gets what that service says when it has nothing (404, an empty hit list, 503),
and a test can take a service away from a complete recording
(`without: [/ensembl\.org/]`) or answer one of the site's own `/api` proxies
(`backend: { pubmed: () => … }`). Each scenario was confirmed to fail on the
commit before its fix.

The one remaining intentional re-composition is `buildChangeForms` in
`tests/clinvar-match.test.js`, which mirrors logic composed inline in the
ClinVar card rather than a named production function.
