# Variant Search (Beta)

## Accepted input formats

Coordinates are interpreted as **GRCh37/hg19** — that is the assembly every
resolution endpoint in this app targets (`grch37.rest.ensembl.org`), so an hg38
position will either miss or silently resolve to the wrong locus.

HGVS is accepted in genomic, coding and protein form:

```
chrX:g.20148675dup      NM_001412.4:c.388dup      EIF1AX:p.Gln130ProfsTer3
```

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

## Upstream reliability and the lookup progress panel

Every nomenclature lookup depends on `grch37.rest.ensembl.org`. That host
intermittently returns 500/503 and has been measured taking **20-40 seconds** to
answer a cold `variant_recoder` or `vep` query. The app previously made a single
attempt behind a 6-7 second deadline, so a slow or unwell upstream produced:

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
- Available model choices are kept in sync between the frontend selector and backend allowlist: `openai/gpt-5-mini`, `openai/gpt-5.4-nano`, `openai/gpt-oss-120b`, `openai/gpt-4o-mini`, `openai/gpt-5-nano`, `google/gemini-2.5-flash-lite`, `google/gemini-3-flash-preview`, `google/gemini-3.1-flash-lite`, `anthropic/claude-3-haiku`, `deepseek/deepseek-v4-flash`, `deepseek/deepseek-v4-pro`, `x-ai/grok-4.3`, `minimax/minimax-m2.7`, `nvidia/nemotron-3-super-120b-a12b:free`, and `qwen/qwen3-32b`.
- JSON payload:

```json
{
  "model": "openai/gpt-5-mini",
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

- **Per-request caps (always on):** the completion is capped at `max_tokens` 8000, the context is clamped to 1,500,000 chars, user notes to 4,000 chars, and request bodies over ~4 MB (just under Vercel's serverless request-body ceiling) are rejected with `413`. These bound the cost of any single call. The context clamp is generous because rich genes assemble very large payloads (e.g. HER2/ERBB2 amplification ≈ 957 KB, almost all openFDA label text). To keep such payloads under the ceiling, the frontend caps the openFDA **record count** (first 40) for the AI query while keeping each record's **full** `indications_and_usage` text; when records are dropped, an `openfda_records_truncated` note is added to the context so the model knows the list was capped. Large contexts may still exceed some models' windows — prefer a large-window model.
- **Origin allowlist:** browser requests are restricted to the live site plus test hosts. Defaults: `https://drdoubleb.com`, `https://www.drdoubleb.com`, `https://variant-search-for-chat-gpt.vercel.app`, and this project's own previews (`variant-search-for-chat-gpt-*.vercel.app`) — a bare `*.vercel.app` was previously allowed, which admitted *anyone's* Vercel deployment. Override with `AI_ALLOWED_ORIGINS` (comma-separated; entries may embed `*` wildcards that match within one DNS label; a bare `*` disables the check). Requests with no `Origin` header (curl, server-to-server) are allowed through to the other layers.
- **Cloudflare Turnstile:** set `TURNSTILE_SECRET_KEY` to require a bot-challenge token (`turnstile_token` in the body) on owner-key requests. Unset → skipped.
- **Rate limiting (Upstash Redis):** set `UPSTASH_REDIS_REST_URL` and `UPSTASH_REDIS_REST_TOKEN` to enable per-IP and a global daily circuit breaker. Defaults: `AI_RL_IP_PER_MIN` 12, `AI_RL_IP_PER_DAY` 120, `AI_RL_GLOBAL_PER_DAY` 1500. On limit, returns `429` (per-IP) or `503` (global daily) with a `Retry-After` header. Unset → skipped.

### Why AI reviews fail, and the token budget

`max_tokens` is shared with the model's **reasoning** tokens on reasoning models (`gpt-5-*`, `deepseek-v4-*`, `grok-*`). A model that thinks for 1,500 tokens has only the remainder left for the answer, so a low cap makes the JSON stop mid-array — or come back empty when reasoning consumes the whole budget. That is why the cap is 8000 rather than 3000, and why the same variant can succeed on one model and fail on another.

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
