# Variant Search (Beta)

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

- **Per-request caps (always on):** the completion is capped at `max_tokens` 3000, the context is clamped to 1,500,000 chars, user notes to 4,000 chars, and request bodies over ~4 MB (just under Vercel's serverless request-body ceiling) are rejected with `413`. These bound the cost of any single call. The context clamp is generous because rich genes assemble very large payloads (e.g. HER2/ERBB2 amplification ≈ 957 KB, almost all openFDA label text). To keep such payloads under the ceiling, the frontend caps the openFDA **record count** (first 40) for the AI query while keeping each record's **full** `indications_and_usage` text; when records are dropped, an `openfda_records_truncated` note is added to the context so the model knows the list was capped. Large contexts may still exceed some models' windows — prefer a large-window model.
- **Origin allowlist:** browser requests are restricted to the live site plus test hosts. Defaults: `https://drdoubleb.com`, `https://www.drdoubleb.com`, `https://variant-search-for-chat-gpt.vercel.app`, and `*.vercel.app` previews. Override with `AI_ALLOWED_ORIGINS` (comma-separated; `*` disables the check). Requests with no `Origin` header (curl, server-to-server) are allowed through to the other layers.
- **Cloudflare Turnstile:** set `TURNSTILE_SECRET_KEY` to require a bot-challenge token (`turnstile_token` in the body) on owner-key requests. Unset → skipped.
- **Rate limiting (Upstash Redis):** set `UPSTASH_REDIS_REST_URL` and `UPSTASH_REDIS_REST_TOKEN` to enable per-IP and a global daily circuit breaker. Defaults: `AI_RL_IP_PER_MIN` 12, `AI_RL_IP_PER_DAY` 120, `AI_RL_GLOBAL_PER_DAY` 1500. On limit, returns `429` (per-IP) or `503` (global daily) with a `Retry-After` header. Unset → skipped.

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
