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
notation (for example `p.Arg273His` → `R273H`) before scoring, improving matching
against TP53 dataset rows that use single-letter protein codes.

### Dataset source links

The backend now hard-codes TP53 dataset endpoints from `tp53.cancer.gov`:

- `https://tp53.cancer.gov/view_data?bq_view_name=MutationView`
- `https://tp53.cancer.gov/view_data?bq_view_name=MutationViewDownload`
- `https://tp53.cancer.gov/view_data?bq_view_name=TumorVariantDownload`
- `https://tp53.cancer.gov/view_data?bq_view_name=GermlineDownload`

### Optional env var override

Set this in Vercel Project Settings → Environment Variables:

- `TP53_MUTATION_DATASET_URL` = custom direct URL to override the built-in hard-coded list.

If not set, the function uses the built-in hard-coded TP53 URLs above.

### Troubleshooting mode

Pass `debug: true` in the POST body to receive:
- normalized protein query value used for matching
- dataset fetch/download attempt log (`dataset_attempts`)
- sample protein values parsed from rows

## Optional OpenRouter AI review

This repo includes an optional AI interpretation endpoint at `api/ai-review.js`. The webpage adds a bottom "Optional AI Review" card after a variant lookup. The card does not send data automatically; users choose a model and click **Run AI review** to send the retrieved annotation, transcript, FDA companion-diagnostic, and clinical-trial context to the backend proxy.

### Endpoint

- `POST /api/ai-review`
- Requires the Vercel environment variable `OPENROUTER_API_KEY`.
- Available model choices are kept in sync between the frontend selector and backend allowlist: `openai/gpt-4.1-mini`, `openai/gpt-oss-120b`, `openai/gpt-4o-mini`, `openai/gpt-5-mini`, `openai/gpt-5.4-nano`, `openai/gpt-5-nano`, `google/gemini-2.5-flash-lite`, `google/gemini-3-flash-preview`, `google/gemini-3.1-flash-lite`, `anthropic/claude-3-haiku`, `deepseek/deepseek-v4-flash`, `deepseek/deepseek-v4-pro`, `x-ai/grok-4.3`, `minimax/minimax-m2.7`, and `nvidia/nemotron-3-super-120b-a12b:free`.
- JSON payload:

```json
{
  "model": "openai/gpt-4.1-mini",
  "context": {
    "submitted_variant": "BRAF V600E",
    "tumor_type": "melanoma",
    "gene": "BRAF"
  }
}
```

The endpoint prompts OpenRouter for strict JSON containing pathogenicity, AMP tier, FDA-approved therapies, clinical trials, a short interpretation summary, and limitations. AI output is intended for research and education only and must be verified against current FDA labeling, guidelines, curated databases, and trial eligibility criteria before clinical use.

The AI review payload also includes a `supplemental_card_data` object populated from live card lookups when available, including direct ClinVar VCV/nearby-variant results, CIViC API assertions, gnomAD v4, SpliceAI Lookup scores, PubMed article previews, COSMIC extended data, and the TP53 mutation database for TP53 variants.


## SpliceAI Lookup proxy

The SpliceAI card uses `api/spliceai.js` as a lightweight proxy to the Broad Institute SpliceAI Lookup API. The proxy accepts variants in `chr-pos-ref-alt` format and forwards interactive-use requests with hg, distance, mask, and Gencode (`bc`) parameters. Returned scores are displayed in the SpliceAI card and included in the AI review payload under `supplemental_card_data.spliceai_lookup`.
