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

- **Per-request caps (always on):** the completion is capped at `max_tokens` 3000, the context is clamped to 200,000 chars, user notes to 4,000 chars, and request bodies over 512 KB are rejected with `413`. These bound the cost of any single call.
- **Origin allowlist:** browser requests are restricted to the live site plus test hosts. Defaults: `https://drdoubleb.com`, `https://www.drdoubleb.com`, `https://variant-search-for-chat-gpt.vercel.app`, and `*.vercel.app` previews. Override with `AI_ALLOWED_ORIGINS` (comma-separated; `*` disables the check). Requests with no `Origin` header (curl, server-to-server) are allowed through to the other layers.
- **Cloudflare Turnstile:** set `TURNSTILE_SECRET_KEY` to require a bot-challenge token (`turnstile_token` in the body) on owner-key requests. Unset → skipped.
- **Rate limiting (Upstash Redis):** set `UPSTASH_REDIS_REST_URL` and `UPSTASH_REDIS_REST_TOKEN` to enable per-IP and a global daily circuit breaker. Defaults: `AI_RL_IP_PER_MIN` 12, `AI_RL_IP_PER_DAY` 120, `AI_RL_GLOBAL_PER_DAY` 1500. On limit, returns `429` (per-IP) or `503` (global daily) with a `Retry-After` header. Unset → skipped.

### Bring-your-own key (BYO-key)

Callers may supply their own OpenRouter key as `openrouter_api_key` in the POST body (format `sk-...`). BYO-key requests pay their own cost, so they **bypass Turnstile and the rate limiter** (a useful valve when the shared key or a user's own balance runs out). A malformed key returns `400`. The success response includes `used_user_key: true|false`.

### Copy-the-prompt mode

Send `{ "mode": "prompt", "context": { ... } }` to receive the fully assembled prompt (`{ prompt, model }`) without calling any model — this powers the "copy prompt" button so users can paste it into their own LLM. It requires no key and incurs no cost (origin-checked only).

### Optional dependencies

Rate limiting uses [`@upstash/ratelimit`](https://github.com/upstash/ratelimit-js) and `@upstash/redis` (declared in `package.json`). They are imported lazily and only when the `UPSTASH_*` env vars are present, so deployments that don't enable rate limiting don't need them at runtime.

The AI review payload also includes a `supplemental_card_data` object populated from live card lookups when available, including direct ClinVar VCV/nearby-variant results, CIViC API assertions, gnomAD v4, SpliceAI Lookup scores, PubMed article previews, COSMIC extended data, and the TP53 mutation database for TP53 variants.


## SpliceAI Lookup proxy

The SpliceAI card uses `api/spliceai.js` as a lightweight proxy to the Broad Institute SpliceAI Lookup API. The proxy accepts variants in `chr-pos-ref-alt` format and forwards interactive-use requests with hg, distance, mask, and Gencode (`bc`) parameters. Returned scores are displayed in the SpliceAI card and included in the AI review payload under `supplemental_card_data.spliceai_lookup`.
