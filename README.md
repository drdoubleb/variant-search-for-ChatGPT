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
