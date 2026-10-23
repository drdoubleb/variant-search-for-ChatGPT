# Backend dependency and input security

This project has a static frontend served on Hostinger and Node.js serverless
functions deployed separately on Vercel.

## Reproducible dependencies

- `package-lock.json` is committed. Keep it in sync with `package.json`.
- Install with `npm ci` locally or in CI. To update a dependency intentionally,
  run `npm install <package>@<version>` and commit both manifests.
- `npm test` runs the existing Node tests plus the security regression cases.
- GitHub Actions checks each PR and push to `main`, and re-runs a production
  dependency audit weekly. The audit fails for high/critical advisories and
  reports lower-severity ones without blocking the build.

## Serverless input checks

- The gnomAD proxy now accepts only GET and validates human chromosomes,
  positive genomic coordinates, optional GRCh38 hints and DNA alleles, returning
  HTTP 400 for malformed values.
- ClinVar gene, protein and region queries and PubMed searches have bounded
  input lengths and reject malformed coordinate, limit and repeated-query values.
- The AI-review endpoint checks the actual parsed JSON body size, even when the
  Content-Length header is missing or inaccurate.
- Browser Origin matching now requires the correct HTTPS scheme for configured
  production hosts and project Vercel previews. Requests without an Origin
  header remain allowed by design; Origin checks are **not authentication**.

## What to deploy

- This change set touches **only Vercel API handlers and developer tooling**;
  `index.html`, `script.js` and `style.css` are unchanged.
- Merging to the configured production branch should deploy the API through
  Vercel. Verify the build, a few normal lookups, and AI review behavior.
- **No Hostinger upload is necessary for this change set.**
- Production Node.js is configured to 24.x in Vercel. Keep it patched and
  maintain OpenRouter spending limits separately from the Redis rate limiter.

The security checks are preventive and cannot guarantee that every dependency,
upstream service or endpoint is free of vulnerabilities.
