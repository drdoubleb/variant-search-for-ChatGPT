// Serverless proxy for the gnomAD v4 GraphQL API.
// GET /api/gnomad-v4?chrom=7&pos37=140453136&ref=A&alt=T
// Performs GRCh37→GRCh38 liftover via Ensembl, then queries gnomAD v4.1.
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const GNOMAD_API = 'https://gnomad.broadinstitute.org/api';
const ENSEMBL_REST = 'https://rest.ensembl.org';

// Dataset hardcoded inline (not as a variable) to avoid DatasetId enum type issues.
const VARIANT_QUERY = `
query GnomadVariant($variantId: String!) {
  variant(variantId: $variantId, dataset: gnomad_r4) {
    variant_id
    chrom
    pos
    ref
    alt
    genome {
      ac
      an
      af
      populations {
        id
        ac
        an
        af
      }
    }
    exome {
      ac
      an
      af
      populations {
        id
        ac
        an
        af
      }
    }
  }
}
`;

async function liftoverHg19ToHg38(chrom, pos) {
    const c = String(chrom).replace(/^chr/i, '');
    const p = parseInt(pos, 10);
    const url = `${ENSEMBL_REST}/map/human/GRCh37/${c}:${p}..${p}:1/GRCh38?content-type=application/json`;
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 8000);
    try {
        const res = await fetch(url, {
            headers: { 'Accept': 'application/json' },
            signal: controller.signal,
        });
        if (!res.ok) return null;
        const data = await res.json();
        if (data && data.mappings && data.mappings.length > 0) {
            const mapped = data.mappings[0].mapped;
            if (mapped && mapped.start) return mapped.start;
        }
        return null;
    } catch {
        return null;
    } finally {
        clearTimeout(timer);
    }
}

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET, OPTIONS');
    if (req.method === 'OPTIONS') { res.status(204).end(); return; }

    const { chrom, pos37, ref, alt } = req.query;
    if (!chrom || !pos37 || !ref || !alt) {
        return res.status(200).json({ status: 'error', message: 'Missing required parameters: chrom, pos37, ref, alt' });
    }

    // Step 1: liftover GRCh37 → GRCh38
    const pos38 = await liftoverHg19ToHg38(chrom, pos37);
    if (!pos38) {
        return res.status(200).json({ status: 'liftover_failed', message: `Could not map ${chrom}:${pos37} to GRCh38` });
    }

    const c = String(chrom).replace(/^chr/i, '');
    const variantId = `${c}-${pos38}-${ref.toUpperCase()}-${alt.toUpperCase()}`;

    // Step 2: query gnomAD v4 GraphQL API
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 10000);
    try {
        const response = await fetch(GNOMAD_API, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Accept': 'application/json',
                'Origin': 'https://gnomad.broadinstitute.org',
                'Referer': 'https://gnomad.broadinstitute.org/',
                'User-Agent': 'Mozilla/5.0 (compatible; variant-search/1.0)',
            },
            body: JSON.stringify({
                query: VARIANT_QUERY,
                variables: { variantId },
            }),
            signal: controller.signal,
        });

        const text = await response.text();
        let body;
        try { body = JSON.parse(text); } catch { body = {}; }

        if (!response.ok) {
            return res.status(200).json({
                status: 'api_error',
                message: `gnomAD API HTTP ${response.status}`,
                grch38Id: variantId,
                detail: text.slice(0, 300),
            });
        }
        if (body.errors && body.errors.length > 0) {
            return res.status(200).json({
                status: 'api_error',
                message: body.errors.map(e => e.message).join('; '),
                grch38Id: variantId,
            });
        }
        const variant = body.data && body.data.variant;
        if (!variant) {
            return res.status(200).json({ status: 'not_found', grch38Id: variantId });
        }
        return res.status(200).json({ status: 'found', grch38Id: variantId, data: variant });
    } catch (err) {
        return res.status(200).json({
            status: 'error',
            message: err.message || String(err),
            grch38Id: variantId,
        });
    } finally {
        clearTimeout(timer);
    }
}
