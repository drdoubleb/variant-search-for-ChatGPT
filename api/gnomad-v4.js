// Serverless proxy for the gnomAD v4 GraphQL API.
// GET /api/gnomad-v4?variantId=7-140453136-A-T
// variantId must use GRCh38 coordinates in the format chrom-pos-ref-alt.
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const GNOMAD_API = 'https://gnomad.broadinstitute.org/api';

const VARIANT_QUERY = `
query GnomadVariant($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
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

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET, OPTIONS');
    if (req.method === 'OPTIONS') { res.status(204).end(); return; }

    const variantId = (req.query.variantId || '').trim();
    if (!variantId) {
        return res.status(200).json({ error: 'Missing variantId parameter' });
    }

    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 12000);

    try {
        const response = await fetch(GNOMAD_API, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Accept': 'application/json',
            },
            body: JSON.stringify({
                query: VARIANT_QUERY,
                variables: { variantId, dataset: 'gnomad_r4' },
            }),
            signal: controller.signal,
        });

        const text = await response.text();
        let body;
        try { body = JSON.parse(text); } catch { body = {}; }

        if (!response.ok) {
            return res.status(200).json({ error: `gnomAD API error ${response.status}: ${text.slice(0, 300)}` });
        }
        if (body.errors && body.errors.length > 0) {
            return res.status(200).json({ error: body.errors.map(e => e.message).join('; '), data: body.data || null });
        }
        return res.status(200).json({ data: body.data || null });
    } catch (err) {
        return res.status(200).json({ error: err.message || String(err) });
    } finally {
        clearTimeout(timer);
    }
}
