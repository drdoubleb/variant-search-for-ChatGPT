// Serverless proxy for the gnomAD v4 GraphQL API.
// GET /api/gnomad-v4?chrom=7&pos37=140453136&ref=A&alt=T
// Performs GRCh37→GRCh38 liftover via Ensembl, then queries gnomAD v4.1.
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const GNOMAD_API = 'https://gnomad.broadinstitute.org/api';
const ENSEMBL_REST = 'https://rest.ensembl.org';

// Only genome is queried for gnomad_r4 (v4 is genome-sequencing-based).
// Exome data (if present) uses a separate dataset. Populations are fetched
// alongside the main AF so the population table can be populated.
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
  }
}
`;

// gnomAD v4 also has joint exome+genome data under a separate dataset ID.
// Query it separately so we can display both if available.
const EXOME_QUERY = `
query GnomadVariantExome($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
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

async function gnomadPost(query, variables) {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 10000);
    try {
        const response = await fetch(GNOMAD_API, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Accept': 'application/json',
                'User-Agent': 'GeneScape/1.0',
            },
            body: JSON.stringify({ query, variables }),
            signal: controller.signal,
        });
        const text = await response.text();
        return { ok: response.ok, status: response.status, text };
    } catch (err) {
        return { ok: false, status: 0, text: err.message || String(err) };
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

    // Step 2: query gnomAD v4 for genome data (primary dataset)
    const genomeFetch = await gnomadPost(VARIANT_QUERY, { variantId, dataset: 'gnomad_r4' });
    if (!genomeFetch.ok) {
        return res.status(200).json({
            status: 'api_error',
            message: `gnomAD API HTTP ${genomeFetch.status}`,
            grch38Id: variantId,
            detail: genomeFetch.text.slice(0, 400),
        });
    }

    let genomeBody;
    try { genomeBody = JSON.parse(genomeFetch.text); } catch { genomeBody = {}; }

    if (genomeBody.errors && genomeBody.errors.length > 0) {
        return res.status(200).json({
            status: 'api_error',
            message: genomeBody.errors.map(e => e.message).join('; '),
            grch38Id: variantId,
        });
    }

    const genomeVariant = genomeBody.data && genomeBody.data.variant;
    if (!genomeVariant) {
        return res.status(200).json({ status: 'not_found', grch38Id: variantId });
    }

    // Step 3: optionally query exome dataset (gnomad_r4_non_ukb has exomes;
    // gnomad_r4 genome-only will return null exome — that's fine).
    let exomeData = null;
    const exomeFetch = await gnomadPost(EXOME_QUERY, { variantId, dataset: 'gnomad_r4' });
    if (exomeFetch.ok) {
        try {
            const exomeBody = JSON.parse(exomeFetch.text);
            exomeData = exomeBody.data && exomeBody.data.variant && exomeBody.data.variant.exome;
        } catch { /* ignore */ }
    }

    return res.status(200).json({
        status: 'found',
        grch38Id: variantId,
        data: {
            variant_id: genomeVariant.variant_id,
            chrom: genomeVariant.chrom,
            pos: genomeVariant.pos,
            ref: genomeVariant.ref,
            alt: genomeVariant.alt,
            genome: genomeVariant.genome || null,
            exome: exomeData || null,
        },
    });
}
