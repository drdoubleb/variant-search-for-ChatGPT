// Serverless proxy for the gnomAD v4 GraphQL API.
// GET /api/gnomad-v4?chrom=7&pos37=140453136&ref=A&alt=T
// Performs GRCh37→GRCh38 liftover via Ensembl, then queries gnomAD v4.1.
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.

const GNOMAD_API = 'https://gnomad.broadinstitute.org/api';
const ENSEMBL_REST = 'https://rest.ensembl.org';

// gnomAD v4 variant query.
// Notes from schema introspection:
//   - variant() does NOT accept a reference_genome argument
//   - VariantPopulation does NOT have an af field; compute af = ac/an client-side
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
      }
    }
  }
}
`;

// Region query used as an indel fallback. The position-only liftover can land a
// multi-base variant a few bases away from how gnomAD indexes it, so when the exact
// ID misses we scan a small window for a record carrying the same alleles.
const REGION_QUERY = `
query GnomadRegion($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      pos
      ref
      alt
      genome { ac an af populations { id ac an } }
      exome { ac an af populations { id ac an } }
    }
  }
}
`;

// Map a single GRCh37 position to GRCh38 via Ensembl. NOTE: this lifts only the start
// coordinate (pos..pos) and the caller reuses the original ref/alt verbatim, so the
// resulting variant ID is reliable for SNVs but not necessarily for indels/MNVs, whose
// ref allele spans multiple bases and whose GRCh38 representation can differ. It also
// takes mappings[0] without validating strand/length.
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

async function gnomadPost(operationName, query, variables) {
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
            body: JSON.stringify({ operationName, query, variables }),
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
    const refU = String(ref).toUpperCase();
    const altU = String(alt).toUpperCase();
    const variantId = `${c}-${pos38}-${refU}-${altU}`;

    // Flag non-SNV lookups: the position-only liftover + verbatim ref/alt can produce an
    // inexact GRCh38 ID for indels/MNVs, so a not-found result is not conclusive for them.
    const isSnv = /^[ACGT]$/.test(refU) && /^[ACGT]$/.test(altU);
    const caveat = isSnv
        ? undefined
        : 'Liftover mapped the start position only and reused the ref/alt alleles. For indels/MNVs the GRCh38 variant ID may be inexact; a nearby-region scan for the same alleles runs when the exact ID misses, but a not-found result is still not fully conclusive.';

    // Step 2: query gnomAD v4
    const result = await gnomadPost('GnomadVariant', VARIANT_QUERY, {
        variantId,
        dataset: 'gnomad_r4',
    });

    if (!result.ok) {
        return res.status(200).json({
            status: 'api_error',
            message: `gnomAD API HTTP ${result.status}`,
            grch38Id: variantId,
            detail: result.text.slice(0, 500),
            ...(caveat ? { caveat } : {}),
        });
    }

    let body;
    try { body = JSON.parse(result.text); } catch { body = {}; }

    if (body.errors && body.errors.length > 0) {
        return res.status(200).json({
            status: 'api_error',
            message: body.errors.map(e => e.message).join('; '),
            grch38Id: variantId,
            ...(caveat ? { caveat } : {}),
        });
    }

    let variant = body.data && body.data.variant;
    let matchedVia;

    // Step 3 (non-SNVs only): the exact ID missed, which for an indel may just mean
    // the liftover shifted the anchor. Scan a window around the lifted position for a
    // record with identical alleles before reporting "not found".
    if (!variant && !isSnv) {
        const span = Math.max(refU.length, altU.length);
        const start = Math.max(1, pos38 - span - 5);
        const stop = pos38 + span + 5;
        const regionResult = await gnomadPost('GnomadRegion', REGION_QUERY, { chrom: c, start, stop });
        if (regionResult.ok) {
            let regionBody;
            try { regionBody = JSON.parse(regionResult.text); } catch { regionBody = {}; }
            const candidates = regionBody?.data?.region?.variants || [];
            const hit = candidates.find(v => String(v.ref).toUpperCase() === refU && String(v.alt).toUpperCase() === altU);
            if (hit) {
                variant = { ...hit, chrom: c };
                matchedVia = `region scan ${c}:${start}-${stop} (exact ID ${variantId} not indexed)`;
            }
        }
    }

    if (!variant) {
        return res.status(200).json({ status: 'not_found', grch38Id: variantId, ...(caveat ? { caveat } : {}) });
    }

    return res.status(200).json({
        status: 'found',
        grch38Id: variant.variant_id || variantId,
        ...(caveat ? { caveat } : {}),
        ...(matchedVia ? { matchedVia } : {}),
        data: {
            variant_id: variant.variant_id,
            chrom: variant.chrom,
            pos: variant.pos,
            ref: variant.ref,
            alt: variant.alt,
            genome: variant.genome || null,
            exome: variant.exome || null,
        },
    });
}
