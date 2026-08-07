// Serverless proxy for the gnomAD v4 GraphQL API.
// GET /api/gnomad-v4?chrom=7&pos37=140453136&ref=A&alt=T[&pos38=140753336]
// Performs GRCh37→GRCh38 liftover via Ensembl, then queries gnomAD v4.1.
// Always returns HTTP 200 — errors go into the body so the card degrades cleanly.
//
// `pos38` is an optional hint: when the caller already knows the GRCh38 start
// (myvariant.info returns one in dbnsfp.hg38 / dbsnp.hg38 for most indexed
// variants) we use it and skip the Ensembl round-trip entirely. That makes the
// common path faster and immune to Ensembl outages.

import { liftoverPosition, normalizeChrom } from './_liftover.js';

const GNOMAD_API = 'https://gnomad.broadinstitute.org/api';

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

// Resolve the GRCh38 start coordinate for a GRCh37 position.
//
// NOTE (unchanged by the liftover work): only the start coordinate is lifted, and
// the caller reuses the original ref/alt verbatim, so the resulting variant ID is
// reliable for SNVs but not necessarily for indels/MNVs, whose ref allele spans
// multiple bases and whose GRCh38 representation can differ. The region-scan
// fallback further down exists to cover that.
//
// Returns the same outcome shape as liftoverPosition().
async function resolvePos38(chrom, pos37, pos38Hint) {
    // Trust a well-formed caller-supplied GRCh38 start and skip Ensembl. This is
    // the same coordinate the liftover would compute, sourced from the annotation
    // the client already holds.
    const hint = Number.parseInt(pos38Hint, 10);
    if (Number.isInteger(hint) && hint > 0) {
        return { ok: true, pos: hint, host: 'caller-supplied' };
    }
    return liftoverPosition(chrom, pos37, { from: 'GRCh37', to: 'GRCh38' });
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

    const { chrom, pos37, ref, alt, pos38: pos38Hint } = req.query;
    if (!chrom || !pos37 || !ref || !alt) {
        return res.status(200).json({ status: 'error', message: 'Missing required parameters: chrom, pos37, ref, alt' });
    }

    // Step 1: resolve the GRCh38 position (caller hint, else liftover).
    // A provider outage and a genuinely unmappable coordinate get distinct
    // statuses — reporting an Ensembl 500 as "liftover failed for <coord>" blames
    // the variant for someone else's downtime and hides the fact that retrying
    // later would work.
    const lift = await resolvePos38(chrom, pos37, pos38Hint);
    if (!lift.ok) {
        return res.status(200).json({
            status: lift.reason === 'unavailable' ? 'liftover_unavailable' : 'liftover_failed',
            message: lift.reason === 'unavailable'
                ? `Ensembl assembly-mapping API is unreachable, so ${chrom}:${pos37} could not be converted to GRCh38.`
                : `Could not map ${chrom}:${pos37} to GRCh38`,
            ...(lift.detail ? { detail: lift.detail } : {}),
        });
    }
    const pos38 = lift.pos;

    const c = normalizeChrom(chrom);
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

    // gnomAD reports an absent variant as a GraphQL *error* ("Variant not found")
    // alongside data.variant: null, not as a plain null. Treating that as an API
    // error was wrong twice over: it showed "gnomAD v4.1 error: Variant not found"
    // where the clean "Not found in gnomAD v4.1" belongs, and — because the check
    // returned early — it made the indel region-scan below unreachable, which is
    // precisely the fallback that covers a liftover landing an indel a few bases
    // off. Absence is a normal answer here; keep going.
    const errors = (body.errors || []).filter((e) => !/variant not found/i.test((e && e.message) || ''));
    if (errors.length > 0) {
        return res.status(200).json({
            status: 'api_error',
            message: errors.map(e => e.message).join('; '),
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
