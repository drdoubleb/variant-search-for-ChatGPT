// Vercel edge-cache headers for the data proxies.
//
// The `_` prefix keeps Vercel from treating this file as a routable serverless
// function — it is imported by the sibling handlers only.
//
// Why this exists: hot variants (BRAF V600E, EGFR L858R, TP53 R175H…) produce
// byte-identical ClinVar/CIViC/gnomAD/PubMed responses on every lookup.
// `s-maxage` lets Vercel's CDN answer those repeats without invoking the
// function or touching the upstream at all — faster cards, far less pressure
// on NCBI's shared rate limit, and cached variants keep working through an
// upstream outage. `stale-while-revalidate` serves the stale copy instantly
// while refreshing in the background once the TTL lapses.
//
// Set the header ONLY on responses that are safe to share: successful lookups
// and deterministic negatives ("no record", "unmappable"). Transient failures
// (upstream 5xx, timeouts) must not be pinned in the CDN — use setNoStore in
// catch blocks (the CDN skips most non-2xx anyway; the explicit no-store also
// keeps intermediaries honest).

export function setEdgeCache(res, seconds, swrSeconds = 86400) {
    res.setHeader('Cache-Control', `public, s-maxage=${seconds}, stale-while-revalidate=${swrSeconds}`);
    // The origin allowlist (api/_origin.js) runs inside the function. Without
    // Vary the CDN would serve a cached 200 to a disallowed origin and bypass
    // it; with Vary each origin gets its own cache entry and disallowed ones
    // still reach the function's 403.
    res.setHeader('Vary', 'Origin');
}

export function setNoStore(res) {
    res.setHeader('Cache-Control', 'no-store');
}
