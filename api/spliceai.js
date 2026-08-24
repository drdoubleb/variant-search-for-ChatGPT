// Serverless proxy for the Broad Institute SpliceAI Lookup API.
// GET /api/spliceai?variant=chr7-140453136-A-T&hg=37&distance=500&mask=0&bc=basic
// Upstream API docs: https://github.com/broadinstitute/SpliceAI-lookup
// The upstream service is intended for interactive use only; do not batch query.

const SPLICEAI_API_BASE_BY_HG = {
    '37': 'https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/',
    '38': 'https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/'
};

function sanitizeVariant(value) {
    const variant = String(value || '').trim();
    if (!/^(?:chr)?[0-9XYMT]{1,2}[-\s:]+[0-9]{1,9}[-\s:]+[ACGT]+[-\s:>]+[ACGT]+$/i.test(variant)) {
        return '';
    }
    const parts = variant.replace(/^chr/i, 'chr').split(/[-\s:>]+/).filter(Boolean);
    if (parts.length !== 4) return '';
    const chrom = parts[0].replace(/^chr/i, 'chr');
    return [chrom, parts[1], parts[2].toUpperCase(), parts[3].toUpperCase()].join('-');
}

function sanitizeInteger(value, fallback, min, max) {
    const n = parseInt(value, 10);
    if (!Number.isFinite(n)) return fallback;
    return Math.min(Math.max(n, min), max);
}

import { rejectDisallowedOrigin } from './_origin.js';
import { setEdgeCache } from './_cache.js';

export default async function handler(req, res) {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

    if (req.method === 'OPTIONS') return res.status(204).end();
    if (req.method !== 'GET') return res.status(405).json({ error: 'Method not allowed' });
    // Block third-party websites from using this deployment as their backend
    // (no-Origin callers pass — see api/_origin.js).
    if (rejectDisallowedOrigin(req, res)) return;

    const variant = sanitizeVariant(req.query.variant);
    if (!variant) return res.status(400).json({ error: 'Missing or invalid variant. Expected chr-pos-ref-alt.' });

    const hg = String(req.query.hg || '37');
    if (!SPLICEAI_API_BASE_BY_HG[hg]) return res.status(400).json({ error: 'Invalid hg value. Use 37 or 38.' });

    const distance = sanitizeInteger(req.query.distance, 500, 1, 10000);
    const mask = String(req.query.mask ?? '0') === '1' ? '1' : '0';
    const bc = String(req.query.bc || 'basic') === 'comprehensive' ? 'comprehensive' : 'basic';

    const params = new URLSearchParams({ hg, variant, distance: String(distance), mask, bc });
    const url = `${SPLICEAI_API_BASE_BY_HG[hg]}?${params}`;
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 20000);

    try {
        const response = await fetch(url, {
            headers: { 'Accept': 'application/json', 'User-Agent': 'variant-search-tool/1.0' },
            signal: controller.signal
        });
        const text = await response.text();
        let data;
        try { data = JSON.parse(text); } catch { data = { raw: text }; }

        if (!response.ok) {
            return res.status(200).json({ error: `SpliceAI Lookup returned ${response.status}`, detail: text.slice(0, 500), variant, hg, distance, mask, bc });
        }

        // SpliceAI scores are deterministic per variant — cache aggressively.
        setEdgeCache(res, 604800);
        return res.status(200).json({ variant, hg, distance, mask, bc, data });
    } catch (err) {
        const message = err?.name === 'AbortError' ? 'SpliceAI Lookup request timed out' : (err.message || 'SpliceAI Lookup failed');
        return res.status(200).json({ error: message, variant, hg, distance, mask, bc });
    } finally {
        clearTimeout(timer);
    }
}
