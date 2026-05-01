const DEFAULT_DATASET_CANDIDATES = [
  // Hard-coded TP53 database endpoints (May 1, 2026).
  // These map to the official dataset views exposed by tp53.cancer.gov.
  'https://tp53.cancer.gov/view_data?bq_view_name=MutationView',
  'https://tp53.cancer.gov/view_data?bq_view_name=MutationViewDownload',
  'https://tp53.cancer.gov/view_data?bq_view_name=TumorVariantDownload',
  'https://tp53.cancer.gov/view_data?bq_view_name=GermlineDownload'
];

const CACHE_TTL_MS = 1000 * 60 * 60 * 6; // 6 hours
const MAX_RESPONSE_MATCHES = 5;

let cache = {
  loadedAt: 0,
  rows: null,
  datasetUrl: null
};

function normalizeKey(k) {
  return String(k || '')
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '');
}

function normalizeText(v) {
  return String(v || '').trim().toLowerCase();
}

function parseCsvLine(line) {
  const out = [];
  let cur = '';
  let inQuotes = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (ch === '"') {
      if (inQuotes && line[i + 1] === '"') {
        cur += '"';
        i++;
      } else {
        inQuotes = !inQuotes;
      }
    } else if (ch === ',' && !inQuotes) {
      out.push(cur);
      cur = '';
    } else {
      cur += ch;
    }
  }
  out.push(cur);
  return out;
}

function parseCsv(text) {
  const lines = String(text || '')
    .replace(/\r\n/g, '\n')
    .replace(/\r/g, '\n')
    .split('\n')
    .filter(Boolean);
  if (lines.length < 2) return [];
  const headers = parseCsvLine(lines[0]).map(h => h.trim());
  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    const cols = parseCsvLine(lines[i]);
    const row = {};
    headers.forEach((h, idx) => {
      row[h] = cols[idx] !== undefined ? cols[idx] : '';
    });
    rows.push(row);
  }
  return rows;
}

function pickValue(row, normalizedKeys) {
  const entries = Object.entries(row);
  for (const [k, v] of entries) {
    if (normalizedKeys.includes(normalizeKey(k))) {
      return v;
    }
  }
  return '';
}

function stripProteinPrefix(p) {
  const s = String(p || '').trim();
  return s.replace(/^p\./i, '');
}

function toProteinCompact(p) {
  return stripProteinPrefix(p).replace(/\s+/g, '').toUpperCase();
}

function parseGenomicPosition(g) {
  const m = String(g || '').match(/:g\.(\d+)/i);
  if (!m) return null;
  const n = Number(m[1]);
  return Number.isFinite(n) ? n : null;
}

async function fetchDatasetText() {
  const envUrl = process.env.TP53_MUTATION_DATASET_URL;
  const candidates = envUrl ? [envUrl, ...DEFAULT_DATASET_CANDIDATES] : DEFAULT_DATASET_CANDIDATES;
  for (const url of candidates) {
    try {
      const res = await fetch(url, {
        headers: { 'User-Agent': 'variant-search-tp53-proxy/1.0' }
      });
      if (!res.ok) continue;
      const text = await res.text();
      if (text && text.includes(',')) {
        return { text, datasetUrl: url };
      }
    } catch {
      // try next URL
    }
  }
  return null;
}

async function ensureDatasetLoaded() {
  const now = Date.now();
  if (cache.rows && now - cache.loadedAt < CACHE_TTL_MS) {
    return cache;
  }
  const fetched = await fetchDatasetText();
  if (!fetched) {
    throw new Error(
      'Unable to fetch TP53 dataset. Set TP53_MUTATION_DATASET_URL to the direct CSV URL from tp53.cancer.gov.'
    );
  }
  const rows = parseCsv(fetched.text);
  cache = {
    loadedAt: now,
    rows,
    datasetUrl: fetched.datasetUrl
  };
  return cache;
}

function buildMatches(rows, { protein, cdna, genomic }) {
  const proteinNeedle = toProteinCompact(protein);
  const cdnaNeedle = normalizeText(cdna);
  const genomicPos = parseGenomicPosition(genomic);

  const matches = [];
  for (const row of rows) {
    const prot = pickValue(row, ['protdescription', 'aachangeinhuman']);
    const cd = pickValue(row, ['cdescription', 'gdescription']);
    const hg19 = pickValue(row, ['hg19chr17coordinates']);
    const hg38 = pickValue(row, ['hg38chr17coordinates']);
    const somaticCount = pickValue(row, ['somaticcount']);
    const germlineCount = pickValue(row, ['germlinecount']);
    const exonIntron = pickValue(row, ['exonintron']);
    const mutId = pickValue(row, ['mutid']);

    const protCompact = toProteinCompact(prot);
    const cdNorm = normalizeText(cd);

    let score = 0;
    if (proteinNeedle && protCompact && protCompact.includes(proteinNeedle)) score += 4;
    if (cdnaNeedle && cdNorm && cdNorm.includes(cdnaNeedle)) score += 5;
    if (genomicPos !== null) {
      if (String(hg19).includes(String(genomicPos))) score += 2;
      if (String(hg38).includes(String(genomicPos))) score += 2;
    }
    if (score === 0) continue;

    matches.push({
      score,
      mut_id: mutId || '',
      protein: prot || '',
      cdna_or_genomic: cd || '',
      hg19_coordinate: hg19 || '',
      hg38_coordinate: hg38 || '',
      somatic_count: somaticCount || '',
      germline_count: germlineCount || '',
      exon_intron: exonIntron || ''
    });
  }
  matches.sort((a, b) => b.score - a.score);
  return matches;
}

export default async function handler(req, res) {
  if (req.method === 'OPTIONS') {
    res.setHeader('Access-Control-Allow-Origin', '*');
    res.setHeader('Access-Control-Allow-Methods', 'POST,GET,OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');
    return res.status(204).end();
  }

  if (req.method === 'GET') {
    return res.status(200).json({
      ok: true,
      service: 'tp53-proxy',
      usage: 'POST JSON: { gene, protein, cdna, genomic }'
    });
  }

  if (req.method !== 'POST') {
    return res.status(405).json({ error: 'Method not allowed' });
  }

  try {
    const body = req.body || {};
    const gene = String(body.gene || '').trim().toUpperCase();
    if (gene && gene !== 'TP53') {
      return res.status(200).json({
        source: 'tp53-database',
        match_count: 0,
        note: 'Gene is not TP53; TP53 database query skipped.'
      });
    }

    const { rows, datasetUrl } = await ensureDatasetLoaded();
    const matches = buildMatches(rows, body);
    const top = matches.slice(0, MAX_RESPONSE_MATCHES);
    const best = top[0] || null;

    return res.status(200).json({
      source: 'tp53-database',
      dataset_url: datasetUrl,
      total_records: rows.length,
      match_count: matches.length,
      classification: best?.exon_intron || '',
      prevalence: best?.somatic_count ? `somatic count: ${best.somatic_count}` : '',
      note: matches.length
        ? 'Matched against TP53 MutationView dataset using protein/cDNA/genomic similarity.'
        : 'No TP53 dataset row matched the supplied variant.',
      matches: top
    });
  } catch (err) {
    return res.status(502).json({
      error: 'TP53 backend query failed',
      detail: err?.message || String(err)
    });
  }
}
