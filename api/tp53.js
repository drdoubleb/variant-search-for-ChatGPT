// Last-resort fallback URLs (release r21, Jan 2025). Discovery (below) normally
// supplies the current release; these are used only if the download page can't be
// reached or parsed. The view_data?bq_view_name=* URLs render HTML, not CSV — the
// actual files live under /static/data/.
const DEFAULT_DATASET_CANDIDATES = [
  'https://tp53.cancer.gov/static/data/MutationView_r21.csv',
  'https://tp53.cancer.gov/static/data/TumorVariantDownload_r21.csv',
  'https://tp53.cancer.gov/static/data/GermlineDownload_r21.csv'
];

// The NCI TP53 database download page lists the current release's CSV files (e.g.
// MutationView_r21.csv). Scraping it lets a new release be picked up automatically
// instead of 404-ing against a hard-coded release number.
const TP53_DOWNLOAD_PAGE = 'https://tp53.cancer.gov/get_tp53data';

// Discover current-release CSV URLs from the download page. Returns absolute URLs for
// the datasets the proxy can parse (MutationView first — it carries the columns
// buildMatches needs), or [] when discovery fails so callers fall back to the
// hard-coded list. The MutationView pattern is anchored on a leading slash so it does
// not also match InducedMutationView.
async function discoverDatasetUrls() {
  try {
    const res = await fetch(TP53_DOWNLOAD_PAGE, {
      headers: {
        'User-Agent': 'Mozilla/5.0 (compatible; variant-search-tp53-proxy/1.0)',
        'Accept': 'text/html,*/*'
      }
    });
    if (!res.ok) return [];
    const html = await res.text();
    const found = [...String(html).matchAll(/static\/data\/[A-Za-z0-9_]+_r\d+\.csv/gi)].map(m => m[0]);
    const pick = (re) => {
      const hit = found.find(h => re.test(h));
      if (!hit) return null;
      try { return new URL(hit, TP53_DOWNLOAD_PAGE).toString(); } catch { return null; }
    };
    const urls = [
      pick(/\/MutationView_r\d+\.csv$/i),
      pick(/\/TumorVariantDownload_r\d+\.csv$/i),
      pick(/\/GermlineDownload_r\d+\.csv$/i)
    ].filter(Boolean);
    return Array.from(new Set(urls));
  } catch {
    return [];
  }
}

const CACHE_TTL_MS = 1000 * 60 * 60 * 6; // 6 hours
const MAX_RESPONSE_MATCHES = 5;

let cache = {
  loadedAt: 0,
  rows: null,
  datasetUrl: null,
  debug: null
};

function normalizeKey(k) {
  return String(k || '')
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '');
}

function normalizeText(v) {
  return String(v || '').trim().toLowerCase();
}

// Single-pass CSV parser that tracks quote state ACROSS newlines, so a quoted field
// containing a newline (common in the NCI TP53 free-text columns) does not split its
// row. Handles doubled-quote ("") escaping. Splitting on lines first (as the previous
// implementation did) corrupted column alignment for any such row.
function parseCsv(text) {
  const s = String(text || '').replace(/\r\n/g, '\n').replace(/\r/g, '\n');
  const rows = [];
  let row = [];
  let cur = '';
  let inQuotes = false;
  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    if (inQuotes) {
      if (ch === '"') {
        if (s[i + 1] === '"') { cur += '"'; i++; }
        else inQuotes = false;
      } else {
        cur += ch;
      }
    } else if (ch === '"') {
      inQuotes = true;
    } else if (ch === ',') {
      row.push(cur);
      cur = '';
    } else if (ch === '\n') {
      row.push(cur);
      cur = '';
      rows.push(row);
      row = [];
    } else {
      cur += ch;
    }
  }
  // Flush the final field/row when the text doesn't end in a newline.
  if (cur !== '' || row.length > 0) {
    row.push(cur);
    rows.push(row);
  }
  if (rows.length < 2) return [];
  const headers = rows[0].map(h => String(h).trim());
  const out = [];
  for (let r = 1; r < rows.length; r++) {
    const cols = rows[r];
    if (cols.length === 1 && cols[0] === '') continue; // skip blank rows
    const obj = {};
    headers.forEach((h, idx) => {
      obj[h] = cols[idx] !== undefined ? cols[idx] : '';
    });
    out.push(obj);
  }
  return out;
}

function looksLikeHtml(text) {
  return /<!doctype html|<html[\s>]|<table[\s>]/i.test(String(text || ''));
}

function extractDownloadLinkFromHtml(html, baseUrl) {
  const s = String(html || '');
  const hrefMatches = [...s.matchAll(/href\s*=\s*["']([^"']+)["']/gi)].map(m => m[1]);
  const candidates = hrefMatches.filter((h) =>
    /csv|download|export|view_data/i.test(h)
  );
  for (const href of candidates) {
    try {
      return new URL(href, baseUrl).toString();
    } catch {
      // continue
    }
  }
  // Common fallback patterns used by data-table download buttons.
  try {
    const u = new URL(baseUrl);
    if (!u.searchParams.has('download')) {
      u.searchParams.set('download', '1');
      return u.toString();
    }
  } catch {
    // ignore URL parse fallback
  }
  return null;
}

function buildColumnResolver(rows) {
  const first = rows && rows[0] ? rows[0] : null;
  const keys = first ? Object.keys(first) : [];
  const normalizedToActual = new Map();
  for (const key of keys) {
    const normalized = normalizeKey(key);
    if (normalized && !normalizedToActual.has(normalized)) {
      normalizedToActual.set(normalized, key);
    }
  }

  return function pickValue(row, normalizedKeys) {
    for (const normalizedKey of normalizedKeys) {
      const actualKey = normalizedToActual.get(normalizedKey);
      if (!actualKey) continue;
      const value = row[actualKey];
      if (value !== undefined && value !== null && value !== '') return value;
    }
    return '';
  };
}

function stripProteinPrefix(p) {
  const s = String(p || '').trim();
  return s.replace(/^p\./i, '');
}

const AA3_TO_1 = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
  HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
  THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V', TER: '*', STOP: '*'
};

function toProteinSingleLetter(value) {
  const raw = stripProteinPrefix(value).replace(/\s+/g, '');
  const m = raw.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|Stop|\*)$/i);
  if (!m) return raw.toUpperCase();
  const ref = AA3_TO_1[String(m[1]).toUpperCase()];
  const alt = AA3_TO_1[String(m[3]).toUpperCase()] || (m[3] === '*' ? '*' : null);
  if (!ref || !alt) return raw.toUpperCase();
  return `${ref}${m[2]}${alt}`.toUpperCase();
}

function toProteinCompact(p) {
  const raw = stripProteinPrefix(p).replace(/\s+/g, '').toUpperCase();
  const single = toProteinSingleLetter(p);
  // Keep both representations available for matching by returning the canonical
  // single-letter code when possible; otherwise preserve raw uppercase text.
  return single || raw;
}

function parseGenomicPosition(g) {
  const m = String(g || '').match(/:g\.(\d+)/i);
  if (!m) return null;
  const n = Number(m[1]);
  return Number.isFinite(n) ? n : null;
}

// The ref amino acid + codon number of a compact protein change (e.g. "R175H" -> "R175").
// Returns null when there is no parseable ref+number.
function proteinCodonKey(compact) {
  const m = String(compact || '').match(/^([A-Z*])(\d+)/);
  return m ? `${m[1]}${m[2]}` : null;
}

// Same codon (same ref residue + position), different/again alt — e.g. R175C vs R175H.
function sameProteinCodon(a, b) {
  const ka = proteinCodonKey(a);
  const kb = proteinCodonKey(b);
  return !!ka && ka === kb;
}

// True when the coordinate string contains `pos` as a whole integer token. Handles a
// single position or a range (e.g. "7578406-7578490"), and — unlike a substring test —
// does not let 7578406 match 17578406.
function coordinateHasPosition(coord, pos) {
  const nums = String(coord || '').match(/\d+/g);
  return Array.isArray(nums) && nums.some(n => Number(n) === pos);
}

async function fetchDatasetText() {
  const envUrl = process.env.TP53_MUTATION_DATASET_URL;
  // Priority: explicit env override → discovered current-release URLs → hard-coded r21.
  // Candidate groups are resolved lazily so a working env override skips the
  // discovery scrape (an extra fetch of the download page) entirely.
  const candidateGroups = [
    async () => (envUrl ? [envUrl] : []),
    discoverDatasetUrls,
    async () => DEFAULT_DATASET_CANDIDATES
  ];
  const debugAttempts = [];
  const tried = new Set();
  for (const group of candidateGroups) {
  let groupUrls = [];
  try { groupUrls = await group(); } catch { groupUrls = []; }
  for (const url of groupUrls) {
    if (!url || tried.has(url)) continue;
    tried.add(url);
    try {
      const res = await fetch(url, {
        headers: {
          'User-Agent': 'Mozilla/5.0 (compatible; variant-search-tp53-proxy/1.0)',
          'Accept': 'text/csv,text/plain,*/*'
        }
      });
      if (!res.ok) {
        debugAttempts.push({ url, ok: false, status: res.status, reason: 'non-200 response' });
        continue;
      }
      const text = await res.text();
      const contentType = res.headers.get('content-type') || '';
      const isHtml = looksLikeHtml(text) || /text\/html/i.test(contentType);
      if (text && text.includes(',') && !isHtml) {
        debugAttempts.push({ url, ok: true, status: res.status, contentType, rowsHint: text.split('\n').length });
        return { text, datasetUrl: url, debugAttempts };
      }

      // Some TP53 links render a webpage with a download button. Try to extract and follow it.
      if (isHtml) {
        const extracted = extractDownloadLinkFromHtml(text, url);
        debugAttempts.push({ url, ok: true, status: res.status, contentType, extractedDownloadUrl: extracted || '' });
        if (extracted) {
          try {
            const res2 = await fetch(extracted, {
              headers: {
                'User-Agent': 'Mozilla/5.0 (compatible; variant-search-tp53-proxy/1.0)',
                'Accept': 'text/csv,text/plain,*/*'
              }
            });
            const text2 = await res2.text();
            const ct2 = res2.headers.get('content-type') || '';
            const html2 = looksLikeHtml(text2) || /text\/html/i.test(ct2);
            if (res2.ok && text2 && text2.includes(',') && !html2) {
              debugAttempts.push({ url: extracted, ok: true, status: res2.status, contentType: ct2, rowsHint: text2.split('\n').length });
              return { text: text2, datasetUrl: extracted, debugAttempts };
            }
            debugAttempts.push({ url: extracted, ok: res2.ok, status: res2.status, contentType: ct2, reason: 'download candidate did not return CSV' });
          } catch (subErr) {
            debugAttempts.push({ url: extracted, ok: false, reason: `download fetch error: ${subErr?.message || String(subErr)}` });
          }
        }
      }
    } catch {
      debugAttempts.push({ url, ok: false, reason: 'fetch exception' });
      // try next URL
    }
  }
  }
  return { text: null, datasetUrl: null, debugAttempts };
}

async function ensureDatasetLoaded() {
  const now = Date.now();
  if (cache.rows && now - cache.loadedAt < CACHE_TTL_MS) {
    return cache;
  }
  let fetched = null;
  try {
    fetched = await fetchDatasetText();
  } catch {
    fetched = null;
  }
  const rows = fetched && fetched.text ? parseCsv(fetched.text) : null;
  if (!rows || rows.length === 0) {
    // Refresh failed (or returned no usable rows). Serve the last-good dataset if we
    // have one, rather than 502-ing the whole TP53 card — important when a release
    // bump temporarily breaks the URLs or the upstream is down.
    if (cache.rows && cache.rows.length) {
      cache.staleSince = cache.staleSince || now;
      return cache;
    }
    throw new Error(
      `Unable to fetch TP53 dataset. Attempts: ${JSON.stringify((fetched && fetched.debugAttempts) || [])}`
    );
  }
  cache = {
    loadedAt: now,
    rows,
    datasetUrl: fetched.datasetUrl,
    debug: fetched.debugAttempts || []
  };
  return cache;
}

function buildMatches(rows, { protein, cdna, genomic }) {
  const pickValue = buildColumnResolver(rows);
  const proteinNeedle = toProteinCompact(protein);
  const cdnaNeedle = normalizeText(cdna);
  const genomicPos = parseGenomicPosition(genomic);

  const matches = [];
  for (const row of rows) {
    // Core identification fields
    const prot = pickValue(row, ['protdescription', 'aachangeinhuman']);
    const cd = pickValue(row, ['cdescription', 'gdescription']);
    const hg19 = pickValue(row, ['hg19chr17coordinates']);
    const hg38 = pickValue(row, ['hg38chr17coordinates']);

    const protCompact = toProteinCompact(prot);
    const cdNorm = normalizeText(cd);

    let score = 0;
    // Protein: an exact token match scores highest. A same-codon match (same ref+
    // position, different alt — e.g. R175C when querying R175H) is surfaced at a lower
    // score so related variants appear without ever outranking the exact hit. The old
    // substring test let a short needle (e.g. "R175") match R175H/C/L indiscriminately.
    let sameCodon = false;
    if (proteinNeedle && protCompact) {
      if (protCompact === proteinNeedle) {
        score += 4;
      } else if (sameProteinCodon(protCompact, proteinNeedle)) {
        score += 1;
        sameCodon = true;
      }
    }
    // cDNA: substring is retained — cDNA strings are specific and formatting varies.
    if (cdnaNeedle && cdNorm && cdNorm.includes(cdnaNeedle)) score += 5;
    // Genomic: exact integer-token compare, not substring (7578406 must not match 17578406).
    if (genomicPos !== null) {
      if (coordinateHasPosition(hg19, genomicPos)) score += 2;
      if (coordinateHasPosition(hg38, genomicPos)) score += 2;
    }
    if (score === 0) continue;

    const somaticCount = pickValue(row, ['somaticcount']);
    const germlineCount = pickValue(row, ['germlinecount']);
    const exonIntron = pickValue(row, ['exonintron']);
    const mutId = pickValue(row, ['mutid']);
    const mutationType = pickValue(row, ['type', 'mutationtype', 'varianttype']);
    const codonNum = pickValue(row, ['codonnumber']);
    const cpgSite = pickValue(row, ['cpgsite']);
    const spliceSite = pickValue(row, ['splicesite']);
    const hotspot = pickValue(row, ['hotspot']);
    const effect = pickValue(row, ['effect']);
    const effectGroup = pickValue(row, ['effectgroup3', 'effectgroup']);
    const taClass = pickValue(row, ['transactivationclass', 'taclass']);
    const lofClass = pickValue(row, ['dnelofclass']);
    const dneClass = pickValue(row, ['dneclass']);
    const structFunc = pickValue(row, ['structurefunctionclass']);
    const residueFunc = pickValue(row, ['residuefunction']);
    const domainFunc = pickValue(row, ['domainfunction']);
    const structMotif = pickValue(row, ['structuralmotif']);
    const agvgd = pickValue(row, ['agvgdclass']);
    const sift = pickValue(row, ['siftclass']);
    const polyphen2 = pickValue(row, ['polyphen2']);
    const bayesDel = pickValue(row, ['bayesdel']);
    const revel = pickValue(row, ['revel']);
    const tcgaCount = pickValue(row, ['tcgaicgcgeniecount']);
    const taWaf1 = pickValue(row, ['waf1nwt']);
    const taMdm2 = pickValue(row, ['mdm2nwt']);
    const taBax = pickValue(row, ['baxnwt']);
    const taH1433s = pickValue(row, ['h1433snwt']);
    const taAip1 = pickValue(row, ['aip1nwt']);
    const taGadd45 = pickValue(row, ['gadd45nwt']);
    const taNoxa = pickValue(row, ['noxanwt']);
    const taP53r2 = pickValue(row, ['p53r2nwt']);
    const spliceAiDsAl = pickValue(row, ['spliceaidssal']);
    const spliceAiDsAg = pickValue(row, ['spliceaidsag']);
    const spliceAiDsDg = pickValue(row, ['spliceaidsdg']);
    const spliceAiDsDl = pickValue(row, ['spliceaidsdl']);

    matches.push({
      score,
      same_codon_match: sameCodon,
      mut_id: mutId || '',
      protein: prot || '',
      cdna_or_genomic: cd || '',
      hg19_coordinate: hg19 || '',
      hg38_coordinate: hg38 || '',
      somatic_count: somaticCount || '',
      germline_count: germlineCount || '',
      tcga_icgc_genie_count: tcgaCount || '',
      exon_intron: exonIntron || '',
      codon_number: codonNum || '',
      mutation_type: mutationType || '',
      cpg_site: cpgSite || '',
      splice_site: spliceSite || '',
      hotspot: hotspot || '',
      effect: effect || '',
      effect_group: effectGroup || '',
      ta_class: taClass || '',
      lof_class: lofClass || '',
      dne_class: dneClass || '',
      structural_class: structFunc || '',
      residue_function: residueFunc || '',
      domain_function: domainFunc || '',
      structural_motif: structMotif || '',
      agvgd_class: agvgd || '',
      sift_class: sift || '',
      polyphen2: polyphen2 || '',
      bayes_del: bayesDel || '',
      revel: revel || '',
      ta_targets: {
        WAF1: taWaf1 || '',
        MDM2: taMdm2 || '',
        BAX: taBax || '',
        '14-3-3σ': taH1433s || '',
        AIP1: taAip1 || '',
        GADD45: taGadd45 || '',
        NOXA: taNoxa || '',
        P53R2: taP53r2 || '',
      },
      splice_ai: {
        DS_AL: spliceAiDsAl || '',
        DS_AG: spliceAiDsAg || '',
        DS_DG: spliceAiDsDg || '',
        DS_DL: spliceAiDsDl || '',
      },
    });
  }
  matches.sort((a, b) => b.score - a.score);
  return matches;
}

import { rejectDisallowedOrigin } from './_origin.js';

export default async function handler(req, res) {
  res.setHeader('Access-Control-Allow-Origin', '*');
  res.setHeader('Access-Control-Allow-Methods', 'POST,GET,OPTIONS');
  res.setHeader('Access-Control-Allow-Headers', 'Content-Type');

  if (req.method === 'OPTIONS') {
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

  // Block third-party websites from using this deployment as their backend
  // (no-Origin callers pass — see api/_origin.js).
  if (rejectDisallowedOrigin(req, res)) return;

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

    const debugEnabled = String(body.debug || '').toLowerCase() === 'true' || body.debug === true;
    const { rows, datasetUrl, debug } = await ensureDatasetLoaded();
    const matches = buildMatches(rows, body);
    const top = matches.slice(0, MAX_RESPONSE_MATCHES);
    // match_count reflects strong matches (exact protein / cDNA / genomic, score >= 2).
    // Same-codon relatives (score 1) are surfaced in `matches` for context but counted
    // separately and never drive the pathogenicity summary.
    const strongMatches = matches.filter(m => m.score >= 2);
    const relatedCodonCount = matches.length - strongMatches.length;
    const best = strongMatches[0] || null;

    // Build a concise pathogenicity summary from the best match for top-level consumption
    const pathogenicitySummary = best ? {
      mutation_type: best.mutation_type || '',
      location: best.exon_intron || '',
      codon_number: best.codon_number || '',
      effect: best.effect || '',
      effect_group: best.effect_group || '',
      ta_class: best.ta_class || '',
      lof_class: best.lof_class || '',
      dne_class: best.dne_class || '',
      hotspot: best.hotspot || '',
      cpg_site: best.cpg_site || '',
      somatic_count: best.somatic_count || '',
      germline_count: best.germline_count || '',
      tcga_icgc_genie_count: best.tcga_icgc_genie_count || '',
      structural_class: best.structural_class || '',
      residue_function: best.residue_function || '',
      domain_function: best.domain_function || '',
      structural_motif: best.structural_motif || '',
      splice_site: best.splice_site || '',
      agvgd_class: best.agvgd_class || '',
      sift_class: best.sift_class || '',
      polyphen2: best.polyphen2 || '',
      bayes_del: best.bayes_del || '',
      revel: best.revel || '',
      ta_targets: best.ta_targets || {},
      splice_ai: best.splice_ai || {},
    } : null;

    const responsePayload = {
      source: 'tp53-database',
      dataset_url: datasetUrl,
      total_records: rows.length,
      match_count: strongMatches.length,
      related_codon_count: relatedCodonCount,
      pathogenicity: pathogenicitySummary,
      classification: best?.exon_intron || '',
      prevalence: best?.somatic_count ? `somatic count: ${best.somatic_count}` : '',
      note: strongMatches.length
        ? 'Matched against TP53 MutationView dataset by exact protein/cDNA/genomic position.'
        : (relatedCodonCount
          ? 'No exact TP53 record; showing same-codon variant(s) at this position for context.'
          : 'No TP53 dataset row matched the supplied variant.'),
      matches: top
    };
    if (debugEnabled) {
      const debugPickValue = buildColumnResolver(rows);
      responsePayload.debug = {
        query: {
          protein_input: body.protein || '',
          protein_normalized: toProteinCompact(body.protein || ''),
          cdna_input: body.cdna || '',
          genomic_input: body.genomic || ''
        },
        dataset_attempts: debug || [],
        sample_protein_values: rows.slice(0, 25).map(r => debugPickValue(r, ['protdescription', 'aachangeinhuman'])).filter(Boolean).slice(0, 10),
        available_columns: rows.length > 0 ? Object.keys(rows[0]) : []
      };
    }
    return res.status(200).json(responsePayload);
  } catch (err) {
    return res.status(502).json({
      error: 'TP53 backend query failed',
      detail: err?.message || String(err)
    });
  }
}
