const DEFAULT_DATASET_CANDIDATES = [
  // Direct static CSV downloads from the NCI TP53 Database (release r21, Jan 2025).
  // The view_data?bq_view_name=* URLs render HTML pages, not CSV files.
  // The actual data files are served from the /static/data/ path.
  'https://tp53.cancer.gov/static/data/MutationView_r21.csv',
  'https://tp53.cancer.gov/static/data/TumorVariantDownload_r21.csv',
  'https://tp53.cancer.gov/static/data/GermlineDownload_r21.csv'
];

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

async function fetchDatasetText() {
  const envUrl = process.env.TP53_MUTATION_DATASET_URL;
  const candidates = envUrl ? [envUrl, ...DEFAULT_DATASET_CANDIDATES] : DEFAULT_DATASET_CANDIDATES;
  const debugAttempts = [];
  for (const url of candidates) {
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
  return { text: null, datasetUrl: null, debugAttempts };
}

async function ensureDatasetLoaded() {
  const now = Date.now();
  if (cache.rows && now - cache.loadedAt < CACHE_TTL_MS) {
    return cache;
  }
  const fetched = await fetchDatasetText();
  if (!fetched || !fetched.text) {
    throw new Error(
      `Unable to fetch TP53 dataset. Attempts: ${JSON.stringify((fetched && fetched.debugAttempts) || [])}`
    );
  }
  const rows = parseCsv(fetched.text);
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
    if (proteinNeedle && protCompact && protCompact.includes(proteinNeedle)) score += 4;
    if (cdnaNeedle && cdNorm && cdNorm.includes(cdnaNeedle)) score += 5;
    if (genomicPos !== null) {
      if (String(hg19).includes(String(genomicPos))) score += 2;
      if (String(hg38).includes(String(genomicPos))) score += 2;
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
    const best = top[0] || null;

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
      match_count: matches.length,
      pathogenicity: pathogenicitySummary,
      classification: best?.exon_intron || '',
      prevalence: best?.somatic_count ? `somatic count: ${best.somatic_count}` : '',
      note: matches.length
        ? 'Matched against TP53 MutationView dataset using protein/cDNA/genomic similarity.'
        : 'No TP53 dataset row matched the supplied variant.',
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
