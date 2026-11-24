import { tp53AssayComparison } from './_tp53-assays.js';

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
// Pure extraction core (no network), split out so the test suite can drive it
// with canned HTML.
function extractDatasetUrls(html, base = TP53_DOWNLOAD_PAGE) {
  const found = [...String(html).matchAll(/static\/data\/[A-Za-z0-9_]+_r\d+\.csv/gi)].map(m => m[0]);
  const pick = (re) => {
    const hit = found.find(h => re.test(h));
    if (!hit) return null;
    try { return new URL(hit, base).toString(); } catch { return null; }
  };
  return [
    pick(/\/MutationView_r\d+\.csv$/i),
    pick(/\/TumorVariantDownload_r\d+\.csv$/i),
    pick(/\/GermlineDownload_r\d+\.csv$/i)
  ].filter(Boolean);
}

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
    return Array.from(new Set(extractDatasetUrls(html)));
  } catch {
    return [];
  }
}

const CACHE_TTL_MS = 1000 * 60 * 60 * 6; // 6 hours
const MAX_RESPONSE_MATCHES = 5;
const GENOMIC_COLUMNS = {
  GRCh37: { accession: 'NC_000017.10', hgvs: 'hgvshg19', description: 'gdescription', coordinate: 'hg19chr17coordinates' },
  GRCh38: { accession: 'NC_000017.11', hgvs: 'hgvshg38', description: 'gdescriptiongrch38', coordinate: 'hg38chr17coordinates' }
};
const MATCH_RANK = { exact_genomic: 5, exact_cdna: 4, same_protein: 3, same_codon: 2, coordinate_only: 1 };
const FUNCTIONAL_SOURCES = {
  transactivation: {
    name: 'Kato et al. (2003)',
    method: 'Yeast reporter assays on eight p53-responsive promoters',
    units: '% of wild-type activity',
    url: 'https://pubmed.ncbi.nlm.nih.gov/12826609/'
  },
  cellular_dne_lof: {
    name: 'Giacomelli et al. (2018)',
    method: 'NCI TP53 Database classification of cellular growth-suppression and dominant-negative activity',
    url: 'https://pubmed.ncbi.nlm.nih.gov/30224644/'
  },
  dominant_negative: {
    name: 'NCI TP53 Database user manual',
    method: 'Dominant-negative transactivation classification curated across studies',
    url: 'https://tp53.cancer.gov/help'
  },
  structure_prediction: {
    name: 'Mathe et al. (2006)',
    method: 'Delaunay tessellation model of 3D structure, trained on the Kato transactivation data',
    url: 'https://tp53.cancer.gov/dt_method'
  }
};

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

function normalizeAssembly(value) {
  const name = normalizeText(value);
  if (name === 'grch37' || name === 'hg19') return 'GRCh37';
  if (name === 'grch38' || name === 'hg38') return 'GRCh38';
  return '';
}

function genomicAllele(assembly, pos, ref, alt) {
  if (!assembly || !Number.isSafeInteger(pos) || pos < 1 || ref === alt
    || !/^[ACGT]+$/.test(ref) || !/^[ACGT]+$/.test(alt)) return null;
  // Remove VCF padding from substitutions. MutationView frequently omits the
  // reference bases of indels, so those require cDNA or stay contextual.
  while (ref.length > 1 && alt.length > 1 && ref[0] === alt[0]) {
    ref = ref.slice(1);
    alt = alt.slice(1);
    pos++;
  }
  while (ref.length > 1 && alt.length > 1 && ref.at(-1) === alt.at(-1)) {
    ref = ref.slice(0, -1);
    alt = alt.slice(0, -1);
  }
  return { assembly, chrom: '17', pos, ref, alt };
}

function genomicHgvsSite(value, assemblyHint, chromHint = '') {
  const match = String(value || '').trim().match(/^(?:(NC_000017\.(?:10|11)|(?:chr)?17):)?g\.(\d+)(?:_(\d+))?(.*)$/i);
  if (!match || (!match[1] && !chromHint)) return null;
  const [, prefix, position, end, change] = match;
  const accessionAssembly = Object.keys(GENOMIC_COLUMNS).find(build => GENOMIC_COLUMNS[build].accession === prefix?.toUpperCase());
  if (accessionAssembly && assemblyHint && assemblyHint !== accessionAssembly) return null;
  const assembly = accessionAssembly || assemblyHint;
  if (!assembly) return null;
  const pos = Number(position);
  const substitution = !end && change.match(/^([ACGT])>([ACGT])$/i);
  const allele = substitution ? genomicAllele(assembly, pos, substitution[1].toUpperCase(), substitution[2].toUpperCase()) : null;
  return { assembly, pos, allele };
}

function queryGenomicSite(query) {
  const { genomic_allele: tuple, genomic } = query;
  const assembly = normalizeAssembly(tuple?.assembly);
  if (assembly && /^(?:chr)?17$/i.test(String(tuple.chrom))) {
    const pos = Number(tuple.pos);
    if (Number.isSafeInteger(pos) && pos > 0) {
      return { assembly, pos, allele: genomicAllele(assembly, pos, String(tuple.ref || '').toUpperCase(), String(tuple.alt || '').toUpperCase()) };
    }
  }
  return genomicHgvsSite(genomic, normalizeAssembly(query.assembly));
}

function sameGenomicAllele(a, b) {
  return a && b && a.assembly === b.assembly && a.pos === b.pos && a.ref === b.ref && a.alt === b.alt;
}

function parseCdna(value) {
  const match = String(value || '').trim().match(/^((?:NM_\d+|ENST\d+)\.\d+):(c\.[*-]?\d+(?:[+-]\d+)?(?:_[*-]?\d+(?:[+-]\d+)?)?(?:[ACGT]+>[ACGT]+|delins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|inv))$/i);
  return match ? { transcript: match[1].toUpperCase(), change: match[2].toUpperCase() } : null;
}

function cdnaMatchBasis(a, b) {
  if (!a || !b || a.change !== b.change) return '';
  if (a.transcript === b.transcript) return 'same_versioned_transcript';
  // These two RefSeq records have identical 1182-base CDS sequences. This
  // equivalence is limited to coding SNVs because UTRs and introns can differ.
  const pair = [a.transcript, b.transcript].sort().join(',');
  const codingSnv = a.change.match(/^C\.(\d+)[ACGT]>[ACGT]$/);
  return pair === 'NM_000546.5,NM_000546.6' && codingSnv && Number(codingSnv[1]) >= 1 && Number(codingSnv[1]) <= 1182
    ? 'identical_refseq_cds_for_coding_snv' : '';
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

function buildMatches(rows, query) {
  const pickValue = buildColumnResolver(rows);
  const proteinNeedle = toProteinCompact(query.protein);
  const cdnaNeedle = parseCdna(query.cdna);
  const site = queryGenomicSite(query);
  const columns = GENOMIC_COLUMNS[site?.assembly];

  const matches = [];
  for (const row of rows) {
    // Core identification fields
    const prot = pickValue(row, ['protdescription', 'aachangeinhuman', 'aachange']);
    const cd = pickValue(row, ['cdescription', 'gdescription']);
    const cdnaHgvs = pickValue(row, ['hgvsnm000546']);
    const hg19 = pickValue(row, ['hg19chr17coordinates']);
    const hg38 = pickValue(row, ['hg38chr17coordinates']);

    const protCompact = toProteinCompact(prot);
    const rowSite = columns ? genomicHgvsSite(pickValue(row, [columns.hgvs, columns.description]), site.assembly, '17') : null;
    const exactGenomic = sameGenomicAllele(site?.allele, rowSite?.allele);
    const genomicConflict = site?.allele && rowSite?.allele && !exactGenomic;
    const cdnaBasis = genomicConflict ? '' : cdnaMatchBasis(cdnaNeedle, parseCdna(cdnaHgvs || cd));
    const sameProtein = proteinNeedle && proteinNeedle !== '?' && proteinNeedle === protCompact;
    const coordinateMatch = columns && (rowSite?.pos === site.pos || coordinateHasPosition(pickValue(row, [columns.coordinate]), site.pos));
    const matchType = exactGenomic ? 'exact_genomic'
      : cdnaBasis ? 'exact_cdna'
      : sameProtein ? 'same_protein'
      : sameProteinCodon(protCompact, proteinNeedle) ? 'same_codon'
      : coordinateMatch ? 'coordinate_only'
      : '';
    if (!matchType) continue;

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
    const spliceAiDsAl = pickValue(row, ['spliceaidsal']);
    const spliceAiDsAg = pickValue(row, ['spliceaidsag']);
    const spliceAiDsDg = pickValue(row, ['spliceaidsdg']);
    const spliceAiDsDl = pickValue(row, ['spliceaidsdl']);

    matches.push({
      match_type: matchType,
      match_basis: cdnaBasis || null,
      matched_genomic_allele: exactGenomic ? rowSite.allele : null,
      same_codon_match: matchType === 'same_codon',
      mut_id: mutId || '',
      protein: prot || '',
      cdna_or_genomic: cd || '',
      transcript_hgvs: cdnaHgvs || '',
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
  matches.sort((a, b) => MATCH_RANK[b.match_type] - MATCH_RANK[a.match_type]);
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
      usage: 'POST JSON: { gene, protein, cdna (versioned HGVS), genomic (HGVS), assembly, genomic_allele: { assembly, chrom, pos, ref, alt } }'
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
    const { rows, datasetUrl, debug, staleSince } = await ensureDatasetLoaded();
    const matches = buildMatches(rows, body);
    const top = matches.slice(0, MAX_RESPONSE_MATCHES);
    const genomicMatches = matches.filter(m => m.match_type === 'exact_genomic');
    const cdnaMatches = matches.filter(m => m.match_type === 'exact_cdna');
    const candidates = genomicMatches.length ? genomicMatches : cdnaMatches;
    const best = candidates.length === 1 ? candidates[0] : null;
    const status = best ? 'exact' : candidates.length > 1 ? 'ambiguous' : top.length ? 'related_only' : 'no_match';
    const notes = {
      exact: best?.match_type === 'exact_genomic'
        ? 'Exact genomic allele matched by assembly, chromosome, position, REF and ALT.'
        : 'Exact cDNA change on a compatible versioned transcript.',
      ambiguous: 'Multiple TP53 dataset rows match the allele. An allele-specific summary could not be selected.',
      related_only: 'TP53 dataset rows provide related protein, codon or coordinate evidence. Exact allele identity could not be established.',
      no_match: 'No TP53 dataset row matched the supplied variant.'
    };

    // Only a uniquely selected exact record supplies the allele-specific summary.
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
      match_count: candidates.length,
      related_codon_count: matches.filter(m => m.same_codon_match).length,
      omitted_match_count: matches.length - top.length,
      selection: { status, match_type: candidates[0]?.match_type ?? null, candidate_count: candidates.length },
      functional_sources: FUNCTIONAL_SOURCES,
      assay_comparison: tp53AssayComparison(toProteinCompact(best ? best.protein : body.protein), best?.transcript_hgvs),
      dataset_stale: Boolean(staleSince),
      pathogenicity: pathogenicitySummary,
      classification: best?.exon_intron || '',
      prevalence: best?.somatic_count ? `somatic count: ${best.somatic_count}` : '',
      note: notes[status],
      matches: top
    };
    if (debugEnabled) {
      const debugPickValue = buildColumnResolver(rows);
      responsePayload.debug = {
        query: {
          protein_input: body.protein || '',
          protein_normalized: toProteinCompact(body.protein || ''),
          cdna_input: body.cdna || '',
          genomic_input: body.genomic || '',
          genomic_allele: body.genomic_allele || null
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

// Named exports for the test suite (tests/tp53-parse.test.js) — the Vercel
// runtime reads only the default export above.
export { parseCsv, proteinCodonKey, sameProteinCodon, coordinateHasPosition, extractDatasetUrls };
