#!/usr/bin/env python3
"""Generate docs/example-output.html from verification1 tier output files.
Run from the repo root: python docs/gen_example_output.py
"""
import csv
import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T3_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier3.tsv')
T2_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier2.tsv')
T1_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier1.maf')
OUT     = os.path.join(BASE, 'docs/example-output.html')


def load_tier(path):
    """Return (headers, data_rows), skipping any leading comment lines and the metadata row."""
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = list(reader)
    start = 0
    while start < len(rows) and rows[start] and rows[start][0].startswith('#'):
        start += 1
    headers = rows[start]
    data = rows[start + 2:]   # skip the source-metadata row
    return headers, data


t3h, t3d = load_tier(T3_PATH)
t2h, t2d = load_tier(T2_PATH)
t1h, t1d = load_tier(T1_PATH)

t3_set   = set(t3h)
t2_vis   = json.dumps([i for i, h in enumerate(t2h) if h in t3_set])
# For T1 (MAF format) some column names differ from T3; include common MAF equivalents
t1_extra = {
    'Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'Chromosome',
    'Start_Position', 'End_Position', 'NCBI_Build', 'Variant_Classification',
    'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'FILTER',
    'Tumor_Sample_Barcode', 'NM_Transcript',
}
t1_vis   = json.dumps([i for i, h in enumerate(t1h) if h in t3_set or h in t1_extra])

# ── Data as JS ─────────────────────────────────────────────────────────────────
T3H_JS   = json.dumps(t3h)
T3D_JS   = json.dumps(t3d)
T2H_JS   = json.dumps(t2h)
T2D_JS   = json.dumps(t2d)
T1H_JS   = json.dumps(t1h)
T1D_JS   = json.dumps(t1d)

# ── HTML template (placeholders replaced below) ─────────────────────────────────
TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>SMART — Example Output | Verification 1</title>

  <!-- DataTables 2.x (no jQuery) -->
  <link rel="stylesheet" href="https://cdn.datatables.net/2.0.8/css/dataTables.dataTables.min.css">

  <style>
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      font-size: 14px;
      background: #f8f9fa;
      color: #333;
    }

    /* ── Header ── */
    header {
      background: #2c3e50;
      color: #fff;
      padding: 16px 28px;
    }
    header h1 { font-size: 1.3em; margin-bottom: 4px; }
    header .sub { font-size: 0.82em; color: #adc8da; margin-top: 2px; }
    header a { color: #adc8da; text-decoration: none; }
    header a:hover { color: #fff; }

    /* ── Layout ── */
    .container { padding: 20px 28px 10px; }

    .info-box {
      background: #e8f4fd;
      border-left: 4px solid #3498db;
      padding: 10px 16px;
      border-radius: 0 5px 5px 0;
      margin-bottom: 20px;
      font-size: 0.88em;
      line-height: 1.6;
    }
    .info-box a { color: #1a6ea8; }

    /* ── Tabs ── */
    .tabs {
      display: flex;
      border-bottom: 2px solid #ddd;
      margin-bottom: 18px;
    }
    .tab-btn {
      padding: 9px 22px;
      border: none;
      background: none;
      cursor: pointer;
      font-size: 0.92em;
      color: #666;
      border-bottom: 3px solid transparent;
      margin-bottom: -2px;
      transition: background 0.12s;
    }
    .tab-btn:hover { color: #2c3e50; background: #f0f0f0; }
    .tab-btn.active {
      color: #2c3e50;
      font-weight: 600;
      border-bottom-color: #2c3e50;
    }
    .badge {
      display: inline-block;
      background: #e0e0e0;
      color: #555;
      border-radius: 10px;
      padding: 0 7px;
      font-size: 0.75em;
      margin-left: 5px;
    }
    .tab-btn.active .badge { background: #2c3e50; color: #fff; }

    /* ── Panels ── */
    .panel { display: none; }
    .panel.active { display: block; }
    .panel-desc {
      margin-bottom: 12px;
      font-size: 0.88em;
      color: #555;
      line-height: 1.6;
    }
    .panel-desc code {
      background: #f0f0f0;
      padding: 1px 5px;
      border-radius: 3px;
      font-size: 0.92em;
    }
    .controls {
      display: flex;
      gap: 10px;
      align-items: center;
      margin-bottom: 10px;
      flex-wrap: wrap;
    }
    .btn {
      padding: 5px 14px;
      border-radius: 4px;
      border: 1px solid #bbb;
      background: #fff;
      cursor: pointer;
      font-size: 0.84em;
      color: #333;
      white-space: nowrap;
      transition: background 0.12s;
    }
    .btn:hover { background: #f0f0f0; }
    .btn.on { background: #2c3e50; color: #fff; border-color: #2c3e50; }
    .hint { font-size: 0.8em; color: #888; }

    /* ── DataTables overrides ── */
    div.dt-container { width: 100%; }
    div.dt-container .dt-search { margin-bottom: 8px; }
    div.dt-container .dt-search input {
      padding: 5px 10px;
      border: 1px solid #ccc;
      border-radius: 4px;
      font-size: 0.88em;
      min-width: 240px;
    }
    table.dataTable { font-size: 0.79em; }
    table.dataTable thead th {
      background: #2c3e50 !important;
      color: #fff !important;
      white-space: nowrap;
      padding: 6px 10px;
      cursor: pointer;
    }
    table.dataTable thead th:hover { background: #3d5368 !important; }
    table.dataTable tbody td {
      padding: 4px 8px;
      vertical-align: top;
      max-width: 280px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }
    table.dataTable tbody tr { background: #fff; }
    table.dataTable tbody tr:nth-child(even) { background: #f8f9fa; }
    table.dataTable tbody tr:hover td { background: #ebf3fb !important; }

    /* ── Colour badges ── */
    .lvl {
      display: inline-block;
      padding: 1px 6px;
      border-radius: 3px;
      font-size: 0.8em;
      font-weight: 700;
      white-space: nowrap;
      letter-spacing: 0.01em;
    }
    .L1   { background: #1a9641; color: #fff; }
    .L2   { background: #0571b0; color: #fff; }
    .L3A  { background: #e66101; color: #fff; }
    .L3B  { background: #fdb863; color: #333; }
    .L4   { background: #92c5de; color: #333; }
    .LR1  { background: #ca0020; color: #fff; }
    .LR2  { background: #f4a582; color: #333; }
    .LDx  { background: #7b2d8b; color: #fff; }
    .LPx  { background: #984ea3; color: #fff; }
    .LFda { background: #5e2d8b; color: #fff; }

    .onc-oc  { display:inline-block; background:#d4edda; color:#155724; border-radius:3px; padding:1px 5px; font-size:0.8em; }
    .onc-lo  { display:inline-block; background:#c3e6cb; color:#0c4128; border-radius:3px; padding:1px 5px; font-size:0.8em; }
    .onc-inc { display:inline-block; background:#fff3cd; color:#856404; border-radius:3px; padding:1px 5px; font-size:0.8em; }
    .onc-vus { display:inline-block; background:#e8e8e8; color:#444;    border-radius:3px; padding:1px 5px; font-size:0.8em; }
    .onc-neu { display:inline-block; background:#f8d7da; color:#721c24; border-radius:3px; padding:1px 5px; font-size:0.8em; }

    .hs   { color: #155724; font-weight: 600; }
    .path { color: #721c24; font-weight: 500; }
    .ben  { color: #155724; font-weight: 500; }

    /* ── Footer ── */
    footer {
      text-align: center;
      padding: 20px;
      color: #999;
      font-size: 0.82em;
      border-top: 1px solid #e8e8e8;
      margin-top: 30px;
    }
    footer a { color: #3498db; }
  </style>
</head>
<body>

<header>
  <h1>SMART &mdash; Example Output</h1>
  <p class="sub">
    Verification 1 &nbsp;&middot;&nbsp; 22 annotated variants &nbsp;&middot;&nbsp;
    1 tumour sample (TUMOR) &nbsp;&middot;&nbsp; SMART v0.2.0 &nbsp;&middot;&nbsp;
    OncoKB v7.0 &nbsp;&middot;&nbsp; VEP 114.2 &nbsp;&middot;&nbsp; GRCh38
    &nbsp;|&nbsp;
    <a href="./">&#8592; Field reference</a>
    &nbsp;&middot;&nbsp;
    <a href="https://github.com/Manuel-DominguezCBG/SMART">GitHub &#8599;</a>
  </p>
</header>

<div class="container">

  <div class="info-box">
    <strong>About this dataset:</strong> Verification 1 contains 18 curated synthetic variants
    (14 SNV/indel &plus; 4 CNA) covering key clinical scenarios across
    NRAS, IDH1, PIK3CA, EGFR, BRAF, GNAQ, PTEN, KRAS, BRCA2, DICER1, TP53,
    ERBB2, MET, CDKN2A (two isoforms), and CDK4.
    The CDKN2A deletion appears as two rows &mdash; one for each preferred isoform
    (p16/INK4a &amp; p14ARF) &mdash; demonstrating SMART&rsquo;s multi-transcript output.
    All three files below come from a single pipeline run;
    column counts differ because each tier applies a different field filter
    (see the <a href="./">field reference</a> for details).
    Click any column header to sort. Use the filter box to search across all fields.
  </div>

  <!-- Tabs -->
  <div class="tabs">
    <button class="tab-btn active" onclick="showTab('t3', this)">
      Tier 3 &mdash; Clinical
      <span class="badge">__T3_NCOLS__ cols</span>
    </button>
    <button class="tab-btn" onclick="showTab('t2', this)">
      Tier 2 &mdash; Bioinformatics
      <span class="badge">__T2_NCOLS__ cols</span>
    </button>
    <button class="tab-btn" onclick="showTab('t1', this)">
      Tier 1 &mdash; Full MAF
      <span class="badge">__T1_NCOLS__ cols</span>
    </button>
  </div>

  <!-- ── Tier 3 ── -->
  <div id="panel-t3" class="panel active">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier3.tsv</code></strong> &nbsp;&mdash;&nbsp;
      Optimised for clinical scientists. All __T3_NCOLS__ columns are shown:
      variant identity, functional consequence, population frequency, pathogenicity scores
      (REVEL, SpliceAI, LOEUF), ClinVar, CancerHotspots, and top-line OncoKB actionability
      (oncogenicity, mutation effect, therapeutic/diagnostic levels).
      Two-row header in the file (field names + source metadata); only field names shown here.
    </p>
    <table id="tbl-t3" class="display" style="width:100%"></table>
  </div>

  <!-- ── Tier 2 ── -->
  <div id="panel-t2" class="panel">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier2.tsv</code></strong> &nbsp;&mdash;&nbsp;
      Optimised for bioinformaticians. Extends Tier 3 with full gnomAD population stratification,
      all SpliceAI delta scores, complete structured CIViC fields, and every OncoKB JSON expansion
      (<code>ONCOKB_TX_*</code>, <code>ONCOKB_DIAG_*</code>, <code>ONCOKB_PROG_*</code>)
      unpacked into indexed columns with drug names, tumour types, and evidence citations.
    </p>
    <div class="controls">
      <button id="btn-t2" class="btn" onclick="toggleExt('t2')">
        Show all __T2_NCOLS__ columns
      </button>
      <span class="hint">Default: core clinical columns (__T3_NCOLS__ of __T2_NCOLS__). Toggle to reveal extended fields.</span>
    </div>
    <table id="tbl-t2" class="display" style="width:100%"></table>
  </div>

  <!-- ── Tier 1 ── -->
  <div id="panel-t1" class="panel">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier1.maf</code></strong> &nbsp;&mdash;&nbsp;
      Standard MAF format for downstream tools (cBioPortal, R/Python pipelines, oncoPrint generators).
      Carries all non-dropped fields including every expanded OncoKB treatment, diagnostic, and
      prognostic entry, plus raw VCF FORMAT fields.
      Column count scales with OncoKB evidence breadth; this run produced __T1_NCOLS__ columns.
      The first line of the file is <code>#SMART_VERSION 0.2.0</code>.
    </p>
    <div class="controls">
      <button id="btn-t1" class="btn" onclick="toggleExt('t1')">
        Show all __T1_NCOLS__ columns
      </button>
      <span class="hint">Default: core columns only. Toggle to reveal all __T1_NCOLS__ MAF fields.</span>
    </div>
    <table id="tbl-t1" class="display" style="width:100%"></table>
  </div>

</div><!-- /.container -->

<footer>
  SMART &middot; Somatic Mutation Annotation and Reporting Tool &middot; University Hospital Southampton<br>
  <a href="./">Field reference</a> &nbsp;&middot;&nbsp;
  <a href="https://github.com/Manuel-DominguezCBG/SMART">GitHub</a> &nbsp;&middot;&nbsp;
  <a href="https://github.com/Manuel-DominguezCBG/SMART/blob/main/README.md">Documentation</a>
</footer>

<!-- ── Embedded data ─────────────────────────────────────────────────────────── -->
<script>
var T3H = __T3_HEADERS__;
var T3D = __T3_DATA__;
var T2H = __T2_HEADERS__;
var T2D = __T2_DATA__;
var T1H = __T1_HEADERS__;
var T1D = __T1_DATA__;
var T2_VIS = new Set(__T2_VIS__);
var T1_VIS = new Set(__T1_VIS__);
</script>

<!-- ── DataTables (vanilla JS, no jQuery) ────────────────────────────────────── -->
<script src="https://cdn.datatables.net/2.0.8/js/dataTables.min.js"></script>

<script>
// ── Level / oncogenicity badge renderer ─────────────────────────────────────────
var LEVEL_CLASS = {
  'LEVEL_1': 'L1', 'LEVEL_2': 'L2', 'LEVEL_3A': 'L3A', 'LEVEL_3B': 'L3B',
  'LEVEL_4': 'L4', 'LEVEL_R1': 'LR1', 'LEVEL_R2': 'LR2'
};
function lvlBadge(v) {
  var cls = LEVEL_CLASS[v];
  if (cls) return '<span class="lvl ' + cls + '">' + v + '</span>';
  if (/^LEVEL_Dx/.test(v)) return '<span class="lvl LDx">' + v + '</span>';
  if (/^LEVEL_Px/.test(v)) return '<span class="lvl LPx">' + v + '</span>';
  if (/^LEVEL_Fda/.test(v)) return '<span class="lvl LFda">' + v + '</span>';
  return null;
}
var ONC_MAP = {
  'Oncogenic': 'onc-oc', 'Likely Oncogenic': 'onc-lo',
  'Inconclusive': 'onc-inc', 'VUS': 'onc-vus',
  'Likely Neutral': 'onc-neu', 'Neutral': 'onc-neu',
  'Unknown': 'onc-vus'
};

function makeRenderer(headers) {
  return function(data, type, row, meta) {
    if (type !== 'display') return data;
    if (!data) return '';
    var col = headers[meta.col];

    // OncoKB levels
    var b = lvlBadge(data);
    if (b) return b;

    // Comma-separated list of levels (e.g. LEVEL_2 treatment fields)
    if (data.indexOf('LEVEL_') === -1 && data.indexOf(',') !== -1) {
      // plain drug list — just show as-is (will be truncated by CSS)
    }

    // Oncogenicity
    var oc = ONC_MAP[data];
    if (oc) return '<span class="' + oc + '">' + data + '</span>';

    // Hotspot flag
    if ((col === 'ONCOKB_HOTSPOT' || col === 'CancerHotspots_HOTSPOT') &&
        (data === 'True' || data === '1')) {
      return '<span class="hs">&#10003; Hotspot</span>';
    }

    // ClinVar significance
    if (col === 'ClinVar_CLNSIG' || col === 'CLNSIG') {
      if (/^Pathogenic/i.test(data)) return '<span class="path">' + data + '</span>';
      if (/^Likely.pathogenic/i.test(data)) return '<span class="path">' + data + '</span>';
      if (/^Benign/i.test(data) || /^Likely.benign/i.test(data))
        return '<span class="ben">' + data + '</span>';
    }

    return data;
  };
}

// ── Table initialisation ──────────────────────────────────────────────────────
var DT = {};
var extOn = { t2: false, t1: false };

function initTable(tblId, headers, data, visSet) {
  var render = makeRenderer(headers);
  var cols = headers.map(function(h) {
    return { title: h, defaultContent: '', render: render };
  });
  var dt;
  try {
    dt = new DataTable('#' + tblId, {
      data: data,
      columns: cols,
      scrollX: true,
      paging: false,
      order: [],
      searching: true,
      info: false
    });
  } catch(e) {
    document.getElementById(tblId).insertAdjacentHTML(
      'afterend',
      '<p style="color:red;padding:8px">Table error: ' + e.message + '</p>'
    );
    return null;
  }
  // Hide extended columns after init using the batch API.
  if (visSet) {
    var extIdx = [];
    for (var i = 0; i < headers.length; i++) {
      if (!visSet.has(i)) extIdx.push(i);
    }
    if (extIdx.length) dt.columns(extIdx).visible(false, false);
    dt.columns.adjust().draw(false);
  }
  return dt;
}

// ── Tab switching (lazy-init T2/T1 on first view) ─────────────────────────────
function showTab(id, btn) {
  document.querySelectorAll('.panel').forEach(function(p) { p.classList.remove('active'); });
  document.querySelectorAll('.tab-btn').forEach(function(b) { b.classList.remove('active'); });
  document.getElementById('panel-' + id).classList.add('active');
  btn.classList.add('active');

  if (!DT[id]) {
    // Small delay so the panel is painted before DataTables measures widths.
    setTimeout(function() {
      if (id === 't2') DT.t2 = initTable('tbl-t2', T2H, T2D, T2_VIS);
      if (id === 't1') DT.t1 = initTable('tbl-t1', T1H, T1D, T1_VIS);
      if (DT[id]) DT[id].columns.adjust().draw(false);
    }, 60);
  } else {
    DT[id].columns.adjust().draw(false);
  }
}

// ── Column visibility toggle for T2 / T1 ─────────────────────────────────────
function toggleExt(tier) {
  if (!DT[tier]) return;
  extOn[tier] = !extOn[tier];
  var vis     = (tier === 't2') ? T2_VIS : T1_VIS;
  var headers = (tier === 't2') ? T2H    : T1H;
  var extIdx  = [];
  for (var i = 0; i < headers.length; i++) {
    if (!vis.has(i)) extIdx.push(i);
  }
  if (extIdx.length) DT[tier].columns(extIdx).visible(extOn[tier], false);
  DT[tier].columns.adjust().draw(false);
  var btn = document.getElementById('btn-' + tier);
  btn.classList.toggle('on', extOn[tier]);
  btn.textContent = extOn[tier]
    ? 'Show core columns only'
    : 'Show all ' + headers.length + ' columns';
}

// ── Boot: Tier 3 is always visible on load ────────────────────────────────────
DT.t3 = initTable('tbl-t3', T3H, T3D, null);
</script>

</body>
</html>
"""

html = TEMPLATE \
    .replace('__T3_HEADERS__', T3H_JS) \
    .replace('__T3_DATA__',    T3D_JS) \
    .replace('__T2_HEADERS__', T2H_JS) \
    .replace('__T2_DATA__',    T2D_JS) \
    .replace('__T1_HEADERS__', T1H_JS) \
    .replace('__T1_DATA__',    T1D_JS) \
    .replace('__T2_VIS__',     t2_vis) \
    .replace('__T1_VIS__',     t1_vis) \
    .replace('__T3_NCOLS__',   str(len(t3h))) \
    .replace('__T2_NCOLS__',   str(len(t2h))) \
    .replace('__T1_NCOLS__',   str(len(t1h)))

with open(OUT, 'w', encoding='utf-8') as f:
    f.write(html)

size_kb = os.path.getsize(OUT) // 1024
print(f'Written: {OUT}')
print(f'T3: {len(t3h)} cols × {len(t3d)} rows')
print(f'T2: {len(t2h)} cols × {len(t2d)} rows')
print(f'T1: {len(t1h)} cols × {len(t1d)} rows')
print(f'File size: {size_kb} KB')
