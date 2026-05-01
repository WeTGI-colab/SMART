#!/usr/bin/env python3
"""
Generate docs/example-output.html from verification1 tier output files.
Pure static HTML approach — no CDN, no DataTables, no external dependencies.
Run from the repo root:  python docs/gen_example_output.py
"""
import csv
import html as HL
import os

BASE   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T3_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier3.tsv')
T2_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier2.tsv')
T1_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier1.maf')
OUT     = os.path.join(BASE, 'docs/example-output.html')

# ── Data loading ───────────────────────────────────────────────────────────────
def load_tier(path):
    with open(path, newline='', encoding='utf-8') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    start = 0
    while start < len(rows) and rows[start] and rows[start][0].startswith('#'):
        start += 1
    return rows[start], rows[start + 2:]   # headers, data (skip metadata row)

t3h, t3d = load_tier(T3_PATH)
t2h, t2d = load_tier(T2_PATH)
t1h, t1d = load_tier(T1_PATH)

# ── Column ordering: core columns first, extended after ────────────────────────
t3_set   = set(t3h)
t3_order = {h: i for i, h in enumerate(t3h)}   # name → position in T3

T1_EXTRA = {
    'Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'Chromosome',
    'Start_Position', 'End_Position', 'NCBI_Build', 'Variant_Classification',
    'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'FILTER', 'NM_Transcript',
}

def reorder(headers, data, core_fn):
    """Return (new_headers, new_data, core_count).
    core_fn(i, h) → True means column is 'core' (always visible)."""
    core_idx = [i for i, h in enumerate(headers) if core_fn(i, h)]
    ext_idx  = [i for i in range(len(headers)) if not core_fn(i, h := headers[i])]
    # Sort T2 core cols in T3 column order so they match T3 exactly
    core_idx.sort(key=lambda i: t3_order.get(headers[i], len(t3h) + i))
    new_order   = core_idx + ext_idx
    new_headers = [headers[i] for i in new_order]
    new_data    = [[row[i] if i < len(row) else '' for i in new_order] for row in data]
    return new_headers, new_data, len(core_idx)

t2h_r, t2d_r, t2_core = reorder(t2h, t2d, lambda i, h: h in t3_set)
t1h_r, t1d_r, t1_core = reorder(t1h, t1d, lambda i, h: h in t3_set or h in T1_EXTRA)

# ── Cell renderer (colour badges) ──────────────────────────────────────────────
_LEVELS = {
    'LEVEL_1': 'L1', 'LEVEL_2': 'L2', 'LEVEL_3A': 'L3A', 'LEVEL_3B': 'L3B',
    'LEVEL_4': 'L4', 'LEVEL_R1': 'LR1', 'LEVEL_R2': 'LR2',
}
_ONC = {
    'Oncogenic': 'onc-oc', 'Likely Oncogenic': 'onc-lo', 'Inconclusive': 'onc-inc',
    'VUS': 'onc-vus', 'Likely Neutral': 'onc-neu', 'Neutral': 'onc-neu',
    'Unknown': 'onc-vus',
}

def render_cell(col, val):
    if not val:
        return ''
    lv = _LEVELS.get(val)
    if lv:
        return f'<span class="lvl {lv}">{HL.escape(val)}</span>'
    if val.startswith('LEVEL_Dx'):
        return f'<span class="lvl LDx">{HL.escape(val)}</span>'
    if val.startswith('LEVEL_Px'):
        return f'<span class="lvl LPx">{HL.escape(val)}</span>'
    if val.startswith('LEVEL_Fda'):
        return f'<span class="lvl LFda">{HL.escape(val)}</span>'
    oc = _ONC.get(val)
    if oc:
        return f'<span class="{oc}">{HL.escape(val)}</span>'
    if col in ('ONCOKB_HOTSPOT', 'CancerHotspots_HOTSPOT') and val in ('True', '1'):
        return '<span class="hs">&#10003; Hotspot</span>'
    if col in ('ClinVar_CLNSIG', 'CLNSIG'):
        vl = val.lower()
        if 'pathogenic' in vl and 'likely' not in vl:
            return f'<span class="path">{HL.escape(val)}</span>'
        if 'pathogenic' in vl:
            return f'<span class="path">{HL.escape(val)}</span>'
        if 'benign' in vl:
            return f'<span class="ben">{HL.escape(val)}</span>'
    return HL.escape(val)

# ── HTML table builder ─────────────────────────────────────────────────────────
def build_table(tbl_id, headers, data):
    p = []
    p.append(f'<div class="tbl-wrap"><table class="dt" id="{tbl_id}">')
    # thead
    p.append('<thead><tr>')
    for i, h in enumerate(headers):
        esc = HL.escape(h)
        p.append(f'<th onclick="sort(\'{tbl_id}\',{i})" title="{esc}">{esc} <span class="si">&#8597;</span></th>')
    p.append('</tr></thead><tbody>')
    # tbody
    for row in data:
        vals = [row[i] if i < len(row) else '' for i in range(len(headers))]
        # pre-compute search string (all values joined, lowercase)
        search = HL.escape('|'.join(v.lower() for v in vals))
        p.append(f'<tr data-s="{search}">')
        for i, val in enumerate(vals):
            cell = render_cell(headers[i], val)
            title = f' title="{HL.escape(val)}"' if len(val) > 50 else ''
            p.append(f'<td{title}>{cell}</td>')
        p.append('</tr>')
    p.append('</tbody></table></div>')
    return ''.join(p)

tbl_t3 = build_table('tbl-t3', t3h,   t3d)
tbl_t2 = build_table('tbl-t2', t2h_r, t2d_r)
tbl_t1 = build_table('tbl-t1', t1h_r, t1d_r)

# ── HTML page ──────────────────────────────────────────────────────────────────
PAGE = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>SMART — Example Output | Verification 1</title>
  <style>
    *,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
    body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;font-size:14px;background:#f8f9fa;color:#333}}

    header{{background:#2c3e50;color:#fff;padding:16px 28px}}
    header h1{{font-size:1.3em;margin-bottom:4px}}
    header .sub{{font-size:.82em;color:#adc8da}}
    header a{{color:#adc8da;text-decoration:none}}
    header a:hover{{color:#fff}}

    .container{{padding:20px 28px 10px}}

    .info-box{{background:#e8f4fd;border-left:4px solid #3498db;padding:10px 16px;
               border-radius:0 5px 5px 0;margin-bottom:20px;font-size:.88em;line-height:1.6}}
    .info-box a{{color:#1a6ea8}}

    .tabs{{display:flex;border-bottom:2px solid #ddd;margin-bottom:18px}}
    .tab-btn{{padding:9px 22px;border:none;background:none;cursor:pointer;font-size:.92em;
              color:#666;border-bottom:3px solid transparent;margin-bottom:-2px}}
    .tab-btn:hover{{color:#2c3e50;background:#f0f0f0}}
    .tab-btn.active{{color:#2c3e50;font-weight:600;border-bottom-color:#2c3e50}}
    .badge{{display:inline-block;background:#e0e0e0;color:#555;border-radius:10px;
            padding:0 7px;font-size:.75em;margin-left:5px}}
    .tab-btn.active .badge{{background:#2c3e50;color:#fff}}

    .panel{{display:none}}.panel.active{{display:block}}
    .panel-desc{{margin-bottom:12px;font-size:.88em;color:#555;line-height:1.6}}
    .panel-desc code{{background:#f0f0f0;padding:1px 5px;border-radius:3px;font-size:.92em}}

    .controls{{display:flex;gap:10px;align-items:center;margin-bottom:10px;flex-wrap:wrap}}
    .search-box{{padding:5px 10px;border:1px solid #ccc;border-radius:4px;
                 font-size:.88em;min-width:240px;flex:1}}
    .btn{{padding:5px 14px;border-radius:4px;border:1px solid #bbb;background:#fff;
          cursor:pointer;font-size:.84em;color:#333;white-space:nowrap}}
    .btn:hover{{background:#f0f0f0}}
    .btn.on{{background:#2c3e50;color:#fff;border-color:#2c3e50}}
    .hint{{font-size:.8em;color:#888}}

    .tbl-wrap{{overflow-x:auto;max-height:65vh;overflow-y:auto;border:1px solid #dee2e6;border-radius:4px}}
    table.dt{{border-collapse:collapse;font-size:.79em;white-space:nowrap}}
    table.dt thead th{{position:sticky;top:0;background:#2c3e50;color:#fff;
                       padding:6px 10px;cursor:pointer;user-select:none;border-right:1px solid #3d5368}}
    table.dt thead th:hover{{background:#3d5368}}
    table.dt thead .si{{opacity:.5;font-size:.75em}}
    table.dt tbody td{{padding:4px 8px;border-bottom:1px solid #eee;
                       max-width:280px;overflow:hidden;text-overflow:ellipsis}}
    table.dt tbody tr:nth-child(even){{background:#f8f9fa}}
    table.dt tbody tr:hover td{{background:#ebf3fb}}
    tr.hidden{{display:none}}

    /* Extended columns hidden by default; shown when table has .show-all */
    #tbl-t2:not(.show-all) th:nth-child(n+{t2_core + 1}),
    #tbl-t2:not(.show-all) td:nth-child(n+{t2_core + 1}){{display:none}}
    #tbl-t1:not(.show-all) th:nth-child(n+{t1_core + 1}),
    #tbl-t1:not(.show-all) td:nth-child(n+{t1_core + 1}){{display:none}}

    /* Badges */
    .lvl{{display:inline-block;padding:1px 6px;border-radius:3px;font-size:.8em;font-weight:700}}
    .L1{{background:#1a9641;color:#fff}} .L2{{background:#0571b0;color:#fff}}
    .L3A{{background:#e66101;color:#fff}} .L3B{{background:#fdb863;color:#333}}
    .L4{{background:#92c5de;color:#333}} .LR1{{background:#ca0020;color:#fff}}
    .LR2{{background:#f4a582;color:#333}} .LDx{{background:#7b2d8b;color:#fff}}
    .LPx{{background:#984ea3;color:#fff}} .LFda{{background:#5e2d8b;color:#fff}}
    .onc-oc{{background:#d4edda;color:#155724;border-radius:3px;padding:1px 5px;font-size:.8em}}
    .onc-lo{{background:#c3e6cb;color:#0c4128;border-radius:3px;padding:1px 5px;font-size:.8em}}
    .onc-inc{{background:#fff3cd;color:#856404;border-radius:3px;padding:1px 5px;font-size:.8em}}
    .onc-vus{{background:#e8e8e8;color:#444;border-radius:3px;padding:1px 5px;font-size:.8em}}
    .onc-neu{{background:#f8d7da;color:#721c24;border-radius:3px;padding:1px 5px;font-size:.8em}}
    .hs{{color:#155724;font-weight:600}}
    .path{{color:#721c24;font-weight:500}}
    .ben{{color:#155724;font-weight:500}}

    footer{{text-align:center;padding:20px;color:#999;font-size:.82em;
            border-top:1px solid #e8e8e8;margin-top:30px}}
    footer a{{color:#3498db}}
  </style>
</head>
<body>

<header>
  <h1>SMART &mdash; Example Output</h1>
  <p class="sub">
    Verification&nbsp;1 &nbsp;&middot;&nbsp; 22&nbsp;annotated variants &nbsp;&middot;&nbsp;
    1&nbsp;tumour sample (TUMOR) &nbsp;&middot;&nbsp; SMART&nbsp;v0.2.0 &nbsp;&middot;&nbsp;
    OncoKB&nbsp;v7.0 &nbsp;&middot;&nbsp; VEP&nbsp;114.2 &nbsp;&middot;&nbsp; GRCh38
    &nbsp;|&nbsp;
    <a href="./">&larr; Field reference</a> &nbsp;&middot;&nbsp;
    <a href="https://github.com/Manuel-DominguezCBG/SMART">GitHub &nearr;</a>
  </p>
</header>

<div class="container">

  <div class="info-box">
    <strong>About this dataset:</strong> Verification&nbsp;1 contains 18 curated synthetic variants
    (14&nbsp;SNV/indel &plus; 4&nbsp;CNA) covering NRAS, IDH1, PIK3CA, EGFR, BRAF, GNAQ, PTEN, KRAS,
    BRCA2, DICER1, TP53, ERBB2, MET, CDKN2A&nbsp;(&times;2&nbsp;isoforms), and CDK4.
    The CDKN2A deletion appears as two rows &mdash; one per isoform (p16/INK4a &amp; p14ARF) &mdash;
    demonstrating multi-transcript output.
    Click any column header to sort &uarr;&darr;. Use the search box to filter rows.
    See the <a href="./">field reference</a> for full column definitions.
  </div>

  <div class="tabs">
    <button class="tab-btn active" onclick="showTab('t3',this)">
      Tier&nbsp;3 &mdash; Clinical <span class="badge">{len(t3h)}&nbsp;cols</span>
    </button>
    <button class="tab-btn" onclick="showTab('t2',this)">
      Tier&nbsp;2 &mdash; Bioinformatics <span class="badge">{len(t2h)}&nbsp;cols</span>
    </button>
    <button class="tab-btn" onclick="showTab('t1',this)">
      Tier&nbsp;1 &mdash; Full&nbsp;MAF <span class="badge">{len(t1h)}&nbsp;cols</span>
    </button>
  </div>

  <!-- Tier 3 -->
  <div id="panel-t3" class="panel active">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier3.tsv</code></strong> &mdash;
      Optimised for clinical scientists. All {len(t3h)} columns shown: variant identity,
      functional consequence, population frequency, pathogenicity scores (REVEL, SpliceAI, LOEUF),
      ClinVar, CancerHotspots, and top-line OncoKB actionability.
    </p>
    <div class="controls">
      <input class="search-box" type="text" placeholder="Search all fields&hellip;"
             oninput="filterRows('tbl-t3',this.value)">
    </div>
    {tbl_t3}
  </div>

  <!-- Tier 2 -->
  <div id="panel-t2" class="panel">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier2.tsv</code></strong> &mdash;
      Optimised for bioinformaticians. Extends Tier&nbsp;3 with full gnomAD population
      stratification, all SpliceAI delta scores, complete structured CIViC fields, and every
      OncoKB JSON expansion (<code>ONCOKB_TX_*</code>, <code>ONCOKB_DIAG_*</code>,
      <code>ONCOKB_PROG_*</code>). Core columns (matching Tier&nbsp;3) are shown first.
    </p>
    <div class="controls">
      <input class="search-box" type="text" placeholder="Search all fields&hellip;"
             oninput="filterRows('tbl-t2',this.value)">
      <button id="btn-t2" class="btn" onclick="toggleExt('t2')">
        Show all {len(t2h)} columns
      </button>
      <span class="hint">Default: {t2_core} core columns. Toggle to reveal all {len(t2h)}.</span>
    </div>
    {tbl_t2}
  </div>

  <!-- Tier 1 -->
  <div id="panel-t1" class="panel">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier1.maf</code></strong> &mdash;
      Standard MAF format for downstream tools (cBioPortal, R/Python pipelines).
      Carries all non-dropped fields including every expanded OncoKB treatment, diagnostic,
      and prognostic entry. This run produced {len(t1h)} columns.
      Core columns (matching Tier&nbsp;3) are shown first.
    </p>
    <div class="controls">
      <input class="search-box" type="text" placeholder="Search all fields&hellip;"
             oninput="filterRows('tbl-t1',this.value)">
      <button id="btn-t1" class="btn" onclick="toggleExt('t1')">
        Show all {len(t1h)} columns
      </button>
      <span class="hint">Default: {t1_core} core columns. Toggle to reveal all {len(t1h)}.</span>
    </div>
    {tbl_t1}
  </div>

</div>

<footer>
  SMART &middot; Somatic Mutation Annotation and Reporting Tool &middot; University Hospital Southampton<br>
  <a href="./">Field reference</a> &nbsp;&middot;&nbsp;
  <a href="https://github.com/Manuel-DominguezCBG/SMART">GitHub</a> &nbsp;&middot;&nbsp;
  <a href="https://github.com/Manuel-DominguezCBG/SMART/blob/main/README.md">Documentation</a>
</footer>

<script>
var extOn = {{t2:false, t1:false}};

function showTab(id, btn) {{
  document.querySelectorAll('.panel').forEach(function(p){{p.classList.remove('active')}});
  document.querySelectorAll('.tab-btn').forEach(function(b){{b.classList.remove('active')}});
  document.getElementById('panel-'+id).classList.add('active');
  btn.classList.add('active');
}}

function filterRows(tblId, q) {{
  q = q.toLowerCase().trim();
  document.querySelectorAll('#'+tblId+' tbody tr').forEach(function(tr) {{
    tr.classList.toggle('hidden', q !== '' && !tr.dataset.s.includes(q));
  }});
}}

function toggleExt(tier) {{
  extOn[tier] = !extOn[tier];
  document.getElementById('tbl-'+tier).classList.toggle('show-all', extOn[tier]);
  var btn = document.getElementById('btn-'+tier);
  btn.classList.toggle('on', extOn[tier]);
  var nc   = tier==='t2' ? {len(t2h)}  : {len(t1h)};
  var nc3  = tier==='t2' ? {t2_core} : {t1_core};
  btn.textContent = extOn[tier] ? 'Show core only ('+nc3+')' : 'Show all '+nc+' columns';
}}

function sort(tblId, col) {{
  var tbl   = document.getElementById(tblId);
  var tbody = tbl.querySelector('tbody');
  var rows  = Array.from(tbody.querySelectorAll('tr'));
  var asc   = tbl.dataset.sc == col && tbl.dataset.sd === 'a';
  rows.sort(function(a,b) {{
    var av = a.cells[col] ? a.cells[col].textContent.trim() : '';
    var bv = b.cells[col] ? b.cells[col].textContent.trim() : '';
    var n  = parseFloat(av), m = parseFloat(bv);
    if (!isNaN(n) && !isNaN(m)) return asc ? m-n : n-m;
    return asc ? bv.localeCompare(av) : av.localeCompare(bv);
  }});
  rows.forEach(function(r){{tbody.appendChild(r)}});
  tbl.dataset.sc = col;
  tbl.dataset.sd = asc ? 'd' : 'a';
}}
</script>

</body>
</html>
"""

with open(OUT, 'w', encoding='utf-8') as f:
    f.write(PAGE)

size_kb = os.path.getsize(OUT) // 1024
print(f'Written: {OUT}')
print(f'T3: {len(t3h)} cols x {len(t3d)} rows')
print(f'T2: {len(t2h)} cols x {len(t2d)} rows  (core visible: {t2_core})')
print(f'T1: {len(t1h)} cols x {len(t1d)} rows  (core visible: {t1_core})')
print(f'File size: {size_kb} KB')
