#!/usr/bin/env python3
"""
Generate docs/example-output.html from verification1 tier output files.
Pure static HTML — no CDN, no DataTables, no external dependencies.
Run from the repo root:  python docs/gen_example_output.py
"""
import csv
import html as HL
import os

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T3_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier3.tsv')
T2_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier2.tsv')
T1_PATH = os.path.join(BASE, 'tests/verification1/output/output/Final_result_tier1.maf')
OUT     = os.path.join(BASE, 'docs/example-output.html')

# ── Data loading ───────────────────────────────────────────────────────────────
def load_tier(path):
    """Return (comment_lines, headers, metadata_row_or_None, data_rows).

    Detects the two-row header format (field names + 'desc | source | version')
    used by Tier 2 and Tier 3 TSV files.  Tier 1 MAF has a leading comment line
    and no metadata row, so metadata_row is None and data starts at the row
    immediately after the header.
    """
    with open(path, newline='', encoding='utf-8') as f:
        rows = list(csv.reader(f, delimiter='\t'))

    # Collect leading comment lines (e.g. '#SMART_VERSION 0.2.0')
    start, comments = 0, []
    while start < len(rows) and rows[start] and rows[start][0].startswith('#'):
        comments.append('\t'.join(rows[start]))
        start += 1

    headers  = rows[start]
    next_row = rows[start + 1] if start + 1 < len(rows) else []

    # The metadata row uses 'description | source | version' format.
    # A plain data value (e.g. 'GRCh38') will never contain ' | '.
    if next_row and ' | ' in next_row[0]:
        return comments, headers, next_row, rows[start + 2:]
    else:
        return comments, headers, None, rows[start + 1:]

t3_comments, t3h, t3_meta, t3d = load_tier(T3_PATH)
t2_comments, t2h, t2_meta, t2d = load_tier(T2_PATH)
t1_comments, t1h, t1_meta, t1d = load_tier(T1_PATH)

# ── Column ordering: core first, extended after ────────────────────────────────
t3_set   = set(t3h)
t3_order = {h: i for i, h in enumerate(t3h)}   # name → T3 position

T1_EXTRA = {
    'Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'Chromosome',
    'Start_Position', 'End_Position', 'NCBI_Build', 'Variant_Classification',
    'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'FILTER', 'NM_Transcript',
}

def reorder(headers, meta, data, core_fn):
    """Return (new_headers, new_meta, new_data, core_count).
    core_fn(header_name) → True keeps the column as always-visible.
    Core columns are sorted to match the T3 column order.
    """
    core_idx = sorted(
        [i for i, h in enumerate(headers) if core_fn(h)],
        key=lambda i: t3_order.get(headers[i], len(t3h) + i)
    )
    ext_idx  = [i for i in range(len(headers)) if not core_fn(headers[i])]
    order    = core_idx + ext_idx

    def pick(row):
        return [row[i] if i < len(row) else '' for i in order]

    new_meta = pick(meta) if meta else None
    return (
        [headers[i] for i in order],
        new_meta,
        [pick(row) for row in data],
        len(core_idx),
    )

t2h_r, t2_meta_r, t2d_r, t2_core = reorder(t2h, t2_meta, t2d, lambda h: h in t3_set)
t1h_r, t1_meta_r, t1d_r, t1_core = reorder(t1h, t1_meta, t1d,
                                              lambda h: h in t3_set or h in T1_EXTRA)

# ── Cell renderer ──────────────────────────────────────────────────────────────
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
        if 'pathogenic' in vl:
            return f'<span class="path">{HL.escape(val)}</span>'
        if 'benign' in vl:
            return f'<span class="ben">{HL.escape(val)}</span>'
    return HL.escape(val)

# ── HTML table builder ─────────────────────────────────────────────────────────
def build_table(tbl_id, headers, meta, data, comments=None):
    p = []

    # Comment / version line above the table (Tier 1 only)
    if comments:
        for c in comments:
            p.append(f'<div class="ver-line"><code>{HL.escape(c)}</code></div>')

    p.append(f'<div class="tbl-wrap"><table class="dt" id="{tbl_id}">')

    # ── thead ──
    p.append('<thead>')

    # Row 1: column names (clickable for sort)
    p.append('<tr class="col-names">')
    for i, h in enumerate(headers):
        esc = HL.escape(h)
        p.append(f'<th onclick="srt(\'{tbl_id}\',{i})" title="{esc}">'
                 f'{esc}<span class="si"> &#8597;</span></th>')
    p.append('</tr>')

    # Row 2: source metadata (Tier 2 and Tier 3 only)
    if meta:
        p.append('<tr class="col-meta">')
        for m in meta:
            esc = HL.escape(m)
            p.append(f'<th title="{esc}">{esc}</th>')
        p.append('</tr>')

    p.append('</thead><tbody>')

    # ── tbody ──
    for row in data:
        vals   = [row[i] if i < len(row) else '' for i in range(len(headers))]
        search = HL.escape('|'.join(v.lower() for v in vals))
        p.append(f'<tr data-s="{search}">')
        for i, val in enumerate(vals):
            cell  = render_cell(headers[i], val)
            title = f' title="{HL.escape(val)}"' if len(val) > 50 else ''
            p.append(f'<td{title}>{cell}</td>')
        p.append('</tr>')

    p.append('</tbody></table></div>')
    return ''.join(p)

tbl_t3 = build_table('tbl-t3', t3h,   t3_meta,   t3d,   t3_comments)
tbl_t2 = build_table('tbl-t2', t2h_r, t2_meta_r, t2d_r, t2_comments)
tbl_t1 = build_table('tbl-t1', t1h_r, t1_meta_r, t1d_r, t1_comments)

# ── Approximate height of col-names row for sticky positioning ─────────────────
# padding 6px top+bottom + ~17px line-height ≈ 29 px; use 32 for safety
META_TOP = 32

# ── Full HTML page ─────────────────────────────────────────────────────────────
PAGE = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>SMART — Example Output | Verification 1</title>
  <style>
    *,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
    body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
         font-size:14px;background:#f8f9fa;color:#333}}

    /* ── Header ── */
    header{{background:#2c3e50;color:#fff;padding:16px 28px}}
    header h1{{font-size:1.3em;margin-bottom:4px}}
    header .sub{{font-size:.82em;color:#adc8da}}
    header a{{color:#adc8da;text-decoration:none}}
    header a:hover{{color:#fff}}

    /* ── Layout ── */
    .container{{padding:20px 28px 10px}}
    .info-box{{background:#e8f4fd;border-left:4px solid #3498db;padding:10px 16px;
               border-radius:0 5px 5px 0;margin-bottom:20px;font-size:.88em;line-height:1.6}}
    .info-box a{{color:#1a6ea8}}

    /* ── Tabs ── */
    .tabs{{display:flex;border-bottom:2px solid #ddd;margin-bottom:18px}}
    .tab-btn{{padding:9px 22px;border:none;background:none;cursor:pointer;font-size:.92em;
              color:#666;border-bottom:3px solid transparent;margin-bottom:-2px}}
    .tab-btn:hover{{color:#2c3e50;background:#f0f0f0}}
    .tab-btn.active{{color:#2c3e50;font-weight:600;border-bottom-color:#2c3e50}}
    .badge{{display:inline-block;background:#e0e0e0;color:#555;border-radius:10px;
            padding:0 7px;font-size:.75em;margin-left:5px}}
    .tab-btn.active .badge{{background:#2c3e50;color:#fff}}

    /* ── Panels ── */
    .panel{{display:none}}.panel.active{{display:block}}
    .panel-desc{{margin-bottom:12px;font-size:.88em;color:#555;line-height:1.6}}
    .panel-desc code{{background:#f0f0f0;padding:1px 5px;border-radius:3px;font-size:.92em}}

    .controls{{display:flex;gap:10px;align-items:center;margin-bottom:8px;flex-wrap:wrap}}
    .search-box{{padding:5px 10px;border:1px solid #ccc;border-radius:4px;
                 font-size:.88em;min-width:240px;flex:1}}
    .btn{{padding:5px 14px;border-radius:4px;border:1px solid #bbb;background:#fff;
          cursor:pointer;font-size:.84em;color:#333;white-space:nowrap}}
    .btn:hover{{background:#f0f0f0}}
    .btn.on{{background:#2c3e50;color:#fff;border-color:#2c3e50}}
    .hint{{font-size:.8em;color:#888}}

    /* ── Version / comment line (Tier 1) ── */
    .ver-line{{font-family:monospace;font-size:.84em;background:#f4f4f4;
               border:1px solid #dee2e6;border-bottom:none;border-radius:4px 4px 0 0;
               padding:5px 10px;color:#555;letter-spacing:.01em}}

    /* ── Table wrapper & table ── */
    .tbl-wrap{{overflow:auto;max-height:65vh;border:1px solid #dee2e6;border-radius:0 0 4px 4px}}
    .ver-line + .tbl-wrap{{border-top:none}}
    table.dt{{border-collapse:collapse;font-size:.79em;white-space:nowrap;width:100%}}

    /* Column names row — sticky at top */
    table.dt thead tr.col-names th{{
      position:sticky;top:0;z-index:3;
      background:#2c3e50;color:#fff;
      padding:6px 10px;cursor:pointer;user-select:none;
      border-right:1px solid #3d5368
    }}
    table.dt thead tr.col-names th:hover{{background:#3d5368}}
    .si{{opacity:.55;font-size:.72em}}

    /* Metadata row — sticky just below the names row */
    table.dt thead tr.col-meta th{{
      position:sticky;top:{META_TOP}px;z-index:2;
      background:#3a5066;color:#b8d0e4;
      padding:4px 10px;font-weight:normal;font-size:.72em;
      cursor:default;border-right:1px solid #4a6278;
      max-width:260px;overflow:hidden;text-overflow:ellipsis
    }}

    /* Data rows */
    table.dt tbody td{{padding:4px 8px;border-bottom:1px solid #eee;
                       max-width:280px;overflow:hidden;text-overflow:ellipsis}}
    table.dt tbody tr:nth-child(even){{background:#f8f9fa}}
    table.dt tbody tr:hover td{{background:#ebf3fb}}
    tr.hidden{{display:none}}

    /* Extended columns hidden by default; .show-all reveals them */
    #tbl-t2:not(.show-all) th:nth-child(n+{t2_core + 1}),
    #tbl-t2:not(.show-all) td:nth-child(n+{t2_core + 1}){{display:none}}
    #tbl-t1:not(.show-all) th:nth-child(n+{t1_core + 1}),
    #tbl-t1:not(.show-all) td:nth-child(n+{t1_core + 1}){{display:none}}

    /* ── Badges ── */
    .lvl{{display:inline-block;padding:1px 6px;border-radius:3px;font-size:.8em;font-weight:700}}
    .L1{{background:#1a9641;color:#fff}}  .L2{{background:#0571b0;color:#fff}}
    .L3A{{background:#e66101;color:#fff}} .L3B{{background:#fdb863;color:#333}}
    .L4{{background:#92c5de;color:#333}}  .LR1{{background:#ca0020;color:#fff}}
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

    /* ── Footer ── */
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
    (14&nbsp;SNV/indel &plus; 4&nbsp;CNA) covering NRAS, IDH1, PIK3CA, EGFR, BRAF, GNAQ, PTEN,
    KRAS, BRCA2, DICER1, TP53, ERBB2, MET, CDKN2A&nbsp;(&times;2&nbsp;isoforms), and CDK4.
    The CDKN2A deletion appears as two rows &mdash; one per isoform (p16/INK4a &amp; p14ARF).
    Tables match the actual file format: Tier&nbsp;2 and Tier&nbsp;3 show both header rows
    (field&nbsp;names&nbsp;+ source&nbsp;metadata); Tier&nbsp;1 shows the
    <code>#SMART_VERSION</code> comment line and MAF column headers.
    Click any column name to sort &uarr;&darr;. Use the search box to filter rows.
  </div>

  <div class="tabs">
    <button class="tab-btn active" onclick="showTab('t3',this)">
      Tier&nbsp;3 &mdash; Clinical <span class="badge">{len(t3h)}&nbsp;cols &middot; {len(t3d)}&nbsp;rows</span>
    </button>
    <button class="tab-btn" onclick="showTab('t2',this)">
      Tier&nbsp;2 &mdash; Bioinformatics <span class="badge">{len(t2h)}&nbsp;cols &middot; {len(t2d)}&nbsp;rows</span>
    </button>
    <button class="tab-btn" onclick="showTab('t1',this)">
      Tier&nbsp;1 &mdash; Full&nbsp;MAF <span class="badge">{len(t1h)}&nbsp;cols &middot; {len(t1d)}&nbsp;rows</span>
    </button>
  </div>

  <!-- ── Tier 3 ── -->
  <div id="panel-t3" class="panel active">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier3.tsv</code></strong> &mdash;
      Clinically focused. {len(t3h)}&nbsp;columns, two-row header
      (field&nbsp;names + <em>description&nbsp;| source&nbsp;| version</em>).
      Includes variant identity, consequence, population frequency, pathogenicity scores,
      ClinVar, CancerHotspots, and top-line OncoKB actionability.
    </p>
    <div class="controls">
      <input class="search-box" type="text" placeholder="Search all fields&hellip;"
             oninput="filterRows('tbl-t3',this.value)">
    </div>
    {tbl_t3}
  </div>

  <!-- ── Tier 2 ── -->
  <div id="panel-t2" class="panel">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier2.tsv</code></strong> &mdash;
      Bioinformatics view. {len(t2h)}&nbsp;columns, two-row header
      (field&nbsp;names + <em>description&nbsp;| source&nbsp;| version</em>).
      Extends Tier&nbsp;3 with full gnomAD stratification, all SpliceAI scores,
      complete CIViC fields, and all OncoKB JSON expansions
      (<code>ONCOKB_TX_*</code>, <code>ONCOKB_DIAG_*</code>, <code>ONCOKB_PROG_*</code>).
      Core columns (matching Tier&nbsp;3) are shown first.
    </p>
    <div class="controls">
      <input class="search-box" type="text" placeholder="Search all fields&hellip;"
             oninput="filterRows('tbl-t2',this.value)">
      <button id="btn-t2" class="btn" onclick="toggleExt('t2')">
        Show all {len(t2h)}&nbsp;columns
      </button>
      <span class="hint">Default: {t2_core}&nbsp;core columns shown.</span>
    </div>
    {tbl_t2}
  </div>

  <!-- ── Tier 1 ── -->
  <div id="panel-t1" class="panel">
    <p class="panel-desc">
      <strong>File: <code>Final_result_tier1.maf</code></strong> &mdash;
      Standard MAF format. {len(t1h)}&nbsp;columns. First line is the
      <code>#SMART_VERSION</code> comment; second line is the MAF column-name header
      (no separate metadata row in MAF format). Designed for cBioPortal, R/Python pipelines,
      and oncoPrint generators.
    </p>
    <div class="controls">
      <input class="search-box" type="text" placeholder="Search all fields&hellip;"
             oninput="filterRows('tbl-t1',this.value)">
      <button id="btn-t1" class="btn" onclick="toggleExt('t1')">
        Show all {len(t1h)}&nbsp;columns
      </button>
      <span class="hint">Default: {t1_core}&nbsp;core columns shown.</span>
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
  var nc  = tier==='t2' ? {len(t2h)}  : {len(t1h)};
  var nc3 = tier==='t2' ? {t2_core} : {t1_core};
  btn.textContent = extOn[tier] ? 'Show core only ('+nc3+')' : 'Show all '+nc+' columns';
}}

function srt(tblId, col) {{
  var tbl   = document.getElementById(tblId);
  var tbody = tbl.querySelector('tbody');
  var rows  = Array.from(tbody.querySelectorAll('tr'));
  var asc   = tbl.dataset.sc == col && tbl.dataset.sd === 'a';
  rows.sort(function(a,b) {{
    var av = a.cells[col] ? a.cells[col].textContent.trim() : '';
    var bv = b.cells[col] ? b.cells[col].textContent.trim() : '';
    var n=parseFloat(av), m=parseFloat(bv);
    if (!isNaN(n) && !isNaN(m)) return asc ? m-n : n-m;
    return asc ? bv.localeCompare(av) : av.localeCompare(bv);
  }});
  rows.forEach(function(r){{tbody.appendChild(r)}});
  tbl.dataset.sc=col; tbl.dataset.sd=asc?'d':'a';
}}
</script>

</body>
</html>
"""

with open(OUT, 'w', encoding='utf-8') as f:
    f.write(PAGE)

size_kb = os.path.getsize(OUT) // 1024
print(f'Written: {OUT}  ({size_kb} KB)')
print(f'T3: {len(t3h)} cols x {len(t3d)} rows  meta={"yes" if t3_meta else "no"}  comment={t3_comments}')
print(f'T2: {len(t2h)} cols x {len(t2d)} rows  meta={"yes" if t2_meta else "no"}  core_visible={t2_core}')
print(f'T1: {len(t1h)} cols x {len(t1d)} rows  meta={"yes" if t1_meta else "no"}  core_visible={t1_core}  comment={t1_comments}')
