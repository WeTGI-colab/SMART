#!/usr/bin/env python3
"""Generate docs/fields.md — the public field reference — from Config.yaml.

Why this script exists
----------------------
`docs/fields.md` is linked from the README as "Full field reference" and is the page a
reporting scientist opens to find out what an output column means. It was maintained by
hand, and hand maintenance failed in the way hand maintenance always fails:

  * it listed 354 concrete columns against Config.yaml's 427 — 73 missing, among them
    every dbNSFP predictor, every COSMIC CMC column, all ten GENIE_* columns, CADD,
    ANNOTATION_FLAGS, SPLICEAI_MAX_DS, TRANSCRIPT_MATCH, TRANSCRIPT_REQUESTED and
    ONCOKB_MAF_QUERY;
  * it carried 10 `ONCOKB_PROG_*_…` names that match no entry and no pattern in
    Config.yaml, so a reader could look up a column the configuration does not describe;
  * LOEUF's source was still a bare citation that the root Config.yaml had already
    corrected (Karczewski 2020 is not where those values come from);
  * the page's own JavaScript filter read an undeclared `tier` variable, so the search
    box threw on every keystroke and filtered nothing.

So the page is generated. Config.yaml is the single thing to edit; this script renders it,
and tests/verification5/test_fields_doc.py fails if the two have drifted apart.

Output is deterministic — no timestamps, no run-dependent anything — so that
`--check` is meaningful and a regeneration produces no spurious diff.

Usage
-----
    python3 utils/gen_fields_doc.py            # rewrite docs/fields.md
    python3 utils/gen_fields_doc.py --check    # exit 1 if docs/fields.md is stale
"""

import argparse
import collections
import fnmatch
import html
import os
import sys

import yaml

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
CONFIG = os.path.join(ROOT, "Config.yaml")
OUT = os.path.join(ROOT, "docs", "fields.md")
SCRIPT = "utils/gen_fields_doc.py"

# --------------------------------------------------------------------------------------
# Source grouping
# --------------------------------------------------------------------------------------
# Config.yaml records `source` as free text, and 53 distinct strings appear in it — among
# them 34 variants of "dbNSFP (<predictor>)" and one 60-word paragraph about the gnomAD
# constraint release. Grouping on the raw string would produce 53 sections, most with one
# row. These rules fold the strings into the resource families a reader actually thinks in,
# while every table still prints the full `source` text unaltered.

GROUP_ORDER = [
    "VCF", "VEP", "Computed by SMART", "OncoKB", "ClinVar", "CIViC",
    "Cancer Hotspots", "COSMIC", "AACR Project GENIE", "SpliceAI",
    "dbNSFP", "REVEL", "CADD", "gnomAD", "Manta",
]

# These notes are rendered inside a block-level <div>, where kramdown does not parse
# Markdown, so they are written in HTML rather than with backticks.
GROUP_NOTES = {
    "VCF": "Columns the variant caller wrote, or that are copied straight out of the "
           "input VCF record. Their descriptions and versions are harvested from the "
           "VCF's own <code>##INFO</code>/<code>##FORMAT</code> headers at run time "
           "where the file provides them (see <code>docs/field_provenance.md</code>), "
           "so a run may report a version more specific than the one below.",
    "VEP": "Columns from Ensembl VEP's CSQ field, including its plugins. The version is "
           "the VEP release pinned in the <code>Dockerfile</code>.",
    "Computed by SMART": "Columns no external resource provides: SMART derives them in "
                         "<code>scripts/post_analysis.py</code> from the columns above.",
    "OncoKB": "Columns from the OncoKB API and the OncoKB MafAnnotator. OncoKB is "
              "versioned by data release rather than by file, so these entries point at "
              "the <code>ONCOKB_DATA_VERSION</code> column of the run itself. The "
              "indexed <code>ONCOKB_TX_*</code>/<code>ONCOKB_DIAG_*</code> columns are "
              "documented by pattern, in "
              "<a href=\"#columns-documented-by-pattern\">the section below</a>.",
    "dbNSFP": "Pathogenicity predictors and conservation scores read from the dbNSFP "
              "release built by <code>utils/build_dbnsfp.py</code>. Each predictor "
              "contributes its own score and, usually, its own rankscore.",
    "gnomAD": "Population allele frequencies and the gene/transcript constraint metrics.",
    "AACR Project GENIE": "Counts from the GENIE somatic table built by "
                          "<code>utils/build_genie_o4.py</code>, used by the SVIG-UK O4 "
                          "code.",
}

# --------------------------------------------------------------------------------------
# Tiers
# --------------------------------------------------------------------------------------
# Read from the routing loop in scripts/post_analysis.py, which is the only thing that
# decides it. Not from a prose description of it elsewhere: create_config/README.md's table
# had Tier 1 and Tier 3 the other way round for as long as it existed. The code:
#
#     if tier == "drop":  continue                  # dropped from every output
#     final_cols.append(col)                        # -> Final_result_tier1.maf
#     if tier == 2:      tier2_cols.append(col)
#     elif tier == 3:    tier3_cols.append(col); tier2_cols.append(col)
#
# so tier 3 reaches both TSVs, tier 2 only the Tier 2 TSV, and anything else — tier 1
# included — only the MAF.

TIER_ROUTING = [
    ("3", "Yes", "Yes", "Yes",
     "The clinical set: the narrowest, most interpretable columns."),
    ("2", "&ndash;", "Yes", "Yes",
     "The bioinformatician set: supporting evidence and the raw scores behind Tier 3."),
    ("1", "&ndash;", "&ndash;", "Yes",
     "Everything else that is kept. <code>tier: 1</code> is also what any unrecognised "
     "value behaves as, and is the default for a column with no entry at all."),
    ("drop", "&ndash;", "&ndash;", "&ndash;",
     "Removed from every output — intermediates, duplicated columns and the JSON blobs "
     "that have already been expanded into columns of their own."),
]

TIER_RANK = {3: 0, 2: 1, 1: 2, "drop": 3}
TIER_BADGE = {3: "t3", 2: "t2", 1: "t1", "drop": "drop"}


def tier_key(tier):
    """Sort key for a tier: most widely published first, dropped last."""
    return TIER_RANK.get(tier, 2.5)


def version_cell(version):
    """The version, or a visible statement that the configuration gives none.

    Three entries (VAF, SPLICEAI_MAX_DS, ANNOTATION_FLAGS) legitimately have no `version`:
    SMART computes them, so there is no external release to name. An empty cell reads like
    an omission, which is exactly the doubt this page exists to remove.
    """
    text = ("" if version is None else str(version)).strip()
    return f"<small>{esc(text)}</small>" if text else "<small><em>n/a</em></small>"


def tier_cell(tier):
    label = "drop" if tier == "drop" else ("&mdash;" if tier is None else str(tier))
    return f'<span class="badge {TIER_BADGE.get(tier, "t1")}">{label}</span>'


# --------------------------------------------------------------------------------------
# Reading Config.yaml
# --------------------------------------------------------------------------------------

def load(path=CONFIG):
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def entries(cfg):
    """The concrete column entries, name -> metadata, in Config.yaml's own order."""
    return dict(cfg.get("fields") or {})


def patterns(cfg):
    """The fnmatch patterns, de-duplicated on the pattern itself.

    Config.yaml currently lists 37 entries of which two are identical
    (`ONCOKB_TX_*_levelAssociatedCancerType.mainType.name`), so 36 distinct patterns are
    documented. A later duplicate that disagreed with the first would be a real
    configuration error, so it is reported rather than silently dropped.
    """
    out, seen, clashes = [], {}, []
    for p in cfg.get("field_patterns") or []:
        key = p.get("pattern", "")
        if key in seen:
            if seen[key] != p:
                clashes.append(key)
            continue
        seen[key] = p
        out.append(p)
    for key in clashes:
        print(f"warning: field_patterns lists {key} twice with different metadata; "
              f"documenting the first", file=sys.stderr)
    return out


# post_analysis.py expands three OncoKB JSON arrays into indexed columns
# (_ONCOKB_JSON_FIELDS: treatments -> ONCOKB_TX, diagnosticImplications -> ONCOKB_DIAG,
# prognosticImplications -> ONCOKB_PROG). Config.yaml carries patterns for the first two
# only, which the page says out loud rather than leaving the reader to wonder -- but it
# says it only while it is true, so that adding the missing patterns removes the notice.
PROG_PROBES = ("ONCOKB_PROG_0_levelOfEvidence", "ONCOKB_PROG_0_description",
               "ONCOKB_PROG_0_tumorType.name")


def undocumented_prog_columns(field_entries, pats):
    """The representative ONCOKB_PROG_* columns that nothing in the config describes."""
    known = [p.get("pattern", "") for p in pats]
    return [c for c in PROG_PROBES
            if c not in field_entries
            and not any(fnmatch.fnmatch(c, k) for k in known)]


def source_group(source):
    """Which section a `source` string belongs in."""
    s = (source or "").strip()
    low = s.lower()
    if low.startswith("dbnsfp"):
        return "dbNSFP"
    if low.startswith("manta") or low.startswith("vcf_manta"):
        return "Manta"
    if low.startswith("gnomad"):
        return "gnomAD"
    if low.startswith("vcf") or low.startswith("input vcf"):
        return "VCF"
    for family in ("VEP", "OncoKB", "ClinVar", "CIViC", "Cancer Hotspots", "COSMIC",
                   "SpliceAI", "REVEL", "CADD"):
        if low.startswith(family.lower()):
            return family
    if low.startswith("aacr project genie") or low.startswith("genie"):
        return "AACR Project GENIE"
    if low.startswith("smart") or low.startswith("computed"):
        return "Computed by SMART"
    # An unrecognised source is never dropped: it gets its own section, named after the
    # text before any parenthetical, and sorts after the known families.
    return s.split(" (")[0] or "Unattributed"


def grouped(field_entries):
    """{group: [(name, meta), ...]} with the groups and rows in their final order."""
    buckets = collections.defaultdict(list)
    for name, meta in field_entries.items():
        buckets[source_group(str(meta.get("source", "")))].append((name, meta))
    for rows in buckets.values():
        rows.sort(key=lambda nm: (tier_key(nm[1].get("tier")), nm[0].lower()))
    known = [g for g in GROUP_ORDER if g in buckets]
    rest = sorted(g for g in buckets if g not in GROUP_ORDER)
    return collections.OrderedDict((g, buckets[g]) for g in known + rest)


# --------------------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------------------

def esc(value):
    """A Config.yaml value as table-cell text: escaped, and on one line.

    Two things have to be dealt with.

    Markup: descriptions carry `<`, `>`, `&` and the pipe-separated CSQ field lists
    verbatim, so everything is escaped; nothing in Config.yaml is intended as markup.

    Line breaks: ten entries -- the SpliceAI block -- are written in YAML styles that keep
    their newlines, and an embedded newline would split one table row across several lines
    of the page. In Markdown a blank line inside a block of HTML ends that block, so a row
    could be published as literal text. `scripts/post_analysis.py` collapses the same
    whitespace for the same reason before writing a MAF record (`_EMBEDDED_WS_RE`).
    """
    text = "" if value is None else str(value)
    return html.escape(" ".join(text.split()), quote=True)


def anchor(text):
    """A stable id for a group heading, so a section can be linked to."""
    keep = [c.lower() if c.isalnum() else "-" for c in text]
    return "source-" + "".join(keep).strip("-")


def plural(n, noun):
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


STYLE = """<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
.controls { margin: 1.2em 0; display: flex; gap: 0.8em; flex-wrap: wrap; align-items: center;
            position: sticky; top: 0; background: #fff; padding: 0.6em 0; z-index: 5; }
.controls input { flex: 1; min-width: 200px; padding: 0.45em 0.7em; border: 1px solid #ccc; border-radius: 4px; font-size: 0.95em; }
.controls select { padding: 0.45em 0.7em; border: 1px solid #ccc; border-radius: 4px; font-size: 0.95em; }
.controls label { font-size: 0.9em; color: #555; }
table { width: 100%; border-collapse: collapse; font-size: 0.87em; }
th { background: #2c3e50; color: #fff; padding: 0.5em 0.7em; text-align: left; }
td { padding: 0.42em 0.7em; border-bottom: 1px solid #eee; vertical-align: top; }
td.src, td.ver { font-size: 0.92em; color: #444; }
tr:hover td { background: #f5f8ff; }
tr.hidden, .group.hidden { display: none; }
code { background: #f0f0f0; padding: 0.1em 0.35em; border-radius: 3px; font-size: 0.9em; word-break: break-all; }
.badge { display: inline-block; padding: 0.15em 0.5em; border-radius: 3px; font-size: 0.8em; font-weight: 600; white-space: nowrap; }
.badge.t3 { background: #d4edda; color: #155724; }
.badge.t2 { background: #cce5ff; color: #004085; }
.badge.t1 { background: #fff3cd; color: #856404; }
.badge.drop { background: #f8d7da; color: #721c24; }
.warn { background: #fff8e1; border-left: 4px solid #f0ad4e; padding: 0.8em 1em; margin: 1.2em 0; border-radius: 0 4px 4px 0; font-size: 0.9em; }
.generated { background: #eef4ff; border-left: 4px solid #2c3e50; padding: 0.8em 1em; margin: 1.2em 0; border-radius: 0 4px 4px 0; font-size: 0.9em; }
.groupnote { font-size: 0.9em; color: #555; margin: 0.4em 0 0.8em; }
#result-count { font-size: 0.9em; color: #555; margin-top: 0.3em; }
</style>"""


def render(cfg, config_name="Config.yaml", script_name=SCRIPT):
    field_entries = entries(cfg)
    pats = patterns(cfg)
    groups = grouped(field_entries)
    per_tier = collections.Counter(e.get("tier") for e in field_entries.values())
    n_pub = sum(n for t, n in per_tier.items() if t != "drop")

    L = []
    w = L.append

    # --- Jekyll front matter (this page is published at /SMART/fields/) -----------------
    w("---")
    w("layout: default")
    w("title: SMART — Field Reference")
    w("permalink: /fields/")
    w("---")
    w("")
    w(f"<!-- GENERATED FILE — DO NOT EDIT BY HAND.")
    w(f"     Written by {script_name} from {config_name}:")
    w(f"     {len(field_entries)} column entries and {len(pats)} field patterns.")
    w(f"     To change this page, edit {config_name} and run:")
    w(f"         python3 {script_name}")
    w(f"     tests/verification5/test_fields_doc.py fails if the two disagree. -->")
    w("")
    w(STYLE)
    w("")
    w("# SMART — Field Reference")
    w("")
    w("SMART (**S**omatic **M**utation **A**nnotation and **R**eporting **T**ool) integrates")
    w("annotations from multiple sources into a single output table. This page describes every")
    w("field the pipeline is configured to produce — what it contains, where it comes from,")
    w("which version of that tool or database produced it, and which of the three output files")
    w("it reaches.")
    w("")
    w('<div class="generated">')
    w(f"<strong>This page is generated.</strong> It is rendered from")
    w(f"<code>{esc(config_name)}</code> — the pipeline's own field configuration — by")
    w(f"<code>{esc(script_name)}</code>. Do not edit it by hand: edit")
    w(f"<code>{esc(config_name)}</code> and run <code>python3 {esc(script_name)}</code>.")
    w("Hand maintenance is what left the previous version of this page 73 columns short of")
    w("the configuration, and <code>tests/verification5/test_fields_doc.py</code> now fails")
    w("if the page and the configuration disagree.")
    w("</div>")
    w("")
    w("&#128200; **[View example output →](example-output.html)**")
    w("Browse all three output files from the Verification 1 test run (22 variants, 1 sample)")
    w("as interactive, searchable tables with colour-coded OncoKB levels.")
    w("")

    # --- Counts -----------------------------------------------------------------------
    w("## What is on this page")
    w("")
    w(f"`{config_name}` describes **{len(field_entries)} columns** by name, grouped below into")
    w(f"{len(groups)} sources, plus **{len(pats)} patterns** covering the indexed columns whose")
    w("names are only known once a variant has been annotated:")
    w("")
    w("| | Count |")
    w("|---|---|")
    w(f"| Columns described by name | {len(field_entries)} |")
    w(f"| &nbsp;&nbsp;of which reach an output file | {n_pub} |")
    w(f"| &nbsp;&nbsp;of which are dropped (`tier: drop`) | {per_tier.get('drop', 0)} |")
    w(f"| &nbsp;&nbsp;in the Tier 3 clinical sheet (`tier: 3`) | {per_tier.get(3, 0)} |")
    w(f"| &nbsp;&nbsp;in the Tier 2 sheet only (`tier: 2`) | {per_tier.get(2, 0)} |")
    w(f"| &nbsp;&nbsp;in the Tier 1 MAF only (`tier: 1`) | {per_tier.get(1, 0)} |")
    w(f"| Patterns covering indexed columns | {len(pats)} |")
    w("")
    w("A single run produces more columns than this: the patterns expand to as many columns as")
    w("the data needs, and columns a VCF documents for itself need no entry here. The Tier 1")
    w("MAF of a real run has carried around 1,000.")
    w("")

    # --- Tiers ------------------------------------------------------------------------
    w("## What the tiers mean")
    w("")
    w("Every entry carries a `tier`, and the tier decides which of the three output files the")
    w("column appears in. The table below is read from the routing loop in")
    w("`scripts/post_analysis.py`, which is the only thing that decides it:")
    w("")
    w("| `tier` in the configuration | `Final_result_tier3.tsv` | `Final_result_tier2.tsv` | `Final_result_tier1.maf` | What it is for |")
    w("|---|---|---|---|---|")
    for value, t3, t2, t1, note in TIER_ROUTING:
        w(f"| `{value}` | {t3} | {t2} | {t1} | {note} |")
    w("")
    w("Two consequences are worth stating plainly, because they are easy to get backwards:")
    w("")
    w("- **Tier 3 is the narrow clinical sheet and Tier 2 contains it.** A `tier: 3` column is")
    w("  written to the Tier 2 sheet as well, so Tier 2 is a superset of Tier 3, and the Tier 1")
    w("  MAF is a superset of both.")
    w("- **Tier 1 is not a curated set, and is not a sheet.** It is where a column lands when")
    w("  nothing else claims it, including a column with no entry in the configuration at all,")
    w("  and such a column appears in neither TSV. To put a column in front of a clinical")
    w("  scientist, give it `tier: 3`; to give it to a bioinformatician only, `tier: 2`.")
    w("")

    # --- Controls ---------------------------------------------------------------------
    w("## Field reference")
    w("")
    w("Within each source, columns are ordered by the tier they reach — Tier 3 first, dropped")
    w("columns last — and alphabetically within a tier.")
    w("")
    w('<div class="controls">')
    w('  <input type="text" id="search" placeholder="Search field name or description…" oninput="filterTable()">')
    w("  <label>Source:")
    w('    <select id="src-filter" onchange="filterTable()">')
    w('      <option value="">All sources</option>')
    for group in groups:
        w(f"      <option>{esc(group)}</option>")
    w("      <option>Patterns</option>")
    w("    </select>")
    w("  </label>")
    w("  <label>Tier:")
    w('    <select id="tier-filter" onchange="filterTable()">')
    w('      <option value="">All tiers</option>')
    w('      <option value="3">3 — clinical sheet</option>')
    w('      <option value="2">2 — bioinformatician sheet</option>')
    w('      <option value="1">1 — MAF only</option>')
    w('      <option value="drop">drop — not output</option>')
    w("    </select>")
    w("  </label>")
    w("</div>")
    w('<div id="result-count"></div>')
    w("")

    for group, rows in groups.items():
        w(f'<div class="group" data-group="{esc(group)}">')
        w("")
        w(f'<h3 id="{anchor(group)}">{esc(group)} '
          f"<small>({plural(len(rows), 'column')})</small></h3>")
        w("")
        if GROUP_NOTES.get(group):
            w(f'<p class="groupnote">{GROUP_NOTES[group]}</p>')
            w("")
        w("<table>")
        w("<thead>")
        w("<tr><th>Field</th><th>Description</th><th>Source</th><th>Version</th><th>Tier</th></tr>")
        w("</thead>")
        w("<tbody>")
        for name, meta in rows:
            tier = meta.get("tier")
            w(f'<tr data-field="{esc(name)}" data-source="{esc(group)}" '
              f'data-tier="{esc("drop" if tier == "drop" else tier)}">'
              f"<td><code>{esc(name)}</code></td>"
              f"<td>{esc(meta.get('description', ''))}</td>"
              f'<td class="src">{esc(meta.get("source", ""))}</td>'
              f'<td class="ver">{version_cell(meta.get("version"))}</td>'
              f"<td>{tier_cell(tier)}</td></tr>")
        w("</tbody>")
        w("</table>")
        w("")
        w("</div>")
        w("")

    # --- Patterns ---------------------------------------------------------------------
    w("## Columns documented by pattern")
    w("")
    w("OncoKB returns its treatment, diagnostic and prognostic implications as JSON arrays, and")
    w("`scripts/post_analysis.py` expands each array element into columns of its own —")
    w("`ONCOKB_TX_0_level`, `ONCOKB_TX_1_level`, and so on for as many implications as the")
    w("variant has. How many there are is a property of the data, so these columns cannot be")
    w("listed by name in advance.")
    w("")
    w("They are documented by `fnmatch` pattern instead. At run time")
    w("`scripts/vcf_metadata.py` (`resolve()`) tries the exact column name against the entries")
    w("above first, and only then walks the patterns in order, taking the metadata of the first")
    w("one that matches. A column matching none of either is reported as undocumented by")
    w("`report_undocumented_columns()` and still appears in the output — it simply carries no")
    w("description in the metadata header row.")
    w("")
    w("The `*` stands for the array index. These are patterns, not column names: the table")
    w("below is reproduced from the configuration exactly, and is not expanded into concrete")
    w("column names, because the expansion depends on the run.")
    w("")
    if undocumented_prog_columns(field_entries, pats):
        w('<div class="warn">')
        w("<strong>One gap, stated rather than hidden.</strong>")
        w("<code>scripts/post_analysis.py</code> expands three OncoKB arrays in this way:")
        w("treatments into <code>ONCOKB_TX_*</code>, diagnostic implications into")
        w("<code>ONCOKB_DIAG_*</code> and prognostic implications into")
        w("<code>ONCOKB_PROG_*</code>. The configuration has patterns for the first two")
        w("only, so the <code>ONCOKB_PROG_*</code> columns reach the Tier 1 MAF with no")
        w("description, source or version, and")
        w("<code>report_undocumented_columns()</code> lists them when a run finishes. They")
        w("are absent from the table below because the configuration does not describe")
        w("them &mdash; not because the columns do not exist. Adding")
        w("<code>ONCOKB_PROG_*</code> patterns to the configuration documents them and")
        w("removes this notice.")
        w("</div>")
        w("")
    w('<div class="group" data-group="Patterns">')
    w("")
    w("<table>")
    w("<thead>")
    w("<tr><th>Pattern</th><th>Description</th><th>Source</th><th>Version</th><th>Tier</th></tr>")
    w("</thead>")
    w("<tbody>")
    for p in pats:
        tier = p.get("tier")
        w(f'<tr data-pattern="{esc(p.get("pattern", ""))}" data-source="Patterns" '
          f'data-tier="{esc("drop" if tier == "drop" else tier)}">'
          f"<td><code>{esc(p.get('pattern', ''))}</code></td>"
          f"<td>{esc(p.get('description', ''))}</td>"
          f'<td class="src">{esc(p.get("source", ""))}</td>'
          f'<td class="ver">{version_cell(p.get("version"))}</td>'
          f"<td>{tier_cell(tier)}</td></tr>")
    w("</tbody>")
    w("</table>")
    w("")
    w("</div>")
    w("")

    # --- Where the metadata comes from -------------------------------------------------
    w("## Where a column's description, source and version come from")
    w("")
    w("The table above is the configuration. A run can be more specific than the")
    w("configuration, and where it is, the run wins: a VCF documents its own fields in its")
    w("`##INFO`/`##FORMAT` headers, for the exact versions that produced them, and")
    w("`scripts/vcf_metadata.py` harvests those headers and prefers them. The OncoKB columns")
    w("likewise take the data version observed in the run. `docs/field_provenance.md` sets out")
    w("the precedence in full.")
    w("")
    w("So the second header row of a Tier 2 or Tier 3 sheet may carry a more precise source or")
    w("version than this page does. Where they differ, the sheet is describing the run and this")
    w("page is describing the configuration.")
    w("")

    # --- Script -----------------------------------------------------------------------
    w("<script>")
    w("function filterTable() {")
    w("  const search = document.getElementById('search').value.toLowerCase();")
    w("  const src    = document.getElementById('src-filter').value;")
    w("  const tier   = document.getElementById('tier-filter').value;")
    w("  let visible = 0;")
    w("  document.querySelectorAll('.group').forEach(group => {")
    w("    let shown = 0;")
    w("    group.querySelectorAll('tbody tr').forEach(row => {")
    w("      const text    = row.textContent.toLowerCase();")
    w("      const rowSrc  = row.dataset.source || '';")
    w("      const rowTier = row.dataset.tier   || '';")
    w("      const match = (!search || text.includes(search))")
    w("                 && (!src    || rowSrc === src)")
    w("                 && (!tier   || rowTier === tier);")
    w("      row.classList.toggle('hidden', !match);")
    w("      if (match) { shown++; visible++; }")
    w("    });")
    w("    group.classList.toggle('hidden', shown === 0);")
    w("  });")
    w("  const total = document.querySelectorAll('.group tbody tr').length;")
    w("  document.getElementById('result-count').textContent =")
    w("    visible + ' of ' + total + ' entries shown';")
    w("}")
    w("filterTable();")
    w("</script>")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=CONFIG, help="field configuration to render")
    ap.add_argument("--out", default=OUT, help="page to write")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if the page on disk is not what this "
                         "script would write")
    args = ap.parse_args(argv)

    cfg = load(args.config)
    text = render(cfg, config_name=os.path.relpath(args.config, ROOT))

    if args.check:
        current = ""
        if os.path.exists(args.out):
            with open(args.out, encoding="utf-8") as fh:
                current = fh.read()
        if current == text:
            print(f"{os.path.relpath(args.out, ROOT)} is up to date")
            return 0
        print(f"{os.path.relpath(args.out, ROOT)} is stale: regenerate with "
              f"`python3 {SCRIPT}`", file=sys.stderr)
        return 1

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"wrote {os.path.relpath(args.out, ROOT)}: "
          f"{len(entries(cfg))} columns, {len(patterns(cfg))} patterns")
    return 0


if __name__ == "__main__":
    sys.exit(main())
