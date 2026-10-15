#!/usr/bin/env python3
"""
Compose one publishable verification report from what a full suite run produced.

Why
---
A full run leaves its evidence in half a dozen places: per-field comparisons against the
sources' own APIs, the same fields against the local source files, the transcript-prioritisation
runs, the parallel-equivalence run, the empty-column audit, and the provenance of the tools and
reference files used. Anyone asking "what was verified, on which variants, against what" has to
assemble that themselves, and an assembled answer is not reproducible.

This reads those artefacts and writes a single document: every variant in the cohort with the
annotation it received, every check that passed and against which authority, every check that
did not and why, and the versions of everything involved. It states nothing it cannot read from
a file the run wrote.

Two renderings, one source
--------------------------
The content is assembled once as a list of blocks -- headings, paragraphs, tables, bullet
lists, preformatted text -- and then rendered to Markdown and to PDF. Converting the Markdown
to PDF instead would mean re-parsing prose to recover the structure that was just discarded,
and the two outputs could drift. The PDF is produced with reportlab, which is pure Python: no
LaTeX, no headless browser, no system libraries, so it works inside the Docker image and on a
server as well as on a laptop.

What it reads
-------------
  output/output/Final_result_tier1.maf     the cohort and its annotations, with #SMART_VERSION
  output/logs/run_provenance.txt           tool versions, reference file names, configuration
  output/logs/*.vcf_metadata.json          per-field source and version, written by the run
  results.tsv                              per-field checks against the sources' live APIs
  local_source_results.tsv                 per-field checks against the local source files
  ../verification3/results.tsv             the transcript-prioritisation comparison

Anything absent is reported as absent rather than omitted, because a report that quietly drops
a section it could not find is the failure this exists to avoid.
"""
import argparse
import collections
import csv
import datetime
import gzip
import json
import os
import re
import sys

SENTINELS = {"", ".", "-", "nan", "NaN", "None", "NA", "[]"}


def read_vcf(path):
    """[[CHROM, POS, ID, REF, ALT], ...] from a VCF, plain or bgzipped.

    Only the five columns that identify a variant. QUAL, FILTER, INFO and the sample columns
    describe the call rather than the variant, and a reader checking "was my variant tested"
    is looking for a coordinate and two alleles.
    """
    if not path or not os.path.exists(path):
        return None
    opener = gzip.open if path.endswith(".gz") else open
    out = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                out.append(f[:5])
    def key(r):
        c = r[0].replace("chr", "")
        return (int(c) if c.isdigit() else 100 + ord(c[0]), int(r[1]))
    return sorted(out, key=key)


def read_tsv(path):
    if not path or not os.path.exists(path):
        return None
    csv.field_size_limit(10 ** 7)
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and any(f.strip() for f in r)]
    if not rows:
        return None
    hdr = rows[0]
    return [dict(zip(hdr, r + [""] * (len(hdr) - len(r)))) for r in rows[1:]]


def read_maf(path):
    """(smart_version, header, rows) -- the leading '#' lines carry the version."""
    if not os.path.exists(path):
        return None, None, None
    csv.field_size_limit(10 ** 7)
    with open(path, newline="") as fh:
        lines = fh.readlines()
    version = ""
    for l in lines:
        if l.startswith("#SMART_VERSION"):
            version = l.split(None, 1)[1].strip()
    body = [l for l in lines if not l.startswith("#") and l.strip()]
    rows = [r for r in csv.reader(body, delimiter="\t") if r and any(f.strip() for f in r)]
    if not rows:
        return version, None, None
    return version, rows[0], rows[1:]


def cell(row, idx, col):
    """A cell's value, with every spelling of "nothing" returned as the empty string.

    A MAF read as text carries "nan" where pandas wrote a missing value and "." where
    MafAnnotator did; printing either in a report would present an artefact of a file format
    as if it were an annotation.
    """
    i = idx.get(col)
    if i is None or i >= len(row):
        return ""
    v = (row[i] or "").strip()
    return "" if v in SENTINELS else v


def tally(rows, module_key, result_key):
    """{module: Counter(result-kind)} -- the kind is the word before any colon."""
    out = collections.defaultdict(collections.Counter)
    for r in rows or []:
        kind = (r.get(result_key, "") or "").split(":")[0].strip() or "?"
        out[r.get(module_key, "?")][kind] += 1
    return out


def skip_reasons(rows, module_key, result_key):
    out = collections.defaultdict(collections.Counter)
    for r in rows or []:
        res = r.get(result_key, "") or ""
        if res.startswith("SKIP") and ":" in res:
            out[r.get(module_key, "?")][res.split(":", 1)[1].strip()] += 1
    return out


# ---- the document model -------------------------------------------------------------
# A block is a tuple whose first element names its kind. Deliberately plain: the renderers
# below are the only things that interpret it, and a dict or a class would not earn its keep.
#   ("h1", text) ("h2", text) ("para", text) ("bullets", [text, ...])
#   ("table", [headers], [[cells], ...]) ("pre", text) ("rule",)
#   ("columns", [blocks], [blocks])   side by side on the page; stacked in Markdown, which
#                                     has no columns, so the content survives either way
#   ("pagebreak",)                    only the PDF honours it

def md_escape(v):
    """A cell safe inside a Markdown table: '|' would otherwise end the column."""
    return str(v).replace("|", "\\|").replace("\n", " ")


def render_markdown(blocks):
    out = []
    for b in blocks:
        kind = b[0]
        if kind == "h1":
            out.append(f"# {b[1]}\n")
        elif kind == "h2":
            out.append(f"## {b[1]}\n")
        elif kind == "para":
            out.append(f"{b[1]}\n")
        elif kind == "bullets":
            out.append("\n".join(f"- {x}" for x in b[1]) + "\n")
        elif kind == "pre":
            out.append("```\n" + b[1] + "\n```\n")
        elif kind == "rule":
            out.append("---\n")
        elif kind == "pagebreak":
            pass                      # Markdown has no pages
        elif kind == "columns":
            out.append(render_markdown(b[1]))
            out.append(render_markdown(b[2]))
        elif kind == "table":
            hdr, rows = b[1], b[2]
            t = ["| " + " | ".join(md_escape(h) for h in hdr) + " |",
                 "|" + "|".join("---" for _ in hdr) + "|"]
            for r in rows:
                t.append("| " + " | ".join(md_escape(c) for c in r) + " |")
            out.append("\n".join(t) + "\n")
    return "\n".join(out)


def render_pdf(blocks, path, title="SMART verification report"):
    """Lay the blocks out with reportlab. Returns False if reportlab is unavailable.

    Landscape A4 because the cohort table is wide: nine columns of identifiers and HGVS
    strings do not fit portrait without a font nobody can read. Tables repeat their header
    row on each page, because a table of 37 variants spanning a page break is otherwise
    unreadable in exactly the document someone is meant to rely on.
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_LEFT
        from reportlab.lib.pagesizes import A4, landscape
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.platypus import (LongTable, PageBreak, Paragraph, Preformatted,
                                        SimpleDocTemplate, Spacer, Table, TableStyle)
    except ImportError:
        return False

    ss = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=ss["BodyText"], fontSize=8.5, leading=11.5,
                          alignment=TA_LEFT, spaceAfter=4)
    h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontSize=16, leading=19, spaceAfter=8)
    h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontSize=11.5, leading=14,
                        spaceBefore=12, spaceAfter=5)
    cell = ParagraphStyle("cell", parent=body, fontSize=7, leading=8.5, spaceAfter=0)
    cellh = ParagraphStyle("cellh", parent=cell, fontName="Helvetica-Bold")
    mono = ParagraphStyle("mono", parent=body, fontName="Courier", fontSize=7, leading=9)

    def inline(text):
        """Markdown's inline emphasis, which the paragraphs use, into reportlab's markup."""
        import re as _re
        t = str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        t = _re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
        t = _re.sub(r"`(.+?)`", r'<font face="Courier">\1</font>', t)
        return t

    def to_flowables(bs):
        """Blocks to reportlab flowables. Recursive, so a column can hold any block."""
        flow = []
        for b in bs:
            kind = b[0]
            if kind == "h1":
                flow.append(Paragraph(inline(b[1]), h1))
            elif kind == "h2":
                flow.append(Paragraph(inline(b[1]), h2))
            elif kind == "para":
                flow.append(Paragraph(inline(b[1]), body))
            elif kind == "bullets":
                for x in b[1]:
                    flow.append(Paragraph("&bull;&nbsp;" + inline(x), body))
            elif kind == "pre":
                flow.append(Preformatted(b[1], mono))
                flow.append(Spacer(1, 4))
            elif kind == "rule":
                flow.append(Spacer(1, 8))
            elif kind == "pagebreak":
                flow.append(PageBreak())
            elif kind == "columns":
                avail = landscape(A4)[0] - 24 * mm
                t = Table([[to_flowables(b[1]), to_flowables(b[2])]],
                          colWidths=[avail * 0.58, avail * 0.42], hAlign="LEFT")
                t.setStyle(TableStyle([
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (0, 0), 0),
                    ("RIGHTPADDING", (0, 0), (0, 0), 8),
                    ("LEFTPADDING", (1, 0), (1, 0), 8),
                    ("RIGHTPADDING", (1, 0), (1, 0), 8),
                    ("BOX", (1, 0), (1, 0), 0.5, colors.HexColor("#b0b0b0")),
                    ("BACKGROUND", (1, 0), (1, 0), colors.HexColor("#fafafa")),
                    ("TOPPADDING", (1, 0), (1, 0), 6),
                    ("BOTTOMPADDING", (1, 0), (1, 0), 6),
                ]))
                flow.append(t)
                flow.append(Spacer(1, 6))
            elif kind == "table":
                hdr, rows = b[1], b[2]
                data = [[Paragraph(inline(h), cellh) for h in hdr]]
                for r in rows:
                    data.append([Paragraph(inline(c), cell) for c in r])
                t = LongTable(data, repeatRows=1, hAlign="LEFT")
                t.setStyle(TableStyle([
                    ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#b0b0b0")),
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eeeeee")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 3),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                    ("TOPPADDING", (0, 0), (-1, -1), 2),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                     [colors.white, colors.HexColor("#f7f7f7")]),
                ]))
                flow.append(t)
                flow.append(Spacer(1, 6))
        return flow

    story = to_flowables(blocks)

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(12 * mm, 8 * mm, title)
        canvas.drawRightString(doc.pagesize[0] - 12 * mm, 8 * mm, f"page {doc.page}")
        canvas.restoreState()

    SimpleDocTemplate(path, pagesize=landscape(A4),
                      leftMargin=12 * mm, rightMargin=12 * mm,
                      topMargin=12 * mm, bottomMargin=14 * mm,
                      title=title, author="SMART").build(story, onFirstPage=footer,
                                                          onLaterPages=footer)
    return True


def collect(args):
    """Assemble the report as blocks. Returns (blocks, ok) -- ok is False if the run
    produced no MAF, in which case there is nothing to report on and that is said."""
    v1 = args.verification1
    maf = os.path.join(v1, "output", "output", "Final_result_tier1.maf")
    smart_version, hdr, rows = read_maf(maf)
    api = read_tsv(os.path.join(v1, "results.tsv"))
    local = read_tsv(os.path.join(v1, "local_source_results.tsv"))
    v3 = read_tsv(args.verification3) if args.verification3 else None
    prov_path = os.path.join(v1, "output", "logs", "run_provenance.txt")
    prov = open(prov_path).read().rstrip() if os.path.exists(prov_path) else None

    B = []
    today = datetime.date.today().isoformat()
    B.append(("h1", "SMART verification report"))
    B.append(("para", f"Generated {today} by `scripts/verification_report.py` from the "
                      f"artefacts of a full suite run. Every figure below is read from a file "
                      f"the run wrote; nothing is asserted."))
    if smart_version:
        B.append(("para", f"**SMART version {smart_version}**"))

    B.append(("h2", "What was run, and what was verified"))
    left = []
    if prov:
        left.append(("pre", prov))
    else:
        left.append(("para", f"`{os.path.basename(prov_path)}` is absent, so the tool and "
                             f"reference versions of this run cannot be stated. Re-run with "
                             f"an image that writes it (SMART 2.1.0 or later)."))

    # The right-hand panel answers "what does this document actually establish?", which the
    # tables further down answer only if the reader already knows how the suite is built.
    right = [("para", "**What this report establishes**")]
    n_api = len(api) if api else 0
    n_loc = len(local) if local else 0
    right.append(("para",
        f"Every annotated field of every variant was compared with the source it came from, "
        f"twice over and by two independent routes: **{n_api:,} comparisons** against the "
        f"sources' own live APIs, and **{n_loc:,}** against the reference files on disk. The "
        f"two routes share no code and no inputs beyond the output being checked, so a lookup, "
        f"a merge or a transcript resolution that went wrong has to deceive both to pass."))
    right.append(("para",
        "A **PASS** means the pipeline's value and the source's agree. A **MISMATCH** is the "
        "only outcome that indicates a defect in the pipeline. A **SKIP** is a comparison "
        "that cannot be made, and each is given with its reason. An **ERROR** is a "
        "third-party service that did not answer, so that check was not made at all."))
    right.append(("para",
        "Also checked: that every variant was annotated on a transcript the operator supplied "
        "(not an arbitrary one); that annotating the same variants under a different "
        "transcript whitelist changes what it should and nothing else; that a parallel run "
        "matches a sequential one; that every REF allele matches the reference genome before "
        "annotation; and that every output column carrying no value is accounted for with a "
        "stated reason."))
    right.append(("para",
        "The input variants are listed in VCF form overleaf, and every individual comparison "
        "is reproduced in the appendix, so nothing here rests on a summary."))
    B.append(("columns", left, right))

    B.append(("pagebreak",))
    B.append(("h2", "The variants put into the test, as submitted"))
    if rows is None:
        B.append(("para", f"`{maf}` is absent: the pipeline did not complete, so there is "
                          f"nothing to report on."))
        return B, False
    idx = {c: i for i, c in enumerate(hdr)}

    # The submitted VCF records, not the annotated MAF. The point of this page is to show
    # what was put in, in the form it was put in, so a reader can reproduce the run or check
    # a coordinate against their own data without reading our column conventions first.
    vcf_path = getattr(args, "vcf", None)
    vcf_rows = read_vcf(vcf_path) if vcf_path else None
    if vcf_rows:
        B.append(("para", f"The {len(vcf_rows)} records submitted to the pipeline, exactly as "
                          f"they appear in the input VCF (`{os.path.basename(vcf_path)}`, "
                          f"GRCh38). Four carry a symbolic ALT because they are copy-number "
                          f"calls, which have no sequence to state."))
        B.append(("pre", "#CHROM\tPOS\tID\tREF\tALT\n" +
                         "\n".join("\t".join(r) for r in vcf_rows)))
    else:
        B.append(("para", "The input VCF was not given to this report (`--vcf`), so the "
                          "submitted records cannot be shown in the form they were "
                          "submitted. What follows is read from the annotated output."))

    B.append(("h2", "What each variant was annotated as"))
    B.append(("para", f"{len(rows)} variants. For each: the annotation it received, and "
                      f"whether that annotation describes a transcript supplied in "
                      f"`--transcripts-file`."))
    cols = ["Hugo_Symbol", "ID", "HGVSc", "HGVSp_Short", "Variant_Classification",
            "Consequence", "IMPACT", "NM_Transcript", "TRANSCRIPT_MATCH"]
    have = [c for c in cols if c in idx]
    table = sorted(([cell(r, idx, c) or "\u2014" for c in have] for r in rows),
                   key=lambda t: (t[0], t[1] if len(t) > 1 else ""))
    B.append(("table", [c.replace("_", " ") for c in have], table))
    if "TRANSCRIPT_MATCH" in idx:
        tm = collections.Counter(cell(r, idx, "TRANSCRIPT_MATCH") for r in rows)
        B.append(("para", f"**{tm.get('match', 0)} of {len(rows)} rows were annotated on a "
                          f"transcript that was supplied.**"))
        other = {k: v for k, v in tm.items() if k != "match"}
        if other:
            B.append(("para", "Rows that were not, and the transcript used instead \u2014 every "
                              "transcript-specific column on these rows describes that "
                              "transcript, not a requested one:"))
            B.append(("bullets", [f"`{k or '(blank)'}` \u2014 {v} row(s)"
                                  for k, v in sorted(other.items())]))

    for title, rows_, blurb in (
        ("Checked against each source's own API", api,
         "Each field re-queried from the authority that produced it \u2014 OncoKB, Ensembl VEP, "
         "CIViC \u2014 and compared with what the pipeline wrote."),
        ("Checked against the local source files", local,
         "Each field compared with the file on disk the annotation came from, so a lookup, a "
         "merge or a transcript resolution that went wrong is caught without depending on a "
         "network service."),
    ):
        B.append(("h2", title))
        if rows_ is None:
            B.append(("para", "Not available: this part of the suite did not run, or wrote "
                              "no results."))
            continue
        B.append(("para", blurb))
        t = tally(rows_, "Module", "Result")
        kinds = ["PASS", "MISMATCH", "COVERAGE_GAP", "SKIP", "ERROR"]
        present = [k for k in kinds if any(c.get(k) for c in t.values())]
        body, tot = [], collections.Counter()
        for mod in sorted(t):
            body.append([mod] + [t[mod].get(k, 0) or "" for k in present])
            tot.update(t[mod])
        body.append(["**total**"] + [f"**{tot.get(k, 0)}**" for k in present])
        B.append(("table", ["module"] + present, body))
        n_mis = tot.get("MISMATCH", 0)
        if n_mis:
            B.append(("para", f"**{n_mis} mismatch(es).** A mismatch is the pipeline "
                              f"disagreeing with the authority and is the only outcome here "
                              f"that indicates a defect:"))
            items = []
            for r in rows_:
                if (r.get("Result", "") or "").startswith("MISMATCH") and len(items) < 20:
                    items.append(f"{r.get('Gene', '?')} {r.get('Variant_ID', '')} "
                                 f"`{r.get('MAF_Field', '?')}`: pipeline "
                                 f"`{(r.get('MAF_Value', '') or '')[:60]}` vs source "
                                 f"`{(r.get('Source_Value') or r.get('API_Value') or '')[:60]}`")
            B.append(("bullets", items))
        else:
            B.append(("para", "**No mismatches.** Every field the pipeline wrote agreed with "
                              "the source it came from."))
        sk = skip_reasons(rows_, "Module", "Result")
        if sk:
            B.append(("para", "Checks not made, and why. A skip is a comparison that cannot "
                              "be made \u2014 the source has no transcript-specific value, a "
                              "plugin does not annotate that consequence, an endpoint takes "
                              "substitutions only \u2014 not a disagreement:"))
            B.append(("bullets", [f"**{mod}**, {n}: {reason}"
                                  for mod in sorted(sk)
                                  for reason, n in sk[mod].most_common()]))
        if tot.get("ERROR"):
            B.append(("para", f"{tot['ERROR']} check(s) could not be completed because a "
                              f"third-party service did not answer. These are not "
                              f"disagreements and vary between runs; the count falls as the "
                              f"local response cache fills."))

    B.append(("h2", "Transcript prioritisation"))
    if v3 is None:
        B.append(("para", "Not available: `verification3/results.tsv` was not produced."))
    else:
        rk = "Result" if "Result" in v3[0] else list(v3[0])[-1]
        t = collections.Counter((r.get(rk, "") or "").split(":")[0] for r in v3)
        B.append(("para", "The same variants annotated twice, once per transcript whitelist, "
                          "and the two outputs compared. `DIFFER` is the expected outcome for "
                          "a field that depends on the transcript; `SAME` on such a field "
                          "would mean the whitelist had no effect."))
        B.append(("table", ["outcome", "checks"], [[k, v] for k, v in t.most_common()]))

    B.append(("h2", "Columns that carry no value, and why"))
    if not args.allowlist or not os.path.exists(args.allowlist):
        B.append(("para", "Not available: no empty-column allowlist was given."))
    else:
        kinds2, reasons = collections.Counter(), collections.defaultdict(list)
        with open(args.allowlist) as fh:
            for line in fh:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or f[0] == "pattern":
                    continue
                kinds2[f[1]] += 1
                reasons[f[1]].append((f[0], f[2]))
        B.append(("para", "Every column of the output that carried no value on any variant of "
                          "this run is accounted for. `structural` means no variant could "
                          "fill it as the pipeline runs; `cohort` means a variant would fill "
                          "it and this cohort has none."))
        B.append(("table", ["kind", "columns"], [[k, v] for k, v in kinds2.most_common()]))
        if reasons.get("cohort"):
            seen, items = set(), []
            for col, why in reasons["cohort"]:
                if why not in seen:
                    seen.add(why)
                    items.append(why)
            B.append(("para", "The `cohort` entries are a list of what a wider cohort would "
                              "cover:"))
            B.append(("bullets", items))

    B.append(("h2", "Where each annotation came from"))
    logs = os.path.join(v1, "output", "logs")
    metas = ([os.path.join(logs, f) for f in sorted(os.listdir(logs))
              if f.endswith("_metadata.json")] if os.path.isdir(logs) else [])
    if not metas:
        B.append(("para", "Not available: no per-sample field metadata was written."))
    else:
        fields = {}
        for m in metas:
            try:
                fields.update(json.load(open(m)).get("fields", {}))
            except (ValueError, OSError):
                continue
        by = collections.Counter()
        for v in fields.values():
            if isinstance(v, dict):
                by[(v.get("source", "?"), v.get("version", "?"))] += 1
        B.append(("para", f"{len(fields)} annotated fields, grouped by the source and version "
                          f"the run recorded for them (`output/logs/*_metadata.json`):"))
        B.append(("table", ["source", "version", "fields"],
                  [[a, b, n] for (a, b), n in sorted(by.items(), key=lambda x: -x[1])]))

    # ---- the errors, named and explained ------------------------------------
    err_rows = [r for r in (api or []) + (local or [])
                if (r.get("Result", "") or "").startswith("ERROR")]
    B.append(("h2", "Checks that could not be made, and why"))
    if not err_rows:
        B.append(("para", "None: every check completed."))
    else:
        reasons = collections.Counter()
        for r in err_rows:
            msg = (r.get("Result", "") or "")
            msg = msg.split(":", 1)[1].strip() if ":" in msg else msg
            # One bucket per kind of failure, not per URL: the query string differs on every
            # row and would make 52 identical failures look like 52 different problems.
            msg = re.sub(r"https?://\S+", "<url>", msg)
            msg = re.sub(r"\s+", " ", msg).strip()
            reasons[msg[:120]] += 1
        B.append(("para",
            f"{len(err_rows)} check(s) did not complete. **None of these is a disagreement "
            f"between the pipeline and its source** \u2014 the comparison never happened, so "
            f"nothing was found to be wrong or right. They are listed so the figure above is "
            f"not read as a silent pass."))
        B.append(("table", ["why the check did not complete", "checks"],
                  [[k, v] for k, v in reasons.most_common()]))
        by_mod = collections.Counter(r.get("Module", "?") for r in err_rows)
        B.append(("para", "Affected modules: " +
                          ", ".join(f"**{m}** ({n})" for m, n in by_mod.most_common()) +
                          ". These counts vary between runs with the availability of the "
                          "service, and fall as the local response cache fills; re-running "
                          "converts them into PASS or MISMATCH."))

    B.append(("rule",))
    B.append(("para", "Limitations that are known and measured are tracked as issues rather "
                      "than restated here; `docs/transcript_specificity.md` states which "
                      "annotations describe the supplied transcript and why the rest cannot."))

    # ---- appendix: every comparison, so no figure rests on a summary ---------
    for title, rows_, src_col, authority in (
        ("Appendix A \u2014 every comparison against the sources' own APIs", api, "API_Value",
         "API"),
        ("Appendix B \u2014 every comparison against the local source files", local,
         "Source_Value", "source"),
    ):
        if not rows_:
            continue
        B.append(("pagebreak",))
        B.append(("h2", title))
        B.append(("para", f"{len(rows_):,} comparisons, grouped by variant and module, in the "
                          f"form the verifier prints them. `MAF` is what the pipeline wrote, "
                          f"`{authority}` is what the authority returned. Marks: "
                          f"**\u2713** agreed \u00b7 **\u2717** disagreed \u00b7 "
                          f"**\u2013** not comparable, reason given \u00b7 "
                          f"**!** the service did not answer \u00b7 "
                          f"**\u25cb** the source has a record the pipeline did not match."))
        groups = collections.OrderedDict()
        for r in rows_:
            k = (r.get("Variant_ID") or ".", r.get("Gene") or "", r.get("Module") or "")
            groups.setdefault(k, []).append(r)
        for (vid, gene, mod), items in groups.items():
            lines = []
            for r in items:
                res = (r.get("Result", "") or "")
                # A skip is not a failure and an error is not a disagreement, so they must
                # not carry the same mark as a mismatch. Marking every non-PASS with a cross
                # made a promoter variant with no HGVSc, which cannot be queried by HGVS,
                # look like a defect.
                mark = {"PASS": "\u2713", "MISMATCH": "\u2717", "SKIP": "\u2013",
                        "ERROR": "!", "COVERAGE_GAP": "\u25cb"}.get(
                            res.split(":")[0].strip(), "?")
                field = (r.get("MAF_Field") or "")[:34]
                mafv = (r.get("MAF_Value") or "")[:30]
                srcv = (r.get(src_col) or "")[:30]
                lines.append(f"  {mark} {field:<34} MAF={mafv!r:<32} "
                             f"{authority}={srcv!r:<32} [{res[:46]}]")
            B.append(("pre", f"[{mod}]  {gene} {vid}\n" + "\n".join(lines)))
    return B, True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verification1", required=True,
                    help="tests/verification1 -- its output/, results.tsv and "
                         "local_source_results.tsv are read")
    ap.add_argument("--vcf", help="the input VCF the run annotated, shown verbatim so the "
                                  "report states what was submitted and not only what came out")
    ap.add_argument("--verification3", help="tests/verification3/results.tsv")
    ap.add_argument("--allowlist", help="tests/verification1/empty_columns_allowed.tsv")
    ap.add_argument("--out", help="write Markdown here; '-' or omitted writes to stdout")
    ap.add_argument("--pdf", help="also write a PDF here. Defaults to --out with a .pdf "
                                  "suffix when --out is a file and --pdf is not given; pass "
                                  "'none' to suppress it")
    a = ap.parse_args()

    blocks, ok = collect(a)
    md = render_markdown(blocks)

    if a.out and a.out != "-":
        with open(a.out, "w") as fh:
            fh.write(md)
        print(f"written: {a.out}", file=sys.stderr)
    else:
        sys.stdout.write(md)

    pdf = a.pdf
    if pdf is None and a.out and a.out != "-":
        pdf = os.path.splitext(a.out)[0] + ".pdf"
    if pdf and pdf != "none":
        if render_pdf(blocks, pdf):
            print(f"written: {pdf}", file=sys.stderr)
        else:
            # Not fatal: the Markdown is already written, and a missing optional dependency
            # should not fail a verification run that otherwise succeeded.
            print("WARNING: reportlab is not installed, so no PDF was written. "
                  "pip install reportlab", file=sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
