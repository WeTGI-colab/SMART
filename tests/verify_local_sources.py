#!/usr/bin/env python3
"""
Verify SMART's annotation against the local source files it was annotated from.

Why this exists
---------------
`tests/verification1/verify.py` compares the output against the REST APIs of OncoKB, VEP,
CIViC, COSMIC and CADD. That covers 219 of the 425 described columns. The rest are annotated
from files on disk through VEP custom annotations and plugins, which have no API to query, so
they were never checked against anything:

    49  dbNSFP          44  ClinVar            10  SpliceAI
    10  Cancer Hotspots 10  AACR Project GENIE  1  REVEL
     1  LOEUF           26  gnomAD (see below)

This script closes that gap for the sources that are a transcription: it looks the variant up
in the same file the pipeline annotated from and compares the value SMART wrote. The
equivalent of the API check, against a local index.

What it does NOT cover
----------------------
**gnomAD (26 columns).** There is no gnomAD file in the pipeline. `gnomADe_AF`, `gnomADg_AF`
and the related columns come out of the VEP cache via `--everything`, so there is nothing on
disk to look up. Verifying them means either querying the VEP REST API, whose gnomAD release
may differ from the cache's and so cannot prove a mismatch wrong, or obtaining the matching
gnomAD release. Left out deliberately rather than silently.

**GENIE (10 columns).** Those are not transcriptions: `GENIE_O4_MATCH_HAEM`,
`GENIE_TOP_TUMOUR_TYPES` and the rest are aggregates computed by `scripts/genie_annotation.py`
from a counts table that holds one row per variant per OncoTree code. Verifying them means
recomputing the aggregation, which checks our arithmetic rather than our transcription. That
is a separate module and a more valuable check -- the GENIE indel left-alignment defect in the
changelog was exactly this class of error -- so it is kept apart from the lookups here.

Matching
--------
Variants are keyed on CHROM:POS:REF:ALT, taken from the MAF's `Chromosome`,
`Start_Position`, `Reference_Allele` and `Tumor_Seq_Allele2`. SMART keeps the VCF's own allele
representation in the MAF rather than converting to the MAF indel convention, so the key works
against the source files directly. Contig naming differs between sources and is declared per
module: ClinVar, SpliceAI, REVEL and dbNSFP use `1`, Cancer Hotspots uses `chr1`.

Where a position carries several rows -- dbNSFP and REVEL hold one per transcript -- a PASS
means the value SMART wrote is present among the rows the source has for that variant. It does
not assert which transcript VEP chose, because the file does not record that choice. Those
cases are noted in the output.

Usage
-----
    python3 tests/verify_local_sources.py \\
        --maf tests/verification1/output/FINAL_Table/verification1.oncoKB.maf \\
        --refs /path/to/refs \\
        --output local_results.tsv

    # one source at a time
    python3 tests/verify_local_sources.py --maf ... --refs ... --modules clinvar spliceai

Exit code 0 if nothing mismatched, 1 on any MISMATCH or ERROR.
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys
from collections import defaultdict

import pandas as pd

# ---------------------------------------------------------------- source files

SOURCES = {
    "clinvar":  dict(path="ClinVar/clinvar.vcf.gz",                        contig="plain"),
    "hotspots": dict(path="CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz", contig="chr"),
    "spliceai": dict(path="SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz",  contig="plain"),
    "spliceai_indel": dict(path="SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz",
                           contig="plain"),
    "revel":    dict(path="REVEL/new_tabbed_revel_grch38.tsv.gz",          contig="plain"),
    "loeuf":    dict(path="gnomAD_constraints/loeuf_dataset_grch38.tsv.gz", contig=None),
    "dbnsfp":   dict(path="dbNSFP/dbNSFP5.0a_grch38.smart.gz",             contig="plain"),
}

# SpliceAI's INFO is one pipe-delimited record; the plugin names the columns after it
SPLICEAI_FIELDS = ["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL",
                   "DP_AG", "DP_AL", "DP_DG", "DP_DL"]

EMPTY = ("", ".", "-", "nan", "NaN", "None", "NA", "[]")


def contig(chrom, style):
    c = str(chrom)
    if style == "chr":
        return c if c.startswith("chr") else "chr" + c
    return c[3:] if c.startswith("chr") else c


def tabix(path, region):
    """Rows at `region`, or [] if the region is empty. Raises on a real tabix failure."""
    r = subprocess.run(["tabix", path, region], capture_output=True, text=True)
    if r.returncode != 0:
        # an absent contig is not an error worth aborting for
        if "could not load" in r.stderr or "not found" in r.stderr.lower():
            return []
        raise RuntimeError(f"tabix {path} {region}: {r.stderr.strip()}")
    return [l.split("\t") for l in r.stdout.splitlines() if l and not l.startswith("#")]


# VEP rewrites the list separators of a custom-annotation value: ',' and '|' both become
# '&', because '|' delimits the CSQ field and ',' separates alleles. So a value that differs
# only in its separator is the same value, and comparing the elements as a set also absorbs
# a difference in order.
LIST_SEPARATORS = "&,|"

# A VCF flag has no value of its own. The table renders its presence as 1, and a float
# column renders that as 1.0, while this script reads it from the source as True.
TRUTHY = {"1", "1.0", "true", "yes", "present"}


def norm(v):
    """The comparison form. Mirrors what post_analysis does to a value, no more."""
    s = str(v).strip()
    if s in EMPTY:
        return ""
    s = s.replace("%3D", "=")
    if s.lower() in TRUTHY:
        return "True"
    try:                                  # 0.10 and 0.1 are the same number
        return f"{float(s):g}"
    except ValueError:
        pass
    if any(sep in s for sep in LIST_SEPARATORS):
        parts = s.replace(",", "&").replace("|", "&").split("&")
        return "&".join(sorted(p.strip() for p in parts if p.strip()))
    return s


def info_dict(info):
    d = {}
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        elif part:
            d[part] = "True"              # a flag, as vcf2table renders it
    return d


def read_maf(path, **kw):
    """Read a SMART MAF without losing values that contain '#'.

    The obvious `comment="#"` is wrong here: pandas honours a comment character anywhere
    in a line, not only at its start, so a row is silently truncated at the first '#' in
    any field and every later column comes back empty. Real annotation values contain
    '#' -- ClinVar's CLNVI holds `UniProtKB:P01111#VAR_006847` -- so that read drops the
    rest of the row. See issue #4. The leading `#SMART_VERSION` lines are skipped by
    counting them instead.
    """
    skip = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                skip += 1
            else:
                break
    return pd.read_csv(path, sep="\t", skiprows=skip, **kw)


# ---------------------------------------------------------------- the modules

def _vcf_rows(refs, source, row):
    """Source rows at this variant's position, matched on REF and ALT."""
    spec = SOURCES[source]
    path = os.path.join(refs, spec["path"])
    if not os.path.exists(path):
        return None, f"source file absent: {spec['path']}"
    c = contig(row["Chromosome"], spec["contig"])
    pos = int(row["Start_Position"])
    hits = [r for r in tabix(path, f"{c}:{pos}-{pos}")
            if len(r) > 4 and r[3] == row["Reference_Allele"]
            and row["Tumor_Seq_Allele2"] in r[4].split(",")]
    return hits, None


def check_vcf_prefixed(refs, row, cols, source, prefix):
    """ClinVar and Cancer Hotspots: column `<prefix>_<INFOKEY>` holds that INFO value."""
    hits, err = _vcf_rows(refs, source, row)
    if err:
        return [(c, "", f"SKIP: {err}") for c in cols]
    merged = {}
    for r in hits:
        merged.update(info_dict(r[7]))
    out = []
    for col in cols:
        key = col[len(prefix) + 1:] if col.startswith(prefix + "_") else None
        if key is None:
            # The bare `<prefix>` column is the source record's own ID -- ClinVar's variant
            # id, Cancer Hotspots' COSM id. Where the source record has no id, VEP falls back
            # to the location, so that spelling is accepted too.
            ids = set()
            for r in hits:
                if len(r) > 2 and r[2] not in (".", ""):
                    ids.add(r[2])
                    # a source record may carry several ';'-separated ids and VEP keeps
                    # the first, so each one on its own is an acceptable answer too
                    ids.update(r[2].split(";"))
            loc = {f"{row['Chromosome']}:{row['Start_Position']}-{row['Start_Position']}"}
            out.append((col, (ids or (loc if hits else {""})), None))
            continue
        out.append((col, merged.get(key, ""), None))
    return [(c, v, n) for c, v, n in out]


def check_spliceai(refs, row, cols):
    """The four deltas and four positions, from the pipe record."""
    source = "spliceai" if (len(row["Reference_Allele"]) == 1
                            and len(row["Tumor_Seq_Allele2"]) == 1) else "spliceai_indel"
    hits, err = _vcf_rows(refs, source, row)
    if err:
        return [(c, "", f"SKIP: {err}") for c in cols]
    records = []
    for r in hits:
        for rec in info_dict(r[7]).get("SpliceAI", "").split(","):
            f = rec.split("|")
            if len(f) == len(SPLICEAI_FIELDS):
                records.append(dict(zip(SPLICEAI_FIELDS, f)))
    out = []
    for col in cols:
        if col == "SpliceAI_cutoff":
            out.append((col, None, "SKIP: a plugin verdict, not a value in the source"))
            continue
        field = col.replace("SpliceAI_pred_", "")
        vals = {r.get(field, "") for r in records}
        out.append((col, vals, None))      # a set: any match is a pass
    return out


def check_tsv(refs, row, cols, source, header_names, pos_col, value_cols):
    """REVEL and dbNSFP: a tabix-indexed TSV keyed on position, ref and alt.

    A position can hold several rows, one per transcript, so each column yields the set of
    values the source offers and a PASS means SMART's value is among them.
    """
    spec = SOURCES[source]
    path = os.path.join(refs, spec["path"])
    if not os.path.exists(path):
        return [(c, "", f"SKIP: source file absent: {spec['path']}") for c in cols]
    c = contig(row["Chromosome"], spec["contig"])
    pos = int(row["Start_Position"])
    idx = {n: i for i, n in enumerate(header_names)}
    hits = []
    for r in tabix(path, f"{c}:{pos}-{pos}"):
        if len(r) <= max(idx.values()):
            continue
        if (r[idx[pos_col]] == str(pos) and r[idx["ref"]] == row["Reference_Allele"]
                and r[idx["alt"]] == row["Tumor_Seq_Allele2"]):
            hits.append(r)
    out = []
    for col in cols:
        if col not in idx:
            out.append((col, None, "SKIP: no column of this name in the source"))
            continue
        out.append((col, {r[idx[col]] for r in hits}, None))
    return out


def check_loeuf(refs, row, cols, table):
    """LOEUF, looked up by the transcript VEP annotated (the plugin uses match_by=transcript)."""
    if table is None:
        return [(c, "", "SKIP: source file absent") for c in cols]
    tx = str(row.get("Feature", "")).split(".")[0]
    if not tx or tx in EMPTY:
        return [(c, None, "SKIP: no Feature/transcript on this row") for c in cols]
    return [(c, table.get(tx, ""), None) for c in cols]


def load_loeuf(refs):
    path = os.path.join(refs, SOURCES["loeuf"]["path"])
    if not os.path.exists(path):
        return None
    out = {}
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
        ti, vi = hdr.index("transcript"), hdr.index("lof.oe_ci.upper")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(ti, vi):
                out[f[ti].split(".")[0]] = f[vi]
    return out


def tsv_header(refs, source):
    path = os.path.join(refs, SOURCES[source]["path"])
    if not os.path.exists(path):
        return None
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").lstrip("#").split("\t")


# ---------------------------------------------------------------- comparison

def verdict(maf_value, source_value):
    """PASS when SMART's value is what the source offers.

    `source_value` may be a set, for a position the source describes more than once.
    """
    m = norm(maf_value)
    if isinstance(source_value, set):
        cands = {norm(v) for v in source_value}
        if not cands or cands == {""}:
            return ("PASS", "both empty") if m == "" else ("MISMATCH", "source has no value")
        if m in cands:
            return ("PASS", f"one of {len(cands)} source rows" if len(cands) > 1 else "")
        return ("MISMATCH", f"source offers {sorted(cands)[:4]}")
    s = norm(source_value)
    if m == s:
        return ("PASS", "")
    return ("MISMATCH", "")


def run(maf_path, refs, modules, out_path):
    df = read_maf(maf_path, low_memory=False)
    cols = set(df.columns)
    import yaml
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    cfg = yaml.safe_load(open(os.path.join(root, "Config.yaml")))["fields"]

    def of(pred):
        return sorted(c for c in cols if pred((cfg.get(c) or {}).get("source", "")))

    groups = {
        "clinvar":  of(lambda s: s == "ClinVar"),
        "hotspots": of(lambda s: "Hotspot" in s),
        "spliceai": of(lambda s: s == "SpliceAI"),
        "revel":    of(lambda s: s == "REVEL"),
        "loeuf":    [c for c in cols if c == "LOEUF"],
        "dbnsfp":   of(lambda s: s.startswith("dbNSFP")),
    }
    # VEP's own fields, which Config.yaml attributes to ClinVar: they come from
    # --everything, not from our --custom ClinVar annotation, so they are not in
    # clinvar.vcf.gz and must not be compared against it. The misattribution is itself
    # provenance drift of the kind issue #16 is about.
    NOT_FROM_CUSTOM = {"CLIN_SIG", "PHENO", "PUBMED", "SOMATIC"}
    groups["clinvar"] = [c for c in groups["clinvar"] if c not in NOT_FROM_CUSTOM]

    loeuf_table = load_loeuf(refs) if "loeuf" in modules else None
    revel_hdr = tsv_header(refs, "revel") if "revel" in modules else None
    dbnsfp_hdr = tsv_header(refs, "dbnsfp") if "dbnsfp" in modules else None

    results = []
    for _, row in df.iterrows():
        if str(row.get("Tumor_Seq_Allele2", "")).startswith("<"):
            continue                       # structural records are in no SNV/indel source
        vid = row.get("Variant_ID") or row.get("dbSNP_RS") or "."
        base = dict(Sample=row.get("Tumor_Sample_Barcode", ""),
                    Variant_ID=vid, Gene=row.get("Hugo_Symbol", ""))
        plan = []
        if "clinvar" in modules and groups["clinvar"]:
            plan.append(("ClinVar", check_vcf_prefixed(
                refs, row, groups["clinvar"], "clinvar", "ClinVar")))
        if "hotspots" in modules and groups["hotspots"]:
            plan.append(("CancerHotspots", check_vcf_prefixed(
                refs, row, groups["hotspots"], "hotspots", "CancerHotspots")))
        if "spliceai" in modules and groups["spliceai"]:
            plan.append(("SpliceAI", check_spliceai(refs, row, groups["spliceai"])))
        if "revel" in modules and groups["revel"] and revel_hdr:
            plan.append(("REVEL", check_tsv(refs, row, groups["revel"], "revel",
                                            revel_hdr, "grch38_pos", groups["revel"])))
        if "loeuf" in modules and groups["loeuf"]:
            plan.append(("LOEUF", check_loeuf(refs, row, groups["loeuf"], loeuf_table)))
        if "dbnsfp" in modules and groups["dbnsfp"] and dbnsfp_hdr:
            plan.append(("dbNSFP", check_tsv(refs, row, groups["dbnsfp"], "dbnsfp",
                                             dbnsfp_hdr, "pos(1-based)", groups["dbnsfp"])))
        for module, checks in plan:
            for col, src, note in checks:
                if note and note.startswith("SKIP"):
                    results.append(dict(base, Module=module, MAF_Field=col,
                                        MAF_Value=row.get(col, ""), Source_Value="",
                                        Result=note, Note=""))
                    continue
                res, why = verdict(row.get(col, ""), src)
                shown = (";".join(sorted({str(x) for x in src}))[:200]
                         if isinstance(src, set) else str(src)[:200])
                results.append(dict(base, Module=module, MAF_Field=col,
                                    MAF_Value=str(row.get(col, ""))[:200],
                                    Source_Value=shown, Result=res, Note=why))

    fields = ["Module", "Sample", "Variant_ID", "Gene", "MAF_Field", "MAF_Value",
              "Source_Value", "Result", "Note"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)

    n = len(results)
    counts = defaultdict(int)
    for r in results:
        counts[r["Result"].split(":")[0]] += 1
    print(f"\nLocal source verification — {os.path.basename(maf_path)}")
    print(f"  checks: {n}")
    for k in sorted(counts):
        print(f"    {k:10} {counts[k]}")
    by_mod = defaultdict(lambda: defaultdict(int))
    for r in results:
        by_mod[r["Module"]][r["Result"].split(":")[0]] += 1
    print(f"\n  {'module':16} {'PASS':>6} {'MISMATCH':>9} {'SKIP':>6}")
    for m in sorted(by_mod):
        d = by_mod[m]
        print(f"  {m:16} {d.get('PASS', 0):>6} {d.get('MISMATCH', 0):>9} {d.get('SKIP', 0):>6}")
    bad = [r for r in results if r["Result"] == "MISMATCH"]
    if bad:
        print(f"\n  {len(bad)} mismatch(es), first 15:")
        for r in bad[:15]:
            print(f"    {r['Module']:12} {r['Variant_ID']:16} {r['MAF_Field']:26} "
                  f"MAF={r['MAF_Value'][:24]!r} source={r['Source_Value'][:30]!r} {r['Note']}")
    print(f"\n  written: {out_path}")
    return 1 if bad or counts.get("ERROR") else 0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--maf", required=True)
    p.add_argument("--refs", required=True, help="the reference directory the run used")
    p.add_argument("--output", default="local_source_results.tsv")
    p.add_argument("--modules", nargs="+",
                   default=["clinvar", "hotspots", "spliceai", "revel", "loeuf", "dbnsfp"])
    a = p.parse_args()
    return run(a.maf, a.refs, set(a.modules), a.output)


if __name__ == "__main__":
    sys.exit(main())
