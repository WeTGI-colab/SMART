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
# ';' is added because dbNSFP packs one value per transcript into a single field with it,
# and VEP rewrites that to '&' too. A ';' cannot occur inside a VCF INFO value (it separates
# the fields) and in dbNSFP it is only ever the list separator, so treating it as one here
# cannot flatten a value that legitimately contains it.
LIST_SEPARATORS = "&,|;"

# A VCF flag has no value of its own. The source gives its presence as True while the table
# writes 1, or 1.0 in a float column. That equivalence cannot live here: a genuine score of
# 1.0 -- phastCons reaches it often -- is a number and must stay one. So only the unambiguous
# spellings are folded here, and 1 == True is decided in verdict() where both sides are known.
TRUTHY = {"true", "yes", "present"}


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
        parts = s.replace(",", "&").replace("|", "&").replace(";", "&").split("&")
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


def check_revel(refs, row, cols, header):
    """REVEL, with the transcript question made explicit.

    The source file records one transcript per score and those are older Ensembl
    identifiers, so the transcript we report is often absent from it. REVEL scores the
    amino-acid substitution rather than the transcript structure, so a score computed on
    another transcript is still the right answer **provided the substitution is the same**
    -- which the file states in its aaref/aaalt columns. The plugin runs with no_match=1
    for that reason, and this check is where the proviso is tested rather than assumed.
    """
    spec = SOURCES["revel"]
    path = os.path.join(refs, spec["path"])
    if not os.path.exists(path):
        return [(c, "", f"SKIP: source file absent: {spec['path']}") for c in cols]
    idx = {n: i for i, n in enumerate(header)}
    c = contig(row["Chromosome"], spec["contig"])
    pos = int(row["Start_Position"])
    hits = [r for r in tabix(path, f"{c}:{pos}-{pos}")
            if len(r) > max(idx.values())
            and r[idx["grch38_pos"]] == str(pos)
            and r[idx["ref"]] == row["Reference_Allele"]
            and r[idx["alt"]] == row["Tumor_Seq_Allele2"]]
    if not hits:
        return [(c2, set(), None) for c2 in cols]

    # the amino-acid change the table reports, as VEP writes it: "V/E"
    aa = str(row.get("Amino_acids", "") or "")
    want = tuple(aa.split("/")) if "/" in aa else None
    tx = str(row.get("Feature", "")).split(".")[0]

    agree = [r for r in hits
             if want is None or (r[idx["aaref"]], r[idx["aaalt"]]) == want]
    same_tx = [r for r in hits
               if tx and tx in {i.split(".")[0]
                                for i in r[idx["Ensembl_transcriptid"]].split(";")}]

    note = None
    if not same_tx and agree:
        note = (f"score computed on another transcript "
                f"({agree[0][idx['Ensembl_transcriptid']][:40]}), same substitution {aa}")
    elif want and not agree:
        return [(c2, {r[idx[c2]] for r in hits} if c2 in idx else set(),
                 None) for c2 in cols]      # let the value comparison call it

    rows_to_use = agree or hits
    return [(c2, {r[idx[c2]] for r in rows_to_use} if c2 in idx else set(), note)
            for c2 in cols]


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


# ---------------------------------------------------------------- gnomAD, via the API

# gnomAD is the one source in the unverified set with no file on disk: `gnomADe_AF`,
# `gnomADg_AF` and the population columns come out of the VEP cache via --everything, so
# there is nothing to look up locally and the Ensembl REST API is the only route. It is
# included here because it belongs to the same gap, not because it is a local source.
#
# verification1/verify.py excluded gnomAD on the grounds that the API's gnomAD release may
# differ from the cache's, which would make a mismatch unprovable. On the variants checked
# so far the two agree exactly, so the exclusion looks unnecessary -- but a disagreement
# here should be read as "these may be different releases" before it is read as an error.

VEP_REST = "https://rest.ensembl.org/vep/human/region"

# MAF column -> key in colocated_variants[].frequencies[<alt>]
#   gnomADe_AF -> gnomade        gnomADg_NFE_AF -> gnomadg_nfe        AFR_AF -> afr
def gnomad_key(col):
    stem = col[:-3] if col.endswith("_AF") else col      # drop the _AF suffix
    if stem in ("gnomADe", "gnomADg"):
        return stem.lower()
    for pre in ("gnomADe_", "gnomADg_"):
        if stem.startswith(pre):
            return pre[:-1].lower() + "_" + stem[len(pre):].lower()
    return stem.lower()                                   # the 1000 Genomes populations


def fetch_frequencies(chrom, pos, ref, alt, cache_dir, sleep=0.2):
    """`frequencies` for this ALT from the Ensembl VEP REST API, or (None, reason)."""
    import json
    import time
    import urllib.error
    import urllib.request
    if len(ref) != 1 or len(alt) != 1 or alt.startswith("<"):
        return None, "SKIP: the region endpoint takes substitutions only"
    c = chrom[3:] if chrom.startswith("chr") else chrom
    url = (f"{VEP_REST}/{c}:{pos}-{pos}/{alt}"
           f"?content-type=application/json&af=1&af_gnomade=1&af_gnomadg=1")
    cache = os.path.join(cache_dir, f"{c}_{pos}_{ref}_{alt}.json") if cache_dir else None
    if cache and os.path.exists(cache):
        payload = json.load(open(cache))
    else:
        payload, last = None, None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    payload = json.loads(r.read().decode())
                break
            except urllib.error.HTTPError as e:
                last = f"HTTP {e.code}"
                # 4xx other than 429 is a permanent refusal of this request, not congestion.
                # A 400 here usually means the variant itself is malformed -- the API rejects
                # a substitution whose ALT equals the reference base, which is how a
                # strand-flipped allele presents (see issue #19).
                if 400 <= e.code < 500 and e.code != 429:
                    return None, f"ERROR: VEP API rejected the request (HTTP {e.code})"
                # 429 means we are asking too fast; honour Retry-After when it is given
                wait = float(e.headers.get("Retry-After", 0) or 0) or (2 ** attempt)
                time.sleep(wait)
            except Exception as e:
                last = type(e).__name__
                time.sleep(2 ** attempt)
        if payload is None:
            return None, f"ERROR: VEP API unreachable after 4 tries ({last})"
        if cache:
            os.makedirs(cache_dir, exist_ok=True)
            json.dump(payload, open(cache, "w"))
        time.sleep(sleep)                   # the public endpoint is rate limited
    if not isinstance(payload, list) or not payload:
        return None, "SKIP: no VEP response"
    for col in payload[0].get("colocated_variants") or []:
        freqs = col.get("frequencies") or {}
        if alt in freqs:
            return freqs[alt], None
    return {}, None                         # answered, and the variant has no frequencies


def check_gnomad(row, cols, cache_dir, sleep):
    freqs, err = fetch_frequencies(str(row["Chromosome"]), int(row["Start_Position"]),
                                   row["Reference_Allele"], row["Tumor_Seq_Allele2"],
                                   cache_dir, sleep)
    if err:
        return [(c, "", err) for c in cols]
    return [(c, freqs.get(gnomad_key(c), ""), None) for c in cols]


# ---------------------------------------------------------------- GENIE, by recomputation

# GENIE's columns are not transcriptions: genie_annotation.py aggregates a counts table that
# holds one row per variant per set of tumour types. So verifying them means recomputing the
# aggregate, which checks the arithmetic rather than the copying -- a stronger check, and the
# one that matters, since the GENIE indel left-alignment defect in the changelog was exactly
# an arithmetic error.
#
# What is recomputed independently and what is only checked for consistency:
#
#   GENIE_SAME_NUCLEOTIDE_{ALL,HAEM,SOLID}   recomputed from the table for substitutions.
#   GENIE_O4_MATCH_*, GENIE_TOP_TUMOUR_TYPES  invariants only. Reproducing them means
#       reimplementing the O4 matching rules (same amino-acid change, truncating at or
#       downstream of a residue, in-frame containment), and a verifier that re-ran the
#       pipeline's own logic would agree with it by construction and prove nothing.
#   indels   skipped for the recomputation: matching them needs the same left-alignment the
#       pipeline applies, so an independent attempt would differ for reasons of its own.

GENIE_HAEM = {"Myeloid", "Lymphoid"}


def load_genie_counts(refs):
    """The counts table as {(chrom, start, end, ref, alt): [(n, [lineages]), ...]}."""
    path = os.path.join(refs, "GENIE", "genie_o4_counts.tsv.gz")
    if not os.path.exists(path):
        return None
    out = {}
    with gzip.open(path, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            key = (r["chrom"], r["start"], r["end"], r["ref"], r["alt"])
            lineages = [x for x in (r.get("lineage") or "").split("|") if x]
            out.setdefault(key, []).append((int(r["n_patients"]), lineages))
    return out


def check_genie(row, cols, counts):
    """Recompute the same-nucleotide counts; check the rest for internal consistency."""
    if counts is None:
        return [(c, "", "SKIP: GENIE counts table absent") for c in cols]
    ref, alt = row["Reference_Allele"], row["Tumor_Seq_Allele2"]
    out = []

    if len(ref) == 1 and len(alt) == 1 and not alt.startswith("<"):
        chrom = str(row["Chromosome"])
        chrom = chrom if chrom.startswith("chr") else "chr" + chrom
        pos = str(int(row["Start_Position"]))
        entries = counts.get((chrom, pos, pos, ref, alt), [])
        all_n = sum(n for n, _ in entries)
        haem = sum(n for n, ls in entries if any(l in GENIE_HAEM for l in ls))
        solid = sum(n for n, ls in entries if any(l == "Solid" for l in ls))
        recomputed = {"GENIE_SAME_NUCLEOTIDE_ALL": all_n,
                      "GENIE_SAME_NUCLEOTIDE_HAEM": haem,
                      "GENIE_SAME_NUCLEOTIDE_SOLID": solid}
    else:
        recomputed = {}

    def num(col):
        """The column as a number, or None when it holds no value.

        'nan' must not become float('nan'): every comparison against it is False, so an
        unannotated row would fail an invariant it was never subject to.
        """
        s = str(row.get(col, "")).strip()
        if s in EMPTY:
            return None
        try:
            v = float(s)
        except ValueError:
            return None
        return None if v != v else v

    for col in cols:
        if col in recomputed:
            out.append((col, str(recomputed[col]), None))
            continue
        if col in ("GENIE_SAME_NUCLEOTIDE_ALL", "GENIE_SAME_NUCLEOTIDE_HAEM",
                   "GENIE_SAME_NUCLEOTIDE_SOLID"):
            out.append((col, None, "SKIP: indel, needs the pipeline's own left-alignment"))
            continue
        if col in ("GENIE_O4_MATCH_HAEM", "GENIE_O4_MATCH_SOLID"):
            # an invariant rather than a value: a subset cannot exceed the whole
            whole, part = num("GENIE_O4_MATCH_ALL"), num(col)
            if whole is None or part is None:
                out.append((col, None, "SKIP: no O4 match on this row, nothing to bound"))
            else:
                out.append((col, str(part) if part <= whole else f"<= {whole}", None))
            continue
        out.append((col, None, "SKIP: derived by the O4 matching rules, not recomputed here"))
    return out


# ---------------------------------------------------------------- comparison

# The VEP cache writes an allele frequency to three significant figures while the REST API
# gives four, so 1.16e-05 and 1.159e-05 are the same number recorded twice. Compare numbers
# within a relative tolerance rather than as strings.
REL_TOLERANCE = 5e-3


def _same_number(a, b):
    try:
        x, y = float(a), float(b)
    except (TypeError, ValueError):
        return False
    if x == y:
        return True
    scale = max(abs(x), abs(y))
    return scale > 0 and abs(x - y) / scale <= REL_TOLERANCE


def _flag_equal(a, b):
    """One side is a flag's presence, the other the 1 the table writes for it."""
    pair = {a, b}
    return "True" in pair and ("1" in pair or "1.0" in pair)


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
        if any(_same_number(m, c) for c in cands):
            return ("PASS", "equal within rounding")
        if any(_flag_equal(m, c) for c in cands):
            return ("PASS", "flag presence")
        return ("MISMATCH", f"source offers {sorted(cands)[:4]}")
    s = norm(source_value)
    if m == s:
        return ("PASS", "")
    if _same_number(m, s):
        return ("PASS", "equal within rounding")
    if _flag_equal(m, s):
        return ("PASS", "flag presence")
    return ("MISMATCH", "")


def run(maf_path, refs, modules, out_path, cache_dir=None, sleep=0.2):
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
        "gnomad":   of(lambda s: s.startswith("gnomAD")),
        "genie":    sorted(c for c in cols if c.startswith("GENIE_")),
    }
    # VEP's own fields, which Config.yaml attributes to ClinVar: they come from
    # --everything, not from our --custom ClinVar annotation, so they are not in
    # clinvar.vcf.gz and must not be compared against it. The misattribution is itself
    # provenance drift of the kind issue #16 is about.
    NOT_FROM_CUSTOM = {"CLIN_SIG", "PHENO", "PUBMED", "SOMATIC"}
    groups["clinvar"] = [c for c in groups["clinvar"] if c not in NOT_FROM_CUSTOM]

    loeuf_table = load_loeuf(refs) if "loeuf" in modules else None
    genie_counts = load_genie_counts(refs) if "genie" in modules else None
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
            plan.append(("REVEL", check_revel(refs, row, groups["revel"], revel_hdr)))
        if "loeuf" in modules and groups["loeuf"]:
            plan.append(("LOEUF", check_loeuf(refs, row, groups["loeuf"], loeuf_table)))
        if "genie" in modules and groups["genie"]:
            plan.append(("GENIE", check_genie(row, groups["genie"], genie_counts)))
        if "gnomad" in modules and groups["gnomad"]:
            plan.append(("gnomAD", check_gnomad(row, groups["gnomad"], cache_dir, sleep)))
        if "dbnsfp" in modules and groups["dbnsfp"] and dbnsfp_hdr:
            plan.append(("dbNSFP", check_tsv(refs, row, groups["dbnsfp"], "dbnsfp",
                                             dbnsfp_hdr, "pos(1-based)", groups["dbnsfp"])))
        for module, checks in plan:
            for col, src, note in checks:
                if note and (note.startswith("SKIP") or note.startswith("ERROR")):
                    results.append(dict(base, Module=module, MAF_Field=col,
                                        MAF_Value=row.get(col, ""), Source_Value="",
                                        Result=note, Note=""))
                    continue
                res, why = verdict(row.get(col, ""), src)
                if note:
                    why = f"{note}; {why}" if why else note
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
                   default=["clinvar", "hotspots", "spliceai", "revel", "loeuf", "dbnsfp",
                            "gnomad", "genie"])
    p.add_argument("--cache", default=None,
                   help="directory for cached VEP API responses, so a re-run costs nothing")
    p.add_argument("--sleep", type=float, default=0.2,
                   help="pause between VEP API calls (the public endpoint is rate limited)")
    a = p.parse_args()
    return run(a.maf, a.refs, set(a.modules), a.output, a.cache, a.sleep)


if __name__ == "__main__":
    sys.exit(main())
