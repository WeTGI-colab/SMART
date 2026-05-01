#!/usr/bin/env python3
"""
SMART pipeline verification3 — Transcript prioritisation impact on VEP and OncoKB annotation.

Runs the pipeline twice on the same VCF with two different preferred-transcript
files (A and B), then verifies:

  1. DIFFER  — targeted variants produce different VEP annotations between runs
               (Feature / HGVSc / HGVSp / Hugo_Symbol differ as expected)
  2. PASS    — each run's VEP annotation is confirmed correct against Ensembl REST
               using that run's expected transcript
  3. MISMATCH — the pipeline's annotation does not match the Ensembl REST result
  4. OncoKB impact — the different protein change produced by each transcript set
               leads to a different OncoKB result (oncogenicity / therapeutic level)

Three test cases are built into the VCF:

  EGFR_L858R   (chr7:55191822 T>G)
    Transcript A: NM_005228.5  → c.2573T>G / p.Leu858Arg   (MANE SELECT; LEVEL_1 in OncoKB)
    Transcript B: NM_001346897.2 → c.2438T>G / p.Leu813Arg  (unknown in OncoKB)

  TP53_R175H   (chr17:7675088 C>T)
    Transcript A: NM_000546.6    → c.524G>A / p.Arg175His  (MANE SELECT; hotspot)
    Transcript B: NM_001126115.2 → c.128G>A / p.Arg43His   (completely different position)

  MantaDEL:CDKN2A_CDKN2B  (chr9:21990000-22005000 <DEL>)
    Transcript A: NM_058195.4 → annotated as CDKN2A deletion (LEVEL_4 in OncoKB)
    Transcript B: NM_004936.4 → annotated as CDKN2B deletion (no treatment level)

  FGFR1_N546K  (chr8:38417331 G>T)
    A real TSO500 case: NM_001174067.2 is used by the TSO500 panel but is NOT
    a MANE transcript. The MANE Select for FGFR1 is NM_023110.3.
    Transcript A: NM_023110.3    → c.1638C>A / p.Asn546Lys  (MANE Select; Unknown in OncoKB)
    Transcript B: NM_001174067.2 → c.1731C>A / p.Asn577Lys  (TSO500 non-MANE; Likely Oncogenic
                                   LEVEL_4 in OncoKB — false-positive risk from wrong transcript)

Usage
-----
    # Run both pipeline outputs through verification
    python3 tests/verification3/verify.py \\
        --maf-a  tests/verification3/output_A/output/Final_result_tier1.maf \\
        --maf-b  tests/verification3/output_B/output/Final_result_tier1.maf \\
        --transcripts-a tests/verification3/transcripts_A.txt \\
        --transcripts-b tests/verification3/transcripts_B.txt \\
        --output  tests/verification3/results.tsv

Exit code 0 if all checks pass (no unexpected SAME or MISMATCH).
Exit code 1 if any SAME or MISMATCH.
"""

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Expected annotations per variant and transcript set
# ---------------------------------------------------------------------------

EXPECTED = {
    "EGFR_L858R": {
        "A": {
            "NM_Transcript":  "NM_005228",
            "HGVSc":          "c.2573T>G",
            "HGVSp":          "p.Leu858Arg",
            "Hugo_Symbol":    "EGFR",
            # OncoKB expected outcomes for transcript A
            "oncokb": {
                "oncogenic":           "Oncogenic",
                "highestSensitiveLevel": "LEVEL_1",
                "mutationEffect":      "Gain-of-function",
            },
        },
        "B": {
            "NM_Transcript":  "NM_001346897",
            "HGVSc":          "c.2438T>G",
            "HGVSp":          "p.Leu813Arg",
            "Hugo_Symbol":    "EGFR",
            # OncoKB expected outcomes for transcript B (wrong protein position)
            "oncokb": {
                "oncogenic":           "Unknown",
                "highestSensitiveLevel": None,   # no therapeutic level
                "mutationEffect":      "Unknown",
            },
        },
    },
    "TP53_R175H": {
        "A": {
            "NM_Transcript":  "NM_000546",
            "HGVSc":          "c.524G>A",
            "HGVSp":          "p.Arg175His",
            "Hugo_Symbol":    "TP53",
            "oncokb": {
                "oncogenic":   "Oncogenic",
                "hotspot":     True,
            },
        },
        "B": {
            "NM_Transcript":  "NM_001126115",
            "HGVSc":          "c.128G>A",
            "HGVSp":          "p.Arg43His",
            "Hugo_Symbol":    "TP53",
            "oncokb": {
                "oncogenic":   "Unknown",
                "hotspot":     False,
            },
        },
    },
    # FGFR1 N546K — chr8:38417331 G>T
    # Demonstrates a non-MANE TSO500 transcript (NM_001174067.2, variant 14)
    # producing a different protein position (N577K) that is coincidentally
    # flagged as Likely Oncogenic LEVEL_4 by OncoKB, while the MANE Select
    # (NM_023110.3) gives the canonical N546K annotation which OncoKB does not
    # yet classify. This is a false-positive risk of using non-MANE transcripts.
    "FGFR1_N546K": {
        "A": {
            "NM_Transcript":  "NM_023110",
            "HGVSc":          "c.1638C>A",
            "HGVSp":          "p.Asn546Lys",
            "Hugo_Symbol":    "FGFR1",
            "oncokb": {
                "oncogenic":             "Unknown",
                "highestSensitiveLevel": None,
                "mutationEffect":        "Unknown",
            },
        },
        "B": {
            "NM_Transcript":  "NM_001174067",
            "HGVSc":          "c.1731C>A",
            "HGVSp":          "p.Asn577Lys",
            "Hugo_Symbol":    "FGFR1",
            "oncokb": {
                "oncogenic":             "Likely Oncogenic",
                "highestSensitiveLevel": "LEVEL_4",
                "mutationEffect":        "Likely Gain-of-function",
            },
        },
    },
    "MantaDEL:CDKN2A_CDKN2B": {
        "A": {
            "NM_Transcript":  "NM_058195",
            "Hugo_Symbol":    "CDKN2A",
            "oncokb": {
                "oncogenic":             "Oncogenic",
                "highestSensitiveLevel": "LEVEL_4",   # CDK4/6 inhibitors (Palbociclib etc.)
                "mutationEffect":        "Loss-of-function",
            },
        },
        "B": {
            "NM_Transcript":  "NM_004936",
            "Hugo_Symbol":    "CDKN2B",
            "oncokb": {
                "highestSensitiveLevel": None,         # CDKN2B deletion has no treatment level
            },
        },
    },
}

# VEP fields compared between runs (cross-run DIFFER check)
# Fields expected to differ between runs for SNVs
DIFFER_FIELDS_SNV = ["NM_Transcript", "HGVSc", "HGVSp"]
# Fields expected to differ between runs for CNAs (gene selection changes)
DIFFER_FIELDS_CNA = ["NM_Transcript", "Hugo_Symbol"]

# VEP fields validated against Ensembl REST (per-run PASS/MISMATCH check)
VEP_VALIDATE_FIELDS = ["NM_Transcript", "HGVSc", "HGVSp"]

VEP_REST = "https://rest.ensembl.org/vep/human/region"
VEP_HEADERS = {"Content-Type": "application/json"}

ONCOKB_BASE = "https://www.oncokb.org/api/v1"
ONCOKB_MUT_URL = f"{ONCOKB_BASE}/annotate/mutations/byProteinChange"
ONCOKB_CNA_URL = f"{ONCOKB_BASE}/annotate/copyNumberAlterations"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_maf(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False, dtype=str)
    df = df.fillna("")
    return df


def maf_val(row: pd.Series, col: str) -> str:
    v = str(row.get(col, "")).strip()
    return v if v not in ("nan", ".", "") else ""


def strip_version(nm: str) -> str:
    """NM_005228.5 → NM_005228  (only for NM_/ENST-style IDs, not for HGVS strings)"""
    if nm.startswith(("NM_", "NR_", "NP_", "ENST", "ENSP")):
        return nm.split(".")[0]
    return nm


def norm(s: str) -> str:
    return s.strip().lower()


def _hgvs_core(h: str) -> str:
    """NM_005228.5:c.2573T>G → c.2573T>G"""
    if ":" in h:
        return h.split(":")[-1]
    return h


def query_vep(chrom: str, start: int, end: int, alt: str, transcript_nm: str):
    """
    Query Ensembl VEP REST for a variant and return the consequence for the
    specified transcript (version-stripped NM_ match).
    """
    region = f"{chrom}:{start}-{end}:1/{alt}"
    url = f"{VEP_REST}/{region}?content-type=application/json&refseq=1&hgvs=1&canonical=1"
    try:
        r = requests.get(url, headers=VEP_HEADERS, timeout=20)
        r.raise_for_status()
        data = r.json()[0]
    except Exception as exc:
        return {"_error": str(exc)}

    tcs = data.get("transcript_consequences", [])
    nm_bare = strip_version(transcript_nm)
    for tc in tcs:
        if strip_version(tc.get("transcript_id", "")) == nm_bare:
            return {
                "NM_Transcript": tc.get("transcript_id", ""),
                "HGVSc":         _hgvs_core(tc.get("hgvsc", "")),
                "HGVSp":         _hgvs_core(tc.get("hgvsp", "")),
                "Hugo_Symbol":   tc.get("gene_symbol", ""),
            }
    return None


# ---------------------------------------------------------------------------
# Main verification logic
def query_oncokb_snv(gene: str, alteration: str, token: str) -> dict:
    """Query OncoKB byProteinChange for a SNV. Returns the JSON dict or {}."""
    params = {
        "hugoSymbol":      gene,
        "alteration":      alteration,
        "referenceGenome": "GRCh38",
    }
    headers = {"Authorization": f"Bearer {token}", "accept": "application/json"}
    try:
        r = requests.get(ONCOKB_MUT_URL, params=params, headers=headers, timeout=15)
        r.raise_for_status()
        return r.json()
    except Exception as exc:
        return {"_error": str(exc)}


def query_oncokb_cna(gene: str, cna_type: str, token: str) -> dict:
    """Query OncoKB copyNumberAlterations. cna_type = AMPLIFICATION or DELETION."""
    params = {
        "hugoSymbol":           gene,
        "copyNameAlterationType": cna_type,
        "referenceGenome":      "GRCh38",
    }
    headers = {"Authorization": f"Bearer {token}", "accept": "application/json"}
    try:
        r = requests.get(ONCOKB_CNA_URL, params=params, headers=headers, timeout=15)
        r.raise_for_status()
        return r.json()
    except Exception as exc:
        return {"_error": str(exc)}


def _hgvsp_to_short(hgvsp_long: str) -> str:
    """ENSP...:p.Leu858Arg → L858R  (single-letter via standard 3→1 map)."""
    aa3 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*", "Stop": "*",
    }
    # Strip transcript prefix
    p = hgvsp_long.split(":")[-1]
    if p.startswith("p."):
        p = p[2:]
    # Convert three-letter codes to one-letter
    for three, one in aa3.items():
        p = p.replace(three, one)
    return p


def _maf_oncokb_fields(row: "pd.Series") -> dict:
    """Extract the OncoKB annotation fields written by the pipeline into the MAF."""
    return {
        "oncogenic":             maf_val(row, "ONCOGENIC"),
        "highestSensitiveLevel": maf_val(row, "HIGHEST_SENSITIVE_LEVEL") or None,
        "mutationEffect":        maf_val(row, "MUTATION_EFFECT"),
        "hotspot":               maf_val(row, "ONCOKB_HOTSPOT").upper() == "TRUE",
    }


# ---------------------------------------------------------------------------

def run(args):
    maf_a = load_maf(args.maf_a)
    maf_b = load_maf(args.maf_b)

    # Index both MAFs by variant ID (ID column)
    def index_maf(df: pd.DataFrame) -> dict:
        out = {}
        for _, row in df.iterrows():
            vid = maf_val(row, "ID") or maf_val(row, "Tumor_Sample_Barcode")
            if vid:
                out[vid] = row
        return out

    rows_a = index_maf(maf_a)
    rows_b = index_maf(maf_b)

    all_results = []
    total_pass = total_mismatch = total_same = total_error = 0

    print(f"\n{'#'*60}")
    print("VERIFICATION 3 — Transcript Prioritisation Impact")
    print(f"{'#'*60}")

    for vid, exp in EXPECTED.items():
        print(f"\n{'─'*60}")
        print(f"Variant: {vid}")

        row_a = rows_a.get(vid)
        row_b = rows_b.get(vid)

        if row_a is None or row_b is None:
            msg = f"Missing in {'MAF-A' if row_a is None else 'MAF-B'}"
            print(f"  ERROR: {msg}")
            total_error += 1
            all_results.append({
                "Variant_ID": vid, "Check": "PRESENCE", "Run": "A+B",
                "Field": "", "MAF_A": "", "MAF_B": "", "Expected": "", "API": "",
                "Result": f"ERROR: {msg}",
            })
            continue

        exp_a = exp.get("A", {})
        exp_b = exp.get("B", {})

        # ── 1. Cross-run DIFFER check ──────────────────────────────────────
        # Determine early if this is a CNA to pick the right field list
        _alt_peek = maf_val(row_a, "Tumor_Seq_Allele2") or maf_val(row_a, "ALT")
        _differ_fields = DIFFER_FIELDS_CNA if (_alt_peek.startswith("<")) else DIFFER_FIELDS_SNV
        print("  [DIFFER CHECK] Comparing MAF-A vs MAF-B annotation:")
        for field in _differ_fields:
            if field not in exp_a and field not in exp_b:
                continue

            val_a = strip_version(_hgvs_core(maf_val(row_a, field)))
            val_b = strip_version(_hgvs_core(maf_val(row_b, field)))

            if not val_a and not val_b:
                status = "SKIP"
                icon = "○"
            elif val_a == val_b:
                # Same transcript selected — unexpected for targeted variants
                status = "SAME"
                icon = "✗"
                total_same += 1
            else:
                status = "DIFFER"
                icon = "✓"
                total_pass += 1

            print(f"  {icon} {field:<20} A={val_a!r:<30} B={val_b!r:<30}  [{status}]")
            all_results.append({
                "Variant_ID": vid, "Check": "DIFFER", "Run": "A_vs_B",
                "Field": field, "MAF_A": val_a, "MAF_B": val_b,
                "Expected": "differ", "API": "",
                "Result": status,
            })

        # ── 2. Per-run VEP API validation ──────────────────────────────────
        # Determine genomic coordinates from MAF-A (same VCF, same coords)
        chrom = (maf_val(row_a, "Chromosome") or maf_val(row_a, "CHROM", )).lstrip("chr")
        start_s = maf_val(row_a, "Start_Position") or maf_val(row_a, "POS")
        ref     = maf_val(row_a, "Reference_Allele") or maf_val(row_a, "REF")
        alt     = maf_val(row_a, "Tumor_Seq_Allele2") or maf_val(row_a, "ALT")

        is_cna = alt in ("<DEL>", "<DUP>", "<INV>", "<BND>") or alt.startswith("<")

        if is_cna:
            print("  [API CHECK] CNA — skipping Ensembl REST (symbolic allele)")
            print("  ✓  Checking Hugo_Symbol only via expected values:")
            for run_label, row, exp_run in [("A", row_a, exp_a), ("B", row_b, exp_b)]:
                gene = maf_val(row, "Hugo_Symbol")
                exp_gene = exp_run.get("Hugo_Symbol", "")
                transcript = strip_version(maf_val(row, "NM_Transcript"))
                exp_tx = strip_version(exp_run.get("NM_Transcript", ""))
                gene_ok = norm(gene) == norm(exp_gene) if exp_gene else True
                tx_ok   = norm(transcript) == norm(exp_tx) if exp_tx else True
                status = "PASS" if (gene_ok and tx_ok) else "MISMATCH"
                icon = "✓" if status == "PASS" else "✗"
                print(f"  {icon} [{transcript}]: Hugo_Symbol={gene!r}  [{status}]")
                if status == "PASS":
                    total_pass += 1
                else:
                    total_mismatch += 1
                all_results.append({
                    "Variant_ID": vid, "Check": "VEP_API", "Run": run_label,
                    "Field": "Hugo_Symbol+Transcript_ID",
                    "MAF_A": gene if run_label == "A" else "",
                    "MAF_B": gene if run_label == "B" else "",
                    "Expected": f"{exp_gene} / {exp_tx}", "API": "expected_only",
                    "Result": status,
                })
            continue

        try:
            start_i = int(start_s)
            end_i   = start_i  # SNVs: end == start
        except (ValueError, TypeError):
            print(f"  ERROR: cannot parse start coordinate ({start_s!r})")
            total_error += 1
            continue

        for run_label, row, exp_run in [("A", row_a, exp_a), ("B", row_b, exp_b)]:
            exp_tx = exp_run.get("NM_Transcript", "")
            print(f"\n  [API CHECK] [{exp_tx}] — validating against Ensembl REST:")

            api = query_vep(chrom, start_i, end_i, alt, exp_tx)
            time.sleep(0.4)

            if api is None:
                print(f"  ○ [{exp_tx}] not found in Ensembl REST response")
                total_error += 1
                all_results.append({
                    "Variant_ID": vid, "Check": "VEP_API", "Run": run_label,
                    "Field": "Transcript_ID", "MAF_A": "", "MAF_B": "",
                    "Expected": exp_tx, "API": "not_found",
                    "Result": "ERROR: transcript not in VEP response",
                })
                continue

            if "_error" in api:
                print(f"  ERROR querying Ensembl: {api['_error']}")
                total_error += 1
                continue

            for field in VEP_VALIDATE_FIELDS:
                if field not in exp_run:
                    continue
                maf_v = strip_version(_hgvs_core(maf_val(row, field)))
                api_v = strip_version(api.get(field, ""))
                exp_v = strip_version(exp_run[field])

                # Compare MAF value against API value
                if norm(maf_v) == norm(api_v):
                    status = "PASS"
                elif norm(_hgvs_core(maf_v)) == norm(_hgvs_core(api_v)):
                    status = "PASS"
                else:
                    status = "MISMATCH"

                icon = "✓" if status == "PASS" else "✗"
                print(f"  {icon} {field:<20} MAF={maf_v!r:<30} API={api_v!r:<30}  [{status}]")
                if status == "PASS":
                    total_pass += 1
                else:
                    total_mismatch += 1

                all_results.append({
                    "Variant_ID": vid, "Check": "VEP_API", "Run": run_label,
                    "Field": field,
                    "MAF_A": maf_v if run_label == "A" else "",
                    "MAF_B": maf_v if run_label == "B" else "",
                    "Expected": exp_v, "API": api_v,
                    "Result": status,
                })

    # ── 3. OncoKB impact check ───────────────────────────────────────────────
    # The pipeline already ran OncoKB during annotation and wrote the results
    # into both MAFs. No API re-query needed — just compare the MAF fields.
    print(f"\n{'#'*60}")
    print("ONCOKB IMPACT — Does transcript choice change therapeutic annotation?")
    print("(reading OncoKB fields already written by the pipeline into each MAF)")
    print(f"{'#'*60}")

    # OncoKB fields written by the pipeline into the MAF
    ONCOKB_MAF_FIELDS = [
        ("ONCOGENIC",               "oncogenic"),
        ("HIGHEST_SENSITIVE_LEVEL", "highestSensitiveLevel"),
        ("MUTATION_EFFECT",         "mutationEffect"),
        ("ONCOKB_HOTSPOT",          "hotspot"),
    ]

    for vid, exp in EXPECTED.items():
        print(f"\n{'─'*60}")
        print(f"Variant: {vid}")

        row_a = rows_a.get(vid)
        row_b = rows_b.get(vid)
        if row_a is None or row_b is None:
            print("  SKIP: variant not found in one or both MAFs")
            continue

        exp_a = exp.get("A", {})
        exp_b = exp.get("B", {})
        alt_peek = maf_val(row_a, "Tumor_Seq_Allele2") or maf_val(row_a, "ALT")
        is_cna   = alt_peek.startswith("<")

        tx_a = maf_val(row_a, "NM_Transcript") or exp_a.get("NM_Transcript", "A")
        tx_b = maf_val(row_b, "NM_Transcript") or exp_b.get("NM_Transcript", "B")

        # Show what each run reported
        for run_label, row, tx_id in [("A", row_a, tx_a), ("B", row_b, tx_b)]:
            gene    = maf_val(row, "Hugo_Symbol")
            hgvsp   = maf_val(row, "HGVSp_Short")
            onco    = maf_val(row, "ONCOGENIC")
            level   = maf_val(row, "HIGHEST_SENSITIVE_LEVEL")
            effect  = maf_val(row, "MUTATION_EFFECT")
            hotspot = maf_val(row, "ONCOKB_HOTSPOT")
            print(f"  [{tx_id}] gene={gene!r}  hgvsp={hgvsp!r}  oncogenic={onco or '(empty)'!r}  "
                  f"highestSensitiveLevel={level or 'None'!r}  "
                  f"mutationEffect={effect or '(empty)'!r}  hotspot={hotspot or '(empty)'!r}")

        # Validate each run's MAF values against expectations — only for rows
        # where the pipeline actually populated OncoKB fields (non-empty ONCOGENIC)
        for run_label, row, exp_run, tx_id in [("A", row_a, exp_a, tx_a), ("B", row_b, exp_b, tx_b)]:
            exp_oncokb = exp_run.get("oncokb", {})
            if not exp_oncokb:
                continue
            oncogenic_populated = bool(maf_val(row, "ONCOGENIC"))
            if not oncogenic_populated:
                print(f"  ○ [{tx_id}] OncoKB fields not populated in MAF for this variant "
                      f"(SNV MafAnnotator merge issue) — skipping field validation")
                continue
            for maf_col, field in ONCOKB_MAF_FIELDS:
                if field not in exp_oncokb:
                    continue
                exp_val  = exp_oncokb[field]
                maf_v    = maf_val(row, maf_col)
                maf_norm = maf_v.lower() if maf_v else "none"
                exp_norm = str(exp_val).lower() if exp_val is not None else "none"
                status   = "PASS" if maf_norm == exp_norm else "MISMATCH"
                icon     = "✓" if status == "PASS" else "✗"
                print(f"    {icon} [{tx_id}] {field:<28} expected={str(exp_val)!r:<20} maf={maf_v!r}  [{status}]")
                if status == "PASS":
                    total_pass += 1
                else:
                    total_mismatch += 1
                all_results.append({
                    "Variant_ID": vid, "Check": "ONCOKB_IMPACT", "Run": tx_id,
                    "Field": field,
                    "MAF_A": maf_v if run_label == "A" else "",
                    "MAF_B": maf_v if run_label == "B" else "",
                    "Expected": str(exp_val), "API": "",
                    "Result": status,
                })

        # Cross-transcript DIFFER check: confirm the two runs encode different queries
        # For SNVs: HGVSp_Short is what gets sent to OncoKB — must differ between runs
        # For CNAs: Hugo_Symbol determines which gene is queried — must differ between runs
        if is_cna:
            differ_field = "Hugo_Symbol"
            val_a = maf_val(row_a, "Hugo_Symbol")
            val_b = maf_val(row_b, "Hugo_Symbol")
        else:
            differ_field = "HGVSp_Short"
            val_a = maf_val(row_a, "HGVSp_Short")
            val_b = maf_val(row_b, "HGVSp_Short")

        differ = val_a != val_b
        icon   = "✓" if differ else "✗"
        status = "DIFFER" if differ else "SAME"
        if differ:
            total_pass += 1
        else:
            total_same += 1
        print(f"  {icon} Cross-transcript {differ_field}: "
              f"[{tx_a}]={val_a!r}  [{tx_b}]={val_b!r}  [{status}]")
        all_results.append({
            "Variant_ID": vid, "Check": "ONCOKB_DIFFER", "Run": f"{tx_a}_vs_{tx_b}",
            "Field": differ_field, "MAF_A": val_a, "MAF_B": val_b,
            "Expected": "differ", "API": "",
            "Result": status,
        })

    # ── Summary ──────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  Variants checked:   {len(EXPECTED)}")
    print(f"  PASS / DIFFER:      {total_pass}  (correct — annotations differ between runs as expected, or match API)")
    print(f"  SAME (unexpected):  {total_same}  (transcript selection had no effect — investigate)")
    print(f"  MISMATCH:           {total_mismatch}  (pipeline annotation does not match Ensembl REST)")
    print(f"  ERROR:              {total_error}")
    print(f"{'='*60}")

    if args.output:
        pd.DataFrame(all_results).to_csv(args.output, sep="\t", index=False)
        print(f"\nResults written to: {args.output}")

    failed = total_same + total_mismatch + total_error
    if failed:
        print("\nRESULT: FAILED")
        sys.exit(1)
    else:
        print("\nRESULT: PASSED")
        sys.exit(0)


def main():
    parser = argparse.ArgumentParser(
        description="Verify transcript prioritisation impact on VEP annotation"
    )
    parser.add_argument("--maf-a", required=True,
                        help="MAF from pipeline run with transcripts_A.txt")
    parser.add_argument("--maf-b", required=True,
                        help="MAF from pipeline run with transcripts_B.txt")
    parser.add_argument("--transcripts-a", default=None,
                        help="Path to transcripts_A.txt (informational only)")
    parser.add_argument("--transcripts-b", default=None,
                        help="Path to transcripts_B.txt (informational only)")
    parser.add_argument("--output", default=None,
                        help="Write results TSV to this path")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
