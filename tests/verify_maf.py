#!/usr/bin/env python3
"""
SMART pipeline — MAF verification script.

Usage:
    python tests/verify_maf.py <tier1_maf>

Exit code:
    0  all checks passed
    1  one or more checks failed

Add new expected variants to the EXPECTED list at the bottom of this file
when new verification sets are created.
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

try:
    from pymaftools.core.MAF import MAF as PyMAF
    PYMAFTOOLS_AVAILABLE = True
except ImportError:
    PYMAFTOOLS_AVAILABLE = False


# ---------------------------------------------------------------------------
# Expected variants per verification set
# Key: ID column value from the VCF (used to match rows)
# ---------------------------------------------------------------------------
EXPECTED = [
    # verification1 — 10 known oncogenic hotspot variants
    {
        "ID":                    "BRAF_V600E",
        "Hugo_Symbol":           "BRAF",
        "HGVSp_Short":           "p.V600E",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "KRAS_G12D",
        "Hugo_Symbol":           "KRAS",
        "HGVSp_Short":           "p.G12D",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "KRAS_G13D",
        "Hugo_Symbol":           "KRAS",
        "HGVSp_Short":           "p.G13D",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "PIK3CA_H1047R",
        "Hugo_Symbol":           "PIK3CA",
        "HGVSp_Short":           "p.H1047R",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "TP53_R175H",
        "Hugo_Symbol":           "TP53",
        "HGVSp_Short":           "p.R175H",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "EGFR_L858R",
        "Hugo_Symbol":           "EGFR",
        "HGVSp_Short":           "p.L858R",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "IDH1_R132H",
        "Hugo_Symbol":           "IDH1",
        "HGVSp_Short":           "p.R132H",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "NRAS_Q61R",
        "Hugo_Symbol":           "NRAS",
        "HGVSp_Short":           "p.Q61R",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    {
        "ID":                    "PTEN_R130stop",
        "Hugo_Symbol":           "PTEN",
        "HGVSp_Short":           "p.R130*",
        "Variant_Classification":"Nonsense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "False",
    },
    {
        "ID":                    "ERBB2_S310F",
        "Hugo_Symbol":           "ERBB2",
        "HGVSp_Short":           "p.S310F",
        "Variant_Classification":"Missense_Mutation",
        "Variant_Type":          "SNP",
        "VARIANT_IN_ONCOKB":     "True",
        "ONCOKB_HOTSPOT":        "True",
    },
    # verification1 — 4 Manta-format CNVs (well-known oncogenic copy-number alterations)
    {
        "ID":                    "MantaDUP:ERBB2_AMP",
        "Hugo_Symbol":           "ERBB2",
        "ONCOKB_QUERY_TYPE":     "CNA",
        "ONCOKB_ONCOGENIC":      "Oncogenic",
    },
    {
        "ID":                    "MantaDUP:MET_AMP",
        "Hugo_Symbol":           "MET",
        "ONCOKB_QUERY_TYPE":     "CNA",
        "ONCOKB_ONCOGENIC":      "Oncogenic",
    },
    {
        "ID":                    "MantaDUP:CDK4_AMP",
        "Hugo_Symbol":           "CDK4",
        "ONCOKB_QUERY_TYPE":     "CNA",
        "ONCOKB_ONCOGENIC":      "Oncogenic",
    },
    {
        "ID":                    "MantaDEL:CDKN2A_DEL",
        "Hugo_Symbol":           "CDKN2A",
        "ONCOKB_QUERY_TYPE":     "CNA",
        "ONCOKB_ONCOGENIC":      "Oncogenic",
    },
]


# ---------------------------------------------------------------------------
# Structural checks (apply to every MAF regardless of content)
# ---------------------------------------------------------------------------
REQUIRED_COLUMNS = [
    "Hugo_Symbol",
    "Tumor_Sample_Barcode",
    "HGVSp_Short",
    "HGVSp",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Variant_Type",
    "NCBI_Build",
]

VALID_CLASSIFICATIONS = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
    "Splice_Region", "Silent", "Nonstop_Mutation", "Translation_Start_Site",
    "5'UTR", "3'UTR", "Intron", "5'Flank", "3'Flank", "IGR", "RNA",
    # VEP consequences for structural variants (Manta DUP/DEL)
    "copy_number_gain", "copy_number_loss",
    "transcript_amplification", "transcript_ablation",
    "feature_amplification", "feature_truncation",
}

VALID_TYPES = {"SNP", "DNP", "TNP", "ONP", "INS", "DEL"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_maf(path: Path) -> pd.DataFrame:
    """Read MAF, skipping the #version comment line."""
    with open(path) as fh:
        first = fh.readline()
    skip = 1 if first.startswith("#") else 0
    return pd.read_csv(path, sep="\t", skiprows=skip, dtype=str)


def check(condition: bool, msg: str, passed: list, failed: list) -> None:
    if condition:
        passed.append(f"  PASS  {msg}")
    else:
        failed.append(f"  FAIL  {msg}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Verify SMART tier1 MAF output.")
    parser.add_argument("maf", help="Path to Final_result_tier1.maf")
    args = parser.parse_args()

    maf_path = Path(args.maf)
    if not maf_path.exists():
        print(f"ERROR: file not found: {maf_path}")
        sys.exit(1)

    # --- check #version header ---
    passed, failed = [], []
    with open(maf_path) as fh:
        first_line = fh.readline().rstrip()
    check(first_line.startswith("#version"), f"#version header present ({first_line!r})", passed, failed)

    df = load_maf(maf_path)

    # --- structural checks ---
    for col in REQUIRED_COLUMNS:
        check(col in df.columns, f"Required column present: {col}", passed, failed)

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        print("\n[STRUCTURAL CHECKS]")
        for m in passed + failed:
            print(m)
        print(f"\nAborting: missing required columns {missing_cols}")
        sys.exit(1)

    # --- classification/type validity ---
    bad_class = df[~df["Variant_Classification"].isin(VALID_CLASSIFICATIONS)]["Variant_Classification"].unique()
    check(len(bad_class) == 0,
          f"All Variant_Classification values are valid MAF terms (found invalid: {list(bad_class)})",
          passed, failed)

    bad_type = df[~df["Variant_Type"].isin(VALID_TYPES)]["Variant_Type"].unique()
    check(len(bad_type) == 0,
          f"All Variant_Type values are valid MAF terms (found invalid: {list(bad_type)})",
          passed, failed)

    check(df["Hugo_Symbol"].notna().all(),
          "No empty Hugo_Symbol values", passed, failed)

    check(df["Variant_Classification"].notna().all() and (df["Variant_Classification"] != "").all(),
          "No empty Variant_Classification values", passed, failed)

    check((df["NCBI_Build"] == "GRCh38").all(),
          "NCBI_Build is GRCh38 for all rows", passed, failed)

    # --- per-variant content checks ---
    id_col = "ID" if "ID" in df.columns else None
    for exp in EXPECTED:
        variant_id = exp["ID"]
        if id_col is None:
            failed.append(f"  FAIL  [{variant_id}] ID column not found in MAF")
            continue

        rows = df[df[id_col] == variant_id]
        if rows.empty:
            failed.append(f"  FAIL  [{variant_id}] variant not found in MAF")
            continue

        row = rows.iloc[0]
        for field, expected_val in exp.items():
            if field == "ID":
                continue
            actual = str(row.get(field, "")) if field in df.columns else "COLUMN_MISSING"
            check(actual == expected_val,
                  f"[{variant_id}] {field}: expected={expected_val!r}  actual={actual!r}",
                  passed, failed)

    # --- pymaftools compatibility check ---
    if PYMAFTOOLS_AVAILABLE:
        try:
            maf_obj = PyMAF(df)

            # 1. Object is created and has the right number of rows
            check(len(maf_obj) == len(df),
                  f"pymaftools: MAF object created with correct row count ({len(maf_obj)})",
                  passed, failed)

            # 2. Variant_Classification breakdown matches expectations
            vc_counts = maf_obj["Variant_Classification"].value_counts().to_dict()
            expected_missense = sum(1 for e in EXPECTED if e["Variant_Classification"] == "Missense_Mutation")
            expected_nonsense = sum(1 for e in EXPECTED if e["Variant_Classification"] == "Nonsense_Mutation")
            check(vc_counts.get("Missense_Mutation", 0) == expected_missense,
                  f"pymaftools: Missense_Mutation count = {vc_counts.get('Missense_Mutation', 0)} (expected {expected_missense})",
                  passed, failed)
            check(vc_counts.get("Nonsense_Mutation", 0) == expected_nonsense,
                  f"pymaftools: Nonsense_Mutation count = {vc_counts.get('Nonsense_Mutation', 0)} (expected {expected_nonsense})",
                  passed, failed)

            # 3. All Variant_Type values recognised by pymaftools
            vt_unique = set(maf_obj["Variant_Type"].dropna().unique())
            check(vt_unique <= VALID_TYPES,
                  f"pymaftools: all Variant_Type values valid ({vt_unique})",
                  passed, failed)

            # 4. filter_maf works (basic API smoke test)
            missense_only = maf_obj.filter_maf(["Missense_Mutation"])
            check(len(missense_only) == expected_missense,
                  f"pymaftools: filter_maf(Missense_Mutation) returns {len(missense_only)} rows",
                  passed, failed)

        except Exception as exc:
            failed.append(f"  FAIL  pymaftools: unexpected error — {exc}")
    else:
        print("  SKIP  pymaftools not installed — skipping compatibility checks")

    # --- summary ---
    print(f"\n{'='*60}")
    print(f"MAF verification: {maf_path.name}")
    print(f"{'='*60}")
    print(f"  Rows:    {len(df)}")
    print(f"  Columns: {len(df.columns)}")
    print()
    print("[ALL CHECKS]")
    for m in passed + failed:
        print(m)

    print(f"\n  Passed: {len(passed)}  /  Failed: {len(failed)}  /  Total: {len(passed)+len(failed)}")
    print("="*60)

    if failed:
        print("RESULT: FAILED")
        sys.exit(1)
    else:
        print("RESULT: ALL CHECKS PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
