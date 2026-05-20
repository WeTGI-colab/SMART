#!/usr/bin/env python3
"""
build_genie_lookup.py — Pre-process a somatic variant MAF into a compact
lookup table used by SMART's post_analysis.py to apply SVIG-UK O4 scoring.

Run this ONCE after downloading the somatic variant dataset. The resulting
genie_lookup.tsv.gz file is small (~5–15 MB) and is placed in the reference
directory for use by every subsequent SMART pipeline run.

Usage:
    python3 utils/build_genie_lookup.py \\
        --input  /path/to/data_mutations.txt \\
        --output /path/to/refs/GENIE/genie_lookup.tsv.gz

Compatible input formats:
    - AACR Project GENIE   (download from Synapse: syn7222066, file: data_mutations.txt)
    - MSK-IMPACT 2017      (Zehir et al., Nat Med 2017; cBioPortal study msk_impact_2017)
    - Any standard MAF with Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode columns

SVIG-UK O4 thresholds applied during post_analysis.py:
    Missense / splice variants:
        count > 10   →  O4 Strong  [+4]
        count 5–10   →  O4 Moderate [+2]
        count 1–4    →  O4 Supporting [+1]

    Frameshift / nonsense variants:
        count > 50   →  O4 Strong  [+4]
        count 20–50  →  O4 Moderate [+2]
        count 10–19  →  O4 Supporting [+1]

    In-frame indels:
        count > 50   →  O4 Strong  [+4]
        count 20–50  →  O4 Moderate [+2]
        count 10–19  →  O4 Supporting [+1]
"""

import argparse
import gzip
import os
import sys

import pandas as pd


# MAF columns we need — both standard names and common alternatives
GENE_COLS    = ["Hugo_Symbol", "Gene", "SYMBOL"]
PROTEIN_COLS = ["HGVSp_Short", "Protein_Change", "AAChange"]
SAMPLE_COLS  = ["Tumor_Sample_Barcode", "sample_id", "Sample_ID"]
PATIENT_COLS = ["Patient_ID", "patient_id", "PATIENT_ID"]


def _first_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def load_maf(path):
    print(f"Loading MAF: {path}")
    # Skip comment lines (start with #), handle tab-delimited
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        skip = 0
        for line in f:
            if line.startswith("#"):
                skip += 1
            else:
                break
    df = pd.read_csv(
        path,
        sep="\t",
        skiprows=skip,
        dtype=str,
        low_memory=False,
    )
    print(f"  {len(df):,} rows, {len(df.columns)} columns")
    return df


def build_lookup(df):
    """
    Build a lookup: (gene, protein_change) → deduplicated patient count.

    Deduplication is critical — GENIE includes multiple samples per patient
    (primary + metastasis, multiple timepoints). Counting samples instead of
    patients inflates O4 scores. We count unique patients per variant.
    If no Patient_ID column exists, we fall back to Tumor_Sample_Barcode
    (still better than no deduplication).
    """
    gene_col    = _first_col(df, GENE_COLS)
    protein_col = _first_col(df, PROTEIN_COLS)
    sample_col  = _first_col(df, SAMPLE_COLS)
    patient_col = _first_col(df, PATIENT_COLS)

    if not gene_col:
        sys.exit("ERROR: No gene symbol column found. Expected one of: " + str(GENE_COLS))
    if not protein_col:
        sys.exit("ERROR: No protein change column found. Expected one of: " + str(PROTEIN_COLS))
    if not sample_col:
        sys.exit("ERROR: No sample barcode column found. Expected one of: " + str(SAMPLE_COLS))

    id_col = patient_col or sample_col
    id_label = "patient" if patient_col else "sample (no Patient_ID column found)"
    print(f"\n  Deduplicating by {id_label}: '{id_col}'")

    # Keep only rows with a gene and a protein change
    work = df[[gene_col, protein_col, id_col]].copy()
    work.columns = ["gene", "protein_change", "id"]
    work = work.dropna(subset=["gene", "protein_change"])
    work = work[work["protein_change"].str.startswith("p.", na=False)]
    print(f"  {len(work):,} rows with valid gene + protein change")

    # Normalise protein change: strip ENSP/NP prefix if present
    # e.g. "ENSP00000123:p.V600E" → "p.V600E"
    work["protein_change"] = work["protein_change"].str.replace(
        r"^[A-Z0-9_]+:\.", "p.", regex=True
    )

    # Count unique patients per (gene, protein_change)
    grouped = (
        work.groupby(["gene", "protein_change"])["id"]
        .nunique()
        .reset_index()
        .rename(columns={"id": "count"})
    )

    print(f"  {len(grouped):,} unique (gene, protein_change) combinations")
    print(f"  Top 5 by count:")
    for _, row in grouped.nlargest(5, "count").iterrows():
        print(f"    {row['gene']} {row['protein_change']}  →  {row['count']}")

    return grouped


def write_lookup(df, output_path):
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    opener = gzip.open if output_path.endswith(".gz") else open
    with opener(output_path, "wt", encoding="utf-8") as f:
        df.to_csv(f, sep="\t", index=False)
    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"\nLookup written: {output_path}  ({size_mb:.1f} MB)")


def main():
    parser = argparse.ArgumentParser(
        description="Build a (gene, protein_change) → count lookup from a somatic MAF."
    )
    parser.add_argument("--input",  required=True,
                        help="Input MAF file (plain or .gz). Accepts GENIE, MSK-IMPACT, etc.")
    parser.add_argument("--output", required=True,
                        help="Output lookup file path (.tsv or .tsv.gz recommended).")
    args = parser.parse_args()

    df  = load_maf(args.input)
    lut = build_lookup(df)
    write_lookup(lut, args.output)

    print("\nDone. Place the output file at:")
    print("  /path/to/refs/GENIE/genie_lookup.tsv.gz")
    print("and pass it to post_analysis.py via --genie-counts.")


if __name__ == "__main__":
    main()
