#!/usr/bin/env python3
"""
run_azurify.py — Run the Azurify ML classifier on SMART Tier 3 output.

Azurify (faryabiLab/Azurify, bioRxiv 2025) uses a CatBoost gradient-boosted
tree model trained on >15,000 clinically classified variants to predict
variant pathogenicity (Pathogenic / Likely Pathogenic / VUS /
Likely Benign / Benign).

This script:
  1. Reads SMART's Final_result_tier3.tsv
  2. Converts it to Azurify's tab-delimited input format
  3. Runs Azurify (with snpEff annotation + internal feature lookup)
  4. Merges the ML predictions back into the SMART output
  5. Writes Final_result_tier3_ml.tsv

Usage:
    python3 utils/run_azurify.py \\
        --smart-output  /path/to/Output_Results \\
        --azurify-dir   /path/to/Azurify \\
        --snpeff-jar    /path/to/snpEff/snpEff.jar \\
        --output-dir    /path/to/Output_Results/azurify

Requirements (install outside Docker):
    pip install catboost liftover tqdm pandas
    git clone https://github.com/faryabiLab/Azurify /path/to/Azurify

Notes:
  - Requires >16 GB RAM.
  - Azurify works on SNVs and small indels only; CNAs and structural
    variants are skipped and preserved in the merged output with
    empty ML columns.
  - VAF must be a decimal in SMART (e.g. 0.33); it is converted to
    percentage (33.0) for Azurify internally.
  - Coordinates are converted from 1-based (VCF/MAF) to 0-based
    half-open (BED), as expected by Azurify.
"""
import argparse
import os
import subprocess
import sys
import glob

import pandas as pd


# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Run Azurify ML classifier on SMART Tier 3 output."
    )
    p.add_argument(
        "--smart-output", required=True,
        help="Directory containing Final_result_tier3.tsv "
             "(the Output_Results/ folder from a SMART run)."
    )
    p.add_argument(
        "--azurify-dir", required=True,
        help="Path to a cloned Azurify repository "
             "(git clone https://github.com/faryabiLab/Azurify)."
    )
    p.add_argument(
        "--snpeff-jar", required=True,
        help="Full path to snpEff.jar (required by Azurify for "
             "internal variant annotation)."
    )
    p.add_argument(
        "--output-dir", default=None,
        help="Directory for Azurify intermediate files and the merged "
             "output. Defaults to <smart-output>/azurify_ml/."
    )
    p.add_argument(
        "--no-drug-targets", action="store_true",
        help="Use the Azurify model variant that omits drug-target "
             "features (--no_drug_targets flag in Azurify)."
    )
    return p.parse_args()


# ── SMART → Azurify input conversion ─────────────────────────────────────────

def load_smart_tier3(smart_output_dir):
    """Find and load Final_result_tier3.tsv from the SMART output directory."""
    path = os.path.join(smart_output_dir, "Final_result_tier3.tsv")
    if not os.path.isfile(path):
        sys.exit(f"ERROR: Final_result_tier3.tsv not found in {smart_output_dir}")

    # Row 0 = column names, Row 1 = source metadata → skip metadata row
    df = pd.read_csv(path, sep="\t", skiprows=[1], dtype=str)
    print(f"Loaded {len(df)} rows from {path}")
    return df, path


def build_azurify_input(smart_df):
    """
    Convert SMART Tier 3 columns to the Azurify input format.

    Azurify expects:
      CHROM  START  STOP  REF  ALT  VAF   (tab-delimited)

    Coordinate system:
      SMART uses 1-based positions (VCF/MAF standard).
      Azurify uses 0-based half-open intervals (BED standard).
        START = Start_Position - 1
        STOP  = Start_Position  (for SNVs and most indels)

    VAF:
      SMART stores VAF as a decimal fraction (e.g. 0.326).
      Azurify expects a percentage (e.g. 32.6).
    """
    df = smart_df.copy()

    # Keep a merge key (original CHROM + 1-based position + alleles)
    df["_MERGE_KEY"] = (
        df["Chromosome"].astype(str) + ":" +
        df["Start_Position"].astype(str) + ":" +
        df["REF"].astype(str) + ":" +
        df["ALT"].astype(str)
    )

    # Filter to rows suitable for Azurify (SNVs and indels — skip
    # symbolic alleles like <DEL>, <DUP>, <INS> used for CNAs)
    is_snv_indel = (
        df["REF"].notna() & df["ALT"].notna() &
        ~df["ALT"].str.startswith("<") &
        ~df["REF"].str.startswith("<")
    )
    az_df = df[is_snv_indel].copy()
    skipped = len(df) - len(az_df)
    if skipped:
        print(f"  Skipped {skipped} CNA/symbolic rows (not supported by Azurify).")

    # Coordinate conversion: 1-based → 0-based BED
    az_df["START"] = pd.to_numeric(az_df["Start_Position"], errors="coerce") - 1
    az_df["STOP"]  = pd.to_numeric(az_df["Start_Position"], errors="coerce")

    # VAF: decimal → percentage
    az_df["VAF_PCT"] = pd.to_numeric(az_df["VAF"], errors="coerce") * 100

    # Build the 6-column Azurify input table.
    # Pass SAMPLE_ID as an extra column so we can track rows after the run.
    az_input = pd.DataFrame({
        "CHROM": az_df["Chromosome"],
        "START": az_df["START"].astype("Int64"),
        "STOP":  az_df["STOP"].astype("Int64"),
        "REF":   az_df["REF"],
        "ALT":   az_df["ALT"],
        "VAF":   az_df["VAF_PCT"].round(4),
        "SAMPLE": az_df.get("Tumor_Sample_Barcode", "UNKNOWN"),
        "_MERGE_KEY": az_df["_MERGE_KEY"],
    })

    az_input = az_input.dropna(subset=["CHROM", "START", "REF", "ALT"])
    print(f"  Built Azurify input: {len(az_input)} variants.")
    return az_input, az_df["_MERGE_KEY"].tolist()


# ── Run Azurify ───────────────────────────────────────────────────────────────

def run_azurify(azurify_dir, input_file, output_dir, snpeff_jar, no_drug_targets):
    """Call azurify.py as a subprocess."""
    cmd = [
        sys.executable,
        os.path.join(azurify_dir, "azurify.py"),
        "-i", input_file,
        "-o", output_dir,
        "-s", snpeff_jar,
        "-g", "hg38",
    ]
    if no_drug_targets:
        cmd += ["-d", "1"]

    print(f"\nRunning Azurify:")
    print("  " + " ".join(cmd))
    print("  (this may take several minutes — Azurify requires >16 GB RAM)\n")

    result = subprocess.run(cmd, cwd=azurify_dir)
    if result.returncode != 0:
        sys.exit(f"ERROR: Azurify exited with code {result.returncode}.")

    # Find the output TSV produced by Azurify
    matches = glob.glob(os.path.join(output_dir, "*.azurify.tsv"))
    if not matches:
        sys.exit(f"ERROR: No .azurify.tsv file found in {output_dir} after run.")
    return matches[0]


# ── Merge predictions back into SMART output ──────────────────────────────────

ML_COLS = ["Pathogenicity", "BP", "PP", "LBP", "LPP", "VP"]

ML_COL_DESCRIPTIONS = {
    "Pathogenicity": "Azurify ML classification: Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign",
    "BP":            "Azurify probability — Benign",
    "PP":            "Azurify probability — Pathogenic",
    "LBP":           "Azurify probability — Likely Benign",
    "LPP":           "Azurify probability — Likely Pathogenic",
    "VP":            "Azurify probability — VUS (Variant of Uncertain Significance)",
}


def merge_ml_into_smart(smart_df, azurify_output_file):
    """
    Merge Azurify ML predictions into the SMART Tier 3 DataFrame.

    Matching is done on CHROM:START(1-based):REF:ALT.
    Rows that Azurify could not annotate (CNAs, liftover failures)
    receive empty ML columns.
    """
    az = pd.read_csv(azurify_output_file, sep="\t", dtype=str)

    # Reconstruct the merge key from Azurify output coordinates.
    # Azurify START is 0-based → convert back to 1-based for matching.
    az["START_1based"] = (pd.to_numeric(az["START"], errors="coerce") + 1).astype("Int64").astype(str)
    az["_MERGE_KEY"] = (
        az["CHROM"].astype(str) + ":" +
        az["START_1based"] + ":" +
        az["REF"].astype(str) + ":" +
        az["ALT"].astype(str)
    )

    # Keep only the ML result columns + merge key
    az_ml = az[["_MERGE_KEY"] + [c for c in ML_COLS if c in az.columns]].copy()
    az_ml = az_ml.drop_duplicates(subset="_MERGE_KEY")

    # Build merge key on the SMART side
    smart_df["_MERGE_KEY"] = (
        smart_df["Chromosome"].astype(str) + ":" +
        smart_df["Start_Position"].astype(str) + ":" +
        smart_df["REF"].astype(str) + ":" +
        smart_df["ALT"].astype(str)
    )

    merged = smart_df.merge(az_ml, on="_MERGE_KEY", how="left")
    merged = merged.drop(columns=["_MERGE_KEY"])

    matched = merged["Pathogenicity"].notna().sum() if "Pathogenicity" in merged.columns else 0
    print(f"\nMerge complete: {matched}/{len(merged)} rows received ML predictions.")
    return merged


# ── Write merged output ───────────────────────────────────────────────────────

def write_output(merged_df, smart_output_dir, original_tier3_path):
    """
    Write Final_result_tier3_ml.tsv with a two-row header
    (column names + source metadata) matching SMART's convention.
    """
    out_path = os.path.join(smart_output_dir, "Final_result_tier3_ml.tsv")

    # Build the source-metadata row: reuse existing descriptions for
    # SMART columns, add new descriptions for ML columns.
    # Read the original metadata row from tier3.
    original = pd.read_csv(original_tier3_path, sep="\t", header=None, nrows=2)
    meta_map = dict(zip(original.iloc[0], original.iloc[1]))

    for col, desc in ML_COL_DESCRIPTIONS.items():
        if col in merged_df.columns:
            meta_map[col] = f"{desc} | Azurify | faryabiLab/Azurify"

    meta_row = [meta_map.get(c, "") for c in merged_df.columns]

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\t".join(merged_df.columns) + "\n")
        f.write("\t".join(str(v) for v in meta_row) + "\n")
        merged_df.to_csv(f, sep="\t", index=False, header=False)

    print(f"\nMerged output written to: {out_path}")
    print(f"Columns added: {[c for c in ML_COLS if c in merged_df.columns]}")

    # Print a quick summary of ML classifications
    if "Pathogenicity" in merged_df.columns:
        print("\nAzurify classification summary:")
        counts = merged_df["Pathogenicity"].value_counts(dropna=False)
        for label, n in counts.items():
            print(f"  {label}: {n}")

    return out_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # Validate paths
    # Resolve all paths to absolute to avoid cwd issues when Azurify runs
    args.smart_output = os.path.abspath(args.smart_output)
    args.azurify_dir  = os.path.abspath(args.azurify_dir)
    args.snpeff_jar   = os.path.abspath(args.snpeff_jar)

    if not os.path.isdir(args.smart_output):
        sys.exit(f"ERROR: --smart-output directory not found: {args.smart_output}")
    if not os.path.isdir(args.azurify_dir):
        sys.exit(f"ERROR: --azurify-dir not found: {args.azurify_dir}")
    if not os.path.isfile(args.snpeff_jar):
        sys.exit(f"ERROR: --snpeff-jar not found: {args.snpeff_jar}")

    output_dir = os.path.abspath(
        args.output_dir or os.path.join(args.smart_output, "azurify_ml")
    )
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("SMART × Azurify — ML variant classification")
    print("=" * 60)

    # 1. Load SMART Tier 3
    print("\n[1/4] Loading SMART Tier 3 output...")
    smart_df, tier3_path = load_smart_tier3(args.smart_output)

    # 2. Convert to Azurify input
    print("\n[2/4] Converting to Azurify input format...")
    az_input, _ = build_azurify_input(smart_df)
    az_input_file = os.path.abspath(os.path.join(output_dir, "smart_azurify_input.txt"))
    az_input.to_csv(az_input_file, sep="\t", index=False)
    print(f"  Azurify input written to: {az_input_file}")

    # 3. Run Azurify
    print("\n[3/4] Running Azurify classifier...")
    az_output_file = run_azurify(
        azurify_dir=args.azurify_dir,
        input_file=az_input_file,
        output_dir=output_dir,
        snpeff_jar=args.snpeff_jar,
        no_drug_targets=args.no_drug_targets,
    )
    print(f"  Azurify output: {az_output_file}")

    # 4. Merge and write
    print("\n[4/4] Merging ML predictions into SMART output...")
    merged = merge_ml_into_smart(smart_df, az_output_file)
    write_output(merged, args.smart_output, tier3_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
