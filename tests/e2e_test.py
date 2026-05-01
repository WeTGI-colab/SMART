#!/usr/bin/env python3
"""
End-to-end regression test for SMART.

Runs the SMART Docker pipeline on a small test VCF, then compares
selected columns in the output Final_result_tier3.tsv against a
pre-validated reference table.

Usage:
    python tests/e2e_test.py

Prerequisites:
    - Docker Desktop running
    - SMART image built (docker build -t smart:latest .)
    - tests/data/test_input.vcf populated with test variants
    - tests/data/validated_output.tsv populated with expected output rows
    - Fill in the CONFIG section below before running
"""

import os
import sys
import shutil
import subprocess
import tempfile
import pandas as pd

# =============================================================================
# CONFIG — fill these in before running
# =============================================================================

# Path to the SMART Docker image
SMART_IMAGE = "smart:latest"

# OncoKB API token
ONCOKB_TOKEN = ""  # e.g. "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"

# Path to the reference resources directory mounted as /refs inside Docker
REF_DIR = ""  # e.g. "/path/to/refs"

# Path to the TSO500 transcript whitelist
TRANSCRIPTS_FILE = ""  # e.g. "/path/to/TSO500_transcripts_list.txt"

# Path to Config.yaml
CONFIG_FILE = os.path.join(os.path.dirname(__file__), "..", "Config.yaml")

# Columns to compare between the pipeline output and the validated table.
# The tier3 output has a two-row header; row 1 = column names, row 2 = metadata.
# List the exact column names from row 1 of Final_result_tier3.tsv.
COLUMNS_TO_CHECK = [
    # Fill in the columns you want to validate, e.g.:
    # "SYMBOL",
    # "HGVSp_Short",
    # "Consequence",
    # "ONCOGENIC",
    # "HIGHEST_LEVEL",
    # "ClinVar_CLNSIG",
    # "VAF",
]

# Column(s) used to match rows between the output and the validated table.
# Must uniquely identify each variant (used as join key).
ROW_KEY_COLUMNS = ["Chromosome", "Start_Position", "REF", "ALT"]

# Paths to test data files
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_VCF = os.path.join(SCRIPT_DIR, "data", "test_input.vcf")
VALIDATED_TABLE = os.path.join(SCRIPT_DIR, "data", "validated_output.tsv")

# =============================================================================
# HELPERS
# =============================================================================

def check_prerequisites():
    errors = []
    if not ONCOKB_TOKEN:
        errors.append("ONCOKB_TOKEN is not set")
    if not REF_DIR or not os.path.isdir(REF_DIR):
        errors.append(f"REF_DIR not found: {REF_DIR!r}")
    if not TRANSCRIPTS_FILE or not os.path.isfile(TRANSCRIPTS_FILE):
        errors.append(f"TRANSCRIPTS_FILE not found: {TRANSCRIPTS_FILE!r}")
    if not os.path.isfile(CONFIG_FILE):
        errors.append(f"CONFIG_FILE not found: {CONFIG_FILE!r}")
    if not os.path.isfile(TEST_VCF):
        errors.append(f"TEST_VCF not found: {TEST_VCF!r}")
    if not os.path.isfile(VALIDATED_TABLE):
        errors.append(f"VALIDATED_TABLE not found: {VALIDATED_TABLE!r}")
    if not COLUMNS_TO_CHECK:
        errors.append("COLUMNS_TO_CHECK is empty — add columns to validate")
    result = subprocess.run(
        ["docker", "info"], capture_output=True
    )
    if result.returncode != 0:
        errors.append("Docker daemon is not running")
    result = subprocess.run(
        ["docker", "image", "inspect", SMART_IMAGE], capture_output=True
    )
    if result.returncode != 0:
        errors.append(f"Docker image not found: {SMART_IMAGE!r}")
    if errors:
        print("PREREQUISITES FAILED:")
        for e in errors:
            print(f"  ✗ {e}")
        sys.exit(1)
    print("Prerequisites OK")


def run_smart(workdir: str) -> str:
    """Run SMART inside Docker and return the path to Final_result_tier3.tsv."""
    vcf_dir = os.path.join(workdir, "OriginalVcf")
    os.makedirs(vcf_dir)

    sample_name = "e2e_test_sample"
    dest_vcf = os.path.join(vcf_dir, f"{sample_name}.dragen.concat_snv_sv.vcf.gz")

    # bgzip-compress the test VCF (using Docker since bgzip may not be local)
    print(f"Compressing test VCF...")
    result = subprocess.run(
        [
            "docker", "run", "--rm",
            "-v", f"{os.path.dirname(TEST_VCF)}:/input",
            "-v", f"{vcf_dir}:/out",
            "quay.io/biocontainers/htslib:1.21--h566b1c6_1",
            "sh", "-c",
            f"bgzip -c /input/{os.path.basename(TEST_VCF)} > /out/{os.path.basename(dest_vcf)}"
            f" && tabix -p vcf /out/{os.path.basename(dest_vcf)}"
        ],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("Failed to compress test VCF:")
        print(result.stderr)
        sys.exit(1)

    # Copy config
    shutil.copy(CONFIG_FILE, os.path.join(workdir, "Config.yaml"))

    print("Running SMART pipeline...")
    result = subprocess.run(
        [
            "docker", "run", "--rm",
            "-v", f"{workdir}:/data",
            "-v", f"{REF_DIR}:/refs:ro",
            SMART_IMAGE,
            ONCOKB_TOKEN,
            "--transcripts-file", f"/data/{os.path.basename(TRANSCRIPTS_FILE)}",
            "--ref-dir", "/refs",
            "--config", "/data/Config.yaml",
            "--no-pass",       # test VCF is already PASS-only
            "--keep-tmp",
        ],
        capture_output=True, text=True
    )
    print(result.stdout)
    if result.returncode != 0:
        print("SMART pipeline FAILED:")
        print(result.stderr)
        sys.exit(1)

    output_table = os.path.join(workdir, "Final_result_tier3.tsv")
    if not os.path.isfile(output_table):
        print(f"Output table not found: {output_table}")
        sys.exit(1)
    return output_table


def load_tier3(path: str) -> pd.DataFrame:
    """Load a tier3 TSV, skipping the metadata row (row 2)."""
    df = pd.read_csv(path, sep="\t", header=0, skiprows=[1], dtype=str)
    df.columns = df.columns.str.strip()
    return df


def compare_tables(output_path: str, validated_path: str) -> bool:
    output_df = load_tier3(output_path)
    validated_df = load_tier3(validated_path)

    missing_cols = [c for c in COLUMNS_TO_CHECK + ROW_KEY_COLUMNS if c not in output_df.columns]
    if missing_cols:
        print(f"Columns missing from pipeline output: {missing_cols}")
        sys.exit(1)
    missing_cols = [c for c in COLUMNS_TO_CHECK + ROW_KEY_COLUMNS if c not in validated_df.columns]
    if missing_cols:
        print(f"Columns missing from validated table: {missing_cols}")
        sys.exit(1)

    output_df  = output_df.set_index(ROW_KEY_COLUMNS)
    validated_df = validated_df.set_index(ROW_KEY_COLUMNS)

    all_passed = True
    results = []

    for idx in validated_df.index:
        if idx not in output_df.index:
            results.append({
                "variant": str(idx),
                "column": "—",
                "status": "FAIL",
                "expected": "row present",
                "got": "row MISSING from output",
            })
            all_passed = False
            continue

        for col in COLUMNS_TO_CHECK:
            expected = str(validated_df.loc[idx, col]).strip()
            got      = str(output_df.loc[idx, col]).strip()
            passed   = (expected == got)
            results.append({
                "variant": str(idx),
                "column":  col,
                "status":  "PASS" if passed else "FAIL",
                "expected": expected,
                "got":      got,
            })
            if not passed:
                all_passed = False

    # Print report
    col_w = max(len(c) for c in COLUMNS_TO_CHECK + ["column"]) + 2
    print()
    print(f"{'VARIANT':<50} {'COLUMN':<{col_w}} {'STATUS':<6}  EXPECTED  →  GOT")
    print("─" * 120)
    for r in results:
        line = f"{r['variant']:<50} {r['column']:<{col_w}} {r['status']:<6}"
        if r["status"] == "FAIL":
            line += f"  {r['expected']!r}  →  {r['got']!r}"
        print(line)
    print("─" * 120)
    passed_n = sum(1 for r in results if r["status"] == "PASS")
    print(f"\nResult: {passed_n}/{len(results)} checks passed")
    return all_passed


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("SMART end-to-end test")
    print("=" * 60)

    check_prerequisites()

    with tempfile.TemporaryDirectory(prefix="smart_e2e_") as workdir:
        print(f"Working directory: {workdir}")

        # Copy transcripts file into workdir so Docker can mount it
        shutil.copy(TRANSCRIPTS_FILE, workdir)

        output_table = run_smart(workdir)
        passed = compare_tables(output_table, VALIDATED_TABLE)

    print()
    if passed:
        print("ALL CHECKS PASSED")
        sys.exit(0)
    else:
        print("SOME CHECKS FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
