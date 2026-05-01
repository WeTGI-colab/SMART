#!/usr/bin/env python3
"""
SMART verification4 — parallel processing plumbing check.

Checks that running with --nodes 2 produces correct output for every sample:
  1. Per-sample log file exists in logs/ and contains the completion marker.
  2. variant_counts.txt has exactly one data row per sample.
  3. Output MAF exists and is non-empty for every sample.
  4. All MAFs are content-identical when the Tumor_Sample_Barcode column is
     excluded (same input VCF → same annotations, regardless of sample name).

Usage:
    python3 tests/verification4/verify.py \\
        --output-dir tests/verification4/output \\
        --samples par_sample_A par_sample_B \\
        --results  tests/verification4/results.tsv

Exit code 0 if all checks pass, 1 if any check fails.
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--results", required=True, type=Path)
    return p.parse_args()


def check(results: list, name: str, status: str, detail: str = ""):
    results.append({"check": name, "status": status, "detail": detail})
    marker = "PASS" if status == "PASS" else "FAIL"
    suffix = f" — {detail}" if detail else ""
    print(f"  [{marker}] {name}{suffix}")


def maf_rows_without_sample_col(maf_path: Path) -> list[str]:
    """Return data rows with Tumor_Sample_Barcode column blanked out."""
    lines = maf_path.read_text().splitlines()
    # MAF files may have comment lines starting with #
    header_line = next(l for l in lines if not l.startswith("#"))
    cols = header_line.split("\t")
    try:
        tsb_idx = cols.index("Tumor_Sample_Barcode")
    except ValueError:
        tsb_idx = None

    normalised = []
    for line in lines:
        if line.startswith("#"):
            continue
        if line == header_line:
            continue
        if tsb_idx is not None:
            parts = line.split("\t")
            if len(parts) > tsb_idx:
                parts[tsb_idx] = "<SAMPLE>"
            normalised.append("\t".join(parts))
        else:
            normalised.append(line)
    return normalised


def main():
    args = parse_args()
    results = []
    failed = False

    logs_dir = args.output_dir / "logs"
    counts_file = args.output_dir / "variant_counts.txt"
    final_table_dir = args.output_dir / "FINAL_Table"

    print(f"\nOutput directory: {args.output_dir}")
    print(f"Samples:          {', '.join(args.samples)}")
    print()

    # ------------------------------------------------------------------
    # 1. Per-sample log files exist and contain completion marker
    # ------------------------------------------------------------------
    print("--- Log file checks ---")
    for sample in args.samples:
        log_path = logs_dir / f"{sample}.log"
        if not log_path.exists():
            check(results, f"log_exists:{sample}", "FAIL",
                  f"not found: {log_path}")
            failed = True
            continue
        check(results, f"log_exists:{sample}", "PASS")

        content = log_path.read_text()
        marker = f"--- {sample} complete ---"
        if marker in content:
            check(results, f"log_complete:{sample}", "PASS")
        else:
            check(results, f"log_complete:{sample}", "FAIL",
                  f"completion marker '{marker}' not found in log")
            failed = True

    # ------------------------------------------------------------------
    # 2. variant_counts.txt has one row per sample
    # ------------------------------------------------------------------
    print("\n--- Variant counts check ---")
    if not counts_file.exists():
        check(results, "counts_file_exists", "FAIL", str(counts_file))
        failed = True
    else:
        check(results, "counts_file_exists", "PASS")
        lines = [l for l in counts_file.read_text().splitlines() if l.strip()]
        data_rows = lines[1:]  # skip header
        expected = len(args.samples)
        if len(data_rows) == expected:
            check(results, "counts_row_count", "PASS",
                  f"{len(data_rows)} rows as expected")
        else:
            check(results, "counts_row_count", "FAIL",
                  f"expected {expected} data rows, got {len(data_rows)}")
            failed = True

        # Check every sample appears in counts
        for sample in args.samples:
            found = any(row.startswith(sample + "\t") for row in data_rows)
            if found:
                check(results, f"counts_has_sample:{sample}", "PASS")
            else:
                check(results, f"counts_has_sample:{sample}", "FAIL",
                      "sample name not found in variant_counts.txt")
                failed = True

    # ------------------------------------------------------------------
    # 3. Output MAF exists and is non-empty for every sample
    # ------------------------------------------------------------------
    print("\n--- MAF existence checks ---")
    maf_paths = {}
    for sample in args.samples:
        maf_path = final_table_dir / f"{sample}.oncoKB.maf"
        if not maf_path.exists():
            check(results, f"maf_exists:{sample}", "FAIL", str(maf_path))
            failed = True
        elif maf_path.stat().st_size == 0:
            check(results, f"maf_exists:{sample}", "FAIL", "file is empty")
            failed = True
        else:
            check(results, f"maf_exists:{sample}", "PASS")
            maf_paths[sample] = maf_path

    # ------------------------------------------------------------------
    # 4. All MAFs are content-identical (excluding Tumor_Sample_Barcode)
    # ------------------------------------------------------------------
    print("\n--- MAF content identity check ---")
    if len(maf_paths) < 2:
        check(results, "maf_content_identical", "SKIP",
              "fewer than 2 MAFs available — skipping comparison")
    else:
        samples_with_mafs = list(maf_paths.keys())
        reference_sample = samples_with_mafs[0]
        ref_rows = maf_rows_without_sample_col(maf_paths[reference_sample])

        all_identical = True
        for sample in samples_with_mafs[1:]:
            other_rows = maf_rows_without_sample_col(maf_paths[sample])
            if ref_rows == other_rows:
                check(results, f"maf_identical:{reference_sample}_vs_{sample}",
                      "PASS", f"{len(ref_rows)} data rows match")
            else:
                diff_count = sum(
                    1 for a, b in zip(ref_rows, other_rows) if a != b
                )
                row_diff = abs(len(ref_rows) - len(other_rows))
                detail = (
                    f"{diff_count} differing rows, "
                    f"{row_diff} row count difference"
                )
                check(results, f"maf_identical:{reference_sample}_vs_{sample}",
                      "FAIL", detail)
                failed = True
                all_identical = False

        if all_identical and len(samples_with_mafs) >= 2:
            check(results, "maf_content_identical", "PASS",
                  "all sample MAFs match after excluding Tumor_Sample_Barcode")

    # ------------------------------------------------------------------
    # Write results TSV
    # ------------------------------------------------------------------
    args.results.parent.mkdir(parents=True, exist_ok=True)
    with args.results.open("w") as f:
        f.write("Check\tStatus\tDetail\n")
        for r in results:
            f.write(f"{r['check']}\t{r['status']}\t{r['detail']}\n")

    print(f"\nResults written to: {args.results}")

    total = len(results)
    passed = sum(1 for r in results if r["status"] == "PASS")
    print(f"\nSummary: {passed}/{total} checks passed")

    if failed:
        print("VERIFICATION 4: FAIL")
        sys.exit(1)
    else:
        print("VERIFICATION 4: PASS")
        sys.exit(0)


if __name__ == "__main__":
    main()
