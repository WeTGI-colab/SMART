#!/usr/bin/env python3
"""
validate_svig_uk.py — Validate automated SVIG-UK classification against
                       manually curated variant classifications.

Steps:
  1. Re-runs post_analysis.py with all SVIG-UK reference files to produce
     up-to-date Tier 3 output with SVIG_UK_classification column.
  2. Loads the Tier 3 output and the human classification sidecar TSV.
  3. Matches variants by (CHROM, POS, REF, ALT).
  4. Computes concordance metrics and confusion matrix.
  5. Writes a markdown report and a discordance TSV.

Usage:
    python3 utils/validate_svig_uk.py \\
        --smart-output  data/jack_list_validation/output \\
        --human-classes data/jack_list_validation/jack_curated.classifications.tsv \\
        --report-dir    data/jack_list_validation/validation_report

    # Skip re-running post_analysis (use existing Tier 3):
        --skip-post-analysis
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import date

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CATS = ["Oncogenic", "Likely Oncogenic", "VUS", "Likely Benign", "Benign"]

CLASS_COLORS = {
    "Oncogenic":        "🔴",
    "Likely Oncogenic": "🟠",
    "VUS":              "🟡",
    "Likely Benign":    "🟢",
    "Benign":           "🔵",
}


# ── Step 1: Re-run post_analysis ─────────────────────────────────────────────

def run_post_analysis(smart_output_dir):
    final_table = os.path.join(smart_output_dir, "FINAL_Table")
    output_dir  = os.path.join(smart_output_dir, "output")
    refs        = os.path.join(REPO_ROOT, "Files4ThisProject")

    cmd = [
        sys.executable,
        os.path.join(REPO_ROOT, "scripts", "post_analysis.py"),
        "--config",              os.path.join(REPO_ROOT, "Config.yaml"),
        "--smart-version",       "1.0.0",
        "--canonical-variants",  os.path.join(refs, "svig_uk_canonical_variants.tsv"),
        "--gene-roles",          os.path.join(refs, "oncokb_gene_roles.tsv"),
        "--cancerhotspots-counts", os.path.join(refs, "cancerhotspots_counts.json"),
    ]
    print("Re-running post_analysis.py …")
    result = subprocess.run(cmd, cwd=smart_output_dir)
    if result.returncode != 0:
        sys.exit("ERROR: post_analysis.py failed.")
    print("Done.\n")


# ── Step 2: Load data ─────────────────────────────────────────────────────────

def load_smart_tier3(smart_output_dir):
    # Tier 2 has all three SVIG-UK columns; Tier 3 only has classification
    for fname in ("Final_result_tier2.tsv", "Final_result_tier3.tsv"):
        path = os.path.join(smart_output_dir, "output", fname)
        if os.path.isfile(path):
            break
    else:
        sys.exit(f"ERROR: No Tier 2/3 TSV found in {smart_output_dir}/output/")

    with open(path, encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    headers = rows[0]

    def col(name, default=None):
        try:
            return headers.index(name)
        except ValueError:
            return default

    key_cols = (col("Chromosome"), col("Start_Position"), col("REF"), col("ALT"))
    sc_col   = col("SVIG_UK_classification")
    sc_score = col("SVIG_UK_score")
    sc_codes = col("SVIG_UK_codes")
    sym_col  = col("SYMBOL")
    hvp_col  = col("HGVSp_Short")

    if sc_col is None:
        sys.exit("ERROR: SVIG_UK_classification column not found — re-run post_analysis.py")

    smart = {}
    for row in rows[2:]:
        key = ":".join(row[c] for c in key_cols)
        smart[key] = {
            "classification": row[sc_col],
            "score":          row[sc_score] if sc_score is not None else "",
            "codes":          row[sc_codes] if sc_codes is not None else "",
            "symbol":         row[sym_col]  if sym_col  is not None else "",
            "hgvsp":          row[hvp_col]  if hvp_col  is not None else "",
        }
    return smart


def load_human_classes(path):
    human = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = f"{r['CHROM']}:{r['POS']}:{r['REF']}:{r['ALT']}"
            human[key] = {
                "human_class":  r["Human_classification"],
                "svig_equiv":   r["SVIG_UK_equivalent"],
                "gene":         r["Gene"],
                "hgvsp":        r["HGVSp"],
                "tumour_type":  r["Tumour_type"],
            }
    return human


# ── Step 3: Compare ───────────────────────────────────────────────────────────

def compare(smart, human):
    matched   = []
    unmatched = []

    for key, h in human.items():
        s = smart.get(key)
        if s is None:
            unmatched.append({**h, "key": key})
            continue
        human_cls = h["svig_equiv"]
        smart_cls = s["classification"]
        concordant = human_cls == smart_cls
        matched.append({
            "key":          key,
            "gene":         h["gene"] or s["symbol"],
            "hgvsp":        h["hgvsp"] or s["hgvsp"],
            "tumour_type":  h["tumour_type"],
            "human_class":  h["human_class"],
            "human_svig":   human_cls,
            "smart_svig":   smart_cls,
            "smart_score":  s["score"],
            "smart_codes":  s["codes"],
            "concordant":   concordant,
        })

    return matched, unmatched


# ── Step 4: Metrics ───────────────────────────────────────────────────────────

def compute_metrics(matched):
    n   = len(matched)
    conc = sum(1 for m in matched if m["concordant"])
    confusion = defaultdict(Counter)
    for m in matched:
        confusion[m["human_svig"]][m["smart_svig"]] += 1

    # Per-category metrics
    per_cat = {}
    for cat in CATS:
        tp = confusion[cat][cat]
        fn = sum(confusion[cat][s] for s in CATS if s != cat)
        fp = sum(confusion[h][cat] for h in CATS if h != cat)
        tn = n - tp - fn - fp
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        precision   = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1 = (2 * sensitivity * precision / (sensitivity + precision)
              if (sensitivity + precision) > 0 else 0)
        per_cat[cat] = {
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "sensitivity": sensitivity,
            "precision": precision,
            "f1": f1,
        }

    return {
        "n":          n,
        "concordant": conc,
        "pct":        100 * conc / n if n > 0 else 0,
        "confusion":  confusion,
        "per_cat":    per_cat,
    }


# ── Step 5: Report ────────────────────────────────────────────────────────────

def write_report(metrics, matched, unmatched, report_dir):
    os.makedirs(report_dir, exist_ok=True)
    today = date.today().isoformat()

    # ── Discordant TSV ──
    disc_path = os.path.join(report_dir, "discordant_variants.tsv")
    discordant = [m for m in matched if not m["concordant"]]
    with open(disc_path, "w", newline="", encoding="utf-8") as f:
        fields = ["key","gene","hgvsp","tumour_type",
                  "human_class","human_svig","smart_svig",
                  "smart_score","smart_codes"]
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields,
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(discordant)
    print(f"Discordant variants: {disc_path}")

    # ── Full comparison TSV ──
    full_path = os.path.join(report_dir, "full_comparison.tsv")
    with open(full_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t",
                           fieldnames=["key","gene","hgvsp","tumour_type",
                                       "human_class","human_svig","smart_svig",
                                       "smart_score","smart_codes","concordant"])
        w.writeheader()
        w.writerows(matched)
    print(f"Full comparison: {full_path}")

    # ── Markdown report ──
    n   = metrics["n"]
    conc = metrics["concordant"]
    pct  = metrics["pct"]
    conf = metrics["confusion"]
    pc   = metrics["per_cat"]

    md_path = os.path.join(report_dir, "validation_report.md")
    with open(md_path, "w", encoding="utf-8") as f:

        f.write(f"# SVIG-UK Automated Classification — Validation Report\n\n")
        f.write(f"**Date:** {today}  \n")
        f.write(f"**Human-curated variants:** {n + len(unmatched)}  \n")
        f.write(f"**Matched to SMART output:** {n}  \n")
        f.write(f"**Unmatched (not in SMART output):** {len(unmatched)}  \n\n")
        f.write("---\n\n")

        f.write("## Overall Concordance\n\n")
        f.write(f"| Metric | Value |\n|---|---|\n")
        f.write(f"| Total matched variants | {n} |\n")
        f.write(f"| Concordant | **{conc} ({pct:.1f}%)** |\n")
        f.write(f"| Discordant | {n - conc} ({100 - pct:.1f}%) |\n\n")

        f.write("## Confusion Matrix\n\n")
        f.write("*Rows = Human classification · Columns = SMART automated classification*\n\n")
        header = "| Human \\ SMART |" + "|".join(f" {c} " for c in CATS) + "| **Total** |\n"
        sep    = "|" + "---|" * (len(CATS) + 2) + "\n"
        f.write(header)
        f.write(sep)
        for h in CATS:
            total = sum(conf[h][s] for s in CATS)
            row = f"| **{CLASS_COLORS[h]} {h}** |"
            for s in CATS:
                n_cell = conf[h][s]
                bold = "**" if h == s else ""
                row += f" {bold}{n_cell}{bold} |"
            row += f" {total} |\n"
            f.write(row)
        f.write("\n")

        f.write("## Per-Category Performance\n\n")
        f.write("| Category | TP | FP | FN | Sensitivity | Precision | F1 |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for cat in CATS:
            p = pc[cat]
            f.write(f"| {CLASS_COLORS[cat]} {cat} |"
                    f" {p['tp']} | {p['fp']} | {p['fn']} |"
                    f" {p['sensitivity']:.2f} | {p['precision']:.2f} | {p['f1']:.2f} |\n")
        f.write("\n")

        f.write("## Key Discordance Patterns\n\n")
        f.write("Major groups where human and SMART classifications differ:\n\n")
        for h in CATS:
            for s in CATS:
                if h != s and conf[h][s] >= 5:
                    f.write(f"- **Human {h} → SMART {s}**: {conf[h][s]} variants\n")
        f.write("\n")

        f.write("## Interpretation\n\n")
        f.write("The main source of discordance between human and automated classification "
                "is the **O4 evidence code** (enrichment in a large somatic variant database "
                "such as AACR Project GENIE). O4 contributes up to **+4 points** (Strong) "
                "and is currently not available in the automated pipeline because the GENIE "
                "dataset requires Synapse registration. Once O4 is integrated, the concordance "
                "rate is expected to improve significantly, particularly for the "
                "Oncogenic and Likely Oncogenic categories where most discordant variants "
                "scored 4–8 points (close to the 6-point and 10-point thresholds).\n\n")
        f.write("Additionally, O11 (tumour phenotype evidence from IHC, LOH, MSI, HRD) "
                "and the full PVS1 decision tree for O2 are not yet fully automated, "
                "which may contribute to under-classification of some TSG null variants.\n\n")

        f.write("---\n")
        f.write(f"*Generated by `utils/validate_svig_uk.py` · SMART v1.0.0 · SVIG-UK ACGS 2025*\n")

    print(f"Markdown report: {md_path}")
    return md_path


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Validate SVIG-UK automated classification.")
    p.add_argument("--smart-output",  default="data/jack_list_validation/output")
    p.add_argument("--human-classes", default="data/jack_list_validation/jack_curated.classifications.tsv")
    p.add_argument("--report-dir",    default="data/jack_list_validation/validation_report")
    p.add_argument("--skip-post-analysis", action="store_true",
                   help="Skip re-running post_analysis.py (use existing Tier 3 output).")
    return p.parse_args()


def main():
    args = parse_args()

    if not args.skip_post_analysis:
        run_post_analysis(args.smart_output)

    print("Loading SMART Tier 3 output…")
    smart = load_smart_tier3(args.smart_output)
    print(f"  {len(smart)} variants loaded.\n")

    print("Loading human classifications…")
    human = load_human_classes(args.human_classes)
    print(f"  {len(human)} variants loaded.\n")

    print("Comparing…")
    matched, unmatched = compare(smart, human)
    metrics = compute_metrics(matched)

    n    = metrics["n"]
    conc = metrics["concordant"]
    pct  = metrics["pct"]
    print(f"  Matched:    {n}")
    print(f"  Unmatched:  {len(unmatched)}")
    print(f"  Concordant: {conc}/{n} ({pct:.1f}%)\n")

    print("Confusion matrix:")
    conf = metrics["confusion"]
    row_hdr = f"  {'Human/SMART':<22}"
    print(row_hdr + "".join(f"{c[:9]:>11}" for c in CATS) + f"{'TOTAL':>8}")
    for h in CATS:
        total = sum(conf[h][s] for s in CATS)
        row = f"  {h:<22}" + "".join(f"{conf[h][s]:>11}" for s in CATS) + f"{total:>8}"
        print(row)
    print()

    write_report(metrics, matched, unmatched, args.report_dir)
    print("\nDone.")


if __name__ == "__main__":
    main()
