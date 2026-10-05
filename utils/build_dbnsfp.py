#!/usr/bin/env python3

# =============================================================================
# build_dbnsfp.py
# =============================================================================
#
# DESCRIPTION:
#   Builds the tabix-indexed dbNSFP file that SMART passes to VEP's dbNSFP plugin
#   (optional annotation, like CADD). dbNSFP academic releases are licensed
#   CC BY-NC-ND 4.0 and must be downloaded by each user from https://www.dbnsfp.org;
#   the built file must not be redistributed.
#
#   Only in silico predictors that SMART does not already annotate are kept, one
#   entry per tool, so no score is repeated:
#     - excluded because SMART already has them: SIFT, SIFT4G (VEP SIFT), PolyPhen-2
#       HDIV/HVAR (VEP PolyPhen), REVEL (REVEL plugin), CADD (CADD plugin);
#     - one version per tool: VARITY_R (not VARITY_ER or the leave-one-out models),
#       BayesDel without allele frequency (population frequency is annotated
#       separately from gnomAD), Eigen (not Eigen-PC), GERP++ (not GERP 91 mammals),
#       phyloP and phastCons on the 100-way vertebrate alignment only;
#     - the raw score and, where dbNSFP provides one, the categorical prediction;
#       rank scores are rescaled copies of the same score and are left out.
#   The position, allele, amino-acid and Ensembl identifier columns required by the
#   plugin are kept.
#
# INPUT:
#   --dbnsfp-dir  directory with dbNSFP5.0a_variant.chr*.gz (from the release zip)
#
# OUTPUT:
#   --out         bgzipped, tabix-indexed file (e.g. /refs/dbNSFP/dbNSFP5.0a_grch38.smart.gz)
# =============================================================================

import argparse
import glob
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

PLUGIN_COLUMNS = ["#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt", "genename",
                  "Ensembl_geneid", "Ensembl_transcriptid", "Ensembl_proteinid"]

PREDICTOR_COLUMNS = [
    "MutationTaster_score", "MutationTaster_pred",
    "MutationAssessor_score", "MutationAssessor_pred",
    "PROVEAN_score", "PROVEAN_pred",
    "VEST4_score",
    "MetaSVM_score", "MetaSVM_pred",
    "MetaLR_score", "MetaLR_pred",
    "MetaRNN_score", "MetaRNN_pred",
    "M-CAP_score", "M-CAP_pred",
    "MutPred_score",
    "MVP_score",
    "gMVP_score",
    "MPC_score",
    "PrimateAI_score", "PrimateAI_pred",
    "DEOGEN2_score", "DEOGEN2_pred",
    "BayesDel_noAF_score", "BayesDel_noAF_pred",
    "ClinPred_score", "ClinPred_pred",
    "LIST-S2_score", "LIST-S2_pred",
    "VARITY_R_score",
    "ESM1b_score", "ESM1b_pred",
    "AlphaMissense_score", "AlphaMissense_pred",
    "PHACTboost_score",
    "MutFormer_score",
    "MutScore_score",
    "Aloft_pred", "Aloft_Confidence",
    "DANN_score",
    "fathmm-XF_coding_score", "fathmm-XF_coding_pred",
    "Eigen-phred_coding",
    "GERP++_RS",
    "phyloP100way_vertebrate",
    "phastCons100way_vertebrate",
]

# Predictors already annotated by SMART through other sources (never selected).
ALREADY_IN_SMART = re.compile(r"^(SIFT|SIFT4G|Polyphen2|REVEL|CADD)_")

CHROM_ORDER = [str(c) for c in range(1, 23)] + ["X", "Y", "M"]


def select_columns(header):
    """Indices (1-based, for cut) of the plugin and predictor columns in a dbNSFP header."""
    names = header.rstrip("\n").split("\t")
    wanted = PLUGIN_COLUMNS + PREDICTOR_COLUMNS
    missing = [c for c in wanted if c not in names]
    if missing:
        raise ValueError(f"columns not found in this dbNSFP release: {missing}")
    clash = [c for c in PREDICTOR_COLUMNS if ALREADY_IN_SMART.match(c)]
    if clash:
        raise ValueError(f"predictors already annotated by SMART: {clash}")
    return [names.index(c) + 1 for c in wanted]


def chromosome_files(directory):
    files = {}
    for path in glob.glob(os.path.join(directory, "dbNSFP*_variant.chr*.gz")):
        m = re.search(r"_variant\.chr([0-9XYM]+)\.gz$", path)
        if m:
            files[m.group(1)] = path
    return [files[c] for c in CHROM_ORDER if c in files]


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build SMART's dbNSFP plugin file")
    ap.add_argument("--dbnsfp-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    for tool in ("bgzip", "tabix", "cut", "awk"):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} not found on PATH")
    files = chromosome_files(args.dbnsfp_dir)
    if not files:
        sys.exit(f"no dbNSFP variant files in {args.dbnsfp_dir}")
    with gzip.open(files[0], "rt") as fh:
        header = fh.readline()
    fields = ",".join(map(str, select_columns(header)))

    tmp = tempfile.mkdtemp(prefix="dbnsfp_", dir=os.path.dirname(os.path.abspath(args.out)))
    body = os.path.join(tmp, "body.tsv")
    try:
        with open(body, "w") as out:
            out.write("\t".join(PLUGIN_COLUMNS + PREDICTOR_COLUMNS) + "\n")
        for path in files:
            print(f">>> {os.path.basename(path)}", file=sys.stderr, flush=True)
            # Select columns first (dbNSFP lines are very long), drop the header, and check that
            # positions ascend: release files are sorted by position, so no sort is needed.
            check = ("awk -F'\\t' 'NR==1{next} {if ($2+0 < p) exit 3; p=$2+0; print}'")
            cmd = f"gzip -dc '{path}' | cut -f {fields} | {check} >> '{body}'"
            done = subprocess.run(["bash", "-o", "pipefail", "-c", cmd])
            if done.returncode == 3:
                raise SystemExit(f"{path} is not sorted by position; sort it before building")
            done.check_returncode()
        threads = str(max(1, (os.cpu_count() or 2) - 1))
        subprocess.run(["bash", "-o", "pipefail", "-c", f"bgzip -@ {threads} -c '{body}' > '{args.out}'"], check=True)
        subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", args.out], check=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"wrote {args.out} ({len(PREDICTOR_COLUMNS)} predictor columns)", file=sys.stderr)


if __name__ == "__main__":
    main()
