#!/usr/bin/env python3

# =============================================================================
# build_svig_uk_resources.py
# =============================================================================
#
# DESCRIPTION:
#   Builds the two licence-gated lookup tables used by
#   scripts/svig_uk_classifier.py from the user's own COSMIC downloads. COSMIC
#   data cannot be redistributed, so these tables are never shipped with SMART;
#   each user builds them once from their licensed files.
#
#   OUTPUT FILES (written to --out-dir):
#
#     - cgc_gene_roles.tsv
#         Gene mode of action from the COSMIC Cancer Gene Census
#         (ROLE_IN_CANCER). SVIG-UK names the CGC as the resource for mode of
#         action (guideline §7.2; Supplementary Table 2, B2).
#
#     - cosmic_o4_counts.tsv.gz
#         Somatic case counts for evidence code O4, aggregated to the protein
#         level as SVIG-UK requires: missense entries are counted per amino-acid
#         change "regardless of the underlying nucleotide change"
#         (Supplementary Table 1, O4 — missense case counting). The count is
#         COSMIC_SAMPLE_MUTATED across all tissues: the CMC DISEASE breakdown is
#         empty for ~99.5% of mutations (v104), so on-target tissue counting is
#         not possible from this file.
#
# INPUT:
#   --cgc-tar   Cosmic_CancerGeneCensus_Tsv_v*_GRCh38.tar
#   --cmc-tar   CancerMutationCensus_AllData_Tsv_v*_GRCh37.tar
#               (the GRCh37 release is fine: O4 counts are protein-level)
#
# USAGE:
#   python3 utils/build_svig_uk_resources.py \
#       --cgc-tar /refs/COSMIC/Cosmic_CancerGeneCensus_Tsv_v104_GRCh38.tar \
#       --cmc-tar /refs/COSMIC/CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar \
#       --out-dir /refs/SVIG_UK
# =============================================================================

import argparse
import csv
import gzip
import io
import os
import re
import sys
import tarfile
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from svig_uk_classifier import normalise_protein  # noqa: E402

# COSMIC "Mutation Description AA" -> O4 case-counting class
MUT_CLASS = {
    "Substitution - Missense": "missense",
    "Substitution - Nonsense": "truncating",
    "Deletion - Frameshift": "truncating",
    "Insertion - Frameshift": "truncating",
    "Complex - frameshift": "truncating",
    "Frameshift": "truncating",
    "Deletion - In frame": "inframe_del",
    "Complex - deletion inframe": "inframe_del",
    "Insertion - In frame": "inframe_ins",
    "Complex - insertion inframe": "inframe_ins",
}


def _open_member(tar_path, suffix):
    tar = tarfile.open(tar_path)
    member = next(m for m in tar.getmembers() if m.name.endswith(suffix))
    raw = tar.extractfile(member)
    return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8", errors="replace")


def build_cgc(cgc_tar, out_path):
    fh = _open_member(cgc_tar, ".tsv.gz")
    reader = csv.DictReader(fh, delimiter="\t")
    n = 0
    with open(out_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene", "role_in_cancer", "tier", "tumour_types_somatic"])
        for r in reader:
            w.writerow([r["GENE_SYMBOL"], r.get("ROLE_IN_CANCER", ""),
                        r.get("TIER", ""), r.get("TUMOUR_TYPES_SOMATIC", "")])
            n += 1
    print(f"cgc_gene_roles.tsv: {n} genes")


def build_o4(cmc_tar, out_path):
    fh = _open_member(cmc_tar, ".tsv.gz")
    reader = csv.DictReader(fh, delimiter="\t")
    # key: (gene, transcript, protein, class, start, end) -> samples mutated
    agg = defaultdict(int)
    seen = 0
    for r in reader:
        seen += 1
        cls = MUT_CLASS.get(r.get("Mutation Description AA", ""))
        if not cls:
            continue
        protein = normalise_protein(r.get("Mutation AA", ""))
        if not protein:
            continue
        try:
            n = int(r.get("COSMIC_SAMPLE_MUTATED") or 0)
        except ValueError:
            continue
        if n == 0:
            continue
        key = (r["GENE_NAME"].split("_")[0], r.get("ACCESSION_NUMBER", ""), protein, cls,
               r.get("AA_MUT_START", ""), r.get("AA_MUT_STOP", ""))
        agg[key] += n
    with gzip.open(out_path, "wt", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene", "transcript", "protein", "mut_class", "aa_start", "aa_end", "n_samples"])
        for (gene, tx, protein, cls, s, e), n in sorted(agg.items()):
            w.writerow([gene, tx, protein, cls, s, e, n])
    print(f"cosmic_o4_counts.tsv.gz: {len(agg)} protein-level entries from {seen} CMC rows")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cgc-tar", required=True)
    ap.add_argument("--cmc-tar", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    build_cgc(args.cgc_tar, os.path.join(args.out_dir, "cgc_gene_roles.tsv"))
    build_o4(args.cmc_tar, os.path.join(args.out_dir, "cosmic_o4_counts.tsv.gz"))


if __name__ == "__main__":
    main()
