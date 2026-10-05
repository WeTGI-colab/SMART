#!/usr/bin/env python3

# =============================================================================
# genie_annotation.py
# =============================================================================
#
# DESCRIPTION:
#   Adds AACR Project GENIE case counts to a SMART per-sample MAF, as annotation
#   (Tiers 1-3), independently of the Tier 4 classification. Counts are unique
#   patients, matched by genomic change with residues computed on the transcript
#   SMART selected, using the table built by utils/build_genie_o4.py:
#
#     GENIE_SAME_NUCLEOTIDE_*   exact nucleotide change
#     GENIE_O4_MATCH_*          the entries SVIG-UK O4 counts for this variant class
#                               (same amino acid change; frameshift/nonsense at the same
#                               residue or downstream; exact insertion or nested deletions;
#                               exact change for canonical splice variants)
#     *_ALL / *_HAEM / *_SOLID  all tumour types / haematological / solid
#     GENIE_O4_MATCH_TYPE       what was counted for GENIE_O4_MATCH_*
#     GENIE_TOP_TUMOUR_TYPES    OncoTree codes with most patients for the O4 match
#     GENIE_MULTI_TUMOUR_PATIENTS  patients carrying the match in samples of more than one
#                               tumour type (counted once in ALL)
#     GENIE_COHORT              patients in the GENIE release (denominators)
#
#   A patient is counted once per column even if several samples carry the variant.
#   GENIE data are licence-restricted; the table must be built from the user's own release.
#
# USAGE:
#   genie_annotation.py --maf SAMPLE.maf --refs /refs [--out OUT.maf]
# =============================================================================

import argparse
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svig_uk_classifier as svig  # noqa: E402

HAEM = {"Myeloid", "Lymphoid"}
COLUMNS = ["GENIE_SAME_NUCLEOTIDE_ALL", "GENIE_SAME_NUCLEOTIDE_HAEM", "GENIE_SAME_NUCLEOTIDE_SOLID",
           "GENIE_O4_MATCH_TYPE", "GENIE_O4_MATCH_ALL", "GENIE_O4_MATCH_HAEM", "GENIE_O4_MATCH_SOLID",
           "GENIE_TOP_TUMOUR_TYPES", "GENIE_MULTI_TUMOUR_PATIENTS", "GENIE_COHORT"]


def split_counts(entries):
    """Unique-patient counts for all / haematological / solid tumours."""
    all_n = sum(e.n for e in entries)
    haem = sum(e.n for e in entries if any(l in HAEM for _, _, l in e.types))
    solid = sum(e.n for e in entries if any(l == "Solid" for _, _, l in e.types))
    return all_n, haem, solid


def top_tumour_types(entries, k=5):
    codes = Counter()
    for e in entries:
        for code, _, _ in e.types:
            codes[code or "unknown"] += e.n
    return "; ".join(f"{c} {n}" for c, n in codes.most_common(k))


def load_cohort(path):
    sizes = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sizes[(r["level"], r["value"])] = int(r["n_patients"])
    return sizes


def cohort_label(sizes):
    if not sizes:
        return ""
    return (f"all {sizes.get(('all', 'all'), 0):,}; haematological {sizes.get(('lineage', 'Haematological'), 0):,}; "
            f"solid {sizes.get(('lineage', 'Solid'), 0):,} patients")


def annotate_row(row, res, cohort):
    gene = row.get("Hugo_Symbol", "")
    entries = res.genie.by_gene.get(gene, [])
    out = {c: "" for c in COLUMNS}
    out["GENIE_COHORT"] = cohort
    allele = svig.maf_allele(row, res.genome)
    same = [e for e in entries if allele and (e.chrom, e.start, e.end, e.ref, e.alt) == allele]
    out["GENIE_SAME_NUCLEOTIDE_ALL"], out["GENIE_SAME_NUCLEOTIDE_HAEM"], out["GENIE_SAME_NUCLEOTIDE_SOLID"] = \
        split_counts(same)
    protein = svig.normalise_protein(row.get("HGVSp") or row.get("HGVSp_Short"))
    vclass = svig.variant_class(row.get("Consequence"))
    cls, what, matched, _ = svig.genie_matches(row, res, vclass, protein, entries)
    if cls is not None:
        out["GENIE_O4_MATCH_TYPE"] = what
        out["GENIE_O4_MATCH_ALL"], out["GENIE_O4_MATCH_HAEM"], out["GENIE_O4_MATCH_SOLID"] = split_counts(matched)
        out["GENIE_TOP_TUMOUR_TYPES"] = top_tumour_types(matched)
        out["GENIE_MULTI_TUMOUR_PATIENTS"] = sum(e.n for e in matched if len(e.types) > 1)
    else:
        out["GENIE_TOP_TUMOUR_TYPES"] = top_tumour_types(same)
        out["GENIE_MULTI_TUMOUR_PATIENTS"] = sum(e.n for e in same if len(e.types) > 1)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Add GENIE case counts to a SMART MAF")
    ap.add_argument("--maf", required=True)
    ap.add_argument("--out", help="output MAF (default: overwrite --maf)")
    ap.add_argument("--refs", required=True)
    ap.add_argument("--genie-o4", help="default: <refs>/GENIE/genie_o4_counts.tsv.gz")
    args = ap.parse_args(argv)

    genie_path = args.genie_o4 or os.path.join(args.refs, "GENIE", "genie_o4_counts.tsv.gz")
    if not os.path.exists(genie_path):
        print(f"GENIE counts not provided ({genie_path}); GENIE annotation skipped")
        return 0
    with open(args.maf, newline="") as fh:
        lines = fh.readlines()
    preamble = [l for l in lines if l.startswith("#")]
    reader = csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t")
    fieldnames, rows = list(reader.fieldnames or []), list(reader)

    res = svig.Resources()
    genes = {r.get("Hugo_Symbol", "") for r in rows}
    res.genie = svig.GenieO4(genie_path, genes=genes)
    models = os.path.join(args.refs, "TranscriptModels", "ncbiRefSeq.hg38.txt.gz")
    if os.path.exists(models):
        res.load_transcript_models(models, accessions={svig._strip_version(r.get(c)) for r in rows
                                                       for c in ("Feature", "MANE_SELECT", "NM_Transcript") if r.get(c)})
    fasta = os.path.join(args.refs, "liftover", "hg38.fa")
    if os.path.exists(fasta + ".fai"):
        res.genome = svig.Genome(fasta)
    cohort = cohort_label(load_cohort(genie_path[:-len(".tsv.gz")] + ".cohort.tsv"
                                      if genie_path.endswith(".tsv.gz") else genie_path + ".cohort.tsv"))

    fieldnames += [c for c in COLUMNS if c not in fieldnames]
    with open(args.out or args.maf, "w", newline="") as fh:
        fh.writelines(preamble)
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            r.update(annotate_row(r, res, cohort))
            w.writerow(r)
    print(f"GENIE annotation: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
