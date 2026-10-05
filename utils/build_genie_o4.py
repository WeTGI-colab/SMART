#!/usr/bin/env python3

# =============================================================================
# build_genie_o4.py
# =============================================================================
#
# DESCRIPTION:
#   Builds the AACR Project GENIE case-count table used by svig_uk_classifier.py
#   for SVIG-UK O4 ("Enriched in a somatic variant database ... international
#   (e.g. GENIE)") [S T1 O4].
#
#   For every mutation in data_mutations_extended.txt the builder
#     1. records the sample's OncoTree code, tissue and lineage (Myeloid / Lymphoid /
#        Solid / Unknown), so that case counting can be limited to 'on-target' tumour
#        types at the level the user chooses (code and subtypes, tissue, lineage, all);
#     2. lifts the GRCh37 coordinates of the public GENIE release to GRCh38 with
#        the UCSC chain file SMART already uses, checks the reference allele
#        against the GRCh38 FASTA and left-aligns indels, so that variants are
#        matched by genomic change rather than by GENIE's transcript numbering;
#     3. counts unique patients per variant: "multiple samples from the same
#        patient are only counted as a single entry" [S T1 O4]. A patient carrying the
#        variant in samples of different tumour types is one row listing all of them,
#        so the classifier counts the patient once if any of those samples is on-target.
#
#   GENIE is released under a data-use agreement that forbids redistribution.
#   The output is derived from GENIE: keep it with your licensed reference data
#   and never commit or share it.
#
# INPUT:
#   --genie-dir   directory with data_mutations_extended.txt and data_clinical_sample.txt
#   --oncotree    OncoTree JSON release (tumorTypes with code and tissue)
#   --chain       hg19ToHg38.over.chain(.gz)
#   --fasta       faidx-indexed GRCh38 FASTA
#
# OUTPUT (tab-separated, gzip), one row per GRCh38 variant x set of tumour types:
#   gene variant_class chrom start end ref alt ref_aa alt_aa oncotree_code tissue lineage n_patients
#   (oncotree_code, tissue and lineage are '|'-separated, aligned lists when a patient's
#    samples differ; tissue is empty for codes retired from the OncoTree release)
#   and <out>.cohort.tsv: unique patients per level (all, lineage, tissue, OncoTree code),
#   the denominators shown next to the counts in SMART's GENIE annotation columns.
#   (start/end/ref/alt in MAF style: 1-based, '-' for the empty allele;
#    ref_aa/alt_aa from GENIE's HGVSp_Short for single-residue changes)
# =============================================================================

import argparse
import bisect
import csv
import gzip
import json
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from svig_uk_classifier import Genome, left_align_maf  # noqa: E402

KEEP_CLASSES = {
    "Missense_Mutation": "missense", "Nonsense_Mutation": "truncating",
    "Frame_Shift_Del": "truncating", "Frame_Shift_Ins": "truncating",
    "In_Frame_Del": "inframe_del", "In_Frame_Ins": "inframe_ins",
    "Splice_Site": "splice", "Splice_Region": "splice", "Translation_Start_Site": "start",
    "Nonstop_Mutation": "nonstop",
}
HAEM_CANCER_TYPE_OVERRIDE = {  # organ-coded haematological neoplasms (as in the O4/O7 project)
    "Non-Hodgkin Lymphoma": "Lymphoid", "Hodgkin Lymphoma": "Lymphoid",
    "Mature B-Cell Neoplasms": "Lymphoid", "Mature T and NK Neoplasms": "Lymphoid",
}
SINGLE_RESIDUE = re.compile(r"^p\.([A-Z*])(\d+)([A-Z*])$")


class Chain:
    """Point liftover from a UCSC chain file (ungapped blocks, best chain first)."""

    def __init__(self, path):
        self.blocks = defaultdict(list)  # q_chrom -> [(q_start, q_end, t_chrom, t_start, t_strand, t_size)]
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            header = None
            for line in fh:
                f = line.split()
                if not f:
                    continue
                if f[0] == "chain":
                    t_name, t_size, t_strand, t_start = f[2], int(f[3]), f[4], int(f[5])
                    q_name, q_strand, q_start = f[7], f[9], int(f[10])
                    header = [t_name, t_size, t_strand, t_start, q_name, q_strand, q_start]
                    continue
                size = int(f[0])
                t_name, t_size, t_strand, t_pos, q_name, q_strand, q_pos = header
                if q_strand == "+":
                    self.blocks[t_name].append((t_pos, t_pos + size, q_name, q_pos, t_size))
                if len(f) == 3:
                    header[3] = t_pos + size + int(f[1])
                    header[6] = q_pos + size + int(f[2])
        # UCSC naming: the chain's "target" is hg19 (source) and "query" is hg38
        self.starts = {}
        for chrom, blocks in self.blocks.items():
            blocks.sort()
            self.starts[chrom] = [b[0] for b in blocks]

    def lift(self, chrom, pos):
        """1-based hg19 position -> (hg38 chrom, 1-based position) or None."""
        blocks = self.blocks.get(chrom)
        if not blocks:
            return None
        i = bisect.bisect_right(self.starts[chrom], pos - 1) - 1
        if i < 0:
            return None
        s, e, q_chrom, q_start, _ = blocks[i]
        if not s <= pos - 1 < e:
            return None
        return q_chrom, q_start + (pos - 1 - s) + 1


def load_lineages(genie_dir, oncotree_path):
    tissue = {t["code"]: t["tissue"] for t in json.load(open(oncotree_path))["tumorTypes"]}
    patient, sample_type = {}, {}
    with open(os.path.join(genie_dir, "data_clinical_sample.txt")) as fh:
        reader = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in reader:
            patient[r["SAMPLE_ID"]] = r["PATIENT_ID"]
            ct = r.get("CANCER_TYPE") or ""
            t = tissue.get(r.get("ONCOTREE_CODE") or "")
            if ct in HAEM_CANCER_TYPE_OVERRIDE:
                lin = HAEM_CANCER_TYPE_OVERRIDE[ct]
            elif t in ("Myeloid", "Lymphoid"):
                lin = t
            elif t or (ct and ct != "UNKNOWN"):
                lin = "Solid"
            else:
                lin = "Unknown"
            code = r.get("ONCOTREE_CODE") or ""
            sample_type[r["SAMPLE_ID"]] = (code, t or "", lin)
    return patient, sample_type


def cohort_sizes(patient, sample_type):
    """Unique patients per level; a patient with samples of several tumour types counts once in 'all'."""
    per = defaultdict(set)
    for sample, pid in patient.items():
        code, tissue, lineage = sample_type[sample]
        per[("all", "all")].add(pid)
        per[("lineage", lineage)].add(pid)
        if lineage in ("Myeloid", "Lymphoid"):
            per[("lineage", "Haematological")].add(pid)
        if tissue:
            per[("tissue", tissue)].add(pid)
        if code:
            per[("code", code)].add(pid)
    return {k: len(v) for k, v in per.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description="GENIE case counts for SVIG-UK O4")
    ap.add_argument("--genie-dir", required=True)
    ap.add_argument("--oncotree", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    patient, sample_type = load_lineages(args.genie_dir, args.oncotree)
    chain = Chain(args.chain)
    genome = Genome(args.fasta)
    carriers = defaultdict(lambda: defaultdict(set))  # variant key -> {patient: {(code, tissue, lineage)}}
    stats = defaultdict(int)

    with open(os.path.join(args.genie_dir, "data_mutations_extended.txt")) as fh:
        reader = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in reader:
            stats["rows"] += 1
            vclass = KEEP_CLASSES.get(r["Variant_Classification"])
            if vclass is None:
                continue
            if r.get("NCBI_Build") not in ("GRCh37", "37", "hg19"):
                stats["not GRCh37"] += 1
                continue
            sample = r["Tumor_Sample_Barcode"]
            if sample not in patient:
                stats["sample without clinical record"] += 1
                continue
            chrom = "chr" + r["Chromosome"].replace("chr", "").replace("MT", "M")
            start, end = int(r["Start_Position"]), int(r["End_Position"])
            ref, alt = r["Reference_Allele"], r["Tumor_Seq_Allele2"]
            a, b = chain.lift(chrom, start), chain.lift(chrom, end)
            if a is None or b is None or a[0] != b[0] or b[1] - a[1] != end - start:
                stats["liftover failed"] += 1
                continue
            chrom38, start38, end38 = a[0], a[1], b[1]
            if ref != "-":
                if genome.fetch(chrom38, start38 - 1, end38) != ref.upper():
                    stats["reference mismatch after liftover"] += 1
                    continue
            if ref == "-" or alt == "-":
                start38, end38, ref, alt = left_align_maf(genome, chrom38, start38, end38, ref, alt)
            m = SINGLE_RESIDUE.match(r.get("HGVSp_Short") or "")
            ref_aa, alt_aa = (m.group(1), m.group(3)) if m else ("", "")
            key = (r["Hugo_Symbol"], vclass, chrom38, start38, end38, ref, alt, ref_aa, alt_aa)
            carriers[key][patient[sample]].add(sample_type[sample])
            stats["kept"] += 1

    counts = defaultdict(int)
    for key, by_patient in carriers.items():
        for types in by_patient.values():
            types = sorted(types)
            if len(types) > 1:
                stats["patients with the variant in more than one tumour type"] += 1
            counts[key + tuple("|".join(t[i] for t in types) for i in range(3))] += 1
    with gzip.open(args.out, "wt") as out:
        out.write("gene\tvariant_class\tchrom\tstart\tend\tref\talt\tref_aa\talt_aa"
                  "\toncotree_code\ttissue\tlineage\tn_patients\n")
        for key in sorted(counts, key=lambda k: (k[0], k[2], k[3])):
            out.write("\t".join(map(str, key)) + f"\t{counts[key]}\n")
    cohort_path = args.out[:-len(".tsv.gz")] + ".cohort.tsv" if args.out.endswith(".tsv.gz") else args.out + ".cohort.tsv"
    with open(cohort_path, "w") as out:
        out.write("level\tvalue\tn_patients\n")
        for (level, value), n in sorted(cohort_sizes(patient, sample_type).items()):
            out.write(f"{level}\t{value}\t{n}\n")
    for k, v in stats.items():
        print(f"{k}: {v:,}", file=sys.stderr)
    print(f"wrote {len(counts):,} variant x tumour-type rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
