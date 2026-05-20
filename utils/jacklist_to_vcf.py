#!/usr/bin/env python3
"""
jacklist_to_vcf.py — Convert a curated variant list (jack_list.txt format)
to a bgzip-compressed, tabix-indexed VCF ready for SMART.

Usage:
    python3 utils/jacklist_to_vcf.py \
        --input  Files4ThisProject/jack_list.txt \
        --output data/jack_list_validation/jack_curated.vcf.gz \
        --sample CURATED

The output VCF also writes a sidecar TSV (same name, .classification.tsv)
containing the original human classification for each variant, keyed by
the same CHROM:POS:REF:ALT identifier used in the SMART output.

Input columns (tab-delimited, no header):
  0  variant_id      e.g. 1:114716123-C-T
  1  tumour_type     e.g. Myeloid
  2  chrom           e.g. 1
  3  pos             e.g. 114716123
  4  ref             e.g. C
  5  alt             e.g. T
  6  gene            e.g. NRAS
  7  transcript      e.g. NM_002524.5
  8  hgvsc           e.g. c.38G>A
  9  hgvsp           e.g. p.Gly13Asp  (or None)
  10 variant_type    e.g. SNP
  11 consequence     e.g. Missense
  12 classification  ONC | LONC | VUS | VUS(HOT) | LBEN | BEN
  13 date_classified
  14 analyst
  15 date_reviewed
  16 reviewer
  17 status          e.g. Curated
"""

import argparse
import csv
import os
import subprocess
import sys

# Map human classification labels to SVIG-UK equivalent
CLASS_MAP = {
    "ONC":      "Oncogenic",
    "LONC":     "Likely Oncogenic",
    "VUS":      "VUS",
    "VUS(HOT)": "VUS",
    "LBEN":     "Likely Benign",
    "BEN":      "Benign",
}

VCF_HEADER = """\
##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=TUMOUR_TYPE,Number=1,Type=String,Description="Tumour type from curation list">
##INFO=<ID=HUMAN_CLASS,Number=1,Type=String,Description="Human oncogenicity classification">
##INFO=<ID=SVIG_UK_EQUIV,Number=1,Type=String,Description="SVIG-UK equivalent of human classification">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def parse_args():
    p = argparse.ArgumentParser(description="Convert jack_list.txt to SMART-compatible VCF")
    p.add_argument("--input",  required=True, help="Path to jack_list.txt")
    p.add_argument("--output", required=True, help="Output VCF path (e.g. curated.vcf.gz)")
    p.add_argument("--sample", default="CURATED", help="Sample name in VCF (default: CURATED)")
    return p.parse_args()


def main():
    args = parse_args()

    # Read input
    variants = []
    skipped  = 0
    with open(args.input, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 13:
                skipped += 1
                continue
            cls = cols[12].strip()
            if cls not in CLASS_MAP:
                skipped += 1
                continue
            try:
                chrom = "chr" + cols[2].strip().lstrip("chr")
                pos   = int(cols[3].strip())
                ref   = cols[4].strip().upper()
                alt   = cols[5].strip().upper()
            except (ValueError, IndexError):
                skipped += 1
                continue
            variants.append({
                "chrom":       chrom,
                "pos":         pos,
                "ref":         ref,
                "alt":         alt,
                "gene":        cols[6].strip() if len(cols) > 6 else ".",
                "tumour_type": cols[1].strip() if len(cols) > 1 else ".",
                "human_class": cls,
                "svig_equiv":  CLASS_MAP[cls],
                "hgvsp":       cols[9].strip() if len(cols) > 9 and cols[9].strip() not in ("None", "") else ".",
            })

    # Sort by chrom (numeric then X/Y) then pos
    def sort_key(v):
        c = v["chrom"].replace("chr", "")
        try:
            return (0, int(c), v["pos"])
        except ValueError:
            return (1, c, v["pos"])

    variants.sort(key=sort_key)
    print(f"Loaded {len(variants)} variants ({skipped} skipped — unknown class or bad format)")

    # Write plain VCF
    vcf_plain = args.output.replace(".gz", "")
    os.makedirs(os.path.dirname(os.path.abspath(vcf_plain)) or ".", exist_ok=True)

    with open(vcf_plain, "w") as fh:
        fh.write(VCF_HEADER.format(sample=args.sample))
        for v in variants:
            var_id = f"{v['gene']}_{v['hgvsp']}"
            info = (f"TUMOUR_TYPE={v['tumour_type']};"
                    f"HUMAN_CLASS={v['human_class']};"
                    f"SVIG_UK_EQUIV={v['svig_equiv']}")
            # Fake allelic depth: 50 ref, 25 alt → VAF = 0.33
            fmt_val = "0/1:50,25:0.33:75"
            fh.write(f"{v['chrom']}\t{v['pos']}\t{var_id}\t"
                     f"{v['ref']}\t{v['alt']}\t.\tPASS\t{info}\t"
                     f"GT:AD:AF:DP\t{fmt_val}\n")

    print(f"VCF written: {vcf_plain}")

    # bgzip + tabix
    out_gz = args.output if args.output.endswith(".gz") else args.output + ".gz"
    subprocess.run(["bgzip", "-f", vcf_plain], check=True)
    subprocess.run(["tabix", "-p", "vcf", out_gz], check=True)
    print(f"Compressed + indexed: {out_gz}")

    # Write classification sidecar TSV
    tsv_out = out_gz.replace(".vcf.gz", ".classifications.tsv")
    with open(tsv_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["CHROM", "POS", "REF", "ALT", "Gene",
                    "HGVSp", "Tumour_type",
                    "Human_classification", "SVIG_UK_equivalent"])
        for v in variants:
            w.writerow([
                v["chrom"], v["pos"], v["ref"], v["alt"],
                v["gene"], v["hgvsp"], v["tumour_type"],
                v["human_class"], v["svig_equiv"],
            ])
    print(f"Classification sidecar: {tsv_out}")
    print(f"\nDistribution of human classifications:")
    from collections import Counter
    counts = Counter(v["human_class"] for v in variants)
    for cls, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {CLASS_MAP[cls]:20s} ({cls:8s})  {n:4d}")


if __name__ == "__main__":
    main()
