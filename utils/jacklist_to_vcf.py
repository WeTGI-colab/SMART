#!/usr/bin/env python3
"""
jacklist_to_vcf.py — Convert a curated variant list to a bgzip-compressed,
tabix-indexed VCF ready for SMART.

Two export formats are recognised automatically:

  * the original `jack_list.txt`: tab-separated, **no header**, 18 fixed columns;
  * a header-based export (comma- or tab-separated), e.g. `jack2.csv` with
    VariantChr,VariantPosition,VariantREF,VariantALT,Gene,Transcript,HGVSc,HGVSp,
    VariantType,VariantEffect,Classification,CreatedDate,CreatedBy,CuratedDate,
    CuratedBy,VariantStatus

Usage:
    python3 utils/jacklist_to_vcf.py \
        --input  jack2.csv \
        --output data/jack2_validation/jack2.vcf.gz \
        --transcripts data/jack2_validation/transcripts.txt \
        --tumour-type Myeloid

The output VCF also writes a sidecar TSV (same name, .classifications.tsv) with the
human classification for each variant, keyed by the same CHROM:POS:REF:ALT identifier
used in the SMART output, plus the curation provenance (transcript, dates, analyst,
reviewer, status, source file).

Several exports can be given at once: they are concatenated and deduplicated on
CHROM:POS:REF:ALT, keeping the most recently reviewed row. Variants curated twice with
different classifications are reported so they can be resolved before use.

    --tumour-type TYPE    tumour type for exports that do not carry one (affects O4)
    --status Curated      keep only signed-off rows (default: keep everything)
    --since 01/06/2026    keep rows classified on or after this date
    --transcripts FILE    write the transcript list for --transcripts-file
    --keep-duplicates     one VCF record per input row, as in the first version

Required fields are chrom, pos, ref, alt and the classification; everything else is
provenance, except the tumour type and the transcript, which change the result.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import Counter

# Map human classification labels to SVIG-UK equivalent
CLASS_MAP = {
    "ONC":      "Oncogenic",
    "LONC":     "Likely Oncogenic",
    "VUS":      "VUS",
    "VUS(HOT)": "VUS",
    "LBEN":     "Likely Benign",
    "BEN":      "Benign",
}

# Legacy jack_list.txt: field -> column index
LEGACY_COLUMNS = {
    "tumour_type": 1, "chrom": 2, "pos": 3, "ref": 4, "alt": 5, "gene": 6,
    "transcript": 7, "hgvsc": 8, "hgvsp": 9, "vtype": 10, "consequence": 11,
    "human_class": 12, "date_class": 13, "analyst": 14, "date_review": 15,
    "reviewer": 16, "status": 17,
}

# Header-based exports: field -> accepted column names (normalised)
HEADER_ALIASES = {
    "chrom":       ("variantchr", "chrom", "chr", "chromosome"),
    "pos":         ("variantposition", "pos", "position", "start", "startposition"),
    "ref":         ("variantref", "ref", "referenceallele"),
    "alt":         ("variantalt", "alt", "alternateallele"),
    "gene":        ("gene", "genesymbol", "hugosymbol"),
    "transcript":  ("transcript", "nmtranscript", "refseq", "refseqtranscript"),
    "hgvsc":       ("hgvsc", "cdna", "hgvscoding"),
    "hgvsp":       ("hgvsp", "protein", "hgvsprotein", "hgvspshort"),
    "vtype":       ("varianttype", "type"),
    "consequence": ("varianteffect", "consequence", "effect"),
    "human_class": ("classification", "class", "humanclassification"),
    "date_class":  ("createddate", "dateclassified", "dateclassification"),
    "analyst":     ("createdby", "analyst", "scientist"),
    "date_review": ("curateddate", "datereviewed", "reviewdate"),
    "reviewer":    ("curatedby", "reviewer", "checkedby"),
    "status":      ("variantstatus", "status"),
    "tumour_type": ("tumourtype", "tumortype", "diseasetype", "oncotreecode", "oncotree"),
}

NULLS = ("", "none", "na", "n/a", ".", "30/12/1899", "1899-12-30")


def clean(value):
    """Strip a cell; Excel nulls and placeholders become '.'."""
    v = str(value or "").strip()
    return "." if v.lower() in NULLS else v


def as_date(text):
    """DD/MM/YYYY (or YYYY-MM-DD) -> sortable YYYYMMDD; anything else sorts first."""
    m = re.match(r"^(\d{1,2})/(\d{1,2})/(\d{4})$", text or "")
    if m:
        return f"{m.group(3)}{int(m.group(2)):02d}{int(m.group(1)):02d}"
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", text or "")
    return f"{m.group(1)}{m.group(2)}{m.group(3)}" if m else "00000000"


def normalise(name):
    return re.sub(r"[^a-z0-9]", "", str(name or "").lower())


def is_header_format(first_line):
    """A header line is one naming the chromosome column rather than holding a value."""
    cells = {normalise(c) for c in re.split(r"[,\t]", first_line)}
    return bool(cells & set(HEADER_ALIASES["chrom"])) and bool(cells & set(HEADER_ALIASES["human_class"]))


def read_export(path, default_tumour):
    """Read one export in either format. Returns (rows, counters)."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        lines = fh.read().splitlines()
    lines = [l for l in lines if l.strip()]
    if not lines:
        return [], Counter()

    rows, counts = [], Counter()
    source = os.path.basename(path)

    if is_header_format(lines[0]):
        delim = "\t" if lines[0].count("\t") > lines[0].count(",") else ","
        reader = csv.DictReader(lines, delimiter=delim)
        index = {}
        for field, aliases in HEADER_ALIASES.items():
            for column in reader.fieldnames or []:
                if normalise(column) in aliases:
                    index[field] = column
                    break
        missing = [f for f in ("chrom", "pos", "ref", "alt", "human_class") if f not in index]
        if missing:
            sys.exit(f"{source}: required column(s) not found: {', '.join(missing)}\n"
                     f"  columns present: {', '.join(reader.fieldnames or [])}")
        if "tumour_type" not in index:
            counts["no_tumour_column"] = 1
        for record in reader:
            rows.append({f: clean(record.get(c)) for f, c in index.items()})
    else:
        for line in lines:
            cols = line.split("\t")
            rows.append({f: clean(cols[i]) if len(cols) > i else "."
                         for f, i in LEGACY_COLUMNS.items()})

    out = []
    for r in rows:
        if r.get("tumour_type", ".") == ".":
            r["tumour_type"] = default_tumour
        r["source"] = source
        for field in HEADER_ALIASES:
            r.setdefault(field, ".")
        try:
            r["chrom"] = "chr" + r["chrom"].lstrip("chr")
            r["pos"] = int(r["pos"])
            r["ref"] = r["ref"].upper()
            r["alt"] = r["alt"].upper()
        except (ValueError, AttributeError):
            counts["bad_coordinates"] += 1
            continue
        if r["human_class"] not in CLASS_MAP:
            counts[f"unknown class: {r['human_class']}"] += 1
            continue
        r["svig_equiv"] = CLASS_MAP[r["human_class"]]
        out.append(r)
    return out, counts


def deduplicate(variants):
    """One record per CHROM:POS:REF:ALT, keeping the most recently reviewed.

    Returns (kept, conflicts) where conflicts are variants curated more than once
    with different classifications — those need a human decision.
    """
    by_key = {}
    for v in variants:
        by_key.setdefault(f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}", []).append(v)
    kept, conflicts = [], []
    for key, rows in by_key.items():
        rows.sort(key=lambda v: (as_date(v["date_review"]), as_date(v["date_class"])), reverse=True)
        if len({r["human_class"] for r in rows}) > 1:
            conflicts.append((key, rows))
        kept.append(rows[0])
    return kept, conflicts


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
    p = argparse.ArgumentParser(description="Convert a curated variant list to a SMART-compatible VCF")
    p.add_argument("--input", required=True, nargs="+",
                   help="One or more curation exports (jack_list.txt format or a header-based CSV/TSV)")
    p.add_argument("--output", required=True, help="Output VCF path (e.g. curated.vcf.gz)")
    p.add_argument("--sample", default="CURATED", help="Sample name in VCF (default: CURATED)")
    p.add_argument("--tumour-type", default="unknown",
                   help="tumour type for exports without one (default: unknown, which makes O4 "
                        "count every tumour type)")
    p.add_argument("--status", help="keep only rows with this status (e.g. Curated)")
    p.add_argument("--since", help="keep only rows classified on or after DD/MM/YYYY")
    p.add_argument("--transcripts", help="write the transcript list for --transcripts-file")
    p.add_argument("--keep-duplicates", action="store_true",
                   help="one VCF record per input row instead of per variant")
    return p.parse_args()


def main():
    args = parse_args()
    cutoff = as_date(args.since) if args.since else None
    if args.since and cutoff == "00000000":
        sys.exit(f"--since must be DD/MM/YYYY, not {args.since!r}")

    variants, counts = [], Counter()
    for path in args.input:
        rows, file_counts = read_export(path, args.tumour_type)
        counts.update(file_counts)
        if file_counts.get("no_tumour_column"):
            print(f"NOTE: {os.path.basename(path)} has no tumour type column — using "
                  f"{args.tumour_type!r} for every row (this changes O4)")
        for r in rows:
            if args.status and r["status"] != args.status:
                counts["filtered by --status"] += 1
                continue
            if cutoff and as_date(r["date_class"]) < cutoff:
                counts["filtered by --since"] += 1
                continue
            variants.append(r)

    print(f"Loaded {len(variants)} rows from {len(args.input)} file(s)")
    for reason, n in sorted(counts.items()):
        if reason != "no_tumour_column":
            print(f"  skipped {n}: {reason}")
    if not variants:
        sys.exit("nothing left to convert")

    if not args.keep_duplicates:
        variants, conflicts = deduplicate(variants)
        print(f"{len(variants)} unique variants after deduplication")
        if conflicts:
            print(f"\nWARNING: {len(conflicts)} variant(s) curated more than once with "
                  f"different classifications — the most recent review was kept:")
            for key, rows in conflicts:
                calls = ", ".join(f"{r['human_class']} ({r['date_review']}, {r['source']})" for r in rows)
                print(f"  {key}  {rows[0]['gene']} {rows[0]['hgvsp']}: {calls}")

    # Sort by chrom (numeric then X/Y) then pos
    def sort_key(v):
        c = v["chrom"].replace("chr", "")
        try:
            return (0, int(c), v["pos"])
        except ValueError:
            return (1, c, v["pos"])

    variants.sort(key=sort_key)

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
            fh.write(f"{v['chrom']}\t{v['pos']}\t{var_id}\t"
                     f"{v['ref']}\t{v['alt']}\t.\tPASS\t{info}\t"
                     f"GT:AD:AF:DP\t0/1:50,25:0.33:75\n")
    print(f"VCF written: {vcf_plain}")

    out_gz = args.output if args.output.endswith(".gz") else args.output + ".gz"
    subprocess.run(["bgzip", "-f", vcf_plain], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", out_gz], check=True)
    print(f"Compressed + indexed: {out_gz}")

    tsv_out = out_gz.replace(".vcf.gz", ".classifications.tsv")
    with open(tsv_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["CHROM", "POS", "REF", "ALT", "Gene", "HGVSp", "Tumour_type",
                    "Human_classification", "SVIG_UK_equivalent",
                    "Transcript", "HGVSc", "Consequence",
                    "Date_classified", "Analyst", "Date_reviewed", "Reviewer",
                    "Status", "Source_file"])
        for v in variants:
            w.writerow([v["chrom"], v["pos"], v["ref"], v["alt"], v["gene"], v["hgvsp"],
                        v["tumour_type"], v["human_class"], v["svig_equiv"],
                        v["transcript"], v["hgvsc"], v["consequence"],
                        v["date_class"], v["analyst"], v["date_review"], v["reviewer"],
                        v["status"], v["source"]])
    print(f"Classification sidecar: {tsv_out}")

    if args.transcripts:
        wanted = sorted({v["transcript"] for v in variants if v["transcript"] != "."})
        with open(args.transcripts, "w") as fh:
            fh.write("\n".join(wanted) + "\n")
        print(f"Transcript list: {args.transcripts} ({len(wanted)} transcripts, "
              f"{len({v['gene'] for v in variants})} genes)")

    print("\nDistribution of human classifications:")
    for cls, n in Counter(v["human_class"] for v in variants).most_common():
        print(f"  {CLASS_MAP[cls]:20s} ({cls:8s})  {n:4d}")


if __name__ == "__main__":
    main()
