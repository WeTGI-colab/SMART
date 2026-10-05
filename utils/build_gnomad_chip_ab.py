#!/usr/bin/env python3

# =============================================================================
# build_gnomad_chip_ab.py
# =============================================================================
#
# DESCRIPTION:
#   Builds the gnomAD v4.1 allele-balance table used by svig_uk_classifier.py for
#   the SVIG-UK O3 CHIP correction [S T1 O3]: "In such genes, the 'allele balance'
#   chart in gnomAD can be used to assess whether entries in the gnomAD dataset may
#   represent CHIP contamination ... low variant allele frequency (VAF) alleles ...
#   can be excluded and the total allele frequency based on the number of alleles
#   with VAF between 0.3 and 0.65 can be re-calculated."
#
#   For every coding RefSeq transcript of the requested genes, coding exon parts +/-25 bp are
#   queried remotely (tabix over HTTPS) in the gnomAD v4.1 exome and genome sites
#   VCFs. Only variants with AC > 0 are written.
#
# REQUIREMENTS:
#   tabix (htslib) with HTTPS support; internet access to storage.googleapis.com.
#
# OUTPUT (tab-separated, gzip):
#   chrom pos ref alt dataset filter AC AN nhomalt ab_hist_alt_bin_freq
#   (ab_hist bins: 0.00|0.05|...|1.00, heterozygous individuals per bin)
# =============================================================================

import argparse
import gzip
import subprocess
import sys
import tempfile
import time

GNOMAD = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/{ds}/gnomad.{ds}.v4.1.sites.{chrom}.vcf.bgz"
# Genes named in [S T1 O3]; the guideline adds "(but not limited to)".
DEFAULT_GENES = ["DNMT3A", "ASXL1", "TET2", "TP53"]
FLANK = 25
TABIX_ATTEMPTS = 4


def gene_regions(models_path, genes):
    regions = {}
    with gzip.open(models_path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13 or f[12] not in genes or not f[1].startswith("NM_") or "_" in f[2]:
                continue
            starts = [int(x) for x in f[9].rstrip(",").split(",")]
            ends = [int(x) for x in f[10].rstrip(",").split(",")]
            cds_start, cds_end = int(f[6]), int(f[7])
            spans = regions.setdefault(f[12], (f[2], []))[1]
            for s, e in zip(starts, ends):
                s, e = max(s, cds_start), min(e, cds_end)  # coding part of the exon
                if s < e:
                    spans.append((s - FLANK + 1, e + FLANK))  # 1-based closed
    merged = {}
    for gene, (chrom, spans) in regions.items():
        out = []
        for s, e in sorted(spans):
            if out and s <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[gene] = (chrom, out)
    return merged


def info_value(info, key):
    for item in info.split(";"):
        if item.startswith(key + "="):
            return item.split("=", 1)[1]
    return ""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--transcript-models", required=True, help="UCSC ncbiRefSeq.txt.gz (hg38)")
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    regions = gene_regions(args.transcript_models, set(args.genes))
    missing = set(args.genes) - set(regions)
    if missing:
        sys.exit(f"genes without RefSeq models: {sorted(missing)}")
    n = 0
    index_dir = tempfile.mkdtemp(prefix="gnomad_tbi_")  # tabix saves remote indexes in its cwd
    with gzip.open(args.out, "wt") as out:
        out.write("chrom\tpos\tref\talt\tdataset\tfilter\tAC\tAN\tnhomalt\tab_hist_alt_bin_freq\n")
        for gene, (chrom, spans) in sorted(regions.items()):
            for ds in ("exomes", "genomes"):
                url = GNOMAD.format(ds=ds, chrom=chrom)
                cmd = ["tabix", url] + [f"{chrom}:{s}-{e}" for s, e in spans]
                for attempt in range(1, TABIX_ATTEMPTS + 1):
                    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=index_dir)
                    if proc.returncode == 0:
                        break
                    print(f"{gene} {ds}: tabix failed (attempt {attempt}): {proc.stderr.strip()[:200]}",
                          file=sys.stderr)
                    time.sleep(10 * attempt)
                else:
                    sys.exit(f"gnomAD query failed for {gene} {ds}")
                seen = set()
                for line in proc.stdout.splitlines():
                    f = line.split("\t", 8)
                    key = (f[1], f[3], f[4])
                    if key in seen:
                        continue  # overlapping regions return a record more than once
                    seen.add(key)
                    ac = info_value(f[7], "AC")
                    if not ac or int(ac) == 0:
                        continue
                    out.write("\t".join([f[0], f[1], f[3], f[4], ds, f[6], ac, info_value(f[7], "AN"),
                                         info_value(f[7], "nhomalt"),
                                         info_value(f[7], "ab_hist_alt_bin_freq")]) + "\n")
                    n += 1
                print(f"{gene} {ds}: {len(seen)} records", file=sys.stderr)
    print(f"wrote {n} variants -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
