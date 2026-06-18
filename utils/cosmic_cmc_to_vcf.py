"""
cosmic_cmc_to_vcf.py  —  Convert a COSMIC Cancer Mutation Census (CMC) export to a
                          bgzip-compressed, tabix-indexed VCF for use with VEP
                          --custom annotation in SMART.

Why this script exists
----------------------
COSMIC data is NOT redistributable (see https://www.cosmickb.org/terms/), so SMART
cannot ship the COSMIC reference file. Instead, each user downloads the Cancer Mutation
Census themselves with their own (free for non-commercial research) COSMIC account, and
this script turns it into the small VCF the pipeline consumes. No COSMIC data is bundled
with SMART — only this transformation code.

How to get the input
--------------------
1. Register / log in at https://cancer.sanger.ac.uk/cosmic
2. Downloads -> project "Cancer Mutation Census" -> "All Data CMC" -> GRCh37*
   (the CMC is only released on GRCh37, but the file carries a GRCh38 coordinate
   column, which is what this script uses — so the output is GRCh38-native.)
3. You will get: CancerMutationCensus_AllData_Tsv_v<NN>_GRCh37.tar
   (containing CancerMutationCensus_AllData_v<NN>_GRCh37.tsv.gz)

Usage
-----
    python cosmic_cmc_to_vcf.py \\
        --input  CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar \\
        --output cosmic_cmc_grch38.vcf.gz

    # the .tar, the inner .tsv.gz, or a plain .tsv are all accepted as --input

Requirements
------------
    Python >= 3.8  (standard library only)
    bgzip and tabix on PATH (or passed via --bgzip / --tabix)
    Unix `sort` on PATH (used for low-memory coordinate sorting of the big file)

Output
------
    <output>          — bgzip-compressed, GRCh38-coordinate VCF
    <output>.tbi      — tabix index

    INFO fields (extracted by VEP with --custom ...,COSMIC,vcf,exact,0,GENE,CNT,TESTED,TIER,ONC_TSG,CGC_TIER,DNDS):
      GENE     <- GENE_NAME
      CNT      <- COSMIC_SAMPLE_MUTATED      (samples carrying this mutation = recurrence)
      TESTED   <- COSMIC_SAMPLE_TESTED       (samples tested at this locus)
      TIER     <- MUTATION_SIGNIFICANCE_TIER (1 / 2 / 3 / Other)
      ONC_TSG  <- ONC_TSG                    (oncogene / TSG, from Cancer Gene Census)
      CGC_TIER <- CGC_TIER                   (Cancer Gene Census tier 1 / 2)
      DNDS     <- DNDS_DISEASE_QVAL_SIG      (dN/dS driver significance)
    The COSV id (GENOMIC_MUTATION_ID) is written to the VCF ID column.

Scope / limitations
-------------------
    Only substitutions where REF and ALT are equal-length A/C/G/T strings (SNVs and
    MNVs) are emitted: those match VEP `exact` by position+allele with no reference
    genome needed. True indels are skipped (they require left-normalisation against the
    genome FASTA) and counted in the run summary. Rows without a GRCh38 coordinate are
    also skipped.
"""

import argparse
import csv
import gzip
import io
import os
import re
import subprocess
import sys
import tarfile
from datetime import date

# CMC column names we need (matched case-insensitively, non-alnum -> '_')
COL_CHROMPOS = "mutation_genome_position_grch38"
COL_COSV     = "genomic_mutation_id"
COL_WT       = "genomic_wt_allele_seq"
COL_MUT      = "genomic_mut_allele_seq"
COL_GENE     = "gene_name"
COL_CNT      = "cosmic_sample_mutated"
COL_TESTED   = "cosmic_sample_tested"
COL_TIER     = "mutation_significance_tier"
COL_ONC_TSG  = "onc_tsg"
COL_CGC_TIER = "cgc_tier"
COL_DNDS     = "dnds_disease_qval_sig"

_VALID_ALLELE = re.compile(r"^[ACGT]+$")


def normalise(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "_", name.strip().lower())


def open_cmc_text(path: str):
    """Return a text-mode line iterator for the CMC TSV, accepting .tar / .gz / .tsv."""
    if path.endswith(".tar"):
        tar = tarfile.open(path, "r")
        member = next((m for m in tar.getmembers()
                       if m.name.endswith(".tsv.gz") or m.name.endswith(".tsv")), None)
        if member is None:
            sys.exit(f"ERROR: no .tsv/.tsv.gz found inside {path}")
        fh = tar.extractfile(member)
        if member.name.endswith(".gz"):
            fh = gzip.GzipFile(fileobj=fh)
        return io.TextIOWrapper(fh, encoding="utf-8", newline=""), tar
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline=""), None
    return open(path, "rt", encoding="utf-8", newline=""), None


def find_idx(header_norm, col, label):
    if col not in header_norm:
        sys.exit(f"ERROR: required CMC column '{label}' not found.\n"
                 f"Header seen (normalised): {header_norm}")
    return header_norm.index(col)


def info_val(v: str) -> str:
    """Sanitise a value for a VCF INFO field.

    Besides ';' '=' and whitespace, commas must be replaced too: a comma in a VCF
    INFO value is the multi-value separator, so a value like "oncogene, fusion" would
    be split by VEP and silently truncated. Collapse it to a single token.
    """
    v = (v or "").strip()
    if v in ("", ".", "NS", "NA", "N/A", "None"):
        return "."
    return (v.replace(";", "|").replace("=", ":")
             .replace(", ", "/").replace(",", "/").replace(" ", "_"))[:200]


def convert(input_path, output_path, bgzip_bin, tabix_bin, sort_bin):
    fh, tar = open_cmc_text(input_path)
    try:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        hn = [normalise(c) for c in header]
        print(f"[INFO] CMC header: {len(header)} columns")

        i_pos   = find_idx(hn, COL_CHROMPOS, "Mutation genome position GRCh38")
        i_cosv  = find_idx(hn, COL_COSV,     "GENOMIC_MUTATION_ID")
        i_wt    = find_idx(hn, COL_WT,       "GENOMIC_WT_ALLELE_SEQ")
        i_mut   = find_idx(hn, COL_MUT,      "GENOMIC_MUT_ALLELE_SEQ")
        i_gene  = find_idx(hn, COL_GENE,     "GENE_NAME")
        i_cnt   = find_idx(hn, COL_CNT,      "COSMIC_SAMPLE_MUTATED")
        i_test  = find_idx(hn, COL_TESTED,   "COSMIC_SAMPLE_TESTED")
        i_tier  = find_idx(hn, COL_TIER,     "MUTATION_SIGNIFICANCE_TIER")
        i_onc   = find_idx(hn, COL_ONC_TSG,  "ONC_TSG")
        i_cgc   = find_idx(hn, COL_CGC_TIER, "CGC_TIER")
        i_dnds  = find_idx(hn, COL_DNDS,     "DNDS_DISEASE_QVAL_SIG")

        body_path = output_path.replace(".gz", "") + ".body.tmp"
        n_total = n_emit = n_nocoord = n_indel = n_badallele = 0

        with open(body_path, "w") as body:
            for row in reader:
                n_total += 1
                chrompos = row[i_pos].strip()
                if not chrompos or ":" not in chrompos:
                    n_nocoord += 1
                    continue
                chrom, _, span = chrompos.partition(":")
                start = span.split("-", 1)[0]
                chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
                if not start.isdigit():
                    n_nocoord += 1
                    continue

                ref = row[i_wt].strip().upper()
                alt = row[i_mut].strip().upper()
                if not ref or not alt:            # insertion/deletion (one side empty)
                    n_indel += 1
                    continue
                if not _VALID_ALLELE.match(ref) or not _VALID_ALLELE.match(alt):
                    n_badallele += 1              # non-ACGT (e.g. N, symbolic)
                    continue
                if len(ref) != len(alt):          # delins — needs FASTA normalisation
                    n_indel += 1
                    continue

                cosv = (row[i_cosv].strip() or ".")
                info = (
                    f"GENE={info_val(row[i_gene])};"
                    f"CNT={info_val(row[i_cnt])};"
                    f"TESTED={info_val(row[i_test])};"
                    f"TIER={info_val(row[i_tier])};"
                    f"ONC_TSG={info_val(row[i_onc])};"
                    f"CGC_TIER={info_val(row[i_cgc])};"
                    f"DNDS={info_val(row[i_dnds])}"
                )
                body.write(f"{chrom}\t{start}\t{cosv}\t{ref}\t{alt}\t.\t.\t{info}\n")
                n_emit += 1
    finally:
        fh.close()
        if tar is not None:
            tar.close()

    print(f"[INFO] rows read: {n_total}")
    print(f"[INFO] emitted (SNV/MNV): {n_emit}")
    print(f"[INFO] skipped — indels: {n_indel}, no GRCh38 coord: {n_nocoord}, "
          f"bad allele: {n_badallele}")
    if n_emit == 0:
        os.remove(body_path)
        sys.exit("ERROR: no records emitted — check the input file.")

    # coordinate-sort the body with unix sort (low memory for the big CMC file)
    sorted_path = body_path + ".sorted"
    print("[RUN] sort (by chrom, pos)")
    subprocess.run([sort_bin, "-k1,1", "-k2,2n", body_path, "-o", sorted_path], check=True)
    os.remove(body_path)

    # assemble header + sorted body
    plain_vcf = output_path.replace(".gz", "") + ".tmp.vcf"
    today = date.today().isoformat()
    with open(plain_vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={today}\n")
        out.write("##source=COSMIC_Cancer_Mutation_Census\n")
        out.write("##reference=GRCh38\n")
        out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="COSMIC gene name">\n')
        out.write('##INFO=<ID=CNT,Number=1,Type=Integer,Description="COSMIC samples carrying this mutation (COSMIC_SAMPLE_MUTATED)">\n')
        out.write('##INFO=<ID=TESTED,Number=1,Type=Integer,Description="COSMIC samples tested at this locus (COSMIC_SAMPLE_TESTED)">\n')
        out.write('##INFO=<ID=TIER,Number=1,Type=String,Description="CMC Mutation Significance Tier (1/2/3/Other)">\n')
        out.write('##INFO=<ID=ONC_TSG,Number=1,Type=String,Description="Gene role: oncogene/TSG (Cancer Gene Census)">\n')
        out.write('##INFO=<ID=CGC_TIER,Number=1,Type=String,Description="Cancer Gene Census tier (1/2)">\n')
        out.write('##INFO=<ID=DNDS,Number=1,Type=String,Description="dN/dS disease q-value significance (driver signal)">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        with open(sorted_path) as body:
            for line in body:
                out.write(line)
    os.remove(sorted_path)

    final_vcf = output_path if output_path.endswith(".gz") else output_path + ".gz"
    print(f"[RUN] bgzip -> {final_vcf}")
    with open(final_vcf, "wb") as fout:
        subprocess.run([bgzip_bin, "-f", "-c", plain_vcf], stdout=fout, check=True)
    os.remove(plain_vcf)

    print(f"[RUN] tabix -> {final_vcf}.tbi")
    subprocess.run([tabix_bin, "-p", "vcf", final_vcf], check=True)

    print(f"[DONE] {final_vcf}  ({n_emit} variants)")
    print(f"[DONE] {final_vcf}.tbi")


def main():
    p = argparse.ArgumentParser(
        description="Convert a COSMIC Cancer Mutation Census export to a GRCh38 VCF for VEP --custom")
    p.add_argument("-i", "--input", required=True,
                   help="CancerMutationCensus_AllData_*.tar (or the inner .tsv.gz / .tsv)")
    p.add_argument("-o", "--output", required=True,
                   help="Output VCF path, e.g. cosmic_cmc_grch38.vcf.gz")
    p.add_argument("--bgzip", default="bgzip", help="bgzip binary (default: from PATH)")
    p.add_argument("--tabix", default="tabix", help="tabix binary (default: from PATH)")
    p.add_argument("--sort",  default="sort",  help="sort binary (default: from PATH)")
    args = p.parse_args()
    convert(args.input, args.output, args.bgzip, args.tabix, args.sort)


if __name__ == "__main__":
    main()
