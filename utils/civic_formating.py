"""
civic_tsv_to_vcf.py  —  Convert a CIViC nightly-VariantSummaries.tsv to a
                         bgzip-compressed, tabix-indexed VCF for use with VEP
                         --custom annotation.

Usage
-----
    python civic_tsv_to_vcf.py \\
        --input  nightly-VariantSummaries.tsv \\
        --output civic.vcf.gz \\
        [--assembly grch38]          # grch38 (default) or grch37

Requirements
------------
    Python >= 3.8  (standard library only; no third-party packages needed)
    bgzip and tabix must be on PATH (or passed via --bgzip / --tabix)

Output
------
    <output>          — bgzip-compressed VCF
    <output>.tbi      — tabix index

Column mapping (CIViC nightly TSV → VCF)
-----------------------------------------
GRCh38 mode (default):
    CHROM  ← chromosome (stripped of 'chr' prefix if present)
    POS    ← start38   (1-based; CIViC is already 1-based)
    REF    ← reference_bases  (. → N for variants without explicit alleles)
    ALT    ← variant_bases    (. → <CNV> or similar placeholder)

GRCh37 mode:
    POS    ← start      (GRCh37 coordinate)

All remaining columns are written as semicolon-separated key=value pairs
in the INFO field so VEP can pull them with --custom.

Notes on CIViC TSV format
--------------------------
The nightly release contains both GRCh37 and GRCh38 coordinates.
The exact column names have been stable since ~2020:

    variant_id, gene, entrez_id, variant, summary, variant_types,
    hgvs_expressions, last_review_date, civic_variant_evidence_score,
    allele_registry_id, clinvar_ids, variant_aliases,
    chromosome, start, stop, reference_bases, variant_bases,        ← GRCh37
    chromosome2, start2, stop2,                                      ← GRCh37 secondary (fusions)
    representative_transcript, ensembl_version, reference_build,
    chromosome_grch38, start_grch38, stop_grch38,                   ← GRCh38
    chromosome2_grch38, start2_grch38, stop2_grch38,                ← GRCh38 secondary (fusions)
    representative_transcript2, gene_civic_url, variant_civic_url

If the nightly file has changed column names, the script will print the
actual header and raise an informative error rather than silently producing
garbage output.
"""

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
import tempfile
from datetime import date

try:
    from pyliftover import LiftOver
    _PYLIFTOVER_AVAILABLE = True
except ImportError:
    _PYLIFTOVER_AVAILABLE = False


# ---------------------------------------------------------------------------
# Accepted column name sets  (normalised to lowercase, underscores)
# We check both the new (grch38-suffixed) and old column variants.
# ---------------------------------------------------------------------------

CHROM_COLS_38  = ("chromosome_grch38", "chr_grch38", "chrom_grch38")
START_COLS_38  = ("start_grch38",)
STOP_COLS_38   = ("stop_grch38",)

CHROM_COLS_37  = ("chromosome", "chr", "chrom")
START_COLS_37  = ("start",)
STOP_COLS_37   = ("stop",)

REF_COLS       = ("reference_bases", "ref", "ref_bases")
ALT_COLS       = ("variant_bases",   "alt", "alt_bases", "var_bases")

# CSQ sub-field order must match the pipe-delimited description in Config.yaml:
# CIViC_CSQ: Allele|Consequence|SYMBOL|Entrez Gene ID|Feature_type|Feature|HGVSc|HGVSp|
#             CIViC Variant Name|CIViC Variant ID|CIViC Variant Aliases|CIViC Variant URL|
#             CIViC Molecular Profile Name|CIViC Molecular Profile ID|
#             CIViC Molecular Profile Aliases|CIViC Molecular Profile URL|CIViC HGVS|
#             Allele Registry ID|ClinVar IDs|CIViC Molecular Profile Score|
#             CIViC Entity Type|CIViC Entity ID|CIViC Entity URL|CIViC Entity Source|
#             CIViC Entity Variant Origin|CIViC Entity Status|CIViC Entity Significance|
#             CIViC Entity Direction|CIViC Entity Disease|CIViC Entity Therapies|
#             CIViC Entity Therapy Interaction Type|CIViC Evidence Phenotypes|
#             CIViC Evidence Level|CIViC Evidence Rating|CIViC Assertion ACMG Codes|
#             CIViC Assertion AMP Category|CIViC Assertion NCCN Guideline|
#             CIViC Assertion Regulatory Approval|CIViC Assertion FDA Companion Test
#
# The nightly-VariantSummaries TSV provides per-variant (not per-evidence) data,
# so evidence-level fields (Entity Type, Disease, Therapies, etc.) are left empty.
# The pipeline uses --custom ...,GN,VT,CSQ  (3 INFO fields extracted from this VCF).
CSQ_N_FIELDS = 39


def normalise_col(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "_", name.strip().lower())


def find_col(header_norm: list, candidates: tuple, label: str) -> int:
    for cand in candidates:
        if cand in header_norm:
            return header_norm.index(cand)
    raise ValueError(
        f"Cannot find required column '{label}'.\n"
        f"Tried: {candidates}\n"
        f"Available (normalised): {header_norm}"
    )


def sanitise_info_value(v: str) -> str:
    """Escape characters that would break VCF INFO syntax."""
    if not v or v in (".", "N/A", "NA", "None", ""):
        return "."
    # Remove semicolons (INFO field separator) and equals signs
    v = v.replace(";", "|").replace("=", ":").replace(" ", "_")
    # Truncate very long values (e.g. HGVS expressions) to keep VCF readable
    if len(v) > 500:
        v = v[:497] + "..."
    return v


def sanitise_csq_subfield(v: str) -> str:
    """Escape a value for use as a pipe-delimited CSQ sub-field.

    VEP replaces both '|' and ',' in embedded custom CSQ values with '&'
    when writing the main VCF CSQ field. A comma-separated multi-value
    within one sub-field (e.g. variant aliases, HGVS expressions, ClinVar IDs)
    would therefore be silently expanded into extra '&'-delimited fields,
    shifting all subsequent fields and corrupting the column mapping.
    Replace commas and pipes with '/' to keep each sub-field as a single token.
    """
    if not v or v in (".", "N/A", "NA", "None", ""):
        return ""
    return (v.replace("|", "/")
             .replace(",", "/")
             .replace(";", "/")
             .replace("=", ":")
             .replace(" ", "_"))[:300]


def build_csq(row: dict, alt: str) -> str:
    """
    Build the pipe-delimited CSQ INFO value from a CIViC nightly TSV row.

    Field order matches Config.yaml CIViC_CSQ description (39 fields):
      Allele|Consequence|SYMBOL|Entrez Gene ID|Feature_type|Feature|HGVSc|HGVSp|
      CIViC Variant Name|CIViC Variant ID|CIViC Variant Aliases|CIViC Variant URL|
      CIViC Molecular Profile Name|CIViC Molecular Profile ID|
      CIViC Molecular Profile Aliases|CIViC Molecular Profile URL|CIViC HGVS|
      Allele Registry ID|ClinVar IDs|CIViC Molecular Profile Score|
      CIViC Entity Type|CIViC Entity ID|CIViC Entity URL|CIViC Entity Source|
      CIViC Entity Variant Origin|CIViC Entity Status|CIViC Entity Significance|
      CIViC Entity Direction|CIViC Entity Disease|CIViC Entity Therapies|
      CIViC Entity Therapy Interaction Type|CIViC Evidence Phenotypes|
      CIViC Evidence Level|CIViC Evidence Rating|CIViC Assertion ACMG Codes|
      CIViC Assertion AMP Category|CIViC Assertion NCCN Guideline|
      CIViC Assertion Regulatory Approval|CIViC Assertion FDA Companion Test
    """
    s = sanitise_csq_subfield

    gene     = s(row.get("feature_name", "") or row.get("gene", ""))
    variant  = s(row.get("variant", ""))
    vid      = s(row.get("variant_id", ""))
    aliases  = s(row.get("variant_aliases", ""))
    vtype    = s(row.get("variant_types", ""))
    entrez   = s(row.get("entrez_id", ""))
    transcript = s(row.get("representative_transcript", ""))
    hgvs_raw = row.get("hgvs_descriptions", "") or ""
    allele_reg = s(row.get("allele_registry_id", ""))
    clinvar  = s(row.get("clinvar_ids", ""))
    mp_id    = s(row.get("single_variant_molecular_profile_id", ""))
    var_url  = s(row.get("variant_civic_url", ""))
    feat_type = s(row.get("feature_type", ""))

    # Derive HGVSc / HGVSp from the combined hgvs_descriptions field
    hgvsc = hgvsp = ""
    for expr in hgvs_raw.split(","):
        expr = expr.strip()
        if "c." in expr:
            hgvsc = s(expr)
        elif "p." in expr:
            hgvsp = s(expr)

    # Molecular profile name: CIViC uses "{GENE} {VARIANT}" for single-variant profiles
    mp_name = f"{gene}_{variant}" if gene and variant else ""
    mp_url  = ""  # not available in nightly TSV

    # Build the 39 pipe-separated fields (empty string = no data for that field)
    fields = [
        s(alt),        # 1  Allele
        vtype,         # 2  Consequence
        gene,          # 3  SYMBOL
        entrez,        # 4  Entrez Gene ID
        feat_type,     # 5  Feature_type
        transcript,    # 6  Feature
        hgvsc,         # 7  HGVSc
        hgvsp,         # 8  HGVSp
        variant,       # 9  CIViC Variant Name
        vid,           # 10 CIViC Variant ID
        aliases,       # 11 CIViC Variant Aliases
        var_url,       # 12 CIViC Variant URL
        mp_name,       # 13 CIViC Molecular Profile Name
        mp_id,         # 14 CIViC Molecular Profile ID
        "",            # 15 CIViC Molecular Profile Aliases  (not in TSV)
        mp_url,        # 16 CIViC Molecular Profile URL
        s(hgvs_raw),   # 17 CIViC HGVS
        allele_reg,    # 18 Allele Registry ID
        clinvar,       # 19 ClinVar IDs
        "",            # 20 CIViC Molecular Profile Score    (not in TSV)
        "",            # 21 CIViC Entity Type
        "",            # 22 CIViC Entity ID
        "",            # 23 CIViC Entity URL
        "",            # 24 CIViC Entity Source
        "",            # 25 CIViC Entity Variant Origin
        "",            # 26 CIViC Entity Status
        "",            # 27 CIViC Entity Significance
        "",            # 28 CIViC Entity Direction
        "",            # 29 CIViC Entity Disease
        "",            # 30 CIViC Entity Therapies
        "",            # 31 CIViC Entity Therapy Interaction Type
        "",            # 32 CIViC Evidence Phenotypes
        "",            # 33 CIViC Evidence Level
        "",            # 34 CIViC Evidence Rating
        "",            # 35 CIViC Assertion ACMG Codes
        "",            # 36 CIViC Assertion AMP Category
        "",            # 37 CIViC Assertion NCCN Guideline
        "",            # 38 CIViC Assertion Regulatory Approval
        "",            # 39 CIViC Assertion FDA Companion Test
    ]

    assert len(fields) == CSQ_N_FIELDS, f"CSQ field count mismatch: {len(fields)} != {CSQ_N_FIELDS}"
    return "|".join(fields)


def allele_ok(base: str) -> bool:
    """Return True if the base string is a valid VCF allele (not empty/dot)."""
    return bool(base) and base not in (".", "N/A", "NA", "None", "")


def convert(input_path: str, output_path: str, assembly: str,
            bgzip_bin: str, tabix_bin: str,
            chain_file: str = None) -> None:

    assembly = assembly.lower().replace("-", "")
    use_38 = assembly == "grch38"

    # -----------------------------------------------------------------------
    # Pass 1: read entire TSV into memory (file is ~600 KB, easily fits RAM)
    # -----------------------------------------------------------------------
    opener = gzip.open if input_path.endswith(".gz") else open
    with opener(input_path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        raw_header = reader.fieldnames
        if raw_header is None:
            raise ValueError("Input file appears to be empty.")

        header_norm = [normalise_col(c) for c in raw_header]

        print(f"[INFO] Detected {len(raw_header)} columns in TSV header.")
        print(f"[INFO] Assembly mode: {assembly}")

        # Locate required columns
        # If GRCh38 is requested but GRCh38-specific columns are absent,
        # try liftover from GRCh37 using a chain file (pyliftover).
        liftover = None
        if use_38:
            grch38_available = any(c in header_norm for c in CHROM_COLS_38)
            if not grch38_available:
                if chain_file and _PYLIFTOVER_AVAILABLE:
                    print(f"[INFO] GRCh38 columns not found — lifting over GRCh37→GRCh38 using {chain_file}")
                    liftover = LiftOver(chain_file)
                else:
                    if not _PYLIFTOVER_AVAILABLE:
                        print("[WARN] GRCh38 columns not found and pyliftover is not installed.")
                    elif not chain_file:
                        print("[WARN] GRCh38 columns not found and no --chain-file provided.")
                    print("[WARN] Falling back to GRCh37 coordinates (annotation will not work with a GRCh38 pipeline).")
                    use_38 = False

        if use_38 and liftover is None:
            chrom_idx = find_col(header_norm, CHROM_COLS_38, "chromosome_grch38")
            start_idx = find_col(header_norm, START_COLS_38, "start_grch38")
        else:
            chrom_idx = find_col(header_norm, CHROM_COLS_37, "chromosome")
            start_idx = find_col(header_norm, START_COLS_37, "start")

        try:
            ref_idx = find_col(header_norm, REF_COLS, "reference_bases")
        except ValueError:
            ref_idx = None

        try:
            alt_idx = find_col(header_norm, ALT_COLS, "variant_bases")
        except ValueError:
            alt_idx = None

        rows = list(reader)

    print(f"[INFO] Read {len(rows)} data rows from TSV.")

    # -----------------------------------------------------------------------
    # Pass 2: build VCF records
    # -----------------------------------------------------------------------
    vcf_records = []
    skipped = 0

    for row in rows:
        values = list(row.values())

        chrom = values[chrom_idx].strip() if values[chrom_idx] else ""
        start = values[start_idx].strip() if values[start_idx] else ""

        # Skip rows with missing coordinates
        if not chrom or chrom in (".", "N/A", "NA") or \
           not start or start in (".", "N/A", "NA") or not start.isdigit():
            skipped += 1
            continue

        # Normalise chromosome name: ensure 'chr' prefix for liftover, then strip
        chrom_bare = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
        pos = int(start)

        # Liftover GRCh37 → GRCh38 if requested
        if liftover is not None:
            result = liftover.convert_coordinate(f"chr{chrom_bare}", pos - 1)  # 0-based input
            if not result:
                skipped += 1
                continue
            _, pos38_0, strand, _ = result[0]
            chrom = re.sub(r"^chr", "", result[0][0], flags=re.IGNORECASE)
            pos = pos38_0 + 1  # back to 1-based
        else:
            chrom = chrom_bare

        ref = values[ref_idx].strip().upper() if ref_idx is not None else ""
        alt = values[alt_idx].strip().upper() if alt_idx is not None else ""

        # Handle missing alleles (CNVs, expression events, fusions, etc.)
        if not allele_ok(ref):
            ref = "N"
        if not allele_ok(alt):
            vt = row.get("variant_types", "") or "UNK"
            alt = "<" + re.sub(r"[^A-Za-z0-9_]", "_", vt)[:18] + ">"

        # Skip variants with no specific nucleotide allele — VEP exact-match cannot
        # match a symbolic allele (e.g. <missense_variant>) against a real ALT.
        # These are typically compound variants (e.g. S310F/Y) that CIViC represents
        # without a single nucleotide change in the nightly TSV.
        if alt.startswith("<"):
            skipped += 1
            continue

        # Build GN (gene name), VT (variant name), CSQ (structured annotation)
        gn  = sanitise_info_value(row.get("feature_name", "") or row.get("gene", ""))
        vt  = sanitise_info_value(row.get("variant", ""))
        csq = build_csq(row, alt)

        info_str = f"GN={gn};VT={vt};CSQ={csq}"

        vcf_records.append((chrom, pos, ref, alt, info_str))

    print(f"[INFO] {len(vcf_records)} records will be written "
          f"({skipped} skipped — missing {assembly} coordinates).")

    # -----------------------------------------------------------------------
    # Sort by chromosome then position
    # -----------------------------------------------------------------------
    def chrom_sort_key(rec):
        chrom = rec[0]
        # Numeric chromosomes sort as integers; X=23, Y=24, M/MT=25
        mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        if chrom in mapping:
            return (mapping[chrom], rec[1])
        try:
            return (int(chrom), rec[1])
        except ValueError:
            return (99, rec[1])

    vcf_records.sort(key=chrom_sort_key)

    # -----------------------------------------------------------------------
    # Write plain VCF to a temp file, then bgzip + tabix
    # -----------------------------------------------------------------------
    tmp_vcf = output_path.replace(".gz", "") + ".tmp.vcf"
    today = date.today().isoformat()

    with open(tmp_vcf, "w") as out:
        # VCF header
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={today}\n")
        out.write(f"##source=CIViC_nightly_VariantSummaries\n")
        out.write(f"##reference={assembly}\n")

        # INFO meta-lines — must match the fields extracted by VEP --custom ...,GN,VT,CSQ
        out.write('##INFO=<ID=GN,Number=1,Type=String,Description="CIViC gene name (HGNC symbol)">\n')
        out.write('##INFO=<ID=VT,Number=1,Type=String,Description="CIViC variant name">\n')
        csq_format = (
            "Allele|Consequence|SYMBOL|Entrez_Gene_ID|Feature_type|Feature|HGVSc|HGVSp|"
            "CIViC_Variant_Name|CIViC_Variant_ID|CIViC_Variant_Aliases|CIViC_Variant_URL|"
            "CIViC_Molecular_Profile_Name|CIViC_Molecular_Profile_ID|"
            "CIViC_Molecular_Profile_Aliases|CIViC_Molecular_Profile_URL|CIViC_HGVS|"
            "Allele_Registry_ID|ClinVar_IDs|CIViC_Molecular_Profile_Score|"
            "CIViC_Entity_Type|CIViC_Entity_ID|CIViC_Entity_URL|CIViC_Entity_Source|"
            "CIViC_Entity_Variant_Origin|CIViC_Entity_Status|CIViC_Entity_Significance|"
            "CIViC_Entity_Direction|CIViC_Entity_Disease|CIViC_Entity_Therapies|"
            "CIViC_Entity_Therapy_Interaction_Type|CIViC_Evidence_Phenotypes|"
            "CIViC_Evidence_Level|CIViC_Evidence_Rating|CIViC_Assertion_ACMG_Codes|"
            "CIViC_Assertion_AMP_Category|CIViC_Assertion_NCCN_Guideline|"
            "CIViC_Assertion_Regulatory_Approval|CIViC_Assertion_FDA_Companion_Test"
        )
        out.write(
            f'##INFO=<ID=CSQ,Number=.,Type=String,'
            f'Description="CIViC variant annotation. Format: {csq_format}">\n'
        )

        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for chrom, pos, ref, alt, info in vcf_records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n")

    print(f"[INFO] Plain VCF written to: {tmp_vcf}")

    # bgzip
    final_vcf = output_path if output_path.endswith(".gz") else output_path + ".gz"
    print(f"[RUN] bgzip -> {final_vcf}")
    subprocess.run([bgzip_bin, "-f", "-c", tmp_vcf],
                   stdout=open(final_vcf, "wb"), check=True)
    os.remove(tmp_vcf)

    # tabix
    print(f"[RUN] tabix -> {final_vcf}.tbi")
    subprocess.run([tabix_bin, "-p", "vcf", final_vcf], check=True)

    print(f"[DONE] Output: {final_vcf}")
    print(f"[DONE] Index:  {final_vcf}.tbi")


def main():
    parser = argparse.ArgumentParser(
        description="Convert CIViC nightly-VariantSummaries.tsv to a VCF for VEP --custom"
    )
    parser.add_argument("-i", "--input",  required=True,
                        help="Path to nightly-VariantSummaries.tsv (plain or .gz)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output path (e.g. civic.vcf.gz)")
    parser.add_argument("-a", "--assembly", default="grch38",
                        choices=["grch38", "grch37"],
                        help="Reference assembly to use for coordinates (default: grch38)")
    parser.add_argument("--bgzip", default="bgzip",
                        help="Path to bgzip binary (default: bgzip from PATH)")
    parser.add_argument("--tabix", default="tabix",
                        help="Path to tabix binary (default: tabix from PATH)")
    parser.add_argument("--chain-file", default=None,
                        help="hg19ToHg38 chain file for liftover (used when TSV lacks GRCh38 "
                             "columns and pyliftover is installed; e.g. hg19ToHg38.over.chain)")
    args = parser.parse_args()

    convert(
        input_path=args.input,
        output_path=args.output,
        assembly=args.assembly,
        bgzip_bin=args.bgzip,
        tabix_bin=args.tabix,
        chain_file=args.chain_file,
    )


if __name__ == "__main__":
    main()
