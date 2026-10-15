#!/usr/bin/env python3
# =============================================================================
# vcf2table.py
# =============================================================================
#
# DESCRIPTION:
#   Parses a VCF file annotated with VEP and OncoKB and converts it into a
#   structured tabular format for downstream analysis.
#
#   The script selects a single transcript per variant based on a provided
#   whitelist (NM_* transcripts), extracts relevant annotations, and expands
#   complex fields into a flat table suitable for reporting and integration
#   with other tools (e.g. MAF-based workflows).
#
#   It also performs additional transformations to produce standardized fields
#   commonly used in somatic variant analysis pipelines.
#
# INPUT:
#   - VCF file annotated with VEP and OncoKB
#   - Transcript whitelist file (NM_* identifiers)
#
# OUTPUT:
#   - Tab-separated table containing:
#       - Core VCF fields
#       - Selected VEP annotations
#       - Expanded INFO fields
#       - Expanded OncoKB JSON annotations
#       - Standardized fields for downstream use (e.g. MAF-like columns)
#
# PROCESSING STEPS:
#   1. Parse VCF header to extract CSQ annotation structure
#   2. Load transcript whitelist (NM_*)
#   3. Select one transcript per variant based on priority rules
#   4. Expand VEP annotations (CSQ field)
#   5. Parse and expand INFO fields
#   6. Convert data into a pandas DataFrame
#   7. Compute derived fields (End_Position, VAF-related fields)
#   8. Standardize output columns (MAF-like structure)
#   9. Expand OncoKB annotations (JSON fields into columns)
#   10. Write final table output
#
# USAGE:
#   python vcf2table.py \
#       --vcf input.vcf \
#       --transcripts transcripts.txt \
#       --output output.tsv
#
# DEPENDENCIES:
#   - Python 3.8+
#   - pandas
#
# AUTHOR:
#   Variant annotation pipeline — Manuel
# =============================================================================

import csv
import re
import sys
import argparse
import pandas as pd
import json
from pandas import json_normalize


# ------------------------------------------------------------
# Argument parser
# ------------------------------------------------------------
def get_args():
    """
    Parse command-line arguments.

    Defines required inputs for the script, including:
      - Input VCF file
      - Transcript whitelist file
      - Output file path
      - Optional debug flag

    Returns
    -------
    argparse.Namespace
        Parsed arguments object.
    """
    parser = argparse.ArgumentParser(
        description="Extract full VEP + OncoKB annotations for one selected transcript per variant"
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-t", "--transcripts", required=True, help="Transcript whitelist file (NM_*)")
    parser.add_argument("--tumour-column", "--tumor-column", dest="tumour_column",
                        help="name the tumour sample column, for a paired VCF whose columns "
                             "are named something other than TUMOR/NORMAL and which carries "
                             "no '##tumor_sample=' header")
    parser.add_argument("--normal-column", dest="normal_column",
                        help="name the matched normal sample column")
    parser.add_argument("-o", "--output", required=True, help="Final output file")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


# ------------------------------------------------------------
# Load NM transcript whitelist
# ------------------------------------------------------------
FIXED_VCF_COLUMNS = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
TUMOUR_NAMES = ("TUMOR", "TUMOUR", "TUMOR_SAMPLE", "T")
NORMAL_NAMES = ("NORMAL", "NORMAL_SAMPLE", "N", "GERMLINE")


def select_sample_columns(meta_lines, vcf_columns):
    """Which sample column is the tumour, and which the normal.

    SMART reports one sample per run, and every number it derives from a sample column --
    GT, AD, DP, VAF -- describes the tumour. It used to take whichever column came first,
    with no check that there was only one: for a Strelka or SomaticSniper VCF, where the
    normal precedes the tumour, that reported the normal's genotype and allele fraction for
    a somatic variant. A variant at 55% in the tumour came out at 2% with GT 0/0, which
    reads as a subclonal call rather than as the wrong sample.

    Resolved in order of authority:

      1. one sample column -- it is the tumour, which is the tumour-only case and the
         laboratory's routine DRAGEN input;
      2. `##tumor_sample=` / `##normal_sample=` in the header. GATK Mutect2 4.x writes these
         and they are explicit, so they win over any guess from a name;
      3. the column names, against TUMOR/TUMOUR and NORMAL. The three paired callers in
         tests/verification2 -- Mutect2 3.8, Strelka, SomaticSniper -- all use those literal
         names and none of them declares the header fields above;
      4. nothing matched: raise, naming the columns, rather than pick one. A wrong sample is
         a plausible low allele fraction, which is worse than a refusal.

    Returns (tumour, normal_or_None, how) where `how` names the rule that decided, so the
    run can state it.
    """
    samples = [c for c in vcf_columns if c not in FIXED_VCF_COLUMNS]
    if not samples:
        return None, None, "no sample column"
    if len(samples) == 1:
        return samples[0], None, "the only sample column"

    declared_t = declared_n = None
    for line in meta_lines:
        m = re.match(r"##tumor_sample=(.+)$", line.strip())
        if m:
            declared_t = m.group(1).strip()
        m = re.match(r"##normal_sample=(.+)$", line.strip())
        if m:
            declared_n = m.group(1).strip()
    if declared_t and declared_t in samples:
        other = declared_n if declared_n in samples else None
        if other is None:
            other = next((s for s in samples if s != declared_t), None)
        return declared_t, other, "##tumor_sample in the header"

    upper = {s.upper(): s for s in samples}
    t = next((upper[n] for n in TUMOUR_NAMES if n in upper), None)
    n = next((upper[x] for x in NORMAL_NAMES if x in upper), None)
    if t:
        return t, n or next((s for s in samples if s != t), None), "the column name"

    raise RuntimeError(
        "cannot tell which sample column is the tumour: found "
        + ", ".join(repr(s) for s in samples)
        + ". SMART reports one sample, and every value it derives from a sample column "
          "describes the tumour, so guessing would put the normal's allele fraction in the "
          "report. Name the tumour column with --tumour-column, or give the VCF a "
          "'##tumor_sample=' header."
    )


def load_nm_transcripts(filename):
    """
    Load transcript whitelist from file.

    Reads a list of NM_* transcript identifiers and removes version suffixes
    (e.g. NM_000123.4 → NM_000123) for consistent matching.

    Parameters
    ----------
    filename : str
        Path to transcript whitelist file.

    Returns
    -------
    set
        Set of transcript identifiers without version numbers.
    """
    with open(filename) as f:
        return {line.strip().split(".")[0] for line in f if line.strip()}


# ------------------------------------------------------------
# MAIN: parse VCF → return DataFrame
# ------------------------------------------------------------
def main(args):
    """
    Parse a VCF file and extract structured variant annotations.

    Performs:
      - VCF header parsing to extract CSQ annotation schema
      - Selection of a single transcript per variant using NM whitelist
      - Extraction of VEP annotations (CSQ field)
      - Parsing of INFO fields into key-value pairs
      - Assembly of a structured DataFrame

    Transcript selection priority:
      1. MANE_SELECT / MANE / MANE_PLUS_CLINICAL
      2. HGVSc-based NM transcript
      3. Fallback to first available annotation

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing parsed and structured variant data.
    """

    # --------------------------------------------------------
    # Parse VCF header
    # --------------------------------------------------------
    CSQ_FIELDS = None
    VCF_COLUMNS = None
    META_LINES = []

    with open(args.vcf) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'Format: (.+)">', line)
                CSQ_FIELDS = m.group(1).split("|")

            if line.startswith("##"):
                META_LINES.append(line)
            elif line.startswith("#CHROM"):
                VCF_COLUMNS = line.strip().lstrip("#").split("\t")
                break

    if not CSQ_FIELDS:
        raise RuntimeError("CSQ header not found")

    nm_transcripts = load_nm_transcripts(args.transcripts)
    rows = []

    # --------------------------------------------------------
    # Parse variants
    # --------------------------------------------------------
    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            values = line.rstrip("\n").split("\t")
            record = dict(zip(VCF_COLUMNS, values))

            info_raw = record.pop("INFO")

            # Parse INFO safely
            info = {}
            for part in info_raw.split(";"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    info[k] = v
                else:
                    info[part] = True

            csq_raw = info.pop("CSQ", None)

            # oncokb2.0.py pins the transcript it used via this INFO tag.
            # When present it takes priority over the whitelist selection so
            # both scripts always agree on which CSQ entry the row represents.
            # Strip it from info so it doesn't appear twice in the output row.
            pinned_nm = info.pop("ONCOKB_PREFERRED_TRANSCRIPT", "") or ""
            if pinned_nm == ".":
                pinned_nm = ""

            # Initialize row with all non-INFO VCF columns
            row = dict(record)

            # ------------------------------------------------
            # Handle CSQ
            # ------------------------------------------------
            selected_vep = {f: "" for f in CSQ_FIELDS}
            nm_selected = ""

            if csq_raw:
                first_vep = None
                first_nm = ""
                first_feature = None
                for csq in csq_raw.split(","):
                    fields = csq.split("|")
                    if len(fields) < len(CSQ_FIELDS):
                        fields.extend([""] * (len(CSQ_FIELDS) - len(fields)))

                    vep = dict(zip(CSQ_FIELDS, fields))

                    # Keep the first record as fallback before any transcript match
                    if first_vep is None:
                        first_vep = vep

                    nm = ""
                    for key in ("MANE_SELECT", "MANE", "MANE_PLUS_CLINICAL"):
                        if vep.get(key, "").startswith("NM_"):
                            nm = vep[key].split(".")[0]
                            break

                    # VEP annotates RefSeq transcripts directly in Feature
                    # (e.g. NM_001346897.2) without populating MANE fields.
                    if not nm:
                        feat = vep.get("Feature", "")
                        if feat.startswith("NM_"):
                            nm = feat.split(".")[0]

                    if not nm and vep.get("HGVSc", "").startswith("NM_"):
                        nm = vep["HGVSc"].split(":")[0].split(".")[0]

                    # Record what the first block was, so a fallback can name the transcript
                    # it actually used. first_nm was initialised and never assigned, so a
                    # fallback row left NM_Transcript empty -- which reads as "no transcript"
                    # rather than "an arbitrary one", the more dangerous of the two.
                    if first_feature is None:
                        first_feature = vep.get("Feature", "")
                        first_nm = nm

                    if pinned_nm:
                        # Use the transcript pinned by oncokb2.0.py
                        if nm == pinned_nm:
                            selected_vep = vep
                            nm_selected = nm
                            break
                    elif nm in nm_transcripts:
                        # Standard whitelist-based selection
                        selected_vep = vep
                        nm_selected = nm
                        break

                if not nm_selected and first_vep is not None:
                    # No match found — use the first CSQ record
                    selected_vep = first_vep
                    nm_selected = first_nm

            row["NM_Transcript"] = nm_selected
            # Two columns stating whether the annotation describes the transcript that was
            # asked for. Without them a fallback is invisible: VEP annotates every transcript
            # and this script keeps one, so a row whose gene the whitelist does not cover
            # silently carries whichever block VEP emitted first.
            #
            #   TRANSCRIPT_REQUESTED  the supplied accession that was used, or empty when none
            #                         of the supplied transcripts was annotated here
            #   TRANSCRIPT_MATCH      "match"          the selected transcript is one that was
            #                                          supplied
            #                         "fallback:<acc>" none was annotated for this variant, so
            #                                          VEP's first block was used; <acc> names it
            #
            # The test is membership of the supplied list, not which code path got there. A
            # first version reported the OncoKB-pinned path as its own outcome, "pinned:<acc>",
            # which was wrong: oncokb2.0.py pins a transcript only when
            # `candidate in preferred_tx`, so a pin is always one of the supplied transcripts
            # and therefore a match. It made 32 of 37 rows of the verification cohort look as
            # though the transcript had been imposed from elsewhere, and -- worse -- a real
            # fallback would have been lost among them.
            if nm_selected and nm_selected in nm_transcripts:
                requested, match = nm_selected, "match"
            else:
                used = first_feature or first_nm or "unknown"
                requested, match = "", f"fallback:{used}"
            row["TRANSCRIPT_REQUESTED"] = requested
            row["TRANSCRIPT_MATCH"] = match
            row.update(selected_vep)

            # Add remaining INFO fields
            for k, v in info.items():
                row[k] = v

            rows.append(row)

    # --------------------------------------------------------
    # Build output columns
    # --------------------------------------------------------
    output_columns = (
        [c for c in VCF_COLUMNS if c != "INFO"] +
        ["NM_Transcript"] +
        CSQ_FIELDS +
        sorted({k for r in rows for k in r} -
               set(VCF_COLUMNS) -
               set(CSQ_FIELDS) -
               {"NM_Transcript"})
    )

    # First: build the DataFrame
    df = pd.DataFrame(rows, columns=output_columns)
    # The VCF's meta lines travel with the frame so that final_requiretments can decide which
    # sample column is the tumour from '##tumor_sample=' where a caller declares it. Carried
    # on .attrs rather than by changing main's return type, which tests/test_multi_transcript
    # depends on.
    df.attrs["meta_lines"] = META_LINES
    # The VCF's own columns, which is what identifies a sample column. The frame's columns
    # cannot: by this point it also holds NM_Transcript, every CSQ field and the derived
    # TRANSCRIPT_* and END columns, and anything not in the nine fixed VCF columns would look
    # like a sample. Passing df.columns made even a single-sample VCF unidentifiable.
    df.attrs["vcf_columns"] = VCF_COLUMNS

    # Nothing else here — keep FORMAT and sample column as-is
    return df




# ------------------------------------------------------------
# Additional requirements
# ------------------------------------------------------------
def calculate_end(row):
    """
    Compute genomic end position for a variant.

    Uses variant type (VARIANT_CLASS) and allele lengths to determine the
    correct end coordinate for:
      - SNVs
      - Insertions
      - Deletions
      - Structural variants

    Parameters
    ----------
    row : pandas.Series
        Row containing variant information.

    Returns
    -------
    int or str
        Calculated end position, or empty string if not applicable.
    """
    pos = int(row["POS"])
    ref = str(row["REF"])
    alt = str(row["ALT"])
    variant_class = str(row.get("VARIANT_CLASS", "")).lower()

    if variant_class == "snv":
        return pos

    if variant_class == "chromosome_breakpoint":
        return ""

    if variant_class == "insertion":
        svlen = row.get("SVLEN")
        if pd.notna(svlen):
            try:
                svlen = int(str(svlen).split(",")[0])
                return pos + svlen - 1
            except ValueError:
                pass
        return pos + len(alt) - 1

    if variant_class == "deletion":
        return pos + len(ref) - 1

    # No VARIANT_CLASS, or one not listed above. A VCF record spans len(REF) bases by
    # definition, which is also what post_analysis.calculate_end_position and
    # oncokb_maf_query.split compute -- returning `pos` instead gave a one-base end for
    # every multi-base REF, and the three implementations disagreed with one right.
    if ref and set(ref.upper()) <= set("ACGTN"):
        return pos + len(ref) - 1
    return pos


aa_dict = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

AA3_RE = re.compile("|".join(aa_dict))


def convert_hgvsp_short(hgvsp):
    """
    Convert long HGVSp notation to short protein change format.

    Every three-letter amino-acid code is replaced by its one-letter code, so the
    conversion works for all HGVS protein forms, not only substitutions:

        p.Val600Glu               → p.V600E
        p.Arg248Ter               → p.R248*
        p.Thr15GlufsTer55         → p.T15Efs*55
        p.Met1673_Leu1676del      → p.M1673_L1676del
        p.His24_Ser37delinsPro    → p.H24_S37delinsP
        p.Thr477%3D               → p.T477=

    Lower-case HGVS keywords (del, ins, dup, fs, ext) never match a three-letter
    code, which always starts with a capital letter.

    Parameters
    ----------
    hgvsp : str
        Protein change annotation.

    Returns
    -------
    str
        Shortened protein change representation.
    """
    if not isinstance(hgvsp, str) or "p." not in hgvsp:
        return ""
    protein = hgvsp.split(":")[-1].replace("%3D", "=")
    return AA3_RE.sub(lambda m: aa_dict[m.group(0)], protein)


def final_requiretments(df, output_csv, meta_lines=(), vcf_columns=None,
                        tumour_column=None, normal_column=None):
    """
    Apply final transformations and generate output table.

    Performs:
      - End position calculation
      - Standardization of MAF-like fields
      - Reordering of columns for readability
      - Expansion of OncoKB annotations:
          - ONCOKB_JSON (flattened)
          - ONCOKB_DIAGNOSTIC_IMPLICATIONS (array expansion)
          - ONCOKB_TREATMENTS (array expansion with drug flattening)
      - Cleanup of special characters in text fields

    Writes final output as a tab-separated file.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe from main() function.
    output_csv : str
        Path to output file.

    Returns
    -------
    None
    """
    # A multi-allelic record has no single ALT for the row to be about: nothing splits
    # them, so VAF and every allele count come out empty (post_analysis._counts) and the
    # row's HGVS, consequence and OncoKB answer describe an allele nothing selected. The
    # input this pipeline is built for -- DRAGEN somatic -- writes one record per ALT, so
    # this is normally zero; it is announced rather than assumed, because the failure is
    # silent and a blank allele fraction looks like a caller that reported none.
    # See issue 33.
    multiallelic = df["ALT"].astype(str).str.contains(",", na=False) if "ALT" in df.columns else None
    if multiallelic is not None and multiallelic.any():
        n = int(multiallelic.sum())
        sys.stderr.write(
            f"WARNING: {n} multi-allelic record(s) in this VCF. VAF, t_ref_count, "
            f"t_alt_count and t_depth will be empty for them, and the allele their "
            f"annotation describes is unspecified. Split the VCF first "
            f"(bcftools norm -m-any) to annotate them. See issue 33.\n")
        for _, r in df[multiallelic].head(10).iterrows():
            sys.stderr.write(f"         {r.get('CHROM')}:{r.get('POS')} "
                             f"{r.get('REF')}>{r.get('ALT')}\n")
        if n > 10:
            sys.stderr.write(f"         ... and {n - 10} more\n")

    # INFO END is carried by structural variants only, so in a VCF holding even one such
    # record the column exists and is empty for every SNV and indel. Gating on the
    # column's absence therefore skipped the calculation for the whole file as soon as
    # one Manta record was present.
    #
    # This is not only a reporting defect. MafAnnotator's per-row query resolution falls
    # back to Genomic_Change for a row with no usable protein change, and that fallback
    # tests only that the genomic *headers* exist -- so with an empty End_Position it
    # builds the query from a blank end and OncoKB answers nothing, which is reported as
    # ONCOGENIC=Unknown. Every variant with no protein change (promoter, splice-site,
    # regulatory) was affected. See issue 32.
    if "END" not in df.columns:
        df["END"] = df.apply(calculate_end, axis=1)
    else:
        blank = df["END"].isna() | (df["END"].astype(str).str.strip() == "")
        if blank.any():
            df.loc[blank, "END"] = df.loc[blank].apply(calculate_end, axis=1)

    fixed_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"]
    # Which column is the tumour, decided by select_sample_columns rather than by position.
    # The chosen columns are renamed to canonical names so that nothing downstream has to
    # repeat the decision -- post_analysis expands FORMAT_DATA for the tumour and
    # FORMAT_DATA_NORMAL for the matched normal, and a reader of the table can see which
    # sample each came from in Tumor_Sample_Barcode and Matched_Norm_Sample_Barcode.
    if vcf_columns is None:
        vcf_columns = [c for c in df.columns if c in FIXED_VCF_COLUMNS]
    if tumour_column:
        tumour_col, normal_col, how = tumour_column, normal_column, "--tumour-column"
        if tumour_col not in df.columns:
            raise RuntimeError(f"--tumour-column {tumour_col!r} is not a column of this VCF; "
                               f"it has {[c for c in df.columns if c not in fixed_cols]}")
    else:
        tumour_col, normal_col, how = select_sample_columns(meta_lines, vcf_columns)
    print(f"Tumour sample column: {tumour_col!r} ({how})"
          + (f"; matched normal: {normal_col!r}" if normal_col else "; no matched normal"))

    df["NCBI_Build"] = "GRCh38"
    df["Tumor_Sample_Barcode"] = tumour_col or ""
    df["Matched_Norm_Sample_Barcode"] = normal_col or ""
    renames = {}
    if tumour_col:
        renames[tumour_col] = "FORMAT_DATA"
    if normal_col:
        renames[normal_col] = "FORMAT_DATA_NORMAL"
    if renames:
        df = df.rename(columns=renames)
    sample_col = tumour_col

    if "HGVSp" in df.columns:
        df["HGVSp_Short"] = df["HGVSp"].apply(convert_hgvsp_short)
    else:
        df["HGVSp_Short"] = ""

    df["Hugo_Symbol"] = df["SYMBOL"] if "SYMBOL" in df.columns else ""

    df["Chromosome"] = df["CHROM"]
    df["Start_Position"] = df["POS"]
    df["End_Position"] = df["END"]
    df["Reference_Allele"] = df["REF"]
    df["Tumor_Seq_Allele1"] = df["REF"]
    df["Tumor_Seq_Allele2"] = df["ALT"]

    # --- Variant_Classification from VEP Consequence ---
    consequence_to_classification = {
        "missense_variant":                 "Missense_Mutation",
        "stop_gained":                      "Nonsense_Mutation",
        "stop_lost":                        "Nonstop_Mutation",
        "synonymous_variant":               "Silent",
        "splice_donor_variant":             "Splice_Site",
        "splice_acceptor_variant":          "Splice_Site",
        "splice_region_variant":            "Splice_Region",
        "frameshift_variant":               "Frame_Shift_Del",
        "inframe_insertion":                "In_Frame_Ins",
        "inframe_deletion":                 "In_Frame_Del",
        "start_lost":                       "Translation_Start_Site",
        "5_prime_UTR_variant":              "5'UTR",
        "3_prime_UTR_variant":              "3'UTR",
        "intron_variant":                   "Intron",
        "upstream_gene_variant":            "5'Flank",
        "downstream_gene_variant":          "3'Flank",
        "intergenic_variant":               "IGR",
        "non_coding_transcript_exon_variant": "RNA",
    }

    def get_variant_classification(row):
        consequence = row.get("Consequence", "")
        if not isinstance(consequence, str):
            return ""
        ref = str(row.get("REF", ""))
        alt = str(row.get("ALT", ""))
        for term in consequence.split("&"):
            if term == "frameshift_variant":
                return "Frame_Shift_Ins" if len(alt) > len(ref) else "Frame_Shift_Del"
            if term in consequence_to_classification:
                return consequence_to_classification[term]
        return consequence.split("&")[0]

    def get_variant_type(ref, alt):
        ref = str(ref)
        alt = str(alt)
        if len(ref) == 1 and len(alt) == 1:
            return "SNP"
        elif len(ref) > len(alt):
            return "DEL"
        elif len(ref) < len(alt):
            return "INS"
        return "SNP"

    df["Variant_Classification"] = df.apply(get_variant_classification, axis=1)
    df["Variant_Type"] = df.apply(lambda r: get_variant_type(r["REF"], r["ALT"]), axis=1)

    first_cols = [
        "NCBI_Build", "Hugo_Symbol", "Tumor_Sample_Barcode",
        "HGVSp_Short", "HGVSp", "Chromosome", "Start_Position",
        "End_Position", "Reference_Allele", "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type",
    ]

    first_cols = [c for c in first_cols if c in df.columns]
    remaining = [c for c in df.columns if c not in first_cols]

    df = df[first_cols + remaining]

    # --------------------------------------------------------
    # Expand ONCOKB_JSON into separate columns
    # --------------------------------------------------------
    if "ONCOKB_JSON" in df.columns:
        # Parse JSON safely
        parsed = df["ONCOKB_JSON"].apply(
            lambda x: json.loads(x) if isinstance(x, str) and x.strip().startswith("{") else {}
        )

        # Deep flatten JSON (max_level=None allows full expansion)
        expanded = json_normalize(parsed, sep=".", max_level=None)

        # Prefix to avoid collisions
        expanded.columns = [f"ONCOKB_{c}" for c in expanded.columns]

        # Rebuild dataframe: drop raw column, append expanded
        df_no_json = df.drop(columns=["ONCOKB_JSON"])
        df = pd.concat([df_no_json, expanded], axis=1)

    # --------------------------------------------------------
    # Expand ONCOKB_DIAGNOSTIC_IMPLICATIONS into separate columns
    # --------------------------------------------------------
    if "ONCOKB_DIAGNOSTIC_IMPLICATIONS" in df.columns:
        # Parse JSON array safely
        parsed = df["ONCOKB_DIAGNOSTIC_IMPLICATIONS"].apply(
            lambda x: json.loads(x) if isinstance(x, str) and x.strip().startswith("[") else []
        )

        # Find max number of entries across all rows
        max_entries = parsed.apply(len).max()

        expanded_diag = pd.DataFrame(index=df.index)

        for i in range(max_entries):
            # Extract the i-th entry from each row (or empty dict if not present)
            entry_series = parsed.apply(lambda lst: lst[i] if i < len(lst) else {})

            # Flatten nested dicts
            flat = json_normalize(entry_series.tolist(), sep=".", max_level=None)
            flat.index = df.index
            flat.columns = [f"ONCOKB_DIAG_{i}_{c}" for c in flat.columns]

            expanded_diag = pd.concat([expanded_diag, flat], axis=1)

        # Rebuild: drop raw column, append expanded
        df_no_diag = df.drop(columns=["ONCOKB_DIAGNOSTIC_IMPLICATIONS"])
        df = pd.concat([df_no_diag, expanded_diag], axis=1)

    # --------------------------------------------------------
    # Expand ONCOKB_TREATMENTS into separate columns
    # --------------------------------------------------------
    if "ONCOKB_TREATMENTS" in df.columns:

        # Parse JSON array safely
        parsed = df["ONCOKB_TREATMENTS"].apply(
            lambda x: json.loads(x) if isinstance(x, str) and x.strip().startswith("[") else []
        )

        # Find max number of entries across all rows
        max_entries = parsed.apply(len).max()

        expanded_tx = pd.DataFrame(index=df.index)

        for i in range(max_entries):
            entry_series = parsed.apply(lambda lst: lst[i] if i < len(lst) else {})

            # Flatten drugs array into a single "Drug1 + Drug2" string before normalizing
            def flatten_drugs(entry):
                if not entry:
                    return entry
                entry = dict(entry)
                drugs = entry.get("drugs", [])
                if isinstance(drugs, list):
                    entry["drugs"] = " + ".join(d.get("drugName", "") for d in drugs)
                return entry

            entry_series = entry_series.apply(flatten_drugs)

            flat = json_normalize(entry_series.tolist(), sep=".", max_level=None)
            flat.index = df.index
            flat.columns = [f"ONCOKB_TX_{i}_{c}" for c in flat.columns]

            # Replace tab and newline characters with _ in all string columns
            for col in flat.columns:
                if flat[col].dtype == object:
                    flat[col] = flat[col].apply(
                        lambda x: re.sub(r"[\t\n\r]", "_", x) if isinstance(x, str) else x
                    )

            expanded_tx = pd.concat([expanded_tx, flat], axis=1)

        # Rebuild: drop raw column, append expanded
        df_no_tx = df.drop(columns=["ONCOKB_TREATMENTS"])
        df = pd.concat([df_no_tx, expanded_tx], axis=1)

    df.to_csv(output_csv, index=False, sep="\t")

# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    args = get_args()
    df = main(args)
    final_requiretments(df, args.output, meta_lines=df.attrs.get("meta_lines", ()),
                        vcf_columns=df.attrs.get("vcf_columns"),
                        tumour_column=getattr(args, "tumour_column", None),
                        normal_column=getattr(args, "normal_column", None))
