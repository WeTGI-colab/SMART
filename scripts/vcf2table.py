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
    parser.add_argument("-o", "--output", required=True, help="Final output file")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


# ------------------------------------------------------------
# Load NM transcript whitelist
# ------------------------------------------------------------
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

    with open(args.vcf) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'Format: (.+)">', line)
                CSQ_FIELDS = m.group(1).split("|")

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

    return pos


aa_dict = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

def convert_hgvsp_short(hgvsp):
    """
    Convert long HGVSp notation to short protein change format.

    Example:
        p.Val600Glu → p.V600E

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
    match = re.search(r"p\.([A-Za-z]+)(\d+)([A-Za-z\*]+)", hgvsp)
    if match:
        ref, pos, alt = match.groups()
        return f"p.{aa_dict.get(ref, ref)}{pos}{aa_dict.get(alt, alt)}"
    return hgvsp.split(":")[-1]


def final_requiretments(df, output_csv):
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
    if "END" not in df.columns:
        df["END"] = df.apply(calculate_end, axis=1)

    fixed_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT"]
    sample_candidates = [c for c in df.columns if c not in fixed_cols]
    sample_col = sample_candidates[0] if sample_candidates else None

    df["NCBI_Build"] = "GRCh38"
    df["Tumor_Sample_Barcode"] = sample_col if sample_col else ""

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
    final_requiretments(df, args.output)
