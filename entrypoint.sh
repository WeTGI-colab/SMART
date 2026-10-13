#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# SMART — Somatic Mutation Annotation and Reporting Tool
# =============================================================================

usage() {
cat <<EOF
SMART — Somatic Mutation Annotation and Reporting Tool (Dockerised)

Usage:
  docker run -v /path/to/data:/data -v /path/to/refs:/refs smart \\
    <ONCOKB_TOKEN> --transcripts-file FILE --ref-dir /refs --config FILE [OPTIONS]

Required:
  --ref-dir DIR            Reference resources directory (mounted at /refs)
  --output-dir DIR         Directory where all results will be written. Must be a
                           mounted volume (e.g. -v /host/results:/results) so that
                           output is preserved after the container exits.

Optional (full pipeline):
  <ONCOKB_TOKEN>           OncoKB API token. If omitted, the pipeline stops
                           automatically after VEP annotation (VEP-only mode).
  --transcripts-file FILE  Transcript list for prioritisation (required when token given)
  --spliceai compute|lookup  Compute SpliceAI against the reported transcripts (default:
                           compute) or look it up in the precomputed tables. The tables
                           carry no transcript and use a +/-50 bp window; see issue #2
  --spliceai-mask 0|1      Raw (0) or masked (1) scores. Default 1, which is what the
                           model's developers recommend for variant interpretation
  --config FILE            YAML config file for post_analysis (required when token given)

  The ref-dir should contain:
    DIR/liftover/hg19ToHg38.over.chain
    DIR/liftover/hg38.fa
    DIR/Plugins
    DIR/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz
    DIR/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz
    DIR/REVEL/new_tabbed_revel_grch38.tsv.gz
    DIR/ClinVar/clinvar.vcf.gz
    DIR/CIVIC/civic_grch38.vcf.gz
    DIR/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz
    DIR/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz
    Optional: DIR/COSMIC/cosmic_cmc_grch38.vcf.gz, DIR/CADD/*, DIR/dbNSFP/dbNSFP5.0a_grch38.smart.gz

Options:
  --pass / --no-pass               Enable/disable PASS filtering (default: ON)
  --ref-check fail|warn            Stop, or only warn, when a VCF's REF alleles do
                                   not match the reference genome (default: fail)
  --liftover / --no-liftover       Enable/disable liftover (default: ON)
  --clean-tmp / --keep-tmp         Delete/keep intermediate files (default: clean)
  --clean-tables / --keep-tables   Delete/keep table files after post_analysis (default: clean)
  --jobs N                         Number of samples to process in parallel on this
                                   machine (default: 1 = sequential). Each job is one
                                   sample; all jobs share the same CPU and RAM.
  --tumour-type TYPE               Tumour type for the Tier 4 SVIG-UK classification: an
                                   OncoTree code (e.g. AML, MDS, LUAD) or myeloid | lymphoid |
                                   haem | solid | unknown (default: unknown). A per-variant
                                   TUMOUR_TYPE column in the input takes precedence.
  --o4-on-target LEVEL             GENIE entries counted for O4: code | tissue | lineage | all
                                   (default: lineage)
  --no-classification              Do not write the Tier 4 classification
  --vep-only                       Stop after VEP annotation. Skips OncoKB,
                                   vcf2table, MafAnnotator, and post-analysis.
                                   Annotated VCFs are written to AnnotatedVcf/.
                                   Also triggered automatically when no token is given.
  --help                           Show this help

Volumes:
  Mount your data directory (containing OriginalVcf/) at /data
  Mount your reference directory at /refs (or wherever --ref-dir points)

Examples:
  docker run --rm \\
    -v \$(pwd):/data \\
    -v /mnt/data1/vep_test:/refs \\
    smart \\
    "\$ONCOKB_TOKEN" \\
    --transcripts-file /data/TSO500_transcripts_list.txt \\
    --ref-dir /refs \\
    --config /data/Config.yaml

EOF
}

if [[ $# -lt 1 ]] || [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    usage
    exit 0
fi

# Defaults
DO_PASS_FILTER=1
REF_CHECK_MODE=fail
DO_LIFTOVER=1
CLEAN_TMP=1
CLEAN_TABLES=1
JOBS=1
VEP_ONLY=0
TUMOUR_TYPE="unknown"
O4_ON_TARGET="lineage"
DO_CLASSIFICATION=1

# The OncoKB token is the first positional argument.
# If it is absent (first arg is a flag) or an empty string, automatically
# activate VEP-only mode — the pipeline stops after VEP annotation and
# skips OncoKB, vcf2table, MafAnnotator, and post-analysis.
if [[ -z "$1" || "$1" == --* ]]; then
    ONCOKB_TOKEN=""
    VEP_ONLY=1
    echo "INFO: No OncoKB token provided — running in VEP-only mode."
    echo "      To run the full pipeline supply a token as the first argument."
    echo ""
else
    ONCOKB_TOKEN="$1"
    shift
fi
TRANSCRIPTS_FILE=""
# SpliceAI: compute the scores against the reported transcripts (default) or look them up
# in Illumina's precomputed tables. The tables carry no transcript, were generated on
# GENCODE v24 and within a +/-50 bp window; see issue #2. Masked is what the model's
# developers recommend for variant interpretation, and what this laboratory uses.
SPLICEAI_MODE="compute"
SPLICEAI_MASK="1"
REF_DIR=""
CONFIG_FILE=""
OUTPUT_DIR_BASE=""

# Parse flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-pass)
            DO_PASS_FILTER=0
            shift
            ;;
        --pass)
            DO_PASS_FILTER=1
            shift
            ;;
        --ref-check)
            REF_CHECK_MODE="$2"
            if [[ "$REF_CHECK_MODE" != "fail" && "$REF_CHECK_MODE" != "warn" ]]; then
                echo "ERROR: --ref-check takes 'fail' or 'warn', got '$REF_CHECK_MODE'"
                exit 1
            fi
            shift 2
            ;;
        --no-liftover)
            DO_LIFTOVER=0
            shift
            ;;
        --liftover)
            DO_LIFTOVER=1
            shift
            ;;
        --clean-tmp)
            CLEAN_TMP=1
            shift
            ;;
        --keep-tmp)
            CLEAN_TMP=0
            shift
            ;;
        --clean-tables)
            CLEAN_TABLES=1
            shift
            ;;
        --keep-tables)
            CLEAN_TABLES=0
            shift
            ;;
        --jobs)
            [[ $# -lt 2 ]] && { echo "ERROR: --jobs requires a number"; exit 2; }
            JOBS="$2"
            [[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --jobs must be a positive integer"; exit 2; }
            shift 2
            ;;
        --vep-only)
            VEP_ONLY=1
            shift
            ;;
        --tumour-type)
            [[ $# -lt 2 ]] && { echo "ERROR: --tumour-type requires a value"; exit 2; }
            TUMOUR_TYPE="$2"
            shift 2
            ;;
        --o4-on-target)
            [[ $# -lt 2 ]] && { echo "ERROR: --o4-on-target requires a value"; exit 2; }
            O4_ON_TARGET="$2"
            [[ "$O4_ON_TARGET" =~ ^(code|tissue|lineage|all)$ ]] || { echo "ERROR: --o4-on-target must be code, tissue, lineage or all"; exit 2; }
            shift 2
            ;;
        --no-classification)
            DO_CLASSIFICATION=0
            shift
            ;;
        --input-dir)
            [[ $# -lt 2 ]] && { echo "ERROR: --input-dir requires a path"; exit 2; }
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            [[ $# -lt 2 ]] && { echo "ERROR: --output-dir requires a path"; exit 2; }
            OUTPUT_DIR_BASE="$2"
            shift 2
            ;;
        --spliceai)
            [[ $# -lt 2 ]] && { echo "ERROR: --spliceai requires compute|lookup"; exit 2; }
            SPLICEAI_MODE="$2"
            [[ "$SPLICEAI_MODE" != "compute" && "$SPLICEAI_MODE" != "lookup" ]] && {
                echo "ERROR: --spliceai takes 'compute' or 'lookup', got '$SPLICEAI_MODE'"; exit 2; }
            shift 2
            ;;
        --spliceai-mask)
            [[ $# -lt 2 ]] && { echo "ERROR: --spliceai-mask requires 0|1"; exit 2; }
            SPLICEAI_MASK="$2"
            [[ "$SPLICEAI_MASK" != "0" && "$SPLICEAI_MASK" != "1" ]] && {
                echo "ERROR: --spliceai-mask takes '0' (raw) or '1' (masked), got '$SPLICEAI_MASK'"; exit 2; }
            shift 2
            ;;
        --transcripts-file)
            [[ $# -lt 2 ]] && { echo "ERROR: --transcripts-file requires a path"; exit 2; }
            TRANSCRIPTS_FILE="$2"
            shift 2
            ;;
        --ref-dir)
            [[ $# -lt 2 ]] && { echo "ERROR: --ref-dir requires a path"; exit 2; }
            REF_DIR="$2"
            shift 2
            ;;
        --config)
            [[ $# -lt 2 ]] && { echo "ERROR: --config requires a path"; exit 2; }
            CONFIG_FILE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 2
            ;;
    esac
done

# Validate required arguments
[[ -z "$OUTPUT_DIR_BASE" ]] && {
    echo "ERROR: --output-dir is required."
    echo "       Mount a host directory as a volume and pass it here, e.g.:"
    echo "         -v /host/results:/results  --output-dir /results"
    echo "       Without a mounted volume the results will be lost when the container exits."
    usage; exit 2;
}
[[ -z "$REF_DIR" ]] && { echo "ERROR: --ref-dir is required"; usage; exit 2; }
[[ ! -d "$REF_DIR" ]] && { echo "ERROR: Reference directory not found: $REF_DIR"; exit 2; }
if [[ $VEP_ONLY -eq 0 ]]; then
    [[ -z "$TRANSCRIPTS_FILE" ]] && { echo "ERROR: --transcripts-file is required (or use --vep-only)"; usage; exit 2; }
    [[ ! -f "$TRANSCRIPTS_FILE" ]] && { echo "ERROR: Transcript file not found: $TRANSCRIPTS_FILE"; exit 2; }
    [[ -z "$CONFIG_FILE" ]] && { echo "ERROR: --config is required (or use --vep-only)"; usage; exit 2; }
    [[ ! -f "$CONFIG_FILE" ]] && { echo "ERROR: Config file not found: $CONFIG_FILE"; exit 2; }
fi
INPUT_DIR="${INPUT_DIR:-./OriginalVcf}"
[[ ! -d "$INPUT_DIR" ]] && { echo "ERROR: Input directory not found: $INPUT_DIR"; exit 2; }

# SMART version — read from VERSION file bundled in the image
SMART_VERSION="$(cat "$(dirname "${BASH_SOURCE[0]}")/VERSION" 2>/dev/null || echo "unknown")"

# Reference files
# You may have problems with the filename of ref files
# as org. change them frecuenly (e.g., new versions)
CHAIN="$REF_DIR/liftover/hg19ToHg38.over.chain"
REF="$REF_DIR/liftover/hg38.fa"
VEP_DIR="$REF_DIR"
VEP_PLUGIN_DIR="$REF_DIR/Plugins"
SPLICEAI_SNV="$REF_DIR/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$REF_DIR/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
REVEL_FILE="$REF_DIR/REVEL/new_tabbed_revel_grch38.tsv.gz"
CLINVAR_VCF="$REF_DIR/ClinVar/clinvar.vcf.gz"
CIVIC_VCF="$REF_DIR/CIVIC/civic_grch38.vcf.gz"
LOEUF_FILE="$REF_DIR/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
CANCER_HOTSPOTS_VCF="$REF_DIR/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"
COSMIC_VCF="$REF_DIR/COSMIC/cosmic_cmc_grch38.vcf.gz"   # optional, licence-gated (build with utils/cosmic_cmc_to_vcf.py)
CADD_SNV="$REF_DIR/CADD/whole_genome_SNVs.tsv.gz"         # optional, very large (~80GB) — install with INSTALL_CADD=1 utils/get_ref_files.sh
CADD_INDEL="$REF_DIR/CADD/gnomad.genomes.r4.0.indel.tsv.gz"  # optional
DBNSFP_FILE="$REF_DIR/dbNSFP/dbNSFP5.0a_grch38.smart.gz"   # optional, academic licence (build with utils/build_dbnsfp.py)
GENIE_COUNTS="$REF_DIR/GENIE/genie_o4_counts.tsv.gz"       # optional, GENIE data-use agreement (build with utils/build_genie_o4.py)

# Validate reference files
REQUIRED_PATHS=(
    "$CHAIN" "$REF" "$VEP_PLUGIN_DIR"
    "$SPLICEAI_SNV" "$SPLICEAI_INDEL" "$REVEL_FILE"
    "$CLINVAR_VCF" "$CIVIC_VCF" "$LOEUF_FILE" "$CANCER_HOTSPOTS_VCF"
)
# Next loop checks all neccesary files such as ClinVar VCF
for path in "${REQUIRED_PATHS[@]}"; do
    [[ ! -e "$path" ]] && { echo "ERROR: Required resource not found: $path"; exit 2; }
done

# COSMIC Cancer Mutation Census is OPTIONAL — it requires a COSMIC licence so it is
# not shipped. If the user has built and mounted the VCF (utils/cosmic_cmc_to_vcf.py),
# add it as a VEP --custom annotation; otherwise annotation runs unchanged without it.
COSMIC_CUSTOM=""
if [[ -e "$COSMIC_VCF" ]]; then
    COSMIC_CUSTOM="--custom $COSMIC_VCF,COSMIC,vcf,exact,0,GENE,CNT,TESTED,TIER,ONC_TSG,CGC_TIER,DNDS"
fi

# CADD is OPTIONAL — the GRCh38 score files are very large (~80 GB for whole-genome
# SNVs) so they are not shipped. If both the SNV and indel files are mounted, add the
# CADD plugin (adds CADD_PHRED/CADD_RAW); otherwise annotation runs unchanged without it.
CADD_PLUGIN=""
if [[ -e "$CADD_SNV" && -e "$CADD_INDEL" ]]; then
    CADD_PLUGIN="--plugin CADD,$CADD_SNV,$CADD_INDEL"
fi

# dbNSFP is OPTIONAL — academic releases (CC BY-NC-ND 4.0) must be downloaded by each
# user and cannot be shipped. If the file built with utils/build_dbnsfp.py is mounted, add
# the dbNSFP plugin (in silico predictors not otherwise annotated, e.g. AlphaMissense,
# MetaRNN, BayesDel); otherwise annotation runs unchanged without it.
DBNSFP_PLUGIN=""
if [[ -e "$DBNSFP_FILE" && -e "${DBNSFP_FILE}.tbi" ]]; then
    DBNSFP_PLUGIN="--plugin dbNSFP,$DBNSFP_FILE,ALL"
fi

# Fetch OncoKB version info only when a token is available
if [[ $VEP_ONLY -eq 0 ]]; then
    ONCOKB_INFO=$(curl -s https://www.oncokb.org/api/v1/info || echo "")
    ONCOKB_DATA_VERSION=$(echo "$ONCOKB_INFO" | grep -oP '"dataVersion":\{"version":"\K[^"]+' || echo "unknown")
    ONCOKB_DATA_DATE=$(echo "$ONCOKB_INFO" | grep -oP '"dataVersion":\{"version":"[^"]+","date":"\K[^"]+' || echo "unknown")
    ONCOKB_API_VERSION=$(echo "$ONCOKB_INFO" | grep -oP '"apiVersion":\{"version":"\K[^"]+' || echo "unknown")
else
    ONCOKB_DATA_VERSION="N/A (VEP-only mode)"
    ONCOKB_DATA_DATE=""
    ONCOKB_API_VERSION="N/A (VEP-only mode)"
fi

echo "============================================================"
echo "SMART — Somatic Mutation Annotation and Reporting Tool  v${SMART_VERSION}"
echo "============================================================"
echo "Tool versions:"
echo "  VEP:               $(vep --help 2>&1 | grep -oP 'ensembl-vep\s*:\s*\K\S+' || echo 'unknown')"
echo "  GATK:              $(gatk --version 2>&1 | grep -oP '\d+\.\d+\.\d+\.\d+' | head -1 || echo 'unknown')"
echo "  bcftools:          $(bcftools --version | head -1 || echo 'unknown')"
echo "  Python:            $(python3 --version 2>&1 || echo 'unknown')"
echo "  OncoKB API:        $ONCOKB_API_VERSION"
echo ""

echo "Reference versions:"
echo "  ClinVar:           $(basename "$CLINVAR_VCF")"
echo "  CIViC:             $(basename "$CIVIC_VCF")"
echo "  CancerHotSpots:    $(basename "$CANCER_HOTSPOTS_VCF")"
echo "  COSMIC CMC:        $([[ -e "$COSMIC_VCF" ]] && basename "$COSMIC_VCF" || echo '(not provided — optional, licence-gated)')"
echo "  SpliceAI SNV:      $(basename "$SPLICEAI_SNV")"
echo "  SpliceAI INDEL:    $(basename "$SPLICEAI_INDEL")"
echo "  REVEL:             $(basename "$REVEL_FILE")"
echo "  CADD:              $([[ -e "$CADD_SNV" && -e "$CADD_INDEL" ]] && basename "$CADD_SNV" || echo '(not provided — optional, ~80GB)')"
echo "  dbNSFP:            $([[ -e "$DBNSFP_FILE" ]] && basename "$DBNSFP_FILE" || echo '(not provided — optional, academic licence)')"
echo "  GENIE counts:      $([[ -e "$GENIE_COUNTS" ]] && basename "$GENIE_COUNTS" || echo '(not provided — optional, GENIE data-use agreement)')"
echo "  LOEUF:             $(basename "$LOEUF_FILE")"
echo "  OncoKB data:       $ONCOKB_DATA_VERSION ($ONCOKB_DATA_DATE)"

echo ""

echo ""
echo "Config:"
echo "  PASS filter:      $([[ $DO_PASS_FILTER -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Reference check:  $REF_CHECK_MODE"
echo "  SpliceAI:         $SPLICEAI_MODE, mask=$SPLICEAI_MASK"
echo "  Liftover:         $([[ $DO_LIFTOVER -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  VEP only:         $([[ $VEP_ONLY -eq 1 ]] && echo "YES — stopping after VEP $( [[ -z "$ONCOKB_TOKEN" ]] && echo '(no token provided)' || echo '(--vep-only flag)' )" || echo NO)"
echo "  Clean tmp:        $([[ $CLEAN_TMP -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Clean tables:     $([[ $CLEAN_TABLES -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Parallel jobs:    $JOBS"
echo "  Transcript file:  ${TRANSCRIPTS_FILE:-N/A (VEP-only mode)}"
echo "  Config file:      ${CONFIG_FILE:-N/A (VEP-only mode)}"
echo "  Reference dir:    $REF_DIR"
echo "  Input dir:        $INPUT_DIR"
echo "============================================================"
echo

SCRIPT_DIR="/opt/smart/scripts"

# Create output directories
INPUT_DIR="${INPUT_DIR:-./OriginalVcf}"
FILTERED_DIR="${OUTPUT_DIR_BASE}/FilteredVcf"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/LiftOverVcf"
REJECT_DIR="${OUTPUT_DIR_BASE}/LiftOverVcf/Rejected"
ANNOTATED_DIR="${OUTPUT_DIR_BASE}/AnnotatedVcf"
ONCOKB_DIR="${OUTPUT_DIR_BASE}/OncoKB_VCF"
TABLE_DIR="${OUTPUT_DIR_BASE}/Table"
FINAL_TABLE_DIR="${OUTPUT_DIR_BASE}/FINAL_Table"
LOGS_DIR="${OUTPUT_DIR_BASE}/logs"

# We create directories where all tmp files are saved
mkdir -p "$OUTPUT_DIR_BASE" "$FILTERED_DIR" "$OUTPUT_DIR" "$REJECT_DIR" "$ANNOTATED_DIR" "$ONCOKB_DIR" "$TABLE_DIR" "$FINAL_TABLE_DIR" "$LOGS_DIR"

echo "Cleaning previous results..."
rm -f "$OUTPUT_DIR"/*gz "$OUTPUT_DIR"/*tbi "$OUTPUT_DIR"/*.vcf 2>/dev/null || true
rm -f "$REJECT_DIR"/* 2>/dev/null || true
rm -f "$ANNOTATED_DIR"/* 2>/dev/null || true
rm -f "$FILTERED_DIR"/* 2>/dev/null || true
rm -f "$ONCOKB_DIR"/* 2>/dev/null || true
rm -f "$TABLE_DIR"/* 2>/dev/null || true
rm -f "$FINAL_TABLE_DIR"/* 2>/dev/null || true
echo "Cleanup complete."

COUNT_FILE="${OUTPUT_DIR_BASE}/variant_counts.txt"

# --- SpliceAI annotation (compute mode only) ---------------------------------
# One annotation per run, built from the same transcript list used for prioritisation, so a
# computed score describes the transcript the report names. Rows are named GENE--REFSEQ and
# SpliceAI reports its annotation's NAME as the SYMBOL field, which is how a score becomes
# traceable to a transcript -- the thing the precomputed tables never recorded.
SPLICEAI_ANNOTATION=""
if [[ "$SPLICEAI_MODE" == "compute" && -n "$TRANSCRIPTS_FILE" ]]; then
    GENCODE_TABLE="$REF_DIR/SpliceAI-lookup/references_files/annotations/gencode.v43.annotation.txt.gz"
    if [[ ! -f "$GENCODE_TABLE" ]]; then
        echo "ERROR: --spliceai compute needs $GENCODE_TABLE"
        echo "       Either provide it or run with --spliceai lookup."
        exit 1
    fi
    SPLICEAI_ANNOTATION="${OUTPUT_DIR_BASE}/spliceai_annotation.txt"
    echo ">>> Building the SpliceAI annotation for the reported transcripts"
    python3 "$SCRIPT_DIR/build_spliceai_annotation.py" \
        --annotation "$GENCODE_TABLE" \
        --transcripts "$TRANSCRIPTS_FILE" \
        --out "$SPLICEAI_ANNOTATION" || {
            echo "ERROR: could not build the SpliceAI annotation"; exit 1; }
fi
echo -e "Sample\tOriginal_Variants\tPASS_Variants\tLifted_Variants\tVEP_Annotated_Variants\tOncoKB_Annotated_Variants\tVariantsInTable" > "$COUNT_FILE"

echo "------------------------------------------------------------"
echo "Input VCF directory: $INPUT_DIR"

# Count VCF files
VCF_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.vcf.gz" | wc -l)
echo "Number of VCF files found: $VCF_COUNT"

# Optional: list them (useful for debugging)
if [[ $VCF_COUNT -gt 0 ]]; then
    echo "VCF files:"
    find "$INPUT_DIR" -maxdepth 1 -name "*.vcf.gz"
else
    echo "WARNING: No VCF files found in $INPUT_DIR"
fi

echo "------------------------------------------------------------"

# Optional hard stop if no files
if [[ $VCF_COUNT -eq 0 ]]; then
    echo "ERROR: No input VCFs detected. Exiting."
    exit 1
fi

# ==============================================================================
# process_sample VCF LOG_FILE
#   Runs the full per-sample pipeline in a subshell (safe for backgrounding).
#   Writes its variant-count row to COUNT_FILE.<sample> so parallel runs do
#   not race on the shared COUNT_FILE; the main process merges them afterwards.
#   All stdout/stderr is appended to LOG_FILE.
# ==============================================================================
process_sample() {
    local vcf="$1"
    local log_file="$2"
    local basename_
    basename_=$(basename "$vcf")
    local sample="${basename_%.vcf.gz}"

    {
    echo "###################### Processing: $sample ##########################################"

    # --- PASS filtering ---
    local FILTERED_VCF="$FILTERED_DIR/${sample}.PASS_only.vcf.gz"
    local FILTERED_VCF_TBI="${FILTERED_VCF}.tbi"
    local INPUT_FOR_NEXT

    if [[ $DO_PASS_FILTER -eq 1 ]]; then
        echo ">>> PASS filtering: $sample"
        zcat "$vcf" \
            | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $7 == "PASS"' \
            | bgzip > "$FILTERED_VCF"
        tabix -p vcf "$FILTERED_VCF" || { echo "ERROR: tabix indexing failed for $FILTERED_VCF"; exit 1; }
        INPUT_FOR_NEXT="$FILTERED_VCF"
    else
        echo ">>> PASS filtering DISABLED: $sample"
        INPUT_FOR_NEXT="$vcf"
    fi

    # --- LiftOver ---
    local LIFTED_VCF_GZ="$OUTPUT_DIR/${sample}.hg19tohg38_liftover.vcf.gz"
    local LIFTED_VCF="$OUTPUT_DIR/${sample}.hg19tohg38_liftover.vcf"
    local LIFTED_VCF_TBI="${LIFTED_VCF_GZ}.tbi"
    local REJECTED_VCF="$REJECT_DIR/${sample}.hg19tohg38_rejected.vcf.gz"
    local REJECTED_VCF_TBI="${REJECTED_VCF}.tbi"
    local VEP_INPUT

    if [[ $DO_LIFTOVER -eq 1 ]]; then
        echo ">>> LiftOver: $sample"
        gatk LiftoverVcf \
            -C "$CHAIN" \
            -I "$INPUT_FOR_NEXT" \
            -O "$LIFTED_VCF_GZ" \
            -R "$REF" \
            --WRITE_ORIGINAL_POSITION \
            --REJECT "$REJECTED_VCF"
        gunzip -c "$LIFTED_VCF_GZ" > "$LIFTED_VCF"
        VEP_INPUT="$LIFTED_VCF"
        echo ">>> LiftOver complete: $sample"
    else
        echo ">>> LiftOver DISABLED: $sample"
        VEP_INPUT="$INPUT_FOR_NEXT"
    fi

    # --- SpliceAI (compute mode) ---
    # Run on $VEP_INPUT, so the scores are for the file VEP annotates, in the assembly the
    # reference describes. The released model scores single-base substitutions and indels
    # whose REF matches the genome; it silently returns nothing for multi-nucleotide
    # variants, so those are reported rather than left to look like an absent prediction.
    local SPLICEAI_VCF=""
    if [[ "$SPLICEAI_MODE" == "compute" && -n "$SPLICEAI_ANNOTATION" ]]; then
        echo ">>> SpliceAI (compute, mask=$SPLICEAI_MASK): $sample"
        SPLICEAI_VCF="$OUTPUT_DIR/${sample}.spliceai.vcf"
        if ! spliceai -I "$VEP_INPUT" -O "$SPLICEAI_VCF" -R "$REF" \
                -A "$SPLICEAI_ANNOTATION" -D 500 -M "$SPLICEAI_MASK"; then
            echo "ERROR: SpliceAI failed for $sample. Re-run with '--spliceai lookup' to use"
            echo "       the precomputed tables instead."
            exit 1
        fi
        local N_MNV
        N_MNV=$(awk 'BEGIN{FS="\t"} !/^#/ && length($4)>1 && length($4)==length($5) {n++} END{print n+0}' "$VEP_INPUT")
        [[ "$N_MNV" -gt 0 ]] && echo "WARNING: $N_MNV multi-nucleotide variant(s) in $sample receive no SpliceAI score (the released model does not support them)"
    fi

    # --- Reference check ---
    # On $VEP_INPUT rather than the original VCF: after liftover this is the file in the
    # assembly $REF describes, and it is exactly what VEP will read. VEP does not reject a
    # VCF whose REF disagrees -- it takes the reference base and keeps the submitted ALT --
    # so an unchecked mismatch becomes a consequence for a variant nobody submitted.
    echo ">>> Reference check: $sample"
    if ! python3 "$SCRIPT_DIR/check_reference.py" \
            --vcf "$VEP_INPUT" --reference "$REF" --mode "$REF_CHECK_MODE"; then
        echo "ERROR: $sample has REF alleles that do not match the reference genome."
        echo "       Re-run with '--ref-check warn' to annotate anyway."
        exit 1
    fi

    # In compute mode the plugin is left out entirely: the scores come from our own run and
    # are merged into the table after vcf2table, so the column names are unchanged either way.
    local SPLICEAI_PLUGIN=""
    if [[ "$SPLICEAI_MODE" == "lookup" ]]; then
        SPLICEAI_PLUGIN="--plugin SpliceAI,snv=$SPLICEAI_SNV,indel=$SPLICEAI_INDEL,cutoff=0.5"
    fi

    # --- VEP Annotation ---
    echo ">>> VEP annotation: $sample"
    local ANNOTATED_VCF="$ANNOTATED_DIR/${sample}_annotated.vcf"
    vep \
        -i "$VEP_INPUT" \
        -o "$ANNOTATED_VCF" \
        --vcf --everything --offline --merged \
        --dir "$VEP_DIR" \
        --dir_plugins "$VEP_PLUGIN_DIR" \
        --assembly GRCh38 \
        --no_escape \
        --fasta "$REF" \
        $SPLICEAI_PLUGIN \
        --plugin REVEL,file="$REVEL_FILE",no_match=1 \
        --custom "$CLINVAR_VCF",ClinVar,vcf,exact,0,AF_ESP,AF_EXAC,AF_TGP,ALLELEID,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNSIGSCV,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ONCDN,ONCDNINCL,ONCDISDB,ONCDISDBINCL,ONC,ONCINCL,ONCREVSTAT,ONCSCV,ONCCONF,ORIGIN,RS,SCIDN,SCIDNINCL,SCIDISDB,SCIDISDBINCL,SCIREVSTAT,SCI,SCIINCL,SCISCV \
        --custom "$CIVIC_VCF",CIViC,vcf,exact,0,GN,VT,CSQ \
        --plugin LOEUF,file="$LOEUF_FILE",match_by=transcript \
        --custom "$CANCER_HOTSPOTS_VCF",CancerHotspots,vcf,exact,0,HOTSPOT,HOTSPOT_GENE,HOTSPOT_HGVSp,HOTSPOT3D,HOTSPOT3D_GENE,HOTSPOT3D_HGVSp,HOTSPOTNC,HOTSPOTNC_GENE,HOTSPOTNC_HGVSc \
        $COSMIC_CUSTOM \
        $CADD_PLUGIN \
        $DBNSFP_PLUGIN

    # --- Early exit when --vep-only is set ---
    if [[ $VEP_ONLY -eq 1 ]]; then
        echo ">>> VEP-only mode: skipping OncoKB, vcf2table, MafAnnotator for $sample"
        local ORIGINAL_COUNT PASS_COUNT LIFTED_COUNT ANNOTATED_COUNT
        ORIGINAL_COUNT=$(zgrep -vc "^#" "$vcf")
        if [[ $DO_PASS_FILTER -eq 1 ]]; then
            PASS_COUNT=$(zgrep -vc "^#" "$FILTERED_VCF")
        else
            PASS_COUNT="$ORIGINAL_COUNT"
        fi
        if [[ $DO_LIFTOVER -eq 1 ]]; then
            LIFTED_COUNT=$(grep -vc "^#" "$LIFTED_VCF")
        else
            LIFTED_COUNT="$PASS_COUNT"
        fi
        ANNOTATED_COUNT=$(grep -vc "^#" "$ANNOTATED_VCF")
        echo -e "${sample}\t${ORIGINAL_COUNT}\t${PASS_COUNT}\t${LIFTED_COUNT}\t${ANNOTATED_COUNT}\tN/A\tN/A" \
            > "${COUNT_FILE}.${sample}"
        echo "--- $sample complete (VEP only) ---"
        return 0
    fi

    # --- OncoKB annotation ---
    echo ">>> OncoKB annotation: $sample"
    local ONCOKB_VCF="$ONCOKB_DIR/${sample}.oncoKB.vcf"
    python3 "$SCRIPT_DIR/oncokb2.0.py" "$ANNOTATED_VCF" "$ONCOKB_VCF" \
        --oncokb_token "$ONCOKB_TOKEN" \
        --tumor_mode generic \
        --preferred-transcripts "$TRANSCRIPTS_FILE"

    # --- VCF to Table ---
    echo ">>> VCF2TABLE: $sample"
    local TABLE_CSV="$TABLE_DIR/${sample}.oncoKB.csv"
    python3 "$SCRIPT_DIR/vcf2table.py" \
        --vcf "$ONCOKB_VCF" \
        --transcripts "$TRANSCRIPTS_FILE" \
        --output "$TABLE_CSV" \
        --debug

    # --- MafAnnotator ---
    echo ">>> MafAnnotator: $sample"
    local ANNOTATED_MAF="$FINAL_TABLE_DIR/${sample}.oncoKB.maf"
    # Variants whose selected transcript is not OncoKB's curated transcript are queried by
    # genomic change, so different protein numbering does not hide known variants.
    local MAF_BY_PROTEIN="$TABLE_DIR/${sample}.oncokb_by_protein.tsv"
    local MAF_BY_GENOMIC="$TABLE_DIR/${sample}.oncokb_by_genomic.tsv"
    python3 "$SCRIPT_DIR/oncokb_maf_query.py" split \
        --table "$TABLE_CSV" --token "$ONCOKB_TOKEN" \
        --protein "$MAF_BY_PROTEIN" --genomic "$MAF_BY_GENOMIC"
    python3 /opt/oncokb-annotator/MafAnnotator.py \
        -i "$MAF_BY_PROTEIN" \
        -o "${MAF_BY_PROTEIN%.tsv}.maf" \
        -b "$ONCOKB_TOKEN"
    local GENOMIC_FAILED=""
    if [[ $(tail -n +2 "$MAF_BY_GENOMIC" | wc -l) -gt 0 ]]; then
        if ! python3 /opt/oncokb-annotator/MafAnnotator.py \
            -i "$MAF_BY_GENOMIC" \
            -o "${MAF_BY_GENOMIC%.tsv}.maf" \
            -q Genomic_Change -r GRCh38 \
            -b "$ONCOKB_TOKEN"; then
            echo "WARNING: OncoKB genomic-change query failed for $sample — querying those variants by protein change"
            GENOMIC_FAILED="--genomic-failed"
            python3 /opt/oncokb-annotator/MafAnnotator.py \
                -i "$MAF_BY_GENOMIC" \
                -o "${MAF_BY_GENOMIC%.tsv}.maf" \
                -b "$ONCOKB_TOKEN"
        fi
    fi
    python3 "$SCRIPT_DIR/oncokb_maf_query.py" merge \
        --protein "${MAF_BY_PROTEIN%.tsv}.maf" --genomic "${MAF_BY_GENOMIC%.tsv}.maf" \
        $GENOMIC_FAILED --out "$ANNOTATED_MAF" \
        || { echo "ERROR: OncoKB MafAnnotator outputs incomplete for $sample"; return 1; }
    rm -f "$MAF_BY_PROTEIN" "$MAF_BY_GENOMIC" "${MAF_BY_PROTEIN%.tsv}.maf" "${MAF_BY_GENOMIC%.tsv}.maf"

    # --- Merge the computed SpliceAI scores into the table ---
    # The column names are the plugin's, so everything downstream -- the tier split, the
    # ANNOTATION_FLAGS threshold, the classifier's O6 and its cryptic-site check -- is
    # unchanged whichever mode produced them. Rows the model did not score are left empty
    # rather than carrying a precomputed value, so the table never mixes the two sources
    # silently in one column.
    if [[ "$SPLICEAI_MODE" == "compute" && -n "$SPLICEAI_VCF" && -f "$SPLICEAI_VCF" ]]; then
        echo ">>> Merging computed SpliceAI scores: $sample"
        python3 "$SCRIPT_DIR/spliceai_patch_maf.py" \
            --maf "$ANNOTATED_MAF" --spliceai-vcf "$SPLICEAI_VCF" \
            --out "${ANNOTATED_MAF}.spliceai" --on-missing blank \
            && mv "${ANNOTATED_MAF}.spliceai" "$ANNOTATED_MAF" \
            || { echo "ERROR: could not merge SpliceAI scores for $sample"; return 1; }
    fi

    # --- GENIE case counts (optional annotation; user-built, licence-restricted table) ---
    if [[ -e "$GENIE_COUNTS" ]]; then
        echo ">>> GENIE annotation: $sample"
        python3 "$SCRIPT_DIR/genie_annotation.py" --maf "$ANNOTATED_MAF" --refs "$REF_DIR" \
            || { echo "ERROR: GENIE annotation failed for $sample"; return 1; }
    fi

    # --- Variant counts (written to per-sample temp file to avoid concurrent writes) ---
    local ORIGINAL_COUNT PASS_COUNT LIFTED_COUNT ANNOTATED_COUNT ONCOKB_COUNT TABLE_COUNT
    ORIGINAL_COUNT=$(zgrep -vc "^#" "$vcf")
    if [[ $DO_PASS_FILTER -eq 1 ]]; then
        PASS_COUNT=$(zgrep -vc "^#" "$FILTERED_VCF")
    else
        PASS_COUNT="$ORIGINAL_COUNT"
    fi
    if [[ $DO_LIFTOVER -eq 1 ]]; then
        LIFTED_COUNT=$(grep -vc "^#" "$LIFTED_VCF")
    else
        LIFTED_COUNT="$PASS_COUNT"
    fi
    ANNOTATED_COUNT=$(grep -vc "^#" "$ANNOTATED_VCF")
    ONCOKB_COUNT=$(grep -vc "^#" "$ONCOKB_VCF")
    TABLE_COUNT=$(( $(wc -l < "$TABLE_CSV") - 1 ))

    echo -e "${sample}\t${ORIGINAL_COUNT}\t${PASS_COUNT}\t${LIFTED_COUNT}\t${ANNOTATED_COUNT}\t${ONCOKB_COUNT}\t${TABLE_COUNT}" \
        > "${COUNT_FILE}.${sample}"

    # --- Field metadata, harvested while the VCFs still exist ---
    # A VCF documents its own fields, and the cleanup below deletes the annotated and OncoKB
    # VCFs, so the headers are read now and kept as a sidecar for post_analysis to describe
    # the columns from. One per sample: the merge step checks the samples agree.
    python3 "$SCRIPT_DIR/vcf_metadata.py" \
        --input "$vcf" --annotated "$ANNOTATED_VCF" --oncokb "$ONCOKB_VCF" \
        --out "$LOGS_DIR/${sample}.vcf_metadata.json" \
        || echo "WARNING: could not harvest field metadata for $sample; the two-row header will fall back to Config.yaml"

    # --- Cleanup ---
    if [[ $CLEAN_TMP -eq 1 ]]; then
        echo ">>> Cleaning temporary files for: $sample"
        rm -f "$FILTERED_VCF" "$FILTERED_VCF_TBI"
        rm -f "$LIFTED_VCF_GZ" "$LIFTED_VCF" "$LIFTED_VCF_TBI"
        rm -f "$REJECTED_VCF" "$REJECTED_VCF_TBI"
        rm -f "$ANNOTATED_VCF"
        rm -f "$ONCOKB_VCF"
    fi

    echo "--- $sample complete ---"
    } >> "$log_file" 2>&1
}

# ==============================================================================
# Parallel sample processing (--nodes controls concurrency).
# Uses wait -n (bash >= 4.3) to drain one slot before launching the next job.
# Exit statuses are collected via per-PID wait in the post-loop phase.
# ==============================================================================
declare -a ALL_PIDS=()
declare -a ALL_SAMPLES=()
RUNNING=0

for vcf in "$INPUT_DIR"/*.vcf.gz; do
    [[ -e "$vcf" ]] || continue
    basename_=$(basename "$vcf")
    sample="${basename_%.vcf.gz}"
    LOG_FILE="${LOGS_DIR}/${sample}.log"

    # Throttle: wait for one slot to free up before launching the next job
    if (( RUNNING >= JOBS )); then
        wait -n 2>/dev/null || true
        (( RUNNING-- )) || true
    fi

    echo ">>> Launching: $sample (log: $LOG_FILE)"
    process_sample "$vcf" "$LOG_FILE" &
    ALL_PIDS+=($!)
    ALL_SAMPLES+=("$sample")
    (( RUNNING++ )) || true
done

# Collect exit statuses for all launched jobs
PIPELINE_FAILED=0
for i in "${!ALL_PIDS[@]}"; do
    if ! wait "${ALL_PIDS[$i]}"; then
        echo "ERROR: Sample '${ALL_SAMPLES[$i]}' failed — see log: ${LOGS_DIR}/${ALL_SAMPLES[$i]}.log"
        PIPELINE_FAILED=1
    else
        echo ">>> OK: ${ALL_SAMPLES[$i]}"
    fi
done
[[ $PIPELINE_FAILED -eq 1 ]] && exit 1

# Merge per-sample count rows into the shared COUNT_FILE (preserving launch order)
for sample in "${ALL_SAMPLES[@]}"; do
    count_tmp="${COUNT_FILE}.${sample}"
    if [[ -f "$count_tmp" ]]; then
        cat "$count_tmp" >> "$COUNT_FILE"
        rm -f "$count_tmp"
    fi
done

if [[ $VEP_ONLY -eq 1 ]]; then
    echo "============================================================"
    echo "SMART PIPELINE COMPLETE (VEP-only mode)"
    echo "VEP-annotated VCFs are in: $ANNOTATED_DIR"
    echo "============================================================"
    exit 0
fi

echo "###################### Post Analysis ##########################################"
# post_analysis.py uses relative paths (./FINAL_Table/, ./output/) so it must
# run from OUTPUT_DIR_BASE where those directories were created.
cd "$OUTPUT_DIR_BASE"
# The OncoKB data version is a property of the run, not of the tool: it is applied to the
# MafAnnotator columns, which carry no header of their own.
RUN_ONCOKB_VERSION=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ONCOKB_DATA_VERSION") c=i; next}
                                  c && $c != "" {print $c}' FINAL_Table/*.maf 2>/dev/null \
                     | sort -u | head -1)
python3 "$SCRIPT_DIR/post_analysis.py" --config "$CONFIG_FILE" \
    --smart-version "$SMART_VERSION" \
    --vcf-metadata logs \
    ${RUN_ONCOKB_VERSION:+--oncokb-version "$RUN_ONCOKB_VERSION"} \
    $([[ $DO_LIFTOVER -eq 1 ]] && echo --liftover) \
    || { echo "ERROR: post_analysis.py failed — aborting"; exit 1; }

# --- Tier 4: SVIG-UK oncogenicity classification ---
# Separate from annotation (Tiers 1-3): reads the Tier 1 MAF and writes Final_result_tier4.tsv.
# Classifier resources under the reference directory are optional; codes whose resource is
# missing are reported as not evaluated in SVIG_UK_NOTES.
if [[ $DO_CLASSIFICATION -eq 1 ]]; then
    echo "###################### Tier 4 classification (SVIG-UK) ##########################"
    python3 "$SCRIPT_DIR/svig_uk_classifier.py" \
        --maf output/Final_result_tier1.maf \
        --tier4 output/Final_result_tier4.tsv \
        --tier4-json output/Final_result_tier4.json \
        --refs "$REF_DIR" \
        --tumour-type "$TUMOUR_TYPE" \
        --o4-on-target "$O4_ON_TARGET" \
        || { echo "ERROR: SVIG-UK classification failed — aborting"; exit 1; }
fi

if [[ $CLEAN_TABLES -eq 1 ]]; then
    echo ">>> Cleaning table files after post_analysis"
    rm -f "$TABLE_DIR"/*
    rm -f "$FINAL_TABLE_DIR"/*
fi

echo "============================================================"
echo "SMART PIPELINE COMPLETE"
echo "============================================================"