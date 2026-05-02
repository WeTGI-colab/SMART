#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# SMART — Somatic Mutation Annotation and Reporting Tool
# ==============================================================================

usage() {
cat <<EOF
SMART — Somatic Mutation Annotation and Reporting Tool (Dockerised)

Usage:
  docker run -v /path/to/data:/data -v /path/to/refs:/refs smart \\
    <ONCOKB_TOKEN> --transcripts-file FILE --ref-dir /refs --config FILE [OPTIONS]

Required:
  <ONCOKB_TOKEN>           OncoKB API token
  --transcripts-file FILE  Transcript list for prioritisation
  --ref-dir DIR            Reference resources directory (mounted at /refs)
  --config FILE            YAML config file for post_analysis (Config.yaml)

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

Options:
  --pass / --no-pass               Enable/disable PASS filtering (default: ON)
  --liftover / --no-liftover       Enable/disable liftover (default: ON)
  --clean-tmp / --keep-tmp         Delete/keep intermediate files (default: clean)
  --clean-tables / --keep-tables   Delete/keep table files after post_analysis (default: clean)
  --jobs N                         Number of samples to process in parallel on this
                                   machine (default: 1 = sequential). Each job is one
                                   sample; all jobs share the same CPU and RAM.
  --vep-only                       Stop after VEP annotation. Skips OncoKB,
                                   vcf2table, MafAnnotator, and post-analysis.
                                   Annotated VCFs are written to AnnotatedVcf/.
                                   OncoKB token is still required but not used.
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

ONCOKB_TOKEN="$1"
shift

# Defaults
DO_PASS_FILTER=1
DO_LIFTOVER=1
CLEAN_TMP=1
CLEAN_TABLES=1
JOBS=1
VEP_ONLY=0
TRANSCRIPTS_FILE=""
REF_DIR=""
CONFIG_FILE=""
OUTPUT_DIR_BASE="/output"

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

# I get oncoKB versions dinamically, I hope they dont change this frecuently
ONCOKB_INFO=$(curl -s https://www.oncokb.org/api/v1/info || echo "")
ONCOKB_DATA_VERSION=$(echo "$ONCOKB_INFO" | grep -oP '"dataVersion":\{"version":"\K[^"]+' || echo "unknown")
ONCOKB_DATA_DATE=$(echo "$ONCOKB_INFO" | grep -oP '"dataVersion":\{"version":"[^"]+","date":"\K[^"]+' || echo "unknown")
ONCOKB_API_VERSION=$(echo "$ONCOKB_INFO" | grep -oP '"apiVersion":\{"version":"\K[^"]+' || echo "unknown")

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
echo "  SpliceAI SNV:      $(basename "$SPLICEAI_SNV")"
echo "  SpliceAI INDEL:    $(basename "$SPLICEAI_INDEL")"
echo "  REVEL:             $(basename "$REVEL_FILE")"
echo "  LOEUF:             $(basename "$LOEUF_FILE")"
echo "  OncoKB data:       $ONCOKB_DATA_VERSION ($ONCOKB_DATA_DATE)"

echo ""

echo ""
echo "Config:"
echo "  PASS filter:      $([[ $DO_PASS_FILTER -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  Liftover:         $([[ $DO_LIFTOVER -eq 1 ]] && echo ENABLED || echo DISABLED)"
echo "  VEP only:         $([[ $VEP_ONLY -eq 1 ]] && echo YES — stopping after VEP || echo NO)"
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
        --plugin SpliceAI,snv="$SPLICEAI_SNV",indel="$SPLICEAI_INDEL",cutoff=0.5 \
        --plugin REVEL,"$REVEL_FILE" \
        --custom "$CLINVAR_VCF",ClinVar,vcf,exact,0,AF_ESP,AF_EXAC,AF_TGP,ALLELEID,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNSIGSCV,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ONCDN,ONCDNINCL,ONCDISDB,ONCDISDBINCL,ONC,ONCINCL,ONCREVSTAT,ONCSCV,ONCCONF,ORIGIN,RS,SCIDN,SCIDNINCL,SCIDISDB,SCIDISDBINCL,SCIREVSTAT,SCI,SCIINCL,SCISCV \
        --custom "$CIVIC_VCF",CIViC,vcf,exact,0,GN,VT,CSQ \
        --plugin LOEUF,file="$LOEUF_FILE",match_by=transcript \
        --custom "$CANCER_HOTSPOTS_VCF",CancerHotspots,vcf,exact,0,HOTSPOT,HOTSPOT_GENE,HOTSPOT_HGVSp,HOTSPOT3D,HOTSPOT3D_GENE,HOTSPOT3D_HGVSp,HOTSPOTNC,HOTSPOTNC_GENE,HOTSPOTNC_HGVSc

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
    python3 /opt/oncokb-annotator/MafAnnotator.py \
        -i "$TABLE_CSV" \
        -o "$ANNOTATED_MAF" \
        -b "$ONCOKB_TOKEN"

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
python3 "$SCRIPT_DIR/post_analysis.py" --config "$CONFIG_FILE" \
    --smart-version "$SMART_VERSION" \
    || { echo "ERROR: post_analysis.py failed — aborting"; exit 1; }

if [[ $CLEAN_TABLES -eq 1 ]]; then
    echo ">>> Cleaning table files after post_analysis"
    rm -f "$TABLE_DIR"/*
    rm -f "$FINAL_TABLE_DIR"/*
fi

echo "============================================================"
echo "SMART PIPELINE COMPLETE"
echo "============================================================"