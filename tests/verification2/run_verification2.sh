#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_verification2.sh — Run SMART on a selected verification2 VCF
#
# Usage:
#   ./run_verification2.sh                  # interactive menu
#   ./run_verification2.sh pisces           # run Pisces directly
#   ./run_verification2.sh mutect2          # run MuTect2 directly
#   ./run_verification2.sh strelka          # run Strelka directly
#   ./run_verification2.sh somaticsniper    # run SomaticSniper (--no-pass)
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR"
OUTPUT_DIR="$SCRIPT_DIR/output"
REFS_DIR="/Volumes/ExternalSSD/refs"
CONFIG="$SCRIPT_DIR/Config.yaml"
TRANSCRIPTS="$SCRIPT_DIR/verification2_transcripts.txt"
IMAGE="monkiky/smart:1.1.0"
VERIFY_SCRIPT="$(dirname "$SCRIPT_DIR")/verify_maf.py"

if [[ -z "${ONCOKB_TOKEN:-}" ]]; then
    echo "[ERROR] ONCOKB_TOKEN environment variable is not set."
    exit 1
fi

# ---------------------------------------------------------------------------
# Select caller
# ---------------------------------------------------------------------------
if [[ $# -ge 1 ]]; then
    CALLER="$(echo "$1" | tr '[:upper:]' '[:lower:]')"
else
    echo ""
    echo "Select a variant caller to test:"
    echo ""
    echo "  [1] pisces         Pisces 5.2.10 — tumour-only amplicon, single sample, GRCh38"
    echo "  [2] mutect2        GATK MuTect2 — tumour/normal WGS, SEQC2 benchmark, GRCh38"
    echo "  [3] strelka        Illumina Strelka — tumour/normal WGS, SEQC2 benchmark, GRCh38"
    echo "  [4] somaticsniper  SomaticSniper — tumour/normal WGS, no PASS filter (--no-pass)"
    echo ""
    read -rp "Enter number [1-4]: " CHOICE
    case "$CHOICE" in
        1) CALLER="pisces" ;;
        2) CALLER="mutect2" ;;
        3) CALLER="strelka" ;;
        4) CALLER="somaticsniper" ;;
        *) echo "[ERROR] Invalid choice."; exit 1 ;;
    esac
fi

case "$CALLER" in
    pisces)
        VCF="To_test_PISCIS.vcf.gz"
        EXTRA=""
        DESC="Pisces 5.2.10 — tumour-only amplicon, single sample, GRCh38"
        ;;
    mutect2)
        VCF="FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.muTect2.subsample.vcf.gz"
        EXTRA=""
        DESC="GATK MuTect2 — tumour/normal WGS, SEQC2 benchmark, GRCh38"
        ;;
    strelka)
        VCF="FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.strelka.subsample.vcf.gz"
        EXTRA=""
        DESC="Illumina Strelka — tumour/normal WGS, SEQC2 benchmark, GRCh38"
        ;;
    somaticsniper)
        VCF="FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.somaticSniper.subsample.vcf.gz"
        EXTRA="--no-pass"
        DESC="SomaticSniper — tumour/normal WGS, no PASS filter"
        ;;
    *)
        echo "[ERROR] Unknown caller: '$CALLER'. Choose from: pisces, mutect2, strelka, somaticsniper"
        exit 1
        ;;
esac

LOG="$OUTPUT_DIR/smart_${CALLER}.log"

echo ""
echo "============================================================"
echo "  SMART verification2 — $CALLER"
echo "  $DESC"
echo "  VCF: $VCF"
[[ -n "$EXTRA" ]] && echo "  Extra flags: $EXTRA"
echo "============================================================"
echo ""

# Check VCF exists
if [[ ! -f "$DATA_DIR/$VCF" ]]; then
    echo "[ERROR] VCF not found: $DATA_DIR/$VCF"
    exit 1
fi

# Check tabix index
if [[ ! -f "$DATA_DIR/$VCF.tbi" ]]; then
    echo "[INFO] No tabix index found — creating it..."
    docker run --rm \
        -v "$DATA_DIR:/data" \
        quay.io/biocontainers/htslib:1.21--h566b1c6_1 \
        tabix -p vcf "/data/$VCF"
    echo "[INFO] Index created."
fi

mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# Run SMART — mount only the selected VCF via a temp symlink dir
# ---------------------------------------------------------------------------
TMPDIR_INPUT="$(mktemp -d "$SCRIPT_DIR/tmp_input_XXXXXX")"
trap 'rm -rf "$TMPDIR_INPUT"' EXIT

cp "$DATA_DIR/$VCF"         "$TMPDIR_INPUT/$VCF"
cp "$DATA_DIR/$VCF.tbi"     "$TMPDIR_INPUT/$VCF.tbi"
cp "$CONFIG"                   "$TMPDIR_INPUT/Config.yaml"
cp "$TRANSCRIPTS"              "$TMPDIR_INPUT/transcripts.txt"

echo "[INFO] Running SMART... (log: $LOG)"
echo ""

docker run --rm \
    -v "$TMPDIR_INPUT:/data" \
    -v "$OUTPUT_DIR:/output" \
    -v "$REFS_DIR:/refs:ro" \
    "$IMAGE" \
    "$ONCOKB_TOKEN" \
    --transcripts-file /data/transcripts.txt \
    --ref-dir /refs \
    --config /data/Config.yaml \
    --input-dir /data \
    --output-dir /output \
    --no-liftover \
    --keep-tmp \
    $EXTRA \
    2>&1 | tee "$LOG"


