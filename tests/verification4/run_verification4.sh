#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# run_verification4.sh — Parallel processing verification (--nodes)
#
# Copies verification1.vcf.gz twice under different sample names, runs the
# SMART pipeline with --nodes 2, then verify.py checks:
#   - both per-sample log files exist and contain the completion marker
#   - variant_counts.txt has exactly 2 data rows
#   - both output MAFs exist and are content-identical (same input → same output,
#     only Tumor_Sample_Barcode differs)
#
# Usage:
#   export ONCOKB_TOKEN="<your-token>"
#   bash tests/verification4/run_verification4.sh
# ------------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V1_DIR="$(cd "$SCRIPT_DIR/../verification1" && pwd)"

IMAGE="${SMART_IMAGE:-monkiky/smart:latest}"
REFS_DIR="${REFS_DIR:-/Volumes/ExternalSSD/refs}"
TOKEN="${ONCOKB_TOKEN:-}"

if [[ -z "$TOKEN" ]]; then
    echo "ERROR: ONCOKB_TOKEN is not set."
    echo "       Run: export ONCOKB_TOKEN=your_token_here"
    exit 1
fi

VCF_DIR="$SCRIPT_DIR/vcfs"
OUTPUT_DIR="$SCRIPT_DIR/output"

mkdir -p "$VCF_DIR" "$OUTPUT_DIR"

echo "================================================================"
echo " VERIFICATION 4: Parallel processing (--nodes 2)"
echo "================================================================"

echo ""
echo ">>> Preparing two sample copies from verification1 VCF..."
cp "$V1_DIR/verification1.vcf.gz"     "$VCF_DIR/par_sample_A.vcf.gz"
cp "$V1_DIR/verification1.vcf.gz.tbi" "$VCF_DIR/par_sample_A.vcf.gz.tbi"
cp "$V1_DIR/verification1.vcf.gz"     "$VCF_DIR/par_sample_B.vcf.gz"
cp "$V1_DIR/verification1.vcf.gz.tbi" "$VCF_DIR/par_sample_B.vcf.gz.tbi"

echo ">>> Running SMART pipeline with --nodes 2..."
docker run --rm \
    -v "$V1_DIR":/v1data:ro \
    -v "$VCF_DIR":/vcfs \
    -v "$OUTPUT_DIR":/output \
    -v "$REFS_DIR":/refs:ro \
    "$IMAGE" \
    "$TOKEN" \
    --transcripts-file /v1data/verification1_transcripts.txt \
    --config /v1data/Config.yaml \
    --ref-dir /refs \
    --input-dir /vcfs \
    --output-dir /output \
    --no-liftover \
    --keep-tmp \
    --keep-tables \
    --jobs 2

echo ""
echo "================================================================"
echo " VERIFICATION 4: Running verify.py"
echo "================================================================"

python3 "$SCRIPT_DIR/verify.py" \
    --output-dir "$OUTPUT_DIR" \
    --samples par_sample_A par_sample_B \
    --results "$SCRIPT_DIR/results.tsv"
