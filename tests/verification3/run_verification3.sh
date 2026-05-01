#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_verification3.sh — Transcript prioritisation impact verification
#
# Runs the SMART pipeline twice on the same VCF with two different
# preferred-transcript files, then compares the VEP annotations between runs.
#
# Usage:
#   export ONCOKB_TOKEN="<your-token>"
#   bash tests/verification3/run_verification3.sh
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

IMAGE="${SMART_IMAGE:-monkiky/smart:latest}"
REFS_DIR="${REFS_DIR:-/Volumes/ExternalSSD/refs}"
TOKEN="${ONCOKB_TOKEN:-}"

OUTPUT_A="$SCRIPT_DIR/output_A"
OUTPUT_B="$SCRIPT_DIR/output_B"

mkdir -p "$OUTPUT_A" "$OUTPUT_B"

echo "================================================================"
echo "RUN A — Transcript set A (clinically preferred)"
echo "  EGFR: NM_005228.5  → p.Leu858Arg"
echo "  TP53: NM_000546.6  → p.Arg175His"
echo "  CNA:  NM_058195.4  → CDKN2A"
echo "================================================================"

docker run --rm \
  -v "$SCRIPT_DIR":/data \
  -v "$OUTPUT_A":/output \
  -v "$REFS_DIR":/refs:ro \
  "$IMAGE" \
  "$TOKEN" \
  --transcripts-file /data/transcripts_A.txt \
  --ref-dir /refs \
  --config /data/Config.yaml \
  --input-dir /data \
  --output-dir /output \
  --no-liftover \
  --keep-tmp 2>&1 | tee "$OUTPUT_A/smart_A.log"

echo ""
echo "================================================================"
echo "RUN B — Transcript set B (alternative — demonstrates impact)"
echo "  EGFR: NM_001346897.2  → p.Leu813Arg  (wrong position for OncoKB)"
echo "  TP53: NM_001126115.2  → p.Arg43His   (completely different residue)"
echo "  CNA:  NM_004936.4     → CDKN2B       (different gene reported)"
echo "================================================================"

docker run --rm \
  -v "$SCRIPT_DIR":/data \
  -v "$OUTPUT_B":/output \
  -v "$REFS_DIR":/refs:ro \
  "$IMAGE" \
  "$TOKEN" \
  --transcripts-file /data/transcripts_B.txt \
  --ref-dir /refs \
  --config /data/Config.yaml \
  --input-dir /data \
  --output-dir /output \
  --no-liftover \
  --keep-tmp 2>&1 | tee "$OUTPUT_B/smart_B.log"

echo ""
echo "================================================================"
echo "VERIFICATION — Comparing annotations between runs"
echo "================================================================"

MAF_A="$OUTPUT_A/output/Final_result_tier1.maf"
MAF_B="$OUTPUT_B/output/Final_result_tier1.maf"

python3 "$SCRIPT_DIR/verify.py" \
  --maf-a  "$MAF_A" \
  --maf-b  "$MAF_B" \
  --transcripts-a "$SCRIPT_DIR/transcripts_A.txt" \
  --transcripts-b "$SCRIPT_DIR/transcripts_B.txt" \
  --output "$SCRIPT_DIR/results.tsv"
