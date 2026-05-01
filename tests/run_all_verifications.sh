#!/usr/bin/env bash
# ==============================================================================
# run_all_verifications.sh — Run all SMART verification suites
#
# Runs unit tests + verification1, verification3, verification4, and reports results.
# verification2 is a manual compatibility check (multiple callers, no automated
# verify.py) and is excluded from automated runs.
#
# Usage:
#   export ONCOKB_TOKEN="<your-token>"
#   bash tests/run_all_verifications.sh
#
# Options:
#   --image IMAGE    Docker image to use (default: monkiky/smart:latest)
#   --refs  DIR      Reference directory   (default: /Volumes/ExternalSSD/refs)
#   --only  NAME     Run only one suite: unit | verification1 | verification3 | verification4
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

IMAGE="${SMART_IMAGE:-monkiky/smart:latest}"
REFS_DIR="${REFS_DIR:-/Volumes/ExternalSSD/refs}"
TOKEN="${ONCOKB_TOKEN:-}"
ONLY=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --image) IMAGE="$2"; shift 2 ;;
        --refs)  REFS_DIR="$2"; shift 2 ;;
        --only)  ONLY="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$TOKEN" ]]; then
    echo "ERROR: ONCOKB_TOKEN is not set."
    echo "       Run: export ONCOKB_TOKEN=your_token_here"
    exit 1
fi

PASS_COUNT=0
FAIL_COUNT=0
RESULTS=()

run_suite() {
    local name="$1"
    local script="$2"
    local results_tsv="$3"

    echo ""
    echo "████████████████████████████████████████████████████████████"
    echo "  ${name}"
    echo "████████████████████████████████████████████████████████████"

    if bash "$script"; then
        local status
        # Count MISMATCH/ERROR/SAME lines in results TSV (skip header)
        local failures
        failures=$(tail -n +2 "$results_tsv" 2>/dev/null \
            | awk -F'\t' '$2 ~ /MISMATCH|ERROR|SAME/' | wc -l | tr -d ' ')
        if [[ "$failures" -eq 0 ]]; then
            status="PASS"
            PASS_COUNT=$((PASS_COUNT + 1))
        else
            status="FAIL (${failures} issues)"
            FAIL_COUNT=$((FAIL_COUNT + 1))
        fi
    else
        status="FAIL (script error)"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi

    RESULTS+=("  ${name}: ${status}")
}

# ── Unit tests — Multi-transcript (no Docker, no token) ─────────────────────
if [[ -z "$ONLY" || "$ONLY" == "unit" ]]; then
    echo ""
    echo "================================================================"
    echo " UNIT TESTS: Multi-transcript feature (test_multi_transcript.py)"
    echo "================================================================"

    if python3 "$SCRIPT_DIR/test_multi_transcript.py" -v 2>&1; then
        RESULTS+=("  unit tests (multi-transcript): PASS")
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        RESULTS+=("  unit tests (multi-transcript): FAIL")
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
fi

# ── Verification 1 — Field-level API verification ────────────────────────────
if [[ -z "$ONLY" || "$ONLY" == "verification1" ]]; then
    V1_DIR="$SCRIPT_DIR/verification1"
    V1_OUTPUT="$V1_DIR/output"
    V1_RESULTS="$V1_DIR/results.tsv"

    mkdir -p "$V1_OUTPUT"

    echo ""
    echo "================================================================"
    echo " VERIFICATION 1: Running SMART pipeline"
    echo "================================================================"

    rm -rf "$V1_OUTPUT"
    mkdir -p "$V1_OUTPUT"

    docker run --rm \
        -v "$V1_DIR":/data \
        -v "$V1_OUTPUT":/output \
        -v "$REFS_DIR":/refs:ro \
        "$IMAGE" \
        "$TOKEN" \
        --transcripts-file /data/verification1_transcripts.txt \
        --config /data/Config.yaml \
        --ref-dir /refs \
        --input-dir /data \
        --no-liftover \
        --keep-tmp \
        --keep-tables

    echo ""
    echo "================================================================"
    echo " VERIFICATION 1: Running verify.py"
    echo "================================================================"

    V1_MAF="$V1_OUTPUT/output/Final_result_tier1.maf"
    if [[ ! -f "$V1_MAF" ]]; then
        echo "ERROR: MAF not found at $V1_MAF"
        RESULTS+=("  verification1: FAIL (MAF missing)")
        FAIL_COUNT=$((FAIL_COUNT + 1))
    else
        python3 "$V1_DIR/verify.py" \
            --maf "$V1_MAF" \
            --token "$TOKEN" \
            --output "$V1_RESULTS"

        failures=$(tail -n +2 "$V1_RESULTS" 2>/dev/null \
            | awk -F'\t' '$2 ~ /MISMATCH|ERROR|COVERAGE_GAP/' | wc -l | tr -d ' ')
        if [[ "$failures" -eq 0 ]]; then
            RESULTS+=("  verification1: PASS")
            PASS_COUNT=$((PASS_COUNT + 1))
        else
            RESULTS+=("  verification1: FAIL (${failures} issues — see $V1_RESULTS)")
            FAIL_COUNT=$((FAIL_COUNT + 1))
        fi
    fi
fi

# ── Verification 3 — Transcript prioritisation impact ───────────────────────
if [[ -z "$ONLY" || "$ONLY" == "verification3" ]]; then
    V3_DIR="$SCRIPT_DIR/verification3"
    V3_RESULTS="$V3_DIR/results.tsv"

    echo ""
    echo "================================================================"
    echo " VERIFICATION 3: Transcript prioritisation impact"
    echo "================================================================"

    SMART_IMAGE="$IMAGE" REFS_DIR="$REFS_DIR" ONCOKB_TOKEN="$TOKEN" \
        bash "$V3_DIR/run_verification3.sh"

    failures=$(tail -n +2 "$V3_RESULTS" 2>/dev/null \
        | awk -F'\t' '$2 ~ /MISMATCH|ERROR|SAME/' | wc -l | tr -d ' ')
    if [[ "$failures" -eq 0 ]]; then
        RESULTS+=("  verification3: PASS")
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        RESULTS+=("  verification3: FAIL (${failures} issues — see $V3_RESULTS)")
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
fi

# ── Verification 4 — Parallel processing (--jobs) ───────────────────────────
if [[ -z "$ONLY" || "$ONLY" == "verification4" ]]; then
    V4_DIR="$SCRIPT_DIR/verification4"
    V4_RESULTS="$V4_DIR/results.tsv"

    echo ""
    echo "================================================================"
    echo " VERIFICATION 4: Parallel processing (--jobs 2)"
    echo "================================================================"

    SMART_IMAGE="$IMAGE" REFS_DIR="$REFS_DIR" ONCOKB_TOKEN="$TOKEN" \
        bash "$V4_DIR/run_verification4.sh"

    failures=$(tail -n +2 "$V4_RESULTS" 2>/dev/null \
        | awk -F'\t' '$2 == "FAIL"' | wc -l | tr -d ' ')
    if [[ "$failures" -eq 0 ]]; then
        RESULTS+=("  verification4: PASS")
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        RESULTS+=("  verification4: FAIL (${failures} failed checks — see $V4_RESULTS)")
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "████████████████████████████████████████████████████████████"
echo "  VERIFICATION SUMMARY"
echo "████████████████████████████████████████████████████████████"
for r in "${RESULTS[@]}"; do
    echo "$r"
done
echo ""
echo "  Passed: ${PASS_COUNT}   Failed: ${FAIL_COUNT}"
echo "████████████████████████████████████████████████████████████"

[[ "$FAIL_COUNT" -eq 0 ]]
