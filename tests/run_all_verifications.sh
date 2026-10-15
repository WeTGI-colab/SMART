#!/usr/bin/env bash
# ==============================================================================
# run_all_verifications.sh — Run all SMART verification suites
#
# Runs unit tests + verification1, verification3, verification4, verification5,
# localsources, and reports results.
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
#   --only  NAME     Run only one suite: unit | verification1 | verification3 | verification4 |
#                    verification5 | localsources
#                    (unit and verification5 need no Docker image, refs, or token;
#                     localsources needs --refs and verification1's output, but no Docker
#                     and no token except for its gnomAD module, which calls the VEP API)
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

# unit and verification5 are pure-Python suites that need no token; only require
# the token when a Docker-based suite will actually run.
if [[ -z "$TOKEN" && "$ONLY" != "unit" && "$ONLY" != "verification5" ]]; then
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

# ── Verification 5 — Core pipeline unit tests (pytest, no Docker/token/refs) ──
if [[ -z "$ONLY" || "$ONLY" == "verification5" ]]; then
    echo ""
    echo "================================================================"
    echo " VERIFICATION 5: Core pipeline unit tests (pytest)"
    echo "================================================================"

    if bash "$SCRIPT_DIR/verification5/run_verification5.sh"; then
        RESULTS+=("  verification5 (unit, pytest): PASS")
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        RESULTS+=("  verification5 (unit, pytest): FAIL")
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
fi

# These are paths, not results, so they are defined unconditionally. Setting them inside
# the verification1 branch meant that any `--only` run reached the report section at the
# end with V1_DIR unset and aborted under `set -u`, after the summary had already printed
# PASS -- so the run looked successful and then failed on its last line.
V1_DIR="$SCRIPT_DIR/verification1"
V1_OUTPUT="$V1_DIR/output"
V1_RESULTS="$V1_DIR/results.tsv"

# ── Verification 1 — Field-level API verification ────────────────────────────
if [[ -z "$ONLY" || "$ONLY" == "verification1" ]]; then

    mkdir -p "$V1_OUTPUT"

    echo ""
    echo "================================================================"
    echo " VERIFICATION 1: Running SMART pipeline"
    echo "================================================================"

    rm -rf "$V1_OUTPUT"
    mkdir -p "$V1_OUTPUT"

    # The default --ref-check fail, deliberately. Issue #19 is fixed: the three records whose
    # REF did not match GRCh38 have been corrected, so the fixture passes the check and the
    # suite now fails if it ever regresses. This previously carried --ref-check warn; do not
    # put it back to make a failing fixture pass.

    docker run --rm \
        -v "$V1_DIR":/data \
        -v "$V1_OUTPUT":/results \
        -v "$REFS_DIR":/refs:ro \
        "$IMAGE" \
        "$TOKEN" \
        --transcripts-file /data/verification1_transcripts.txt \
        --config /data/Config.yaml \
        --ref-dir /refs \
        --output-dir /results \
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
        # `|| true` because verify.py's exit status is read below rather than aborting the
        # run here: a failure in this suite must not stop verification3, verification4 and
        # the local-source checks from executing, which is what used to happen.
        python3 "$V1_DIR/verify.py" \
            --maf "$V1_MAF" \
            --token "$TOKEN" \
            --output "$V1_RESULTS" || true

        # The Result column, found by name. This previously tested $2 -- the Sample column,
        # which never contains any of these words -- so the count was always zero and
        # verification1 was recorded as PASS unconditionally, whatever verify.py reported.
        if [[ -f "$V1_RESULTS" ]]; then
            mismatches=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Result") c=i; next}
                                      c && $c ~ /^MISMATCH$/ {n++} END{print n+0}' "$V1_RESULTS")
            notes=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Result") c=i; next}
                                 c && $c ~ /^(ERROR|COVERAGE_GAP)/ {n++} END{print n+0}' "$V1_RESULTS")
        else
            mismatches=1; notes=0
        fi

        # Empty-column audit: every column that carries no value on any variant must be
        # listed with a reason. ONCOKB_VUS was empty on every row of every run for years
        # because the code read the API's "vus" field as "VUS", and nothing could notice --
        # an absent key, a null value and an inapplicable field are the same empty cell.
        echo ""
        echo "----------------------------------------------------------------"
        echo " VERIFICATION 1: Empty-column audit"
        echo "----------------------------------------------------------------"
        # Both files, and the second is the one that matters. The audit used to run only on
        # the intermediate FINAL_Table MAF, which is written before post_analysis -- so a
        # column that post_analysis itself emptied passed the audit. End_Position was blank
        # on 44 of 44 rows of the Tier 1 MAF, the file the laboratory reads, and this audit
        # existed to catch exactly that and was pointed one step too early.
        V1_ONCOKB_MAF="$V1_OUTPUT/FINAL_Table/verification1.oncoKB.maf"
        V1_TIER1_MAF="$V1_OUTPUT/output/Final_result_tier1.maf"
        AUDIT_EXIT=0
        # The Tier 1 MAF takes a second allowlist. post_analysis pads the OncoKB indexed
        # columns out to a fixed count, so patterns that are needed there are stale for
        # FINAL_Table, where every slot that exists is populated.
        for maf in "$V1_ONCOKB_MAF" "$V1_TIER1_MAF"; do
            if [[ -f "$maf" ]]; then
                echo "--- $(basename "$(dirname "$maf")")/$(basename "$maf")"
                EXTRA_ALLOW=()
                if [[ "$maf" == "$V1_TIER1_MAF" ]]; then
                    EXTRA_ALLOW=(--allowlist "$V1_DIR/empty_columns_allowed_tier1.tsv")
                fi
                python3 "$PROJECT_DIR/scripts/audit_empty_columns.py" \
                    --maf "$maf" \
                    --allowlist "$V1_DIR/empty_columns_allowed.tsv" \
                    ${EXTRA_ALLOW[@]+"${EXTRA_ALLOW[@]}"} || AUDIT_EXIT=$?
            else
                echo "SKIP: $maf not found"
            fi
        done

        # Transcript provenance: per column, did the value come from the transcript the user
        # supplied? Reported rather than enforced, because the pipeline has columns that are
        # legitimately not transcript-specific and 29 dbNSFP columns that are not resolved to
        # the reported transcript at all (issue #29). The point is that it is stated on every
        # run instead of being invisible.
        echo ""
        echo "----------------------------------------------------------------"
        echo " VERIFICATION 1: Transcript provenance"
        echo "----------------------------------------------------------------"
        if [[ -f "$V1_ONCOKB_MAF" ]]; then
            python3 "$PROJECT_DIR/scripts/transcript_provenance.py" \
                --maf "$V1_ONCOKB_MAF" \
                --transcripts "$V1_DIR/verification1_transcripts.txt" || true
        else
            echo "SKIP: $V1_ONCOKB_MAF not found"
        fi

        if [[ "$mismatches" -eq 0 && "$AUDIT_EXIT" -eq 0 ]]; then
            if [[ "$notes" -gt 0 ]]; then
                RESULTS+=("  verification1: PASS (${notes} unanswered API calls / known coverage gaps)")
            else
                RESULTS+=("  verification1: PASS")
            fi
            PASS_COUNT=$((PASS_COUNT + 1))
        else
            detail=""
            [[ "$mismatches" -gt 0 ]] && detail="${mismatches} mismatch(es)"
            [[ "$AUDIT_EXIT" -ne 0 ]] && detail="${detail:+$detail, }empty-column audit failed"
            RESULTS+=("  verification1: FAIL (${detail} — see $V1_RESULTS)")
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

# ── Local source verification — annotation against the files it came from ────
# verification1 reaches the sources that have a REST API. This one reaches the rest:
# ClinVar, Cancer Hotspots, SpliceAI, REVEL, LOEUF, dbNSFP and GENIE are compared against
# the files on disk the pipeline annotated from, and gnomAD against the VEP API, since it
# has no file (it comes from the VEP cache). Needs verification1 to have run.
if [[ -z "$ONLY" || "$ONLY" == "localsources" ]]; then
    echo ""
    echo "================================================================"
    echo " LOCAL SOURCE VERIFICATION: annotation vs. its own source files"
    echo "================================================================"

    LOCAL_MAF="$SCRIPT_DIR/verification1/output/FINAL_Table/verification1.oncoKB.maf"
    if [[ ! -f "$LOCAL_MAF" ]]; then
        echo "SKIP: $LOCAL_MAF not found — run verification1 first"
        RESULTS+=("  localsources: SKIP (no verification1 output)")
    elif python3 "$SCRIPT_DIR/verify_local_sources.py" \
            --maf "$LOCAL_MAF" --refs "$REFS_DIR" \
            --output "$SCRIPT_DIR/verification1/local_source_results.tsv"; then
        n_err=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Result") c=i; next}
                              c && $c ~ /^ERROR/ {n++} END{print n+0}' \
                "$SCRIPT_DIR/verification1/local_source_results.tsv" 2>/dev/null || echo 0)
        if [[ "$n_err" -gt 0 ]]; then
            RESULTS+=("  localsources: PASS (${n_err} unanswered API calls)")
        else
            RESULTS+=("  localsources: PASS")
        fi
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        RESULTS+=("  localsources: FAIL (mismatches — see local_source_results.tsv)")
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

# ==============================================================================
# One publishable report from everything the run produced. Composed last, so the
# local-source checks and verification3 have written their results first. It reads only
# files the run wrote, and says so where one is missing rather than omitting the section.
# ==============================================================================
if [[ -f "$V1_DIR/output/output/Final_result_tier1.maf" ]]; then
    echo ""
    echo "================================================================"
    echo " VERIFICATION REPORT"
    echo "================================================================"
    python3 "$PROJECT_DIR/scripts/verification_report.py" \
        --verification1 "$V1_DIR" \
        --vcf "$V1_DIR/verification1.vcf" \
        --verification3 "$PROJECT_DIR/tests/verification3/results.tsv" \
        --allowlist "$V1_DIR/empty_columns_allowed.tsv" \
        --out "$PROJECT_DIR/tests/VERIFICATION_REPORT.md" \
        --pdf "$PROJECT_DIR/tests/VERIFICATION_REPORT.pdf" || true
fi

# The suite's verdict, asserted last so it is this script's exit status. It used to be the
# final line; the report block was appended after it, which made the report's status the
# script's and silently turned every failing run into a passing one.
[[ "$FAIL_COUNT" -eq 0 ]]
