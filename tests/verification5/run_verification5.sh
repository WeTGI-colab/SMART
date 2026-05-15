#!/usr/bin/env bash
# Verification 5 — Unit tests for core pipeline logic
# Runs pytest from the verification5 directory, no Docker or external APIs required.
# Usage:  bash tests/verification5/run_verification5.sh [--coverage]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ensure pytest is available
if ! command -v pytest &>/dev/null; then
    echo "ERROR: pytest not found. Install test dependencies:"
    echo "       pip install pytest pytest-cov"
    exit 1
fi

cd "$SCRIPT_DIR"

if [[ "${1:-}" == "--coverage" ]]; then
    echo "Running unit tests with coverage report..."
    pytest . -v \
        --cov=../../scripts \
        --cov-report=term-missing \
        --cov-report=html:../../coverage_html \
        --tb=short
    echo ""
    echo "HTML coverage report written to: coverage_html/index.html"
else
    echo "Running unit tests..."
    pytest . -v --tb=short
fi
