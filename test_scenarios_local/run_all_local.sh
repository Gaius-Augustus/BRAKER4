#!/bin/bash
# Run all local test scenarios sequentially.
# No SLURM, no VARUS, no SRA downloads, no repeat masking.
# 8 threads by default.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Download test data first
echo "Downloading test data if needed..."
source "$PIPELINE_DIR/test_data/download_test_data.sh"
echo ""

SCENARIOS=(
    scenario_01_es
    scenario_03_et_fastq
    scenario_04_etp_bam
    scenario_05_isoseq_bam
    scenario_06_isoseq_fastq
    scenario_07_dual
    scenario_02_ep              # last: runs AUGUSTUS optimization
)

PASSED=0
FAILED=0
FAILED_LIST=""

for scenario in "${SCENARIOS[@]}"; do
    echo ""
    echo "======================================================================"
    echo "  Running: $scenario"
    echo "======================================================================"
    echo ""

    if bash "$SCRIPT_DIR/$scenario/run_test.sh"; then
        echo ""
        echo "PASSED: $scenario"
        PASSED=$((PASSED + 1))
    else
        echo ""
        echo "FAILED: $scenario"
        FAILED=$((FAILED + 1))
        FAILED_LIST="$FAILED_LIST  - $scenario\n"
    fi
done

echo ""
echo "======================================================================"
echo "  LOCAL TEST SUMMARY"
echo "======================================================================"
echo "  Passed: $PASSED / $((PASSED + FAILED))"
echo "  Failed: $FAILED / $((PASSED + FAILED))"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "  Failed scenarios:"
    echo -e "$FAILED_LIST"
    exit 1
else
    echo ""
    echo "  All local tests passed."
fi
