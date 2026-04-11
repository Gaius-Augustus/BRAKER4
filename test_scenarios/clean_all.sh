#!/bin/bash
# Remove BRAKER4 test scenario outputs so tests can be re-run from scratch.
# Does NOT delete config files, samples.csv, run_test.sh, or rulegraph images.
#
# Only touches scenarios that are part of the official BRAKER4 test suite.
# Other scenario_* directories (e.g. user benchmarks) are left alone.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Only these scenarios are cleaned. Add new official scenarios here.
OFFICIAL_SCENARIOS=(
    scenario_08_es
    scenario_09_ep_masking
    scenario_10_et_varus
    scenario_11_etp_sra
    scenario_12_multi_mode
)

echo "Cleaning BRAKER4 test scenario outputs..."

count=0
for name in "${OFFICIAL_SCENARIOS[@]}"; do
    scenario="$SCRIPT_DIR/$name"
    [ -d "$scenario" ] || continue
    removed=""

    for dir in output logs benchmarks augustus_config shared_data .snakemake .ncbi; do
        if [ -d "$scenario/$dir" ]; then
            rm -rf "$scenario/$dir"
            removed="$removed $dir"
        fi
    done

    # Remove AGAT log files
    rm -f "$scenario"/braker.agat.log
    rm -f "$scenario"/.wget-hsts

    if [ -n "$removed" ]; then
        echo "  $name: removed$removed"
        count=$((count + 1))
    fi
done

if [ $count -eq 0 ]; then
    echo "  Nothing to clean."
else
    echo "Cleaned $count scenario(s)."
fi
