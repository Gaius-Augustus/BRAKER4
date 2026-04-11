#!/bin/bash
# Remove all local test scenario outputs so tests can be re-run from scratch.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Cleaning local test scenario outputs..."

count=0
for scenario in "$SCRIPT_DIR"/scenario_*/; do
    name=$(basename "$scenario")
    removed=""

    for dir in output logs benchmarks augustus_config shared_data .snakemake .ncbi; do
        if [ -d "$scenario/$dir" ]; then
            rm -rf "$scenario/$dir"
            removed="$removed $dir"
        fi
    done

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
