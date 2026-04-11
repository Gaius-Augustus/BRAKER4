#!/bin/bash
# Generate rule graph PNGs for all test scenarios (local + HPC)
# Requires: snakemake, graphviz (dot)
#
# Usage: bash generate_rulegraphs.sh
#
# This script creates a rulegraph.png in each scenario directory showing
# the order of rule execution for that scenario's configuration.
# Run this on a machine with graphviz installed (not needed on HPC).

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check for graphviz
if ! command -v dot &> /dev/null; then
    echo "ERROR: graphviz (dot) is not installed."
    echo "Install with: apt install graphviz / conda install graphviz"
    exit 1
fi

PASSED=0
FAILED=0

# Create dummy files for rulegraph generation (removed after)
# These are touched so Snakemake can resolve the DAG without real data.
DUMMY_FILES=()
for f in \
    "$PIPELINE_DIR/test_data/Ostreococcus_tauri.fa" \
    "$PIPELINE_DIR/test_data/Viridiplantae.fa" \
    "$PIPELINE_DIR/test_data/genome.fa" \
    "$PIPELINE_DIR/test_data/proteins.fa" \
    "$PIPELINE_DIR/test_data/RNAseq.bam" \
    "$PIPELINE_DIR/test_data/reads_R1.fastq.gz" \
    "$PIPELINE_DIR/test_data/reads_R2.fastq.gz" \
    "$PIPELINE_DIR/test_data/isoseq.bam" \
    "$PIPELINE_DIR/test_data/isoseq_lib2.bam" \
    "$PIPELINE_DIR/test_data/isoseq.fastq.gz"; do
    if [ ! -e "$f" ]; then
        touch "$f"
        DUMMY_FILES+=("$f")
    fi
done

# FANTASIA-Lite stand-ins. Some scenarios (e.g. scenario_08_es) flip
# BRAKER4_RUN_FANTASIA=1 in their scenario_overrides.sh. The Snakefile then
# parse-time validates that the configured FANTASIA SIF (a file) and HF cache
# (a directory) actually exist on disk, refusing to build the DAG otherwise.
# On a workstation that has never pulled the FANTASIA image, the configured
# NAS path does not exist and the rulegraph generation fails. Provide stand-in
# paths under /tmp and override the env vars accordingly so the parse-time
# checks pass. Snakemake never reads the contents of either path during
# --rulegraph (it only walks the DAG), so empty files / dirs are sufficient.
FANTASIA_DUMMY_DIR=$(mktemp -d -t braker4-rulegraph-fantasia.XXXXXX)
FANTASIA_DUMMY_SIF="$FANTASIA_DUMMY_DIR/fantasia_lite.sif"
FANTASIA_DUMMY_HF="$FANTASIA_DUMMY_DIR/hf_cache"
touch "$FANTASIA_DUMMY_SIF"
mkdir -p "$FANTASIA_DUMMY_HF"
export BRAKER4_FANTASIA_SIF="$FANTASIA_DUMMY_SIF"
export BRAKER4_FANTASIA_HF_CACHE="$FANTASIA_DUMMY_HF"

# Scan both test_scenarios/ and test_scenarios_local/.
#
# Each tier (test_scenarios/ and test_scenarios_local/) has a SHARED biology
# config at $TIER_DIR/config.ini, picked up by the Snakefile via the
# BRAKER4_CONFIG env var. Per-scenario tweaks live in scenario_overrides.sh
# which exports BRAKER4_<KEY>=... env vars that the Snakefile reads on top
# of the shared config. Mirror that contract here so the rulegraph DAG
# matches what `run_test.sh` would actually build.
# Scenarios excluded by name. The athaliana benchmark is a real-data run that
# lives only on brain (the HPC) — its samples.csv references absolute NAS
# paths under /home/nas-hs/projs/braker4/data/Arabidopsis_thaliana/ that do
# not exist on a workstation, so even if we built the DAG it would point at
# files no local user can read. Keep it out of the rulegraph generator.
EXCLUDED_SCENARIOS=(scenario_benchmark_athaliana)

is_excluded() {
    local name=$1
    for excl in "${EXCLUDED_SCENARIOS[@]}"; do
        [ "$name" = "$excl" ] && return 0
    done
    return 1
}

for SCENARIO_DIR in "$PIPELINE_DIR"/test_scenarios/scenario_*/ "$PIPELINE_DIR"/test_scenarios_local/scenario_*/; do
    [ -d "$SCENARIO_DIR" ] || continue

    TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
    SCENARIO_NAME="$(basename "$SCENARIO_DIR")"
    SCENARIO="$(basename "$TIER_DIR")/$SCENARIO_NAME"

    if is_excluded "$SCENARIO_NAME"; then
        echo "SKIP $SCENARIO (excluded by name)"
        continue
    fi
    if [ ! -f "$SCENARIO_DIR/samples.csv" ]; then
        echo "SKIP $SCENARIO (no samples.csv)"
        continue
    fi
    if [ ! -f "$TIER_DIR/config.ini" ]; then
        echo "SKIP $SCENARIO (tier has no config.ini at $TIER_DIR/config.ini)"
        continue
    fi

    echo -n "$SCENARIO: "
    cd "$SCENARIO_DIR"

    # Run snakemake in a subshell so per-scenario BRAKER4_* exports do not
    # leak across iterations.
    DOT_OUTPUT=$(
        export BRAKER4_CONFIG="$TIER_DIR/config.ini"
        if [ -f "$SCENARIO_DIR/scenario_overrides.sh" ]; then
            # shellcheck disable=SC1090
            source "$SCENARIO_DIR/scenario_overrides.sh"
        fi
        snakemake --snakefile "$PIPELINE_DIR/Snakefile" --rulegraph 2>/dev/null \
            | sed -n '/^digraph/,/^}/p' || true
    )

    if [ -n "$DOT_OUTPUT" ]; then
        echo "$DOT_OUTPUT" | dot -Tpng > rulegraph.png
        echo "OK ($(du -h rulegraph.png | cut -f1))"
        PASSED=$((PASSED + 1))
    else
        echo "FAILED (re-run the snakemake call manually with stderr to see why)"
        FAILED=$((FAILED + 1))
    fi
done

# Clean up dummy files
for f in "${DUMMY_FILES[@]}"; do
    rm -f "$f"
done
rm -rf "$FANTASIA_DUMMY_DIR"

echo ""
echo "Done: $PASSED succeeded, $FAILED failed"
