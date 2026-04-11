#!/bin/bash
# Test Scenario 12: Multi-Mode (all 7 BRAKER modes in one samples.csv)
# Runs ES, EP, ET, ETP, IsoSeq-BAM, IsoSeq-FASTA, and Dual modes
# on the small A. thaliana test genome in a single Snakemake invocation.
# This tests multi-sample support and all pipeline branches simultaneously.

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "============================================"
echo "Test Scenario 12: Multi-Mode (7 samples)"
echo "============================================"
echo ""
echo "Configuration:"
echo "  - Genome: A. thaliana chr5 fragment (~1 MB, pre-masked)"
echo "  - Samples: ES, EP, ET(FASTQ), ETP(BAM), IsoSeq-BAM, IsoSeq-FASTA, Dual"
echo ""

# Download test data if not present
source "$PIPELINE_DIR/test_data/download_test_data.sh"
echo ""

# Load shared compute profile, then per-scenario overrides
source "$TIER_DIR/compute_profile.sh"
[ -f "$SCENARIO_DIR/scenario_overrides.sh" ] && source "$SCENARIO_DIR/scenario_overrides.sh"

# Use the shared biology config.ini for this tier
export BRAKER4_CONFIG="$TIER_DIR/config.ini"

DRY_RUN=${DRY_RUN:-false}
if [ "$DRY_RUN" = "true" ]; then
    echo "Running in DRY-RUN mode"
    DRY_RUN_FLAG="-n"
else
    echo "Running FULL EXECUTION"
    DRY_RUN_FLAG=""
fi

cd "$SCENARIO_DIR"

export SINGULARITYENV_PREPEND_PATH=/opt/conda/bin
snakemake \
    --snakefile "$PIPELINE_DIR/Snakefile" \
    --cores "$CORES" --jobs "$CORES" \
    $DRY_RUN_FLAG \
    --printshellcmds \
    --rerun-incomplete \
    --use-singularity \
    --singularity-prefix "$PIPELINE_DIR/.singularity_cache" \
    --singularity-args "-B /home --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS

echo ""
[ "$DRY_RUN" = "true" ] && echo "✓ Dry-run completed!" || echo "✓ Test completed!"
