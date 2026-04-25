#!/bin/bash
# Test Scenario 09: EP mode + masking with O. tauri genome
# BRAKER Mode: EP (Evidence from Proteins)
# Data: O. tauri genome (~12.6 MB) + Viridiplantae proteins
# Tests: RepeatModeler2 + RepeatMasker + full EP pipeline

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "=========================================="
echo "Test Scenario 09: EP Mode + Masking (O. tauri)"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  - Genome: O. tauri (~12.6 MB, needs masking)"
echo "  - Proteins: Viridiplantae ODB12"
echo "  - RNA-Seq: None"
echo "  - BRAKER Mode: EP (proteins only)"
echo "  - Masking: RepeatModeler2 + RepeatMasker"
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
    --latency-wait 120 \
    --restart-times 3 \
    --use-singularity \
    --singularity-prefix "$PIPELINE_DIR/.singularity_cache" \
    --singularity-args "-B /home --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS

echo ""
[ "$DRY_RUN" = "true" ] && echo "✓ Dry-run completed!" || echo "✓ Test completed!"
