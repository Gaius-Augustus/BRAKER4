#!/bin/bash
# Test Scenario 10: ET mode via VARUS (O. tauri)
# BRAKER Mode: ET (RNA-Seq only via VARUS)
# Data: O. tauri genome (pre-masked) + VARUS
#
# VARUS auto-discovers and downloads RNA-Seq data for O. tauri from SRA.
# The genome is used as-is for the "masked" genome (skipping RepeatMasker).

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "=========================================================="
echo "Test Scenario 10: ET Mode via VARUS (O. tauri)"
echo "=========================================================="
echo ""
echo "Configuration:"
echo "  - Genome: O. tauri (~12.6 MB, pre-masked)"
echo "  - Masked Genome: same (skipping RepeatMasker)"
echo "  - RNA-Seq: VARUS (Ostreococcus tauri)"
echo "  - BRAKER Mode: ET (RNA-Seq only, no proteins)"
echo "  - Masking: SKIP (using genome as pre-masked)"
echo ""
echo "NOTE: Requires internet access (VARUS downloads RNA-Seq from SRA)"
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
