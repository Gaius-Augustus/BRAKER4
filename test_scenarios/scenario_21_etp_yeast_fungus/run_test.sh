#!/bin/bash
# Test Scenario 21: ETP mode, S. cerevisiae, fungus = 1 (issue #110)
set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "=============================================================="
echo "Test Scenario 21: ETP, S. cerevisiae, fungus = 1 (issue #110)"
echo "=============================================================="
echo ""
echo "Configuration:"
echo "  - Genome: S. cerevisiae S288C R64 (~12 Mbp, unmasked, masked by the pipeline)"
echo "  - Proteins: OrthoDB v12 Fungi"
echo "  - RNA-Seq: SRA SRR5819434, SRR5819435 (paired)"
echo "  - fungus = 1, GeneMark-ETP with 8 cores"
echo ""

bash "$SCENARIO_DIR/get_data.sh"
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
