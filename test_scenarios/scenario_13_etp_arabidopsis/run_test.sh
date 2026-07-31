#!/bin/bash
# Benchmark Scenario 13: ETP Mode (A. thaliana)
# BRAKER Mode: ETP (genome + proteins + RNA-seq BAM)
# Data: Full A. thaliana genome, OrthoDB proteins, VARUS RNA-seq BAM
# Note: genome_masked is empty — RepeatMasker runs as part of the pipeline.

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "============================================"
echo "Benchmark Scenario 13: ETP Mode (A. thaliana)"
echo "============================================"
echo ""
echo "Configuration:"
echo "  - Genome: A. thaliana (~135 MB, unmasked — masking runs in-pipeline)"
echo "  - Proteins: /home/nas-hs/projs/braker4/data/Arabidopsis_thaliana/proteins.fa"
echo "  - RNA-Seq: VARUS BAM"
echo "  - BRAKER Mode: ETP"
echo "  - skip_optimize_augustus: 1"
echo "  - run_omark: 0"
echo ""

# Ensure singularity is available regardless of which login node we land on.
export PATH=/opt/singularity-3.11.3/bin:${PATH}

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
    --singularity-args "-B /home -B /projects --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS

echo ""
[ "$DRY_RUN" = "true" ] && echo "Dry-run completed!" || echo "Benchmark completed!"
