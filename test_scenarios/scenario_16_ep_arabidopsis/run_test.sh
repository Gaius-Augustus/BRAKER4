#!/bin/bash
# Benchmark Scenario 16: EP Mode (A. thaliana, BRAKER2-style)
# BRAKER Mode: EP (genome + proteins, no RNA-seq)
# Data: Full A. thaliana genome (pre-masked), OrthoDB proteins

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "============================================"
echo "Benchmark Scenario 16: EP Mode (A. thaliana)"
echo "============================================"
echo ""
echo "Configuration:"
echo "  - Genome: A. thaliana (~135 MB, pre-masked)"
echo "  - Proteins: OrthoDB proteins"
echo "  - RNA-Seq: none (EP / BRAKER2-style)"
echo "  - BRAKER Mode: EP"
echo "  - skip_optimize_augustus: 1"
echo "  - run_omark: 0"
echo ""

export PATH=/opt/singularity-3.11.3/bin:${PATH}

source "$TIER_DIR/compute_profile.sh"
[ -f "$SCENARIO_DIR/scenario_overrides.sh" ] && source "$SCENARIO_DIR/scenario_overrides.sh"

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
