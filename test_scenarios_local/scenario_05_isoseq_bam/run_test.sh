#!/bin/bash
# Test Scenario 05: IsoSeq + Proteins + Pre-masked Genome (IsoSeq mode)
# BRAKER Mode: IsoSeq (Long-read transcripts + proteins)
# Data: Pre-masked genome + IsoSeq BAM + Proteins
set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "=================================================================="
echo "Test Scenario 05: IsoSeq + Proteins + Pre-masked (IsoSeq Mode)"
echo "=================================================================="
echo ""
echo "Configuration:"
echo "  - Genome: test_data/genome.fasta"
echo "  - Masked Genome: test_data/genome_masked.fasta"
echo "  - Proteins: test_data/proteins.fasta"
echo "  - IsoSeq: test_data/isoseq.bam (PacBio long-read)"
echo "  - BRAKER Mode: IsoSeq (long-read transcripts + proteins)"
echo "  - Masking: SKIP (using pre-masked genome)"
echo ""
echo "IMPORTANT: IsoSeq ALWAYS requires protein evidence!"
echo "           IsoSeq currently only supports BAM input"
echo "           Uses special IsoSeq container image"
echo ""
echo "NOTE: Demonstrates IsoSeq mode with time-saving pre-masked genome"
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
    --nolock \
    --use-singularity \
    --singularity-prefix "$PIPELINE_DIR/.singularity_cache" \
    --singularity-args "-B /home --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS

echo ""
[ "$DRY_RUN" = "true" ] && echo "✓ Dry-run completed!" || echo "✓ Test completed!"
