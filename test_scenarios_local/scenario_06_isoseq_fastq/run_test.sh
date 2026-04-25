#!/bin/bash
# Test Scenario 06: IsoSeq FASTQ + Proteins (Pre-masked)
# BRAKER Mode: IsoSeq + Proteins (with minimap2 alignment)
# Data: Genome + Pre-masked Genome + IsoSeq FASTQ + Proteins

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "============================================"
echo "Test Scenario 06: IsoSeq FASTQ + Proteins (Pre-masked)"
echo "============================================"
echo ""
echo "Configuration:"
echo "  - Genome: test_data/genome.fa"
echo "  - Masked genome: test_data/genome.fa (pre-masked)"
echo "  - Proteins: test_data/proteins.fa"
echo "  - IsoSeq: test_data/isoseq.fastq.gz (unaligned, needs minimap2)"
echo "  - BRAKER Mode: IsoSeq + Proteins"
echo ""
echo "NOTE: Requires IsoSeq FASTQ test data in test_data/"
echo "      minimap2 aligns IsoSeq reads to genome before BRAKER"
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
