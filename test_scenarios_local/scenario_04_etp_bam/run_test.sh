#!/bin/bash
# Test Scenario 04: BAM + Proteins + Pre-masked Genome (ETP mode)
# BRAKER Mode: ETP (Evidence from Transcripts and Proteins)
# Data: Pre-masked genome + BAM + Proteins
# This is STANDARD BRAKER3 with pre-masked genome!

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "============================================================="
echo "Test Scenario 04: BAM + Proteins + Pre-masked (ETP Mode)"
echo "============================================================="
echo ""
echo "Configuration:"
echo "  - Genome: test_data/genome.fasta"
echo "  - Masked Genome: test_data/genome_masked.fasta"
echo "  - Proteins: test_data/proteins.fasta"
echo "  - RNA-Seq: test_data/RNAseq.bam (pre-aligned)"
echo "  - BRAKER Mode: ETP (standard BRAKER3)"
echo "  - Masking: SKIP (using pre-masked genome)"
echo ""
echo "NOTE: Standard BRAKER3 workflow with time-saving pre-masked genome"
echo ""

# Download test BAM file if not present
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
    --nolock \
    --use-singularity \
    --singularity-prefix "$PIPELINE_DIR/.singularity_cache" \
    --singularity-args "-B /home --env PREPEND_PATH=/opt/conda/bin" \
    $EXECUTOR_ARGS

echo ""
[ "$DRY_RUN" = "true" ] && echo "✓ Dry-run completed!" || echo "✓ Test completed!"
