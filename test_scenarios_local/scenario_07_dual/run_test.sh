#!/bin/bash
# Test Scenario 07: BAM + IsoSeq + Proteins + Pre-masked (Dual-Evidence ETP)
# BRAKER Mode: ETP with dual transcriptomics evidence
# Data: Genome + short-read BAM + IsoSeq BAM + Proteins
# Pipeline: Two separate GeneMark-ETP runs (short-read + isoseq),
#           merge HC genes, single AUGUSTUS training + prediction

set -euo pipefail

SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIER_DIR="$(cd "$SCENARIO_DIR/.." && pwd)"
PIPELINE_DIR="$(cd "$TIER_DIR/.." && pwd)"

echo "================================================================="
echo "Test Scenario 07: BAM + IsoSeq + Proteins + Pre-masked (Dual-Evidence ETP)"
echo "================================================================="
echo ""
echo "Configuration:"
echo "  - Genome: test_data/genome.fa"
echo "  - Proteins: test_data/proteins.fa"
echo "  - RNA-Seq: test_data/RNAseq.bam (short-read, HISAT2-aligned)"
echo "  - IsoSeq: test_data/isoseq.bam (long-read, minimap2-aligned)"
echo "  - BRAKER Mode: dual ETP (separate GeneMark runs, merged training)"
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
[ "$DRY_RUN" = "true" ] && echo "Dry-run completed!" || echo "Test completed!"
