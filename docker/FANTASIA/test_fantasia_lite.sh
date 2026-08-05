#!/bin/bash
#SBATCH --job-name=fantasia-lite-test
#SBATCH --output=fantasia-lite-test_%j.out
#SBATCH --error=fantasia-lite-test_%j.err
#SBATCH --gpus=1
#SBATCH --partition=vision,vision-fast
#SBATCH --cpus-per-task=30
#SBATCH --mem=150G
#SBATCH --time=2:00:00
#SBATCH --get-user-env

# ------------------------------------------------------------------
# Test script: pull the FANTASIA-Lite Singularity image and run it
# on the bundled test.fa with GPU support.
#
# Usage:
#   cd Nextflow_Tiberius_FANTASIA/dockerfiles/FANTASIA
#   sbatch test_fantasia_lite.sh
# ------------------------------------------------------------------

set -euo pipefail

echo "[$(date)] Loading modules"
module load singularity

WORK_DIR="$HOME/git/Nextflow_Tiberius_FANTASIA/dockerfiles/FANTASIA"
SIF="${WORK_DIR}/fantasia-lite.sif"

# Singularity needs a writable temp dir for SIF creation; the default
# SLURM tmpdir may not be writable.
export SINGULARITY_TMPDIR="${WORK_DIR}/tmp_singularity"
export TMPDIR="${SINGULARITY_TMPDIR}"
mkdir -p "$SINGULARITY_TMPDIR"

# --- Step 1: Pull the SIF from Docker Hub if it doesn't exist yet ---
if [ ! -f "$SIF" ]; then
    echo "[$(date)] Pulling fantasia-lite image from Docker Hub"
    singularity pull "$SIF" docker://katharinahoff/fantasia_for_brain:lite.v1.0.0
    rm -rf "$SINGULARITY_TMPDIR"
fi

echo "[$(date)] Using SIF: $SIF"

# --- Step 2: Run the smoke test with GPU ---------------------------
OUTDIR="${WORK_DIR}/test_output_$$"
mkdir -p "$OUTDIR"

echo "[$(date)] Running FANTASIA-Lite on test.fa (GPU mode)"
LOOKUP="${FANTASIA_LOOKUP_DIR:-/home/nas-hs/data/db/EukBinDB/fantasia_v1_lookup}"
singularity exec --nv \
    --bind "$OUTDIR":"$OUTDIR" \
    --bind "$LOOKUP":"$LOOKUP" \
    "$SIF" \
    python3 /opt/fantasia-lite/src/fantasia_pipeline.py \
        --serial-models \
        --embed-models prot_t5 \
        --device cuda \
        --venv-dir /opt/venv \
        --lookup-npz "$LOOKUP/lookup_table.npz" \
        --annotations-json "$LOOKUP/annotations.json" \
        --accessions-json "$LOOKUP/accessions.json" \
        --embeddings-npz "$OUTDIR/query_embeddings.npz" \
        --config-yaml "$OUTDIR/fantasia_config.yaml" \
        --results-csv "$OUTDIR/results.csv" \
        --topgo \
        --topgo-dir "$OUTDIR/topgo" \
        --chunk-dir "$OUTDIR/tmp/fasta_chunks" \
        --chunk-embed-dir "$OUTDIR/tmp/chunk_embeddings" \
        --chunk-results-dir "$OUTDIR/tmp/chunk_results" \
        --chunk-config-dir "$OUTDIR/tmp/chunk_configs" \
        --chunk-failure-dir "$OUTDIR/tmp/failures" \
        --failure-report "$OUTDIR/failed_sequences.csv" \
        /opt/fantasia-lite/test.fa

echo ""
echo "[$(date)] DONE. Check results in: $OUTDIR"
ls -lh "$OUTDIR"
