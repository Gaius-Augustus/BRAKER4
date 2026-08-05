#!/bin/bash
#SBATCH --job-name=braker3_etp_ath
#SBATCH --output=/home/hoffk83/git/BRAKER4/test_scenarios/scenario_19_braker3_etp_arabidopsis/braker3_etp.log
#SBATCH --error=/home/hoffk83/git/BRAKER4/test_scenarios/scenario_19_braker3_etp_arabidopsis/braker3_etp.err
#SBATCH --partition=batch,snowball,pinky
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=4320

set -euo pipefail

SCENARIO_DIR=/home/hoffk83/git/BRAKER4/test_scenarios/scenario_19_braker3_etp_arabidopsis
SIF_CACHE=/home/hoffk83/git/BRAKER4/.singularity_cache
BRAKER3_SIF=$SIF_CACHE/0ee9893ca6784f3cd184129689650c92.simg
GFFCMP_SIF=$SIF_CACHE/2b70909c3b56a6698bf067440b597123.simg

GENOME=/projects/AI-GUSTUS/tiberius_mesangiospermae_training_data/test/genomes/AthalianacolumbiaAraport11.genome.fa
BAM=/home/hoffk83/git/BRAKER4/test_scenarios/scenario_13_etp_arabidopsis/output/ath_etp/varus/varus.sorted.bam
PROTEINS=/home/hoffk83/GALBA2_test_data/proteins.fa
REFERENCE=/projects/AI-GUSTUS/tiberius_mesangiospermae_training_data/test/annotations/AthalianacolumbiaAraport11.gff3
WORKDIR=$SCENARIO_DIR/output
GFFCMP_DIR=$SCENARIO_DIR/gffcompare
THREADS=48

source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity
export PATH=/opt/singularity-3.11.3/bin:${PATH}

if [ ! -f "$BAM" ]; then
    echo "[ERROR] VARUS BAM not found: $BAM"
    exit 1
fi

mkdir -p "$WORKDIR" "$GFFCMP_DIR"

# Copy Augustus config from container (braker.pl requires a writable config)
AUGUSTUS_CONFIG=$SCENARIO_DIR/augustus_config
if [ ! -d "$AUGUSTUS_CONFIG" ]; then
    echo "[INFO] Copying Augustus config from container..."
    singularity exec \
        -B /home -B /projects \
        "$BRAKER3_SIF" \
        cp -r /opt/Augustus/config "$AUGUSTUS_CONFIG"
    chmod -R u+w "$AUGUSTUS_CONFIG"
    echo "[INFO] Augustus config ready."
fi

echo "[INFO] Running braker.pl ETP mode (BAM + proteins, skipOptimize)..."
singularity exec \
    -B /home -B /projects \
    -B "${AUGUSTUS_CONFIG}:/opt/Augustus/config" \
    "$BRAKER3_SIF" \
    /opt/BRAKER/scripts/braker.pl \
        --genome="$GENOME" \
        --bam="$BAM" \
        --prot_seq="$PROTEINS" \
        --species=ath_braker3_etp \
        --workingdir="$WORKDIR" \
        --threads=$THREADS \
        --skipOptimize \
        --etpmode

echo "[INFO] braker.pl finished. Evaluating with gffcompare..."

awk '$3 == "CDS"' "$REFERENCE"          > "$GFFCMP_DIR/reference.CDS.gtf"
awk '$3 == "CDS"' "$WORKDIR/braker.gtf" > "$GFFCMP_DIR/prediction.CDS.gtf"
echo "Reference CDS lines : $(wc -l < $GFFCMP_DIR/reference.CDS.gtf)"
echo "Prediction CDS lines: $(wc -l < $GFFCMP_DIR/prediction.CDS.gtf)"

singularity exec \
    -B /home -B /projects \
    "$GFFCMP_SIF" \
    gffcompare \
        --strict-match -e 3 -T \
        -r "$GFFCMP_DIR/reference.CDS.gtf" \
        -o "$GFFCMP_DIR/gffcompare" \
        "$GFFCMP_DIR/prediction.CDS.gtf"

echo ""
echo "=== gffcompare results ==="
cat "$GFFCMP_DIR/gffcompare.stats"
echo "[INFO] Done."
