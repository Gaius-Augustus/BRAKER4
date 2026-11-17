#!/bin/bash
#SBATCH --job-name=braker4_s16
#SBATCH --output=snakemake_mother.log
#SBATCH --error=snakemake_mother.err
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2880

# Snakemake mother process of this scenario as a SLURM job; the rules are
# submitted as jobs of their own (compute_profile.sh).
# Submit from this scenario directory:  sbatch submit_snakemake.sh
# Tool locations of your cluster (snakemake, singularity) belong in
# ../site_profile.sh, which compute_profile.sh sources.

source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity 2>/dev/null || true

cd "${SLURM_SUBMIT_DIR:?submit this script with sbatch from the scenario directory}"
if [ ! -f run_test.sh ]; then
    echo "ERROR: run_test.sh not found in $PWD; submit from the scenario directory" >&2
    exit 1
fi

export DRY_RUN=false
bash run_test.sh
