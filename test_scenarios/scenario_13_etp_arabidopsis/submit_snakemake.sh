#!/bin/bash
#SBATCH --job-name=braker4_s13
#SBATCH --output=/home/hoffk83/git/BRAKER4/test_scenarios/scenario_13_etp_arabidopsis/snakemake_mother.log
#SBATCH --error=/home/hoffk83/git/BRAKER4/test_scenarios/scenario_13_etp_arabidopsis/snakemake_mother.err
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2880

source /etc/profile.d/modules.sh 2>/dev/null || true
module load singularity

export PATH=/home/hoffk83/miniconda3/bin:/opt/singularity-3.11.3/bin:${PATH}

cd /home/hoffk83/git/BRAKER4/test_scenarios/scenario_13_etp_arabidopsis

export DRY_RUN=false
bash /home/hoffk83/git/BRAKER4/test_scenarios/scenario_13_etp_arabidopsis/run_test.sh
