#!/bin/bash
# Per-scenario overrides for scenario_21_etp_yeast_fungus (issue #110).
# fungus = 1 and 8 cores, as in the report (gmes_petap.pl ran with --cores 8).
export BRAKER4_FUNGUS=1
export CORES="${CORES:-48}"
export PARTITION="${PARTITION:-batch,snowball,pinky}"
export USE_SLURM="${USE_SLURM:-true}"
export BRAKER4_CPUS_PER_TASK="${BRAKER4_CPUS_PER_TASK:-8}"
export BRAKER4_MEM_OF_NODE="${BRAKER4_MEM_OF_NODE:-64000}"
export BRAKER4_MAX_RUNTIME="${BRAKER4_MAX_RUNTIME:-2880}"
export DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-64000}"
if [ "$USE_SLURM" = "true" ]; then
    export EXECUTOR_ARGS="--executor slurm --default-resources slurm_partition=$PARTITION mem_mb=$DEFAULT_MEM_MB"
else
    export EXECUTOR_ARGS=""
fi
