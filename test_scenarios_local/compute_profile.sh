#!/bin/bash
# Shared compute profile for local (no-SLURM) test scenarios.
#
# This file is the SINGLE SOURCE OF TRUTH for the cores/memory configuration
# used by every scenario in test_scenarios_local/. Edit this file once to
# scale the local tests up or down — there is no need to touch the individual
# scenario_*/run_test.sh files.
#
# Every variable below uses the ${VAR:-default} idiom so that an external
# environment variable still wins. That means a single test invocation can
# still override anything ad-hoc, e.g.:
#
#   CORES=4 bash scenario_01_es/run_test.sh
#
# Per-scenario overrides (e.g. scenarios that need no_cleanup=1 for debugging)
# live in each scenario's optional scenario_overrides.sh, which is sourced
# AFTER this file by run_test.sh.

# How many cores snakemake itself uses (snakemake --cores / --jobs).
export CORES="${CORES:-8}"

# SLURM partition list — only consulted if USE_SLURM=true.
export PARTITION="${PARTITION:-batch}"

# Local mode: do NOT submit jobs to SLURM. Set to "true" if you want to run
# the local scenarios on a SLURM cluster anyway.
export USE_SLURM="${USE_SLURM:-false}"

# Per-job resource requests. These are still consulted by the Snakefile even
# in local mode (e.g. for the runtime hint on long-running rules).
export BRAKER4_CPUS_PER_TASK="${BRAKER4_CPUS_PER_TASK:-8}"
export BRAKER4_MEM_OF_NODE="${BRAKER4_MEM_OF_NODE:-32000}"
export BRAKER4_MAX_RUNTIME="${BRAKER4_MAX_RUNTIME:-1440}"

# Maximum mem_mb that snakemake will request as a default-resource per job.
export DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-32000}"

# Build the snakemake --executor / --default-resources arguments.
if [ "$USE_SLURM" = "true" ]; then
    export EXECUTOR_ARGS="--executor slurm --default-resources slurm_partition=$PARTITION mem_mb=$DEFAULT_MEM_MB"
else
    export EXECUTOR_ARGS=""
fi
