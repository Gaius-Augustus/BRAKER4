#!/bin/bash
# Shared compute profile for HPC test scenarios.
#
# This file is the SINGLE SOURCE OF TRUTH for the SLURM/cores configuration
# used by every scenario in test_scenarios/. Edit this file once to point the
# tests at your own cluster — there is no need to touch the individual
# scenario_*/run_test.sh files.
#
# Every variable below uses the ${VAR:-default} idiom so that an external
# environment variable still wins. That means a single test invocation can
# still override anything ad-hoc, e.g.:
#
#   PARTITION=highmem CORES=64 bash scenario_08_es/run_test.sh
#
# Per-scenario overrides (e.g. scenario_09 needs a longer max_runtime than the
# default) live in each scenario's optional scenario_overrides.sh, which is
# sourced AFTER this file by run_test.sh.

# How many cores snakemake itself uses (snakemake --cores / --jobs).
export CORES="${CORES:-48}"

# SLURM partition list (comma-separated). Adjust to match your cluster.
export PARTITION="${PARTITION:-batch,snowball,pinky}"

# Whether to submit each rule as a SLURM job (true) or run all rules locally
# in the snakemake process (false). Set to "false" on machines without SLURM.
export USE_SLURM="${USE_SLURM:-true}"

# Per-job resource requests passed to the Snakefile via env var overrides.
# These are picked up by Snakefile's BRAKER4_* lookup table and override
# whatever is in the shared config.ini.
export BRAKER4_CPUS_PER_TASK="${BRAKER4_CPUS_PER_TASK:-48}"
export BRAKER4_MEM_OF_NODE="${BRAKER4_MEM_OF_NODE:-120000}"
export BRAKER4_MAX_RUNTIME="${BRAKER4_MAX_RUNTIME:-120}"

# Maximum mem_mb that snakemake will request as a default-resource per job.
# This bounds the SLURM --mem request when a rule does not specify mem_mb.
export DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-120000}"

# Build the snakemake --executor / --default-resources arguments.
if [ "$USE_SLURM" = "true" ]; then
    export EXECUTOR_ARGS="--executor slurm --default-resources slurm_partition=$PARTITION mem_mb=$DEFAULT_MEM_MB"
else
    export EXECUTOR_ARGS=""
fi
