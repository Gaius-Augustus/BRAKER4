#!/bin/bash
# Start scenario 21 locally (no SLURM) only once the machine is idle:
# the 15-min load average must stay below IDLE_LOAD for IDLE_CHECKS
# consecutive checks, CHECK_INTERVAL seconds apart.
set -euo pipefail
SCENARIO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IDLE_LOAD="${IDLE_LOAD:-1.5}"
IDLE_CHECKS="${IDLE_CHECKS:-3}"
CHECK_INTERVAL="${CHECK_INTERVAL:-300}"

ok=0
while [ "$ok" -lt "$IDLE_CHECKS" ]; do
    load15=$(awk '{print $3}' /proc/loadavg)
    if awk -v l="$load15" -v t="$IDLE_LOAD" 'BEGIN{exit !(l < t)}'; then
        ok=$((ok + 1))
    else
        ok=0
    fi
    echo "[$(date '+%F %T')] load15=$load15 idle_checks=$ok/$IDLE_CHECKS"
    [ "$ok" -ge "$IDLE_CHECKS" ] && break
    sleep "$CHECK_INTERVAL"
done

echo "[$(date '+%F %T')] machine idle, starting scenario 21"
export USE_SLURM=false CORES=8 BRAKER4_CPUS_PER_TASK=8
export BRAKER4_MEM_OF_NODE=28000 DEFAULT_MEM_MB=28000
# Local snakemake lives in the pixi env (pixi binary from the miniforge3 env "pixi")
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate pixi
pixi run --manifest-path "$SCENARIO_DIR/../../pixi.toml" bash "$SCENARIO_DIR/run_test.sh"
