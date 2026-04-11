#!/bin/bash
# Master test runner for all BRAKER HPC test scenarios (O. tauri genome)
#
# Each scenario's Snakemake runs under nohup so it survives shell/screen issues.
# A PID tracking file lets you monitor and kill processes if needed.
#
# Usage:
#   bash run_all_tests.sh                    # Run all scenarios
#   bash run_all_tests.sh 08 11              # Run specific scenarios by number
#   DRY_RUN=true bash run_all_tests.sh       # Dry-run only
#   bash run_all_tests.sh status             # Show status of running tests
#   bash run_all_tests.sh kill               # Kill all running test processes

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Configuration
# Compute defaults (cores, partitions, mem, runtime) live in compute_profile.sh
# in this directory. Edit that file to change the cluster configuration.
# You can still override the partition list ad-hoc by exporting PARTITION
# (singular) — it propagates through compute_profile.sh into each scenario.
DRY_RUN=${DRY_RUN:-false}
RESULTS_DIR="$SCRIPT_DIR/test_results_$(date +%Y%m%d_%H%M%S)"
PIDFILE="$SCRIPT_DIR/.running_tests.pid"

ALL_SCENARIOS=(
    scenario_08_es
    scenario_09_ep_masking
    scenario_10_et_varus
    scenario_11_etp_sra
    scenario_12_multi_mode
)

# Default to launching all scenarios in parallel — no need to throttle since
# Snakemake itself manages SLURM job concurrency. Override with MAX_CONCURRENT env.
MAX_CONCURRENT=${MAX_CONCURRENT:-${#ALL_SCENARIOS[@]}}

# Handle special commands
if [ "${1:-}" = "status" ]; then
    echo "=== Running BRAKER test processes ==="
    if [ -f "$PIDFILE" ]; then
        while IFS=$'\t' read -r _pid _scen _log; do
            if kill -0 "$_pid" 2>/dev/null; then
                steps=$(grep -c "steps.*done" "$_log" 2>/dev/null | tail -1 || echo "?")
                echo "  RUNNING  PID=$_pid  $_scen  (last: $steps steps)"
            else
                if grep -q "Test completed" "$_log" 2>/dev/null; then
                    echo "  DONE     PID=$_pid  $_scen"
                else
                    echo "  DEAD     PID=$_pid  $_scen"
                fi
            fi
        done < "$PIDFILE"
    else
        echo "  No PID file found. No tests running?"
    fi
    exit 0
fi

is_official_scenario_dir() {
    # Return 0 if $1 is exactly $SCRIPT_DIR/<official_scenario>, else 1.
    # We compare against the explicit ALL_SCENARIOS allow-list rather than
    # globbing scenario_* so that user benchmark directories like
    # scenario_benchmark_athaliana are NOT matched.
    local dir="$1"
    local s
    for s in "${ALL_SCENARIOS[@]}"; do
        [ "$dir" = "$SCRIPT_DIR/$s" ] && return 0
    done
    return 1
}

cancel_slurm_jobs_from_scenarios() {
    # Cancel only SLURM jobs whose WorkDir is exactly one of the official
    # scenario dirs in ALL_SCENARIOS. Jobs from other pipelines (EukAssembly,
    # etc.) and from user benchmarks (scenario_benchmark_*) are left alone.
    command -v squeue >/dev/null 2>&1 || return 0
    local jobids
    jobids=$(squeue -u "$USER" -h -o '%i' 2>/dev/null) || return 0
    [ -z "$jobids" ] && return 0
    local to_cancel=""
    for jid in $jobids; do
        local wd
        wd=$(scontrol show job "$jid" 2>/dev/null | grep -oP 'WorkDir=\K\S+') || continue
        if is_official_scenario_dir "$wd"; then
            to_cancel="$to_cancel $jid"
        fi
    done
    if [ -n "$to_cancel" ]; then
        echo "  Cancelling SLURM jobs from this runner's scenarios:$to_cancel"
        scancel $to_cancel 2>/dev/null || true
    fi
}

kill_snakemake_from_scenarios() {
    # Kill only snakemake processes whose cwd is exactly one of the official
    # scenario dirs. User benchmark scenarios are left alone.
    for pid in $(pgrep -u "$USER" -f snakemake 2>/dev/null); do
        local pwd_link
        pwd_link=$(readlink "/proc/$pid/cwd" 2>/dev/null) || continue
        if is_official_scenario_dir "$pwd_link"; then
            kill "$pid" 2>/dev/null || true
        fi
    done
}

if [ "${1:-}" = "kill" ]; then
    echo "=== Killing all BRAKER test processes ==="
    # Kill tracked PIDs
    if [ -f "$PIDFILE" ]; then
        while IFS=$'\t' read -r _pid _scen _log; do
            if kill -0 "$_pid" 2>/dev/null; then
                echo "  Killing PID=$_pid ($_scen)..."
                kill -- -"$_pid" 2>/dev/null || kill "$_pid" 2>/dev/null || true
            fi
        done < "$PIDFILE"
        rm -f "$PIDFILE"
    fi
    # Kill any orphaned snakemake processes whose cwd is within our scenario dirs
    kill_snakemake_from_scenarios
    # Cancel only SLURM jobs from our scenarios (NOT scancel -u $USER!)
    cancel_slurm_jobs_from_scenarios
    echo "  Done. BRAKER test scenarios stopped (other pipelines untouched)."
    exit 0
fi

# Select scenarios: from args or all
if [ $# -gt 0 ]; then
    SCENARIOS=()
    for num in "$@"; do
        for s in "${ALL_SCENARIOS[@]}"; do
            if [[ "$s" == scenario_${num}_* ]]; then
                SCENARIOS+=("$s")
                break
            fi
        done
    done
else
    SCENARIOS=("${ALL_SCENARIOS[@]}")
fi

echo "============================================================"
echo "  BRAKER Pipeline - HPC Test Runner (O. tauri)"
echo "============================================================"
echo ""
echo "  Scenarios:       ${#SCENARIOS[@]}"
echo "  Max concurrent:  $MAX_CONCURRENT"
echo "  SLURM profile:   $SCRIPT_DIR/compute_profile.sh"
echo "  Dry-run:         $DRY_RUN"
echo "  Results dir:     $RESULTS_DIR"
echo ""

mkdir -p "$RESULTS_DIR"

# Download test data first
echo "Downloading test data if needed..."
if [ -f "$PIPELINE_DIR/test_data/download_test_data.sh" ]; then
    bash "$PIPELINE_DIR/test_data/download_test_data.sh" 2>/dev/null || true
fi
echo ""

# Clean up orphaned snakemake processes + SLURM jobs from previous BRAKER test runs.
# Only touches processes/jobs whose cwd/WorkDir is within our scenario dirs —
# other pipelines (EukAssembly, etc.) are left alone.
kill_snakemake_from_scenarios
cancel_slurm_jobs_from_scenarios
sleep 1

# Clean up stale SLURM log files older than 10 days, but only inside the
# official scenario directories. User benchmarks (e.g. scenario_benchmark_*)
# are left alone so their logs survive across runs.
echo "Cleaning up SLURM log files older than 10 day(s)."
for s in "${ALL_SCENARIOS[@]}"; do
    if [ -d "$SCRIPT_DIR/$s/.snakemake/slurm_logs" ]; then
        find "$SCRIPT_DIR/$s/.snakemake/slurm_logs" -type f -name "*.log" -mtime +10 -delete 2>/dev/null || true
    fi
done

# Initialize PID tracking
> "$PIDFILE"

launch_scenario() {
    local scenario="$1"
    local scenario_dir="$SCRIPT_DIR/$scenario"
    local logfile="$RESULTS_DIR/${scenario}.log"

    if [ ! -f "$scenario_dir/run_test.sh" ]; then
        echo "  SKIP $scenario (no run_test.sh)"
        return
    fi

    # Clean stale state
    rm -rf "$scenario_dir/.snakemake/incomplete/" "$scenario_dir/.snakemake/locks/" 2>/dev/null || true

    # Launch under nohup — survives shell/screen issues.
    # PARTITION / CORES / etc. come from compute_profile.sh, sourced inside
    # run_test.sh; DRY_RUN propagates from this runner's environment.
    cd "$scenario_dir"
    nohup bash -c "DRY_RUN=$DRY_RUN bash run_test.sh" \
        > "$logfile" 2>&1 &
    local pid=$!
    cd "$SCRIPT_DIR"

    # Record PID
    echo -e "${pid}\t${scenario}\t${logfile}" >> "$PIDFILE"
    echo "  LAUNCHED  PID=$pid  $scenario  -> $logfile"
}

wait_for_slots() {
    # Wait until fewer than MAX_CONCURRENT processes are running
    while true; do
        local running=0
        while IFS=$'\t' read -r _pid _scen _log; do
            if kill -0 "$_pid" 2>/dev/null; then
                running=$((running + 1))
            fi
        done < "$PIDFILE"
        if [ "$running" -lt "$MAX_CONCURRENT" ]; then
            break
        fi
        sleep 30
    done
}

# === Launch all scenarios ===
echo "=== Launching ${#SCENARIOS[@]} scenarios ==="
for scenario in "${SCENARIOS[@]}"; do
    wait_for_slots
    launch_scenario "$scenario"
    sleep 5  # Brief stagger
done
echo ""
echo "All scenarios launched. Waiting for completion..."
echo ""

# Wait for all to finish
while true; do
    running=0
    while IFS=$'\t' read -r _pid _scen _log; do
        if kill -0 "$_pid" 2>/dev/null; then
            running=$((running + 1))
        fi
    done < "$PIDFILE"
    if [ "$running" -eq 0 ]; then
        break
    fi
    echo "  $(date +%H:%M:%S) — $running scenarios still running..."
    sleep 60
done
echo "All scenarios complete."
echo ""

# === Results Summary ===
echo "============================================================"
echo "  Results Summary"
echo "============================================================"
echo ""

PASS=0; FAIL=0; SKIP=0
for scenario in "${SCENARIOS[@]}"; do
    logfile="$RESULTS_DIR/${scenario}.log"
    if [ ! -f "$logfile" ]; then
        echo "  ?    $scenario (no log)"
        continue
    fi

    sample_name=$(awk -F',' 'NR==2{print $1}' "$SCRIPT_DIR/$scenario/samples.csv" 2>/dev/null)
    gtf="$SCRIPT_DIR/$scenario/output/${sample_name}/braker.gtf"

    if grep -q "SKIP" "$logfile" 2>/dev/null; then
        echo "  -    $scenario (SKIPPED)"
        SKIP=$((SKIP + 1))
    elif [ "$DRY_RUN" = "true" ] && ! grep -qi "error\|failed" "$logfile" 2>/dev/null; then
        echo "  PASS $scenario (dry-run OK)"
        PASS=$((PASS + 1))
    elif [ -f "$gtf" ] && [ -s "$gtf" ]; then
        genes=$(grep -cP '\tgene\t' "$gtf" 2>/dev/null || echo 0)
        echo "  PASS $scenario ($genes genes)"
        PASS=$((PASS + 1))
    else
        last_error=$(grep -i "error\|failed" "$logfile" 2>/dev/null | tail -1 | head -c 120)
        echo "  FAIL $scenario: $last_error"
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "  PASS: $PASS  FAIL: $FAIL  SKIP: $SKIP  Total: ${#SCENARIOS[@]}"
echo ""
echo "  Logs: $RESULTS_DIR/"
echo "  PIDs: $PIDFILE"
echo "============================================================"

rm -f "$PIDFILE"
[ "$FAIL" -eq 0 ] && exit 0 || exit 1
