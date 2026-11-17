# Sourced by run_genemark_et and run_genemark_ep (automatic mode fallback,
# [PARAMS] auto_mode_fallback).
#
# genemark_fallback_marker TOOL EXIT_CODE MARKER LOG...
#
# Called when gmes_petap.pl ended without genemark.gtf. Decides whether the
# sample may continue in ES mode and, if so, writes the reason to MARKER
# (read by check_et_result / check_ep_result in check_evidence.smk) and
# returns 0. Returns 1 when the failure must stay a failure:
#   - the program was killed by a signal (exit code >= 128: out of memory
#     kill, time limit, segmentation fault), or
#   - a log shows a full disk, exhausted memory, or the crash or kill of a
#     program that GeneMark called (segmentation fault, "error code 11").
# gmes_petap.pl has no error message of its own for too little evidence; it
# stops somewhere in training (division by zero in parse_ET.pl, missing
# model files, no valid sequences). Every error exit of the program itself
# therefore counts as "the hints were not enough to train on". In EP mode
# GeneMark-ES has already run on the same genome in the same environment,
# so the hints are the only difference; in ET mode a technical problem shows
# up again in the GeneMark-ES run that follows.
genemark_fallback_marker() {
    local tool=$1 code=$2 marker=$3
    shift 3
    if [ "$code" -ge 128 ]; then
        return 1
    fi
    if grep -qiE 'No space left on device|Disk quota exceeded|Cannot allocate memory|Out of memory|Segmentation fault|core dumped|error code (9|11|134|137|139)([^0-9]|$)' "$@" 2>/dev/null; then
        return 1
    fi
    local last
    # no match is not an error (the caller may run under set -e -o pipefail)
    last=$(grep -hiE 'error|illegal' "$@" 2>/dev/null | grep -v '^ERROR: GeneMark' | tail -n 1 | tr '\t' ' ' | cut -c1-200 || true)
    mkdir -p "$(dirname "$marker")"
    {
        printf '%s stopped with an error during training (exit code %s) and wrote no genemark.gtf' "$tool" "$code"
        if [ -n "$last" ]; then
            printf '; last error message: %s' "$last"
        fi
        printf '\n'
    } > "$marker"
    return 0
}
