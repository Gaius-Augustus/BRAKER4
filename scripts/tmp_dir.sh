# shellcheck shell=bash
# Scratch directories for rules that write many or large temporary files.
# Source it in a rule's shell block:  source scripts/tmp_dir.sh
#
# Network file systems with slow synchronous writes (brain's CephFS: ~0.4 s
# per file creation on bad days, 2026-10-06) turn tools that write hundreds
# of small files into hour-long jobs.  Such rules work in a private
# directory on the node-local disk and copy only the kept results back to
# the run directory (README.md, "HPC scratch").
#
# make_tmp PREFIX [ROOT] [NEED_GB]
#   Create a private directory PREFIX_XXXXXX under ROOT and print its path.
#   ROOT is [paths] tmp_dir (TMP_ROOT in the Snakefile, which also binds
#   it into containers); empty -> $TMPDIR -> /tmp.  Fails when ROOT is not
#   writable (e.g. a host path not visible inside the container) instead of
#   silently writing the scratch data to the project filesystem, and when
#   NEED_GB is given and ROOT has less free space than that
#   ($BRAKER4_TMP_NEED_GB, when set, replaces NEED_GB: tests, small machines).
make_tmp() {
    local prefix=$1 root=${2:-} need=${3:-} d free
    root=${root:-${TMPDIR:-/tmp}}
    [ -n "${BRAKER4_TMP_NEED_GB:-}" ] && [ -n "$need" ] && need=$BRAKER4_TMP_NEED_GB
    if [ -n "$need" ]; then
        free=$(df -Pk "$root" 2>/dev/null | awk 'NR == 2 { printf "%d", $4 / 1048576 }')
        if [ -n "$free" ] && [ "$free" -lt "$need" ]; then
            echo "ERROR: $root on $(hostname) has $free GB free, $need GB needed" >&2
            return 1
        fi
    fi
    if d=$(mktemp -d "$root/${prefix}_XXXXXX" 2>/dev/null); then
        echo "$d"
        return 0
    fi
    echo "ERROR: cannot create a temporary directory under $root; set" \
         "[paths] tmp_dir in config.ini to a writable path on the" \
         "compute nodes" >&2
    return 1
}

# scratch_dir VAR PREFIX ROOT NEED_GB FALLBACK
#   Set the variable VAR to a make_tmp directory, or to FALLBACK (a
#   directory in the run directory, created) with a warning when no scratch
#   directory with NEED_GB free can be made.  SCRATCH is set to the scratch
#   directory (empty when working in FALLBACK) so the rule can remove it on
#   exit.  Not for command substitution: $(...) runs in a subshell and the
#   variables would be lost.
#       scratch_dir outDir busco_genome_X "{params.tmp_root}" 20 output/X/busco
#       trap 'rm -rf -- "$SCRATCH"' EXIT
#   The job then copies the kept files from $outDir to FALLBACK (copy_back).
SCRATCH=""
scratch_dir() {
    local var=$1 prefix=$2 root=$3 need=$4 fallback=$5
    if SCRATCH=$(make_tmp "$prefix" "$root" "$need"); then
        echo "scratch directory: $SCRATCH" >&2
        printf -v "$var" '%s' "$SCRATCH"
    else
        SCRATCH=""
        echo "WARNING: working in $fallback instead (slow on network file systems)" >&2
        mkdir -p "$fallback"
        printf -v "$var" '%s' "$fallback"
    fi
}

# need_gb FACTOR FILE...
#   Print FACTOR x the total size of FILE... in GB, rounded up, plus 5 GB
#   headroom: the scratch space a job needs when its output scales with its
#   input (minimap2's SAM is ~3x the gzipped reads, a BAM ~1/3 of its SAM).
#   Missing files count as 0.  FACTOR is an integer.
need_gb() {
    local factor=$1 bytes=0 f s
    shift
    for f in "$@"; do
        s=$(stat -L -c %s "$f" 2>/dev/null || echo 0)
        bytes=$((bytes + s))
    done
    echo $(( (bytes * factor + 1073741823) / 1073741824 + 5 ))
}

# copy_back SRC DST [NAME...]
#   Copy the entries NAME... (files or directories, relative to SRC) from
#   SRC to DST; without names, everything in SRC.  Does nothing when SRC
#   and DST are the same directory (scratch_dir fell back to the run
#   directory).  Missing names are skipped: the rule checks its outputs.
copy_back() {
    local src=$1 dst=$2 n
    shift 2
    [ "$(cd "$src" && pwd -P)" = "$(mkdir -p "$dst" && cd "$dst" && pwd -P)" ] && return 0
    if [ $# -eq 0 ]; then
        cp -r "$src/." "$dst/"
        return
    fi
    for n in "$@"; do
        [ -e "$src/$n" ] || continue
        mkdir -p "$dst/$(dirname "$n")"
        rm -rf -- "${dst:?}/$n"
        cp -r "$src/$n" "$dst/$n"
    done
}
