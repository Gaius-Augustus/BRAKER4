#!/usr/bin/env bash
# Build a symlinked copy of the container's GeneMark-ETP bin tree in which
# ProtHint's run_spliced_alignment.pl is replaced by spaln_dispatcher.py
# (issue #98: the Perl-threads dispatcher can hang after the last batch is
# enqueued). Everything else points back to the original files, and the
# scripts that matter locate their helpers from the path they were called
# with (gmetp.pl: FindBin $Bin, prothint.py: abspath(__file__),
# spalnBatch.sh: readlink -f of its directory), so they pick up the copy.
#
# Usage: make_prothint_shadow.sh <dest_dir> [etp_bin_dir]
#   etp_bin_dir defaults to the directory of gmetp.pl on PATH.
# Afterwards use <dest_dir>/gmetp.pl or
# <dest_dir>/gmes/ProtHint/bin/prothint.py instead of the ones on PATH.
set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <dest_dir> [etp_bin_dir]" >&2
    exit 1
fi

script_dir=$(dirname "$(readlink -f "$0")")
dispatcher=$script_dir/spaln_dispatcher.py

if [ $# -eq 2 ]; then
    src=$(readlink -f "$2")
else
    src=$(dirname "$(readlink -f "$(command -v gmetp.pl)")")
fi
prothint_bin=$src/gmes/ProtHint/bin
if [ ! -f "$prothint_bin/run_spliced_alignment.pl" ]; then
    echo "ERROR: $prothint_bin/run_spliced_alignment.pl not found" >&2
    exit 1
fi

# Only ever replace a directory this script created.
if [ -e "$1" ]; then
    if [ ! -f "$1/.prothint_shadow" ]; then
        echo "ERROR: $1 exists and was not created by $0" >&2
        exit 1
    fi
    rm -rf "$1"
fi
mkdir -p "$1"
dest=$(readlink -f "$1")
touch "$dest/.prothint_shadow"

# link_all_but <src_dir> <dest_dir> <entry>: real <dest_dir> holding
# symlinks to everything in <src_dir> except <entry>.
link_all_but() {
    mkdir -p "$2"
    local entry name
    for entry in "$1"/* "$1"/.[!.]*; do
        [ -e "$entry" ] || [ -L "$entry" ] || continue
        name=$(basename "$entry")
        [ "$name" = "$3" ] && continue
        ln -s "$entry" "$2/$name"
    done
}

link_all_but "$src" "$dest" gmes
link_all_but "$src/gmes" "$dest/gmes" ProtHint
link_all_but "$src/gmes/ProtHint" "$dest/gmes/ProtHint" bin
link_all_but "$prothint_bin" "$dest/gmes/ProtHint/bin" run_spliced_alignment.pl
cp "$dispatcher" "$dest/gmes/ProtHint/bin/run_spliced_alignment.pl"
chmod +x "$dest/gmes/ProtHint/bin/run_spliced_alignment.pl"
