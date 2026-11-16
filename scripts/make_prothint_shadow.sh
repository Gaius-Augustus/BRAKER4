#!/usr/bin/env bash
# Build a symlinked copy of the container's GeneMark-ETP bin tree in which
# ProtHint's run_spliced_alignment.pl is replaced by spaln_dispatcher.py
# (issue #98: the Perl-threads dispatcher can hang after the last batch is
# enqueued). Everything else points back to the original files, and the
# scripts that matter locate their helpers from the path they were called
# with (gmetp.pl: FindBin $Bin, prothint.py: abspath(__file__),
# spalnBatch.sh: readlink -f of its directory), so they pick up the copy.
#
# super.pl is also replaced, by a patched copy (#110): it writes the number
# of introns into the third column of set.out, where parse_set.pl expects
# the number of exons ("## gene_id cds_length exons_in_gene", the gmhmme3
# training format). Single-exon genes then count as -1 introns, and on
# genomes with few multi-intron genes (S. cerevisiae) the fallback estimate
# gives $ToInternalExon > 1 and $ToTerminalExon < 0; gmhmme3 stops with
# "negative argument in LOG function" in the masking-penalty step. The
# copy is made with sed at run time, so no GeneMark code ships with BRAKER4.
#
# Usage: make_prothint_shadow.sh <dest_dir> [etp_bin_dir]
#   etp_bin_dir defaults to the directory of gmetp.pl on PATH, else /opt/ETP/bin.
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
    # gmetp.pl is on PATH in the BRAKER3 container; without it, fall back to
    # the container's install location instead of resolving to "." (#75).
    gmetp=$(command -v gmetp.pl || true)
    if [ -n "$gmetp" ]; then
        src=$(dirname "$(readlink -f "$gmetp")")
    elif [ -f /opt/ETP/bin/gmetp.pl ]; then
        src=/opt/ETP/bin
    else
        echo "ERROR: gmetp.pl is not on PATH and /opt/ETP/bin/gmetp.pl does not exist; pass the GeneMark-ETP bin directory as second argument" >&2
        exit 1
    fi
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

# super.pl: exon count instead of intron count in set.out (#110). Patched
# only if the original line is found, so a fixed upstream version stays as
# it is.
super_pat='(scalar @introns) ."\n"'
if [ -f "$src/super.pl" ] && grep -qF "$super_pat" "$src/super.pl"; then
    rm "$dest/super.pl"
    sed 's/(scalar @introns) \."\\n"/(scalar @introns + 1) ."\\n"/' "$src/super.pl" > "$dest/super.pl"
    chmod +x "$dest/super.pl"
    if grep -qF '(scalar @introns + 1)' "$dest/super.pl"; then
        echo "[INFO] super.pl patched: exons_in_gene in set.out (#110)"
    else
        echo "WARNING: super.pl patch did not apply, using the original" >&2
        rm "$dest/super.pl"
        ln -s "$src/super.pl" "$dest/super.pl"
    fi
fi
