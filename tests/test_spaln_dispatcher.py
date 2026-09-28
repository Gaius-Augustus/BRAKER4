"""
Tests for scripts/spaln_dispatcher.py (issue #98), the process-based
replacement for ProtHint's bin/run_spliced_alignment.pl.

The script is imported and driven in-process via main(argv) so that
capsys can capture its stderr logging. It shells out to two stub
scripts standing in for spalnBatch.sh and gff_from_region_to_contig.pl,
installed in a fake --prothint_bin directory under tmp_path. All input
FASTA/pair data is synthetic and generated in-memory.

The stub spalnBatch.sh writes the content of each pair's nuc/prot file
plus one "ARGS:<batch>:<min_exon_score>:<nonCanonical>:<longGene>:<longProtein>"
line to its output file, and deletes the nuc/prot files -- mirroring the
real script closely enough to verify sequence cleaning, argument
plumbing and batch ordering. Two env vars make it misbehave on demand:

  STALL_BATCH / STALL_ONCE  -- the named batch execs "sleep 600" (0 CPU)
                                instead of running; with STALL_ONCE=1 it
                                only does this on the first attempt (a
                                marker file records the attempt).
  PLAIN_SLEEP_BATCH / PLAIN_SLEEP_SECONDS -- the named batch sleeps for a
                                few seconds (consuming no CPU either) and
                                then completes normally; used to check
                                that --stall_timeout 0 disables the
                                watchdog.

The stub gff_from_region_to_contig.pl just copies --in_gff to --out_gff.
"""

import os
import re
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, SCRIPTS)

import spaln_dispatcher as sd  # noqa: E402


STUB_SPALNBATCH = r"""#!/bin/bash
batch="$1"
out="$2"
min_exon_score="$3"
noncanonical="$4"
longgene="$5"
longprotein="$6"

if [ -n "$STALL_BATCH" ] && [ "$batch" = "$STALL_BATCH" ]; then
    if [ "$STALL_ONCE" = "1" ] && [ -f "stalled_once_${batch}" ]; then
        :  # second attempt: fall through and complete normally
    else
        if [ "$STALL_ONCE" = "1" ]; then
            touch "stalled_once_${batch}"
        fi
        echo $$ > "stall_pid_${batch}"
        exec sleep 600
    fi
fi

if [ -n "$PLAIN_SLEEP_BATCH" ] && [ "$batch" = "$PLAIN_SLEEP_BATCH" ]; then
    sleep "$PLAIN_SLEEP_SECONDS"
fi

: > "$out"
while IFS=$'\t' read -r nuc prot; do
    cat "$nuc" >> "$out"
    cat "$prot" >> "$out"
    rm -f "$nuc" "$prot"
done < "$batch"
echo "ARGS:${batch}:${min_exon_score}:${noncanonical}:${longgene}:${longprotein}" >> "$out"
"""

STUB_TO_CONTIG = r"""#!/bin/bash
in_gff=""
out_gff=""
while [ $# -gt 0 ]; do
    case "$1" in
        --in_gff) in_gff="$2"; shift 2 ;;
        --out_gff) out_gff="$2"; shift 2 ;;
        *) shift ;;
    esac
done
cp "$in_gff" "$out_gff"
"""


def make_prothint_bin(tmp_path, missing_spalnbatch=False):
    bin_dir = tmp_path / "prothint_bin"
    bin_dir.mkdir()
    if not missing_spalnbatch:
        spalnbatch = bin_dir / "spalnBatch.sh"
        spalnbatch.write_text(STUB_SPALNBATCH)
        spalnbatch.chmod(0o755)
    to_contig = bin_dir / "gff_from_region_to_contig.pl"
    to_contig.write_text(STUB_TO_CONTIG)
    to_contig.chmod(0o755)
    return str(bin_dir)


def make_pairs(tmp_path, n):
    """n unique nuc/prot pairs. Sequences carry lowercase letters, digits
    and '*' that read_sequences() must uppercase/strip to "ACGTT"/"MNPQQ"."""
    list_lines, nuc_lines, prot_lines = [], [], []
    for i in range(1, n + 1):
        nuc_id, prot_id = "n" + str(i), "p" + str(i)
        list_lines.append(nuc_id + "\t" + prot_id + "\n")
        nuc_lines.append(">" + nuc_id + " desc\nacgt123*t\n")
        prot_lines.append(">" + prot_id + " desc\nmnpq99*q\n")
    (tmp_path / "pairs.lst").write_text("".join(list_lines))
    (tmp_path / "nuc.fa").write_text("".join(nuc_lines))
    (tmp_path / "prot.fa").write_text("".join(prot_lines))
    return (str(tmp_path / "nuc.fa"), str(tmp_path / "prot.fa"),
            str(tmp_path / "pairs.lst"))


def leftover_temp_files(tmp_path):
    return sorted(p.name for p in tmp_path.iterdir()
                  if p.name.startswith(("nuc_", "prot_", "batch_")))


def stalled_pid(tmp_path, batch):
    """PID the stub wrote to stall_pid_<batch> just before it exec'd
    "sleep 600", or None if that batch never stalled."""
    pid_file = tmp_path / ("stall_pid_" + batch)
    if not pid_file.exists():
        return None
    return int(pid_file.read_text().strip())


def process_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


# ---------------------------------------------------------------------------
# 1. Normal run
# ---------------------------------------------------------------------------

def test_normal_run_250_pairs_two_cores(tmp_path, monkeypatch, capsys):
    nuc, prot, lst = make_pairs(tmp_path, 250)
    bin_dir = make_prothint_bin(tmp_path)
    monkeypatch.chdir(tmp_path)

    sd.main(["--nuc", nuc, "--prot", prot, "--list", lst,
             "--cores", "2", "--aligner", "spaln", "--verbose",
             "--prothint_bin", bin_dir])

    err = capsys.readouterr().err
    assert "100.0%" in err
    assert "250/250 pairs aligned" in err

    gff = (tmp_path / "spaln.gff").read_text()
    ids = re.findall(r">n(\d+)", gff)
    assert ids == [str(i) for i in range(1, 251)]
    assert gff.count("ACGTT") == 250
    assert gff.count("MNPQQ") == 250
    assert gff.count("ARGS:batch_") == 3
    assert gff.index("ARGS:batch_1") < gff.index("ARGS:batch_2") < gff.index("ARGS:batch_3")

    assert leftover_temp_files(tmp_path) == []
    assert not (tmp_path / "spaln.regions.gff").exists()


# ---------------------------------------------------------------------------
# 2. prothint.py-style argument abbreviation / types
# ---------------------------------------------------------------------------

def test_prothint_style_arguments(tmp_path, monkeypatch):
    nuc, prot, lst = make_pairs(tmp_path, 2)
    bin_dir = make_prothint_bin(tmp_path)
    monkeypatch.chdir(tmp_path)

    sd.main(["--cores", "1", "--nuc", nuc, "--list", lst, "--prot", prot,
             "--v", "--aligner", "spaln", "--min_exon_score", "25",
             "--nonCanonical", "--longGene", "30000", "--longProtein", "15000",
             "--prothint_bin", bin_dir])

    gff = (tmp_path / "spaln.gff").read_text()
    # min_exon_score must render as "25", not "25.0"; nonCanonical as "1"
    assert "ARGS:batch_1:25:1:30000:15000" in gff


# ---------------------------------------------------------------------------
# 3. Stall once, then succeed on retry
# ---------------------------------------------------------------------------

def test_stall_once_then_succeeds(tmp_path, monkeypatch, capsys):
    nuc, prot, lst = make_pairs(tmp_path, 5)
    bin_dir = make_prothint_bin(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sd, "CPU_CHECK_SECONDS", 0.3)
    monkeypatch.setattr(sd, "POLL_SECONDS", 0.05)
    monkeypatch.setattr(sd, "KILL_WAIT_SECONDS", 5)
    monkeypatch.setenv("STALL_BATCH", "batch_1")
    monkeypatch.setenv("STALL_ONCE", "1")

    sd.main(["--nuc", nuc, "--prot", prot, "--list", lst,
             "--cores", "1", "--aligner", "spaln",
             "--stall_timeout", "1", "--prothint_bin", bin_dir])

    err = capsys.readouterr().err
    assert "WARNING" in err and "batch_1" in err and "retrying" in err

    gff = (tmp_path / "spaln.gff").read_text()
    ids = re.findall(r">n(\d+)", gff)
    assert ids == [str(i) for i in range(1, 6)]
    assert gff.count("ACGTT") == 5
    assert gff.count("MNPQQ") == 5

    assert leftover_temp_files(tmp_path) == []


# ---------------------------------------------------------------------------
# 4. Stall on every attempt: batch is skipped, others still succeed
# ---------------------------------------------------------------------------

def test_stall_always_batch_is_skipped(tmp_path, monkeypatch, capsys):
    nuc, prot, lst = make_pairs(tmp_path, 150)  # batch_1 (100) + batch_2 (50)
    bin_dir = make_prothint_bin(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sd, "CPU_CHECK_SECONDS", 0.3)
    monkeypatch.setattr(sd, "POLL_SECONDS", 0.05)
    monkeypatch.setattr(sd, "KILL_WAIT_SECONDS", 5)
    monkeypatch.setenv("STALL_BATCH", "batch_2")

    sd.main(["--nuc", nuc, "--prot", prot, "--list", lst,
             "--cores", "2", "--aligner", "spaln",
             "--stall_timeout", "1", "--prothint_bin", bin_dir])

    err = capsys.readouterr().err
    assert "WARNING" in err and "batch_2" in err

    gff = (tmp_path / "spaln.gff").read_text()
    assert "ARGS:batch_1" in gff
    assert "ARGS:batch_2" not in gff
    ids = set(re.findall(r">n(\d+)", gff))
    assert ids == {str(i) for i in range(1, 101)}  # batch_2's n101..n150 skipped

    assert leftover_temp_files(tmp_path) == []
    assert (tmp_path / "spaln.gff").exists()

    pid = stalled_pid(tmp_path, "batch_2")
    assert pid is not None
    assert not process_alive(pid)


# ---------------------------------------------------------------------------
# 5. Missing spalnBatch.sh
# ---------------------------------------------------------------------------

def test_missing_spalnbatch_exits_nonzero(tmp_path, monkeypatch):
    nuc, prot, lst = make_pairs(tmp_path, 2)
    bin_dir = make_prothint_bin(tmp_path, missing_spalnbatch=True)
    monkeypatch.chdir(tmp_path)

    with pytest.raises(SystemExit) as exc_info:
        sd.main(["--nuc", nuc, "--prot", prot, "--list", lst,
                 "--cores", "1", "--aligner", "spaln",
                 "--prothint_bin", bin_dir])
    assert exc_info.value.code != 0


# ---------------------------------------------------------------------------
# 6. --stall_timeout 0 disables the watchdog
# ---------------------------------------------------------------------------

def test_stall_timeout_zero_disables_watchdog(tmp_path, monkeypatch):
    nuc, prot, lst = make_pairs(tmp_path, 3)
    bin_dir = make_prothint_bin(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PLAIN_SLEEP_BATCH", "batch_1")
    monkeypatch.setenv("PLAIN_SLEEP_SECONDS", "2")

    sd.main(["--nuc", nuc, "--prot", prot, "--list", lst,
             "--cores", "1", "--aligner", "spaln",
             "--stall_timeout", "0", "--prothint_bin", bin_dir])

    gff = (tmp_path / "spaln.gff").read_text()
    assert "ARGS:batch_1" in gff
    ids = re.findall(r">n(\d+)", gff)
    assert ids == ["1", "2", "3"]
    assert leftover_temp_files(tmp_path) == []
