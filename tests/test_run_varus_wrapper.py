"""
Tests for scripts/run_varus_wrapper.sh with fake varus and samtools on PATH
that record their arguments. No network, no pyVARUS.
"""

import os
import stat
import subprocess

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
WRAPPER = os.path.join(ROOT, "scripts", "run_varus_wrapper.sh")

FAKE_VARUS = r"""#!/bin/bash
echo "varus $*" >> "$FAKE_DIR/calls.txt"
sub="$1"
outdir=""
args=("$@")
for ((i=0;i<${#args[@]};i++)); do
  [ "${args[$i]}" = "--outdir" ] && outdir="${args[$((i+1))]}"
done
case "$sub" in
  runlist)
    mkdir -p "$outdir"
    [ "${FAKE_RUNLIST:-}" = nohits ] && { echo "RuntimeError: No SRA RNA-seq runs found for species 'x'"; exit 1; }
    [ "${FAKE_RUNLIST:-}" = neterror ] && { echo "URLError: connection reset"; exit 1; }
    printf 'Run\tSpecies\nSRR1\tx\nSRR2\tx\nSRR3\tx\n' > "$outdir/Runlist.tsv" ;;
  index)
    mkdir -p "$outdir" ;;
  run)
    mkdir -p "$outdir"
    echo bam > "$outdir/VARUS.bam"
    [ -n "${FAKE_RUN_NO_BAM:-}" ] && rm -f "$outdir/VARUS.bam"
    [ -n "${FAKE_RUN_RC:-}" ] && { [ "$FAKE_RUN_RC" = 3 ] && rm -f "$outdir/VARUS.bam"; exit "$FAKE_RUN_RC"; }
    ;;
esac
exit 0
"""

FAKE_SAMTOOLS = r"""#!/bin/bash
echo "samtools $*" >> "$FAKE_DIR/calls.txt"
case "$1" in
  sort) while [ $# -gt 0 ]; do [ "$1" = -o ] && out="$2"; shift; done
        echo sorted > "$out" ;;
  index) for a; do last="$a"; done; echo csi > "$last.csi" ;;
  view) echo 42 ;;
esac
exit 0
"""


@pytest.fixture
def env(tmp_path):
    fake = tmp_path / "fake"
    fake.mkdir()
    for name, body in (("varus", FAKE_VARUS), ("samtools", FAKE_SAMTOOLS)):
        p = fake / name
        p.write_text(body)
        p.chmod(p.stat().st_mode | stat.S_IEXEC)
    work = tmp_path / "work"
    work.mkdir()
    (work / "genome.fa").write_text(">c\nACGT\n")
    for n in ("a_R1.fastq.gz", "a_R2.fastq.gz", "b_1.fq.gz", "b_2.fq.gz"):
        (work / n).write_text("x")
    e = {k: v for k, v in os.environ.items() if k != "NCBI_EMAIL"}
    e["PATH"] = f"{fake}:{e['PATH']}"
    e["FAKE_DIR"] = str(fake)
    return {"env": e, "work": work, "fake": fake}


def run(env, species, libs=(), max_batches="1000", use_logan="0",
        extra_env=None, scratch=False):
    work = env["work"]
    e = dict(env["env"], **(extra_env or {}))
    varus_dir = (work / ("scratch" if scratch else "out/varus"))
    cmd = ["bash", WRAPPER, str(varus_dir), "genome.fa", species, "2",
           "out/varus/varus.sorted.bam", "log.txt", use_logan, max_batches, *libs]
    p = subprocess.run(cmd, cwd=work, env=e, capture_output=True, text=True)
    calls = env["fake"] / "calls.txt"
    lines = calls.read_text().splitlines() if calls.exists() else []
    log = (work / "log.txt").read_text() if (work / "log.txt").exists() else ""
    return p, lines, log


def varus_run_line(lines):
    return [l for l in lines if l.startswith("varus run ") and not l.startswith("varus runlist")][0]


def test_local_only(env):
    p, lines, log = run(env, "", ["a_R1.fastq.gz,a_R2.fastq.gz", "b_1.fq.gz,b_2.fq.gz"],
                        max_batches="500")
    assert p.returncode == 0, p.stderr + log
    assert not any(l.startswith("varus runlist") for l in lines)
    run_line = varus_run_line(lines)
    w = env["work"]
    assert "--runlist" not in run_line
    assert f"--fastq {w}/a_R1.fastq.gz,{w}/a_R2.fastq.gz" in run_line
    assert f"--fastq {w}/b_1.fq.gz,{w}/b_2.fq.gz" in run_line
    assert "--max-batches 500" in run_line and "--keep-bam" in run_line
    assert "--no-logan" in run_line
    assert run_line.startswith("varus run local FASTQ libraries ")
    out = w / "out" / "varus"
    assert not (out / "varus_runlist.tsv").exists()
    stats = (out / "varus_stats.txt").read_text()
    assert "Local FASTQ libraries: 2" in stats
    assert "SRA runs in runlist" not in stats
    assert "Max batches: 500" in stats and "Total reads in BAM: 42" in stats
    assert (out / "varus.sorted.bam").exists() and (out / "varus.sorted.bam.csi").exists()


def test_sra_only(env):
    p, lines, log = run(env, "Homo sapiens", extra_env={"NCBI_EMAIL": "a@b.c"},
                        use_logan="1")
    assert p.returncode == 0, p.stderr + log
    rl = [l for l in lines if l.startswith("varus runlist")][0]
    assert rl.startswith("varus runlist Homo sapiens ") and "--email a@b.c" in rl
    run_line = varus_run_line(lines)
    assert "--runlist" in run_line and "--fastq" not in run_line
    assert "--no-logan" not in run_line
    assert run_line.startswith("varus run Homo sapiens ")
    assert "--max-batches 1000" in run_line
    out = env["work"] / "out" / "varus"
    assert (out / "varus_runlist.tsv").exists()
    stats = (out / "varus_stats.txt").read_text()
    assert "Species: Homo sapiens" in stats
    assert "SRA runs in runlist: 3" in stats
    assert "Local FASTQ libraries" not in stats


def test_both_with_scratch_dir(env):
    p, lines, log = run(env, "Homo sapiens", ["a_R1.fastq.gz,a_R2.fastq.gz"], scratch=True)
    assert p.returncode == 0, p.stderr + log
    run_line = varus_run_line(lines)
    assert "--runlist" in run_line and "--fastq" in run_line
    out = env["work"] / "out" / "varus"
    stats = (out / "varus_stats.txt").read_text()
    assert "SRA runs in runlist: 3" in stats and "Local FASTQ libraries: 1" in stats
    assert (out / "varus_runlist.tsv").exists()
    assert (out / "varus.sorted.bam").exists()  # copied from scratch


def test_symlink_not_resolved(env):
    w = env["work"]
    (w / "real").mkdir()
    (w / "real" / "x.fastq.gz").write_text("x")
    (w / "link_R1.fastq.gz").symlink_to(w / "real" / "x.fastq.gz")
    p, lines, log = run(env, "", ["link_R1.fastq.gz,a_R2.fastq.gz"])
    assert p.returncode == 0, log
    assert f"--fastq {w}/link_R1.fastq.gz," in varus_run_line(lines)


def test_exit_4_with_bam_is_warning(env):
    p, lines, log = run(env, "", ["a_R1.fastq.gz,a_R2.fastq.gz"], extra_env={"FAKE_RUN_RC": "4"})
    assert p.returncode == 0, log
    assert "[WARNING]" in log and "exit status 4" in log
    assert (env["work"] / "out" / "varus" / "varus_stats.txt").exists()


def test_exit_3_is_error(env):
    p, lines, log = run(env, "", ["a_R1.fastq.gz,a_R2.fastq.gz"], extra_env={"FAKE_RUN_RC": "3"})
    assert p.returncode == 1
    assert "[ERROR]" in log and "5%" in log
    assert not any(l.startswith("samtools") for l in lines)
    assert not (env["work"] / "out" / "varus" / "varus_stats.txt").exists()


def test_other_exit_status_propagates(env):
    p, lines, log = run(env, "", ["a_R1.fastq.gz,a_R2.fastq.gz"], extra_env={"FAKE_RUN_RC": "7"})
    assert p.returncode == 7
    assert "exit status 7" in log


def test_empty_email_sends_none(env):
    p, lines, log = run(env, "Homo sapiens", extra_env={"NCBI_EMAIL": ""})
    assert p.returncode == 0, log
    rl = [l for l in lines if l.startswith("varus runlist")][0]
    assert "--email" not in rl


def test_unset_email_uses_default(env):
    p, lines, log = run(env, "Homo sapiens")
    rl = [l for l in lines if l.startswith("varus runlist")][0]
    assert "--email katharina" in rl


def test_missing_fastq(env):
    p, lines, log = run(env, "", ["a_R1.fastq.gz,nope_R2.fastq.gz"])
    assert p.returncode == 1
    assert "nope_R2.fastq.gz is missing or empty" in log
    assert lines == []


def test_empty_fastq(env):
    (env["work"] / "a_R2.fastq.gz").write_text("")
    p, lines, log = run(env, "", ["a_R1.fastq.gz,a_R2.fastq.gz"])
    assert p.returncode == 1 and "missing or empty" in log


def test_neither_species_nor_libs(env):
    p, lines, log = run(env, "")
    assert p.returncode == 1
    assert "neither a species name nor" in p.stderr
    assert lines == []


def test_no_sra_runs_continues_with_local_libraries(env):
    p, lines, log = run(env, "Genus species", ["a_R1.fastq.gz,a_R2.fastq.gz"],
                        extra_env={"FAKE_RUNLIST": "nohits"})
    assert p.returncode == 0, p.stderr + log
    assert "Sampling from the local FASTQ libraries only" in log
    run_line = varus_run_line(lines)
    assert "--runlist" not in run_line and "--fastq" in run_line
    out = env["work"] / "out" / "varus"
    assert not (out / "varus_runlist.tsv").exists()
    stats = (out / "varus_stats.txt").read_text()
    assert "SRA runs in runlist" not in stats
    assert "SRA: no RNA-Seq runs found" in stats and "Local FASTQ libraries: 1" in stats


def test_no_sra_runs_without_local_libraries_fails(env):
    p, lines, log = run(env, "Genus species", extra_env={"FAKE_RUNLIST": "nohits"})
    assert p.returncode != 0
    assert not any(l.startswith("varus run ") for l in lines)


def test_other_runlist_error_fails_even_with_local_libraries(env):
    p, lines, log = run(env, "Genus species", ["a_R1.fastq.gz,a_R2.fastq.gz"],
                        extra_env={"FAKE_RUNLIST": "neterror"})
    assert p.returncode != 0
    assert "varus runlist ended with exit status 1" in log
    assert not any(l.startswith("varus run ") for l in lines)


def test_missing_bam_after_successful_run_fails(env):
    p, lines, log = run(env, "", ["a_R1.fastq.gz,a_R2.fastq.gz"],
                        extra_env={"FAKE_RUN_NO_BAM": "1"})
    assert p.returncode != 0
    assert "no usable alignments" in log
    assert not any(l.startswith("samtools sort") for l in lines)
