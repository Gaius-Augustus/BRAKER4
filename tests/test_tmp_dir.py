"""Tests for scripts/tmp_dir.sh: make_tmp PREFIX [ROOT] [NEED_GB], scratch_dir
and copy_back.

Sourced into a bash subshell for each case (as rules do via
`source {script_dir}/tmp_dir.sh`).  Rules that write many small files work in a
scratch directory on the node-local disk and copy only the kept results
back (brain's CephFS needs ~0.4 s per file creation on bad days).
"""

import os
import shutil
import subprocess
from pathlib import Path

TMP_DIR_SH = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "tmp_dir.sh")
)
HELPER = Path(TMP_DIR_SH)


def _run_make_tmp(prefix, root=None, cwd=None, env=None):
    """Run `source scripts/tmp_dir.sh; make_tmp PREFIX [ROOT]` and return the
    completed process (stdout is the printed directory path)."""
    cmd = f"source {TMP_DIR_SH!r}; make_tmp {prefix!r}"
    if root is not None:
        cmd += f" {root!r}"
    return subprocess.run(
        ["bash", "-c", cmd],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
    )


def test_explicit_root(tmp_path):
    root = tmp_path / "scratch"
    root.mkdir()
    proc = _run_make_tmp("myprefix", root=str(root))
    assert proc.returncode == 0, proc.stderr
    out = proc.stdout.strip()
    assert out, proc.stderr
    assert os.path.isdir(out)
    assert os.path.dirname(out) == str(root)
    assert os.path.basename(out).startswith("myprefix_")


def test_empty_root_uses_tmpdir_env(tmp_path):
    tmpdir_root = tmp_path / "tmpdir_env"
    tmpdir_root.mkdir()
    env = dict(os.environ)
    env["TMPDIR"] = str(tmpdir_root)
    proc = _run_make_tmp("envprefix", root="", env=env)
    assert proc.returncode == 0, proc.stderr
    out = proc.stdout.strip()
    assert os.path.isdir(out)
    assert os.path.dirname(out) == str(tmpdir_root)


def test_empty_root_and_no_tmpdir_falls_back_to_slash_tmp():
    env = dict(os.environ)
    env.pop("TMPDIR", None)
    proc = _run_make_tmp("notmpdirprefix", root=None, env=env)
    assert proc.returncode == 0, proc.stderr
    out = proc.stdout.strip()
    try:
        assert os.path.isdir(out)
        assert out.startswith("/tmp/")
    finally:
        shutil.rmtree(out, ignore_errors=True)


def test_unwritable_root_fails_without_writing_to_cwd(tmp_path):
    # root's parent doesn't exist -> mktemp under it fails; make_tmp must
    # fail loudly instead of putting scratch data into the working directory.
    bad_root = tmp_path / "does_not_exist" / "nested"
    workdir = tmp_path / "work"
    workdir.mkdir()
    proc = _run_make_tmp("fallbackprefix", root=str(bad_root), cwd=str(workdir))
    assert proc.returncode != 0
    assert "ERROR" in proc.stderr
    assert "tmp_dir" in proc.stderr
    assert proc.stdout.strip() == ""
    assert os.listdir(workdir) == []


def test_two_calls_give_distinct_dirs(tmp_path):
    root = tmp_path / "scratch2"
    root.mkdir()
    proc1 = _run_make_tmp("dup", root=str(root))
    proc2 = _run_make_tmp("dup", root=str(root))
    assert proc1.returncode == 0 and proc2.returncode == 0
    out1 = proc1.stdout.strip()
    out2 = proc2.stdout.strip()
    assert out1 != out2
    assert os.path.isdir(out1)
    assert os.path.isdir(out2)


# ── NEED_GB, scratch_dir, copy_back ───────────────────────────────────────

def _bash(script, **env):
    return subprocess.run(
        ["bash", "-c", f"set -euo pipefail\nsource {HELPER}\n{script}"],
        capture_output=True, text=True, env={**os.environ, **env})


def test_make_tmp_fails_when_space_is_short(tmp_path):
    # No disk has 10 million GB free.
    proc = _bash(f'make_tmp x "{tmp_path}" 10000000')
    assert proc.returncode == 1
    assert "GB free" in proc.stderr and "needed" in proc.stderr
    assert list(tmp_path.iterdir()) == []


def test_need_gb_env_override_replaces_need(tmp_path):
    # BRAKER4_TMP_NEED_GB forces the space check to fail (fallback tests).
    proc = _bash(f'make_tmp x "{tmp_path}" 1', BRAKER4_TMP_NEED_GB="999999999")
    assert proc.returncode == 1
    assert "ERROR" in proc.stderr
    assert list(tmp_path.iterdir()) == []


def test_scratch_dir_prints_scratch_and_sets_SCRATCH(tmp_path):
    fallback = tmp_path / "wd" / "output" / "S1" / "busco"
    proc = _bash(f'scratch_dir d busco_genome_S1 "{tmp_path}" 0 "{fallback}"; '
                 f'echo "$d"; echo "SCRATCH=$SCRATCH"')
    assert proc.returncode == 0, proc.stderr
    lines = proc.stdout.strip().splitlines()
    d = Path(lines[0])
    assert d.is_dir() and d.parent == tmp_path and d.name.startswith("busco_genome_S1_")
    assert lines[1] == f"SCRATCH={d}"
    assert "scratch directory:" in proc.stderr
    assert not fallback.exists()


def test_scratch_dir_falls_back_to_run_dir_with_warning(tmp_path):
    fallback = tmp_path / "wd" / "output" / "S1" / "busco"
    proc = _bash(f'scratch_dir d busco_genome_S1 "{tmp_path}/none" 0 "{fallback}"; '
                 f'echo "$d"; echo "SCRATCH=$SCRATCH"')
    assert proc.returncode == 0, proc.stderr
    lines = proc.stdout.strip().splitlines()
    assert lines == [str(fallback), "SCRATCH="]
    assert fallback.is_dir()
    assert "WARNING: working in" in proc.stderr


def test_scratch_dir_trap_removes_scratch_on_exit(tmp_path):
    fallback = tmp_path / "wd" / "x"
    proc = _bash(f'scratch_dir d s "{tmp_path}" 0 "{fallback}"; '
                 'trap \'rm -rf -- "$SCRATCH"\' EXIT; echo "$d"; touch "$d/f"')
    assert proc.returncode == 0, proc.stderr
    assert not Path(proc.stdout.strip()).exists()


def test_need_gb_one_mib_doubled_is_six(tmp_path):
    f = tmp_path / "genome.fa"
    f.write_bytes(b"A" * (1 << 20))
    proc = _bash(f'need_gb 2 "{f}"')
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "6"


def test_need_gb_scales_input_size_plus_headroom(tmp_path):
    a = tmp_path / "a.fq.gz"
    a.write_bytes(b"x" * 10)
    proc = _bash(f'need_gb 3 "{a}" "{tmp_path}/missing"; need_gb 3')
    assert proc.returncode == 0, proc.stderr
    # 30 bytes round up to 1 GB, + 5 GB headroom; no files: 5
    assert proc.stdout.split() == ["6", "5"]


def test_copy_back_selected_entries(tmp_path):
    src = tmp_path / "scratch"
    (src / "rnaseq" / "stringtie").mkdir(parents=True)
    (src / "rnaseq" / "stringtie" / "transcripts_merged.gff").write_text("gff\n")
    (src / "rnaseq" / "hisat2").mkdir()
    (src / "rnaseq" / "hisat2" / "big.bam").write_text("x" * 10)
    (src / "data").mkdir()
    (src / "data" / "dna.fna").write_text(">c\nA\n")
    (src / "genemark.gtf").write_text("gtf\n")
    dst = tmp_path / "wd" / "output" / "S1" / "busco"
    proc = _bash(f'copy_back "{src}" "{dst}" rnaseq/stringtie/transcripts_merged.gff genemark.gtf missing.txt')
    assert proc.returncode == 0, proc.stderr
    assert (dst / "rnaseq" / "stringtie" / "transcripts_merged.gff").read_text() == "gff\n"
    assert (dst / "genemark.gtf").read_text() == "gtf\n"
    assert not (dst / "rnaseq" / "hisat2").exists()
    assert not (dst / "data").exists()
    assert not (dst / "missing.txt").exists()


def test_copy_back_replaces_stale_directory(tmp_path):
    src = tmp_path / "scratch"
    (src / "output_bins").mkdir(parents=True)
    (src / "output_bins" / "a.fa").write_text(">c\nA\n")
    dst = tmp_path / "wd" / "output" / "S1" / "compleasm"
    (dst / "output_bins").mkdir(parents=True)
    (dst / "output_bins" / "stale.fa").write_text("old\n")
    proc = _bash(f'copy_back "{src}" "{dst}" output_bins')
    assert proc.returncode == 0, proc.stderr
    assert sorted(p.name for p in (dst / "output_bins").iterdir()) == ["a.fa"]


def test_copy_back_everything_without_names(tmp_path):
    src = tmp_path / "scratch"
    (src / "sub").mkdir(parents=True)
    (src / "sub" / "a.json").write_text("{}")
    (src / "meta.json").write_text("{}")
    dst = tmp_path / "wd" / "output" / "S1" / "varus"
    proc = _bash(f'copy_back "{src}" "{dst}"')
    assert proc.returncode == 0, proc.stderr
    assert (dst / "meta.json").exists() and (dst / "sub" / "a.json").exists()


def test_copy_back_is_a_noop_for_the_same_directory(tmp_path):
    d = tmp_path / "wd" / "x"
    d.mkdir(parents=True)
    (d / "f").write_text("1\n")
    proc = _bash(f'copy_back "{d}" "{d}" f; copy_back "{d}" "{tmp_path}/wd/../wd/x"')
    assert proc.returncode == 0, proc.stderr
    assert (d / "f").read_text() == "1\n"


# ── Rules converted to node-local scratch (SCRATCH_PLAN.md, Tier 1) ───────

import re

import pytest

REPO = Path(TMP_DIR_SH).parent.parent
SCRATCH_RULES = {
    "busco_genome": "rules/quality_control/run_busco.smk",
    "busco_proteins": "rules/quality_control/run_busco.smk",
    "run_genemark_etp": "rules/genemark/run_genemark_etp.smk",
    "run_genemark_etp_isoseq": "rules/genemark/run_genemark_etp_isoseq.smk",
    "run_prothint": "rules/genemark/run_prothint.smk",
    "run_prothint_iter2": "rules/genemark/run_prothint_iter2.smk",
    "run_genemark_es": "rules/genemark/run_genemark_es.smk",
    "run_genemark_et": "rules/genemark/run_genemark_et.smk",
    "run_genemark_ep": "rules/genemark/run_genemark_ep.smk",
    "run_augustus_hints": "rules/augustus_predict/run_augustus_hints.smk",
    "run_augustus_hints_iter2": "rules/augustus_predict/run_augustus_hints_iter2.smk",
    "optimize_augustus": "rules/augustus_training/optimize_augustus.smk",
    "run_compleasm": "rules/quality_control/run_compleasm.smk",
    "run_varus": "rules/preprocessing/run_varus.smk",
}


def _rule_text(name, path):
    text = (REPO / path).read_text()
    start = text.index(f"rule {name}:")
    nxt = re.search(r"^rule \w+:", text[start + 1:], re.M)
    return text[start:start + 1 + nxt.start()] if nxt else text[start:]


@pytest.mark.parametrize("name", sorted(SCRATCH_RULES))
def test_rule_uses_scratch_dir_safely(name):
    rule = _rule_text(name, SCRATCH_RULES[name])
    shell = rule[rule.index("shell:"):]
    assert "tmp_root" in rule[:rule.index("shell:")], "params.tmp_root missing"
    assert "Scratch:" in rule, "docstring lacks a Scratch: paragraph"
    src = shell.index("source {script_dir}/tmp_dir.sh")
    call = shell.index("scratch_dir outDir ")
    assert src < call
    # scratch_dir sets variables: never in a subshell or command substitution
    assert not re.search(r"\$\(\s*scratch_dir|\|\s*scratch_dir|\(\s*scratch_dir", shell)
    # the next statement after the (possibly continued) call is the one trap
    lines = shell[call:].splitlines()
    i = 0
    while lines[i].rstrip().endswith("\\"):
        i += 1
    assert lines[i + 1].strip() == "trap 'rm -rf -- \"$SCRATCH\"' EXIT"
    assert len(re.findall(r"^\s*trap ", shell, re.M)) == 1
    assert "{resources.tmpdir}" not in shell
