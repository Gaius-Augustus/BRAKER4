"""
Tests for the run status line of best_by_compleasm.py and its use in the report.

best_by_compleasm used to be called with `|| true`, and a compleasm crash left
the log saying "original was already best" (seen in #91). The script now writes
one line to --status_file: "OK: ...", "COMPLEASM_FAILED: <reason>" or
"FAILED: <reason>". These tests run the script end to end with fake
compleasm.py, getAnnoFastaFromJoingenes.py and tsebra.py executables.
"""

import os
import stat
import subprocess
import sys
import importlib.util

import pytest

_SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")
_BBC = os.path.join(_SCRIPTS, "best_by_compleasm.py")

_FAKE_COMPLEASM = r'''#!/usr/bin/env python3
import os, sys
mode = os.environ.get("FAKE_COMPLEASM_MODE", "ok")
if mode == "crash":
    sys.stderr.write("running hmmsearch\n")
    sys.stderr.write("Error: Target sequence length > 100K, over comparison pipeline limit\n")
    sys.exit(1)
args = sys.argv[1:]
out = args[args.index("-o") + 1]
prot = os.path.basename(args[args.index("-p") + 1])
os.makedirs(out, exist_ok=True)
if mode == "nosummary":
    sys.exit(0)
missing = {"braker": 1.0}.get(prot.split(".")[0], 5.0)
with open(os.path.join(out, "summary.txt"), "w") as f:
    f.write("## lineage: eukaryota_odb12\nS:90.00%, 90\nD:0.00%, 0\n"
            "F:0.00%, 0\nI:0.00%, 0\nM:" + format(missing, ".2f") + "%, 5\nN:100\n")
'''

_FAKE_GETANNO = r'''#!/usr/bin/env python3
import sys
args = sys.argv[1:]
stem = args[args.index("-o") + 1]
with open(stem + ".aa", "w") as f:
    f.write(">g1.t1\nMAAA\n")
'''

_FAKE_TSEBRA = "#!/bin/sh\nexit 0\n"


def _write_exe(path, content):
    with open(path, "w") as f:
        f.write(content)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def setup(tmp_path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    _write_exe(bindir / "compleasm.py", _FAKE_COMPLEASM)
    _write_exe(bindir / "getAnnoFastaFromJoingenes.py", _FAKE_GETANNO)
    _write_exe(bindir / "tsebra.py", _FAKE_TSEBRA)

    indir = tmp_path / "stage"
    (indir / "GeneMark-ETP").mkdir(parents=True)
    gtf = 'chr1\tAUGUSTUS\tCDS\t1\t12\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n'
    for name in ("braker.gtf", "augustus.hints.gtf", "GeneMark-ETP/genemark.gtf"):
        (indir / name).write_text(gtf)
    for name in ("braker.aa", "augustus.hints.aa"):
        (indir / name).write_text(">g1.t1\nMAAA\n")
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATGGCCGCCGCCTAA\n")

    lineage = tmp_path / "lib" / "eukaryota_odb12"
    (lineage / "hmms").mkdir(parents=True)
    (lineage / "hmms" / "1at2759.hmm").write_text("")
    (lineage / "scores_cutoff").write_text("1at2759 10.0\n")

    def run(mode):
        env = dict(os.environ)
        env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
        env["FAKE_COMPLEASM_MODE"] = mode
        status_file = tmp_path / "bbc_status.txt"
        proc = subprocess.run(
            [sys.executable, _BBC, "-m", str(tmp_path / "work"), "-d", str(indir),
             "-g", str(genome), "-p", "eukaryota_odb12",
             "-L", str(tmp_path / "lib"), "-s", str(status_file)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, env=env)
        status = status_file.read_text().strip() if status_file.exists() else None
        return proc, status

    return run


def test_ok_when_original_is_best(setup):
    proc, status = setup("ok")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert status.startswith("OK: original gene set kept")
    assert "STATUS: " + status in proc.stdout


def test_compleasm_crash_is_reported(setup):
    proc, status = setup("crash")
    assert proc.returncode != 0
    assert status.startswith("COMPLEASM_FAILED: compleasm protein on ")
    assert "Target sequence length > 100K" in status
    # the crash must not be reported as a regular outcome
    assert "is the best one" not in proc.stdout


def test_missing_summary_is_reported(setup):
    proc, status = setup("nosummary")
    assert proc.returncode != 0
    assert status.startswith("COMPLEASM_FAILED: compleasm wrote no summary.txt")


def test_missing_lineage_is_reported(setup, tmp_path):
    (tmp_path / "lib" / "eukaryota_odb12" / "scores_cutoff").unlink()
    import shutil
    shutil.rmtree(tmp_path / "lib" / "eukaryota_odb12")
    proc, status = setup("ok")
    assert proc.returncode != 0
    assert status.startswith("COMPLEASM_FAILED: pre-downloaded BUSCO lineage not found")


# ---------------------------------------------------------------------------
# generate_report.py
# ---------------------------------------------------------------------------

def _import_report():
    spec = importlib.util.spec_from_file_location(
        "generate_report", os.path.join(_SCRIPTS, "generate_report.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_report = _import_report()


def test_report_reads_status(tmp_path):
    log = tmp_path / "best_by_compleasm.log"
    log.write_text("STATUS: COMPLEASM_FAILED: compleasm exited with code 1\n")
    assert _report.read_bbc_status(str(log)) == "COMPLEASM_FAILED: compleasm exited with code 1"
    log.write_text("BRAKER is missing 1.0 BUSCOs.\n")
    assert _report.read_bbc_status(str(log)) is None


def test_report_warns_on_failures(tmp_path):
    work = tmp_path / "work"
    out = tmp_path / "out"
    (work / "compleasm_proteins").mkdir(parents=True)
    (work / "best_by_compleasm.log").write_text(
        "STATUS: COMPLEASM_FAILED: compleasm wrote no summary.txt\n")
    (work / "compleasm_proteins" / "summary.txt").write_text(
        "COMPLEASM_FAILED: compleasm exited with code 1 and wrote no summary.txt\n")
    warnings = _report.collect_run_warnings(str(work), str(out))
    assert len(warnings) == 2
    assert "best_by_compleasm failed" in warnings[0]
    assert "proteome failed" in warnings[1]

    methods = _report.generate_methods_text(str(work), "ES (ab initio)")
    assert "best_by_compleasm" in methods and "failed" in methods


def test_report_no_warnings_on_success(tmp_path):
    work = tmp_path / "work"
    (work / "compleasm_proteins").mkdir(parents=True)
    (work / "best_by_compleasm.log").write_text(
        "STATUS: OK: original gene set kept, it has the fewest missing BUSCOs (1.0%)\n"
        "BRAKER is missing 1.0 BUSCOs.\n")
    (work / "compleasm_proteins" / "summary.txt").write_text("M:1.00%, 1\n")
    assert _report.collect_run_warnings(str(work), str(tmp_path / "out")) == []
