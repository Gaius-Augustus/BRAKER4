"""Tests for scripts/genemark_fallback.sh (automatic mode fallback of
GeneMark-ET and GeneMark-EP to ES mode).

genemark_fallback_marker TOOL EXIT_CODE MARKER LOG... decides whether a
gmes_petap.pl run that ended without genemark.gtf may continue in ES mode.
"""

import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "genemark_fallback.sh"


def call(tmp_path, code, logs, tool="GeneMark-ET"):
    """Run the function; return (return code, marker text or None)."""
    marker = tmp_path / "mode_selection" / "marker.txt"
    paths = []
    for i, text in enumerate(logs):
        p = tmp_path / f"log{i}.txt"
        p.write_text(text)
        paths.append(str(p))
    cmd = 'source "$0"; genemark_fallback_marker "$@"'
    proc = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", cmd, str(SCRIPT), tool, str(code), str(marker)] + paths,
        capture_output=True, text=True,
    )
    return proc.returncode, (marker.read_text() if marker.exists() else None)


def test_error_exit_of_the_program_allows_fallback(tmp_path):
    rc, text = call(tmp_path, 2, [
        "some output\n",
        "Use of uninitialized value $ph0 in division (/) at parse_ET.pl line 105.\n"
        "Illegal division by zero at parse_ET.pl line 105.\n",
    ])
    assert rc == 0
    assert text.count("\n") == 1
    assert text.startswith("GeneMark-ET stopped with an error during training (exit code 2)")
    assert "last error message: Illegal division by zero at parse_ET.pl line 105." in text


def test_marker_names_the_tool_and_works_without_error_line(tmp_path):
    rc, text = call(tmp_path, 1, ["nothing useful here\n"], tool="GeneMark-EP")
    assert rc == 0
    assert text == ("GeneMark-EP stopped with an error during training (exit code 1) "
                    "and wrote no genemark.gtf\n")


def test_own_error_line_of_the_rule_is_not_the_reason(tmp_path):
    rc, text = call(tmp_path, 1, [
        "error on call: make_nt_freq_mat.pl --section donor_GT\n"
        "ERROR: GeneMark-ET failed to produce genemark.gtf (exit=1)\n",
    ])
    assert rc == 0
    assert "last error message: error on call: make_nt_freq_mat.pl" in text
    assert "failed to produce" not in text


@pytest.mark.parametrize("code", [137, 139, 143])
def test_killed_by_signal_is_not_a_fallback(tmp_path, code):
    rc, text = call(tmp_path, code, ["Illegal division by zero\n"])
    assert rc == 1
    assert text is None


@pytest.mark.parametrize("line", [
    "write failed: No space left on device",
    "sort: write failed: Disk quota exceeded",
    "perl: Cannot allocate memory",
    "Out of memory!",
])
def test_full_disk_or_memory_is_not_a_fallback(tmp_path, line):
    rc, text = call(tmp_path, 1, ["Illegal division by zero\n", line + "\n"])
    assert rc == 1
    assert text is None


def test_tabs_are_removed_and_long_lines_cut(tmp_path):
    rc, text = call(tmp_path, 1, ["error\tin\tfile " + "x" * 400 + "\n"])
    assert rc == 0
    assert "\t" not in text
    reason = text.split("last error message: ")[1].rstrip("\n")
    assert len(reason) == 200


def test_missing_log_file_does_not_break_the_decision(tmp_path):
    marker = tmp_path / "m.txt"
    cmd = 'source "$0"; genemark_fallback_marker "$@"'
    proc = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", cmd, str(SCRIPT), "GeneMark-ET", "1",
         str(marker), str(tmp_path / "does_not_exist.log")],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert marker.read_text().startswith("GeneMark-ET stopped with an error")


@pytest.mark.parametrize("line", [
    "Segmentation fault (core dumped)",
    "Error on last system call, error code 11",
    "Error on last system call, error code 137",
    "sh: line 1: 12345 Segmentation fault      gmhmme3 -m model.mod",
])
def test_crash_of_called_program_is_no_fallback(tmp_path, line):
    """A crash or kill of a program that GeneMark called is a technical
    failure: no marker, return code 1."""
    log = tmp_path / "gm.log"
    log.write_text("some output\n" + line + "\nAbort program!!!\n")
    marker = tmp_path / "mode_selection" / "x_insufficient.txt"
    script = (f'source "{SCRIPT}"; genemark_fallback_marker GeneMark-ET 2 "{marker}" "{log}"')
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    assert result.returncode == 1
    assert not marker.exists()


def test_error_code_of_training_failure_still_falls_back(tmp_path):
    """Other numbers are not crashes: 'error code 1' and the division by zero
    of parse_ET.pl stay a thin-evidence failure."""
    log = tmp_path / "gm.log"
    log.write_text("Illegal division by zero at /opt/ETP/bin/gmes/parse_ET.pl line 105.\n"
                   "Error on last system call, error code 1\nerror code 110 seen\n")
    marker = tmp_path / "mode_selection" / "x_insufficient.txt"
    script = (f'source "{SCRIPT}"; genemark_fallback_marker GeneMark-ET 2 "{marker}" "{log}"')
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    assert result.returncode == 0
    assert marker.exists()
