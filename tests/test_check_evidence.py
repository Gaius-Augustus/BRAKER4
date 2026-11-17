"""
Tests for scripts/check_evidence.py (automatic mode fallback).

All data is synthetic and written into tmp_path.
"""

import importlib.util
import os

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "check_evidence.py")
_spec = importlib.util.spec_from_file_location("check_evidence", SCRIPT)
ce = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ce)


def hint(feature="intron", score=".", attr="src=E", seq="chr1", start=10, end=50):
    return "\t".join([seq, "b2h", feature, str(start), str(end), str(score), "+", ".", attr])


def write_hints(path, lines):
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


# ---------------------------------------------------------------- hint_multiplicity

@pytest.mark.parametrize(
    "score, attr, expected",
    [
        ("12", "pri=4;src=E", 12),      # bam2hints style, column 6 only
        (".", "mult=7;src=P", 7),       # attribute only
        ("0", "mult=7", 7),             # score 0 does not hide the attribute
        ("12", "mult=5;src=P", 12),     # larger of both: column 6
        ("3", "pri=4;mult=20", 20),     # larger of both: attribute
        (".", "src=E", 1),              # neither: default 1
        ("0", "pri=4", 1),
        (".", "xmult=50", 1),           # not matched inside another key
        ("2", "src=E;xmult=50", 2),
        (".", "src=E;mult=9", 9),
        (".", "mult=9", 9),             # at start of the attribute column
    ],
)
def test_hint_multiplicity(score, attr, expected):
    fields = hint(score=score, attr=attr).split("\t")
    assert ce.hint_multiplicity(fields) == expected


# ---------------------------------------------------------------- count_introns

def test_count_introns_filters_lines(tmp_path):
    lines = [
        "# a comment",
        hint("intron", 10),                 # supported
        hint("Intron", ".", "mult=4"),      # ProtHint style, case-insensitive
        hint("INTRON", 1),                  # counted, not supported
        hint("start", 50),
        hint("stop", 50),
        hint("CDSpart", 50),
        "chr1\tb2h\tintron\t1\t2\t50\t+\t.",          # 8 columns, skipped
        "chr1\tb2h\tintron\t1\t2\t50\t+\t.\ta\tb",    # 10 columns, skipped
        "",
    ]
    path = write_hints(tmp_path / "h.gff", lines)
    assert ce.count_introns(path, 4) == (3, 2)
    assert ce.count_introns(path, 10) == (3, 1)
    assert ce.count_introns(path, 1) == (3, 3)


def test_count_introns_empty_file(tmp_path):
    path = tmp_path / "empty.gff"
    path.write_text("")
    assert ce.count_introns(str(path), 4) == (0, 0)


# ---------------------------------------------------------------- introns subcommand

def run_introns(tmp_path, n_supported, n_weak, min_introns, min_mult=10):
    lines = [hint(score=min_mult, start=i * 100, end=i * 100 + 50) for i in range(n_supported)]
    lines += [hint(score=min_mult - 1, start=i * 100, end=i * 100 + 50) for i in range(n_weak)]
    hints = write_hints(tmp_path / "hints.gff", lines)
    out = tmp_path / "sub" / "decision.tsv"
    rc = ce.main([
        "introns", "--hints", hints, "--evidence", "RNA-Seq",
        "--min-mult", str(min_mult), "--min-introns", str(min_introns),
        "--out", str(out),
    ])
    return rc, out


def test_introns_exactly_at_threshold_is_sufficient(tmp_path):
    rc, out = run_introns(tmp_path, n_supported=5, n_weak=3, min_introns=5)
    assert rc == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "1"
    assert d["check"] == "introns"
    assert d["evidence"] == "RNA-Seq"
    assert d["introns"] == "8"
    assert d["introns_supported"] == "5"
    assert d["min_mult"] == "10"
    assert d["min_introns"] == "5"
    assert "5 RNA-Seq intron hints" in d["reason"]


def test_introns_one_below_threshold_is_insufficient(tmp_path):
    rc, out = run_introns(tmp_path, n_supported=4, n_weak=3, min_introns=5)
    assert rc == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "0"
    assert d["introns_supported"] == "4"
    assert d["introns"] == "7"


def test_introns_empty_hints(tmp_path):
    rc, out = run_introns(tmp_path, n_supported=0, n_weak=0, min_introns=5)
    assert rc == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "0"
    assert d["introns"] == "0"
    assert d["introns_supported"] == "0"


def test_introns_missing_hints_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ce.main(["introns", "--hints", str(tmp_path / "nope.gff"), "--evidence", "x",
                 "--min-mult", "4", "--min-introns", "1", "--out", str(tmp_path / "o.tsv")])


# ------------------------------------------- result subcommand (alias: etp)

def run_etp(tmp_path, marker_text=None, gtf_text=None):
    gtf = tmp_path / "genemark.gtf"
    marker = tmp_path / "ETP_NOT_ENOUGH"
    out = tmp_path / "etp_decision.tsv"
    if gtf_text is not None:
        gtf.write_text(gtf_text)
    if marker_text is not None:
        marker.write_text(marker_text)
    rc = ce.main(["etp", "--gtf", str(gtf), "--marker", str(marker), "--out", str(out)])
    return rc, out


def test_etp_marker_means_insufficient_with_marker_text_as_reason(tmp_path):
    rc, out = run_etp(tmp_path, marker_text="no spliced transcripts\n\n  too few genes  \n",
                      gtf_text="chr1\tx\tCDS\t1\t3\t.\t+\t0\tg;\n")
    assert rc == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "0"
    assert d["check"] == "result"
    assert d["evidence"] == "GeneMark-ETP"
    assert d["reason"] == "no spliced transcripts too few genes"


def test_etp_empty_marker_gets_default_reason(tmp_path):
    rc, out = run_etp(tmp_path, marker_text="")
    assert rc == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "0"
    assert d["reason"] == "GeneMark-ETP stopped for lack of training data"


def test_etp_nonempty_gtf_without_marker_is_sufficient(tmp_path):
    rc, out = run_etp(tmp_path, gtf_text="chr1\tx\tCDS\t1\t3\t.\t+\t0\tg;\n")
    assert rc == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "1"
    assert "finished" in d["reason"]


@pytest.mark.parametrize("gtf_text", [None, ""])
def test_etp_empty_or_missing_gtf_without_marker_fails(tmp_path, gtf_text, capsys):
    rc, out = run_etp(tmp_path, gtf_text=gtf_text)
    assert rc == 1
    assert not out.exists()
    assert "did not finish" in capsys.readouterr().err


@pytest.mark.parametrize("tool", ["GeneMark-ET", "GeneMark-EP"])
def test_result_names_the_tool(tmp_path, tool, capsys):
    gtf = tmp_path / "genemark.gtf"
    marker = tmp_path / "marker.txt"
    out = tmp_path / "decision.tsv"
    # empty marker: default reason names the tool
    marker.write_text("")
    assert ce.main(["result", "--tool", tool, "--gtf", str(gtf),
                    "--marker", str(marker), "--out", str(out)]) == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "0"
    assert d["evidence"] == tool
    assert d["reason"] == f"{tool} stopped for lack of training data"
    # no marker, non-empty gtf: enough
    marker.unlink()
    gtf.write_text("chr1\tx\tCDS\t1\t3\t.\t+\t0\tg;\n")
    assert ce.main(["result", "--tool", tool, "--gtf", str(gtf),
                    "--marker", str(marker), "--out", str(out)]) == 0
    d = ce.read_decision(str(out))
    assert d["sufficient"] == "1"
    assert d["reason"] == f"{tool} finished and wrote genemark.gtf"
    # no marker, empty gtf: the run did not finish
    gtf.write_text("")
    out.unlink()
    assert ce.main(["result", "--tool", tool, "--gtf", str(gtf),
                    "--marker", str(marker), "--out", str(out)]) == 1
    assert not out.exists()
    assert f"{tool} did not finish" in capsys.readouterr().err


# ---------------------------------------------------------------- write_decision

def test_write_decision_sanitises_tabs_and_newlines(tmp_path):
    path = tmp_path / "deep" / "d.tsv"
    ce.write_decision(str(path), [("a", "x\ty\nz"), ("b", 3), ("reason", "ok")])
    lines = path.read_text().split("\n")
    assert lines == ["a\tx y z", "b\t3", "reason\tok", ""]
    assert all(l.count("\t") == 1 for l in lines[:-1])
    assert ce.read_decision(str(path)) == {"a": "x y z", "b": "3", "reason": "ok"}
