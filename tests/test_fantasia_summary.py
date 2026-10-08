"""
Tests for scripts/fantasia_summary.py: fantasia_go_categories.png is a
required output of rule fantasia_summarize and must exist even when no GO
term passes --min-score (a placeholder figure, not a missing file).

All test data is synthetic and generated in-memory.
"""

import os
import subprocess
import sys

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "fantasia_summary.py")

FANTASIA_HEADER = "query_accession,go_id,reliability_index,distance,go_description,category\n"


@pytest.mark.parametrize("rows", [
    "",                                                          # results.csv with a header only
    "g1.t1,GO:0005524,0.21,0.9,ATP binding,molecular_function\n"  # every term below the cutoff
    "g2.t1,GO:0005634,0.40,0.8,nucleus,cellular_component\n",
], ids=["no-rows", "all-below-cutoff"])
def test_placeholder_png_when_no_go_term_passes(tmp_path, rows):
    results = tmp_path / "results.csv"
    results.write_text(FANTASIA_HEADER + rows)
    out = tmp_path / "out"
    proc = subprocess.run([sys.executable, SCRIPT, "--results", str(results), "--out-dir", str(out),
                           "--min-score", "0.5"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    png = out / "fantasia_go_categories.png"
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert (out / "fantasia_summary.txt").is_file()
    assert (out / "fantasia_go_terms.tsv").read_text() == \
        "transcript_id\tgo_id\tgo_name\tgo_namespace\treliability_index\n"
