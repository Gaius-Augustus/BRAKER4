"""
Unit tests for load_score_cutoff and load_length_cutoff in best_by_compleasm.py.

Two bug fixes are exercised:
  1. load_score_cutoff: taxid key must be the part before 'at', e.g.
     '71915at6073' -> '71915', so hmmsearch query names (which carry no
     lineage suffix) match the dict.
  2. load_length_cutoff: odb12 lengths_cutoff has 3 columns (taxid, length,
     sigma), not 4; reading line[1] for length (was line[3]) must not raise
     IndexError.
"""

import sys
import os
import importlib
import importlib.util
import types

import pytest

# ---------------------------------------------------------------------------
# Import helpers
# ---------------------------------------------------------------------------
# best_by_compleasm.py calls argparser.parse_args() at module level, which
# would fail with "required arguments missing" during a plain import.
# We patch sys.argv with the minimum required positional flags so argparse
# is satisfied, then import the module and restore argv.

_SCRIPT = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "best_by_compleasm.py"
)

_FAKE_ARGV = [
    "best_by_compleasm.py",
    "-m", "/tmp/fake_tmp",
    "-d", "/tmp/fake_input",
    "-g", "/tmp/fake_genome.fa",
    "-p", "insecta_odb12",
]


def _import_module():
    """Import best_by_compleasm under a controlled sys.argv."""
    orig_argv = sys.argv[:]
    sys.argv = _FAKE_ARGV[:]
    try:
        spec = importlib.util.spec_from_file_location(
            "best_by_compleasm", os.path.abspath(_SCRIPT)
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.argv = orig_argv
    return mod


_mod = _import_module()
load_score_cutoff = _mod.load_score_cutoff
load_length_cutoff = _mod.load_length_cutoff


# ---------------------------------------------------------------------------
# Tests for load_score_cutoff
# ---------------------------------------------------------------------------

class TestLoadScoreCutoff:
    def test_odb12_key_strips_at_suffix(self, tmp_path):
        """odb12 entry '71915at6073\t190.0' must yield key '71915'."""
        f = tmp_path / "scores_cutoff"
        f.write_text("71915at6073\t190.0\n")
        result = load_score_cutoff(str(f))
        assert "71915" in result, "key should be '71915' without 'at6073' suffix"
        assert result["71915"] == pytest.approx(190.0)

    def test_odb12_value_is_float(self, tmp_path):
        """Score value must be stored as a float."""
        f = tmp_path / "scores_cutoff"
        f.write_text("12345at9999\t42.5\n")
        result = load_score_cutoff(str(f))
        assert isinstance(result["12345"], float)

    def test_multiple_entries(self, tmp_path):
        """All entries in a multi-line file must be present."""
        f = tmp_path / "scores_cutoff"
        f.write_text("100at6073\t10.0\n200at6073\t20.0\n300at6073\t30.0\n")
        result = load_score_cutoff(str(f))
        assert set(result.keys()) == {"100", "200", "300"}
        assert result["200"] == pytest.approx(20.0)

    def test_entry_without_at_suffix(self, tmp_path):
        """An entry with no 'at' in the taxid must not raise and key is unchanged."""
        f = tmp_path / "scores_cutoff"
        f.write_text("71915\t190.0\n")
        result = load_score_cutoff(str(f))
        # '71915'.split('at')[0] == '71915'
        assert "71915" in result
        assert result["71915"] == pytest.approx(190.0)

    def test_missing_file_raises(self, tmp_path):
        """A non-existent file must raise RuntimeError (not a bare IOError)."""
        with pytest.raises(RuntimeError, match="Impossible to read"):
            load_score_cutoff(str(tmp_path / "nonexistent_scores_cutoff"))


# ---------------------------------------------------------------------------
# Tests for load_length_cutoff
# ---------------------------------------------------------------------------

class TestLoadLengthCutoff:
    def test_odb12_three_column_no_index_error(self, tmp_path):
        """3-column odb12 format must parse without IndexError."""
        f = tmp_path / "lengths_cutoff"
        f.write_text("10002at6073\t181.0\t0.00\n")
        # Would raise IndexError with old line[3]
        result = load_length_cutoff(str(f))
        assert "10002at6073" in result

    def test_sigma_zero_replaced_with_one(self, tmp_path):
        """sigma == 0.0 in source must be stored as 1 (sentinel for 'no sigma')."""
        f = tmp_path / "lengths_cutoff"
        f.write_text("10002at6073\t181.0\t0.00\n")
        result = load_length_cutoff(str(f))
        assert result["10002at6073"]["length"] == pytest.approx(181.0)
        assert result["10002at6073"]["sigma"] == 1

    def test_nonzero_sigma_stored_correctly(self, tmp_path):
        """sigma != 0.0 must be stored as-is."""
        f = tmp_path / "lengths_cutoff"
        f.write_text("20003at6073\t250.0\t12.5\n")
        result = load_length_cutoff(str(f))
        assert result["20003at6073"]["length"] == pytest.approx(250.0)
        assert result["20003at6073"]["sigma"] == pytest.approx(12.5)

    def test_multiple_entries(self, tmp_path):
        """Multiple rows must all be loaded."""
        lines = "aaa\t100.0\t5.0\nbbb\t200.0\t0.0\nccc\t300.0\t15.0\n"
        f = tmp_path / "lengths_cutoff"
        f.write_text(lines)
        result = load_length_cutoff(str(f))
        assert set(result.keys()) == {"aaa", "bbb", "ccc"}
        assert result["bbb"]["sigma"] == 1      # was 0.0, replaced
        assert result["ccc"]["sigma"] == pytest.approx(15.0)

    def test_missing_file_returns_empty_dict(self, tmp_path):
        """A missing lengths_cutoff file must return {} (odb12 lineages omit it)."""
        result = load_length_cutoff(str(tmp_path / "nonexistent_lengths_cutoff"))
        assert result == {}
