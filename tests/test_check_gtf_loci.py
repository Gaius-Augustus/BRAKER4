"""
Tests for scripts/check_gtf_loci.py (issue #62).

Every gene and transcript must sit on one sequence and strand, and gene and
transcript lines must span exactly their children. Violations must fail the
script with exit code 1.

All test data is synthetic and generated in-memory.
"""

import os
import subprocess
import sys

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
CHECK = os.path.join(SCRIPTS, "check_gtf_loci.py")


def _line(seq, feature, start, end, strand, attr, source="AUGUSTUS"):
    frame = "0" if feature == "CDS" else "."
    return "\t".join([seq, source, feature, str(start), str(end), ".",
                      strand, frame, attr]) + "\n"


def _feat(seq, feature, start, end, strand, tx="g1.t1", gene="g1"):
    return _line(seq, feature, start, end, strand,
                 f'transcript_id "{tx}"; gene_id "{gene}";')


# Bare TSEBRA/AUGUSTUS gene and transcript lines, as in braker.gtf
GOOD_BARE = (
    _line("X1", "gene", 1000, 3000, "+", "g1")
    + _line("X1", "transcript", 1000, 3000, "+", "g1.t1")
    + _feat("X1", "exon", 1000, 1500, "+")
    + _feat("X1", "CDS", 1000, 1500, "+")
    + _feat("X1", "exon", 2000, 3000, "+")
    + _feat("X1", "CDS", 2000, 3000, "+")
)

# Attribute-style gene and transcript lines, as written by stringtie2utr.py
GOOD_ATTR = (
    _line("X1", "gene", 900, 3200, "-", 'gene_id "g1";')
    + _line("X1", "transcript", 900, 3200, "-", 'transcript_id "g1.t1"; gene_id "g1";')
    + _feat("X1", "three_prime_UTR", 900, 999, "-")
    + _feat("X1", "CDS", 1000, 3000, "-")
    + _feat("X1", "five_prime_UTR", 3001, 3200, "-")
    + _line("X1", "transcript", 1000, 2500, "-", 'transcript_id "g1.t2"; gene_id "g1";')
    + _feat("X1", "CDS", 1000, 2500, "-", tx="g1.t2")
)


def _run(tmp_path, text):
    gtf = tmp_path / "in.gtf"
    gtf.write_text(text)
    return subprocess.run([sys.executable, CHECK, str(gtf)],
                          capture_output=True, text=True)


def test_consistent_bare_lines_pass(tmp_path):
    res = _run(tmp_path, GOOD_BARE)
    assert res.returncode == 0, res.stderr


def test_consistent_attribute_lines_pass(tmp_path):
    res = _run(tmp_path, GOOD_ATTR)
    assert res.returncode == 0, res.stderr


def test_transcript_on_two_sequences_fails(tmp_path):
    text = GOOD_BARE + _feat("X2", "five_prime_UTR", 500, 800, "+")
    res = _run(tmp_path, text)
    assert res.returncode == 1
    assert "g1.t1" in res.stderr and "several sequences" in res.stderr


def test_transcript_on_two_strands_fails(tmp_path):
    text = GOOD_BARE + _feat("X1", "three_prime_UTR", 3001, 3100, "-")
    res = _run(tmp_path, text)
    assert res.returncode == 1
    assert "g1.t1" in res.stderr


def test_transcript_line_span_mismatch_fails(tmp_path):
    text = GOOD_BARE.replace("transcript\t1000\t3000", "transcript\t1000\t9000")
    res = _run(tmp_path, text)
    assert res.returncode == 1
    assert "transcript g1.t1" in res.stderr


def test_gene_span_mismatch_fails(tmp_path):
    # The #62 symptom: gene stretched far beyond its only transcript
    text = GOOD_BARE.replace("gene\t1000\t3000", "gene\t1000\t880000")
    res = _run(tmp_path, text)
    assert res.returncode == 1
    assert "gene g1" in res.stderr


def test_gene_on_other_sequence_than_transcripts_fails(tmp_path):
    text = GOOD_BARE.replace("X1\tAUGUSTUS\tgene", "X2\tAUGUSTUS\tgene")
    res = _run(tmp_path, text)
    assert res.returncode == 1
    assert "gene g1" in res.stderr


def test_duplicate_gene_line_fails(tmp_path):
    text = GOOD_BARE + _line("X2", "gene", 1000, 3000, "+", "g1")
    res = _run(tmp_path, text)
    assert res.returncode == 1
    assert "2 gene lines" in res.stderr
