"""
Regression and positive tests for the single-exon guard added to
merge_features() in scripts/stringtie2utr.py.

Bug (fixed at lines 272-282):
  When a single-exon BRAKER transcript was matched by positional overlap to a
  multi-exon StringTie transcript, distant StringTie exons were incorrectly
  pulled into the BRAKER transcript, extending its coordinates.

All test data is synthetic and generated in-memory.
"""

import sys
import os

# Allow importing the script under test from scripts/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from stringtie2utr import merge_features  # noqa: E402


# ---------------------------------------------------------------------------
# Helper to build a minimal GTF feature line
# ---------------------------------------------------------------------------

def _gtf_line(seqname, source, feature, start, end, score, strand, frame, attrs):
    return "\t".join([seqname, source, feature,
                      str(start), str(end),
                      score, strand, frame, attrs])


# ---------------------------------------------------------------------------
# Test 1 – REGRESSION: single-exon BRAKER gene + multi-exon StringTie
#   transcript whose second exon is far away.
#   Expected: the far exon must NOT be merged; gene bounds stay at 1172-2511.
# ---------------------------------------------------------------------------

def test_single_exon_braker_distant_stringtie_exon_not_merged():
    """Distant second exon of a multi-exon StringTie transcript must not be
    merged into a single-exon BRAKER transcript."""

    braker_tx_id = "g50.t1"
    stringtie_tx_id = "STRG.1.1"

    # BRAKER transcript features (CDS + exon, single-exon gene, Chr1:1172-2511)
    braker_features = [
        _gtf_line("Chr1", "AUGUSTUS", "CDS",   1172, 2508, "0.9", "+", "0",
                  'transcript_id "g50.t1"; gene_id "g50";'),
        _gtf_line("Chr1", "AUGUSTUS", "exon",  1172, 2511, ".",  "+", ".",
                  'transcript_id "g50.t1"; gene_id "g50";'),
        _gtf_line("Chr1", "AUGUSTUS", "stop_codon", 2509, 2511, ".", "+", "0",
                  'transcript_id "g50.t1"; gene_id "g50";'),
    ]

    # StringTie transcript: exon1 overlaps BRAKER exon, exon2 is far away
    stringtie_features = [
        _gtf_line("Chr1", "StringTie", "exon", 1172, 2511, "1000", "+", ".",
                  'gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "5";'),
        _gtf_line("Chr1", "StringTie", "exon", 3001, 6165, "1000", "+", ".",
                  'gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "5";'),
    ]

    tsebra_gtf = {braker_tx_id: list(braker_features)}
    stringtie_gtf = {stringtie_tx_id: list(stringtie_features)}
    selected = {braker_tx_id: stringtie_tx_id}

    result = merge_features(tsebra_gtf, stringtie_gtf, selected)
    merged = result[braker_tx_id]

    # Collect all end coordinates present after merging
    ends = [int(line.split("\t")[4]) for line in merged]
    starts = [int(line.split("\t")[3]) for line in merged]

    assert 6165 not in ends, (
        "Distant StringTie exon end (6165) must NOT be merged into "
        "single-exon BRAKER transcript"
    )
    assert max(ends) == 2511, (
        f"Gene end must remain 2511 after merge, got max end = {max(ends)}"
    )
    assert min(starts) == 1172, (
        f"Gene start must remain 1172 after merge, got min start = {min(starts)}"
    )


# ---------------------------------------------------------------------------
# Test 2 – POSITIVE: single-exon BRAKER gene + single-exon StringTie
#   transcript that overlaps → StringTie exon IS merged.
# ---------------------------------------------------------------------------

def test_single_exon_braker_overlapping_stringtie_exon_is_merged():
    """A StringTie exon that overlaps the BRAKER exon must still be merged."""

    braker_tx_id = "g50.t1"
    stringtie_tx_id = "STRG.2.1"

    braker_features = [
        _gtf_line("Chr1", "AUGUSTUS", "CDS",  1172, 2508, "0.9", "+", "0",
                  'transcript_id "g50.t1"; gene_id "g50";'),
        _gtf_line("Chr1", "AUGUSTUS", "exon", 1172, 2511, ".",  "+", ".",
                  'transcript_id "g50.t1"; gene_id "g50";'),
        _gtf_line("Chr1", "AUGUSTUS", "stop_codon", 2509, 2511, ".", "+", "0",
                  'transcript_id "g50.t1"; gene_id "g50";'),
    ]

    # StringTie exon extends the 5'-UTR region (upstream overlap)
    stringtie_features = [
        _gtf_line("Chr1", "StringTie", "exon", 900, 2511, "1000", "+", ".",
                  'gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "5";'),
    ]

    tsebra_gtf = {braker_tx_id: list(braker_features)}
    stringtie_gtf = {stringtie_tx_id: list(stringtie_features)}
    selected = {braker_tx_id: stringtie_tx_id}

    result = merge_features(tsebra_gtf, stringtie_gtf, selected)
    merged = result[braker_tx_id]

    feature_types = [line.split("\t")[2] for line in merged]
    stringtie_exons = [line for line in merged if "StringTie" in line.split("\t")[1]]

    assert stringtie_exons, (
        "Overlapping StringTie exon must be merged into the single-exon "
        "BRAKER transcript"
    )
    # The merged set must now contain the StringTie exon starting at 900
    starts = [int(line.split("\t")[3]) for line in merged]
    assert 900 in starts, (
        "StringTie exon start (900) must appear in merged features"
    )


# ---------------------------------------------------------------------------
# Test 3 – POSITIVE: multi-exon BRAKER gene matched by intron pattern →
#   ALL StringTie exons are still merged (fix must not regress multi-exon case).
# ---------------------------------------------------------------------------

def test_multi_exon_braker_all_stringtie_exons_merged():
    """For a multi-exon BRAKER transcript the single-exon guard must not
    activate; all StringTie exons are merged as before the fix."""

    braker_tx_id = "g60.t1"
    stringtie_tx_id = "STRG.3.1"

    # Two-exon BRAKER gene: exon1 100-500, exon2 700-1200
    braker_features = [
        _gtf_line("Chr1", "AUGUSTUS", "CDS",  100,  498, "0.9", "+", "0",
                  'transcript_id "g60.t1"; gene_id "g60";'),
        _gtf_line("Chr1", "AUGUSTUS", "exon", 100,  500, ".",  "+", ".",
                  'transcript_id "g60.t1"; gene_id "g60";'),
        _gtf_line("Chr1", "AUGUSTUS", "CDS",  700, 1198, "0.9", "+", "0",
                  'transcript_id "g60.t1"; gene_id "g60";'),
        _gtf_line("Chr1", "AUGUSTUS", "exon", 700, 1200, ".",  "+", ".",
                  'transcript_id "g60.t1"; gene_id "g60";'),
    ]

    # StringTie transcript: same two exons plus a 5'-UTR extension on exon1
    # and a 3'-UTR extension on exon2.
    stringtie_features = [
        _gtf_line("Chr1", "StringTie", "exon",  50,  500, "1000", "+", ".",
                  'gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; cov "5";'),
        _gtf_line("Chr1", "StringTie", "exon", 700, 1400, "1000", "+", ".",
                  'gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "2"; cov "5";'),
    ]

    tsebra_gtf = {braker_tx_id: list(braker_features)}
    stringtie_gtf = {stringtie_tx_id: list(stringtie_features)}
    selected = {braker_tx_id: stringtie_tx_id}

    result = merge_features(tsebra_gtf, stringtie_gtf, selected)
    merged = result[braker_tx_id]

    stringtie_lines = [line for line in merged if "StringTie" in line.split("\t")[1]]

    # Both StringTie exons must have been merged
    assert len(stringtie_lines) == 2, (
        f"Both StringTie exons must be merged for multi-exon BRAKER transcript, "
        f"got {len(stringtie_lines)}"
    )
    ends = [int(line.split("\t")[4]) for line in stringtie_lines]
    assert 500  in ends, "StringTie exon1 (end=500) must be present"
    assert 1400 in ends, "StringTie exon2 (end=1400) must be present"
