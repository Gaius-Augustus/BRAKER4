"""
Tests for scripts/eval_test_set_accuracy.py.

Covers GenBank location parsing (single, join, complement, complement(join),
line-wrapped, partial markers, underscore-containing seqnames), the LOCUS
length sanity check, stop_codon merging, strand-aware nucleotide scoring,
isoform collapsing at the gene level, clipping of predictions to the test
regions, and that format_report() is parseable by
training_summary.parse_accuracy_file().

All test data is synthetic and written to tmp_path only.
"""

import os
import subprocess
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, SCRIPTS)

from eval_test_set_accuracy import (  # noqa: E402
    evaluate,
    format_report,
    parse_genbank_test_set,
    parse_gtf_cds,
    parse_location,
)
from training_summary import parse_accuracy_file  # noqa: E402

SCRIPT = os.path.join(SCRIPTS, "eval_test_set_accuracy.py")


# ---------------------------------------------------------------------------
# GenBank file builder helpers (fixed-width feature table, like
# gff2gbSmallDNA.pl produces)
# ---------------------------------------------------------------------------

def _locus(seq, start, end):
    length = end - start + 1
    return f"LOCUS       {seq}_{start}-{end}       {length} bp    DNA     linear   UNK\n"


def _feat(key, loc):
    return f"     {key:<16}{loc}\n"


def _cont(text):
    return " " * 21 + text + "\n"


def gb_region(seq, start, end, cds_locations):
    """One GenBank LOCUS entry.

    cds_locations: list of location strings (one CDS feature each) or, for a
    line-wrapped location, a list of the raw pieces that appear on
    successive lines of the same CDS feature.
    """
    length = end - start + 1
    out = [
        _locus(seq, start, end),
        "DEFINITION  synthetic test region.\n",
        "FEATURES             Location/Qualifiers\n",
        _feat("source", f"1..{length}"),
        _cont('/organism="synthetic"'),
    ]
    for loc in cds_locations:
        pieces = loc if isinstance(loc, list) else [loc]
        out.append(_feat("CDS", pieces[0]))
        for piece in pieces[1:]:
            out.append(_cont(piece))
        out.append(_cont('/gene="g"'))
        out.append(_cont('/product="hypothetical protein"'))
    out.append("ORIGIN\n")
    out.append("//\n")
    return "".join(out)


def gb_file(tmp_path, entries, name="test.gb"):
    """entries: list of (seq, start, end, cds_locations)."""
    text = "".join(gb_region(*e) for e in entries)
    path = tmp_path / name
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
# GTF builder helpers
# ---------------------------------------------------------------------------

def _attr(gid, tid):
    return f'gene_id "{gid}"; transcript_id "{tid}";'


def cds_line(seq, start, end, strand, gid, tid, feature="CDS", frame="0", source="AUGUSTUS"):
    return f"{seq}\t{source}\t{feature}\t{start}\t{end}\t.\t{strand}\t{frame}\t{_attr(gid, tid)}\n"


def gtf_file(tmp_path, lines, name="pred.gtf"):
    path = tmp_path / name
    path.write_text("".join(lines))
    return path


# ---------------------------------------------------------------------------
# parse_location: single, join, complement, complement(join), partial markers
# ---------------------------------------------------------------------------

def test_parse_location_single():
    assert parse_location("101..200") == ("+", [(101, 200)])


def test_parse_location_join():
    assert parse_location("join(101..200,301..400)") == ("+", [(101, 200), (301, 400)])


def test_parse_location_complement():
    assert parse_location("complement(101..200)") == ("-", [(101, 200)])


def test_parse_location_complement_join():
    assert parse_location("complement(join(101..200,301..400))") == (
        "-",
        [(101, 200), (301, 400)],
    )


def test_parse_location_partial_markers_stripped():
    # < and > partial-feature markers must not end up in the coordinates
    assert parse_location("join(<101..200,301..>400)") == ("+", [(101, 200), (301, 400)])
    assert parse_location("complement(join(<1..50,>951..1000))") == (
        "-",
        [(1, 50), (951, 1000)],
    )


def test_parse_location_sorts_pieces_regardless_of_file_order():
    # join() lists pieces out of order in some tools; the parser must sort
    assert parse_location("join(301..400,101..200)") == ("+", [(101, 200), (301, 400)])


# ---------------------------------------------------------------------------
# parse_genbank_test_set: seqnames with underscores, line-wrapped locations,
# LOCUS length mismatch
# ---------------------------------------------------------------------------

def test_seqname_with_underscores(tmp_path):
    path = gb_file(tmp_path, [("scaffold_007", 501, 800, ["101..200"])])
    regions = parse_genbank_test_set(path)
    assert len(regions) == 1
    assert regions[0]["seq"] == "scaffold_007"
    assert regions[0]["start"] == 501
    assert regions[0]["end"] == 800
    # 101..200 relative, offset = start - 1 = 500
    assert regions[0]["genes"] == [("+", ((601, 700),))]


def test_multiline_wrapped_location_matches_single_line_equivalent(tmp_path):
    wrapped = gb_file(
        tmp_path,
        [("chrW", 1001, 2000, [["join(101..200,", "301..403)"]])],
        name="wrapped.gb",
    )
    flat = gb_file(
        tmp_path,
        [("chrW", 1001, 2000, ["join(101..200,301..403)"])],
        name="flat.gb",
    )
    assert parse_genbank_test_set(wrapped) == parse_genbank_test_set(flat)
    regions = parse_genbank_test_set(wrapped)
    assert regions[0]["genes"] == [("+", ((1101, 1200), (1301, 1403)))]


def test_locus_length_mismatch_raises(tmp_path):
    path = tmp_path / "bad.gb"
    # claims 999 bp, actual end - start + 1 = 1000
    text = (
        "LOCUS       chr1_1001-2000       999 bp    DNA     linear   UNK\n"
        "DEFINITION  synthetic\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..1000\n"
        "ORIGIN\n"
        "//\n"
    )
    path.write_text(text)
    with pytest.raises(ValueError):
        parse_genbank_test_set(path)


# ---------------------------------------------------------------------------
# Scoring behaviour via evaluate()/parse_gtf_cds()
# ---------------------------------------------------------------------------

def test_stop_codon_merged_gives_exact_gene_match(tmp_path):
    # GenBank CDS includes the stop codon: join(101..200,301..403) relative,
    # offset 1000 -> genome (1101,1200),(1301,1403)
    gb = gb_file(tmp_path, [("chr1", 1001, 2000, ["join(101..200,301..403)"])])
    regions = parse_genbank_test_set(gb)

    # AUGUSTUS-style GTF keeps the stop codon as a separate feature
    gtf = gtf_file(
        tmp_path,
        [
            cds_line("chr1", 1101, 1200, "+", "g1", "g1.t1"),
            cds_line("chr1", 1301, 1400, "+", "g1", "g1.t1"),
            cds_line("chr1", 1401, 1403, "+", "g1", "g1.t1", feature="stop_codon"),
        ],
    )
    predictions = parse_gtf_cds(gtf)
    c = evaluate(regions, predictions)

    assert c["nu_tp"] == c["nu_anno"] == c["nu_pred"] == 203
    assert c["ex_tp"] == c["ex_anno"] == c["ex_pred"] == 2
    assert c["gene_anno"] == c["gene_tp_anno"] == 1
    assert c["gene_pred"] == c["gene_tp_pred"] == 1


def test_wrong_strand_prediction_gives_zero_nucleotide_tp(tmp_path):
    gb = gb_file(tmp_path, [("chrW", 1001, 1500, ["101..400"])])
    regions = parse_genbank_test_set(gb)
    assert regions[0]["genes"] == [("+", ((1101, 1400),))]

    # identical coordinates, but predicted on the opposite strand
    gtf = gtf_file(tmp_path, [cds_line("chrW", 1101, 1400, "-", "p1", "p1.t1")])
    predictions = parse_gtf_cds(gtf)
    c = evaluate(regions, predictions)

    assert c["nu_tp"] == 0
    assert c["nu_anno"] == 300
    assert c["nu_pred"] == 300
    assert c["gene_tp_anno"] == 0
    assert c["gene_tp_pred"] == 0


def test_two_isoforms_one_matching_counts_as_one_correct_predicted_gene(tmp_path):
    gb = gb_file(tmp_path, [("chrG", 2001, 3000, ["join(101..200,301..403)"])])
    regions = parse_genbank_test_set(gb)
    assert regions[0]["genes"] == [("+", ((2101, 2200), (2301, 2403)))]

    gtf = gtf_file(
        tmp_path,
        [
            # t1: exact match (stop codon merges back to 2301..2403)
            cds_line("chrG", 2101, 2200, "+", "g1", "g1.t1"),
            cds_line("chrG", 2301, 2400, "+", "g1", "g1.t1"),
            cds_line("chrG", 2401, 2403, "+", "g1", "g1.t1", feature="stop_codon"),
            # t2: same gene, different (non-matching) single-exon structure
            cds_line("chrG", 2101, 2403, "+", "g1", "g1.t2"),
        ],
    )
    predictions = parse_gtf_cds(gtf)
    c = evaluate(regions, predictions)

    assert c["gene_pred"] == 1  # one predicted gene (g1), not one per isoform
    assert c["gene_tp_pred"] == 1  # matched, because t1 matches exactly
    assert c["gene_anno"] == 1
    assert c["gene_tp_anno"] == 1


def test_predictions_outside_all_regions_are_ignored(tmp_path):
    gb = gb_file(tmp_path, [("chrQ", 1001, 1500, ["101..400"])])
    regions = parse_genbank_test_set(gb)

    gtf = gtf_file(
        tmp_path,
        [
            # different sequence entirely, not among the test regions
            cds_line("chrZ", 1101, 1400, "+", "z1", "z1.t1"),
            # same sequence as the region, but entirely outside its bounds
            cds_line("chrQ", 2000, 2100, "+", "q1", "q1.t1"),
        ],
    )
    predictions = parse_gtf_cds(gtf)
    c = evaluate(regions, predictions)

    assert c["gene_pred"] == 0
    assert c["gene_tp_pred"] == 0
    assert c["nu_pred"] == 0
    assert c["ex_pred"] == 0


def test_partially_overlapping_transcript_is_clipped_and_counts_as_predicted_not_matching(tmp_path):
    gb = gb_file(tmp_path, [("chrP", 1001, 1500, ["101..500"])])
    regions = parse_genbank_test_set(gb)
    assert regions[0]["genes"] == [("+", ((1101, 1500),))]

    # exon runs from 901 to 1200: half outside the region, half overlapping
    # the annotated CDS
    gtf = gtf_file(tmp_path, [cds_line("chrP", 901, 1200, "+", "p1", "p1.t1")])
    predictions = parse_gtf_cds(gtf)
    c = evaluate(regions, predictions)

    # clipped to (1001, 1200): counted as a predicted gene/exon...
    assert c["gene_pred"] == 1
    assert c["ex_pred"] == 1
    assert c["nu_pred"] == 200  # 1001..1200
    # ...but it does not reproduce the annotated (1101, 1500) exon
    assert c["gene_tp_pred"] == 0
    assert c["ex_tp"] == 0
    # overlap of clipped (1001,1200) with annotated (1101,1500) is 1101..1200
    assert c["nu_tp"] == 100
    assert c["nu_anno"] == 400


# ---------------------------------------------------------------------------
# format_report() output must be parseable by training_summary
# ---------------------------------------------------------------------------

def test_format_report_parseable_by_training_summary(tmp_path):
    gb = gb_file(tmp_path, [("chr1", 1001, 2000, ["join(101..200,301..403)"])])
    regions = parse_genbank_test_set(gb)
    gtf = gtf_file(
        tmp_path,
        [
            cds_line("chr1", 1101, 1200, "+", "g1", "g1.t1"),
            cds_line("chr1", 1301, 1400, "+", "g1", "g1.t1"),
            cds_line("chr1", 1401, 1403, "+", "g1", "g1.t1", feature="stop_codon"),
        ],
    )
    predictions = parse_gtf_cds(gtf)
    c = evaluate(regions, predictions)
    report = format_report(c, "braker.gtf")

    out = tmp_path / "accuracy_final_gene_set.txt"
    out.write_text(report)

    acc = parse_accuracy_file(str(out))
    for key in ("nu_sen", "nu_sp", "ex_sen", "ex_sp", "gene_sen", "gene_sp", "weighted"):
        assert key in acc
    # exact match scenario -> everything is 100%
    for key in acc:
        assert acc[key] == pytest.approx(100.0)


def test_cli_end_to_end(tmp_path):
    gb = gb_file(tmp_path, [("chr1", 1001, 2000, ["join(101..200,301..403)"])])
    gtf = gtf_file(
        tmp_path,
        [
            cds_line("chr1", 1101, 1200, "+", "g1", "g1.t1"),
            cds_line("chr1", 1301, 1400, "+", "g1", "g1.t1"),
            cds_line("chr1", 1401, 1403, "+", "g1", "g1.t1", feature="stop_codon"),
        ],
    )
    out = tmp_path / "accuracy_final_gene_set.txt"

    result = subprocess.run(
        [sys.executable, SCRIPT, "-t", str(gb), "-g", str(gtf), "-o", str(out)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.exists()

    acc = parse_accuracy_file(str(out))
    assert acc["gene_sen"] == pytest.approx(100.0)
    assert acc["gene_sp"] == pytest.approx(100.0)
