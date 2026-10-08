"""
Regression test for merge_features() in scripts/stringtie2utr.py: UTRs on a
side of the CDS are added only if a kept StringTie exon contains that end of
the CDS.

Bug: a StringTie exon that overlaps a CDS segment but is shorter than it is
dropped. When that was the exon over the first (or last) CDS segment, the
StringTie exons further out still became UTRs, joined to the CDS by an intron
StringTie does not have.

All test data is synthetic and generated in-memory.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from stringtie2utr import compute_utr_features, merge_features  # noqa: E402

CDS = [(1000, 1200), (1500, 1700), (2000, 2100)]


def _line(feature, start, end, strand, source, attrs):
    return "\t".join(["Chr1", source, feature, str(start), str(end), ".", strand, ".", attrs])


def _annotate(strand, stringtie_exons):
    attrs = 'transcript_id "g1.t1"; gene_id "g1";'
    braker = []
    for s, e in CDS:
        braker.append(_line("CDS", s, e, strand, "AUGUSTUS", attrs))
        braker.append(_line("exon", s, e, strand, "AUGUSTUS", attrs))
    st_attrs = 'gene_id "STRG.1"; transcript_id "STRG.1.1";'
    stringtie = [_line("exon", s, e, strand, "StringTie", st_attrs) for s, e in stringtie_exons]
    gtf = merge_features({"g1.t1": braker}, {"STRG.1.1": stringtie}, {"g1.t1": "STRG.1.1"})
    gtf = compute_utr_features(gtf)

    def spans(feature):
        return sorted((int(f.split("\t")[3]), int(f.split("\t")[4]))
                      for f in gtf["g1.t1"] if f.split("\t")[2] == feature)
    return spans


def test_no_utr_beyond_a_dropped_first_exon_plus_strand():
    # 1100-1200 is shorter than the CDS segment 1000-1200 and dropped; 500-700
    # would be joined to the CDS by an intron 701-999 StringTie does not have
    spans = _annotate("+", [(500, 700), (1100, 1200), (1500, 1700), (2000, 2400)])
    assert spans("five_prime_UTR") == []
    assert spans("three_prime_UTR") == [(2101, 2400)]


def test_no_utr_beyond_a_dropped_last_exon_minus_strand():
    # 2000-2050 is shorter than the CDS segment 2000-2100 and dropped; on the
    # minus strand 2300-2400 would have become a 5' UTR
    spans = _annotate("-", [(500, 900), (1000, 1200), (1500, 1700), (2000, 2050), (2300, 2400)])
    assert spans("five_prime_UTR") == []
    assert spans("three_prime_UTR") == [(500, 900)]


def test_utrs_on_both_sides_when_the_end_exons_are_kept():
    spans = _annotate("+", [(500, 700), (900, 1200), (1500, 1700), (2000, 2400)])
    assert spans("five_prime_UTR") == [(500, 700), (900, 999)]
    assert spans("three_prime_UTR") == [(2101, 2400)]
