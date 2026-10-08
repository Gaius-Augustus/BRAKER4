"""
Regression tests for the stop codon handling in compute_utr_features() of
scripts/stringtie2utr.py.

Bug: BRAKER GTF follows GTF2.2, CDS lines exclude the stop codon, which has
its own stop_codon line. The UTR computation took the CDS span as the coding
region, so the stop codon (inside the exon, outside the CDS) became a 3-bp
three_prime_UTR on every transcript, identical to the stop_codon feature.

The stop codon is part of the coding region: a 3' UTR starts after it. The
output keeps the GTF2.2 convention (CDS without stop codon, stop_codon line).

All test data is synthetic and generated in-memory.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from stringtie2utr import compute_utr_features, merge_features  # noqa: E402


def _line(feature, start, end, strand, tx, source="AUGUSTUS"):
    return "\t".join(["Chr1", source, feature, str(start), str(end), ".",
                      strand, ".", f'transcript_id "{tx}"; gene_id "{tx}_g";'])


def _st_exon(start, end, strand, st_tx="STRG.1.1"):
    return "\t".join(["Chr1", "StringTie", "exon", str(start), str(end), "1000",
                      strand, ".",
                      f'gene_id "STRG.1"; transcript_id "{st_tx}"; cov "5";'])


def _coords(features, ftype):
    return sorted((int(f.split("\t")[3]), int(f.split("\t")[4]))
                  for f in features if f.split("\t")[2] == ftype)


def _run(braker, stringtie=None, tx="t1"):
    gtf = {tx: list(braker)}
    if stringtie:
        gtf = merge_features(gtf, {"STRG.1.1": stringtie}, {tx: "STRG.1.1"})
    return compute_utr_features(gtf)[tx]


def _no_utr_on_coding(features):
    coding = _coords(features, "CDS") + _coords(features, "stop_codon")
    for ftype in ("five_prime_UTR", "three_prime_UTR"):
        for us, ue in _coords(features, ftype):
            for cs, ce in coding:
                assert not (us <= ce and ue >= cs), (
                    f"{ftype} {us}-{ue} overlaps coding {cs}-{ce}")


def test_plus_single_exon_stop_codon_is_not_utr():
    f = _run([
        _line("exon", 101, 200, "+", "t1"),
        _line("CDS", 101, 197, "+", "t1"),
        _line("start_codon", 101, 103, "+", "t1"),
        _line("stop_codon", 198, 200, "+", "t1"),
    ])
    assert _coords(f, "three_prime_UTR") == []
    assert _coords(f, "five_prime_UTR") == []
    assert _coords(f, "stop_codon") == [(198, 200)]
    assert _coords(f, "CDS") == [(101, 197)]


def test_minus_single_exon_stop_codon_is_not_utr():
    # the g1.t1 pattern from the scenario 11 report
    f = _run([
        _line("exon", 1, 1224, "-", "t1"),
        _line("CDS", 4, 1224, "-", "t1"),
        _line("stop_codon", 1, 3, "-", "t1"),
        _line("start_codon", 1222, 1224, "-", "t1"),
    ])
    assert _coords(f, "three_prime_UTR") == []
    assert _coords(f, "five_prime_UTR") == []


def test_minus_strand_real_utr_starts_after_stop_codon():
    braker = [
        _line("exon", 1001, 1224, "-", "t1"),
        _line("CDS", 1004, 1224, "-", "t1"),
        _line("stop_codon", 1001, 1003, "-", "t1"),
        _line("start_codon", 1222, 1224, "-", "t1"),
    ]
    f = _run(braker, [_st_exon(901, 1300, "-")])
    assert _coords(f, "three_prime_UTR") == [(901, 1000)]
    assert _coords(f, "five_prime_UTR") == [(1225, 1300)]
    _no_utr_on_coding(f)


def test_plus_strand_real_utr_starts_after_stop_codon():
    braker = [
        _line("exon", 101, 200, "+", "t1"),
        _line("CDS", 101, 197, "+", "t1"),
        _line("start_codon", 101, 103, "+", "t1"),
        _line("stop_codon", 198, 200, "+", "t1"),
    ]
    f = _run(braker, [_st_exon(51, 260, "+")])
    assert _coords(f, "five_prime_UTR") == [(51, 100)]
    assert _coords(f, "three_prime_UTR") == [(201, 260)]
    _no_utr_on_coding(f)


def test_plus_multi_exon_stop_codon_split_by_intron():
    # CDS ends at 198, stop codon = 199-200 + 301, last exon 301-400
    braker = [
        _line("exon", 101, 200, "+", "t1"),
        _line("exon", 301, 400, "+", "t1"),
        _line("CDS", 101, 198, "+", "t1"),
        _line("start_codon", 101, 103, "+", "t1"),
        _line("stop_codon", 199, 200, "+", "t1"),
        _line("stop_codon", 301, 301, "+", "t1"),
    ]
    f = _run(braker)
    # the exon remainder after the stop codon is 3' UTR
    assert _coords(f, "three_prime_UTR") == [(302, 400)]
    assert _coords(f, "five_prime_UTR") == []
    _no_utr_on_coding(f)


def test_minus_multi_exon_stop_codon_split_by_intron():
    # stop codon = 301 (upstream exon) + 199-200 (downstream exon 101-200)
    braker = [
        _line("exon", 101, 200, "-", "t1"),
        _line("exon", 301, 400, "-", "t1"),
        _line("CDS", 302, 400, "-", "t1"),
        _line("start_codon", 398, 400, "-", "t1"),
        _line("stop_codon", 301, 301, "-", "t1"),
        _line("stop_codon", 199, 200, "-", "t1"),
    ]
    f = _run(braker)
    assert _coords(f, "three_prime_UTR") == [(101, 198)]
    assert _coords(f, "five_prime_UTR") == []
    _no_utr_on_coding(f)


def test_multi_exon_no_extension_has_no_utr():
    braker = [
        _line("exon", 101, 200, "+", "t1"),
        _line("exon", 301, 400, "+", "t1"),
        _line("CDS", 101, 200, "+", "t1"),
        _line("CDS", 301, 397, "+", "t1"),
        _line("start_codon", 101, 103, "+", "t1"),
        _line("stop_codon", 398, 400, "+", "t1"),
    ]
    f = _run(braker)
    assert _coords(f, "three_prime_UTR") == []
    assert _coords(f, "five_prime_UTR") == []
