"""
Regression tests for scripts/stringtie2utr.py:

1. A StringTie exon that starts inside the CDS and ends behind it (single-exon
   genes) was kept and renamed three_prime_UTR as a whole, next to the correct
   three_prime_UTR cut from it: a second 3' UTR running from inside the CDS to
   the transcript end (10 transcripts in a scenario 11 run).
2. Exon lines ended at the coding region; the UTRs sat outside every exon.
   In GTF the exon covers the UTR: rebuild_exons() extends the exon lines and
   gives a UTR exon behind an intron its own exon and intron line.
3. GeneMark transcripts in braker.gtf without a stop_codon line, whose exon
   ends 3 bp behind the CDS: those 3 bp are the stop codon, not a UTR.
4. Intron chain matching kept one transcript per intron, so introns shared
   by alternative transcripts were lost; a transcript then matched a StringTie
   model lacking some of its introns and got UTRs inside its introns.

All test data is synthetic and generated in-memory.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from stringtie2utr import (compute_utr_features, merge_features, rebuild_exons,  # noqa: E402
                           coding_intervals, create_introns_hash,
                           find_matching_transcripts, print_gtf)


def _line(feature, start, end, strand, tx, source="AUGUSTUS", score="."):
    return "\t".join(["Chr1", source, feature, str(start), str(end), score,
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
    return rebuild_exons(compute_utr_features(gtf))[tx]


def _utrs_inside_exons(features):
    exons = _coords(features, "exon")
    for ftype in ("five_prime_UTR", "three_prime_UTR"):
        for us, ue in _coords(features, ftype):
            assert any(es <= us and ue <= ee for es, ee in exons), (
                f"{ftype} {us}-{ue} outside every exon {exons}")


def _no_utr_on_coding(features):
    coding = _coords(features, "CDS") + _coords(features, "stop_codon")
    for ftype in ("five_prime_UTR", "three_prime_UTR"):
        for us, ue in _coords(features, ftype):
            for cs, ce in coding:
                assert not (us <= ce and ue >= cs), (
                    f"{ftype} {us}-{ue} overlaps coding {cs}-{ce}")


PLUS_SINGLE = [
    _line("exon", 1001, 2200, "+", "t1"),
    _line("CDS", 1001, 2197, "+", "t1"),
    _line("start_codon", 1001, 1003, "+", "t1"),
    _line("stop_codon", 2198, 2200, "+", "t1"),
]

MINUS_SINGLE = [
    _line("exon", 1001, 2200, "-", "t1"),
    _line("CDS", 1004, 2200, "-", "t1"),
    _line("stop_codon", 1001, 1003, "-", "t1"),
    _line("start_codon", 2198, 2200, "-", "t1"),
]


# --- 1. StringTie exon starting inside the CDS ------------------------------

def test_plus_stringtie_exon_starting_inside_cds_gives_one_3prime_utr():
    # the g3877.t2 pattern: StringTie exon 17 bp into the CDS, ends behind it
    f = _run(PLUS_SINGLE, [_st_exon(1018, 2300, "+")])
    assert _coords(f, "three_prime_UTR") == [(2201, 2300)]
    assert _coords(f, "five_prime_UTR") == []
    assert not [x for x in f if "StringTie" in x.split("\t")[1]
                and x.split("\t")[2] == "exon"]
    _no_utr_on_coding(f)


def test_minus_stringtie_exon_ending_inside_cds_gives_one_5prime_utr():
    f = _run(MINUS_SINGLE, [_st_exon(1020, 2300, "-")])
    assert _coords(f, "five_prime_UTR") == [(2201, 2300)]
    assert _coords(f, "three_prime_UTR") == []
    _no_utr_on_coding(f)


def test_minus_stringtie_exon_starting_inside_cds_gives_one_3prime_utr():
    f = _run(MINUS_SINGLE, [_st_exon(901, 2190, "-")])
    assert _coords(f, "three_prime_UTR") == [(901, 1000)]
    assert _coords(f, "five_prime_UTR") == []
    _no_utr_on_coding(f)


# --- 2. exon lines cover the UTRs -------------------------------------------

def test_single_exon_utr_both_sides_extends_exon():
    f = _run(PLUS_SINGLE, [_st_exon(901, 2300, "+")])
    assert _coords(f, "five_prime_UTR") == [(901, 1000)]
    assert _coords(f, "three_prime_UTR") == [(2201, 2300)]
    assert _coords(f, "exon") == [(901, 2300)]
    assert _coords(f, "CDS") == [(1001, 2197)]
    assert _coords(f, "stop_codon") == [(2198, 2200)]
    _utrs_inside_exons(f)
    # the exon line keeps the predictor's source and attributes
    exon = [x for x in f if x.split("\t")[2] == "exon"][0].split("\t")
    assert exon[1] == "AUGUSTUS"
    assert 'transcript_id "t1"' in exon[8]


def test_utr_exon_behind_intron_gets_exon_and_intron_line():
    braker = [
        _line("exon", 1001, 1200, "+", "t1"),
        _line("exon", 1301, 1500, "+", "t1"),
        _line("CDS", 1001, 1200, "+", "t1"),
        _line("CDS", 1301, 1497, "+", "t1"),
        _line("intron", 1201, 1300, "+", "t1", score="0.9"),
        _line("start_codon", 1001, 1003, "+", "t1"),
        _line("stop_codon", 1498, 1500, "+", "t1"),
    ]
    # 5' UTR exon 701-800 behind an intron, 5' UTR 901-1000 joined to exon 1
    st = [_st_exon(701, 800, "+"), _st_exon(901, 1200, "+"), _st_exon(1301, 1600, "+")]
    f = _run(braker, st)
    assert _coords(f, "five_prime_UTR") == [(701, 800), (901, 1000)]
    assert _coords(f, "three_prime_UTR") == [(1501, 1600)]
    assert _coords(f, "exon") == [(701, 800), (901, 1200), (1301, 1600)]
    assert _coords(f, "intron") == [(801, 900), (1201, 1300)]
    _utrs_inside_exons(f)
    introns = {(int(x.split("\t")[3]), x.split("\t")[5]) for x in f
               if x.split("\t")[2] == "intron"}
    assert (1201, "0.9") in introns  # the predictor's intron line is unchanged
    assert (801, ".") in introns


def test_no_intron_line_added_when_transcript_has_none():
    # rebuild_exons on a transcript whose predictor wrote no intron lines
    feats = [
        _line("exon", 1001, 1200, "+", "t1", source="GeneMark.hmm3"),
        _line("exon", 1301, 1500, "+", "t1", source="GeneMark.hmm3"),
        _line("CDS", 1001, 1200, "+", "t1", source="GeneMark.hmm3"),
        _line("CDS", 1301, 1497, "+", "t1", source="GeneMark.hmm3"),
        _line("stop_codon", 1498, 1500, "+", "t1", source="GeneMark.hmm3"),
        _line("five_prime_UTR", 701, 800, "+", "t1", source="StringTie"),
        _line("five_prime_UTR", 901, 1000, "+", "t1", source="StringTie"),
    ]
    f = rebuild_exons({"t1": feats})["t1"]
    assert _coords(f, "exon") == [(701, 800), (901, 1200), (1301, 1500)]
    assert _coords(f, "intron") == []
    assert all(x.split("\t")[1] == "GeneMark.hmm3" for x in f if x.split("\t")[2] == "exon")


def test_transcript_without_utr_is_unchanged():
    f = _run(PLUS_SINGLE)
    assert sorted(f) == sorted(PLUS_SINGLE)


def test_minus_strand_exons_cover_utrs():
    braker = [
        _line("exon", 1001, 1200, "-", "t1"),
        _line("exon", 1301, 1500, "-", "t1"),
        _line("CDS", 1004, 1200, "-", "t1"),
        _line("CDS", 1301, 1500, "-", "t1"),
        _line("intron", 1201, 1300, "-", "t1"),
        _line("stop_codon", 1001, 1003, "-", "t1"),
        _line("start_codon", 1498, 1500, "-", "t1"),
    ]
    f = _run(braker, [_st_exon(801, 1200, "-"), _st_exon(1301, 1700, "-")])
    assert _coords(f, "three_prime_UTR") == [(801, 1000)]
    assert _coords(f, "five_prime_UTR") == [(1501, 1700)]
    assert _coords(f, "exon") == [(801, 1200), (1301, 1700)]
    _utrs_inside_exons(f)
    _no_utr_on_coding(f)


# --- 3. stop codon without a stop_codon line ---------------------------------

def test_genemark_exon_without_stop_codon_line_plus():
    braker = [
        _line("exon", 1001, 2200, "+", "t1", source="GeneMark.hmm3"),
        _line("CDS", 1001, 2197, "+", "t1", source="GeneMark.hmm3"),
    ]
    assert sorted(coding_intervals(braker)) == [(1001, 2197), (2198, 2200)]
    f = _run(braker)
    assert _coords(f, "three_prime_UTR") == []
    f = _run(braker, [_st_exon(1001, 2300, "+")])
    assert _coords(f, "three_prime_UTR") == [(2201, 2300)]


def test_genemark_exon_without_stop_codon_line_minus():
    braker = [
        _line("exon", 1001, 2200, "-", "t1", source="GeneMark.hmm3"),
        _line("CDS", 1004, 2200, "-", "t1", source="GeneMark.hmm3"),
    ]
    assert sorted(coding_intervals(braker)) == [(1001, 1003), (1004, 2200)]
    f = _run(braker, [_st_exon(901, 2200, "-")])
    assert _coords(f, "three_prime_UTR") == [(901, 1000)]


def test_incomplete_transcript_exon_equal_to_cds_has_no_implicit_stop():
    braker = [
        _line("exon", 1001, 2200, "+", "t1"),
        _line("CDS", 1001, 2200, "+", "t1"),
    ]
    assert coding_intervals(braker) == [(1001, 2200)]


# --- 4. intron chain matching with shared introns ---------------------------

def test_shared_introns_are_kept_for_every_transcript():
    def exons_to_tx(tx, exons, strand="+"):
        feats = []
        for i, (s, e) in enumerate(exons):
            feats.append(_line("exon", s, e, strand, tx))
            feats.append(_line("CDS", s, e, strand, tx))
            if i:
                feats.append(_line("intron", exons[i - 1][1] + 1, s - 1, strand, tx))
        return feats

    braker = {
        "t1": exons_to_tx("t1", [(100, 200), (300, 400), (500, 600), (700, 800)]),
        # t2 skips exon 2: shares introns 201-299? no: 201-499 is its own,
        # 601-699 is shared with t1
        "t2": exons_to_tx("t2", [(100, 200), (500, 600), (700, 800)]),
    }
    bi = create_introns_hash(braker)
    assert bi["t1"] == {"Chr1_201_299_+", "Chr1_401_499_+", "Chr1_601_699_+"}
    assert bi["t2"] == {"Chr1_201_499_+", "Chr1_601_699_+"}

    stringtie = {
        # has t1's chain
        "S1": exons_to_tx("S1", [(50, 200), (300, 400), (500, 600), (700, 900)]),
        # has only the last intron of t1: must not match t1 (old code matched
        # when t1's shared introns were credited to t2)
        "S2": exons_to_tx("S2", [(550, 600), (700, 900)]),
        # has t2's chain
        "S3": exons_to_tx("S3", [(100, 200), (500, 600), (700, 850)]),
    }
    si = create_introns_hash(stringtie)
    m = find_matching_transcripts(bi, si)
    assert m["t1"] == ["S1"]
    assert m["t2"] == ["S3"]


def test_stringtie_alternative_transcripts_sharing_an_intron_both_match():
    braker = {"t1": [
        _line("exon", 100, 200, "+", "t1"), _line("exon", 300, 400, "+", "t1"),
        _line("CDS", 100, 200, "+", "t1"), _line("CDS", 300, 400, "+", "t1"),
        _line("intron", 201, 299, "+", "t1"),
    ]}
    stringtie = {
        "S1": [_st_exon(50, 200, "+", "S1"), _st_exon(300, 450, "+", "S1"),
               "\t".join(["Chr1", "StringTie", "intron", "201", "299", ".", "+", ".", 'transcript_id "S1";'])],
        "S2": [_st_exon(90, 200, "+", "S2"), _st_exon(300, 420, "+", "S2"),
               "\t".join(["Chr1", "StringTie", "intron", "201", "299", ".", "+", ".", 'transcript_id "S2";'])],
    }
    m = find_matching_transcripts(create_introns_hash(braker), create_introns_hash(stringtie))
    assert sorted(m["t1"]) == ["S1", "S2"]


# --- print_gtf: exon lines are not rewritten as UTR lines --------------------

def test_print_gtf_rewrites_only_utr_lines(tmp_path):
    f = _run(PLUS_SINGLE, [_st_exon(901, 2300, "+")])
    gtf = {"t1": f}
    out = tmp_path / "out.gtf"
    print_gtf(str(out), gtf,
              {"t1_g": "\t".join(["Chr1", "AUGUSTUS", "gene", "901", "2300", ".", "+", ".", "t1_g"])},
              {"t1": "t1_g"},
              {"t1": "\t".join(["Chr1", "AUGUSTUS", "transcript", "901", "2300", ".", "+", ".", "t1"])})
    lines = [l.rstrip("\n").split("\t") for l in open(out)]
    sources = {l[2]: l[1] for l in lines}
    assert sources["exon"] == "AUGUSTUS"
    assert sources["five_prime_UTR"] == "stringtie2utr"
    assert sources["three_prime_UTR"] == "stringtie2utr"
    assert [l[2] for l in lines][:2] == ["gene", "transcript"]
