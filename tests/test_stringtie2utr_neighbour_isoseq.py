"""
Tests for the issue #29 follow-up in stringtie2utr.py.

1. UTRs stop at the nearest same-strand neighbour gene that has a StringTie
   match of its own (short-read StringTie often joins adjacent genes into one
   read-through transcript). Ab initio neighbours do not stop a UTR.
2. With --long-read-stringtie (dual mode), a BRAKER transcript matched in the
   IsoSeq-only assembly takes its UTRs from there; the short-read assembly is
   used for all other transcripts.
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from stringtie2utr import merge_features, neighbour_bounds  # noqa: E402

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "stringtie2utr.py")


def _line(seq, source, feature, start, end, strand, attrs, score=".", frame="."):
    return "\t".join([seq, source, feature, str(start), str(end), score, strand, frame, attrs])


def _braker_tx(gene, tx, exons, strand="+", seq="chr1"):
    """BRAKER transcript features: exon + CDS per (start, end), introns between."""
    attrs = f'transcript_id "{tx}"; gene_id "{gene}";'
    feats = []
    for i, (s, e) in enumerate(exons):
        if i:
            feats.append(_line(seq, "AUGUSTUS", "intron", exons[i - 1][1] + 1, s - 1,
                               strand, attrs))
        feats.append(_line(seq, "AUGUSTUS", "exon", s, e, strand, attrs))
        feats.append(_line(seq, "AUGUSTUS", "CDS", s, e, strand, attrs, frame="0"))
    return feats


def _braker_gtf(genes):
    """genes: list of (gene, tx, exons, strand). Returns GTF text."""
    out = []
    for gene, tx, exons, strand in genes:
        start = min(s for s, _ in exons)
        end = max(e for _, e in exons)
        out.append(_line("chr1", "AUGUSTUS", "gene", start, end, strand, gene))
        out.append(_line("chr1", "AUGUSTUS", "transcript", start, end, strand, tx))
        out.extend(_braker_tx(gene, tx, exons, strand))
    return "\n".join(out) + "\n"


def _stringtie_gtf(transcripts):
    """transcripts: list of (tx_id, exons, strand). Returns StringTie GTF text."""
    out = []
    for tx, exons, strand in transcripts:
        gene = tx.rsplit(".", 1)[0]
        start = min(s for s, _ in exons)
        end = max(e for _, e in exons)
        out.append(_line("chr1", "StringTie", "transcript", start, end, strand,
                         f'gene_id "{gene}"; transcript_id "{tx}";', score="1000"))
        for i, (s, e) in enumerate(exons, 1):
            out.append(_line("chr1", "StringTie", "exon", s, e, strand,
                             f'gene_id "{gene}"; transcript_id "{tx}"; exon_number "{i}";',
                             score="1000"))
    return "\n".join(out) + "\n"


def _utrs(gtf_path, tx):
    res = []
    with open(gtf_path) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if "UTR" in x[2] and f'transcript_id "{tx}"' in x[8]:
                res.append((x[2], int(x[3]), int(x[4])))
    return sorted(res, key=lambda r: r[1])


# ---------------------------------------------------------------------------
# neighbour_bounds
# ---------------------------------------------------------------------------

def _gtf_dict():
    return {
        "g1.t1": _braker_tx("g1", "g1.t1", [(1000, 1500), (1600, 2000)]),
        "g2.t1": _braker_tx("g2", "g2.t1", [(3000, 3500), (3600, 4000)]),
        "g3.t1": _braker_tx("g3", "g3.t1", [(5000, 5500)]),
        "g4.t1": _braker_tx("g4", "g4.t1", [(2500, 2800)], strand="-"),
    }


TX2GENE = {"g1.t1": "g1", "g2.t1": "g2", "g3.t1": "g3", "g4.t1": "g4"}


def test_supported_neighbours_bound_each_other():
    bounds = neighbour_bounds(_gtf_dict(), TX2GENE, {"g1.t1", "g2.t1"})
    assert bounds["g1.t1"] == (None, 2999)
    assert bounds["g2.t1"] == (2001, None)


def test_unsupported_neighbour_does_not_bound():
    # g3 has no StringTie match: g2 may extend its 3' UTR into it
    bounds = neighbour_bounds(_gtf_dict(), TX2GENE, {"g2.t1"})
    assert "g2.t1" not in bounds
    # g3 itself is bounded by the supported g2 on its left
    assert bounds["g3.t1"] == (4001, None)


def test_opposite_strand_neighbour_does_not_bound():
    bounds = neighbour_bounds(_gtf_dict(), TX2GENE, {"g4.t1"})
    assert "g1.t1" not in bounds and "g2.t1" not in bounds


def test_isoform_of_same_gene_is_skipped():
    gtf = _gtf_dict()
    # a second, non-overlapping g1 isoform must not bound g1.t1; g2 still does
    gtf["g1.t2"] = _braker_tx("g1", "g1.t2", [(2200, 2400)])
    tx2gene = dict(TX2GENE, **{"g1.t2": "g1"})
    bounds = neighbour_bounds(gtf, tx2gene, {"g1.t2", "g2.t1"})
    assert bounds["g1.t1"] == (None, 2999)


def test_merge_features_stops_readthrough_at_neighbour():
    gtf = _gtf_dict()
    bounds = neighbour_bounds(gtf, TX2GENE, {"g1.t1", "g2.t1"})
    # read-through transcript covering g1 and g2
    st = {"MSTRG.1.1": [
        _line("chr1", "StringTie", "exon", 800, 1500, "+", 'transcript_id "MSTRG.1.1";'),
        _line("chr1", "StringTie", "exon", 1600, 2300, "+", 'transcript_id "MSTRG.1.1";'),
        _line("chr1", "StringTie", "exon", 2900, 3500, "+", 'transcript_id "MSTRG.1.1";'),
        _line("chr1", "StringTie", "exon", 3600, 4200, "+", 'transcript_id "MSTRG.1.1";'),
    ]}
    merged = merge_features(gtf, st, {"g1.t1": "MSTRG.1.1"}, bounds=bounds)
    st_exons = [(int(f.split("\t")[3]), int(f.split("\t")[4]))
                for f in merged["g1.t1"] if f.split("\t")[1] == "StringTie"]
    assert max(e for _, e in st_exons) <= 2999
    assert (800, 1500) in st_exons and (1600, 2300) in st_exons
    assert (2900, 2999) in st_exons


# ---------------------------------------------------------------------------
# End to end: neighbour stop and IsoSeq priority
# ---------------------------------------------------------------------------

def _run(tmp_path, braker, short, long=None):
    g = tmp_path / "braker.gtf"
    s = tmp_path / "short.gtf"
    o = tmp_path / "out.gtf"
    g.write_text(braker)
    s.write_text(short)
    cmd = [sys.executable, SCRIPT, "-g", str(g), "-s", str(s), "-o", str(o)]
    if long is not None:
        lr = tmp_path / "long.gtf"
        lr.write_text(long)
        cmd += ["-l", str(lr)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return o


BRAKER = _braker_gtf([
    ("g1", "g1.t1", [(1000, 1500), (1600, 2000)], "+"),
    ("g2", "g2.t1", [(3000, 3500), (3600, 4000)], "+"),
])

# Short reads: one read-through transcript over g1 and g2.
SHORT = _stringtie_gtf([
    ("MSTRG.1.1", [(800, 1500), (1600, 3500), (3600, 4300)], "+"),
])


def test_readthrough_utr_stops_before_neighbour(tmp_path):
    # the read-through transcript is the only (and longest) match of both genes
    short = _stringtie_gtf([
        ("MSTRG.1.1", [(800, 1500), (1600, 2400), (2900, 3500), (3600, 4300)], "+"),
    ])
    out = _run(tmp_path, BRAKER, short)
    assert _utrs(out, "g1.t1") == [
        ("five_prime_UTR", 800, 999),
        ("three_prime_UTR", 2001, 2400),
        ("three_prime_UTR", 2900, 2999),
    ]
    assert _utrs(out, "g2.t1") == [
        ("five_prime_UTR", 2001, 2400),
        ("five_prime_UTR", 2900, 2999),
        ("three_prime_UTR", 4001, 4300),
    ]


def test_isoseq_utrs_take_priority(tmp_path):
    long = _stringtie_gtf([
        ("MSTRG.1.1", [(950, 1500), (1600, 2100)], "+"),
    ])
    out = _run(tmp_path, BRAKER, SHORT, long)
    # g1 matched in both assemblies: IsoSeq UTRs, although shorter
    assert _utrs(out, "g1.t1") == [
        ("five_prime_UTR", 950, 999),
        ("three_prime_UTR", 2001, 2100),
    ]
    # g2 not matched by IsoSeq: short-read UTRs from the read-through match,
    # stopped at the supported g1 on its left)
    assert _utrs(out, "g2.t1") == [
        ("five_prime_UTR", 2001, 2999),
        ("three_prime_UTR", 4001, 4300),
    ]


def test_without_long_reads_unchanged_for_isolated_gene(tmp_path):
    braker = _braker_gtf([("g2", "g2.t1", [(3000, 3500), (3600, 4000)], "+")])
    out = _run(tmp_path, braker, SHORT)
    assert _utrs(out, "g2.t1") == [
        ("five_prime_UTR", 800, 1500),
        ("five_prime_UTR", 1600, 2999),
        ("three_prime_UTR", 4001, 4300),
    ]
