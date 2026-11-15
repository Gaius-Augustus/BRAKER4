"""
Tests for start and stop codons split by an intron in scripts/normalize_gtf.py
(legacy BRAKER issue #751).

GeneMark (gmst, GeneMark.hmm3) puts the stop codon into the CDS and, for a
codon split by an intron, writes a start_codon/stop_codon line for one part
only. normalize_gtf.py must read the codon across the intron, trim a split
stop codon from the CDS (dropping a CDS that is only stop codon, keeping its
exon) and write one codon line per genomic part with the GTF frame.

All sequences are synthetic. The plus-strand gene is also placed on the minus
strand of a second sequence (reverse complement) to test both strands.
"""

import os
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "normalize_gtf.py")

INTRON = "GT" + "C" * 16 + "AG"          # 20 bp
PAD = "C" * 10


def _revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _lines(chrom, strand, tid, gid, exons, cds, codons):
    """GeneMark-style lines: transcript (bare ID), exons, CDS, codon lines."""
    attr = f'transcript_id "{tid}"; gene_id "{gid}";'
    start = min(s for s, _ in exons)
    end = max(e for _, e in exons)
    out = [f"{chrom}\tAUGUSTUS\tgene\t{start}\t{end}\t.\t{strand}\t.\t{gid}",
           f"{chrom}\tGeneMark.hmm3\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{tid}"]
    for s, e in exons:
        out.append(f"{chrom}\tGeneMark.hmm3\texon\t{s}\t{e}\t.\t{strand}\t.\t{attr}")
    for s, e, fr in cds:
        out.append(f"{chrom}\tGeneMark.hmm3\tCDS\t{s}\t{e}\t.\t{strand}\t{fr}\t{attr}")
    for feat, s, e, fr in codons:
        out.append(f"{chrom}\tGeneMark.hmm3\t{feat}\t{s}\t{e}\t.\t{strand}\t{fr}\t{attr}")
    return "\n".join(out) + "\n"


def _run(tmp_path, genome, gtf):
    g, inp = tmp_path / "genome.fa", tmp_path / "in.gtf"
    out, log = tmp_path / "out.gtf", tmp_path / "norm.log"
    g.write_text("".join(f">{k}\n{v}\n" for k, v in genome.items()))
    inp.write_text(gtf)
    subprocess.run([sys.executable, SCRIPT, "-g", str(g), "-f", str(inp),
                    "-o", str(out), "-l", str(log)],
                   capture_output=True, text=True, check=True)
    rows = [line.split("\t") for line in out.read_text().splitlines() if line]
    return rows, log.read_text()


def _feats(rows, feature):
    return sorted((int(r[3]), int(r[4]), r[7]) for r in rows if r[2] == feature)


# Plus strand, split stop: ATG GCC GCC TA | intron | G  (TAG split 2 + 1)
STOP_SEQ = PAD + "ATGGCCGCCTA" + INTRON + "G" + PAD    # CDS 11-21, intron 22-41, 42
STOP_L = len(STOP_SEQ)


def _minus(p):
    return STOP_L - p + 1


def test_split_stop_plus_strand(tmp_path):
    gtf = _lines("X1", "+", "g1.t1", "g1", [(11, 21), (42, 42)],
                 [(11, 21, 0), (42, 42, 1)], [("start_codon", 11, 13, 0),
                                              ("stop_codon", 42, 42, 1)])
    rows, log = _run(tmp_path, {"X1": STOP_SEQ}, gtf)
    assert _feats(rows, "CDS") == [(11, 19, "0")]
    assert _feats(rows, "stop_codon") == [(20, 21, "0"), (42, 42, "1")]
    assert [(s, e) for s, e, _ in _feats(rows, "exon")] == [(11, 21), (42, 42)]
    assert "split by an intron" in log


def test_split_stop_minus_strand(tmp_path):
    m = _minus
    gtf = _lines("X2", "-", "g2.t1", "g2", [(m(21), m(11)), (m(42), m(42))],
                 [(m(21), m(11), 0), (m(42), m(42), 1)],
                 [("start_codon", m(13), m(11), 0), ("stop_codon", m(42), m(42), 1)])
    rows, log = _run(tmp_path, {"X2": _revcomp(STOP_SEQ)}, gtf)
    assert _feats(rows, "CDS") == [(m(19), m(11), "0")]
    assert _feats(rows, "stop_codon") == [(m(42), m(42), "1"), (m(21), m(20), "0")]
    assert "split by an intron" in log


# Plus strand, split start: A | intron | TG GCC GCC TAA  (ATG split 1 + 2)
START_SEQ = PAD + "A" + INTRON + "TGGCCGCCTAA" + PAD     # 11, intron 12-31, CDS 32-42


def test_split_start_plus_strand(tmp_path):
    gtf = _lines("X1", "+", "g3.t1", "g3", [(11, 11), (32, 42)],
                 [(11, 11, 0), (32, 42, 2)], [("start_codon", 11, 11, 0),
                                              ("stop_codon", 40, 42, 0)])
    rows, log = _run(tmp_path, {"X1": START_SEQ}, gtf)
    assert _feats(rows, "start_codon") == [(11, 11, "0"), (32, 33, "2")]
    assert _feats(rows, "CDS") == [(11, 11, "0"), (32, 39, "2")]
    assert _feats(rows, "stop_codon") == [(40, 42, "0")]
    assert "non-ATG" not in log


def test_unsplit_codon_lines_unchanged(tmp_path):
    # ATG GCC TAA in one exon: stop trimmed, existing codon lines kept as written
    seq = PAD + "ATGGCCTAA" + PAD
    gtf = _lines("X1", "+", "g4.t1", "g4", [(11, 19)], [(11, 19, 0)],
                 [("start_codon", 11, 13, 0), ("stop_codon", 17, 19, 0)])
    rows, log = _run(tmp_path, {"X1": seq}, gtf)
    assert _feats(rows, "CDS") == [(11, 16, "0")]
    assert _feats(rows, "stop_codon") == [(17, 19, "0")]
    assert _feats(rows, "start_codon") == [(11, 13, "0")]
    assert "split by an intron" not in log
