"""
Tests for the ID safety net in scripts/normalize_gtf.py.

- Transcripts with an identical CDS chain are collapsed to one (#71), both
  within a gene and across genes.
- A gene ID used on several loci is split into one gene per locus (#97).

Input mimics TSEBRA output: bare IDs in column 9 of gene and transcript
lines. The genome is poly-C, so no stop codon is trimmed and CDS
coordinates stay as written. All test data is synthetic.
"""

import os
import re
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "normalize_gtf.py")


def _gene(chrom, start, end, strand, gid):
    return f"{chrom}\tAUGUSTUS\tgene\t{start}\t{end}\t.\t{strand}\t.\t{gid}\n"


def _tx(chrom, strand, gid, tid, cds, frames=None):
    """Transcript line plus one exon and one CDS line per (start, end)."""
    frames = frames or ["0"] * len(cds)
    start, end = cds[0][0], cds[-1][1]
    attr = f'transcript_id "{tid}"; gene_id "{gid}";'
    lines = [f"{chrom}\tAUGUSTUS\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{tid}\n"]
    for (s, e), fr in zip(cds, frames):
        lines.append(f"{chrom}\tAUGUSTUS\texon\t{s}\t{e}\t.\t{strand}\t.\t{attr}\n")
        lines.append(f"{chrom}\tAUGUSTUS\tCDS\t{s}\t{e}\t.\t{strand}\t{fr}\t{attr}\n")
    return "".join(lines)


def _run(tmp_path, gtf):
    genome = tmp_path / "genome.fa"
    genome.write_text(">X1\n" + "C" * 5000 + "\n>X2\n" + "C" * 5000 + "\n")
    inp, out, log = tmp_path / "in.gtf", tmp_path / "out.gtf", tmp_path / "norm.log"
    inp.write_text(gtf)
    subprocess.run([sys.executable, SCRIPT, "-g", str(genome), "-f", str(inp),
                    "-o", str(out), "-l", str(log)],
                   capture_output=True, text=True, check=True)
    rows = [line.split("\t") for line in out.read_text().splitlines() if line]
    return rows, log.read_text()


def _ids(col9, key):
    m = re.search(rf'{key} "([^"]+)"', col9)
    return m.group(1) if m else None


def _tx_by_gene(rows):
    """gene_id -> set of transcript_ids, taken from feature lines."""
    out = {}
    for r in rows:
        if r[2] in ("CDS", "exon"):
            out.setdefault(_ids(r[8], "gene_id"), set()).add(_ids(r[8], "transcript_id"))
    return out


CDS_A = [(101, 200), (301, 400), (501, 601)]


def test_identical_isoform_removed(tmp_path):
    gtf = (_gene("X1", 101, 601, "+", "g1")
           + _tx("X1", "+", "g1", "g1.t1", CDS_A)
           + _tx("X1", "+", "g1", "g1.t2", CDS_A)
           + _tx("X1", "+", "g1", "g1.t3", [(101, 200), (301, 601)]))
    rows, log = _run(tmp_path, gtf)
    assert _tx_by_gene(rows) == {"g1": {"g1.t1", "g1.t3"}}
    assert not any("g1.t2" in r[8] for r in rows)
    assert "Duplicate transcripts removed (identical CDS): 1" in log
    assert "DUPLICATE g1.t2: CDS identical to g1.t1" in log


def test_identical_gene_removed(tmp_path):
    gtf = (_gene("X1", 101, 601, "+", "g1") + _tx("X1", "+", "g1", "g1.t1", CDS_A)
           + _gene("X1", 101, 601, "+", "g2") + _tx("X1", "+", "g2", "g2.t1", CDS_A))
    rows, log = _run(tmp_path, gtf)
    assert _tx_by_gene(rows) == {"g1": {"g1.t1"}}
    assert [r[8] for r in rows if r[2] == "gene"] == ["g1"]
    assert "DUPLICATE gene g2" in log


def test_different_strand_or_frame_is_not_a_duplicate(tmp_path):
    gtf = (_gene("X1", 101, 601, "+", "g1") + _tx("X1", "+", "g1", "g1.t1", CDS_A)
           + _gene("X1", 101, 601, "-", "g2") + _tx("X1", "-", "g2", "g2.t1", CDS_A)
           + _gene("X1", 101, 601, "+", "g3")
           + _tx("X1", "+", "g3", "g3.t1", CDS_A, frames=["1", "0", "0"]))
    rows, log = _run(tmp_path, gtf)
    assert set(_tx_by_gene(rows)) == {"g1", "g2", "g3"}
    assert "Duplicate transcripts removed (identical CDS): 0" in log


def test_gene_id_on_two_sequences_is_split(tmp_path):
    # Same gene and transcript IDs on X1 and X2, as reported in #97
    gtf = (_gene("X1", 101, 601, "+", "g1") + _tx("X1", "+", "g1", "g1.t1", CDS_A)
           + _gene("X2", 1101, 1601, "-", "g1")
           + _tx("X2", "-", "g1", "g1.t1", [(1101, 1200), (1301, 1601)]))
    rows, log = _run(tmp_path, gtf)
    genes = {r[8]: (r[0], r[6], r[3], r[4]) for r in rows if r[2] == "gene"}
    assert genes == {"g1": ("X1", "+", "101", "601"),
                     "g1_2": ("X2", "-", "1101", "1601")}
    txs = {r[8]: (r[0], r[6]) for r in rows if r[2] == "transcript"}
    assert txs == {"g1.t1": ("X1", "+"), "g1_2.t1": ("X2", "-")}
    for r in rows:
        if r[2] == "CDS":
            expected = ("g1", "g1.t1") if r[0] == "X1" else ("g1_2", "g1_2.t1")
            assert (_ids(r[8], "gene_id"), _ids(r[8], "transcript_id")) == expected
    assert "Gene IDs split (several loci): 1" in log
    assert "SPLIT g1: 2 loci -> g1, g1_2" in log


def test_non_overlapping_transcripts_are_split(tmp_path):
    gtf = (_gene("X1", 101, 3600, "+", "g1")
           + _tx("X1", "+", "g1", "g1.t1", CDS_A)
           + _tx("X1", "+", "g1", "g1.t2", [(3001, 3300), (3401, 3600)]))
    rows, log = _run(tmp_path, gtf)
    assert _tx_by_gene(rows) == {"g1": {"g1.t1"}, "g1_2": {"g1_2.t2"}}
    genes = {r[8]: (r[3], r[4]) for r in rows if r[2] == "gene"}
    assert genes == {"g1": ("101", "601"), "g1_2": ("3001", "3600")}
    assert "SPLIT g1: 2 loci -> g1, g1_2" in log


def test_overlapping_isoforms_stay_one_gene(tmp_path):
    gtf = (_gene("X1", 101, 900, "+", "g1")
           + _tx("X1", "+", "g1", "g1.t1", CDS_A)
           + _tx("X1", "+", "g1", "g1.t2", [(501, 700), (801, 900)]))
    rows, log = _run(tmp_path, gtf)
    assert _tx_by_gene(rows) == {"g1": {"g1.t1", "g1.t2"}}
    assert "Gene IDs split (several loci): 0" in log
