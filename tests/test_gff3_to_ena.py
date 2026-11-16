"""
Tests for scripts/gff3_to_ena.py (BRAKER4 GFF3 -> GFF3 for ENA decoupled
annotation, https://ena-docs.readthedocs.io/en/latest/submit/gff3.html).

The expectations follow ENA's validator gff3tools (enasequence/gff3tools,
main of 2026-11-16): only INSDC qualifiers, one ID per multi-row feature,
codon_start for the frame, partial in genomic orientation, locus_tag
PREFIX_number, tRNA/rRNA overlap and length rules.

All test data is synthetic: a random genome with planted gene models.
"""

import os
import random
import subprocess
import sys

import pytest

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "gff3_to_ena.py")
sys.path.insert(0, os.path.dirname(SCRIPT))

from gff3_to_ena import Code, merge_intervals, phases, subtract, encode  # noqa: E402

STOPS = {"TAA", "TAG", "TGA"}
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in STOPS and a + b + c != "ATG"]
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def codons(rng, n):
    return "".join(rng.choice(SENSE) for _ in range(n))


class Genome:
    """Random background with stop-free planted sequences."""

    def __init__(self, length, seed=1):
        self.rng = random.Random(seed)
        self.seq = list("".join(self.rng.choice("ACGT") for _ in range(length)))

    def put(self, pos, s):
        self.seq[pos - 1:pos - 1 + len(s)] = list(s)

    def text(self):
        return "".join(self.seq)


def run(tmp_path, gff, fasta, *extra, prefix="BRKTEST", check=True):
    g = tmp_path / "in.gff3"
    g.write_text(gff)
    fa = tmp_path / "genome.fa"
    fa.write_text(fasta)
    out = tmp_path / "out.gff3"
    rep = tmp_path / "report.txt"
    ids = tmp_path / "ids.tsv"
    cmd = [sys.executable, SCRIPT, "--gff3", str(g), "--genome", str(fa),
           "--locus-tag-prefix", prefix, "-o", str(out), "--report", str(rep),
           "--id-map", str(ids), *extra]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if check:
        assert res.returncode == 0, res.stderr
    return res, (out.read_text() if out.exists() else ""), \
        (rep.read_text() if rep.exists() else "")


def rows(text):
    out = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split("\t")
        attrs = dict(kv.split("=", 1) for kv in c[8].split(";"))
        out.append((c[0], c[2], int(c[3]), int(c[4]), c[6], c[7], attrs))
    return out


ALLOWED = {"ID", "Parent", "locus_tag", "product", "transl_table", "codon_start",
           "inference", "partial", "pseudo", "note", "Dbxref"}


def assert_ena_shape(text):
    """Rules every output must meet, independent of the case."""
    lines = text.splitlines()
    assert lines[0] == "##gff-version 3"
    declared = set()
    for line in lines[1:]:
        if line.startswith("##sequence-region"):
            declared.add(line.split()[1])
            continue
        if line.startswith("#"):
            continue
        c = line.split("\t")
        assert len(c) == 9
        assert c[0] in declared
        assert c[6] in "+-"
        assert int(c[3]) <= int(c[4])
        assert (c[7] in "012") if c[2] == "CDS" else c[7] == "."
        assert c[2] not in ("start_codon", "stop_codon", "intron", "transcript")
        for kv in c[8].split(";"):
            key, value = kv.split("=", 1)
            assert key in ALLOWED, key
            assert "=" not in value and ";" not in value


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def test_merge_subtract_phases():
    assert merge_intervals([(5, 9), (1, 3), (4, 4)]) == [(1, 9)]
    assert subtract([(1, 10)], [(4, 6)]) == [(1, 3), (7, 10)]
    assert subtract([(1, 3)], [(1, 3)]) == []
    # 10 nt at phase 0 leaves 1 nt: next segment phase 2
    assert phases([(1, 10), (20, 40)], 0) == [0, 2]
    assert phases([(1, 10), (20, 40)], 1) == [1, 0]


def test_codes():
    std, cil = Code(1), Code(6)
    assert std.is_stop("TAA") and not cil.is_stop("TAA") and cil.is_stop("TGA")
    assert std.is_start("TTG") and std.is_start("ATG") and not cil.is_start("TTG")
    with pytest.raises(Exception):
        Code(2)


def test_encode():
    assert encode("a;b=c,d&e%f+g") == "a%3Bb%3Dc%2Cd%26e%25f%2Bg"


# ---------------------------------------------------------------------------
# protein-coding models
# ---------------------------------------------------------------------------

def coding_case():
    """Plus-strand complete gene in GTF style (stop outside the CDS) and a
    minus-strand gene without start codon whose 5' end has phase 1."""
    gen = Genome(6000)
    rng = gen.rng
    ex1 = "ATG" + codons(rng, 30)                 # 93 nt
    intron = "GT" + "".join(rng.choice("ACGT") for _ in range(96)) + "AG"
    ex2 = codons(rng, 20) + "TAA"                 # 63 nt
    gen.put(1001, ex1)
    gen.put(1094, intron)
    gen.put(1194, ex2)
    # minus-strand gene, biological sequence: 1 nt + 40 codons + TGA
    bio = "C" + codons(rng, 40) + "TGA"           # 124 nt
    gen.put(3001, rc(bio))                        # 3001..3124
    gff = (
        "##gff-version 3\n"
        "X1\tAUGUSTUS\tgene\t1001\t1300\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
        "X1\tAUGUSTUS\tmRNA\t1001\t1300\t0.9\t+\t.\tID=g1.t1;Parent=g1;gene_id=g1\n"
        "X1\tAUGUSTUS\texon\t1001\t1093\t.\t+\t.\tID=e1;Parent=g1.t1\n"
        "X1\tAUGUSTUS\texon\t1194\t1300\t.\t+\t.\tID=e2;Parent=g1.t1\n"
        "X1\tAUGUSTUS\tCDS\t1001\t1093\t.\t+\t0\tID=c1;Parent=g1.t1\n"
        "X1\tAUGUSTUS\tCDS\t1194\t1253\t.\t+\t0\tID=c2;Parent=g1.t1\n"
        "X1\tAUGUSTUS\tintron\t1094\t1193\t.\t+\t.\tID=i1;Parent=g1.t1\n"
        "X1\tAUGUSTUS\tstart_codon\t1001\t1003\t.\t+\t0\tID=sc;Parent=g1.t1\n"
        "X1\tAUGUSTUS\tstop_codon\t1254\t1256\t.\t+\t0\tID=st;Parent=g1.t1\n"
        "X1\tstringtie2utr\tthree_prime_UTR\t1257\t1300\t.\t+\t.\tID=u;Parent=g1.t1\n"
        "X1\tAUGUSTUS\tgene\t3001\t3124\t.\t-\t.\tID=g2\n"
        "X1\tAUGUSTUS\tmRNA\t3001\t3124\t.\t-\t.\tID=g2.t1;Parent=g2\n"
        "X1\tAUGUSTUS\texon\t3001\t3124\t.\t-\t.\tID=e3;Parent=g2.t1\n"
        "X1\tAUGUSTUS\tCDS\t3004\t3124\t.\t-\t1\tID=c3;Parent=g2.t1\n"
        "X1\tAUGUSTUS\tstop_codon\t3001\t3003\t.\t-\t0\tID=st2;Parent=g2.t1\n"
    )
    return gff, ">X1 test\n" + gen.text() + "\n"


def test_complete_and_5prime_partial(tmp_path):
    gff, fa = coding_case()
    versions = tmp_path / "software_versions.tsv"
    versions.write_text("AUGUSTUS\t3.5.0 (commit 220e5a6)\nGeneMark-ET\t4.* (commit a)\n")
    _, out, rep = run(tmp_path, gff, fa, "--taxid", "3702",
                      "--software-versions", str(versions))
    assert_ena_shape(out)
    assert "##species https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/3702" in out
    assert "##sequence-region X1 1 6000" in out
    r = rows(out)
    genes = [x for x in r if x[1] == "gene"]
    assert [g[6]["locus_tag"] for g in genes] == ["BRKTEST_000010", "BRKTEST_000020"]
    # every feature of g1 carries g1's locus tag
    assert all(x[6]["locus_tag"] == "BRKTEST_000010" for x in r if x[2] < 2000)

    cds1 = [x for x in r if x[1] == "CDS" and x[2] < 2000]
    # CDS includes the stop codon, one ID, phases recomputed
    assert [(x[2], x[3], x[5]) for x in cds1] == [(1001, 1093, "0"), (1194, 1256, "0")]
    assert {x[6]["ID"] for x in cds1} == {"cds-g1.t1"}
    a = cds1[0][6]
    assert a["product"] == "hypothetical protein" and a["transl_table"] == "1"
    assert a["inference"] == "ab initio prediction:AUGUSTUS:3.5.0"
    assert "partial" not in a and "codon_start" not in a
    mrna1 = [x for x in r if x[1] == "mRNA" and x[2] < 2000][0]
    assert mrna1[6]["Parent"] == "g1" and "partial" not in mrna1[6]
    utr = [x for x in r if x[1] == "three_prime_UTR"]
    assert [(x[2], x[3]) for x in utr] == [(1257, 1300)]

    # minus strand, 5' end (upper coordinate) partial, frame from phase 1
    cds2 = [x for x in r if x[1] == "CDS" and x[2] > 2000]
    assert [(x[2], x[3], x[5]) for x in cds2] == [(3001, 3124, "1")]
    a = cds2[0][6]
    assert a["codon_start"] == "2" and a["partial"] == "end"
    for t in ("mRNA", "gene", "exon"):
        assert [x for x in r if x[1] == t and x[2] > 2000][0][6]["partial"] == "end"
    assert "partial_5prime\t1" in rep


def test_3prime_partial_drops_utr_and_clips_exon(tmp_path):
    gen = Genome(3000, seed=2)
    cds = "ATG" + codons(gen.rng, 50) + "GA"     # ends without stop, 155 nt
    gen.put(501, cds)                            # 501..655
    gff = (
        "X1\tAUGUSTUS\tgene\t501\t700\t.\t+\t.\tID=g1\n"
        "X1\tAUGUSTUS\tmRNA\t501\t700\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "X1\tAUGUSTUS\texon\t501\t700\t.\t+\t.\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tCDS\t501\t655\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tstart_codon\t501\t503\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAGAT\tthree_prime_UTR\t656\t700\t.\t+\t.\tParent=g1.t1\n"
    )
    _, out, rep = run(tmp_path, gff, ">X1\n" + gen.text() + "\n")
    assert_ena_shape(out)
    r = rows(out)
    assert not [x for x in r if x[1] == "three_prime_UTR"]
    exon = [x for x in r if x[1] == "exon"][0]
    assert (exon[2], exon[3], exon[6]["partial"]) == (501, 655, "end")
    assert [x for x in r if x[1] == "gene"][0][3] == 655
    assert [x for x in r if x[1] == "CDS"][0][6]["partial"] == "end"
    assert "utr_removed_at_partial_end\t1" in rep


def test_missing_codon_rows_stop_taken_from_sequence(tmp_path):
    """GeneMark transcript without start/stop rows; CDS stops before TGA, which
    stringtie2utr labelled as 3' UTR."""
    gen = Genome(3000, seed=3)
    gen.put(1001, "ATG" + codons(gen.rng, 40) + "TGA")   # 1001..1126
    gff = (
        "X1\tGeneMark.hmm3\tgene\t1001\t1126\t.\t+\t.\tID=g1\n"
        "X1\tGeneMark.hmm3\tmRNA\t1001\t1126\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "X1\tGeneMark.hmm3\texon\t1001\t1126\t.\t+\t0\tParent=g1.t1\n"
        "X1\tGeneMark.hmm3\tCDS\t1001\t1123\t.\t+\t0\tParent=g1.t1\n"
        "X1\tstringtie2utr\tthree_prime_UTR\t1124\t1126\t.\t+\t0\tParent=g1.t1\n"
    )
    _, out, rep = run(tmp_path, gff, ">X1\n" + gen.text() + "\n")
    assert_ena_shape(out)
    r = rows(out)
    cds = [x for x in r if x[1] == "CDS"][0]
    assert (cds[2], cds[3]) == (1001, 1126) and "partial" not in cds[6]
    assert not [x for x in r if "UTR" in x[1]]
    assert "stop_codon_taken_from_sequence\t1" in rep


def test_internal_stop_is_pseudo(tmp_path):
    gen = Genome(3000, seed=4)
    gen.put(1001, "ATG" + codons(gen.rng, 10) + "TAG" + codons(gen.rng, 10) + "TAA")
    gff = (
        "X1\tAUGUSTUS\tgene\t1001\t1069\t.\t+\t.\tID=g1\n"
        "X1\tAUGUSTUS\tmRNA\t1001\t1069\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "X1\tAUGUSTUS\tCDS\t1001\t1069\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tstart_codon\t1001\t1003\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tstop_codon\t1067\t1069\t.\t+\t0\tParent=g1.t1\n"
    )
    _, out, rep = run(tmp_path, gff, ">X1\n" + gen.text() + "\n")
    assert_ena_shape(out)
    a = [x for x in rows(out) if x[1] == "CDS"][0][6]
    assert a["pseudo"] == "true"
    assert "product" not in a and "transl_table" not in a
    # exons are rebuilt from the CDS when the input has none
    assert [x for x in rows(out) if x[1] == "exon"]
    assert "pseudo_internal_stop\t1" in rep


def test_ciliate_code(tmp_path):
    """Table 6: TAA is Gln, only TGA stops."""
    gen = Genome(3000, seed=5)
    gen.put(1001, "ATG" + codons(gen.rng, 10) + "TAA" + codons(gen.rng, 10) + "TGA")
    gff = (
        "X1\tAUGUSTUS\tgene\t1001\t1069\t.\t+\t.\tID=g1\n"
        "X1\tAUGUSTUS\tmRNA\t1001\t1069\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "X1\tAUGUSTUS\tCDS\t1001\t1069\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tstart_codon\t1001\t1003\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tstop_codon\t1067\t1069\t.\t+\t0\tParent=g1.t1\n"
    )
    _, out, _ = run(tmp_path, gff, ">X1\n" + gen.text() + "\n", "--translation-table", "6")
    a = [x for x in rows(out) if x[1] == "CDS"][0][6]
    assert a["transl_table"] == "6" and "pseudo" not in a and "partial" not in a


# ---------------------------------------------------------------------------
# ncRNA
# ---------------------------------------------------------------------------

def ncrna_case():
    gen = Genome(8000, seed=6)
    gen.put(1001, "ATG" + codons(gen.rng, 300) + "TAA")   # CDS 1001..1906
    gff = (
        "X1\tAUGUSTUS\tgene\t1001\t1906\t.\t+\t.\tID=g1\n"
        "X1\tAUGUSTUS\tmRNA\t1001\t1906\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "X1\tAUGUSTUS\tCDS\t1001\t1906\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tstart_codon\t1001\t1003\t.\t+\t0\tParent=g1.t1\n"
        "X1\tAUGUSTUS\tstop_codon\t1904\t1906\t.\t+\t0\tParent=g1.t1\n"
        # tRNA with intron, minus strand: joined over its exons
        "X1\ttRNAscan-SE\tgene\t5000\t5084\t.\t-\t.\tID=gene-t1;gene_biotype=tRNA;Name=X1.tRNA1-MetCAT\n"
        "X1\ttRNAscan-SE\ttRNA\t5000\t5084\t64.2\t-\t.\tID=t1;Parent=gene-t1;Name=X1.tRNA1-MetCAT;isotype=Met;anticodon=CAT\n"
        "X1\ttRNAscan-SE\texon\t5047\t5084\t.\t-\t.\tID=t1.exon1;Parent=t1\n"
        "X1\ttRNAscan-SE\texon\t5000\t5035\t.\t-\t.\tID=t1.exon2;Parent=t1\n"
        # tRNA of 40 bp: outside ENA's 50-150 bp
        "X1\ttRNAscan-SE\tgene\t5200\t5239\t.\t+\t.\tID=gene-t2\n"
        "X1\ttRNAscan-SE\ttRNA\t5200\t5239\t30\t+\t.\tID=t2;Parent=gene-t2;isotype=Gly\n"
        # tRNA inside the CDS (other strand): ENA rejects it
        "X1\ttRNAscan-SE\tgene\t1200\t1272\t.\t-\t.\tID=gene-t3\n"
        "X1\ttRNAscan-SE\ttRNA\t1200\t1272\t30\t-\t.\tID=t3;Parent=gene-t3;isotype=Arg\n"
        # pseudo tRNA
        "X1\ttRNAscan-SE\tpseudogene\t5400\t5472\t.\t+\t.\tID=gene-t4\n"
        "X1\ttRNAscan-SE\tpseudogenic_tRNA\t5400\t5472\t20\t+\t.\tID=t4;Parent=gene-t4;isotype=Undet\n"
        # rRNA overlapping the CDS on the same strand: dropped
        "X1\tbarrnap\tgene\t1800\t1950\t.\t+\t.\tID=gene-r1\n"
        "X1\tbarrnap\trRNA\t1800\t1950\t.\t+\t.\tID=r1;Parent=gene-r1;Name=s-rRNA_1_5S_rRNA;product=5S ribosomal RNA\n"
        # partial 18S rRNA on its own
        "X1\tbarrnap\tgene\t6000\t6500\t.\t+\t.\tID=gene-r2\n"
        "X1\tbarrnap\trRNA\t6000\t6500\t.\t+\t.\tID=r2;Parent=gene-r2;product=18S ribosomal RNA (partial);note=aligned only 30 percent of the 18S ribosomal RNA\n"
        # Infernal snoRNA and miRNA hairpin
        "X1\tInfernal\tgene\t7000\t7104\t.\t-\t.\tID=gene-n1;gene_biotype=snoRNA\n"
        "X1\tInfernal\tsnoRNA\t7000\t7104\t104\t-\t.\tID=n1;Parent=gene-n1;Name=snoR126;Dbxref=RFAM:RF01426;evalue=2e-19;note=Infernal cmscan hit to snoR126 (RF01426);rfam_type=Gene,snRNA\n"
        "X1\tInfernal\texon\t7000\t7104\t.\t-\t.\tID=n1.exon1;Parent=n1\n"
        "X1\tInfernal\tgene\t7300\t7385\t.\t-\t.\tID=gene-n2\n"
        "X1\tInfernal\tmiRNA\t7300\t7385\t66\t-\t.\tID=n2;Parent=gene-n2;Name=MIR171_1;Dbxref=RFAM:RF00643\n"
        # cis-regulatory element without gene
        "X1\tInfernal\triboswitch\t7600\t7660\t30\t+\t.\tID=n3;Name=TPP;Dbxref=RFAM:RF00059\n"
        # lncRNA from FEELnc, two exons
        "X1\tFEELnc\tgene\t2500\t3500\t.\t+\t.\tID=gene-l1;gene_biotype=lncRNA\n"
        "X1\tFEELnc\tlnc_RNA\t2500\t3500\t.\t+\t.\tID=l1;Parent=gene-l1\n"
        "X1\tFEELnc\texon\t2500\t2800\t.\t+\t.\tParent=l1\n"
        "X1\tFEELnc\texon\t3200\t3500\t.\t+\t.\tParent=l1\n"
    )
    return gff, ">X1\n" + gen.text() + "\n"


def test_ncrna(tmp_path):
    gff, fa = ncrna_case()
    _, out, rep = run(tmp_path, gff, fa)
    assert_ena_shape(out)
    r = rows(out)
    t1 = [x for x in r if x[6].get("ID") == "t1"]
    assert [(x[2], x[3]) for x in t1] == [(5000, 5035), (5047, 5084)]
    assert t1[0][1] == "tRNA" and t1[0][6]["product"] == "tRNA-Met"
    assert t1[0][6]["Parent"] == "gene-t1"
    assert not [x for x in r if x[1] == "exon" and x[2] >= 5000]
    ids = {x[6].get("ID") for x in r}
    assert "t2" not in ids and "t3" not in ids and "r1" not in ids
    t4 = [x for x in r if x[6].get("ID") == "t4"][0]
    assert t4[1] == "tRNA" and t4[6]["pseudo"] == "true" and "product" not in t4[6]
    assert [x for x in r if x[6].get("ID") == "gene-t4"][0][6]["pseudo"] == "true"
    r2 = [x for x in r if x[6].get("ID") == "r2"][0]
    assert r2[6]["product"] == "18S ribosomal RNA"
    assert "aligned only 30 percent" in r2[6]["note"]
    n1 = [x for x in r if x[6].get("ID") == "n1"][0]
    assert n1[1] == "snoRNA" and n1[6]["Dbxref"] == "RFAM:RF01426"
    assert "rfam_type" not in n1[6] and "evalue" not in n1[6]
    assert [x for x in r if x[6].get("ID") == "n2"][0][1] == "pre_miRNA"
    rs = [x for x in r if x[1] == "riboswitch"][0]
    assert "locus_tag" not in rs[6] and rs[6]["Dbxref"] == "RFAM:RF00059"
    l1 = [x for x in r if x[6].get("ID") == "l1"]
    assert [(x[1], x[2], x[3]) for x in l1] == [("lnc_RNA", 2500, 2800), ("lnc_RNA", 3200, 3500)]
    for key in ("dropped_trna_length\t1", "dropped_trna_in_cds\t1",
                "dropped_rrna_overlap_cds\t1"):
        assert key in rep
    # locus tags only on genes and their parts, unique per gene
    tags = [x[6]["locus_tag"] for x in r if x[1] == "gene"]
    assert len(tags) == len(set(tags))


# ---------------------------------------------------------------------------
# sequence names, prefix, errors
# ---------------------------------------------------------------------------

SIMPLE = (
    "{sid}\tAUGUSTUS\tgene\t1001\t1072\t.\t+\t.\tID=g1\n"
    "{sid}\tAUGUSTUS\tmRNA\t1001\t1072\t.\t+\t.\tID=g1.t1;Parent=g1\n"
    "{sid}\tAUGUSTUS\tCDS\t1001\t1072\t.\t+\t0\tParent=g1.t1\n"
    "{sid}\tAUGUSTUS\tstart_codon\t1001\t1003\t.\t+\t0\tParent=g1.t1\n"
    "{sid}\tAUGUSTUS\tstop_codon\t1070\t1072\t.\t+\t0\tParent=g1.t1\n"
)


def simple_genome(sid):
    gen = Genome(2000, seed=7)
    gen.put(1001, "ATG" + codons(gen.rng, 22) + "TAA")
    return f">{sid}\n" + gen.text() + "\n"


def test_invalid_prefix(tmp_path):
    res, _, _ = run(tmp_path, SIMPLE.format(sid="X1"), simple_genome("X1"),
                    prefix="br", check=False)
    assert res.returncode != 0 and "locus_tag prefix" in res.stderr


def test_seqid_with_dot_needs_map(tmp_path):
    res, _, _ = run(tmp_path, SIMPLE.format(sid="ctg.a"), simple_genome("ctg.a"),
                    check=False)
    assert res.returncode != 0 and "seqid map" in res.stderr
    m = tmp_path / "map.tsv"
    m.write_text("ctg.a\tCAXYZ010000001.1\n")
    _, out, _ = run(tmp_path, SIMPLE.format(sid="ctg.a"), simple_genome("ctg.a"),
                    "--seqid-map", str(m))
    assert "##sequence-region CAXYZ010000001.1 1 2000" in out
    assert all(x[0] == "CAXYZ010000001.1" for x in rows(out))
    assert "CAXYZ010000001.1\tg1\tg1.t1\tBRKTEST_000010" in (tmp_path / "ids.tsv").read_text()


def test_unmapped_sequence_and_missing_sequence(tmp_path):
    m = tmp_path / "map.tsv"
    m.write_text("other\tOTHER.1\n")
    res, _, _ = run(tmp_path, SIMPLE.format(sid="X1"), simple_genome("X1"),
                    "--seqid-map", str(m), check=False)
    assert res.returncode != 0 and "not in the seqid map" in res.stderr
    res, _, _ = run(tmp_path, SIMPLE.format(sid="X9"), simple_genome("X1"), check=False)
    assert res.returncode != 0 and "missing from the genome" in res.stderr


def test_short_complete_cds_has_inference(tmp_path):
    """Complete CDS under 25 aa need /inference (gff3tools CDS_LENGTH)."""
    _, out, rep = run(tmp_path, SIMPLE.format(sid="X1"), simple_genome("X1"))
    a = [x for x in rows(out) if x[1] == "CDS"][0][6]
    assert a["inference"].startswith("ab initio prediction:AUGUSTUS")
    assert "complete_cds_lt_25aa\t1" in rep
