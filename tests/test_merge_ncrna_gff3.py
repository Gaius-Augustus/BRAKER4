"""
Tests for scripts/merge_ncrna_gff3.py and the Rfam type mapping in
scripts/infernal_to_gff3.py.

Every ncRNA in the merged GFF3 must be gene -> RNA -> exon with gene_biotype
on the gene, as in NCBI/Ensembl GFF3 (Annotrieve counts genes by biotype and
ignores RNAs without a gene). Overlapping ncRNAs are resolved by priority.

All test data is synthetic and generated in-memory.
"""

import os
import subprocess
import sys

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, SCRIPTS)

from infernal_to_gff3 import load_family_types, parse_tblout, rna_class  # noqa: E402

MERGE = os.path.join(SCRIPTS, "merge_ncrna_gff3.py")

CODING = (
    "##gff-version 3\n"
    "X1\tAUGUSTUS\tgene\t1000\t2000\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n"
    "X1\tAUGUSTUS\tmRNA\t1000\t2000\t.\t+\t.\tID=g1.t1;Parent=g1\n"
    "X1\tAUGUSTUS\texon\t1000\t2000\t.\t+\t.\tParent=g1.t1\n"
    "X1\tAUGUSTUS\tCDS\t1000\t2000\t.\t+\t0\tParent=g1.t1\n"
)
TRNA = (
    "##gff-version 3\n"
    "X1\ttRNAscan-SE\ttRNA\t5000\t5084\t64.2\t-\t.\t"
    "ID=s-X1.trna1;Name=X1.tRNA1-MetCAT;gene_biotype=tRNA;\n"
    "X1\ttRNAscan-SE\texon\t5047\t5084\t.\t-\t.\tID=s-X1.trna1.exon1;Parent=s-X1.trna1;\n"
    "X1\ttRNAscan-SE\texon\t5000\t5035\t.\t-\t.\tID=s-X1.trna1.exon2;Parent=s-X1.trna1;\n"
)
INFERNAL = (
    "##gff-version 3\n"
    # same locus as the tRNAscan-SE tRNA: must be dropped
    "X1\tInfernal\ttRNA\t5001\t5084\t48.7\t-\t.\tID=s-ncRNA_1;Name=tRNA;gene_biotype=tRNA\n"
    # new snoRNA: must become gene -> snoRNA -> exon
    "X1\tInfernal\tsnoRNA\t8000\t8100\t71.1\t+\t.\tID=s-ncRNA_2;Name=snoZ105;gene_biotype=snoRNA\n"
    # inside the coding exon on the same strand: must be dropped
    "X1\tInfernal\tsnRNA\t1200\t1300\t50.0\t+\t.\tID=s-ncRNA_3;Name=U2;gene_biotype=snRNA\n"
    # inside the coding exon on the other strand: kept
    "X1\tInfernal\tmiRNA\t1400\t1480\t40.0\t-\t.\tID=s-ncRNA_4;Name=MIR171;gene_biotype=miRNA\n"
    # cis-regulatory element: no gene, passed through unchanged
    "X1\tInfernal\triboswitch\t1500\t1560\t30.0\t+\t.\tID=s-ncRNA_5;Name=TPP\n"
)
LNCRNA = (
    "##gff-version 3\n"
    "X2\tStringTie\tlnc_RNA\t100\t900\t.\t+\t.\tID=s-lncRNA_1;Name=s-lncRNA_1;biotype=lncRNA\n"
    "X2\tStringTie\texon\t100\t300\t.\t+\t.\tID=s-lncRNA_1.exon1;Parent=s-lncRNA_1\n"
    "X2\tStringTie\texon\t500\t900\t.\t+\t.\tID=s-lncRNA_1.exon2;Parent=s-lncRNA_1\n"
    # orphan exon (Parent not in file): dropped with a warning
    "X2\tStringTie\texon\t2000\t2100\t.\t+\t.\tParent=missing\n"
)


def _attrs(col9):
    return dict(p.split("=", 1) for p in col9.strip(";").split(";") if "=" in p)


def _merge(tmp_path, files):
    paths = []
    for name, text in [("coding.gff3", CODING)] + files:
        path = tmp_path / name
        path.write_text(text)
        paths.append(str(path))
    out = tmp_path / "merged.gff3"
    res = subprocess.run([sys.executable, MERGE, "--coding", paths[0],
                          "--ncrna", *paths[1:], str(tmp_path / "absent.gff3"),
                          "-o", str(out)],
                         capture_output=True, text=True, check=True)
    rows = [line.split("\t") for line in out.read_text().splitlines()
            if line and not line.startswith("#")]
    return rows, res.stderr


def _merged(tmp_path):
    return _merge(tmp_path, [("tRNAs.gff3", TRNA), ("infernal.gff3", INFERNAL),
                             ("lncRNAs.gff3", LNCRNA)])


def test_coding_genes_are_copied_unchanged(tmp_path):
    rows, _ = _merged(tmp_path)
    assert ["\t".join(r) for r in rows[:4]] == CODING.splitlines()[1:]


def test_every_ncrna_has_a_gene_with_biotype(tmp_path):
    rows, _ = _merged(tmp_path)
    ids = {_attrs(r[8]).get("ID"): r for r in rows}
    rna_rows = [r for r in rows if r[2] in {"tRNA", "snoRNA", "miRNA", "lnc_RNA"}]
    assert {r[2] for r in rna_rows} == {"tRNA", "snoRNA", "miRNA", "lnc_RNA"}
    for rna in rna_rows:
        gene = ids[_attrs(rna[8])["Parent"]]
        assert gene[2] == "gene"
        assert (gene[3], gene[4], gene[6]) == (rna[3], rna[4], rna[6])
    biotypes = {_attrs(r[8])["gene_biotype"] for r in rows
                if r[2] == "gene" and r[1] != "AUGUSTUS"}
    assert biotypes == {"tRNA", "snoRNA", "miRNA", "lncRNA"}


def test_exons_point_to_their_rna(tmp_path):
    rows, _ = _merged(tmp_path)
    ids = {_attrs(r[8]).get("ID") for r in rows}
    for r in rows:
        parent = _attrs(r[8]).get("Parent")
        if parent:
            assert parent in ids, r
    trna_exons = [r for r in rows if r[1] == "tRNAscan-SE" and r[2] == "exon"]
    assert len(trna_exons) == 2          # intron-containing tRNA keeps both exons
    sno_exons = [r for r in rows if r[2] == "exon"
                 and _attrs(r[8])["Parent"] == "s-ncRNA_2"]
    assert [(e[3], e[4]) for e in sno_exons] == [("8000", "8100")]  # created


def test_overlaps_are_resolved(tmp_path):
    rows, log = _merged(tmp_path)
    names = {_attrs(r[8]).get("ID") for r in rows}
    assert "s-ncRNA_1" not in names      # Infernal tRNA at the tRNAscan-SE locus
    assert "s-ncRNA_3" not in names      # same strand as the coding exon
    assert "s-ncRNA_4" in names          # opposite strand
    assert "1 dropped (overlap with coding exons), 1 dropped (overlap with " \
           "higher-priority ncRNA)" in log


def test_non_rna_features_pass_through_without_gene(tmp_path):
    rows, _ = _merged(tmp_path)
    riboswitch = [r for r in rows if r[2] == "riboswitch"]
    assert len(riboswitch) == 1
    assert "Parent" not in _attrs(riboswitch[0][8])


def test_orphan_features_are_dropped_with_warning(tmp_path):
    rows, log = _merged(tmp_path)
    assert not [r for r in rows if r[3] == "2000"]
    assert "Parent=missing not in file" in log


def test_rfam_types_map_to_ncbi_feature_types():
    assert rna_class("Gene; tRNA;") == ("tRNA", "tRNA")
    assert rna_class("Gene; rRNA;") == ("rRNA", "rRNA")
    assert rna_class("Gene; snRNA; snoRNA; CD-box;") == ("snoRNA", "snoRNA")
    assert rna_class("Gene; snRNA; snoRNA; scaRNA;") == ("scaRNA", "scaRNA")
    assert rna_class("Gene; snRNA; splicing;") == ("snRNA", "snRNA")
    assert rna_class("Gene; miRNA;") == ("miRNA", "miRNA")
    assert rna_class("Gene; lncRNA;") == ("lnc_RNA", "lncRNA")
    assert rna_class("Gene; sRNA;") == ("ncRNA", "ncRNA")
    assert rna_class("Cis-reg; riboswitch;") == ("riboswitch", None)
    assert rna_class("Cis-reg; leader;") == ("regulatory_region", None)
    assert rna_class("Intron;") == ("autocatalytically_spliced_intron", None)
    assert rna_class("") == ("ncRNA", "ncRNA")


def _tblout_row(idx, model, acc, olp, seq_from, seq_to, strand, score, evalue, inc="!"):
    """One cmscan --fmt 2 tblout line, the 29 fields of Infernal 1.1.5."""
    return " ".join(str(f) for f in (
        idx, model, acc, "X1", "-", "CL00112", "cm", 1, 2900,
        seq_from, seq_to, strand, "no", 1, "0.52", "0.0", score, evalue, inc,
        olp, "-", "-", "-", "-", "-", "-", 2900, 1000000,
        model.replace("_", " ")))


def test_overlap_losers_are_dropped_but_winners_kept(tmp_path):
    """Of an overlapping set only '=' is dropped, after the LSU rRNA example of
    the Infernal guide: the best hit is marked '^', not '*', and keeping only
    '*' lost the locus entirely."""
    tblout = tmp_path / "infernal.tblout"
    tblout.write_text(
        "#idx target name  accession  query name ...\n"
        # the winner of the overlapping set: '^', must be kept
        + _tblout_row(1, "LSU_rRNA_archaea", "RF02540", "^", 762872, 765862, "+", "2763.9", "0") + "\n"
        # overlaps hit 1 and scores worse: '=', must be dropped
        + _tblout_row(2, "LSU_rRNA_bacteria", "RF02541", "=", 762874, 765862, "+", "1872.5", "0") + "\n"
        + _tblout_row(3, "LSU_rRNA_eukarya", "RF02543", "=", 762947, 765766, "+", "1361.0", "0") + "\n"
        # overlaps a worse-scoring loser only: '$', must be kept
        + _tblout_row(4, "LSU_rRNA_eukarya", "RF02543", "$", 765700, 766400, "+", "140.2", "1e-30") + "\n"
        # no overlap at all, on the minus strand with seq from > seq to
        + _tblout_row(5, "tRNA", "RF00005", "*", 5084, 5001, "-", "48.7", "1e-10") + "\n"
        # below the inclusion threshold: dropped whatever its olp says
        + _tblout_row(6, "U2", "RF00004", "^", 9000, 9100, "+", "20.0", "3.0", inc="?") + "\n"
    )
    hits = parse_tblout(str(tblout))
    assert [h['rfam_acc'] for h in hits] == ["RF02540", "RF02543", "RF00005"]
    assert [(h['start'], h['end'], h['strand']) for h in hits] == [
        (762872, 765862, '+'), (765700, 766400, '+'), (5001, 5084, '-')]

def test_shipped_rfam_table_covers_all_families():
    types = load_family_types(os.path.join(SCRIPTS, "rfam_family_types.tsv"))
    assert len(types) > 4000
    assert types["RF00005"] == "Gene; tRNA;"
