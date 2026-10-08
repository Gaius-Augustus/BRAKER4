"""
Tests for scripts/feelnc_to_gff3.py: the lncRNA GTF of FEELnc_codpot.pl (exon
lines only) becomes lnc_RNA -> exon GFF3 with IDs <stem>-lncRNA_<n> in genome
order, the StringTie transcript ID as Alias, and the comment lines of the GTF
(the notes of rule run_feelnc) in the header. Synthetic data only.
"""

import os
import subprocess
import sys

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
SCRIPT = os.path.join(SCRIPTS, "feelnc_to_gff3.py")

# FEELnc writes exon lines only, not in genome order
LNC_GTF = (
    'X2\tStringTie\texon\t90543\t90692\t1000\t-\t.\tgene_id "STRG.44"; transcript_id "STRG.44.1"; exon_number "2";\n'
    'X2\tStringTie\texon\t90061\t90452\t1000\t-\t.\tgene_id "STRG.44"; transcript_id "STRG.44.1"; exon_number "1";\n'
    'X1\tStringTie\texon\t5000\t5100\t1000\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.2"; exon_number "1";\n'
    'X1\tStringTie\texon\t5300\t5400\t1000\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.2"; exon_number "2";\n'
    'X1\tStringTie\texon\t1000\t1200\t1000\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.1"; exon_number "1";\n'
    'X1\tStringTie\texon\t1500\t1600\t1000\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.1"; exon_number "2";\n'
)


def run(tmp_path, text):
    gtf = tmp_path / "lncRNAs.gtf"
    gtf.write_text(text)
    out = tmp_path / "lncRNAs.gff3"
    proc = subprocess.run([sys.executable, SCRIPT, "--stem", "s", str(gtf), "-o", str(out)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return out.read_text().splitlines()


def test_transcripts_in_genome_order_with_exons(tmp_path):
    lines = run(tmp_path, LNC_GTF)
    assert lines[0] == "##gff-version 3"
    rows = [l.split("\t") for l in lines[1:]]
    assert [r[2] for r in rows] == ["lnc_RNA", "exon", "exon"] * 3
    rnas = [r for r in rows if r[2] == "lnc_RNA"]
    assert [(r[0], r[3], r[4], r[6]) for r in rnas] == [
        ("X1", "1000", "1600", "+"), ("X1", "5000", "5400", "+"), ("X2", "90061", "90692", "-")]
    assert rnas[0][8] == "ID=s-lncRNA_1;Name=s-lncRNA_1;Alias=STRG.9.1;biotype=lncRNA"
    assert rnas[2][8] == "ID=s-lncRNA_3;Name=s-lncRNA_3;Alias=STRG.44.1;biotype=lncRNA"
    # exons follow their transcript, sorted by start, IDs <id>.exon<k>
    assert rows[1][3:5] == ["1000", "1200"] and rows[1][8] == "ID=s-lncRNA_1.exon1;Parent=s-lncRNA_1"
    assert rows[2][3:5] == ["1500", "1600"] and rows[2][8] == "ID=s-lncRNA_1.exon2;Parent=s-lncRNA_1"
    assert rows[7][3:5] == ["90061", "90452"] and rows[7][8] == "ID=s-lncRNA_3.exon1;Parent=s-lncRNA_3"


def test_comment_lines_of_an_empty_gtf_become_the_header(tmp_path):
    note = "# FEELnc: 3 candidate transcripts, 250 annotated transcripts; FEELnc_codpot.pl needs at least 100 of each to train, no lncRNA called"
    lines = run(tmp_path, note + "\n")
    assert lines == ["##gff-version 3", note]


def test_empty_file(tmp_path):
    assert run(tmp_path, "") == ["##gff-version 3"]


def test_exon_without_transcript_id_fails(tmp_path):
    gtf = tmp_path / "bad.gtf"
    gtf.write_text('X1\tStringTie\texon\t1\t100\t.\t+\t.\tgene_id "STRG.1";\n')
    proc = subprocess.run([sys.executable, SCRIPT, "--stem", "s", str(gtf), "-o", str(tmp_path / "o.gff3")],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert "transcript_id" in proc.stderr
