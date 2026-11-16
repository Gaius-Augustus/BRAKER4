"""
Tests for scripts/restore_genemark_codons.py.

filterGenemark.pl keeps only the CDS line of single-exon genes; the script
must put the start_codon/stop_codon rows back from the GeneMark GTF, and only
those rows. All test data is synthetic.
"""

import os
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "restore_genemark_codons.py")

SOURCE = (
    # single-exon gene, plus strand
    'X1\tGeneMark.hmm3\tstart_codon\t100\t102\t.\t+\t0\tgene_id "1_g"; transcript_id "1_t";\n'
    'X1\tGeneMark.hmm3\tCDS\t100\t399\t.\t+\t0\tgene_id "1_g"; transcript_id "1_t"; cds_type "Single";\n'
    'X1\tGeneMark.hmm3\tstop_codon\t397\t399\t.\t+\t0\tgene_id "1_g"; transcript_id "1_t";\n'
    # single-exon gene, minus strand
    'X1\tGeneMark.hmm3\tstop_codon\t1000\t1002\t.\t-\t0\tgene_id "2_g"; transcript_id "2_t";\n'
    'X1\tGeneMark.hmm3\tCDS\t1000\t1299\t.\t-\t0\tgene_id "2_g"; transcript_id "2_t"; cds_type "Single";\n'
    'X1\tGeneMark.hmm3\tstart_codon\t1297\t1299\t.\t-\t0\tgene_id "2_g"; transcript_id "2_t";\n'
    # multi-exon gene, already complete in the filtered file
    'X1\tGeneMark.hmm3\tstart_codon\t2000\t2002\t.\t+\t0\tgene_id "3_g"; transcript_id "3_t";\n'
    'X1\tGeneMark.hmm3\tCDS\t2000\t2100\t.\t+\t0\tgene_id "3_g"; transcript_id "3_t";\n'
    'X1\tGeneMark.hmm3\tCDS\t2200\t2301\t.\t+\t1\tgene_id "3_g"; transcript_id "3_t";\n'
    'X1\tGeneMark.hmm3\tstop_codon\t2299\t2301\t.\t+\t0\tgene_id "3_g"; transcript_id "3_t";\n'
    # gene whose ID is reused in the filtered file on another sequence
    'X2\tGeneMark.hmm3\tstart_codon\t100\t102\t.\t+\t0\tgene_id "4_g"; transcript_id "4_t";\n'
)

FILTERED = (
    'X1\tGeneMark.hmm3\tCDS\t100\t399\t.\t+\t0\tgene_id "1_g"; transcript_id "1_t"; cds_type "Single";\n'
    'X1\tGeneMark.hmm3\tCDS\t1000\t1299\t.\t-\t0\tgene_id "2_g"; transcript_id "2_t"; cds_type "Single";\n'
    'X1\tGeneMark.hmm3\tstart_codon\t2000\t2002\t.\t+\t0\tgene_id "3_g"; transcript_id "3_t";\n'
    'X1\tGeneMark.hmm3\tCDS\t2000\t2100\t.\t+\t0\tgene_id "3_g"; transcript_id "3_t";\n'
    'X1\tGeneMark.hmm3\tCDS\t2200\t2301\t.\t+\t1\tgene_id "3_g"; transcript_id "3_t";\n'
    'X1\tGeneMark.hmm3\tstop_codon\t2299\t2301\t.\t+\t0\tgene_id "3_g"; transcript_id "3_t";\n'
    'X1\tGeneMark.hmm3\tCDS\t5000\t5299\t.\t+\t0\tgene_id "4_g"; transcript_id "4_t";\n'
)


def run(tmp_path, filtered):
    src = tmp_path / "genemark.sorted.gtf"
    src.write_text(SOURCE)
    good = tmp_path / "genemark.f.good.gtf"
    good.write_text(filtered)
    empty = tmp_path / "genemark.f.bad.gtf"
    empty.write_text("")
    res = subprocess.run([sys.executable, SCRIPT, "--source", str(src), str(good), str(empty)],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    return good.read_text(), res.stderr


def test_single_exon_codons_restored(tmp_path):
    out, err = run(tmp_path, FILTERED)
    lines = out.splitlines()
    # original lines unchanged and in order, only codon rows added
    assert [l for l in lines if '"4_t"' in l or "\tCDS\t" in l] == \
        [l for l in FILTERED.splitlines() if '"4_t"' in l or "\tCDS\t" in l]
    added = [l for l in lines if l + "\n" not in FILTERED]
    assert len(added) == 4
    t1 = [l.split("\t")[2:5] for l in lines if '"1_t"' in l]
    assert t1 == [["CDS", "100", "399"], ["start_codon", "100", "102"], ["stop_codon", "397", "399"]]
    t2 = {l.split("\t")[2] for l in lines if '"2_t"' in l}
    assert t2 == {"CDS", "start_codon", "stop_codon"}
    # multi-exon gene keeps exactly its two codon rows
    assert sum(1 for l in lines if '"3_t"' in l and "_codon" in l) == 2
    # codon of a different sequence is not copied
    assert not [l for l in lines if '"4_t"' in l and "_codon" in l]
    assert "4 start/stop codon rows restored" in err


def test_nothing_to_do_leaves_file_unchanged(tmp_path):
    complete = SOURCE.split("X2")[0]
    out, err = run(tmp_path, complete)
    assert out == complete
    assert "0 start/stop codon rows restored" in err
