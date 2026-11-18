"""
Tests for the VARUS sentence of generate_methods_text() and the local-FASTQ
section of generate_html() in scripts/generate_report.py. Fixture: a tmp
directory with only varus/varus_stats.txt.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import generate_report as gr  # noqa: E402


def methods(tmp_path, stats):
    (tmp_path / "varus").mkdir()
    if stats is not None:
        (tmp_path / "varus" / "varus_stats.txt").write_text(stats)
    return gr.generate_methods_text(str(tmp_path), "BRAKER1")


LOCAL = "pyVARUS Run Summary\n===\nLocal FASTQ libraries: {n}\nMax batches: 1000\nTotal reads in BAM: 42\n"
SRA = "pyVARUS Run Summary\n===\nSpecies: Homo sapiens\nSRA runs in runlist: 3\nMax batches: 1000\n"


def test_local_only_singular(tmp_path):
    t = methods(tmp_path, LOCAL.format(n=1))
    assert "sampled from **1 local FASTQ library** using VARUS" in t
    assert "Sequence Read Archive using VARUS" not in t


def test_local_only_plural(tmp_path):
    t = methods(tmp_path, LOCAL.format(n=3))
    assert "**3 local FASTQ libraries**" in t and "aligned with HISAT2" in t


def test_sra_plus_local(tmp_path):
    t = methods(tmp_path, SRA + "Local FASTQ libraries: 2\n")
    assert "from runs of the NCBI Sequence Read Archive and from **2 local FASTQ libraries**" in t


def test_sra_only(tmp_path):
    t = methods(tmp_path, SRA)
    assert "automatically selected and downloaded from the NCBI Sequence Read Archive" in t
    assert "local FASTQ" not in t


def test_no_stats_file(tmp_path):
    t = methods(tmp_path, None)
    assert "automatically selected and downloaded" in t
    assert "VARUS (Stanke et al., 2019)." in t


def test_no_varus_dir(tmp_path):
    assert "VARUS" not in gr.generate_methods_text(str(tmp_path), "BRAKER1")


# --- varus_sampling_html(): table of the runs and libraries with sampled reads

RUN_STATS = (
    "accesion-id;timesDownloaded;avgUmrPercent;badQuality;maxNumOfBatches\n"
    "SRR1;3;88.5%;False;100\n"
    "SRR2;0;0.0%;False;50\n"
    "SRR3;0;1.2%;True;70\n"
    "liver;5;91.0%;False;8\n"
)
RUNLIST = "@Run_acc\ttotal_spots\ttotal_bases\nSRR1\t5000000\t1\nSRR2\t2500000\t1\nSRR3\t3500000\t1\n"


def results_dir(tmp_path):
    results = tmp_path / "output" / "s" / "results"
    results.mkdir(parents=True)
    return results


def test_sampling_table_sra_and_local(tmp_path):
    results = results_dir(tmp_path)
    (results / "varus_run_statistics.csv").write_text(RUN_STATS)
    (results / "varus_runlist.tsv").write_text(RUNLIST)
    (results / "varus_stats.txt").write_text(SRA)
    html = gr.varus_sampling_html(str(results))
    assert "could sample from 3 SRA runs and 1 local FASTQ library" in html
    assert "took 8 batches of reads from 2 of them" in html
    assert "1 were rejected" in html
    # only the two sources with sampled reads are rows, the library first (5 batches)
    assert html.count("<tr><td>") == 2
    assert html.index("<td>liver</td><td>local FASTQ</td><td>5</td>") < html.index("SRR1</a>")
    assert "sra/SRR1" in html and "SRR2" not in html and "SRR3" not in html
    assert "SRA runs in runlist: 3" in html


def test_sampling_table_local_only_from_varus_directory(tmp_path):
    results = results_dir(tmp_path)
    varus = results.parent / "varus"
    varus.mkdir()
    (varus / "RunStatistics.csv").write_text(
        "accesion-id;timesDownloaded;avgUmrPercent;badQuality;maxNumOfBatches\nreads;7;93.1%;False;8\n")
    html = gr.varus_sampling_html(str(results))
    assert "could sample from 1 local FASTQ library" in html
    assert "<td>reads</td><td>local FASTQ</td><td>7</td><td>8</td><td>93.1%</td>" in html
    assert "ncbi.nlm.nih.gov" not in html and "rejected" not in html


def test_one_batch_is_singular(tmp_path):
    results = results_dir(tmp_path)
    (results / "varus_run_statistics.csv").write_text(
        "accesion-id;timesDownloaded;avgUmrPercent;badQuality;maxNumOfBatches\nmap;1;32.7%;False;1\n")
    html = gr.varus_sampling_html(str(results))
    assert "took 1 batch of reads from 1 of them" in html


def test_run_list_without_statistics_is_not_called_selected(tmp_path):
    results = results_dir(tmp_path)
    (results / "varus_runlist.tsv").write_text(RUNLIST)
    html = gr.varus_sampling_html(str(results))
    assert "SRA Runs Available for Sampling" in html
    assert "selected" not in html
    assert "5,000,000" in html


def test_no_varus_files_no_section(tmp_path):
    assert gr.varus_sampling_html(str(results_dir(tmp_path))) == ""
