"""
Tests for [PARAMS] varus_fastq: helpers of rules/common.smk (data types, IDs,
mode, validation rule 6b, mode fingerprint) and the get_varus_* helpers of
rules/preprocessing/run_varus.smk.

.smk files cannot be imported. The functions are cut out as text (a top-level
`def name` up to the next top-level statement) and exec'd with a fake
`config` and `samples_df`. Synthetic data only.
"""

import os
import re

import pandas as pd
import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
COMMON = os.path.join(ROOT, "rules", "common.smk")
RUN_VARUS = os.path.join(ROOT, "rules", "preprocessing", "run_varus.smk")


def _cut(path, names):
    """Source of the named top-level functions and of the `_MODE_*` constants."""
    lines = open(path).read().split("\n")
    out = []
    i = 0
    while i < len(lines):
        m = re.match(r"def (\w+)\(", lines[i])
        if m and m.group(1) in names:
            out.append(lines[i])
            i += 1
            while i < len(lines) and (lines[i] == "" or lines[i][0] in " \t)"
                                      or lines[i].startswith("#")):
                out.append(lines[i])
                i += 1
            continue
        i += 1
    return "\n".join(out)


COMMON_FUNCS = ["detect_data_types", "get_braker_mode", "get_ena_settings",
                "validate_samples", "get_bam_files", "get_bam_ids",
                "get_sra_ids", "get_fastq_pairs", "get_fastq_ids",
                "get_varus_fastq_pairs", "get_varus_ids", "get_all_rnaseq_ids",
                "mode_evidence_fingerprint"]


def _constants():
    text = open(COMMON).read()
    start = text.index("_MODE_FINGERPRINT_KEYS = (")
    end = text.index("_MODE_MARKER_FILES")
    return text[start:end]


def load(rows, config=None):
    """Namespace with the helpers of common.smk for samples.csv rows."""
    from pathlib import Path
    df = pd.DataFrame(rows)
    ns = {"pd": pd, "os": os, "re": re, "Path": Path, "samples_df": df,
          "config": {} if config is None else config}
    exec(compile(_constants() + "\n" + _cut(COMMON, COMMON_FUNCS), COMMON, "exec"), ns)
    return ns


def row(**kw):
    base = {"sample_name": "s", "genome": "g.fa", "busco_lineage": "x_odb12"}
    base.update(kw)
    return base


FQ = {"fastq_r1": "a_R1.fastq.gz:b_R1.fq.gz", "fastq_r2": "a_R2.fastq.gz:b_R2.fq.gz"}
SPECIES = {"varus_genus": "Homo", "varus_species": "sapiens"}
PROT = {"protein_fasta": "p.fa"}


# ------------------------------------------------------------ four combinations

def test_fastq_off_hisat2_path():
    ns = load([row(**FQ)], {"varus_fastq": False})
    t = ns["detect_data_types"]("s")
    assert t["has_fastq"] and not t["has_varus"]
    assert not t["has_varus_fastq"] and not t["has_varus_sra"]
    assert ns["get_fastq_ids"]("s") == ["a_R1", "b_R1"] or len(ns["get_fastq_ids"]("s")) == 2
    assert ns["get_varus_ids"]("s") == []
    assert ns["get_varus_fastq_pairs"]("s") == []
    assert len(ns["get_all_rnaseq_ids"]("s")) == 2


def test_fastq_key_absent_is_off():
    ns = load([row(**FQ)])
    assert ns["detect_data_types"]("s")["has_fastq"]


def test_fastq_on_goes_through_varus():
    ns = load([row(**FQ)], {"varus_fastq": True})
    t = ns["detect_data_types"]("s")
    assert not t["has_fastq"]
    assert t["has_varus"] and t["has_varus_fastq"] and not t["has_varus_sra"]
    assert ns["get_fastq_ids"]("s") == []
    assert ns["get_varus_ids"]("s") == ["varus"]
    assert ns["get_varus_fastq_pairs"]("s") == [
        ("a_R1.fastq.gz", "a_R2.fastq.gz"), ("b_R1.fq.gz", "b_R2.fq.gz")]
    assert ns["get_all_rnaseq_ids"]("s") == ["varus"]


def test_fastq_on_with_species_both_sources():
    ns = load([row(**FQ, **SPECIES)], {"varus_fastq": True})
    t = ns["detect_data_types"]("s")
    assert t["has_varus_sra"] and t["has_varus_fastq"] and t["has_varus"]
    assert ns["get_all_rnaseq_ids"]("s") == ["varus"]


def test_fastq_off_with_species_sra_only():
    ns = load([row(**FQ, **SPECIES)], {"varus_fastq": False})
    t = ns["detect_data_types"]("s")
    assert t["has_varus_sra"] and not t["has_varus_fastq"] and t["has_fastq"]
    assert len(ns["get_fastq_ids"]("s")) == 2
    assert ns["get_all_rnaseq_ids"]("s")[-1] == "varus"


def test_no_fastq_flag_on_is_inert():
    ns = load([row(**SPECIES)], {"varus_fastq": True})
    t = ns["detect_data_types"]("s")
    assert t["has_varus_sra"] and not t["has_varus_fastq"] and not t["has_fastq"]
    assert ns["get_varus_fastq_pairs"]("s") == []
    assert ns["get_varus_ids"]("s") == ["varus"]


def test_no_fastq_no_species_flag_on():
    ns = load([row()], {"varus_fastq": True})
    t = ns["detect_data_types"]("s")
    assert not t["has_varus"] and not t["has_fastq"]
    assert ns["get_all_rnaseq_ids"]("s") == []
    assert ns["get_braker_mode"]("s") == "es"


def test_flags_are_bool():
    ns = load([row(**FQ)], {"varus_fastq": 1})
    t = ns["detect_data_types"]("s")
    assert t["has_varus_fastq"] is True
    assert t["has_fastq"] is False


@pytest.mark.parametrize("flag", [True, False])
def test_mode_et_and_etp(flag):
    ns = load([row(sample_name="a", **FQ), row(sample_name="b", **FQ, **PROT)],
              {"varus_fastq": flag})
    assert ns["get_braker_mode"]("a") == "et"
    assert ns["get_braker_mode"]("b") == "etp"


def test_mode_varus_fastq_with_isoseq_is_dual():
    ns = load([row(**FQ, **PROT, isoseq_bam="i.bam")], {"varus_fastq": True})
    assert ns["get_braker_mode"]("s") == "dual"


def test_pairs_mismatch_still_errors():
    ns = load([row(fastq_r1="a.fastq.gz:b.fastq.gz", fastq_r2="a2.fastq.gz")],
              {"varus_fastq": True})
    with pytest.raises(ValueError):
        ns["get_varus_fastq_pairs"]("s")


# ------------------------------------------------------------ rule 6b

def validate(config, **cols):
    ns = load([row(**cols)], config)
    ns["validate_samples"]()


@pytest.mark.parametrize("r1,r2", [
    ("a_R1.fastq.gz", "a_R2.fastq.gz"),
    ("a_1.fq.gz", "a_2.fq.gz"),
    ("a_R1.fastq.gz:b_1.fq.gz", "a_R2.fastq.gz:b_2.fq.gz"),
])
def test_rule6b_accepts_gz(r1, r2):
    validate({"varus_fastq": True}, fastq_r1=r1, fastq_r2=r2)


@pytest.mark.parametrize("r1,r2,bad", [
    ("a_R1.fastq", "a_R2.fastq.gz", "a_R1.fastq"),
    ("a_1.fq", "a_2.fq", "a_1.fq"),
    ("a_R1.fastq.gz", "a_R2.fastq", "a_R2.fastq"),
    ("a_R1.fastq.gz:b_R1.fq", "a_R2.fastq.gz:b_R2.fq.gz", "b_R1.fq"),
    ("a_R1.fastq.bz2", "a_R2.fastq.bz2", "a_R1.fastq.bz2"),
])
def test_rule6b_rejects_plain(r1, r2, bad):
    with pytest.raises(ValueError, match=re.escape(bad)) as e:
        validate({"varus_fastq": True}, fastq_r1=r1, fastq_r2=r2)
    assert "gzip" in str(e.value)


@pytest.mark.parametrize("config", [{}, {"varus_fastq": False}])
def test_rule6b_off_accepts_plain(config):
    validate(config, fastq_r1="a.fastq", fastq_r2="a2.fq")


def test_rule6b_without_fastq_ok():
    validate({"varus_fastq": True}, protein_fasta="p.fa")


# ------------------------------------------------------------ fingerprint

def fp(config, **cols):
    return load([row(**cols)], config)["mode_evidence_fingerprint"]("s")


def test_fingerprint_defaults_identical():
    base = {"fallback_min_introns": 0, "fungus": False}
    a = fp(base, **FQ)
    assert a == fp({**base, "varus_fastq": False, "varus_max_batches": 1000}, **FQ)
    assert "varus" not in a.replace("varus_", "x") or "config\tvarus" not in a
    assert "config\tvarus_fastq" not in a and "config\tvarus_max_batches" not in a
    # the old keys are still there
    assert "config\tfungus\tFalse" in a


def test_fingerprint_varus_fastq_changes():
    base = {"fungus": False}
    a = fp(base, **FQ)
    b = fp({**base, "varus_fastq": True}, **FQ)
    assert a != b
    assert "config\tvarus_fastq\tTrue" in b
    assert b.startswith(a.split("config\t")[0])


def test_fingerprint_max_batches_changes():
    a = fp({}, **FQ)
    b = fp({"varus_max_batches": 500}, **FQ)
    assert a != b and "config\tvarus_max_batches\t500" in b
    assert fp({"varus_max_batches": 500}, **FQ) == b


def test_fingerprint_old_text_unchanged_order():
    # keys in the old order, nothing inserted between them
    cfg = {"fallback_min_introns": 5, "fungus": True, "min_contig": 1,
           "gm_max_intergenic": 2, "translation_table": 1}
    text = fp(cfg)
    tail = [l for l in text.split("\n") if l.startswith("config\t")]
    assert [l.split("\t")[1] for l in tail] == [
        "fallback_min_introns", "fungus", "min_contig", "gm_max_intergenic",
        "translation_table"]


# ------------------------------------------------------------ run_varus.smk

def varus_ns(config, **cols):
    ns = load([row(**cols)], config)
    ns["shlex_"] = None
    exec(compile(_cut(RUN_VARUS, ["get_varus_species_name", "get_varus_fastq_libs"]),
                 RUN_VARUS, "exec"), ns)
    return ns


def test_species_name():
    ns = varus_ns({"varus_fastq": True}, **SPECIES)
    assert ns["get_varus_species_name"]("s") == "Homo sapiens"


def test_species_name_empty_without_sra():
    ns = varus_ns({"varus_fastq": True}, **FQ)
    assert ns["get_varus_species_name"]("s") == ""


def test_species_name_strips():
    ns = varus_ns({}, varus_genus=" Homo ", varus_species="sapiens ")
    assert ns["get_varus_species_name"]("s") == "Homo sapiens"


def test_fastq_libs_words():
    ns = varus_ns({"varus_fastq": True}, **FQ)
    assert ns["get_varus_fastq_libs"]("s") == (
        "a_R1.fastq.gz,a_R2.fastq.gz b_R1.fq.gz,b_R2.fq.gz")


def test_fastq_libs_quoted_and_empty():
    ns = varus_ns({"varus_fastq": True}, fastq_r1="my dir/a.fastq.gz",
                  fastq_r2="my dir/b.fastq.gz")
    assert ns["get_varus_fastq_libs"]("s") == "'my dir/a.fastq.gz,my dir/b.fastq.gz'"
    ns = varus_ns({"varus_fastq": False}, **FQ)
    assert ns["get_varus_fastq_libs"]("s") == ""
