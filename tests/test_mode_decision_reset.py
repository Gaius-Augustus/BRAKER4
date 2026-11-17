"""
Tests for the handling of mode decision files across Snakemake calls
(rules/common.smk): mode_evidence_fingerprint, reset_stale_mode_decisions and
the real-run path of _mode_and_decision_files with pending and missing
decisions.

The helper block of rules/common.smk is exec'd as in
test_mode_fallback_helpers.py, here with a pandas samples table and a stub
for Snakemake's `checkpoints` object. All data is synthetic (tmp_path).
"""

import os
import types

import pandas as pd
import pytest
from snakemake.exceptions import IncompleteCheckpointException
from snakemake.io import checkpoint_target

ROOT = os.path.join(os.path.dirname(__file__), "..")
COMMON = os.path.join(ROOT, "rules", "common.smk")

STAGE_FILE = {"rnaseq": "rnaseq_evidence.tsv", "etp": "etp_result.tsv",
              "et": "et_result.tsv", "protein": "protein_evidence.tsv",
              "ep": "ep_result.tsv"}
STAGE_RULE = {"rnaseq": "check_rnaseq_evidence", "etp": "check_etp_result",
              "et": "check_et_result", "protein": "check_protein_evidence",
              "ep": "check_ep_result"}


def _helper_block():
    with open(COMMON) as fh:
        text = fh.read()
    start = text.index("_MODE_STAGE_FILES = {")
    end = text.index("# Input Validation", start)
    return text[start:end]


class FakeCheckpoint:
    """get() raises like Snakemake does for a checkpoint that has not run in
    this call or will run again."""

    def __init__(self, name, pending):
        self.rule = types.SimpleNamespace(name=name)
        self._name = name
        self._pending = pending

    def get(self, **wildcards):
        if self._name in self._pending:
            raise IncompleteCheckpointException(self.rule, checkpoint_target("x"))
        return object()


def make_env(mode, row, dryrun=False, main=True, config=None):
    pending = set()
    cps = types.SimpleNamespace(**{name: FakeCheckpoint(name, pending)
                                   for name in STAGE_RULE.values()})
    ns = {
        "os": os,
        "pd": pd,
        "samples_df": pd.DataFrame([dict(sample_name="s", **row)]),
        "get_braker_mode": lambda s: mode,
        "config": {"auto_mode_fallback": True, "fallback_min_introns": 150} if config is None else config,
        "workflow": types.SimpleNamespace(dryrun=dryrun, is_main_process=main),
        "checkpoints": cps,
    }
    exec(compile(_helper_block(), COMMON, "exec"), ns)
    env = types.SimpleNamespace(**{k: v for k, v in ns.items() if not k.startswith("__")})
    env.pending = pending
    return env


@pytest.fixture
def cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def write_decision(stage, sufficient):
    os.makedirs("output/s/mode_selection", exist_ok=True)
    with open(f"output/s/mode_selection/{STAGE_FILE[stage]}", "w") as fh:
        fh.write(f"sufficient\t{sufficient}\n")


def decisions_on_disk():
    d = "output/s/mode_selection"
    return sorted(f for f in os.listdir(d) if f.endswith(".tsv")) if os.path.isdir(d) else []


def evidence(tmp_path, bam_text="reads"):
    bam = tmp_path / "rnaseq.bam"
    bam.write_text(bam_text)
    prot = tmp_path / "proteins.fa"
    if not prot.exists():
        prot.write_text(">p\nM\n")
    return {"bam_files": str(bam), "protein_fasta": str(prot), "sra_ids": float("nan")}


# ------------------------------------------------------------- fingerprint

def test_fingerprint_changes_with_file_row_and_setting(cwd):
    row = evidence(cwd)
    env = make_env("etp", row)
    first = env.mode_evidence_fingerprint("s")
    assert env.mode_evidence_fingerprint("s") == first
    assert str(row["bam_files"]) in first

    # file content (size) changed
    (cwd / "rnaseq.bam").write_text("more reads than before")
    assert env.mode_evidence_fingerprint("s") != first

    # another library in the row
    other = cwd / "second.bam"
    other.write_text("x")
    env2 = make_env("etp", dict(row, bam_files=f"{row['bam_files']},{other}"))
    assert str(other) in env2.mode_evidence_fingerprint("s")

    # a setting that changes what GeneMark can do
    (cwd / "rnaseq.bam").write_text("reads")
    a = make_env("etp", row, config={"auto_mode_fallback": True, "fallback_min_introns": 150})
    b = make_env("etp", row, config={"auto_mode_fallback": True, "fallback_min_introns": 10})
    assert a.mode_evidence_fingerprint("s") != b.mode_evidence_fingerprint("s")


def test_fingerprint_ignores_other_settings(cwd):
    row = evidence(cwd)
    a = make_env("etp", row, config={"auto_mode_fallback": True, "cpus": 8})
    b = make_env("etp", row, config={"auto_mode_fallback": True, "cpus": 48})
    assert a.mode_evidence_fingerprint("s") == b.mode_evidence_fingerprint("s")


# ------------------------------------------------------------------- reset

def test_first_call_writes_fingerprint_and_keeps_nothing_to_reset(cwd):
    env = make_env("etp", evidence(cwd))
    assert env.reset_stale_mode_decisions("s") is False
    assert os.path.exists("output/s/mode_selection/evidence.fingerprint")


def test_unchanged_evidence_keeps_decisions(cwd):
    env = make_env("etp", evidence(cwd))
    env.reset_stale_mode_decisions("s")
    write_decision("rnaseq", 1)
    write_decision("etp", 0)
    os.makedirs("output/s/results", exist_ok=True)
    open("output/s/results/.done", "w").close()
    assert env.reset_stale_mode_decisions("s") is False
    assert decisions_on_disk() == ["etp_result.tsv", "rnaseq_evidence.tsv"]
    assert os.path.exists("output/s/results/.done")


def test_changed_evidence_removes_decisions_markers_and_done(cwd):
    row = evidence(cwd)
    env = make_env("etp", row)
    env.reset_stale_mode_decisions("s")
    for stage, value in (("rnaseq", 1), ("etp", 0), ("protein", 1), ("ep", 1)):
        write_decision(stage, value)
    open("output/s/mode_selection/etp_insufficient.txt", "w").write("thin\n")
    os.makedirs("output/s/results", exist_ok=True)
    open("output/s/results/.done", "w").close()

    (cwd / "rnaseq.bam").write_text("a better library, much larger")
    assert env.reset_stale_mode_decisions("s") is True
    assert decisions_on_disk() == []
    assert not os.path.exists("output/s/mode_selection/etp_insufficient.txt")
    assert not os.path.exists("output/s/results/.done")
    # the new fingerprint is stored: a second call changes nothing
    assert env.reset_stale_mode_decisions("s") is False


def test_decisions_without_fingerprint_are_kept(cwd):
    """Decision files made before the fingerprint existed are trusted."""
    env = make_env("etp", evidence(cwd))
    write_decision("rnaseq", 1)
    write_decision("etp", 0)
    assert env.reset_stale_mode_decisions("s") is False
    assert decisions_on_disk() == ["etp_result.tsv", "rnaseq_evidence.tsv"]
    assert os.path.exists("output/s/mode_selection/evidence.fingerprint")


def test_dry_run_reports_but_changes_nothing(cwd):
    row = evidence(cwd)
    make_env("etp", row).reset_stale_mode_decisions("s")
    write_decision("rnaseq", 0)
    before = open("output/s/mode_selection/evidence.fingerprint").read()
    (cwd / "rnaseq.bam").write_text("a better library, much larger")
    dry = make_env("etp", row, dryrun=True)
    assert dry.reset_stale_mode_decisions("s") is True
    assert decisions_on_disk() == ["rnaseq_evidence.tsv"]
    assert open("output/s/mode_selection/evidence.fingerprint").read() == before


def test_job_process_and_inactive_samples_do_nothing(cwd):
    row = evidence(cwd)
    # a job process of a cluster executor must not reset anything
    job = make_env("etp", row, main=False)
    assert job.reset_stale_mode_decisions("s") is False
    assert not os.path.exists("output/s/mode_selection")
    # ES, IsoSeq and dual samples have no decisions
    for mode in ("es", "isoseq", "dual"):
        assert make_env(mode, row).reset_stale_mode_decisions("s") is False
    assert not os.path.exists("output/s/mode_selection")
    # switched off
    off = make_env("etp", row, config={"auto_mode_fallback": False})
    assert off.reset_stale_mode_decisions("s") is False
    assert not os.path.exists("output/s/mode_selection")


def test_run_finished_before_the_feature_is_left_alone(cwd):
    os.makedirs("output/s/results")
    open("output/s/results/.done", "w").close()
    env = make_env("etp", evidence(cwd))
    assert env.reset_stale_mode_decisions("s") is False
    assert not os.path.exists("output/s/mode_selection")
    assert os.path.exists("output/s/results/.done")


# ------------------------------------- real-run path: pending and missing

def test_fresh_run_waits_for_each_decision(cwd):
    """No decision file: the input function asks Snakemake to run the check."""
    env = make_env("etp", evidence(cwd))
    env.pending.update(STAGE_RULE.values())
    with pytest.raises(IncompleteCheckpointException):
        env.get_effective_mode("s")
    # rnaseq decided in this call, the next stage has no file yet: wait again
    env.pending.discard("check_rnaseq_evidence")
    write_decision("rnaseq", 1)
    with pytest.raises(IncompleteCheckpointException):
        env.get_effective_mode("s")
    env.pending.discard("check_etp_result")
    write_decision("etp", 1)
    assert env.get_effective_mode("s") == "etp"


def test_missing_later_decision_raises_even_when_get_returns(cwd):
    """A later stage is not in the DAG yet, so get() returns although the
    check never ran: the helper raises itself."""
    env = make_env("etp", evidence(cwd))
    write_decision("rnaseq", 0)        # ETP -> EP, protein decision missing
    with pytest.raises(IncompleteCheckpointException):
        env.get_effective_mode("s")


def test_pending_check_with_old_decision_file_does_not_raise(cwd):
    """Decision files of an earlier call and a check that will run again:
    plan with the old decision, the file is among the inputs."""
    env = make_env("etp", evidence(cwd))
    for stage, value in (("rnaseq", 1), ("etp", 0), ("protein", 1), ("ep", 1)):
        write_decision(stage, value)
    env.pending.update(STAGE_RULE.values())
    assert env.get_effective_mode("s") == "ep"
    files = env.mode_decision_files("s")
    assert [os.path.basename(f) for f in files] == [
        "rnaseq_evidence.tsv", "etp_result.tsv", "protein_evidence.tsv", "ep_result.tsv"]


def test_planned_sample_is_not_unplanned_by_a_changed_decision(cwd):
    """After the decisions were complete once in this process, a decision
    that changes the path must not raise: the missing next decision counts
    as enough for now and its file is an input."""
    env = make_env("etp", evidence(cwd))
    write_decision("rnaseq", 1)
    write_decision("etp", 1)
    assert env.get_effective_mode("s") == "etp"          # planned
    # GeneMark-ETP ran again and now has too little training data
    write_decision("etp", 0)
    assert env.get_effective_mode("s") == "ep"
    files = [os.path.basename(f) for f in env.mode_decision_files("s")]
    assert files == ["rnaseq_evidence.tsv", "etp_result.tsv",
                     "protein_evidence.tsv", "ep_result.tsv"]
    # the protein check ran and says no: ES, and the EP result is not needed
    write_decision("protein", 0)
    assert env.get_effective_mode("s") == "es"
    assert [os.path.basename(f) for f in env.mode_decision_files("s")] == [
        "rnaseq_evidence.tsv", "etp_result.tsv", "protein_evidence.tsv"]


def test_unplanned_sample_still_raises_after_another_sample_was_planned(cwd):
    """The 'planned once' state is per sample."""
    env = make_env("etp", evidence(cwd))
    write_decision("rnaseq", 1)
    write_decision("etp", 1)
    assert env.get_effective_mode("s") == "etp"
    assert "s" in env._MODE_CHAIN_SETTLED
    assert "other" not in env._MODE_CHAIN_SETTLED
