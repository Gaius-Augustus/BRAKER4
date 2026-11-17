"""
Tests for the automatic mode fallback helpers of rules/common.smk and for
read_mode_fallback() of scripts/generate_report.py.

rules/common.smk is Snakemake DSL and cannot be imported. The helper block
(from `_MODE_STAGE_FILES = {` up to the `# Input Validation` banner) is plain
Python, so it is cut out as text and exec'd with stubs for get_braker_mode,
config, workflow and os. The real-run path (checkpoints) is not exercised.
All data is synthetic and written into tmp_path.
"""

import ast
import glob
import os
import re
import types

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
COMMON = os.path.join(ROOT, "rules", "common.smk")


def _helper_block():
    with open(COMMON) as fh:
        text = fh.read()
    start = text.index("_MODE_STAGE_FILES = {")
    end = text.index("# Input Validation", start)
    return text[start:end]  # trailing '# ====' banner line is only a comment


class Env:
    """Namespace with the exec'd helper block."""

    def __init__(self, modes, config=None, dryrun=True):
        self.modes = modes
        ns = {
            "os": os,
            "get_braker_mode": lambda s: modes[s],
            # fallback_min_introns > 0: the two checks before GeneMark are on, so
            # that every decision point is covered; the default 0 has its own
            # tests further down
            "config": {"auto_mode_fallback": True, "fallback_min_introns": 150} if config is None else config,
            "workflow": types.SimpleNamespace(dryrun=dryrun),
        }
        exec(compile(_helper_block(), COMMON, "exec"), ns)
        self.ns = ns

    def __getattr__(self, name):
        return self.ns[name]


@pytest.fixture
def cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def make(mode, **kw):
    return Env({"s": mode}, **kw)


def write_decision(sample, stage, sufficient, reason=None):
    d = f"output/{sample}/mode_selection"
    os.makedirs(d, exist_ok=True)
    name = {"rnaseq": "rnaseq_evidence.tsv", "etp": "etp_result.tsv",
            "et": "et_result.tsv", "protein": "protein_evidence.tsv",
            "ep": "ep_result.tsv"}[stage]
    with open(f"{d}/{name}", "w") as fh:
        fh.write(f"sufficient\t{sufficient}\n")
        if reason:
            fh.write(f"reason\t{reason}\n")


# ------------------------------------------------------------ _walk_mode_chain

def failing(*stages):
    return lambda stage: stage not in stages


@pytest.mark.parametrize("planned, bad, expected", [
    ("etp", [], "etp"),
    ("etp", ["rnaseq"], "ep"),
    ("etp", ["etp"], "ep"),
    ("etp", ["rnaseq", "protein"], "es"),
    ("etp", ["etp", "protein"], "es"),
    ("etp", ["rnaseq", "ep"], "es"),
    ("etp", ["etp", "ep"], "es"),
    ("et", [], "et"),
    ("et", ["rnaseq"], "es"),
    ("et", ["et"], "es"),
    ("ep", [], "ep"),
    ("ep", ["protein"], "es"),
    ("ep", ["ep"], "es"),
])
def test_walk_chain(planned, bad, expected):
    env = make(planned)
    assert env._walk_mode_chain("s", failing(*bad)) == expected


@pytest.mark.parametrize("planned", ["es", "isoseq", "dual"])
def test_walk_chain_unchanged_modes(planned):
    env = make(planned)
    asked = []

    def enough(stage):
        asked.append(stage)
        return False

    assert env._walk_mode_chain("s", enough) == planned
    assert asked == []


# ------------------------------------------------------- effective_mode_on_disk

def test_missing_decisions_count_as_sufficient(cwd):
    env = make("etp")
    assert env.effective_mode_on_disk("s") == "etp"


def test_decision_zero_switches(cwd):
    env = make("etp")
    write_decision("s", "rnaseq", 0)
    assert env.effective_mode_on_disk("s") == "ep"
    write_decision("s", "protein", 0)
    assert env.effective_mode_on_disk("s") == "es"


def test_decision_one_keeps_mode(cwd):
    env = make("et")
    write_decision("s", "rnaseq", 1)
    write_decision("s", "et", 1)
    assert env.effective_mode_on_disk("s") == "et"


def test_decision_without_sufficient_key_is_sufficient(cwd):
    env = make("ep")
    os.makedirs("output/s/mode_selection")
    with open("output/s/mode_selection/protein_evidence.tsv", "w") as fh:
        fh.write("reason\tfoo\n")
    assert env.effective_mode_on_disk("s") == "ep"


def test_et_result_zero_only_after_rnaseq_ok(cwd):
    env = make("et")
    write_decision("s", "et", 0)
    assert env.effective_mode_on_disk("s") == "es"


# --------------------------------------------------------- mode_decision_files

def test_decision_files_dry_run_lists_all_consulted_stages(cwd):
    """In a dry run the list holds the decision files of every stage that the
    walk asked about, whether or not they exist (a missing one is 'enough', so
    the walk continues to the next stage). Stages after a failed one are not
    listed. Planned ETP, nothing on disk: the mode stays etp, so only rnaseq
    and etp are asked."""
    env = make("etp")
    wc = types.SimpleNamespace(sample="s")
    assert env.mode_decision_files(wc) == [
        "output/s/mode_selection/rnaseq_evidence.tsv",
        "output/s/mode_selection/etp_result.tsv",
    ]


def test_decision_files_follow_fallback_path(cwd):
    """After etp -> ep (rnaseq insufficient) the walk asks rnaseq, protein, ep."""
    env = make("etp")
    write_decision("s", "rnaseq", 0)
    wc = types.SimpleNamespace(sample="s")
    assert env.mode_decision_files(wc) == [
        "output/s/mode_selection/rnaseq_evidence.tsv",
        "output/s/mode_selection/protein_evidence.tsv",
        "output/s/mode_selection/ep_result.tsv",
    ]


def test_decision_files_et_and_ep(cwd):
    wc = types.SimpleNamespace(sample="s")
    assert [os.path.basename(p) for p in make("et").mode_decision_files(wc)] == [
        "rnaseq_evidence.tsv", "et_result.tsv"]
    assert [os.path.basename(p) for p in make("ep").mode_decision_files(wc)] == [
        "protein_evidence.tsv", "ep_result.tsv"]


def test_decision_files_accept_plain_sample_string(cwd):
    assert make("ep").mode_decision_files("s")


@pytest.mark.parametrize("planned", ["es", "isoseq", "dual"])
def test_decision_files_empty_for_unaffected_modes(cwd, planned):
    assert make(planned).mode_decision_files("s") == []


# ------------------------------------------------------------ mode_fallback_active

def test_fallback_inactive_when_config_off(cwd):
    env = make("etp", config={"auto_mode_fallback": False})
    assert env.mode_fallback_active("s") is False
    assert env.mode_decision_files("s") == []
    assert env.effective_mode_on_disk("s") == "etp"


def test_fallback_default_is_on(cwd):
    assert make("etp", config={}).mode_fallback_active("s") is True


@pytest.mark.parametrize("planned", ["es", "isoseq", "dual"])
def test_fallback_inactive_for_other_modes(cwd, planned):
    assert make(planned).mode_fallback_active("s") is False


@pytest.mark.parametrize("planned", ["etp", "et", "ep"])
def test_fallback_active_for_planned_modes(cwd, planned):
    assert make(planned).mode_fallback_active("s") is True


def test_fallback_inactive_for_finished_old_run(cwd):
    """results/.done exists and there is no mode_selection directory."""
    os.makedirs("output/s/results")
    open("output/s/results/.done", "w").close()
    env = make("etp")
    assert env.mode_fallback_active("s") is False
    assert env.mode_decision_files("s") == []


def test_fallback_active_for_finished_run_with_decisions(cwd):
    os.makedirs("output/s/results")
    open("output/s/results/.done", "w").close()
    os.makedirs("output/s/mode_selection")
    assert make("etp").mode_fallback_active("s") is True


def test_fallback_active_for_unfinished_run_without_decisions(cwd):
    os.makedirs("output/s/results")
    assert make("etp").mode_fallback_active("s") is True


# ---------------------------------------------------------------- with_mode_decisions

def test_with_mode_decisions_adds_key_without_mutating(cwd):
    env = make("etp")
    inputs = {"a": "x"}
    out = env.with_mode_decisions(types.SimpleNamespace(sample="s"), inputs)
    assert inputs == {"a": "x"}
    assert out is not inputs
    assert out["a"] == "x"
    assert out["mode_decisions"] == env.mode_decision_files("s")


def test_with_mode_decisions_no_decisions_returns_inputs(cwd):
    env = make("es")
    inputs = {"a": "x"}
    out = env.with_mode_decisions(types.SimpleNamespace(sample="s"), inputs)
    assert out == {"a": "x"}
    assert "mode_decisions" not in out


# ------------------------------------------------------- effective_has_transcripts

@pytest.mark.parametrize("planned", ["et", "etp", "isoseq", "dual"])
def test_has_transcripts_true(cwd, planned):
    assert make(planned).effective_has_transcripts("s") is True


@pytest.mark.parametrize("planned", ["ep", "es"])
def test_has_transcripts_false_without_rna(cwd, planned):
    assert make(planned).effective_has_transcripts("s") is False


def test_has_transcripts_false_after_fallback(cwd):
    env = make("etp", dryrun=True)
    write_decision("s", "rnaseq", 0)
    assert env.effective_has_transcripts("s") is False
    env2 = make("et", dryrun=True)
    write_decision("s", "et", 0)
    assert env2.effective_has_transcripts("s") is False


def test_get_effective_mode_dry_run_shows_planned_without_decisions(cwd):
    assert make("etp").get_effective_mode("s") == "etp"


def test_is_dry_run_false_on_broken_workflow(cwd):
    env = make("etp")
    env.ns["workflow"] = object()
    assert env._is_dry_run() is False


# ------------------------------------------------ generate_report.read_mode_fallback

def _load_read_mode_fallback():
    """Load the function by source (ast), so the heavy module imports of
    scripts/generate_report.py are not needed."""
    path = os.path.join(ROOT, "scripts", "generate_report.py")
    with open(path) as fh:
        src = fh.read()
    tree = ast.parse(src)
    node = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "read_mode_fallback")
    ns = {"os": os}
    exec(compile(ast.Module([node], []), path, "exec"), ns)
    return ns["read_mode_fallback"]


rmf = _load_read_mode_fallback()


def _wd(tmp_path, **files):
    d = tmp_path / "mode_selection"
    d.mkdir()
    for name, body in files.items():
        (d / name).write_text(body)
    return str(tmp_path)


def test_report_empty_when_same_mode(tmp_path):
    assert rmf(str(tmp_path), "etp", "etp") == ""


@pytest.mark.parametrize("planned, mode", [("", "es"), ("et", ""), (None, None)])
def test_report_empty_when_unknown(tmp_path, planned, mode):
    assert rmf(str(tmp_path), planned, mode) == ""


def test_report_names_modes_and_reasons(tmp_path):
    wd = _wd(tmp_path,
             **{"rnaseq_evidence.tsv": "sufficient\t0\nreason\tonly 3 RNA-Seq introns\n",
                "protein_evidence.tsv": "sufficient\t1\nreason\tshould not appear\n",
                "ep_result.tsv": "sufficient\t0\nreason\tGeneMark-EP found too few genes\n"})
    note = rmf(wd, "etp", "es")
    assert "ETP" in note and "ES" in note
    assert "only 3 RNA-Seq introns" in note
    assert "GeneMark-EP found too few genes" in note
    assert "should not appear" not in note
    assert note.endswith(".")
    assert "GeneMark-ES trained AUGUSTUS" not in note


def test_report_default_reason_and_no_decisions(tmp_path):
    wd = _wd(tmp_path, **{"et_result.tsv": "sufficient\t0\n"})
    assert "not enough evidence" in rmf(wd, "etp", "ep")
    note = rmf(str(tmp_path / "nowhere"), "ep", "es")
    assert "EP" in note and ":" not in note


def test_report_et_to_es_extra_sentence(tmp_path):
    wd = _wd(tmp_path, **{"rnaseq_evidence.tsv": "sufficient\t0\nreason\tfew introns\n"})
    note = rmf(wd, "et", "es")
    assert "RNA-Seq intron hints were still given to AUGUSTUS" in note
    assert "few introns" in note


# ------------------------------------------------------------ static guard

HELPERS = ("get_effective_mode(", "effective_has_transcripts(", "get_genemark_output",
           "get_genemark_hints_for_filter", "get_genemark_training_gtf",
           "get_augustus_gtf", "get_augustus_hintsfile", "get_extrinsic_cfg")
# the checkpoint file only mentions the helpers in comments
GUARD_EXEMPT = {os.path.join("rules", "genemark", "check_evidence.smk")}


def _code(text):
    """Source without comment lines."""
    return "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))


def _rule_files():
    files = sorted(glob.glob(os.path.join(ROOT, "rules", "**", "*.smk"), recursive=True))
    return [f for f in files if os.path.relpath(f, ROOT) != os.path.join("rules", "common.smk")]


def test_static_guard_scans_files():
    assert len(_rule_files()) > 10


@pytest.mark.parametrize("path", _rule_files(), ids=lambda p: os.path.relpath(p, ROOT))
def test_rules_using_effective_mode_have_decision_inputs(path):
    """A rule that depends on the effective mode needs the decision files as
    input, otherwise Snakemake does not re-evaluate it when a decision
    changes. Checked per file: it must contain mode_decision_files or
    with_mode_decisions( if it uses any mode-dependent helper.

    Intended exception inside merge_hints.smk: for a sample planned as ET the
    decisions are not added (hintsfile.gff is the input of GeneMark-ET, which
    would make a cycle); see test_merge_hints_et_exception."""
    rel = os.path.relpath(path, ROOT)
    if rel in GUARD_EXEMPT:
        pytest.skip("only mentions helpers in comments")
    with open(path) as fh:
        code = _code(fh.read())
    uses = [h for h in HELPERS if h in code]
    if not uses:
        pytest.skip("no mode-dependent helper used")
    assert "mode_decision_files" in code or "with_mode_decisions(" in code, (
        f"{rel} uses {uses} but never adds the mode decision files to its inputs")


def test_guard_covers_expected_files():
    """Make sure the scan really sees the rules that are known to depend on the mode."""
    hits = set()
    for p in _rule_files():
        with open(p) as fh:
            code = _code(fh.read())
        if any(h in code for h in HELPERS):
            hits.add(os.path.basename(p))
    assert {"run_tsebra.smk", "add_utr.smk", "convert_to_gff3.smk",
            "run_augustus_hints.smk", "merge_hints.smk"} <= hits


def test_merge_hints_et_exception():
    path = os.path.join(ROOT, "rules", "postprocessing", "merge_hints.smk")
    with open(path) as fh:
        code = _code(fh.read())
    assert re.search(r"if get_braker_mode\(sample\) != 'et':\s*\n\s*.*\n?.*with_mode_decisions\(", code) \
        or ("get_braker_mode(sample) != 'et'" in code and "with_mode_decisions(" in code)


def test_merge_hints_shell_uses_params_hint_files():
    path = os.path.join(ROOT, "rules", "postprocessing", "merge_hints.smk")
    with open(path) as fh:
        code = _code(fh.read())
    assert "hint_files=" in code
    assert "mode_decisions" in code
    assert "{params.hint_files}" in code
    # no loop over the bare {input}, which would include the decision files
    assert not re.search(r"for\s+\w+\s+in\s+\{input\}", code)
    assert not re.search(r"\{input\}\s*;?\s*do", code)


# --------------------------- default fallback_min_introns = 0: no pre-checks

DEFAULT = {"auto_mode_fallback": True, "fallback_min_introns": 0}


@pytest.mark.parametrize("planned, bad, expected, asked", [
    ("etp", [], "etp", ["etp"]),
    ("etp", ["rnaseq", "protein"], "etp", ["etp"]),      # the checks before GeneMark are not asked
    ("etp", ["etp"], "ep", ["etp", "ep"]),
    ("etp", ["etp", "ep"], "es", ["etp", "ep"]),
    ("et", ["rnaseq"], "et", ["et"]),
    ("et", ["et"], "es", ["et"]),
    ("ep", ["protein"], "ep", ["ep"]),
    ("ep", ["ep"], "es", ["ep"]),
])
def test_walk_chain_without_prechecks(planned, bad, expected, asked):
    """With the default 0 only the results of the GeneMark runs decide."""
    for config in (DEFAULT, {"auto_mode_fallback": True}):
        env = make(planned, config=config)
        seen = []

        def enough(stage):
            seen.append(stage)
            return stage not in bad

        assert env._walk_mode_chain("s", enough) == expected
        assert seen == asked


def test_decision_files_without_prechecks(cwd):
    env = make("etp", config=DEFAULT, dryrun=True)
    assert [os.path.basename(f) for f in env.mode_decision_files("s")] == ["etp_result.tsv"]
    # a pre-check file of an earlier call with another setting is ignored
    write_decision("s", "rnaseq", 0)
    assert env.effective_mode_on_disk("s") == "etp"
    write_decision("s", "etp", 0)
    assert env.effective_mode_on_disk("s") == "ep"
    assert [os.path.basename(f) for f in env.mode_decision_files("s")] == ["etp_result.tsv", "ep_result.tsv"]
