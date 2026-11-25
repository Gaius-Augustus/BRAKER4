"""
Common functions and configuration for the BRAKER4 Snakemake workflow.

This module provides:
- Sample parsing from samples.csv
- Data type detection (BAM, FASTQ, proteins, etc.)
- Input validation
- Helper functions for dynamic input resolution
"""

import pandas as pd
import os
import re
from pathlib import Path

# ==============================================================================
# Configuration and Samples
# ==============================================================================

# Read samples.csv
samples_csv = config.get("samples_file", "samples.csv")
if not os.path.isfile(samples_csv):
    raise FileNotFoundError(
        f"Samples file '{samples_csv}' not found. "
        "Please create it or set samples_file in config.ini [paths]. "
        "See the README for the expected CSV format."
    )
samples_df = pd.read_csv(samples_csv)

# Get list of all samples
SAMPLES = samples_df["sample_name"].tolist()

# ==============================================================================
# Data Type Detection
# ==============================================================================

def detect_data_types(sample):
    """
    Detect which data types are present for a sample.

    Returns dict with flags for each data type.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]

    # [PARAMS] varus_fastq = 1 sends the local FASTQ libraries through pyVARUS
    # (sampling) instead of HISAT2 (all reads). has_fastq then is False, so
    # that everything keyed on it (HISAT2 rules, hint IDs) leaves them alone.
    fastq = pd.notna(row.get("fastq_r1")) and pd.notna(row.get("fastq_r2"))
    varus_sra = pd.notna(row.get("varus_genus")) and pd.notna(row.get("varus_species"))
    varus_fastq = fastq and bool(config.get("varus_fastq", False))

    return {
        "has_proteins": pd.notna(row.get("protein_fasta")),
        "has_bam": pd.notna(row.get("bam_files")),
        "has_fastq": fastq and not varus_fastq,
        "has_sra": pd.notna(row.get("sra_ids")),
        "has_varus": varus_sra or varus_fastq,
        "has_varus_sra": varus_sra,
        "has_varus_fastq": varus_fastq,
        "has_isoseq": pd.notna(row.get("isoseq_bam")) or pd.notna(row.get("isoseq_fastq")),
        "has_isoseq_fastq": pd.notna(row.get("isoseq_fastq")),
        "has_masked_genome": pd.notna(row.get("genome_masked")),
        "needs_masking": pd.isna(row.get("genome_masked")),
        "has_reference_gtf": pd.notna(row.get("reference_gtf"))
    }

def get_braker_mode(sample):
    """
    Determine BRAKER mode for a sample.

    Returns one of es, et, ep, etp, isoseq, dual:

    - ES: Genome only (no RNA-Seq, no proteins, no IsoSeq)
    - ET: RNA-Seq only (no proteins)
    - EP: Proteins only (no RNA-Seq)
    - ETP: RNA-Seq + Proteins
    - IsoSeq: PacBio long-read + Proteins (no short-read RNA-Seq)
    - dual: Short-read RNA-Seq + IsoSeq + Proteins (one combined ETP run with both)
    """
    types = detect_data_types(sample)

    has_rnaseq = any([types["has_bam"], types["has_fastq"],
                      types["has_sra"], types["has_varus"]])

    if types["has_isoseq"]:
        if not types["has_proteins"]:
            raise ValueError(f"Sample {sample}: IsoSeq REQUIRES protein evidence!")
        if has_rnaseq:
            return "dual"
        return "isoseq"
    elif has_rnaseq and types["has_proteins"]:
        return "etp"
    elif has_rnaseq:
        return "et"
    elif types["has_proteins"]:
        return "ep"
    else:
        return "es"


# ==============================================================================
# Automatic mode fallback (ETP -> EP -> ES, ET -> ES)
# ==============================================================================
#
# get_braker_mode() is the mode PLANNED from samples.csv; it is known when
# the Snakefile is parsed and decides which rule files are included. The
# mode a sample really runs in (the EFFECTIVE mode) can be a simpler one
# when the evidence turns out to be too thin ([PARAMS] auto_mode_fallback,
# default on). Checkpoints (rules/genemark/check_evidence.smk) write one
# decision file each to output/<sample>/mode_selection/:
#
#   etp      etp_result.tsv        GeneMark-ETP found training genes (ETP)
#   et       et_result.tsv         GeneMark-ET could train (ET)
#   ep       ep_result.tsv         GeneMark-EP could train (EP, also after
#                                  a fallback from ETP)
#   and only with fallback_min_introns > 0 (default 0), before GeneMark:
#   rnaseq   rnaseq_evidence.tsv   RNA-Seq intron hints (planned ET, ETP)
#   protein  protein_evidence.tsv  ProtHint intron hints (EP)
#
# hintsfile.gff of a planned ET sample is the input of GeneMark-ET, so
# merge_hints cannot depend on the decisions for such a sample: after a
# fallback from ET to ES, AUGUSTUS still gets the RNA-Seq intron hints.
#
# Rules downstream of GeneMark must use get_effective_mode() in their input
# functions and effective_mode_on_disk() in params (with
# mode_decision_files as an input, so that the decisions exist when the
# params are evaluated). Code that runs while the Snakefile is parsed, and
# rules upstream of the checkpoints, must keep get_braker_mode().

_MODE_STAGE_FILES = {
    "rnaseq":  "rnaseq_evidence.tsv",
    "etp":     "etp_result.tsv",
    "et":      "et_result.tsv",
    "protein": "protein_evidence.tsv",
    "ep":      "ep_result.tsv",
}
_MODE_STAGE_CHECKPOINTS = {
    "rnaseq":  "check_rnaseq_evidence",
    "etp":     "check_etp_result",
    "et":      "check_et_result",
    "protein": "check_protein_evidence",
    "ep":      "check_ep_result",
}


def mode_stage_file(sample, stage):
    """Path of the decision file of one stage."""
    return f"output/{sample}/mode_selection/{_MODE_STAGE_FILES[stage]}"


def read_mode_decision(path):
    """Read a decision file (key<TAB>value lines) into a dict."""
    entries = {}
    with open(path) as handle:
        for line in handle:
            key, _, value = line.rstrip("\n").partition("\t")
            if key:
                entries[key] = value
    return entries


def mode_fallback_active(sample):
    """True when the effective mode of the sample is decided at run time."""
    if not config.get("auto_mode_fallback", True):
        return False
    if get_braker_mode(sample) not in ("etp", "et", "ep"):
        return False
    # A run that finished before the fallback existed has no decision
    # files and its intermediate files are gone; asking for the decisions
    # would start it again from the alignments.
    if (os.path.exists(f"output/{sample}/results/.done")
            and not os.path.isdir(f"output/{sample}/mode_selection")):
        return False
    return True


# Settings that change what GeneMark can do with the evidence of a sample
_MODE_FINGERPRINT_KEYS = ("fallback_min_introns", "fungus", "min_contig",
                          "gm_max_intergenic", "translation_table",
                          "varus_fastq", "varus_max_batches")
# Settings added after the fingerprint existed: left out of it at their
# default, so that the fingerprint of an older finished run stays the same.
_MODE_FINGERPRINT_DEFAULTS = {"varus_fastq": False, "varus_max_batches": 1000}
_MODE_MARKER_FILES = ("etp_insufficient.txt", "et_insufficient.txt", "ep_insufficient.txt")
# samples whose decisions were complete at least once in this process
_MODE_CHAIN_SETTLED = set()


def mode_evidence_fingerprint(sample):
    """What the mode decisions of a sample rest on, as text: its row of
    samples.csv, size and modification time of every file named in it, and
    the settings in _MODE_FINGERPRINT_KEYS."""
    import re
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    lines = []
    for column in samples_df.columns:
        value = row.get(column)
        value = "" if pd.isna(value) else str(value)
        lines.append(f"{column}\t{value}")
        for token in re.split(r"[,:;\s]+", value):
            if token and os.path.isfile(token):
                st = os.stat(token)
                lines.append(f"file\t{token}\t{st.st_size}\t{st.st_mtime_ns}")
    for key in _MODE_FINGERPRINT_KEYS:
        if key in _MODE_FINGERPRINT_DEFAULTS and \
                config.get(key, _MODE_FINGERPRINT_DEFAULTS[key]) == _MODE_FINGERPRINT_DEFAULTS[key]:
            continue
        lines.append(f"config\t{key}\t{config.get(key, '')}")
    return "\n".join(lines) + "\n"


def reset_stale_mode_decisions(sample):
    """Remove the decision files of a sample when its evidence changed.

    Called once per sample when the workflow starts. Decision files of an
    earlier Snakemake call are only valid for the evidence they were made
    from. When a file named in samples.csv, the row itself or one of the
    settings in _MODE_FINGERPRINT_KEYS changed (a library was added, the
    alignments were redone), Snakemake computes the sample again, and the
    old decisions must not steer that run: the decision files are removed,
    and all decisions are made again as in a new run. results/.done of the
    sample is removed as well, because Snakemake takes a job whose inputs
    wait for a checkpoint as up to date when its output exists. The fingerprint is
    kept in output/<sample>/mode_selection/evidence.fingerprint. Decision
    files without a fingerprint (made before it existed) are kept.

    Only the main Snakemake process does this, not the job processes of a
    cluster executor, and a dry run only reports it. Returns True when
    decisions were (or, in a dry run, would be) removed."""
    import sys
    if not mode_fallback_active(sample):
        return False
    if not getattr(workflow, "is_main_process", True):
        return False
    directory = f"output/{sample}/mode_selection"
    fingerprint_path = f"{directory}/evidence.fingerprint"
    current = mode_evidence_fingerprint(sample)
    previous = None
    if os.path.exists(fingerprint_path):
        with open(fingerprint_path) as handle:
            previous = handle.read()
    if previous == current:
        return False
    old_files = [os.path.join(directory, name)
                 for name in list(_MODE_STAGE_FILES.values()) + list(_MODE_MARKER_FILES)
                 if os.path.exists(os.path.join(directory, name))]
    reset = previous is not None and bool(old_files)
    if _is_dry_run():
        if reset:
            print(f"NOTE: the evidence of sample {sample} changed since its mode decisions "
                  f"were made; a real run makes them again.", file=sys.stderr)
        return reset
    if reset:
        for path in old_files + [f"output/{sample}/results/.done"]:
            if os.path.exists(path):
                os.remove(path)
        print(f"NOTE: the evidence of sample {sample} changed since its mode decisions "
              f"were made ({directory}/); they are made again in this run.", file=sys.stderr)
    os.makedirs(directory, exist_ok=True)
    with open(fingerprint_path, "w") as handle:
        handle.write(current)
    return reset


def _is_dry_run():
    try:
        return bool(workflow.dryrun)
    except Exception:
        return False


def _walk_mode_chain(sample, enough):
    """Effective mode of a sample; enough(stage) answers one decision.

    The results of GeneMark-ETP, GeneMark-ET and GeneMark-EP are always
    asked. The two checks before GeneMark (RNA-Seq and protein intron
    hints) are only asked with fallback_min_introns > 0; the default 0
    means that only a GeneMark run that really failed changes the mode."""
    planned = get_braker_mode(sample)
    mode = planned
    precheck = int(config.get("fallback_min_introns", 0) or 0) > 0
    if precheck and mode in ("etp", "et") and not enough("rnaseq"):
        mode = "ep" if planned == "etp" else "es"
    if mode == "etp" and not enough("etp"):
        mode = "ep"
    if mode == "et" and not enough("et"):
        mode = "es"
    if precheck and mode == "ep" and not enough("protein"):
        mode = "es"
    if mode == "ep" and not enough("ep"):
        mode = "es"
    return mode


def _mode_and_decision_files(sample, from_disk):
    """(effective mode, decision files it rests on).

    from_disk=False asks the checkpoints, which makes Snakemake run them
    and evaluate the calling input function again afterwards.
    from_disk=True only reads decision files that exist and takes a missing
    one as 'enough evidence': used for params, and for dry runs, which
    therefore show the planned mode until the decisions are made.
    """
    if not mode_fallback_active(sample):
        return get_braker_mode(sample), []
    files = []
    # Once the decisions of a sample were complete in this process, its
    # jobs are planned. From then on an open decision (a check that runs
    # again and changes the path, so that the next decision file does not
    # exist yet) must not raise IncompleteCheckpointException any more:
    # that would turn planned jobs into jobs with unknown inputs, which
    # Snakemake 8 keeps in its job cache and later runs in that state, and
    # it treats such a job as up to date when its old output exists.
    # Instead the missing decision counts as 'enough' for now and its file
    # is among the inputs, so the job waits for the check and is evaluated
    # again afterwards.
    settled = (not from_disk) and sample in _MODE_CHAIN_SETTLED

    def enough(stage):
        path = mode_stage_file(sample, stage)
        files.append(path)
        if from_disk:
            if not os.path.exists(path):
                return True
        else:
            from snakemake.exceptions import IncompleteCheckpointException
            from snakemake.io import checkpoint_target
            stage_cp = getattr(checkpoints, _MODE_STAGE_CHECKPOINTS[stage])
            pending = None
            try:
                stage_cp.get(sample=sample)
            except IncompleteCheckpointException as error:
                # the check has not run in this Snakemake call, or it will
                # run again
                pending = error
            if not os.path.exists(path):
                if settled:
                    # see above: stay planned, the decision file is an input
                    return True
                if pending is not None:
                    raise pending
                # get() only knows the checkpoint jobs that are in the DAG
                # already. A later stage enters the DAG when an earlier
                # decision asks for it, and get() then returns although
                # the stage has not run: ask for it in the same way.
                raise IncompleteCheckpointException(stage_cp.rule, checkpoint_target(path))
            # The file exists. If the check is pending (decision file of an
            # earlier call, check runs again), plan with the old decision:
            # the file is an input of every rule that depends on the mode
            # (mode_decision_files), so such a job cannot start before the
            # check has run again, and Snakemake evaluates it again then.
        return read_mode_decision(path).get("sufficient", "1") == "1"

    mode = _walk_mode_chain(sample, enough)
    if not from_disk:
        _MODE_CHAIN_SETTLED.add(sample)
    return mode, files


def get_effective_mode(sample):
    """Mode the sample really runs in. Only for input functions."""
    return _mode_and_decision_files(sample, from_disk=_is_dry_run())[0]


def effective_mode_on_disk(sample):
    """Effective mode for params: decision files only, no checkpoint call.
    The rule needs mode_decision_files among its inputs."""
    return _mode_and_decision_files(sample, from_disk=True)[0]


def mode_decision_files(wildcards):
    """Input function: the decision files the effective mode rests on."""
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    return _mode_and_decision_files(sample, from_disk=_is_dry_run())[1]


def with_mode_decisions(wildcards, inputs):
    """inputs (a dict for unpack()) plus the decision files as
    'mode_decisions'. Every rule whose inputs or params depend on the
    effective mode needs them as inputs: only then does Snakemake evaluate
    the rule again when a decision changes, e.g. on a second run in the same
    directory after the evidence was improved."""
    decisions = mode_decision_files(wildcards)
    if decisions:
        inputs = dict(inputs)
        inputs["mode_decisions"] = decisions
    return inputs


def effective_has_transcripts(sample):
    """Transcript evidence in use (UTRs, lncRNA): not after a fallback to
    EP or ES mode. Only for input functions."""
    return get_effective_mode(sample) in ("et", "etp", "isoseq", "dual")


# ==============================================================================
# Input Validation
# ==============================================================================

def get_ena_settings(sample):
    """ENA export columns of a sample (optional): locus_tag prefix, taxid,
    seqid map. Empty strings when not set. pandas reads a taxid column with
    gaps as float, so 3702.0 is turned back into '3702'."""
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]

    def cell(name):
        value = row.get(name)
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return ""
        if isinstance(value, float) and value.is_integer():
            value = int(value)
        return str(value).strip()

    return {
        "locus_tag_prefix": cell("ena_locus_tag_prefix"),
        "taxid": cell("ena_taxid"),
        "seqid_map": cell("ena_seqid_map"),
    }


def validate_samples():
    """Validate sample configurations before workflow execution."""
    for idx, row in samples_df.iterrows():
        sample = row["sample_name"]

        # Rule 1: Genome is required
        if pd.isna(row.get("genome")):
            raise ValueError(f"{sample}: genome column is required")

        # Rule 2: No evidence = ES mode (ab initio), which is valid but warn
        has_evidence = any([
            pd.notna(row.get("protein_fasta")),
            pd.notna(row.get("bam_files")),
            pd.notna(row.get("fastq_r1")),
            pd.notna(row.get("sra_ids")),
            pd.notna(row.get("varus_genus")),
            pd.notna(row.get("isoseq_bam")),
            pd.notna(row.get("isoseq_fastq"))
        ])
        if not has_evidence:
            print(f"  WARNING: {sample} has no evidence — will run in ES mode (ab initio only)")

        # Rule 3: IsoSeq requires proteins
        has_isoseq = pd.notna(row.get("isoseq_bam")) or pd.notna(row.get("isoseq_fastq"))
        if has_isoseq and pd.isna(row.get("protein_fasta")):
            raise ValueError(f"{sample}: IsoSeq REQUIRES protein evidence")

        # Rule 5: FASTQ requires both R1 and R2
        if pd.notna(row.get("fastq_r1")) != pd.notna(row.get("fastq_r2")):
            raise ValueError(f"{sample}: FASTQ requires both R1 and R2")

        # Rule 6: VARUS requires genus and species
        if pd.notna(row.get("varus_genus")) != pd.notna(row.get("varus_species")):
            raise ValueError(f"{sample}: VARUS requires both genus and species")

        # Rule 6b: pyVARUS samples from gzip-compressed FASTQ only
        if config.get("varus_fastq", False) and pd.notna(row.get("fastq_r1")) and pd.notna(row.get("fastq_r2")):
            for fq in str(row["fastq_r1"]).split(":") + str(row["fastq_r2"]).split(":"):
                if not fq.endswith((".fastq.gz", ".fq.gz")):
                    raise ValueError(
                        f"{sample}: varus_fastq = 1 needs gzip-compressed FASTQ files "
                        f"named *.fastq.gz or *.fq.gz, got {fq}. Compress the file "
                        f"with gzip or set varus_fastq = 0.")

        # Rule 7: busco_lineage is required
        if pd.isna(row.get("busco_lineage")) or not str(row.get("busco_lineage")).strip():
            raise ValueError(f"{sample}: busco_lineage column is required in samples.csv")

        # Rule 9: ENA export settings (optional columns)
        ena = get_ena_settings(sample)
        if ena["locus_tag_prefix"]:
            # same rule as gff3tools LocusTagExistsValidation
            if not re.match(r"^[A-Z][A-Z0-9]{2,11}$", ena["locus_tag_prefix"]):
                raise ValueError(
                    f"{sample}: ena_locus_tag_prefix '{ena['locus_tag_prefix']}' is invalid. "
                    f"ENA needs the registered prefix: 3-12 upper-case letters or digits, "
                    f"starting with a letter.")
            if ena["taxid"] and not ena["taxid"].isdigit():
                raise ValueError(f"{sample}: ena_taxid must be an NCBI taxonomy ID (a number), "
                                 f"got '{ena['taxid']}'")
            if ena["seqid_map"] and not os.path.isfile(ena["seqid_map"]):
                raise ValueError(f"{sample}: ena_seqid_map file not found: {ena['seqid_map']}")
        elif ena["taxid"] or ena["seqid_map"]:
            print(f"  WARNING: {sample} has ena_taxid/ena_seqid_map but no ena_locus_tag_prefix; "
                  f"no ENA GFF3 is written")

        # Rule 8: Alternative translation tables only supported in ES/ET/EP modes
        tt = config.get('translation_table', 1)
        if tt != 1:
            mode = get_braker_mode(sample)
            if mode in ('etp', 'isoseq', 'dual'):
                raise ValueError(
                    f"{sample}: translation_table={tt} is not supported in {mode.upper()} mode. "
                    f"GeneMark-ETP and the GeneMarkS-T it runs on the assembled transcripts "
                    f"always use the standard code (gmst.pl supports codes 1, 4 and 11 only). "
                    f"Use ES, ET, or EP mode instead, or set translation_table = 1."
                )

# Run validation at workflow start
validate_samples()

# ==============================================================================
# Helper Functions
# ==============================================================================

def get_masked_genome(sample):
    """
    Get the masked genome file path for a sample.

    All paths point to header-cleaned files in the output directory.
    Returns:
    - Pre-masked genome (cleaned): output/{sample}/genome_masked.fa
    - RepeatMasker output: output/{sample}/preprocessing/genome.fa.masked
    - Cleaned unmasked genome: output/{sample}/genome.fa
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]

    if pd.notna(row.get("genome_masked")):
        return f"output/{sample}/genome_masked.fa"

    if GLOBAL_DATA_TYPES["needs_masking"]:
        return f"output/{sample}/preprocessing/genome.fa.masked"

    return f"output/{sample}/genome.fa"

def get_genome(sample):
    """Get the cleaned genome file path (headers simplified to accession only)."""
    return f"output/{sample}/genome.fa"

def get_raw_genome(sample):
    """Get the original user-provided genome path (before header cleaning)."""
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    return row["genome"]

def get_bam_files(sample):
    """
    Get list of BAM file paths for a sample.

    Returns empty list if no BAM files.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]

    if pd.notna(row.get("bam_files")):
        # BAM files are colon-separated
        return row["bam_files"].split(":")
    else:
        return []

def get_bam_ids(sample):
    """
    Get list of BAM IDs for a sample (for wildcards).

    Returns base names without extension.
    """
    bam_files = get_bam_files(sample)
    return [Path(bam).stem for bam in bam_files]

def get_sra_ids(sample):
    """
    Get list of SRA accession IDs for a sample.

    Returns empty list if no SRA IDs specified.
    SRA IDs are colon-separated in samples.csv.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    if pd.notna(row.get("sra_ids")):
        return str(row["sra_ids"]).split(":")
    return []

def get_fastq_pairs(sample):
    """
    Get list of (r1, r2) FASTQ file path tuples for a sample.

    Supports colon-separated paths for multiple pairs.
    Returns empty list if no FASTQ files specified.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    if pd.notna(row.get("fastq_r1")) and pd.notna(row.get("fastq_r2")):
        r1_files = str(row["fastq_r1"]).split(":")
        r2_files = str(row["fastq_r2"]).split(":")
        if len(r1_files) != len(r2_files):
            raise ValueError(
                f"Sample {sample}: fastq_r1 and fastq_r2 must have "
                f"the same number of colon-separated files"
            )
        return list(zip(r1_files, r2_files))
    return []

def get_fastq_ids(sample):
    """
    Get IDs for FASTQ pairs (derived from R1 filenames).

    Strips extensions and _R1/_1 suffixes to produce a clean ID.
    E.g. reads_R1.fastq.gz -> reads

    These are the libraries that HISAT2 aligns in full. With
    [PARAMS] varus_fastq = 1 the list is empty: the libraries are sampled by
    pyVARUS and arrive as the one VARUS BAM (get_varus_ids,
    get_varus_fastq_pairs).
    """
    if detect_data_types(sample)["has_varus_fastq"]:
        return []
    pairs = get_fastq_pairs(sample)
    ids = []
    for r1, _ in pairs:
        name = Path(r1).name
        # Strip all extensions (.fastq.gz, .fq.gz, etc.)
        stem = name.split('.')[0]
        # Remove R1/1 suffix
        for suffix in ['_R1', '_1', '_r1']:
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]
                break
        ids.append(stem)
    return ids

def get_fastq_r1(sample, fastq_id):
    """Get R1 FASTQ path for a given fastq_id."""
    pairs = get_fastq_pairs(sample)
    ids = get_fastq_ids(sample)
    for fid, (r1, _) in zip(ids, pairs):
        if fid == fastq_id:
            return r1
    raise ValueError(f"No FASTQ R1 found for fastq_id {fastq_id} in sample {sample}")

def get_fastq_r2(sample, fastq_id):
    """Get R2 FASTQ path for a given fastq_id."""
    pairs = get_fastq_pairs(sample)
    ids = get_fastq_ids(sample)
    for fid, (_, r2) in zip(ids, pairs):
        if fid == fastq_id:
            return r2
    raise ValueError(f"No FASTQ R2 found for fastq_id {fastq_id} in sample {sample}")

def get_isoseq_fastq_files(sample):
    """Get list of IsoSeq FASTQ/FASTA file paths for a sample (unaligned reads).

    Supports colon-separated paths for multiple files.
    Returns empty list if no IsoSeq FASTQ files specified.
    The column is called isoseq_fastq for consistency with the test data
    (PacBio HiFi reads are typically distributed as FASTQ); minimap2
    accepts FASTA too, so the same column path may be a .fa file.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    if pd.notna(row.get("isoseq_fastq")):
        return str(row["isoseq_fastq"]).split(":")
    return []

def get_isoseq_fastq_ids(sample):
    """Get IDs for IsoSeq FASTQ files (derived from filenames).

    E.g. isoseq_lib1.fastq.gz -> isoseq_lib1
    """
    files = get_isoseq_fastq_files(sample)
    ids = []
    for f in files:
        name = Path(f).name.split('.')[0]
        ids.append(name)
    return ids

def get_isoseq_fastq(sample):
    """Get single IsoSeq FASTQ/FASTA file path (backwards compat).

    For rules that need a single path. Returns first file or None.
    """
    files = get_isoseq_fastq_files(sample)
    return files[0] if files else None

def get_isoseq_bam_files(sample):
    """Get list of pre-aligned IsoSeq BAM file paths for a sample.

    Supports colon-separated paths for multiple files.
    Returns empty list if no IsoSeq BAM files specified.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    if pd.notna(row.get("isoseq_bam")):
        return str(row["isoseq_bam"]).split(":")
    return []

def get_isoseq_bam_ids(sample):
    """Get IDs for IsoSeq BAM files (derived from filenames).

    E.g. isoseq_lib1.bam -> isoseq_lib1
    """
    files = get_isoseq_bam_files(sample)
    return [Path(f).stem for f in files]

def get_isoseq_sorted_bams(sample):
    """Get list of individually sorted IsoSeq BAM paths (before merging).

    Used as inputs to the merge rule.
    """
    bams = []
    for bid in get_isoseq_bam_ids(sample):
        bams.append(f"output/{sample}/isoseq_sorted/{bid}.sorted.bam")
    for fid in get_isoseq_fastq_ids(sample):
        bams.append(f"output/{sample}/minimap2_aligned/{fid}.sorted.bam")
    return bams

def get_isoseq_bam_for_etp(sample):
    """Get the single (possibly merged) IsoSeq BAM path for GeneMark-ETP input.

    IsoSeq BAMs must be merged before passing to ETP.
    If multiple IsoSeq inputs, returns the merged BAM path.
    If single input, returns the sorted BAM directly.
    """
    sorted_bams = get_isoseq_sorted_bams(sample)
    if len(sorted_bams) > 1:
        return f"output/{sample}/isoseq_merged/isoseq.merged.bam"
    elif len(sorted_bams) == 1:
        return sorted_bams[0]
    return None

def get_varus_fastq_pairs(sample):
    """(r1, r2) tuples of the local FASTQ libraries that pyVARUS samples
    from ([PARAMS] varus_fastq = 1); empty list otherwise."""
    if detect_data_types(sample)["has_varus_fastq"]:
        return get_fastq_pairs(sample)
    return []

def get_varus_ids(sample):
    """
    Get VARUS RNA-seq ID for a sample.

    Returns ["varus"] if pyVARUS runs for the sample (SRA sampling with
    varus_genus and varus_species, and/or sampling of the local FASTQ
    libraries with varus_fastq = 1), empty list otherwise.
    """
    types = detect_data_types(sample)
    if types["has_varus"]:
        return ["varus"]
    return []

def get_isoseq_ids(sample):
    """
    Get IsoSeq BAM ID for a sample.

    Returns ["isoseq"] if IsoSeq data is present, empty list otherwise.
    IsoSeq is NOT part of get_all_rnaseq_ids: it does not go through
    bam2hints / filter_introns / GeneMark-ET but through GeneMark-ETP
    (see get_isoseq_bam_for_etp).
    """
    types = detect_data_types(sample)
    if types["has_isoseq"]:
        return ["isoseq"]
    return []

def get_all_rnaseq_ids(sample):
    """
    Get all short-read RNA-Seq evidence IDs for a sample.

    Combines IDs from pre-existing BAM files, SRA downloads, FASTQ files,
    and VARUS. These IDs are used to expand bam2hints and join_hints rules.
    IsoSeq is NOT included here — it goes through GeneMark-ETP, not bam2hints.
    """
    return get_bam_ids(sample) + get_sra_ids(sample) + get_fastq_ids(sample) + get_varus_ids(sample)

def get_rnaseq_bam(wildcards):
    """
    Get sorted BAM + index for any RNA-seq evidence ID.

    Routes to the correct BAM path depending on data source:
    - Pre-existing BAM -> output/{sample}/bam_sorted/{id}.sorted.bam
    - SRA or FASTQ -> output/{sample}/hisat2_aligned/{id}.sorted.bam
    - VARUS -> output/{sample}/varus/varus.sorted.bam
    """
    sample = wildcards.sample
    bam_id = wildcards.bam_id

    if bam_id in get_bam_ids(sample):
        return {
            "bam": f"output/{sample}/bam_sorted/{bam_id}.sorted.bam",
            "csi": f"output/{sample}/bam_sorted/{bam_id}.sorted.bam.csi"
        }
    elif bam_id in get_sra_ids(sample) or bam_id in get_fastq_ids(sample):
        return {
            "bam": f"output/{sample}/hisat2_aligned/{bam_id}.sorted.bam",
            "csi": f"output/{sample}/hisat2_aligned/{bam_id}.sorted.bam.csi"
        }
    elif bam_id in get_varus_ids(sample):
        return {
            "bam": f"output/{sample}/varus/varus.sorted.bam",
            "csi": f"output/{sample}/varus/varus.sorted.bam.csi"
        }
    else:
        raise ValueError(f"Unknown RNA-seq ID: {bam_id} for sample {sample}")

def get_protein_fasta_files(sample):
    """Get list of protein FASTA file paths for a sample.

    Supports colon-separated paths for multiple files.
    Returns empty list if no protein files specified.
    """
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    if pd.notna(row.get("protein_fasta")):
        return str(row["protein_fasta"]).split(":")
    return []

def get_protein_fasta(sample):
    """Get the protein FASTA path for pipeline consumption.

    If multiple files are specified (colon-separated), or any file is
    gzipped (.gz, as OrthoDB partitions are downloaded), returns the path
    to the merged and decompressed output of merge_proteins. A single
    uncompressed file is returned directly. Returns None if no proteins.
    """
    files = get_protein_fasta_files(sample)
    if len(files) > 1 or any(f.endswith(".gz") for f in files):
        return f"output/{sample}/preprocessing/proteins_merged.fa"
    elif len(files) == 1:
        return files[0]
    return None

def get_reference_gtf(sample):
    """Get reference annotation GTF path for a sample, or None."""
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    return row.get("reference_gtf") if pd.notna(row.get("reference_gtf")) else None

def get_busco_lineage(sample):
    """Get BUSCO lineage for a sample — must be set in samples.csv."""
    sample = sample.sample if hasattr(sample, 'sample') else sample
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    lineage = row.get("busco_lineage")
    if not lineage or (isinstance(lineage, float) and pd.isna(lineage)):
        raise ValueError(
            f"Sample '{sample}' is missing required 'busco_lineage' column in samples.csv"
        )
    return lineage

def get_output_dir(wildcards):
    """Get the output directory for a sample."""
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    return f"output/{sample}"

def get_braker_dir(wildcards):
    """Get the BRAKER directory (same as output directory in this workflow)."""
    return get_output_dir(wildcards)

def get_species_name(wildcards):
    """Get the species name for AUGUSTUS (derived from sample name)."""
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    # Replace special characters with underscores for AUGUSTUS compatibility
    return sample.replace("-", "_").replace(".", "_")

# ==============================================================================
# Mode-Dependent Path Routing
# ==============================================================================

def get_genemark_output(wildcards):
    """
    Get the GeneMark training genes GTF for AUGUSTUS training.

    This is the gene set used for flanking region calculation,
    GenBank conversion, and AUGUSTUS parameter training.

    - ES mode: GeneMark-ES ab initio predictions
    - ET mode: GeneMark-ET predictions (genemark.gtf)
    - EP mode: GeneMark-EP predictions (genemark.gtf)
    - ETP mode: GeneMark-ETP curated training genes (training.gtf)
    - IsoSeq mode: GeneMark-ETP training genes (same as ETP, GeneMark-ETP/training.gtf)
    - dual mode: training.gtf of the combined short-read + IsoSeq GeneMark-ETP
      run (GeneMark-ETP-isoseq/training.gtf)
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = get_effective_mode(sample)
    if mode == 'es':
        return f"output/{sample}/GeneMark-ES/genemark.gtf"
    elif mode == 'ep':
        return f"output/{sample}/GeneMark-EP/genemark.gtf"
    elif mode == 'dual':
        return f"output/{sample}/GeneMark-ETP-isoseq/training.gtf"
    elif mode in ('etp', 'isoseq'):
        return f"output/{sample}/GeneMark-ETP/training.gtf"
    else:  # et
        return f"output/{sample}/genemark/genemark.gtf"

def get_genemark_hints_for_filter(wildcards):
    """
    Get the hints file appropriate for filtering GeneMark predictions.

    In ES mode: no hints (empty file)
    In ET mode: RNA-seq intron hints (from bam2hints)
    In EP mode: protein intron hints (from ProtHint)
    In ETP mode: combined hints (from get_etp_hints.py)
    In IsoSeq mode: combined hints from GeneMark-ETP (same as ETP)
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = get_effective_mode(sample)
    if mode == 'es':
        return f"output/{sample}/empty_hints.gff"
    elif mode == 'ep':
        return f"output/{sample}/genemark_hints_ep.gff"
    elif mode in ('etp', 'isoseq'):
        return f"output/{sample}/etp_hints.gff"
    else:  # et
        return f"output/{sample}/bam2hints.gff"

def get_genemark_training_gtf(wildcards):
    """
    Get the training genes GTF for AUGUSTUS training.

    In ET/EP mode: filtered GeneMark genes (genemark.f.good.gtf)
    In ETP/IsoSeq mode: GeneMark-ETP training genes (training.gtf from ETP)
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = get_effective_mode(sample)
    if mode == 'dual':
        return f"output/{sample}/GeneMark-ETP-isoseq/training.gtf"
    elif mode in ('etp', 'isoseq'):
        return f"output/{sample}/GeneMark-ETP/training.gtf"
    else:
        return f"output/{sample}/genemark/genemark.f.good.gtf"

def get_augustus_gtf(wildcards):
    """
    Get the AUGUSTUS predictions GTF for downstream processing.

    In EP mode: uses iteration 2 output (after ProtHint refinement)
    In all other modes: uses iteration 1 output
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = get_effective_mode(sample)
    if mode == 'ep':
        return f"output/{sample}/augustus.hints_iter2.gtf"
    else:
        return f"output/{sample}/augustus.hints.gtf"

def get_augustus_hintsfile(wildcards):
    """
    Get the hintsfile used for the final AUGUSTUS prediction.

    In EP mode: iteration 2 hintsfile (after ProtHint refinement)
    In all other modes: standard hintsfile
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = get_effective_mode(sample)
    if mode == 'ep':
        return f"output/{sample}/hintsfile_iter2.gff"
    else:
        return f"output/{sample}/hintsfile.gff"

def get_extrinsic_cfg(wildcards):
    """
    Get the AUGUSTUS extrinsic configuration file for the sample's mode.

    - ET mode: rnaseq.cfg (RNA-Seq evidence only)
    - EP mode: ep.cfg (protein evidence only)
    - ETP/IsoSeq mode: etp.cfg (combined transcript + protein evidence)

    Used in params: the rule needs mode_decision_files among its inputs.
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = effective_mode_on_disk(sample)
    if mode in ('etp', 'isoseq', 'dual'):
        return "/opt/BRAKER/scripts/cfg/etp.cfg"
    elif mode == 'ep':
        return "/opt/BRAKER/scripts/cfg/ep.cfg"
    else:
        return "/opt/BRAKER/scripts/cfg/rnaseq.cfg"

# ==============================================================================
# Container Configuration
# ==============================================================================

BRAKER3_CONTAINER = config.get("braker3_image", "docker://teambraker/braker3:v3.1.3")
MINIMAP2_CONTAINER = config.get("minimap2_image", "docker://katharinahoff/minimap-minisplice:v0.1")
MINISPLICE_CONTAINER = config.get("minisplice_image", "docker://katharinahoff/minimap-minisplice:v0.1")
RED_CONTAINER = config.get("red_image", "docker://quay.io/biocontainers/red:2018.09.10--h9948957_3")
GFFCOMPARE_CONTAINER = config.get("gffcompare_image", "docker://quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_1")
AGAT_CONTAINER = config.get("agat_image", "docker://quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0")
BARRNAP_CONTAINER = config.get("pybarrnap_image", "docker://quay.io/biocontainers/pybarrnap:0.5.1--pyhdfd78af_0")
BUSCO_CONTAINER = config.get("busco_image", "docker://ezlabgva/busco:v6.1.0_cv2")
OMARK_CONTAINER = config.get("omark_image", "docker://quay.io/biocontainers/omark:0.4.1--pyh7e72e81_0")
TETOOLS_CONTAINER = config.get("tetools_image", "docker://dfam/tetools:latest")
VARUS_CONTAINER = config.get("varus_image", "docker://gaiusaugustus/pyvarus:v2.1.3")
TRNASCAN_CONTAINER = config.get("trnascan_image", "docker://quay.io/biocontainers/trnascan-se:2.0.12--pl5321h031d066_0")
FEELNC_CONTAINER = config.get("feelnc_image", "docker://quay.io/biocontainers/feelnc:0.2--pl526_0")
INFERNAL_CONTAINER = config.get("infernal_image", "docker://quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_2")

# GeneMark-ETP container: the separate braker3:isoseq image is no longer used;
# all GeneMark-ETP rules (etp, isoseq and dual mode) use BRAKER3_CONTAINER.
GENEMARK_ETP_CONTAINER = BRAKER3_CONTAINER

# ==============================================================================
# Workflow Summary
# ==============================================================================

def print_workflow_summary():
    """Print summary of detected samples and their configurations."""
    print("\n" + "="*70)
    print("BRAKER4 Snakemake Workflow")
    print("="*70)
    print(f"\nTotal samples: {len(SAMPLES)}")

    for sample in SAMPLES:
        mode = get_braker_mode(sample)
        types = detect_data_types(sample)
        print(f"\n  {sample}:")
        print(f"    Mode: {mode.upper()}")
        print(f"    Genome: {'pre-masked' if types['has_masked_genome'] else 'needs masking'}")

        evidence = []
        if types["has_bam"]:
            evidence.append(f"BAM ({len(get_bam_files(sample))} files)")
        if types["has_proteins"]:
            evidence.append("Proteins")
        if types["has_fastq"]:
            evidence.append("FASTQ")
        if types["has_sra"]:
            evidence.append("SRA")
        if types["has_varus"]:
            sources = [name for name, used in (("SRA", types["has_varus_sra"]),
                                               ("local FASTQ", types["has_varus_fastq"])) if used]
            evidence.append(f"VARUS ({' + '.join(sources)})")
        if types["has_isoseq"]:
            if types.get("has_isoseq_fastq"):
                evidence.append("IsoSeq (FASTQ/FASTA, needs alignment)")
            else:
                evidence.append("IsoSeq (pre-aligned BAM)")

        print(f"    Evidence: {', '.join(evidence)}")

    print("\n" + "="*70 + "\n")

# Print summary when workflow loads
print_workflow_summary()
