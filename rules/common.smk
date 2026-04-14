"""
Common functions and configuration for BRAKER3 Snakemake workflow.

This module provides:
- Sample parsing from samples.csv
- Data type detection (BAM, FASTQ, proteins, etc.)
- Input validation
- Helper functions for dynamic input resolution
"""

import pandas as pd
import os
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

    return {
        "has_proteins": pd.notna(row.get("protein_fasta")),
        "has_bam": pd.notna(row.get("bam_files")),
        "has_fastq": pd.notna(row.get("fastq_r1")) and pd.notna(row.get("fastq_r2")),
        "has_sra": pd.notna(row.get("sra_ids")),
        "has_varus": pd.notna(row.get("varus_genus")) and pd.notna(row.get("varus_species")),
        "has_isoseq": pd.notna(row.get("isoseq_bam")) or pd.notna(row.get("isoseq_fastq")),
        "has_isoseq_fastq": pd.notna(row.get("isoseq_fastq")),
        "has_masked_genome": pd.notna(row.get("genome_masked")),
        "needs_masking": pd.isna(row.get("genome_masked")),
        "has_reference_gtf": pd.notna(row.get("reference_gtf"))
    }

def get_braker_mode(sample):
    """
    Determine BRAKER mode for a sample.

    - ET: RNA-Seq only (no proteins)
    - EP: Proteins only (no RNA-Seq)
    - ETP: RNA-Seq + Proteins (standard BRAKER3)
    - IsoSeq: PacBio long-read + Proteins (no short-read RNA-Seq)
    - dual: Short-read RNA-Seq + IsoSeq + Proteins (two ETP runs, then merge)
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
# Input Validation
# ==============================================================================

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

        # Rule 7: busco_lineage is required
        if pd.isna(row.get("busco_lineage")) or not str(row.get("busco_lineage")).strip():
            raise ValueError(f"{sample}: busco_lineage column is required in samples.csv")

        # Rule 8: Alternative translation tables only supported in ES/ET/EP modes
        tt = config.get('translation_table', 1)
        if tt != 1:
            mode = get_braker_mode(sample)
            if mode in ('etp', 'isoseq', 'dual'):
                raise ValueError(
                    f"{sample}: translation_table={tt} is not supported in {mode.upper()} mode. "
                    f"GeneMark-ETP does not support alternative genetic codes. "
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
    """
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

def get_varus_ids(sample):
    """
    Get VARUS RNA-seq ID for a sample.

    Returns ["varus"] if VARUS is configured, empty list otherwise.
    """
    types = detect_data_types(sample)
    if types["has_varus"]:
        return ["varus"]
    return []

def get_isoseq_ids(sample):
    """
    Get IsoSeq BAM ID for a sample.

    Returns ["isoseq"] if IsoSeq data is present, empty list otherwise.
    IsoSeq BAMs provide transcript evidence like RNA-Seq and flow through
    bam2hints → filter_introns → GeneMark-ET.
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
            "bai": f"output/{sample}/bam_sorted/{bam_id}.sorted.bam.bai"
        }
    elif bam_id in get_sra_ids(sample) or bam_id in get_fastq_ids(sample):
        return {
            "bam": f"output/{sample}/hisat2_aligned/{bam_id}.sorted.bam",
            "bai": f"output/{sample}/hisat2_aligned/{bam_id}.sorted.bam.bai"
        }
    elif bam_id in get_varus_ids(sample):
        return {
            "bam": f"output/{sample}/varus/varus.sorted.bam",
            "bai": f"output/{sample}/varus/varus.sorted.bam.bai"
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

    If multiple files are specified (colon-separated), returns the path
    to the merged output. If single file, returns it directly.
    Returns None if no proteins.
    """
    files = get_protein_fasta_files(sample)
    if len(files) > 1:
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
    - IsoSeq mode: GeneMark-ETP training genes (same as ETP, with isoseq container)
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = get_braker_mode(sample)
    if mode == 'es':
        return f"output/{sample}/GeneMark-ES/genemark.gtf"
    elif mode == 'ep':
        return f"output/{sample}/GeneMark-EP/genemark.gtf"
    elif mode == 'dual':
        return f"output/{sample}/dual_etp_merged/training.gtf"
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
    mode = get_braker_mode(sample)
    if mode == 'es':
        return f"output/{sample}/empty_hints.gff"
    elif mode == 'ep':
        return f"output/{sample}/genemark_hints_ep.gff"
    elif mode == 'dual':
        return f"output/{sample}/dual_etp_merged/etp_hints.gff"
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
    mode = get_braker_mode(sample)
    if mode == 'dual':
        return f"output/{sample}/dual_etp_merged/training.gtf"
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
    mode = get_braker_mode(sample)
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
    mode = get_braker_mode(sample)
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
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    mode = get_braker_mode(sample)
    if mode in ('etp', 'isoseq', 'dual'):
        return "/opt/BRAKER/scripts/cfg/etp.cfg"
    elif mode == 'ep':
        return "/opt/BRAKER/scripts/cfg/ep.cfg"
    else:
        return "/opt/BRAKER/scripts/cfg/rnaseq.cfg"

# ==============================================================================
# Container Configuration
# ==============================================================================

BRAKER3_CONTAINER = config.get("braker3_image", "docker://teambraker/braker3:v3.0.10")
ISOSEQ_CONTAINER = config.get("isoseq_image", "docker://teambraker/braker3:isoseq")
MINIMAP2_CONTAINER = config.get("minimap2_image", "docker://katharinahoff/minimap-minisplice:v0.1")
MINISPLICE_CONTAINER = config.get("minisplice_image", "docker://katharinahoff/minimap-minisplice:v0.1")
GFFCOMPARE_CONTAINER = config.get("gffcompare_image", "docker://quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_1")
AGAT_CONTAINER = config.get("agat_image", "docker://quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0")
BARRNAP_CONTAINER = config.get("barrnap_image", "docker://quay.io/biocontainers/barrnap:0.9--hdfd78af_4")
BUSCO_CONTAINER = config.get("busco_image", "docker://ezlabgva/busco:v6.0.0_cv1")
OMARK_CONTAINER = config.get("omark_image", "docker://quay.io/biocontainers/omark:0.4.1--pyh7e72e81_0")

# GeneMark-ETP container: uses IsoSeq container when any sample is in IsoSeq mode
# (the IsoSeq container has a GeneMark-ETP build that handles long-read evidence)
# All other rules always use BRAKER3_CONTAINER
# GeneMark-ETP container selection:
# - 'isoseq' mode (IsoSeq only): use ISOSEQ_CONTAINER for the single ETP run
# - 'dual' mode: short-read ETP uses BRAKER3_CONTAINER, IsoSeq ETP uses ISOSEQ_CONTAINER
# - 'etp' mode: use BRAKER3_CONTAINER
HAS_ISOSEQ_ONLY = any(get_braker_mode(s) == 'isoseq' for s in SAMPLES)
GENEMARK_ETP_CONTAINER = ISOSEQ_CONTAINER if HAS_ISOSEQ_ONLY else BRAKER3_CONTAINER

# ==============================================================================
# Workflow Summary
# ==============================================================================

def print_workflow_summary():
    """Print summary of detected samples and their configurations."""
    print("\n" + "="*70)
    print("BRAKER3 Snakemake Workflow")
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
            evidence.append("VARUS")
        if types["has_isoseq"]:
            if types.get("has_isoseq_fastq"):
                evidence.append("IsoSeq (FASTQ/FASTA, needs alignment)")
            else:
                evidence.append("IsoSeq (pre-aligned BAM)")

        print(f"    Evidence: {', '.join(evidence)}")

    print("\n" + "="*70 + "\n")

# Print summary when workflow loads
print_workflow_summary()
