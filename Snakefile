"""
BRAKER3 Snakemake Workflow

Dynamic workflow supporting all BRAKER modes:
- ET: RNA-Seq only (Evidence from Transcripts)
- EP: Proteins only (Evidence from Proteins)
- ETP: RNA-Seq + Proteins (BRAKER3)
- IsoSeq: PacBio long-read + Proteins

Authors: Henning Krall, Katharina J. Hoff
Version: 0.3.0-beta
"""

__author__ = "Henning Krall,Katharina J. Hoff"
__version__ = "0.3.0-beta"

import pandas as pd
import configparser
import os

# ============================================================================
# Configuration
# ============================================================================

# Read config.ini
# The path can be overridden via the BRAKER4_CONFIG environment variable so
# multiple test scenarios can share a single config.ini.
config_ini_path = os.environ.get('BRAKER4_CONFIG', 'config.ini')
# Allow trailing # / ; comments on value lines (ConfigParser default is None).
config_parser = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
if not os.path.isfile(config_ini_path):
    raise FileNotFoundError(
        f"config.ini not found at '{config_ini_path}'. "
        "Set BRAKER4_CONFIG to a valid path or create config.ini in the "
        "working directory. See the README for an example configuration."
    )
config_parser.read(config_ini_path)

# Overlay environment-variable overrides on top of the parsed config so that
# users can change individual PARAMS or SLURM_ARGS values without editing
# config.ini. Variables follow the pattern BRAKER4_<SECTION>_<KEY>, except
# for the well-known SLURM_ARGS keys which use the short BRAKER4_<KEY> form.
# Anything set in the environment wins over the file.
for section in ('PARAMS', 'SLURM_ARGS', 'fantasia'):
    if not config_parser.has_section(section):
        config_parser.add_section(section)

_env_overrides = {
    # SLURM_ARGS
    'BRAKER4_CPUS_PER_TASK':                  ('SLURM_ARGS', 'cpus_per_task'),
    'BRAKER4_MEM_OF_NODE':                    ('SLURM_ARGS', 'mem_of_node'),
    'BRAKER4_MAX_RUNTIME':                    ('SLURM_ARGS', 'max_runtime'),
    # PARAMS
    'BRAKER4_FUNGUS':                         ('PARAMS', 'fungus'),
    'BRAKER4_MIN_CONTIG':                     ('PARAMS', 'min_contig'),
    'BRAKER4_SKIP_OPTIMIZE_AUGUSTUS':         ('PARAMS', 'skip_optimize_augustus'),
    'BRAKER4_USE_DEV_SHM':                    ('PARAMS', 'use_dev_shm'),
    'BRAKER4_GM_MAX_INTERGENIC':              ('PARAMS', 'gm_max_intergenic'),
    'BRAKER4_USE_COMPLEASM_HINTS':            ('PARAMS', 'use_compleasm_hints'),
    'BRAKER4_SKIP_BUSCO':                     ('PARAMS', 'skip_busco'),
    'BRAKER4_RUN_OMARK':                      ('PARAMS', 'run_omark'),
    'BRAKER4_NO_CLEANUP':                     ('PARAMS', 'no_cleanup'),
    'BRAKER4_TRANSLATION_TABLE':              ('PARAMS', 'translation_table'),
    'BRAKER4_GC_DONOR':                       ('PARAMS', 'gc_donor'),
    'BRAKER4_ALLOW_HINTED_SPLICESITES':       ('PARAMS', 'allow_hinted_splicesites'),
    'BRAKER4_RUN_NCRNA':                      ('PARAMS', 'run_ncrna'),
    'BRAKER4_RUN_BEST_BY_COMPLEASM':          ('PARAMS', 'run_best_by_compleasm'),
    'BRAKER4_MASKING_TOOL':                   ('PARAMS', 'masking_tool'),
    'BRAKER4_USE_MINISPLICE':                 ('PARAMS', 'use_minisplice'),
    'BRAKER4_USE_VARUS':                      ('PARAMS', 'use_varus'),
    'BRAKER4_SKIP_SINGLE_EXON_DOWNSAMPLING':  ('PARAMS', 'skip_single_exon_downsampling'),
    'BRAKER4_DOWNSAMPLING_LAMBDA':            ('PARAMS', 'downsampling_lambda'),
    'BRAKER4_AUGUSTUS_CHUNKSIZE':             ('PARAMS', 'augustus_chunksize'),
    'BRAKER4_AUGUSTUS_OVERLAP':               ('PARAMS', 'augustus_overlap'),
    # FANTASIA-Lite (optional, GPU-only functional annotation)
    'BRAKER4_RUN_FANTASIA':                   ('fantasia', 'enable'),
    'BRAKER4_FANTASIA_SIF':                   ('fantasia', 'sif'),
    'BRAKER4_FANTASIA_HF_CACHE':              ('fantasia', 'hf_cache_dir'),
    'BRAKER4_FANTASIA_PARTITION':             ('fantasia', 'partition'),
    'BRAKER4_FANTASIA_GPUS':                  ('fantasia', 'gpus'),
    'BRAKER4_FANTASIA_MEM_MB':                ('fantasia', 'mem_mb'),
    'BRAKER4_FANTASIA_CPUS':                  ('fantasia', 'cpus_per_task'),
    'BRAKER4_FANTASIA_MAX_RUNTIME':           ('fantasia', 'max_runtime'),
    'BRAKER4_FANTASIA_MIN_SCORE':             ('fantasia', 'min_score'),
    'BRAKER4_FANTASIA_ADDITIONAL_PARAMS':     ('fantasia', 'additional_params'),
}
for _env_name, (_section, _key) in _env_overrides.items():
    if _env_name in os.environ:
        config_parser.set(_section, _key, os.environ[_env_name])

augustus_config_path = os.path.abspath(config_parser['paths'].get('augustus_config_path', 'augustus_config'))

# Pass config to Snakemake config dict
config['samples_file'] = config_parser['paths'].get('samples_file', 'samples.csv')
config['fungus'] = config_parser.getboolean('PARAMS', 'fungus', fallback=False)
config['min_contig'] = config_parser.getint('PARAMS', 'min_contig', fallback=10000)
config['braker3_image'] = config_parser.get('containers', 'braker3_image',
                                             fallback='docker://teambraker/braker3:v3.0.10')

config['slurm_args'] = {
    'cpus_per_task': config_parser.getint('SLURM_ARGS', 'cpus_per_task', fallback=1),
    'mem_of_node': config_parser.getint('SLURM_ARGS', 'mem_of_node', fallback=16000),
    'max_runtime': config_parser.getint('SLURM_ARGS', 'max_runtime', fallback=60)
}

config['skip_optimize_augustus'] = config_parser.getboolean(
    'PARAMS',
    'skip_optimize_augustus',
    fallback=False
)

config['use_dev_shm'] = config_parser.getboolean(
    'PARAMS',
    'use_dev_shm',
    fallback=False
)

config['gm_max_intergenic'] = config_parser.getint(
    'PARAMS',
    'gm_max_intergenic',
    fallback=None
)

# Compleasm always runs in BRAKER4. There is intentionally no kill-switch
# for it because best_by_compleasm needs compleasm in protein mode to do
# BUSCO-rescue during the TSEBRA merge step. The compleasm rule also
# produces the genome-completeness QC summary used in the report.
#
# What IS configurable is whether the compleasm-derived CDSpart hints
# (from run_compleasm in genome mode) are fed into the AUGUSTUS hintsfile.
# When True (default), they are. When False, the same compleasm run still
# happens but its hints stay out of the AUGUSTUS hintsfile.
#
# Why this matters: BUSCO-derived c2h CDSpart hints can shift TSEBRA's
# per-transcript scoring on borderline AUGUSTUS predictions, occasionally
# costing a few percentage points of locus-level sensitivity. Disabling
# them produces a configuration closer to native braker.pl, which never
# uses BUSCO hints in the AUGUSTUS run.
config['use_compleasm_hints'] = config_parser.getboolean(
    'PARAMS', 'use_compleasm_hints', fallback=True
)

config['skip_busco'] = config_parser.getboolean(
    'PARAMS', 'skip_busco', fallback=False
)

config['run_omark'] = config_parser.getboolean(
    'PARAMS', 'run_omark', fallback=False
)

config['omamer_db'] = config_parser.get(
    'OMARK', 'omamer_db', fallback=False
)

config['no_cleanup'] = config_parser.getboolean(
    'PARAMS', 'no_cleanup', fallback=False
)

config['translation_table'] = config_parser.getint(
    'PARAMS', 'translation_table', fallback=1
)
# GeneMark-ES (gmes_petap.pl) inside the teambraker/braker3 container
# hard-rejects any genetic code except 1, 6, and 29. Verified at gmes_petap.pl
# sub CheckGcode — it exits with "error, input genetic code X is not supported"
# for anything else. Reject other codes at config-parse time so the failure
# is immediate with a clear explanation, instead of wasting compute on
# upstream rules that would eventually crash at GeneMark.
#
# Note: braker.pl itself nominally accepts tables 6, 10, 12, 25, 26, 27, 28,
# 29, 30, 31 at its top-level check, but the downstream gmes_petap.pl
# invocation would crash for any of 10/12/25/26/27/28/30/31 — so braker.pl's
# broader acceptance list is effectively a braker.pl bug we deliberately do
# not replicate here. See BRAKER_PARITY_AUDIT.md finding #15.
_SUPPORTED_GCODES = {1, 6, 29}
if config['translation_table'] not in _SUPPORTED_GCODES:
    raise ValueError(
        f"Invalid translation_table={config['translation_table']}. "
        f"BRAKER4 only supports {sorted(_SUPPORTED_GCODES)} because "
        "GeneMark-ES (gmes_petap.pl) inside the teambraker/braker3 "
        "container does not accept any other genetic code. "
        "Table 1 is the standard code; tables 6 and 29 are the "
        "ciliate/Mesodinium nuclear codes where TGA codes for Trp. "
        "See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"
    )

config['gc_donor'] = config_parser.getfloat(
    'PARAMS', 'gc_donor', fallback=0.001
)

config['allow_hinted_splicesites'] = config_parser.get(
    'PARAMS', 'allow_hinted_splicesites', fallback='gcag,atac'
)

# AUGUSTUS chunking for parallel hint-based prediction. createAugustusJoblist.pl
# splits the genome into chunks of `augustus_chunksize` bp with `augustus_overlap`
# bp shared between adjacent chunks for de-duplication by join_aug_pred.pl. If a
# gene is longer than the overlap, the two halves on adjacent chunks cannot be
# reconciled — so the overlap must exceed the longest expected gene.
#
# Defaults:
#   augustus_chunksize = 3000000  (BRAKER4 deviation; braker.pl uses 2500000)
#   augustus_overlap   = 500000   (matches braker.pl for hint-based predictions)
config['augustus_chunksize'] = config_parser.getint(
    'PARAMS', 'augustus_chunksize', fallback=3000000
)
config['augustus_overlap'] = config_parser.getint(
    'PARAMS', 'augustus_overlap', fallback=500000
)

config['run_ncrna'] = config_parser.getboolean(
    'PARAMS', 'run_ncrna', fallback=False
)

config['run_best_by_compleasm'] = config_parser.getboolean(
    'PARAMS', 'run_best_by_compleasm', fallback=True
)

config['masking_tool'] = config_parser.get(
    'PARAMS', 'masking_tool', fallback='repeatmasker'
).strip().lower()
if config['masking_tool'] not in ('repeatmasker', 'red'):
    raise ValueError(
        f"masking_tool must be 'repeatmasker' or 'red', got '{config['masking_tool']}'"
    )

config['use_minisplice'] = config_parser.getboolean(
    'PARAMS', 'use_minisplice', fallback=False
)

# FANTASIA-Lite functional annotation (optional, GPU-only).
# Off by default. The container hard-requires an NVIDIA GPU and has only been
# validated on an A100 in the Hoff lab; see README "run_fantasia" section.
config['run_fantasia'] = config_parser.getboolean(
    'fantasia', 'enable', fallback=False
)
config['fantasia'] = {
    'sif':              config_parser.get('fantasia', 'sif', fallback=''),
    'hf_cache_dir':     config_parser.get('fantasia', 'hf_cache_dir', fallback=''),
    'additional_params': config_parser.get('fantasia', 'additional_params', fallback=''),
    'min_score':        config_parser.get('fantasia', 'min_score', fallback='0.5'),
    'partition':        config_parser.get('fantasia', 'partition', fallback=''),
    'gpus':             config_parser.get('fantasia', 'gpus', fallback='1'),
    'mem_mb':           config_parser.get('fantasia', 'mem_mb',
                                          fallback=str(config['slurm_args']['mem_of_node'])),
    'cpus_per_task':    config_parser.get('fantasia', 'cpus_per_task',
                                          fallback=str(config['slurm_args']['cpus_per_task'])),
    'max_runtime':      config_parser.get('fantasia', 'max_runtime',
                                          fallback=str(config['slurm_args']['max_runtime'])),
}
if config['run_fantasia']:
    if not config['fantasia']['sif']:
        raise ValueError(
            "fantasia.enable=1 but fantasia.sif is empty. Set the path to the "
            "FANTASIA-Lite Singularity image (pre-pulled from "
            "docker://katharinahoff/fantasia_for_brain:lite.v0.0.2) in config.ini "
            "[fantasia] sif=, or via BRAKER4_FANTASIA_SIF. You can stage both "
            "the SIF and the ProtT5 cache by running:\n"
            "    BRAKER4_DOWNLOAD_FANTASIA=1 bash test_data/download_test_data.sh"
        )
    if not os.path.isfile(config['fantasia']['sif']):
        raise FileNotFoundError(
            f"fantasia.enable=1 but the configured SIF does not exist on disk: "
            f"{config['fantasia']['sif']}\n"
            "BRAKER4 will not start a run that would only fail later inside the "
            "FANTASIA rule. Either pre-pull the image yourself with "
            "`singularity pull <path> docker://katharinahoff/fantasia_for_brain:lite.v0.0.2`, "
            "or run:\n"
            "    BRAKER4_DOWNLOAD_FANTASIA=1 bash test_data/download_test_data.sh"
        )
    if not config['fantasia']['hf_cache_dir']:
        raise ValueError(
            "fantasia.enable=1 but fantasia.hf_cache_dir is empty. Pre-cache "
            "the Rostlab/prot_t5_xl_uniref50 HuggingFace model and set "
            "fantasia.hf_cache_dir to the cache directory in config.ini, or via "
            "BRAKER4_FANTASIA_HF_CACHE. The same downloader script "
            "(BRAKER4_DOWNLOAD_FANTASIA=1 bash test_data/download_test_data.sh) "
            "stages the cache for you."
        )
    if not os.path.isdir(config['fantasia']['hf_cache_dir']):
        raise FileNotFoundError(
            f"fantasia.enable=1 but the configured HuggingFace cache directory "
            f"does not exist: {config['fantasia']['hf_cache_dir']}\n"
            "FANTASIA-Lite runs in HF_HUB_OFFLINE mode, so the ProtT5 weights "
            "must already be on disk before the rule is invoked. Run:\n"
            "    BRAKER4_DOWNLOAD_FANTASIA=1 bash test_data/download_test_data.sh"
        )

config['username'] = os.environ.get("USER", "unknown")
script_dir = os.path.join(os.path.dirname(workflow.main_snakefile), "scripts")

# Shared data paths (avoid re-downloading per scenario)
_shared_data = os.path.join(os.path.dirname(workflow.main_snakefile), 'shared_data')
config['busco_download_path'] = os.path.abspath(
    config_parser.get('paths', 'busco_download_path',
                      fallback=os.path.join(_shared_data, 'busco_downloads'))
)
# compleasm uses a flat layout ({path}/{lineage}/) while BUSCO uses
# {busco_download_path}/lineages/{lineage}/. By defaulting
# compleasm_download_path to the BUSCO "lineages" subdirectory, a single
# extracted lineage tarball serves both tools and we avoid downloading the
# same data twice. Users can still override this in [paths] if they want
# the caches kept apart.
config['compleasm_download_path'] = os.path.abspath(
    config_parser.get('paths', 'compleasm_download_path',
                      fallback=os.path.join(config['busco_download_path'], 'lineages'))
)
config['rfam_path'] = os.path.abspath(
    config_parser.get('paths', 'rfam_path',
                      fallback=os.path.join(_shared_data, 'rfam'))
)

config['pipeline_version'] = __version__

# ============================================================================
# Include Common Functions and Sample Parsing
# ============================================================================

include: "rules/common.smk"

# ============================================================================
# Conditional Rule Inclusion Based on Data Types
# ============================================================================

# Detect which data types are present across all samples
GLOBAL_DATA_TYPES = {
    'needs_bam_sorting': False,
    'has_bam': False,
    'has_fastq': False,
    'has_sra': False,
    'has_varus': False,
    'has_isoseq': False,
    'has_isoseq_bam': False,
    'has_isoseq_fastq': False,
    'has_proteins': False,
    'needs_masking': False,
    'has_reference_gtf': False
}

for sample in SAMPLES:
    types = detect_data_types(sample)
    if types['has_bam']:
        GLOBAL_DATA_TYPES['needs_bam_sorting'] = True
        GLOBAL_DATA_TYPES['has_bam'] = True
    if types['has_fastq']:
        GLOBAL_DATA_TYPES['has_fastq'] = True
    if types['has_sra']:
        GLOBAL_DATA_TYPES['has_sra'] = True
    if types['has_varus']:
        GLOBAL_DATA_TYPES['has_varus'] = True
    if types['has_isoseq']:
        GLOBAL_DATA_TYPES['has_isoseq'] = True
        if types.get('has_isoseq_fastq'):
            GLOBAL_DATA_TYPES['has_isoseq_fastq'] = True
        else:
            GLOBAL_DATA_TYPES['has_isoseq_bam'] = True
    if types['has_proteins']:
        GLOBAL_DATA_TYPES['has_proteins'] = True
    if types['needs_masking']:
        GLOBAL_DATA_TYPES['needs_masking'] = True
    if types.get('has_reference_gtf'):
        GLOBAL_DATA_TYPES['has_reference_gtf'] = True

# Print what will be included
# Compute per-sample modes for accurate reporting
_MODES = {get_braker_mode(s) for s in SAMPLES}

print("\nIncluding rules based on detected data types:")
if GLOBAL_DATA_TYPES['has_sra']:
    print("  ✓ SRA download (prefetch + fastq-dump)")
if GLOBAL_DATA_TYPES['has_sra'] or GLOBAL_DATA_TYPES['has_fastq']:
    print("  ✓ HISAT2 alignment")
if GLOBAL_DATA_TYPES['needs_bam_sorting']:
    print("  ✓ BAM sorting and validation")
if 'et' in _MODES:
    print("  ✓ RNA-Seq hints generation (bam2hints)")
    print("  ✓ GeneMark-ET training")
if _MODES & {'etp', 'isoseq', 'dual'}:
    print("  ✓ GeneMark-ETP (combined transcript + protein training)")
if 'dual' in _MODES:
    print("  ✓ Dual-ETP mode (separate short-read + IsoSeq GeneMark-ETP runs)")
if 'ep' in _MODES:
    print("  ✓ Protein hints generation (ProtHint)")
    print("  ✓ GeneMark-ES (ab initio seeds)")
    print("  ✓ GeneMark-EP (protein evidence training)")
if 'es' in _MODES:
    print("  ✓ GeneMark-ES (ab initio, no evidence)")
if GLOBAL_DATA_TYPES['has_varus']:
    print("  ✓ VARUS auto-download and alignment")
if GLOBAL_DATA_TYPES['has_isoseq_fastq']:
    print("  ✓ minimap2 IsoSeq FASTQ/FASTA alignment")
if GLOBAL_DATA_TYPES['needs_masking']:
    if config['masking_tool'] == 'red':
        print("  ✓ Red repeat masking (fast, no repeat library)")
    else:
        print("  ✓ RepeatModeler2 + RepeatMasker (genome masking)")
if GLOBAL_DATA_TYPES['has_reference_gtf']:
    print("  ✓ gffcompare evaluation against reference annotation")
print()

# Include SRA download rules
if GLOBAL_DATA_TYPES['has_sra']:
    include: "rules/preprocessing/download_sra.smk"

# Include HISAT2 alignment rules (for SRA downloads and user-provided FASTQs)
if GLOBAL_DATA_TYPES['has_sra'] or GLOBAL_DATA_TYPES['has_fastq']:
    include: "rules/preprocessing/hisat2_align.smk"

# Include minimap2 alignment rules (for unaligned IsoSeq FASTA/FASTQ)
if GLOBAL_DATA_TYPES['has_isoseq_fastq']:
    include: "rules/preprocessing/minimap2_isoseq_align.smk"
    # Optional minisplice splice-site scoring (improves minimap2 junction detection)
    if config.get('use_minisplice', False):
        print("  ✓ minisplice splice-site scoring for IsoSeq alignment")
        include: "rules/preprocessing/run_minisplice.smk"

# Include VARUS auto-download rules
if GLOBAL_DATA_TYPES['has_varus']:
    include: "rules/preprocessing/run_varus.smk"

# Include BAM sorting rules (only for pre-existing BAM files)
if GLOBAL_DATA_TYPES['has_bam']:
    include: "rules/genemark/check_bam_sorted.smk"

# Include IsoSeq BAM sorting/merging rules (for pre-aligned BAM or minimap2 output)
if GLOBAL_DATA_TYPES['has_isoseq_bam'] or GLOBAL_DATA_TYPES['has_isoseq_fastq']:
    include: "rules/preprocessing/check_isoseq_bam.smk"

# Include RNA-Seq hints and GeneMark-ET rules (for short-read RNA-Seq evidence)
if GLOBAL_DATA_TYPES['has_bam'] or GLOBAL_DATA_TYPES['has_varus'] or GLOBAL_DATA_TYPES['has_sra'] or GLOBAL_DATA_TYPES['has_fastq']:
    include: "rules/genemark/bam2hints.smk"
    include: "rules/genemark/join_hints.smk"
    include: "rules/genemark/filter_introns.smk"
    include: "rules/genemark/run_genemark_et.smk"

# Include GeneMark-ES (needed for EP mode seeding and ES mode)
HAS_ES = any(get_braker_mode(s) == 'es' for s in SAMPLES)
if GLOBAL_DATA_TYPES['has_proteins'] or HAS_ES:
    include: "rules/genemark/run_genemark_es.smk"

# Include EP mode rules (protein-only evidence)
if GLOBAL_DATA_TYPES['has_proteins']:
    include: "rules/genemark/run_prothint.smk"
    include: "rules/genemark/prepare_genemark_hints_ep.smk"
    include: "rules/genemark/run_genemark_ep.smk"
    include: "rules/genemark/run_prothint_iter2.smk"

# Include ETP mode rules (RNA-Seq + proteins combined, also used for IsoSeq and dual mode)
HAS_ETP = any(get_braker_mode(s) in ('etp', 'isoseq', 'dual') for s in SAMPLES)
if HAS_ETP:
    include: "rules/genemark/run_genemark_etp.smk"

# Include dual-ETP mode rules (IsoSeq + short-read RNA-Seq + proteins)
HAS_DUAL = any(get_braker_mode(s) == 'dual' for s in SAMPLES)
if HAS_DUAL:
    include: "rules/genemark/run_genemark_etp_isoseq.smk"
    include: "rules/genemark/merge_dual_etp.smk"

# GeneMark filtering (always needed, uses mode-dependent input functions)
include: "rules/genemark/filter_genemark.smk"

# Include AUGUSTUS training rules (always needed)
include: "rules/augustus_training/copy_augustus_config.smk"
include: "rules/augustus_training/create_new_species.smk"
include: "rules/augustus_training/convert_to_genbank.smk"
include: "rules/augustus_training/split_training_set.smk"
include: "rules/augustus_training/train_augustus.smk"
include: "rules/augustus_training/optimize_augustus.smk"

include: "rules/postprocessing/redundancy_removal.smk"
include: "rules/augustus_training/filter_genes_etraining.smk"
include: "rules/augustus_training/downsample_training_genes.smk"
include: "rules/augustus_training/compute_flanking_region.smk"

# Include AUGUSTUS prediction rules
include: "rules/postprocessing/merge_hints.smk"
include: "rules/augustus_predict/run_augustus_hints.smk"
include: "rules/postprocessing/fix_in_frame_stop_codons.smk"

# Include AUGUSTUS iteration 2 for EP mode (ProtHint refinement)
HAS_EP = any(get_braker_mode(s) == 'ep' for s in SAMPLES)
if HAS_EP:
    include: "rules/augustus_predict/run_augustus_hints_iter2.smk"

# TSEBRA merging (always included, mode-aware per sample)
include: "rules/postprocessing/run_tsebra.smk"

# best_by_compleasm: BUSCO-driven gene rescue step (enabled by default).
# Always included so the DAG can flow braker.tsebra.raw.gtf -> braker.tsebra.gtf;
# when run_best_by_compleasm = 0 the rule just copies the file through.
include: "rules/postprocessing/best_by_compleasm.smk"

# HC gene extraction (for ETP and IsoSeq mode samples)
if HAS_ETP:
    include: "rules/postprocessing/extract_hc_training_genes.smk"

# Rescue multi-exon genes where both AUGUSTUS and GeneMark agree on the CDS
# intron chain but TSEBRA dropped them (all introns must be hint-supported).
include: "rules/postprocessing/both_agree_rescue.smk"

include: "rules/postprocessing/filter_stop_codons.smk"
include: "rules/postprocessing/normalize_cds.smk"
include: "rules/postprocessing/postprocess_augustus.smk"

# Include UTR decoration rules (for any mode with transcript evidence)
HAS_TRANSCRIPTS = GLOBAL_DATA_TYPES['has_bam'] or GLOBAL_DATA_TYPES['has_fastq'] or GLOBAL_DATA_TYPES['has_sra'] or GLOBAL_DATA_TYPES['has_varus'] or GLOBAL_DATA_TYPES['has_isoseq']
if HAS_TRANSCRIPTS:
    include: "rules/postprocessing/add_utr.smk"

# Include GFF3 conversion
include: "rules/postprocessing/convert_to_gff3.smk"

# Include QC rules
include: "rules/quality_control/run_compleasm.smk"
include: "rules/quality_control/gene_support.smk"

# BUSCO is the slow QC step (full pipeline including HMMER searches against the
# lineage dataset). Skip it via skip_busco = 1 in config.ini for fast iteration.
if not config.get('skip_busco', False):
    include: "rules/quality_control/run_busco.smk"

# Include OMArk rules (optional, requires LUCA.h5 database)
RUN_OMARK = config_parser.getboolean('PARAMS', 'run_omark', fallback=False)
if RUN_OMARK:
    include: "rules/quality_control/run_omark.smk"

# Include evaluation rules (optional, requires reference annotation)
if GLOBAL_DATA_TYPES['has_reference_gtf']:
    include: "rules/quality_control/run_gffcompare.smk"

# Include preprocessing rules
include: "rules/preprocessing/prepare_genome.smk"

# Include protein merge rule (for multiple protein files)
if GLOBAL_DATA_TYPES['has_proteins']:
    include: "rules/preprocessing/merge_proteins.smk"

# Include masking rules (conditional)
if GLOBAL_DATA_TYPES['needs_masking']:
    if config['masking_tool'] == 'red':
        include: "rules/preprocessing/run_red_masking.smk"
    else:
        include: "rules/preprocessing/run_masking.smk"

# Include ncRNA rules (optional). When run_ncrna=1 the pipeline annotates:
#   - rRNA via barrnap (always part of run_ncrna)
#   - tRNA via tRNAscan-SE
#   - snoRNA, snRNA, miRNA, ribozymes etc. via Infernal/cmscan against Rfam
#   - lncRNA via FEELnc (only when transcript evidence is available)
# All ncRNA predictions are merged into braker_with_ncRNA.gff3.
RUN_NCRNA = config.get('run_ncrna', False)
if RUN_NCRNA:
    include: "rules/ncrna/run_barrnap.smk"
    include: "rules/ncrna/run_trnascan.smk"
    include: "rules/ncrna/run_infernal.smk"
    if HAS_TRANSCRIPTS:
        include: "rules/ncrna/run_feelnc.smk"

# Include FANTASIA-Lite functional annotation rules (optional, GPU-only).
# This is the most fragile component of BRAKER4 -- the upstream container
# requires an NVIDIA GPU (validated on A100). See README "run_fantasia" section.
if config.get('run_fantasia', False):
    print("  ✓ FANTASIA-Lite functional annotation (GPU)")
    include: "rules/postprocessing/run_fantasia.smk"

# Include results collection rule
include: "rules/postprocessing/collect_results.smk"

# ============================================================================
# Target Rule
# ============================================================================

def get_final_outputs():
    """
    Determine final outputs based on sample configurations.

    The collect_results rule gathers all important output files into
    output/{sample}/results/ and removes intermediate files. Its inputs
    are determined dynamically by _get_collect_inputs() in collect_results.smk.
    """
    outputs = []
    for sample in SAMPLES:
        outputs.append(f"output/{sample}/results/.done")
    return outputs

rule all:
    input:
        get_final_outputs()
    default_target: True

# ============================================================================
# Wildcard Constraints
# ============================================================================

wildcard_constraints:
    sample="|".join(SAMPLES)

# ============================================================================
# Workflow Summary
# ============================================================================

print("\nTarget outputs:")
for output in get_final_outputs():
    print(f"  → {output}")
print(f"\n{'='*70}\n")
