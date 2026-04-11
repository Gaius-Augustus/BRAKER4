# BRAKER4 Pipeline - Test Scenarios

Test scenarios covering all unique pipeline paths for the BRAKER4 Snakemake pipeline.

**HPC scenarios** (`test_scenarios/`, 08-11): 4 scenarios on *O. tauri* genome, run on HPC with SLURM.
**Local scenarios** (`test_scenarios_local/`, 01-07): 7 scenarios on small *A. thaliana* fragment, run locally with 8 cores.

Together they cover every input type, every BRAKER mode, and masking — without redundancy.

## Quick Start

```bash
# Run all HPC scenarios (two waves: fast first, then slow)
cd test_scenarios
bash run_all_tests.sh

# Run specific scenarios by number
bash run_all_tests.sh 08 11

# Dry-run only (validate DAGs without executing)
DRY_RUN=true bash run_all_tests.sh

# Run a single scenario
cd test_scenarios/scenario_11_etp_sra
bash run_test.sh
```

For local testing (no SLURM, 8 cores, cheap scenarios only), see `test_scenarios_local/`.

## Test Scenarios

### HPC scenarios (`test_scenarios/`)

Each covers a unique combination of input source and BRAKER mode on the *O. tauri* genome (~12.6 MB):

| #  | Name          | Mode | Evidence       | Masking       | Unique steps tested                      |
|----|---------------|------|----------------|---------------|------------------------------------------|
| 08 | es            | ES   | (none)         | Pre-masked    | GeneMark-ES (ab initio)                  |
| 09 | ep\_masking   | EP   | Proteins       | RepeatMasker  | Masking + ProtHint + GeneMark-EP + iter2 |
| 10 | et\_varus     | ET   | VARUS          | Pre-masked    | VARUS + GeneMark-ET                      |
| 11 | etp\_sra      | ETP  | SRA + proteins | Pre-masked    | SRA download + HISAT2 + GeneMark-ETP     |

### Local scenarios (`test_scenarios_local/`)

Fast scenarios on the small *A. thaliana* chr5 fragment (~1 MB, pre-masked):

| #  | Name            | Mode   | Evidence                     | Unique steps tested                           |
|----|-----------------|--------|------------------------------|-----------------------------------------------|
| 01 | es              | ES     | (none)                       | GeneMark-ES (ab initio)                       |
| 02 | ep              | EP     | Proteins                     | ProtHint + GeneMark-EP + iter2                |
| 03 | et\_fastq       | ET     | FASTQ                        | HISAT2 + GeneMark-ET                          |
| 04 | etp\_bam        | ETP    | BAM + proteins               | BAM check + GeneMark-ETP                      |
| 05 | isoseq\_bam     | IsoSeq | IsoSeq BAM (multi) + prot.  | IsoSeq BAM check + merge + IsoSeq ETP        |
| 06 | isoseq\_fasta   | IsoSeq | IsoSeq FASTA + prot.        | minimap2 + IsoSeq ETP                         |
| 07 | dual            | Dual   | BAM + IsoSeq BAM + prot.    | BAM + IsoSeq + dual ETP + merge               |

### Naming Convention

`scenario_<number>_<mode>[_<source>]`

- `es` = ES mode (ab initio, no evidence)
- `ep` = EP mode (proteins only)
- `et_<source>` = ET mode (RNA-Seq only)
- `etp_<source>` = ETP mode (RNA-Seq + proteins)
- `isoseq_<type>` = IsoSeq mode
- `dual` = Dual mode (short-read + IsoSeq + proteins)

## Test Data

**Included in git** (for local scenarios, *A. thaliana* chr5 fragment):

- `test_data/genome.fa` (993 KB) -- small A. thaliana chr5 fragment (8 contigs)
- `test_data/proteins.fa` (2.6 MB) -- OrthoDB protein subset
- `test_data/isoseq.bam` -- PacBio IsoSeq reads aligned to test genome
- `test_data/isoseq.fastq.gz` -- unaligned IsoSeq reads (for minimap2 scenarios)
- `test_data/reference.gtf` -- reference annotation for gffcompare

**Downloaded on demand** (`bash test_data/download_test_data.sh`):

- `test_data/Ostreococcus_tauri.fa` (13 MB) -- *O. tauri* genome (for HPC scenarios)
- `test_data/Viridiplantae.fa` (5.3 GB) -- ODB12 proteins (for protein scenarios)
- `test_data/RNAseq.bam` (164 MB) -- *A. thaliana* RNA-Seq alignments (local scenarios)
- `test_data/reads_R1.fastq.gz` (77 MB), `reads_R2.fastq.gz` (78 MB) -- paired FASTQ (from RNAseq.bam)
- `test_data/LUCA.h5` (8.8 GB) -- OMAmer database for OMArk (optional, only needed if `run_omark=1`)
- `shared_data/rfam/` (650 MB) -- Rfam covariance models (for `run_ncrna=1`)

## How the Test Runner Works

The `run_all_tests.sh` script launches each scenario's Snakemake process under `nohup`, so it survives shell disconnects and screen detachments. Each Snakemake process then submits its own SLURM jobs.

All 4 HPC scenarios are launched concurrently. Scenarios 08 and 11 (pre-masked) finish faster; 09 (masking) and 10 (VARUS) take longer.

## Monitoring Running Tests

```bash
# Check status of all running test processes
bash run_all_tests.sh status
```

Output shows PID, scenario name, and whether each is RUNNING, DONE, or DEAD:
```
=== Running BRAKER test processes ===
  RUNNING  PID=12345  scenario_11_etp_sra  (last: 15 steps)
  RUNNING  PID=12346  scenario_09_ep_masking  (last: 8 steps)
  DONE     PID=12347  scenario_08_es
```

### Checking individual scenario logs

```bash
# Tail the log of a specific scenario
tail -f test_results_*/scenario_11_etp_sra.log

# Check how many steps completed
grep "steps.*done" test_results_*/scenario_11_etp_sra.log | tail -1
```

### Checking SLURM jobs

```bash
squeue -u $USER
```

## Killing Tests

```bash
# Kill all running test Snakemake processes
bash run_all_tests.sh kill
```

This kills the Snakemake controller processes. SLURM jobs that were already submitted may still be running. To cancel those too:

```bash
scancel -u $USER
```

## Running a Single Scenario

```bash
cd test_scenarios/scenario_11_etp_sra
bash run_test.sh
```

With options:

```bash
DRY_RUN=true bash run_test.sh          # Dry-run
PARTITION=snowball bash run_test.sh     # Override partition list ad-hoc
USE_SLURM=false bash run_test.sh       # Run locally (no SLURM)
CORES=8 bash run_test.sh               # Limit cores (local mode)
```

## Configuration

All scenarios in this directory share **two** configuration files:

| File | What it controls |
|------|------------------|
| `compute_profile.sh` | Cluster / compute knobs: cores, partition list, memory, per-job runtime, SLURM on/off |
| `config.ini` | Biology knobs shared by every scenario: optimisation flags, masking, ncRNA, etc. |

Edit `compute_profile.sh` once to point the test suite at your own cluster — no per-scenario edits needed. The same applies to `test_scenarios_local/compute_profile.sh` for the local (no-SLURM) tier.

### compute_profile.sh

The variables defined in `compute_profile.sh` use the `${VAR:-default}` idiom, so any of them can still be overridden ad-hoc by exporting an environment variable before invoking a scenario:

```bash
PARTITION=highmem CORES=64 bash scenario_08_es/run_test.sh
BRAKER4_MAX_RUNTIME=240 bash scenario_11_etp_sra/run_test.sh
```

Knobs you'll typically tune:

| Variable | Default (HPC) | Default (local) | Effect |
|----------|---------------|-----------------|--------|
| `CORES` | `48` | `8` | Cores snakemake itself uses (`--cores`/`--jobs`) |
| `PARTITION` | `batch,snowball,pinky` | `batch` | SLURM partition list |
| `USE_SLURM` | `true` | `false` | Submit each rule as a SLURM job |
| `BRAKER4_CPUS_PER_TASK` | `48` | `8` | CPUs requested per SLURM job |
| `BRAKER4_MEM_OF_NODE` | `120000` | `32000` | MB requested per SLURM job |
| `BRAKER4_MAX_RUNTIME` | `120` | `1440` | Per-job wall-time limit (minutes) |
| `DEFAULT_MEM_MB` | `120000` | `32000` | Snakemake `--default-resources mem_mb` cap |

### config.ini (shared biology defaults)

The shared `config.ini` carries the biology knobs that every scenario uses. The Snakefile finds it via the `BRAKER4_CONFIG` environment variable, which `run_test.sh` sets to `$TIER_DIR/config.ini`. Example:

```ini
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config

[PARAMS]
run_best_by_compleasm = 1
run_omark = 0
skip_optimize_augustus = 1
skip_single_exon_downsampling = 0
downsampling_lambda = 2
use_varus = 0
use_compleasm_hints = 1
use_dev_shm = 0
gm_max_intergenic = 10000
run_ncrna = 0
```

Note that `[SLURM_ARGS]` is no longer in the shared config — those values are exported by `compute_profile.sh` as `BRAKER4_*` environment variables and the Snakefile picks them up automatically.

### Per-scenario overrides (scenario_overrides.sh)

A handful of scenarios deviate from the shared defaults — for example, `scenario_08_es` is the OMArk QC test and needs `run_omark=1`, while `scenario_09_ep_masking` and `scenario_10_et_varus` need a longer per-job runtime than the 120-minute shared default. Those deviations live in an optional `scenario_overrides.sh` file inside the scenario directory:

```bash
# scenario_08_es/scenario_overrides.sh
export BRAKER4_RUN_OMARK=1
```

`run_test.sh` sources `compute_profile.sh` first, then `scenario_overrides.sh` (if present), so anything exported in the override file wins. Most scenarios do not need an override file at all.

### Environment-variable override reference

Any value in `[SLURM_ARGS]` or `[PARAMS]` can be overridden via a `BRAKER4_<KEY_UPPER>` environment variable. The full list of recognised names lives in `Snakefile` (search for `_env_overrides`). Common ones:

| Env var | Section / key |
|---------|---------------|
| `BRAKER4_CPUS_PER_TASK` | `[SLURM_ARGS] cpus_per_task` |
| `BRAKER4_MEM_OF_NODE` | `[SLURM_ARGS] mem_of_node` |
| `BRAKER4_MAX_RUNTIME` | `[SLURM_ARGS] max_runtime` |
| `BRAKER4_RUN_OMARK` | `[PARAMS] run_omark` |
| `BRAKER4_RUN_NCRNA` | `[PARAMS] run_ncrna` |
| `BRAKER4_USE_VARUS` | `[PARAMS] use_varus` |
| `BRAKER4_SKIP_OPTIMIZE_AUGUSTUS` | `[PARAMS] skip_optimize_augustus` |
| `BRAKER4_NO_CLEANUP` | `[PARAMS] no_cleanup` |
| `BRAKER4_SKIP_BUSCO` | `[PARAMS] skip_busco` |

### Environment variables for run\_all\_tests.sh

| Variable | Default | Description |
|----------|---------|-------------|
| `DRY_RUN` | `false` | Set to `true` for DAG validation only |
| `MAX_CONCURRENT` | (all) | Maximum concurrent Snakemake processes |

## Output

Each scenario produces final outputs in:

```
scenario_XX_name/output/<sample_name>/results/
├── braker.gtf                              # Gene predictions
├── braker.gff3                             # GFF3 format
├── braker.aa                               # Protein sequences
├── braker.codingseq                        # CDS sequences
├── quality_control/
│   ├── training_summary.txt                # Training gene counts + accuracy
│   ├── training_summary.pdf                # Publication-quality training plot
│   ├── busco_summary.txt                   # BUSCO completeness
│   ├── compleasm_summary.txt               # Compleasm completeness
│   └── omark_summary.txt                   # OMArk quality (if enabled)
└── gene_support.tsv                        # Per-gene evidence support
```

All intermediate files are removed automatically after results collection.

## Troubleshooting

### "IncompleteFilesException"

A previous run crashed mid-step. All `run_test.sh` scripts include `--rerun-incomplete` which handles this automatically. If it persists:

```bash
cd scenario_XX_name
snakemake --snakefile ../../Snakefile --unlock
rm -rf .snakemake/incomplete/
bash run_test.sh
```

### "LockException"

Another Snakemake process is (or was) running in the same directory:

```bash
cd scenario_XX_name
snakemake --snakefile ../../Snakefile --unlock
bash run_test.sh
```

### Snakemake dies but SLURM jobs keep running

The `nohup` wrapper prevents most of these. If the Snakemake process dies anyway, the SLURM jobs become orphans. Cancel them and re-run:

```bash
scancel -u $USER
bash run_test.sh   # Snakemake picks up from where it left off
```

### VARUS scenarios fail at GeneMark-ETP

VARUS needs enough data. With `--maxBatches 600` and `--batchSize 75000`, VARUS downloads up to 45M reads. If this is still insufficient (mapping rate < 5%), the genome may need proper repeat masking first, or VARUS may not find suitable RNA-Seq data for this species.

### Container python issues

All `run_test.sh` scripts export `SINGULARITYENV_PREPEND_PATH=/opt/conda/bin` and shell blocks use `export PATH=/opt/conda/bin:$PATH` to ensure the container's python (with biopython, intervaltree, matplotlib) is used instead of the host's. If you see `ModuleNotFoundError`, check that the v3.0.8 container is being used.

## Rule Graphs

Each scenario contains a `rulegraph.png` showing the execution DAG. To regenerate:

```bash
bash test_scenarios/generate_rulegraphs.sh
```

Requires graphviz. Run on a workstation, not on HPC.
