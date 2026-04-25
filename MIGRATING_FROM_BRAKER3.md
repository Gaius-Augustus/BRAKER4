# Migrating from BRAKER3 (`braker.pl`) to BRAKER4

This is a step-by-step guide for users who already know how to run the original BRAKER3 Perl pipeline (`braker.pl`) and want to switch to BRAKER4 (the Snakemake rewrite). It assumes you can already run a successful `braker.pl` annotation and that you have a working Singularity (or Apptainer) installation.

If you have never used BRAKER before, read the main [README.md](README.md) first — this guide will not re-explain ET / EP / ETP modes or what BRAKER does conceptually.

## What changes (and what does not)

| Concept | BRAKER3 (`braker.pl`) | BRAKER4 |
|---|---|---|
| Main entry point | `braker.pl --genome=... --bam=...` | `snakemake --snakefile path/to/Snakefile ...` |
| Input specification | Command-line flags | Two text files: `samples.csv` and `config.ini` |
| Mode selection | Implicit from flags (or `--esmode`) | Implicit from which `samples.csv` columns you fill in |
| Multiple genomes | One run per genome | Many genomes in one `samples.csv`, all run from one `snakemake` invocation |
| Output location | `--workingdir` | `output/{sample_name}/results/` |
| Resume after failure | (very buggy) | Re-run the same `snakemake` command — picks up where it left off |
| Repeat masking | You provide a masked genome | Optional: leave the `genome_masked` column empty and BRAKER4 runs RepeatModeler2 + RepeatMasker for you |
| HPC submission | You wrap `braker.pl` in your own SLURM script | `snakemake --executor slurm` submits each rule as a separate cluster job |
| Container | One container, you call it manually | All containers managed by Snakemake's `--use-singularity` |
| QC outputs | None integrated | BUSCO, compleasm, OMArk, gffcompare integrated; HTML report generated automatically |

What does **not** change:

- The workflow. GeneMark trains on extrinsic evidence, AUGUSTUS trains on the GeneMark output, TSEBRA merges the two. The same algorithms with the same parameters as BRAKER3.
- Container image. The default `[containers] braker3_image` is `docker://teambraker/braker3:v3.0.10` — it is an updated version of the original BRAKER3 container that is no longer compatible with braker.pl and therefore note tagged as latest.
- Output formats. The final files are still `braker.gtf`, `braker.gff3`, `braker.aa`, `braker.codingseq`, `hintsfile.gff` — just collected into the `results/` subdirectory.

## Step 0 — prerequisites

Make sure you have:

- **Snakemake 8.x** (`snakemake --version` should report ≥ 8). Do not use Snakemake 7 — the SLURM executor architecture changed in 8 and BRAKER4 only tests against 8.
- **Singularity (or Apptainer) 3.x**. The Snakemake `--use-singularity` flag will pull the BRAKER3 container and several auxiliary containers (minimap2, gffcompare, AGAT, BUSCO, OMArk, barrnap, dfam/tetools, …) automatically on first run.
- **`snakemake-executor-plugin-slurm`** if you want to run on a SLURM cluster.
- **About 13 GB of disk for the Singularity image cache.** The BRAKER3 image alone is ~3 GB, the dfam/tetools repeat-masking image is ~2.5 GB, BUSCO is ~2.5 GB, AGAT is ~1 GB, plus several smaller helpers. Plus enough disk for your genome FASTAs, intermediates, and outputs (`output/` can grow to 10–50 GB per genome on real data — set `no_cleanup = 0` in `config.ini` to have it pruned to a few hundred MB after each run).

If you currently run `braker.pl` directly inside a Singularity shell, you already have steps 1 and 2 covered.

## Step 1 — clone the repository

```bash
git clone https://github.com/Gaius-Augustus/BRAKER4.git
cd BRAKER4
```

You will run snakemake from a working directory of your choice and point it at the `Snakefile` here. You do not need to install BRAKER4 — there is nothing to install. The Snakefile is the entry point.

## Step 2 — translate your `braker.pl` command

Take the command you currently use, and write down each flag. Then map each flag to either a `samples.csv` column or a `config.ini` parameter using the table below.

### Common BRAKER3 flag → BRAKER4 destination

| `braker.pl` flag | BRAKER4 destination |
|---|---|
| `--genome=genome.fa` | `samples.csv` column `genome` |
| `--bam=rnaseq.bam` (one or several) | `samples.csv` column `bam_files` (colon-separated for multiple) |
| `--prot_seq=proteins.fa` | `samples.csv` column `protein_fasta` (colon-separated for multiple) |
| `--rnaseq_sets_ids=SRR1,SRR2` | `samples.csv` column `sra_ids` |
| `--species=my_species` | `samples.csv` column `sample_name` (becomes the AUGUSTUS species AND the output directory name) |
| `--workingdir=/path/to/wd` | The directory you cd into before running `snakemake` (output goes to `output/{sample_name}/`) |
| `--threads=N` | `snakemake --cores N` (local) or `--cores N --jobs N` (SLURM) |
| `--gff3` | Always on. Both `.gtf` and `.gff3` are produced. |
| `--esmode` | Leave all evidence columns empty in `samples.csv`. ES mode is auto-detected. |
| `--fungus` | `config.ini` `[PARAMS]` `fungus = 1` |
| `--useexisting` | Not needed. Automatic resume is built in. Just re-run the same `snakemake` command. |
| `--softmasking` | Default in the container. No flag needed. |
| `--AUGUSTUS_CONFIG_PATH=...` | `config.ini` `[paths]` `augustus_config_path = ...` (default `augustus_config` works for most users) |
| `--PROTHINT_PATH`, `--GENEMARK_PATH`, `--TSEBRA_PATH`, `--PYTHON3_PATH`, etc. | Not needed. Everything lives in the container. |

### What about IsoSeq?

`braker.pl` does not support PacBio IsoSeq directly. If you were using a separate workflow for that, BRAKER4 has it built in:

- Pre-aligned IsoSeq BAM → `samples.csv` column `isoseq_bam`
- Unaligned IsoSeq FASTQ/FASTA → `samples.csv` column `isoseq_fastq` (BRAKER4 aligns it with minimap2)

Both require `protein_fasta` to also be set (this is enforced).

## Step 3 — create your `samples.csv`

Pick a working directory (it does not have to be inside the cloned repo):

```bash
mkdir ~/my_braker4_run
cd ~/my_braker4_run
```

Write a `samples.csv`. The header is fixed — copy it exactly:

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
```

Then add **one row per genome**. Leave a column empty if you do not have that input. Below are concrete translations of common BRAKER3 commands.

### ETP example (RNA-Seq BAM + proteins, the most common BRAKER3 setup)

If you used to run:

```bash
braker.pl --genome=genome.fa --bam=rnaseq.bam --prot_seq=proteins.fa \
          --species=my_species --workingdir=output --threads=48 --gff3
```

Your `samples.csv` is:

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_species,/path/to/genome.fa,,/path/to/proteins.fa,/path/to/rnaseq.bam,,,,,,,,eukaryota_odb12,
```

Notes:

- **`busco_lineage` is required.** Pick a clade-specific lineage like `arthropoda_odb12`, `viridiplantae_odb12`, `fungi_odb12`, etc. The pipeline uses it for compleasm and BUSCO QC. If you want the most generic option, use `eukaryota_odb12`.
- **Paths can be absolute** (recommended) **or relative** to the working directory you run `snakemake` from.
- The `sample_name` becomes the AUGUSTUS species name and the output directory name. Pick something stable — re-running with the same name resumes; changing it forces a fresh run.

### ET example (RNA-Seq only — BRAKER1)

```bash
braker.pl --genome=genome.fa --bam=rnaseq.bam --species=my_species --threads=48
```

becomes

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_species,/path/to/genome.fa,,,/path/to/rnaseq.bam,,,,,,,,eukaryota_odb12,
```

(Note: `protein_fasta` is empty, that is what triggers ET mode.)

### EP example (proteins only — BRAKER2)

```bash
braker.pl --genome=genome.fa --prot_seq=proteins.fa --species=my_species --threads=48
```

becomes

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_species,/path/to/genome.fa,,/path/to/proteins.fa,,,,,,,,,eukaryota_odb12,
```

### ES example (genome only — ab initio)

```bash
braker.pl --esmode --genome=genome.fa --species=my_species --threads=48
```

becomes

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_species,/path/to/genome.fa,,,,,,,,,,,eukaryota_odb12,
```

(All evidence columns empty → ES mode is detected automatically.)

### Multi-sample example

You can run several genomes in the same invocation:

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
fly,/data/fly.fa,,/data/orthodb_arthropoda.fa,/data/fly_rnaseq.bam,,,,,,,,arthropoda_odb12,
worm,/data/worm.fa,/data/worm_masked.fa,/data/orthodb_metazoa.fa,,/data/worm_R1.fq.gz,/data/worm_R2.fq.gz,,,,,,metazoa_odb12,
plant,/data/plant.fa,,/data/orthodb_viridiplantae.fa,,,,SRR123456:SRR789012,,,,,viridiplantae_odb12,
```

This single `samples.csv` annotates three different organisms in three different ways in one Snakemake invocation. You cannot do that with `braker.pl`.

## Step 4 — create your `config.ini`

In the same working directory, create `config.ini`:

```ini
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config

[PARAMS]
fungus = 0
min_contig = 10000
skip_optimize_augustus = 0
use_compleasm_hints = 1
skip_busco = 0
run_omark = 0
run_ncrna = 0

[SLURM_ARGS]
cpus_per_task = 48
mem_of_node = 120000
max_runtime = 4320
```

The defaults above are appropriate for most ETP runs on a 48-core SLURM node with 120 GB RAM and a 72-hour wall time. Adjust to your cluster.

A few flags you might want to flip:

| Flag | When to set it |
|---|---|
| `fungus = 1` | If your organism is a fungus. Sets the GeneMark branch-point model. |
| `skip_busco = 1` | BUSCO is the slowest QC step (~20 min on plants). Skip if you do not need it. |
| `run_omark = 1` | Run OMArk on the predicted proteins. Requires the LUCA.h5 database (~8.8 GB). |
| `run_ncrna = 1` | Predict ncRNAs (rRNA, tRNA, snoRNA, miRNA, lncRNA) in addition to coding genes. |
| `use_compleasm_hints = 0` | Keep BUSCO compleasm CDSpart hints out of the AUGUSTUS hintsfile. Closer to native braker.pl behavior. |
| `skip_optimize_augustus = 1` | Skip the slow AUGUSTUS parameter optimization step. Saves ~30 min on small test runs but slightly reduces accuracy. Recommended only for testing. |

For the full list see the README's [Description of selected configuration options](README.md#description-of-selected-configuration-options).

## Step 5 — run BRAKER4

Two cases — local workstation or HPC cluster.

### Local (workstation, no SLURM)

From your working directory (the one with `samples.csv` and `config.ini`):

```bash
snakemake \
    --snakefile /path/to/BRAKER4/Snakefile \
    --cores 8 \
    --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home" \
    --latency-wait 120 \
    --restart-times 3
```

- `--cores 8` → use 8 CPU threads.
- `--use-singularity` → run all rules inside their declared containers.
- `--singularity-prefix .singularity_cache` → cache pulled images here. Without this flag, every working directory gets its own copy of the (multi-GB) container image.
- `--singularity-args "-B /home"` → make `/home` visible inside the container. **If your data lives outside `/home`** (e.g. `/scratch`, `/data`, `/gpfs`), add those bind paths: `--singularity-args "-B /home -B /scratch -B /data"`. The container can only see directories that are explicitly bound.

The first run will pull all the containers BRAKER4 needs (~13 GB total: BRAKER3 ~3 GB, dfam/tetools repeat-masking ~2.5 GB, BUSCO ~2.5 GB, AGAT ~1 GB, plus smaller helpers). Give it 10-20 minutes depending on your network. Subsequent runs reuse the cached images. If you set `--singularity-prefix` to a stable shared path, multiple working directories on the same machine share one cache.

### HPC cluster (SLURM)

From your working directory:

```bash
snakemake \
    --snakefile /path/to/BRAKER4/Snakefile \
    --executor slurm \
    --default-resources slurm_partition=batch mem_mb=120000 \
    --cores 48 \
    --jobs 48 \
    --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home -B /scratch" \
    --latency-wait 120 \
    --restart-times 3
```

- `--executor slurm` → submit each rule as a separate `sbatch` job. Set `slurm_partition=` to your cluster's partition name.
- `--cores 48` → cores per individual SLURM job. Match your `cpus_per_task` in `config.ini`.
- `--jobs 48` → maximum number of SLURM jobs in flight at once.
- `mem_mb=120000` → memory per SLURM job (MB). Match `mem_of_node` in `config.ini`.

**Run snakemake itself in `screen` or `tmux`**, or submit it as a long-running SLURM job. The Snakemake driver process must stay alive for the entire pipeline duration — it submits the children, monitors them, and assembles the DAG. If you log out of the head node and the driver dies, the SLURM children keep going but no new jobs get submitted.

### Resume after a failure

If anything fails (out of memory, time limit, network glitch), just **re-run the same command**. Snakemake reads the existing partial outputs, figures out what is missing, and only re-runs what is needed. There is no `--continue` or `--useexisting` flag — resume is the default behavior.

## Step 6 — find your outputs

Your final files are in:

```
output/{sample_name}/results/
├── braker.gtf.gz          # final gene set, GTF
├── braker.gff3.gz         # same set, GFF3 (AGAT-formatted)
├── braker.aa.gz           # protein sequences
├── braker.codingseq.gz    # coding sequences
├── braker_utr.gtf.gz      # UTR-decorated gene models (if RNA-Seq evidence was given)
├── hintsfile.gff.gz       # all extrinsic hints
├── braker_report.html     # self-contained HTML report
├── braker_citations.bib   # BibTeX for all tools used
├── software_versions.tsv  # versions of all software used
├── gene_support.tsv       # per-gene evidence support
└── quality_control/
    ├── busco_summary.txt
    ├── compleasm_summary.txt
    ├── completeness.png
    ├── gene_set_statistics.txt
    └── (more plots and tables)
```

The mapping from BRAKER3 output to BRAKER4 output is one-to-one — `braker.gtf` is still `braker.gtf`. The main differences are:

1. **Outputs are gzipped** by default in `results/`. The `collect_results` rule gzips them on the way out. The unzipped versions also exist in `output/{sample_name}/` itself before cleanup.
2. **Outputs are in a `results/` subdirectory.** Everything else in `output/{sample_name}/` is intermediate state.
3. **The HTML report is new.** Open `braker_report.html` in a browser — it embeds all QC plots, gene-set statistics, software versions, and citations into a single self-contained file you can email to a collaborator.
4. **The BibTeX file is new.** All software cited in your run, ready to drop into a bibliography manager.

## Step 7 — sanity-check the run

Open `braker_report.html` in a browser. It will tell you:

- How many genes were predicted at each pipeline stage (training set → AUGUSTUS → TSEBRA-merged → final).
- BUSCO and compleasm completeness scores against your `busco_lineage`.
- Hint-support statistics (what fraction of predicted genes are backed by RNA-Seq / protein hints).
- Gene-set statistics (CDS lengths, exons per gene, mono-exonic vs multi-exonic ratio, etc.).
- The exact `braker.pl`-equivalent command you would have had to run, plus all software versions and citations.

If a step looks wrong (e.g. far fewer genes than expected, or BUSCO completeness is unexpectedly low), look at the per-rule logs in `logs/{sample_name}/<rule_name>/` to find what happened. Each rule writes its own log file there. On SLURM, the captured stdout/stderr of each job lives in `.snakemake/slurm_logs/rule_<rule_name>/<sample_name>/<jobid>.log`.

## Common gotchas when migrating

### "Singularity cannot find my data"

By default `--singularity-args "-B /home"` only binds `/home`. If your genome FASTA is in `/scratch/myproject/genome.fa`, the container cannot see it. Add the bind path:

```
--singularity-args "-B /home -B /scratch"
```

This is the most common first-run failure for users coming from BRAKER3.

### "I don't have a `busco_lineage` for my organism"

Use the most specific BUSCO odb12 lineage that contains your organism. If you really do not know which to pick, `eukaryota_odb12` works for any eukaryote (it is just less informative for QC than a clade-specific one). **You cannot leave the column empty** — the pipeline requires it for compleasm and best_by_compleasm. If you really do not want any BUSCO at all, you can set both `skip_busco = 1` and `use_compleasm_hints = 0` in `config.ini`, but compleasm itself still runs in protein mode for the `best_by_compleasm` rescue step (it is part of the gene-set merging logic, not optional QC).

### "I don't want my genome masked by RepeatModeler"

Provide your already-masked genome in the `genome_masked` column. The `genome` column should still point at the unmasked version — both are used at different stages. If `genome_masked` is filled, the RepeatModeler/RepeatMasker rules are skipped.

### "I want to override one config value without editing config.ini"

Every value in `[PARAMS]` and `[SLURM_ARGS]` can be overridden via an environment variable named `BRAKER4_<KEY_UPPER>`. For example:

```bash
BRAKER4_FUNGUS=1 BRAKER4_MAX_RUNTIME=240 snakemake --latency-wait 120 --restart-times 3 ...
```

Useful when you want to share one `config.ini` across multiple runs but tweak a few values per-run. The test suite uses this pattern: one shared `config.ini`, plus a small `scenario_overrides.sh` per scenario that exports the differences.

### "Snakemake says my SLURM job has no `mem_mb` resource"

Add `mem_mb=...` to `--default-resources`:

```
--default-resources slurm_partition=batch mem_mb=120000
```

`braker.pl` did not need this because you set memory yourself in the wrapping `sbatch` script. In BRAKER4, Snakemake submits each rule individually and needs to know how much memory to request.

### "I want to keep all intermediate files for debugging"

Set `no_cleanup = 1` in `[PARAMS]`. By default the `collect_results` rule moves the important outputs into `results/` and deletes the intermediate state. With `no_cleanup = 1`, everything in `output/{sample_name}/` is preserved. Disk-hungry — only enable for debugging.

### "I had a `--workingdir` separate from my source code"

In BRAKER4, your "working directory" is whatever directory you `cd` into before running snakemake. It contains `samples.csv` and `config.ini`, and Snakemake creates `output/`, `logs/`, `benchmarks/`, `.snakemake/` underneath it. Pick any directory you like — it does not have to be inside the cloned BRAKER4 repository. Reference the Snakefile via its absolute path (`--snakefile /path/to/BRAKER4/Snakefile`).

### "My runs hung / SLURM showed COMPLETING jobs forever"

If you killed a snakemake driver mid-run, leftover SLURM jobs may still be in the queue and stale lock files may be in `.snakemake/locks/`. Clean both before re-running:

```bash
# scancel any stale jobs (use squeue -u $USER to identify)
scancel <jobid> ...
# remove stale locks
rm -f .snakemake/locks/*
# now re-run snakemake
snakemake --latency-wait 120 --restart-times 3 ...
```

## Where to read more

- [README.md](README.md) — full reference: every config option, every mode, every output file.
- [test_scenarios/](test_scenarios/) — runnable HPC test scenarios for ET, EP, ETP, IsoSeq, dual, ES, the `athaliana` benchmark.
- [test_scenarios_local/](test_scenarios_local/) — same scenarios but for a local workstation (no SLURM).
- [TESTING_GUIDE.md](TESTING_GUIDE.md) — how to run the test suite end-to-end.
- The [multi-mode rulegraph](test_scenarios/scenario_12_multi_mode/rulegraph.png) — visualises the entire DAG including all 7 modes plus repeat masking and QC.

If you get stuck, the per-rule logs in `logs/{sample_name}/` and the SLURM job logs in `.snakemake/slurm_logs/` are usually enough to identify what went wrong. Open an issue at <https://github.com/Gaius-Augustus/BRAKER4/issues> if you cannot.
