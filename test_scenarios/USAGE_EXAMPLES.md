# Usage Examples - BRAKER4 Pipeline

## Quick Reference

| Data You Have | Scenario | Directory |
|---|---|---|
| Genome only (ab initio) | 01 / 08 | `scenario_01_es` / `scenario_08_es` |
| Genome + Proteins | 02 / 09 | `scenario_02_ep` / `scenario_09_ep_masking` |
| Genome + FASTQ | 03 | `scenario_03_et_fastq` |
| Genome + BAM + Proteins | 04 | `scenario_04_etp_bam` |
| Genome + SRA IDs + Proteins | 11 | `scenario_11_etp_sra` |
| Genome + VARUS | 10 | `scenario_10_et_varus` |
| Genome + IsoSeq BAM + Proteins | 05 | `scenario_05_isoseq_bam` |
| Genome + IsoSeq FASTQ + Proteins | 06 | `scenario_06_isoseq_fastq` |
| Genome + BAM + IsoSeq + Proteins | 07 | `scenario_07_dual` |

Then edit `samples.csv` and run:

```bash
DRY_RUN=false bash run_test.sh
```

---

## Example 1: Genome + Proteins (EP mode)

```bash
cd test_scenarios_local/scenario_02_ep
nano samples.csv
```

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,,/path/to/proteins.fasta,,,,,,,,,fungi_odb12,
```

```bash
DRY_RUN=false CORES=16 bash run_test.sh
```

---

## Example 2: Genome + BAM + Proteins (ETP mode)

```bash
cd test_scenarios_local/scenario_04_etp_bam
nano samples.csv
```

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,,/path/to/proteins.fasta,/path/to/rna1.bam:/path/to/rna2.bam,,,,,,,,,metazoa_odb12,
```

**Note**: Multiple BAM files are separated by colons (`:`)

```bash
DRY_RUN=false CORES=32 bash run_test.sh
```

---

## Example 3: Genome + FASTQ reads (ET mode)

```bash
cd test_scenarios_local/scenario_03_et_fastq
nano samples.csv
```

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,,,,,/path/to/reads_R1.fastq.gz,/path/to/reads_R2.fastq.gz,,,,,,eukaryota_odb12,
```

The pipeline will automatically align FASTQ reads with HISAT2 and run BRAKER.

```bash
DRY_RUN=false CORES=24 bash run_test.sh
```

---

## Example 4: Download RNA-Seq from SRA (ETP mode)

Find your SRA accessions at [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra).

```bash
cd test_scenarios/scenario_11_etp_sra
nano samples.csv
```

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,,/path/to/proteins.fasta,,,,SRR1234567:SRR1234568,,,,,eukaryota_odb12,
```

The pipeline will download SRA files, convert to FASTQ, align with HISAT2, and run BRAKER.

```bash
DRY_RUN=false CORES=32 bash run_test.sh
```

---

## Example 5: VARUS auto-discovery (ET mode)

```bash
cd test_scenarios/scenario_10_et_varus
nano samples.csv
```

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,,,,,,,,Drosophila,melanogaster,,,diptera_odb12,
```

VARUS will automatically search NCBI SRA, select high-quality datasets, download and align them.

```bash
# Requires internet access
DRY_RUN=false CORES=48 bash run_test.sh
```

---

## Example 6: IsoSeq BAM + Proteins

```bash
cd test_scenarios_local/scenario_05_isoseq_bam
nano samples.csv
```

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,,/path/to/proteins.fasta,,,,,,,,/path/to/isoseq.bam,,mammalia_odb12,
```

**Note**: IsoSeq BAM must already be aligned to the genome. Proteins are always required for IsoSeq mode.

```bash
DRY_RUN=false CORES=24 bash run_test.sh
```

---

## Example 7: Unaligned IsoSeq FASTA/FASTQ + Proteins

```bash
cd test_scenarios_local/scenario_06_isoseq_fastq
nano samples.csv
```

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,,/path/to/proteins.fasta,,,,,,,,,/path/to/isoseq.fastq.gz,mammalia_odb12,
```

The pipeline will align IsoSeq reads with `minimap2 -ax splice:hq` before running BRAKER.

```bash
DRY_RUN=false CORES=24 bash run_test.sh
```

---

## Example 8: Pre-masked genome

If you already ran RepeatMasker, provide both genome and masked genome to skip masking:

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
my_genome,/path/to/genome.fasta,/path/to/genome.masked.fasta,/path/to/proteins.fasta,,,,,,,,,,fungi_odb12,
```

This saves hours by skipping RepeatModeler2 + RepeatMasker.

---

## Multi-Sample Example

Create a single `samples.csv` with multiple rows:

```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf
species1,/data/sp1/genome.fa,,/data/proteins.fa,/data/sp1/rna.bam,,,,,,,,,fungi_odb12,
species2,/data/sp2/genome.fa,,/data/proteins.fa,,/data/sp2/R1.fq,/data/sp2/R2.fq,,,,,,eukaryota_odb12,
species3,/data/sp3/genome.fa,,/data/proteins.fa,,,,,,Genus,species,,,metazoa_odb12,
```

All three will be processed in parallel.

---

## Common Configuration Adjustments

### Increase timeout for large genomes

Edit `config.ini`. Use a memory value that fits your cluster's largest standard node — 120 GB is the recommended default and runs on most HPC partitions.

```ini
[SLURM_ARGS]
cpus_per_task = 48
mem_of_node = 120000    # 120 GB — fits most HPC nodes
max_runtime = 7200      # 5 days
```

Only raise `mem_of_node` above 120 GB if you have verified that your target SLURM partitions have nodes with that much RAM. Requesting more RAM than any node provides will cause SLURM to reject the job with `BadConstraints`.

### Enable AUGUSTUS optimization

Edit `config.ini`:

```ini
[PARAMS]
skip_optimize_augustus = 0    # Enable optimization (slower but better)
```

---

## Checking Results

```bash
# View final gene predictions
less output/my_genome/braker/braker.gtf

# Check protein sequences
less output/my_genome/braker/braker.aa

# View quality metrics
cat output/my_genome/qc/compleasm_output/summary.txt

# Count genes predicted
grep -c "gene" output/my_genome/braker/braker.gtf
```

---

## Troubleshooting

### Validation (dry-run) before execution

Always test first:

```bash
bash run_test.sh    # Dry-run (validates config)
```

If successful:

```bash
DRY_RUN=false bash run_test.sh    # Full execution
```

### Resume failed run

Snakemake automatically resumes from where it failed:

```bash
# Just run again
DRY_RUN=false bash run_test.sh
```
