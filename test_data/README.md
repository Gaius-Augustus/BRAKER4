# Test Data

Small test dataset for the BRAKER4 pipeline. The test genome corresponds to the last 1,000,000 nucleotides of *Arabidopsis thaliana* chromosome Chr5, split into 8 artificial contigs. RNA-Seq alignments were obtained by VARUS. Protein sequences are a subset of OrthoDB.

This test data is not compiled for optimal prediction accuracy but for quickly testing pipeline components.

## Files included in git

| File | Size | Description |
|------|------|-------------|
| `genome.fa` | 993 KB | Small genome fragment (8 contigs) |
| `proteins.fa` | 2.6 MB | 5,733 protein sequences (OrthoDB subset) |
| `isoseq.bam` | 14 MB | PacBio IsoSeq reads aligned to test genome |
| `isoseq.fastq.gz` | 672 KB | Unaligned IsoSeq reads (for minimap2 scenarios) |
| `isoseq_lib2.bam` | symlink | Symlink to `isoseq.bam` (for multi-library testing) |
| `reference.gtf` | 570 KB | Reference annotation for gffcompare evaluation |
| `RNAseq.hints` | 51 KB | Pre-computed intron hints |

## Files downloaded on demand

Run `bash test_data/download_test_data.sh` to download these files. They are too large for git.

| File | Size | Description |
|------|------|-------------|
| `RNAseq.bam` | 164 MB | RNA-Seq alignments (from Georgia Tech server) |
| `reads_R1.fastq.gz` | 77 MB | Paired-end FASTQ R1 (extracted from RNAseq.bam) |
| `reads_R2.fastq.gz` | 78 MB | Paired-end FASTQ R2 (extracted from RNAseq.bam) |
| `Arabidopsis_thaliana.fa` | ~120 MB | Full A. thaliana TAIR10.1 genome (for repeat masking scenarios) |
| `Viridiplantae.fa` | ~60 MB | Viridiplantae ODB12 proteins (for repeat masking scenarios) |
| `LUCA.h5` | ~15 GB | OMAmer database for OMArk evaluation (optional) |

Note: The `genome_masked.fasta` files referenced in test scenario `samples.csv` files point to `genome.fa` because the small test genome is already effectively masked. On real data, use a properly soft-masked genome.

## Download script

```bash
bash test_data/download_test_data.sh
```

This script downloads all files listed above and creates the `isoseq_lib2.bam` symlink. It is safe to re-run (skips files that already exist).

## Source

Test data originates from the [BRAKER example data](https://github.com/Gaius-Augustus/BRAKER/tree/master/example), extended with IsoSeq and reference annotation files for BRAKER4 testing.
