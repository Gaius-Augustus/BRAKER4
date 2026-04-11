#!/bin/bash
# VARUS wrapper script for Snakemake
#
# VARUS internally changes its working directory (chdir), which conflicts
# with Snakemake's directory management. This wrapper:
# 1. Resolves all paths to absolute before any chdir
# 2. cd into the VARUS working directory
# 3. Creates a complete VARUSparameters.txt
# 4. Runs runVARUS.pl
# 5. Sorts and indexes the output BAM
#
# Usage: run_varus_wrapper.sh <varus_dir> <genome> <genus> <species> <threads> <output_bam> <log>

set -euo pipefail

VARUS_DIR="$1"
GENOME="$2"
GENUS="$3"
SPECIES="$4"
THREADS="$5"
OUTPUT_BAM="$6"
LOGFILE="$7"

# Resolve all paths to absolute (before any chdir)
STARTDIR="$PWD"
GENOME_ABS=$(readlink -f "$GENOME")

# Build absolute paths from relative ones
if [[ "$VARUS_DIR" != /* ]]; then
    VARUS_DIR_ABS="$STARTDIR/$VARUS_DIR"
else
    VARUS_DIR_ABS="$VARUS_DIR"
fi

if [[ "$OUTPUT_BAM" != /* ]]; then
    OUTPUT_BAM_ABS="$STARTDIR/$OUTPUT_BAM"
else
    OUTPUT_BAM_ABS="$OUTPUT_BAM"
fi

if [[ "$LOGFILE" != /* ]]; then
    LOGFILE_ABS="$STARTDIR/$LOGFILE"
else
    LOGFILE_ABS="$LOGFILE"
fi

# Create directories
mkdir -p "$VARUS_DIR_ABS"
mkdir -p "$(dirname "$OUTPUT_BAM_ABS")"
mkdir -p "$(dirname "$LOGFILE_ABS")"

echo "[INFO] Starting VARUS wrapper" > "$LOGFILE_ABS"
echo "[INFO] VARUS directory: $VARUS_DIR_ABS" >> "$LOGFILE_ABS"
echo "[INFO] Genome (absolute): $GENOME_ABS" >> "$LOGFILE_ABS"
echo "[INFO] Species: $GENUS $SPECIES" >> "$LOGFILE_ABS"
echo "[INFO] Threads: $THREADS" >> "$LOGFILE_ABS"
echo "[INFO] Output BAM: $OUTPUT_BAM_ABS" >> "$LOGFILE_ABS"

SPECIES_DIR="${GENUS}_${SPECIES}"

# Change into VARUS working directory (required because VARUS does chdir internally)
cd "$VARUS_DIR_ABS"

# Create complete VARUSparameters.txt with all required parameters.
# runVARUS.pl copies this to the species subdir and adjusts paths.
# The VARUS binary reads the species-level copy - ALL parameters must be present.
cat > VARUSparameters.txt <<VARUSEOF
--aligner HISAT
--batchSize 75000
--blockSize 5000
--components 1
--cost 0.001
--deleteLater 0
--estimator 2
--exportObservationsToFile 0
--exportParametersToFile 0
--fastqDumpCall fastq-dump
--genomeDir .
--genomeFaFile ./genome.fa
--lambda 2
--loadAllOnce 0
--maxBatches 600
--mergeThreshold 2
--outFileNamePrefix ./
--pathToParameters ./VARUSparameters.txt
--pathToRuns ./
--pathToVARUS /opt/VARUS
--profitCondition 0
--pseudoCount 1
--qualityThreshold 5
--randomSeed 1
--readParametersFromFile 1
--runThreadN $THREADS
--verbosityDebug 0
VARUSEOF

echo "[INFO] Created VARUSparameters.txt with all parameters" >> "$LOGFILE_ABS"

# Symlink genome into VARUS directory so runVARUS.pl finds it
if [ ! -e genome.fa ]; then
    ln -sf "$GENOME_ABS" genome.fa
fi

echo "[INFO] Running runVARUS.pl..." >> "$LOGFILE_ABS"

# Run runVARUS.pl which creates the species subdir, builds index,
# adjusts VARUSparameters.txt, then calls the VARUS binary.
# runVARUS.pl does chdir into the species subdir before calling VARUS.
runVARUS.pl \
    --aligner=HISAT \
    --runThreadN="$THREADS" \
    --speciesGenome=genome.fa \
    --readFromTable=0 \
    --createindex=1 \
    --verbosity=5 \
    --latinGenus="$GENUS" \
    --latinSpecies="$SPECIES" \
    --varusParameters=VARUSparameters.txt \
    --outFileDir=. \
    >> "$LOGFILE_ABS" 2>&1 || {
        echo "[WARNING] runVARUS.pl exited with non-zero status, checking for output..." >> "$LOGFILE_ABS"
    }

# runVARUS.pl adjusts VARUSparameters.txt with paths relative to outFileDir,
# but then chdir's into the species subdir before calling the VARUS binary.
# This means the paths in VARUSparameters.txt are wrong from the binary's CWD.
# Fix: overwrite with paths relative to the species subdir (current dir of binary).
if [ -d "$VARUS_DIR_ABS/$SPECIES_DIR" ]; then
    echo "[INFO] Fixing VARUSparameters.txt paths for species subdir..." >> "$LOGFILE_ABS"
    cat > "$VARUS_DIR_ABS/$SPECIES_DIR/VARUSparameters.txt" <<FIXEOF
--aligner HISAT
--batchSize 75000
--blockSize 5000
--components 1
--cost 0.001
--deleteLater 0
--estimator 2
--exportObservationsToFile 0
--exportParametersToFile 0
--fastqDumpCall fastq-dump
--genomeDir ./genome/
--genomeFaFile ../genome.fa
--lambda 2
--loadAllOnce 0
--maxBatches 600
--mergeThreshold 2
--outFileNamePrefix ./
--pathToParameters ./VARUSparameters.txt
--pathToRuns ./
--pathToVARUS /opt/VARUS
--profitCondition 0
--pseudoCount 1
--qualityThreshold 5
--randomSeed 1
--readParametersFromFile 1
--runThreadN $THREADS
--verbosityDebug 0
FIXEOF

    # Now run the VARUS binary directly from the species subdir
    echo "[INFO] Running VARUS binary directly from species subdir..." >> "$LOGFILE_ABS"
    cd "$VARUS_DIR_ABS/$SPECIES_DIR"
    /opt/VARUS/Implementation/./VARUS >> "$LOGFILE_ABS" 2>&1 || {
        echo "[WARNING] VARUS binary exited with non-zero status" >> "$LOGFILE_ABS"
    }
    cd "$VARUS_DIR_ABS"
fi

# Find the VARUS BAM output — it's inside the species subdirectory
VARUS_BAM=""
for candidate in \
    "$VARUS_DIR_ABS/$SPECIES_DIR/VARUS.bam" \
    "$VARUS_DIR_ABS/$SPECIES_DIR/accepted_hits.bam" \
    "$VARUS_DIR_ABS/accepted_hits.bam" \
    "$VARUS_DIR_ABS/VARUS.bam"; do
    if [ -f "$candidate" ]; then
        VARUS_BAM="$candidate"
        break
    fi
done

if [ -z "$VARUS_BAM" ]; then
    echo "[ERROR] VARUS did not produce a BAM file" >> "$LOGFILE_ABS"
    echo "[ERROR] Searched for VARUS.bam and accepted_hits.bam in:" >> "$LOGFILE_ABS"
    echo "[ERROR]   $VARUS_DIR_ABS/$SPECIES_DIR/" >> "$LOGFILE_ABS"
    echo "[ERROR]   $VARUS_DIR_ABS/" >> "$LOGFILE_ABS"
    echo "[ERROR] Files found:" >> "$LOGFILE_ABS"
    find "$VARUS_DIR_ABS" -name "*.bam" -o -name "*.sam" >> "$LOGFILE_ABS" 2>&1 || true
    ls -laR "$VARUS_DIR_ABS/$SPECIES_DIR/" >> "$LOGFILE_ABS" 2>&1 || true
    exit 1
fi

echo "[INFO] VARUS produced: $VARUS_BAM" >> "$LOGFILE_ABS"

# Sort the BAM (VARUS output may not be coordinate-sorted)
echo "[INFO] Sorting BAM..." >> "$LOGFILE_ABS"
samtools sort -@ "$THREADS" -o "$OUTPUT_BAM_ABS" "$VARUS_BAM" 2>> "$LOGFILE_ABS"

# Index the sorted BAM
echo "[INFO] Indexing BAM..." >> "$LOGFILE_ABS"
samtools index -@ "$THREADS" "$OUTPUT_BAM_ABS" 2>> "$LOGFILE_ABS"

# Extract VARUS run list: which SRA experiments were downloaded and how many reads
RUNLIST_FILE="$(dirname "$OUTPUT_BAM_ABS")/varus_runlist.tsv"
echo -e "SRA_accession\treads_downloaded\tstatus" > "$RUNLIST_FILE"
if [ -d "$VARUS_DIR_ABS/$SPECIES_DIR" ]; then
    for sra_dir in "$VARUS_DIR_ABS/$SPECIES_DIR"/SRR* "$VARUS_DIR_ABS/$SPECIES_DIR"/ERR* "$VARUS_DIR_ABS/$SPECIES_DIR"/DRR*; do
        if [ -d "$sra_dir" ]; then
            accession=$(basename "$sra_dir")
            # Count FASTQ reads if available
            n_reads=0
            for fq in "$sra_dir"/*.fastq "$sra_dir"/*.fq "$sra_dir"/*.fastq.gz; do
                if [ -f "$fq" ]; then
                    if [[ "$fq" == *.gz ]]; then
                        n_reads=$((n_reads + $(zcat "$fq" 2>/dev/null | wc -l) / 4))
                    else
                        n_reads=$((n_reads + $(wc -l < "$fq") / 4))
                    fi
                fi
            done
            echo -e "${accession}\t${n_reads}\tdownloaded" >> "$RUNLIST_FILE"
        fi
    done
fi
N_SRA=$(tail -n +2 "$RUNLIST_FILE" | wc -l)
echo "[INFO] VARUS downloaded $N_SRA SRA experiments" >> "$LOGFILE_ABS"
echo "[INFO] Run list: $RUNLIST_FILE" >> "$LOGFILE_ABS"

TOTAL_MAPPED=$(samtools view -c "$OUTPUT_BAM_ABS" 2>/dev/null || echo 0)
TOTAL_READS=$(samtools view -c -F 4 "$OUTPUT_BAM_ABS" 2>/dev/null || echo 0)

# Save VARUS summary stats
VARUS_STATS="$(dirname "$OUTPUT_BAM_ABS")/varus_stats.txt"
cat > "$VARUS_STATS" << STATSEOF
VARUS Run Summary
=================
Species: $GENUS $SPECIES
SRA experiments downloaded: $N_SRA
Total reads in BAM: $TOTAL_MAPPED
Mapped reads: $TOTAL_READS
Batch size: 75000
Max batches: 600
STATSEOF

echo "[INFO] VARUS wrapper complete" >> "$LOGFILE_ABS"
echo "[INFO] Output: $OUTPUT_BAM_ABS" >> "$LOGFILE_ABS"
echo "[INFO] Reads mapped: $TOTAL_MAPPED" >> "$LOGFILE_ABS"
