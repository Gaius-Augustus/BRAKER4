#!/bin/bash
# pyVARUS wrapper for Snakemake.
#
# Usage: run_varus_wrapper.sh <varus_dir> <genome> <species_name> <threads>
#                             <output_bam> <log> <use_logan> <max_batches>
#                             [R1.fq.gz,R2.fq.gz ...]
#
# species_name: "Genus species" for sampling from SRA; "" for no SRA sampling
#               (then at least one local library is needed).
# use_logan:    1 runs the pyVARUS Logan pre-screen, 0 passes --no-logan.
# max_batches:  upper bound on the number of batches (varus run --max-batches).
# R1,R2:        local gzipped FASTQ libraries (*.fastq.gz or *.fq.gz) that
#               pyVARUS samples from (varus run --fastq, pyVARUS >= 2.1.0).
#               With a species name, SRA runs and local libraries compete
#               for the batches.
#
# <varus_dir> may be a node-local scratch directory: the BAM is sorted and
# indexed there and then copied to <output_bam>; the run list and the stats
# are written next to <output_bam>. The log is appended to (the rule
# truncates it).
#
# NCBI Entrez email: set NCBI_EMAIL in the environment or it falls back to
# the lab address below (used only for rate-limit courtesy, not authentication).
# An empty NCBI_EMAIL (export NCBI_EMAIL=) sends no address.

set -euo pipefail

VARUS_DIR="$1"
GENOME="$2"
SPECIES_NAME="$3"
THREADS="$4"
OUTPUT_BAM="$5"
LOGFILE="$6"
USE_LOGAN="${7:-0}"
MAX_BATCHES="${8:-1000}"
shift $(( $# < 8 ? $# : 8 ))
FASTQ_LIBS=("$@")

EMAIL="${NCBI_EMAIL-katharina.hoff@uni-greifswald.de}"

if [ -z "$SPECIES_NAME" ] && [ "${#FASTQ_LIBS[@]}" -eq 0 ]; then
    echo "[ERROR] run_varus_wrapper.sh: neither a species name nor a local FASTQ library given" >&2
    exit 1
fi

# Resolve all paths to absolute
STARTDIR="$PWD"
GENOME_ABS=$(readlink -f "$GENOME")

[[ "$VARUS_DIR"  != /* ]] && VARUS_DIR_ABS="$STARTDIR/$VARUS_DIR"   || VARUS_DIR_ABS="$VARUS_DIR"
[[ "$OUTPUT_BAM" != /* ]] && OUTPUT_BAM_ABS="$STARTDIR/$OUTPUT_BAM" || OUTPUT_BAM_ABS="$OUTPUT_BAM"
[[ "$LOGFILE"    != /* ]] && LOGFILE_ABS="$STARTDIR/$LOGFILE"       || LOGFILE_ABS="$LOGFILE"

mkdir -p "$VARUS_DIR_ABS" "$(dirname "$OUTPUT_BAM_ABS")" "$(dirname "$LOGFILE_ABS")"

INDEX_DIR="$VARUS_DIR_ABS/genome"

echo "[INFO] pyVARUS wrapper start"       >> "$LOGFILE_ABS"
echo "[INFO] Species:  ${SPECIES_NAME:-none (no SRA sampling)}" >> "$LOGFILE_ABS"
echo "[INFO] Local FASTQ libraries: ${#FASTQ_LIBS[@]}" >> "$LOGFILE_ABS"
echo "[INFO] Max batches: $MAX_BATCHES"   >> "$LOGFILE_ABS"
echo "[INFO] Genome:   $GENOME_ABS"       >> "$LOGFILE_ABS"
echo "[INFO] Threads:  $THREADS"          >> "$LOGFILE_ABS"
echo "[INFO] Out BAM:  $OUTPUT_BAM_ABS"   >> "$LOGFILE_ABS"
echo "[INFO] Logan:    $USE_LOGAN"        >> "$LOGFILE_ABS"

LOGAN_ARGS=()
[[ "$USE_LOGAN" == "1" ]] || LOGAN_ARGS=(--no-logan)

# Local libraries: absolute paths without resolving symbolic links, because
# pyVARUS takes the library name and the file type from the file name.
FASTQ_ARGS=()
for lib in "${FASTQ_LIBS[@]}"; do
    r1="${lib%%,*}"
    r2="${lib#*,}"
    for f in "$r1" "$r2"; do
        if [ ! -s "$f" ]; then
            echo "[ERROR] FASTQ file $f is missing or empty" >> "$LOGFILE_ABS"
            exit 1
        fi
    done
    FASTQ_ARGS+=(--fastq "$(realpath -s "$r1"),$(realpath -s "$r2")")
done

# Step 1: fetch SRA run list from NCBI (only with a species name).
# A species without RNA-Seq runs in SRA ends the job, unless there are local
# libraries: then pyVARUS samples from those alone.
RUNLIST_ARGS=()
SRA_USED=0
if [ -n "$SPECIES_NAME" ]; then
    echo "[INFO] Fetching run list..." >> "$LOGFILE_ABS"
    EMAIL_ARGS=()
    [ -z "$EMAIL" ] || EMAIL_ARGS=(--email "$EMAIL")
    RUNLIST_RC=0
    varus runlist "$SPECIES_NAME" \
        --outdir "$VARUS_DIR_ABS" \
        "${EMAIL_ARGS[@]}" \
        >> "$LOGFILE_ABS" 2>&1 || RUNLIST_RC=$?
    if [ "$RUNLIST_RC" -eq 0 ]; then
        RUNLIST_ARGS=(--runlist "$VARUS_DIR_ABS/Runlist.tsv")
        SRA_USED=1
    elif [ "${#FASTQ_LIBS[@]}" -gt 0 ] && grep -q "No SRA RNA-seq runs found" "$LOGFILE_ABS"; then
        echo "[WARNING] SRA has no RNA-Seq runs for $SPECIES_NAME. Sampling from the local FASTQ libraries only." >> "$LOGFILE_ABS"
    else
        echo "[ERROR] varus runlist ended with exit status $RUNLIST_RC" >> "$LOGFILE_ABS"
        exit "$RUNLIST_RC"
    fi
fi

# Step 2: build HISAT2 index
echo "[INFO] Building HISAT2 index..." >> "$LOGFILE_ABS"
varus index "$GENOME_ABS" \
    --outdir  "$INDEX_DIR" \
    --threads "$THREADS" \
    >> "$LOGFILE_ABS" 2>&1

# Step 3: online sampling loop (download or extract, align, score).
# --keep-bam: pyVARUS >= 2.0.0 deletes VARUS.bam after its own StringTie
# assembly; BRAKER4 needs the alignments. Exit status 4 means that this
# assembly failed (the BAM is kept); BRAKER4 does not use it. Exit status 3
# means that no batch passed the quality gate.
echo "[INFO] Running pyVARUS sampling loop..." >> "$LOGFILE_ABS"
VARUS_RC=0
varus run "${SPECIES_NAME:-local FASTQ libraries}" "$GENOME_ABS" \
    "${RUNLIST_ARGS[@]}" \
    "${FASTQ_ARGS[@]}" \
    --index   "$INDEX_DIR/hisatidx" \
    --outdir  "$VARUS_DIR_ABS" \
    --threads "$THREADS" \
    --max-batches "$MAX_BATCHES" \
    --keep-bam \
    "${LOGAN_ARGS[@]}" \
    >> "$LOGFILE_ABS" 2>&1 || VARUS_RC=$?
if [ "$VARUS_RC" -eq 4 ] && [ -s "$VARUS_DIR_ABS/VARUS.bam" ]; then
    echo "[WARNING] The StringTie assembly inside pyVARUS failed (exit status 4). BRAKER4 does not use it and continues with VARUS.bam." >> "$LOGFILE_ABS"
elif [ "$VARUS_RC" -eq 3 ]; then
    echo "[ERROR] pyVARUS accepted no batch: in every batch fewer than 5% of the reads mapped uniquely to the genome. Check that the reads belong to this genome." >> "$LOGFILE_ABS"
    exit 1
elif [ "$VARUS_RC" -ne 0 ]; then
    echo "[ERROR] varus run ended with exit status $VARUS_RC" >> "$LOGFILE_ABS"
    exit "$VARUS_RC"
fi
if [ ! -s "$VARUS_DIR_ABS/VARUS.bam" ]; then
    echo "[ERROR] pyVARUS wrote no VARUS.bam: no usable alignments" >> "$LOGFILE_ABS"
    exit 1
fi

# Step 4: sort and index in <varus_dir>, then copy to <output_bam>
SORTED_BAM="$VARUS_DIR_ABS/varus.sorted.bam"
echo "[INFO] Sorting BAM..." >> "$LOGFILE_ABS"
samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$VARUS_DIR_ABS/VARUS.bam" \
    >> "$LOGFILE_ABS" 2>&1

echo "[INFO] Indexing BAM (CSI)..." >> "$LOGFILE_ABS"
samtools index -c -@ "$THREADS" "$SORTED_BAM" >> "$LOGFILE_ABS" 2>&1

# Same file when <varus_dir> is the output directory (no scratch)
if [ ! "$SORTED_BAM" -ef "$OUTPUT_BAM_ABS" ]; then
    cp "$SORTED_BAM" "$OUTPUT_BAM_ABS"
    cp "$SORTED_BAM.csi" "$OUTPUT_BAM_ABS.csi"
fi

# Step 5: copy run list to standard output location for reports (SRA only)
if [ "$SRA_USED" -eq 1 ]; then
    cp "$VARUS_DIR_ABS/Runlist.tsv" "$(dirname "$OUTPUT_BAM_ABS")/varus_runlist.tsv"
fi

# Step 6: write summary stats
TOTAL_READS=$(samtools view -c "$SORTED_BAM" 2>/dev/null || echo 0)
{
    echo "pyVARUS Run Summary"
    echo "==================="
    if [ -n "$SPECIES_NAME" ]; then
        echo "Species: $SPECIES_NAME"
    fi
    if [ "$SRA_USED" -eq 1 ]; then
        echo "SRA runs in runlist: $(tail -n +2 "$VARUS_DIR_ABS/Runlist.tsv" | wc -l)"
    elif [ -n "$SPECIES_NAME" ]; then
        echo "SRA: no RNA-Seq runs found for this species"
    fi
    if [ "${#FASTQ_LIBS[@]}" -gt 0 ]; then
        echo "Local FASTQ libraries: ${#FASTQ_LIBS[@]}"
    fi
    echo "Max batches: $MAX_BATCHES"
    echo "Total reads in BAM: $TOTAL_READS"
} > "$(dirname "$OUTPUT_BAM_ABS")/varus_stats.txt"

echo "[INFO] pyVARUS wrapper complete — $TOTAL_READS reads in $OUTPUT_BAM_ABS" >> "$LOGFILE_ABS"
