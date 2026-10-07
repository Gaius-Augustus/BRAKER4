#!/bin/bash
# pyVARUS wrapper for Snakemake.
#
# Usage: run_varus_wrapper.sh <varus_dir> <genome> <genus> <species>
#                             <threads> <output_bam> <log> [use_logan]
#
# use_logan: 1 runs the pyVARUS Logan pre-screen, 0 (default) passes --no-logan.
#
# <varus_dir> may be a node-local scratch directory: the BAM is sorted and
# indexed there and then copied to <output_bam>; the run list and the stats
# are written next to <output_bam>. The log is appended to (the rule
# truncates it).
#
# NCBI Entrez email: set NCBI_EMAIL in the environment or it falls back to
# the lab address below (used only for rate-limit courtesy, not authentication).

set -euo pipefail

VARUS_DIR="$1"
GENOME="$2"
GENUS="$3"
SPECIES="$4"
THREADS="$5"
OUTPUT_BAM="$6"
LOGFILE="$7"
USE_LOGAN="${8:-0}"

EMAIL="${NCBI_EMAIL:-katharina.hoff@uni-greifswald.de}"
SPECIES_NAME="$GENUS $SPECIES"

# Resolve all paths to absolute
STARTDIR="$PWD"
GENOME_ABS=$(readlink -f "$GENOME")

[[ "$VARUS_DIR"  != /* ]] && VARUS_DIR_ABS="$STARTDIR/$VARUS_DIR"   || VARUS_DIR_ABS="$VARUS_DIR"
[[ "$OUTPUT_BAM" != /* ]] && OUTPUT_BAM_ABS="$STARTDIR/$OUTPUT_BAM" || OUTPUT_BAM_ABS="$OUTPUT_BAM"
[[ "$LOGFILE"    != /* ]] && LOGFILE_ABS="$STARTDIR/$LOGFILE"       || LOGFILE_ABS="$LOGFILE"

mkdir -p "$VARUS_DIR_ABS" "$(dirname "$OUTPUT_BAM_ABS")" "$(dirname "$LOGFILE_ABS")"

INDEX_DIR="$VARUS_DIR_ABS/genome"

echo "[INFO] pyVARUS wrapper start"       >> "$LOGFILE_ABS"
echo "[INFO] Species:  $SPECIES_NAME"     >> "$LOGFILE_ABS"
echo "[INFO] Genome:   $GENOME_ABS"       >> "$LOGFILE_ABS"
echo "[INFO] Threads:  $THREADS"          >> "$LOGFILE_ABS"
echo "[INFO] Out BAM:  $OUTPUT_BAM_ABS"   >> "$LOGFILE_ABS"
echo "[INFO] Logan:    $USE_LOGAN"        >> "$LOGFILE_ABS"

LOGAN_ARGS=()
[[ "$USE_LOGAN" == "1" ]] || LOGAN_ARGS=(--no-logan)

# Step 1: fetch SRA run list from NCBI
echo "[INFO] Fetching run list..." >> "$LOGFILE_ABS"
varus runlist "$SPECIES_NAME" \
    --outdir "$VARUS_DIR_ABS" \
    --email  "$EMAIL" \
    >> "$LOGFILE_ABS" 2>&1

# Step 2: build HISAT2 index
echo "[INFO] Building HISAT2 index..." >> "$LOGFILE_ABS"
varus index "$GENOME_ABS" \
    --outdir  "$INDEX_DIR" \
    --threads "$THREADS" \
    >> "$LOGFILE_ABS" 2>&1

# Step 3: online sampling loop — download, align, score
echo "[INFO] Running pyVARUS sampling loop..." >> "$LOGFILE_ABS"
varus run "$SPECIES_NAME" "$GENOME_ABS" \
    --runlist "$VARUS_DIR_ABS/Runlist.tsv" \
    --index   "$INDEX_DIR/hisatidx" \
    --outdir  "$VARUS_DIR_ABS" \
    --threads "$THREADS" \
    "${LOGAN_ARGS[@]}" \
    >> "$LOGFILE_ABS" 2>&1

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

# Step 5: copy run list to standard output location for reports
cp "$VARUS_DIR_ABS/Runlist.tsv" "$(dirname "$OUTPUT_BAM_ABS")/varus_runlist.tsv"

# Step 6: write summary stats
N_SRA=$(tail -n +2 "$VARUS_DIR_ABS/Runlist.tsv" | wc -l)
TOTAL_READS=$(samtools view -c "$SORTED_BAM" 2>/dev/null || echo 0)
cat > "$(dirname "$OUTPUT_BAM_ABS")/varus_stats.txt" <<STATSEOF
pyVARUS Run Summary
===================
Species: $SPECIES_NAME
SRA runs in runlist: $N_SRA
Total reads in BAM: $TOTAL_READS
STATSEOF

echo "[INFO] pyVARUS wrapper complete — $TOTAL_READS reads in $OUTPUT_BAM_ABS" >> "$LOGFILE_ABS"
