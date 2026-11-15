#!/usr/bin/env bash
# check_bam_refs.sh BAM GENOME_FAI
#
# Compare the @SQ reference names of BAM with the sequence names of the
# genome (first column of GENOME_FAI, the headers kept by prepare_genome).
# Exits 1 with an ERROR on stderr when none of the BAM references are in
# the genome; prints WARNINGs when only some match or lengths differ.
# Without this check GeneMark-ETP's StringTie step silently produces an
# empty transcripts_merged.fasta and fails much later
# (BRAKER4 #76; Gaius-Augustus/BRAKER#877 for long-read BAMs).
# Used by check_bam_sorted (short reads) and check_isoseq_bam (long reads).
set -euo pipefail
bam=$1
fai=$2
samtools view -H "$bam" | awk -F'\t' -v fai="$fai" -v bam="$bam" '
    BEGIN { while ((getline l < fai) > 0) { split(l, f, "\t"); len[f[1]] = f[2] } }
    $1 == "@SQ" {
        sn = ""; ln = ""
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SN:/) sn = substr($i, 4)
            if ($i ~ /^LN:/) ln = substr($i, 4)
        }
        n++
        if (sn in len) { ok++; if (len[sn] != ln) diff++ }
        else if (miss++ < 5) ex = ex " " sn
    }
    END {
        if (n > 0 && ok == 0) {
            print "ERROR: none of the " n " reference sequences in " bam " match the genome FASTA headers." > "/dev/stderr"
            print "  BAM references, e.g.:" ex > "/dev/stderr"
            print "  Align the reads against the same genome FASTA (identical sequence names) that you pass to BRAKER4." > "/dev/stderr"
            exit 1
        }
        if (miss > 0) printf "WARNING: %d of %d BAM reference sequences are not in the genome, e.g.:%s\n", miss, n, ex > "/dev/stderr"
        if (diff > 0) printf "WARNING: %d BAM reference sequences have a different length than in the genome. Was the BAM made against another assembly version?\n", diff > "/dev/stderr"
    }'
