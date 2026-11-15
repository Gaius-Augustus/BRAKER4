"""
Check if BAM file is coordinate-sorted and sort if necessary.

This rule fixes a common BRAKER issue where unsorted BAM files cause
cryptic failures hours into the pipeline. It also parallelizes sorting
(original BRAKER uses single-threaded sorting).

Input:
    - BAM file (may or may not be sorted)

Output:
    - Coordinate-sorted BAM file
    - Index file (.csi)

Container: teambraker/braker3:latest (contains samtools)
"""

def get_input_bam(wildcards):
    """Get the input BAM file path for a given bam_id."""
    bam_files = get_bam_files(wildcards.sample)
    bam_ids = get_bam_ids(wildcards.sample)
    # Find the BAM file corresponding to this bam_id
    for bam_path, bam_id in zip(bam_files, bam_ids):
        if bam_id == wildcards.bam_id:
            return bam_path
    raise ValueError(f"No BAM file found for bam_id {wildcards.bam_id}")

rule check_bam_sorted:
    """Link a coordinate-sorted BAM, or sort it.

    Scratch: when the BAM must be sorted, the sorted BAM, its .csi and the
    sort chunks are written to a private directory on the node-local disk
    (scripts/tmp_dir.sh, [paths] tmp_dir) and copied to their output
    paths. NEED 2 x BAM + 5 GB; with less free the job writes in
    output/<sample>/bam_sorted/ directly. The symlink branch is unchanged.
    """
    input:
        bam=get_input_bam,
        genome_fai="output/{sample}/genome.fa.fai"
    output:
        bam=temp("output/{sample}/bam_sorted/{bam_id}.sorted.bam"),
        csi=temp("output/{sample}/bam_sorted/{bam_id}.sorted.bam.csi")
    log:
        "logs/{sample}/check_bam_sorted/{bam_id}.log"
    benchmark:
        "benchmarks/{sample}/check_bam_sorted/{bam_id}.txt"
    params:
        tmp_root=TMP_ROOT,
        mode=lambda w: get_braker_mode(w.sample)
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Reference names in the BAM must match the genome FASTA headers
        # (first word, as kept by prepare_genome). Otherwise GeneMark-ETP's
        # StringTie step silently produces an empty transcripts_merged.fasta
        # and fails much later (issue #76).
        bash {script_dir}/check_bam_refs.sh {input.bam} {input.genome_fai} 2> {log} || {{ cat {log} >&2; exit 1; }}

        # Check if BAM is already coordinate-sorted and indexable
        # Even if header says sorted, unmapped reads might be in wrong position
        IS_SORTED=false

        if samtools view -H {input.bam} | grep -q '@HD.*SO:coordinate'; then
            echo "BAM file {input.bam} header indicates coordinate-sorted" >> {log}

            # Try to create a symlink and index it
            ln -sf $(readlink -f {input.bam}) {output.bam} 2>> {log}

            # Test if we can index it
            if samtools index -c -@ {threads} {output.bam} 2>> {log}; then
                echo "BAM file is properly sorted and indexable" >> {log}
                IS_SORTED=true
            else
                echo "BAM file claims to be sorted but cannot be indexed (unmapped reads in wrong position)" >> {log}
                echo "Will re-sort to fix the issue" >> {log}
                rm -f {output.bam} {output.csi}
                IS_SORTED=false
            fi
        else
            echo "BAM file {input.bam} is not coordinate-sorted" >> {log}
            IS_SORTED=false
        fi

        # If not properly sorted, re-sort
        if [ "$IS_SORTED" = "false" ]; then
            echo "Sorting BAM file..." >> {log}

            # Sort chunks and the sorted BAM go to the node-local disk
            source {script_dir}/tmp_dir.sh
            BAM_ABS=$PWD/{output.bam}
            CSI_ABS=$PWD/{output.csi}
            scratch_dir outDir "bamsort_{wildcards.sample}_{wildcards.bam_id}" "{params.tmp_root}" \
                "$(need_gb 2 {input.bam})" "$(dirname "$BAM_ABS")" 2>> {log}
            trap 'rm -rf -- "$SCRATCH"' EXIT
            SORTED_BAM="$outDir/{wildcards.bam_id}.sorted.bam"

            # Sort with parallel threads (fixes BRAKER issue: was -@ 0)
            samtools sort \
                -@ {threads} \
                -T "$outDir/sort" \
                -o "$SORTED_BAM" \
                {input.bam} \
                2>> {log}

            echo "Sorting complete" >> {log}

            # Create index
            samtools index -c -@ {threads} "$SORTED_BAM" 2>> {log}
            echo "Indexing complete" >> {log}
            if [ -n "$SCRATCH" ]; then
                cp "$SORTED_BAM" "$BAM_ABS.tmp"
                mv "$BAM_ABS.tmp" "$BAM_ABS"
                cp "$SORTED_BAM.csi" "$CSI_ABS.tmp"
                mv "$CSI_ABS.tmp" "$CSI_ABS"
            fi
            if [ ! -s "$BAM_ABS" ] || [ ! -s "$CSI_ABS" ]; then
                echo "ERROR: {output.bam} or its .csi missing after the copy back" >> {log}
                exit 1
            fi
        fi

        # A BAM without mapped reads, or a STAR BAM whose spliced reads lack
        # the XS strand tag, passes the checks above, and GeneMark-ETP then
        # crashes hours later in parse_set.pl (Gaius-Augustus/BRAKER#719).
        # GeneMark-ETP runs StringTie, which needs XS; ET mode (bam2hints)
        # does not, so there it is only a warning.
        MAPPED=$(samtools idxstats {output.bam} | awk '{{s += $3}} END {{print s + 0}}')
        if [ "$MAPPED" -eq 0 ]; then
            echo "ERROR: {input.bam} contains no mapped reads. Check the aligner log and that the reads were aligned to this genome." | tee -a {log} >&2
            exit 1
        fi
        echo "Mapped reads: $MAPPED" >> {log}
        set +o pipefail    # awk stops reading after 1M alignments
        read -r N_SPL N_XS < <(samtools view -F 4 {output.bam} | awk -F'\t' '$6 ~ /N/ {{n++; if ($0 ~ /\tXS:A:[+-]/) x++}} NR >= 1000000 {{exit}} END {{print n + 0, x + 0}}')
        set -o pipefail
        echo "Spliced alignments in the first 1M: $N_SPL, with XS tag: $N_XS" >> {log}
        if [ "$N_SPL" -gt 0 ] && [ "$N_XS" -eq 0 ]; then
            MSG="none of the $N_SPL spliced reads in the first 1M alignments of {input.bam} has an XS strand tag. GeneMark-ETP's StringTie step needs it (STAR: --outSAMstrandField intronMotif; HISAT2 writes it by default)"
            if [ "{params.mode}" = "etp" ] || [ "{params.mode}" = "dual" ]; then
                echo "ERROR: $MSG. Realign, or run ET mode (leave protein_fasta empty)." | tee -a {log} >&2
                exit 1
            else
                echo "WARNING: $MSG" >> {log}
            fi
        fi

        echo "Final BAM file: {output.bam}" >> {log}
        echo "Final index file: {output.csi}" >> {log}
        """
