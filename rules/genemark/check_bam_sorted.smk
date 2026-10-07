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
        tmp_root=TMP_ROOT
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
        samtools view -H {input.bam} | awk -F'\t' -v fai={input.genome_fai} '
            BEGIN {{ while ((getline l < fai) > 0) {{ split(l, f, "\t"); len[f[1]] = f[2] }} }}
            $1 == "@SQ" {{
                sn = ""; ln = ""
                for (i = 2; i <= NF; i++) {{
                    if ($i ~ /^SN:/) sn = substr($i, 4)
                    if ($i ~ /^LN:/) ln = substr($i, 4)
                }}
                n++
                if (sn in len) {{ ok++; if (len[sn] != ln) diff++ }}
                else if (miss++ < 5) ex = ex " " sn
            }}
            END {{
                if (n > 0 && ok == 0) {{
                    print "ERROR: none of the " n " reference sequences in {input.bam} match the genome FASTA headers." > "/dev/stderr"
                    print "  BAM references, e.g.:" ex > "/dev/stderr"
                    print "  Align the reads against the same genome FASTA (identical sequence names) that you pass to BRAKER4." > "/dev/stderr"
                    exit 1
                }}
                if (miss > 0) printf "WARNING: %d of %d BAM reference sequences are not in the genome, e.g.:%s\n", miss, n, ex > "/dev/stderr"
                if (diff > 0) printf "WARNING: %d BAM reference sequences have a different length than in the genome. Was the BAM made against another assembly version?\n", diff > "/dev/stderr"
            }}' 2> {log} || {{ cat {log} >&2; exit 1; }}

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

        echo "Final BAM file: {output.bam}" >> {log}
        echo "Final index file: {output.csi}" >> {log}
        """
