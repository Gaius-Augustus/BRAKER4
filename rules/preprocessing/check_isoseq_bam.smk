"""
Sort and index pre-aligned IsoSeq BAM files.

Supports multiple IsoSeq BAMs (colon-separated in samples.csv).
Each BAM is sorted and indexed independently. If multiple BAMs exist,
they are merged downstream by merge_isoseq_bams.

Container: teambraker/braker3:latest (contains samtools)
"""

def get_input_isoseq_bam_by_id(wildcards):
    """Get the pre-aligned IsoSeq BAM path for a given isoseq_id."""
    bam_files = get_isoseq_bam_files(wildcards.sample)
    bam_ids = get_isoseq_bam_ids(wildcards.sample)
    for bam_path, bid in zip(bam_files, bam_ids):
        if bid == wildcards.isoseq_id:
            return bam_path
    raise ValueError(f"No IsoSeq BAM found for id {wildcards.isoseq_id} in sample {wildcards.sample}")

rule check_isoseq_bam:
    """Sort and index an IsoSeq BAM.

    Scratch: the sorted BAM, its .csi and the sort chunks are written to a
    private directory on the node-local disk (scripts/tmp_dir.sh, [paths]
    tmp_dir) and copied to their output paths. NEED 2 x BAM + 5 GB; with
    less free the job writes in output/<sample>/isoseq_sorted/ directly.
    """
    input:
        bam=get_input_isoseq_bam_by_id
    output:
        bam=temp("output/{sample}/isoseq_sorted/{isoseq_id}.sorted.bam"),
        csi=temp("output/{sample}/isoseq_sorted/{isoseq_id}.sorted.bam.csi")
    log:
        "logs/{sample}/check_isoseq_bam/{isoseq_id}.log"
    benchmark:
        "benchmarks/{sample}/check_isoseq_bam/{isoseq_id}.txt"
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
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})

        BAM_ABS=$(readlink -f {input.bam})

        echo "Sorting IsoSeq BAM file {wildcards.isoseq_id}..." > {log}

        # Sort chunks and the sorted BAM go to the node-local disk
        source {script_dir}/tmp_dir.sh
        OUT_BAM_ABS=$PWD/{output.bam}
        OUT_CSI_ABS=$PWD/{output.csi}
        scratch_dir outDir "bamsort_{wildcards.sample}_{wildcards.isoseq_id}" "{params.tmp_root}" \
            "$(need_gb 2 "$BAM_ABS")" "$(dirname "$OUT_BAM_ABS")" 2>> {log}
        trap 'rm -rf -- "$SCRATCH"' EXIT
        SORTED_BAM="$outDir/{wildcards.isoseq_id}.sorted.bam"

        samtools sort -@ {threads} -T "$outDir/sort" -o "$SORTED_BAM" "$BAM_ABS" 2>> {log}
        echo "Sorting complete" >> {log}

        samtools index -c -@ {threads} "$SORTED_BAM" 2>> {log}
        echo "Indexing complete" >> {log}
        if [ -n "$SCRATCH" ]; then
            cp "$SORTED_BAM" "$OUT_BAM_ABS.tmp"
            mv "$OUT_BAM_ABS.tmp" "$OUT_BAM_ABS"
            cp "$SORTED_BAM.csi" "$OUT_CSI_ABS.tmp"
            mv "$OUT_CSI_ABS.tmp" "$OUT_CSI_ABS"
        fi
        if [ ! -s "$OUT_BAM_ABS" ] || [ ! -s "$OUT_CSI_ABS" ]; then
            echo "ERROR: {output.bam} or its .csi missing after the copy back" >> {log}
            exit 1
        fi
        """


rule merge_isoseq_bams:
    """Merge multiple sorted IsoSeq BAMs into one for GeneMark-ETP."""
    input:
        bams=lambda wildcards: get_isoseq_sorted_bams(wildcards.sample)
    output:
        bam="output/{sample}/isoseq_merged/isoseq.merged.bam",
        csi="output/{sample}/isoseq_merged/isoseq.merged.bam.csi"
    log:
        "logs/{sample}/merge_isoseq_bams/merge.log"
    benchmark:
        "benchmarks/{sample}/merge_isoseq_bams/merge.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam})

        echo "Merging $(echo {input.bams} | wc -w) IsoSeq BAMs..." > {log}
        samtools merge -@ {threads} {output.bam} {input.bams} 2>> {log}
        samtools index -c -@ {threads} {output.bam} 2>> {log}
        echo "Merged IsoSeq BAM: $(samtools view -c {output.bam}) reads" >> {log}

        # Remove per-FASTQ sorted BAMs once merged (can be 1-20 GB each).
        for _bam in {input.bams}; do
            rm -f "$_bam" "$_bam.bai" 2>/dev/null || true
        done
        """
