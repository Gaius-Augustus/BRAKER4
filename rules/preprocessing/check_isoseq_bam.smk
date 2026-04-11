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
    input:
        bam=get_input_isoseq_bam_by_id
    output:
        bam=temp("output/{sample}/isoseq_sorted/{isoseq_id}.sorted.bam"),
        bai=temp("output/{sample}/isoseq_sorted/{isoseq_id}.sorted.bam.bai")
    log:
        "logs/{sample}/check_isoseq_bam/{isoseq_id}.log"
    benchmark:
        "benchmarks/{sample}/check_isoseq_bam/{isoseq_id}.txt"
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
        samtools sort -@ {threads} -o {output.bam} "$BAM_ABS" 2>> {log}
        echo "Sorting complete" >> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}
        echo "Indexing complete" >> {log}
        """


rule merge_isoseq_bams:
    """Merge multiple sorted IsoSeq BAMs into one for GeneMark-ETP."""
    input:
        bams=lambda wildcards: get_isoseq_sorted_bams(wildcards.sample)
    output:
        bam="output/{sample}/isoseq_merged/isoseq.merged.bam",
        bai="output/{sample}/isoseq_merged/isoseq.merged.bam.bai"
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
        samtools index -@ {threads} {output.bam} 2>> {log}
        echo "Merged IsoSeq BAM: $(samtools view -c {output.bam}) reads" >> {log}
        """
