"""
minimap2 alignment of IsoSeq long reads to genome.

Supports multiple IsoSeq FASTQ/FASTA files (colon-separated in samples.csv).
Each file is aligned independently, then merged downstream if needed.

When use_minisplice=1, splice site scores from minisplice are passed to
minimap2 via --spsc to improve junction detection (requires minimap2 >= 2.29).

Split into two rules per file:
1. minimap2_isoseq_align: Align with minimap2 -> SAM (minimap2 container)
2. sort_isoseq_sam: Convert SAM -> sorted BAM + index (braker3 container)
"""

USE_MINISPLICE = config.get('use_minisplice', False)


def get_isoseq_fastq_by_id(wildcards):
    """Get the IsoSeq FASTQ/FASTA path for a given isoseq_fastq_id."""
    files = get_isoseq_fastq_files(wildcards.sample)
    ids = get_isoseq_fastq_ids(wildcards.sample)
    for fpath, fid in zip(files, ids):
        if fid == wildcards.isoseq_fastq_id:
            return fpath
    raise ValueError(f"No IsoSeq FASTQ found for id {wildcards.isoseq_fastq_id} in sample {wildcards.sample}")


def _minimap2_inputs(wildcards):
    """Return minimap2 inputs, optionally including minisplice scores."""
    inputs = {
        "genome": get_masked_genome(wildcards.sample),
        "isoseq_reads": get_isoseq_fastq_by_id(wildcards),
    }
    if USE_MINISPLICE:
        inputs["splice_scores"] = f"output/{wildcards.sample}/minisplice/splice_scores.tsv"
    return inputs


rule minimap2_isoseq_align:
    """Align IsoSeq reads with minimap2, output SAM."""
    input:
        unpack(_minimap2_inputs)
    output:
        sam=temp("output/{sample}/minimap2_aligned/{isoseq_fastq_id}.sam")
    log:
        "logs/{sample}/minimap2/{isoseq_fastq_id}_align.log"
    benchmark:
        "benchmarks/{sample}/minimap2/{isoseq_fastq_id}_align.txt"
    params:
        spsc_flag=lambda wildcards, input: f"--spsc={input.splice_scores}" if hasattr(input, 'splice_scores') else "",
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        MINIMAP2_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p output/{wildcards.sample}/minimap2_aligned
        mkdir -p $(dirname {log})

        echo "Aligning IsoSeq reads {wildcards.isoseq_fastq_id} with minimap2 splice:hq preset..." > {log}
        if [ -n "{params.spsc_flag}" ]; then
            echo "Using minisplice splice scores: {params.spsc_flag}" >> {log}
        fi

        minimap2 -ax splice:hq -uf \
            -t {threads} \
            {params.spsc_flag} \
            {input.genome} \
            {input.isoseq_reads} \
            > {output.sam} \
            2>> {log}

        echo "minimap2 alignment complete" >> {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        MM2_VER=$(minimap2 --version 2>&1 || echo "unknown")
        ( flock 9; printf "minimap2\t%s\n" "$MM2_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite minimap2 "$REPORT_DIR"
        """


rule sort_isoseq_sam:
    """Convert SAM to sorted BAM and index."""
    input:
        sam="output/{sample}/minimap2_aligned/{isoseq_fastq_id}.sam"
    output:
        bam="output/{sample}/minimap2_aligned/{isoseq_fastq_id}.sorted.bam",
        bai="output/{sample}/minimap2_aligned/{isoseq_fastq_id}.sorted.bam.bai"
    log:
        "logs/{sample}/minimap2/{isoseq_fastq_id}_sort.log"
    benchmark:
        "benchmarks/{sample}/minimap2/{isoseq_fastq_id}_sort.txt"
    params:
        sort_threads=lambda wildcards, threads: max(1, threads - 1)
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        echo "Converting SAM to sorted BAM..." > {log}

        samtools view -bS --threads 1 {input.sam} | \
            samtools sort -@ {params.sort_threads} -o {output.bam} 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}

        echo "IsoSeq BAM: $(samtools view -c {output.bam}) reads" >> {log}
        """
