"""
HISAT2 alignment of RNA-Seq reads to genome.

Builds a HISAT2 index from the genome, then aligns FASTQ reads
(from SRA download or user-provided) to produce sorted BAM files.

Mirrors braker.pl's make_bam_file() function:
- hisat2-build for indexing
- hisat2 --dta for spliced alignment
- samtools sort for coordinate sorting

Input:
    - Genome FASTA (masked)
    - FASTQ files (paired or unpaired)

Output:
    - Coordinate-sorted BAM file with index

Container: teambraker/braker3:latest (contains hisat2, samtools)
"""


rule hisat2_index:
    """Build HISAT2 genome index (once per sample)."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        index="output/{sample}/hisat2/genome.1.ht2"
    log:
        "logs/{sample}/hisat2/hisat2_build.log"
    benchmark:
        "benchmarks/{sample}/hisat2/hisat2_build.txt"
    params:
        prefix=lambda wildcards: f"output/{wildcards.sample}/hisat2/genome"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p output/{wildcards.sample}/hisat2

        hisat2-build \
            -p {threads} \
            {input.genome} \
            {params.prefix} \
            > {log} 2>&1

        if [ ! -f {output.index} ]; then
            echo "ERROR: hisat2-build failed to create index" >> {log}
            exit 1
        fi

        echo "HISAT2 index built successfully" >> {log}
        """


def _get_align_deps(wildcards):
    """Get input dependencies for HISAT2 alignment based on data source."""
    sample = wildcards.sample
    align_id = wildcards.align_id

    if align_id in get_sra_ids(sample):
        # SRA: depend on download completion marker
        return f"output/{sample}/sra_fastq/{align_id}/.download_complete"
    elif align_id in get_fastq_ids(sample):
        # User-provided FASTQ: depend on actual files
        r1 = get_fastq_r1(sample, align_id)
        r2 = get_fastq_r2(sample, align_id)
        return [r1, r2]
    else:
        raise ValueError(f"Unknown alignment ID: {align_id} for sample {sample}")


rule hisat2_align:
    """Align FASTQ reads with HISAT2 and produce sorted BAM."""
    input:
        index="output/{sample}/hisat2/genome.1.ht2",
        deps=_get_align_deps
    output:
        bam="output/{sample}/hisat2_aligned/{align_id}.sorted.bam",
        bai="output/{sample}/hisat2_aligned/{align_id}.sorted.bam.bai"
    log:
        "logs/{sample}/hisat2/{align_id}.log"
    benchmark:
        "benchmarks/{sample}/hisat2/{align_id}.txt"
    params:
        source=lambda wildcards: "sra" if wildcards.align_id in get_sra_ids(wildcards.sample) else "fastq",
        sra_dir=lambda wildcards: f"output/{wildcards.sample}/sra_fastq",
        r1=lambda wildcards: get_fastq_r1(wildcards.sample, wildcards.align_id) if wildcards.align_id in get_fastq_ids(wildcards.sample) else "",
        r2=lambda wildcards: get_fastq_r2(wildcards.sample, wildcards.align_id) if wildcards.align_id in get_fastq_ids(wildcards.sample) else "",
        index_prefix=lambda wildcards: f"output/{wildcards.sample}/hisat2/genome"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "Aligning {wildcards.align_id} (source: {params.source})..." > {log}

        if [ "{params.source}" = "sra" ]; then
            # SRA-derived FASTQs: check for paired vs unpaired
            if [ -f "{params.sra_dir}/{wildcards.align_id}_1.fastq" ] && \
               [ -f "{params.sra_dir}/{wildcards.align_id}_2.fastq" ]; then
                echo "Paired-end SRA alignment" >> {log}
                hisat2 -x {params.index_prefix} \
                    -1 {params.sra_dir}/{wildcards.align_id}_1.fastq \
                    -2 {params.sra_dir}/{wildcards.align_id}_2.fastq \
                    --dta -p {threads} \
                    2>> {log} | \
                    samtools sort -@ {threads} -o {output.bam}
            elif [ -f "{params.sra_dir}/{wildcards.align_id}.fastq" ]; then
                echo "Single-end SRA alignment" >> {log}
                hisat2 -x {params.index_prefix} \
                    -U {params.sra_dir}/{wildcards.align_id}.fastq \
                    --dta -p {threads} \
                    2>> {log} | \
                    samtools sort -@ {threads} -o {output.bam}
            else
                echo "ERROR: No FASTQ files found for SRA ID {wildcards.align_id}" >> {log}
                ls -la {params.sra_dir}/ >> {log} 2>&1 || true
                exit 1
            fi
        else
            # User-provided paired-end FASTQs
            echo "Paired-end FASTQ alignment: {params.r1} {params.r2}" >> {log}
            hisat2 -x {params.index_prefix} \
                -1 {params.r1} \
                -2 {params.r2} \
                --dta -p {threads} \
                2>> {log} | \
                samtools sort -@ {threads} -o {output.bam}
        fi

        # Index the BAM
        samtools index -@ {threads} {output.bam} 2>> {log}

        N_READS=$(samtools view -c {output.bam})
        echo "Alignment complete: $N_READS reads mapped" >> {log}

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        HISAT2_VER=$(hisat2 --version 2>&1 | head -1 | grep -oP 'version \K\S+' || echo "unknown")
        SAM_VER=$(samtools --version 2>&1 | head -1 | awk '{{print $2}}' || echo "unknown")
        ( flock 9
          printf "HISAT2\t%s\n" "$HISAT2_VER" >> "$VERSIONS_FILE"
          printf "SAMtools\t%s\n" "$SAM_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite hisat2 "$REPORT_DIR"
        cite samtools "$REPORT_DIR"
        """
