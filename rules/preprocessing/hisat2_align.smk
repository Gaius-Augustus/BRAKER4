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
        sentinel=touch("output/{sample}/hisat2/.index_complete")
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

        # Large genomes (>4 GB) produce .ht2l instead of .ht2 — accept either
        if [ ! -f {params.prefix}.1.ht2 ] && [ ! -f {params.prefix}.1.ht2l ]; then
            echo "ERROR: hisat2-build failed to create index" >> {log}
            exit 1
        fi

        echo "HISAT2 index built successfully" >> {log}
        touch {output.sentinel}
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
    """Align FASTQ reads with HISAT2 and produce sorted BAM.

    Scratch: the sorted BAM, its .csi and the sort chunks are written to a
    private directory on the node-local disk (scripts/tmp_dir.sh, [paths]
    tmp_dir); BAM and .csi are then copied to their output paths. NEED
    2 x the FASTQ files + 5 GB; with less free the job writes in
    output/<sample>/hisat2_aligned/ directly.
    """
    input:
        index="output/{sample}/hisat2/.index_complete",
        deps=_get_align_deps
    output:
        bam="output/{sample}/hisat2_aligned/{align_id}.sorted.bam",
        csi="output/{sample}/hisat2_aligned/{align_id}.sorted.bam.csi"
    log:
        "logs/{sample}/hisat2/{align_id}.log"
    benchmark:
        "benchmarks/{sample}/hisat2/{align_id}.txt"
    params:
        source=lambda wildcards: "sra" if wildcards.align_id in get_sra_ids(wildcards.sample) else "fastq",
        sra_dir=lambda wildcards: f"output/{wildcards.sample}/sra_fastq",
        r1=lambda wildcards: get_fastq_r1(wildcards.sample, wildcards.align_id) if wildcards.align_id in get_fastq_ids(wildcards.sample) else "",
        r2=lambda wildcards: get_fastq_r2(wildcards.sample, wildcards.align_id) if wildcards.align_id in get_fastq_ids(wildcards.sample) else "",
        index_prefix=lambda wildcards: f"output/{wildcards.sample}/hisat2/genome",
        hisat2_threads=lambda wildcards, threads: max(1, threads // 2),
        sort_threads=lambda wildcards, threads: max(1, threads - threads // 2),
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

        echo "Aligning {wildcards.align_id} (source: {params.source})..." > {log}

        # Sort chunks and the sorted BAM go to the node-local disk
        source {script_dir}/tmp_dir.sh
        BAM_ABS=$PWD/{output.bam}
        CSI_ABS=$PWD/{output.csi}
        scratch_dir outDir "hisat2_{wildcards.sample}_{wildcards.align_id}" "{params.tmp_root}" \
            "$(need_gb 2 {params.r1} {params.r2} {params.sra_dir}/{wildcards.align_id}_1.fastq {params.sra_dir}/{wildcards.align_id}_2.fastq {params.sra_dir}/{wildcards.align_id}.fastq)" \
            "$(dirname "$BAM_ABS")" 2>> {log}
        trap 'rm -rf -- "$SCRATCH"' EXIT
        SORTED_BAM="$outDir/{wildcards.align_id}.sorted.bam"

        if [ "{params.source}" = "sra" ]; then
            # SRA-derived FASTQs: check for paired vs unpaired
            if [ -f "{params.sra_dir}/{wildcards.align_id}_1.fastq" ] && \
               [ -f "{params.sra_dir}/{wildcards.align_id}_2.fastq" ]; then
                echo "Paired-end SRA alignment" >> {log}
                hisat2 -x {params.index_prefix} \
                    -1 {params.sra_dir}/{wildcards.align_id}_1.fastq \
                    -2 {params.sra_dir}/{wildcards.align_id}_2.fastq \
                    --dta -p {params.hisat2_threads} \
                    2>> {log} | \
                    samtools sort -@ {params.sort_threads} -T "$outDir/sort" -o "$SORTED_BAM"
            elif [ -f "{params.sra_dir}/{wildcards.align_id}.fastq" ]; then
                echo "Single-end SRA alignment" >> {log}
                hisat2 -x {params.index_prefix} \
                    -U {params.sra_dir}/{wildcards.align_id}.fastq \
                    --dta -p {params.hisat2_threads} \
                    2>> {log} | \
                    samtools sort -@ {params.sort_threads} -T "$outDir/sort" -o "$SORTED_BAM"
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
                --dta -p {params.hisat2_threads} \
                2>> {log} | \
                samtools sort -@ {params.sort_threads} -T "$outDir/sort" -o "$SORTED_BAM"
        fi

        # Index the BAM
        samtools index -c -@ {threads} "$SORTED_BAM" 2>> {log}

        N_READS=$(samtools view -c "$SORTED_BAM")
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
        echo "Alignment complete: $N_READS reads mapped" >> {log}

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        HISAT2_VER=$(hisat2 --version 2>&1 | head -1 | grep -oP 'version \K\S+' || echo "unknown")
        SAM_VER=$(samtools --version 2>&1 | head -1 | awk '{{print $2}}' || echo "unknown")
        ( flock 9
          printf "HISAT2\t%s\n" "$HISAT2_VER" >> "$VERSIONS_FILE"
          printf "SAMtools\t%s\n" "$SAM_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Remove SRA FASTQ files after alignment (no longer needed downstream)
        if [ "{params.source}" = "sra" ]; then
            rm -f {params.sra_dir}/{wildcards.align_id}_1.fastq \
                  {params.sra_dir}/{wildcards.align_id}_2.fastq \
                  {params.sra_dir}/{wildcards.align_id}.fastq 2>/dev/null || true
            echo "Removed SRA FASTQ files for {wildcards.align_id}" >> {log}
        fi

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite hisat2 "$REPORT_DIR"
        cite samtools "$REPORT_DIR"
        """
