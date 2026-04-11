"""
Download RNA-Seq data from NCBI SRA.

Uses SRA Toolkit (prefetch + fastq-dump) to download and convert
SRA accessions to FASTQ files, matching braker.pl's download_rna_libs().

prefetch downloads the .sra file, then fastq-dump --split-3 converts it
to paired (_1.fastq/_2.fastq) or unpaired (.fastq) FASTQ files.

Input:
    - SRA accession ID (from samples.csv sra_ids column)

Output:
    - Marker file indicating download is complete
    - FASTQ files in output/{sample}/sra_fastq/

Container: teambraker/braker3:latest (contains SRA Toolkit)
"""

rule download_sra:
    output:
        marker="output/{sample}/sra_fastq/{sra_id}/.download_complete"
    log:
        "logs/{sample}/download_sra/{sra_id}.log"
    benchmark:
        "benchmarks/{sample}/download_sra/{sra_id}.txt"
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/sra_fastq"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}

        echo "Downloading SRA accession {wildcards.sra_id}..." > {log}

        # Step 1: prefetch the .sra file
        prefetch \
            --max-size 35G \
            {wildcards.sra_id} \
            --output-directory {params.outdir} \
            >> {log} 2>&1

        # Verify download
        if [ ! -f "{params.outdir}/{wildcards.sra_id}/{wildcards.sra_id}.sra" ]; then
            echo "ERROR: prefetch failed - .sra file not found" >> {log}
            exit 1
        fi

        echo "prefetch complete, converting to FASTQ..." >> {log}

        # Step 2: fastq-dump --split-3 to produce FASTQ files
        # --split-3 produces _1.fastq/_2.fastq for paired, .fastq for unpaired
        # --force: overwrite existing FASTQ files from previous runs
        rm -f {params.outdir}/{wildcards.sra_id}_1.fastq {params.outdir}/{wildcards.sra_id}_2.fastq {params.outdir}/{wildcards.sra_id}.fastq
        fastq-dump \
            --split-3 \
            {params.outdir}/{wildcards.sra_id}/{wildcards.sra_id}.sra \
            --outdir {params.outdir} \
            >> {log} 2>&1

        # Verify FASTQ output
        if [ -f "{params.outdir}/{wildcards.sra_id}_1.fastq" ] && \
           [ -f "{params.outdir}/{wildcards.sra_id}_2.fastq" ]; then
            echo "Paired-end FASTQ files created" >> {log}
        elif [ -f "{params.outdir}/{wildcards.sra_id}.fastq" ]; then
            echo "Single-end FASTQ file created" >> {log}
        else
            echo "ERROR: fastq-dump produced no FASTQ files" >> {log}
            exit 1
        fi

        # Create marker before cleanup (marker is inside the SRA subdir)
        mkdir -p $(dirname {output.marker})
        touch {output.marker}

        # Clean up .sra file to save space (keep marker dir)
        rm -f {params.outdir}/{wildcards.sra_id}/{wildcards.sra_id}.sra
        echo "SRA download complete for {wildcards.sra_id}" >> {log}

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        SRATOOLKIT_VER=$(prefetch --version 2>&1 | grep -oP '[\d.]+' | head -1 || true)
        ( flock 9; printf "SRA Toolkit\t%s\n" "$SRATOOLKIT_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite sratoolkit "$REPORT_DIR"
        """
