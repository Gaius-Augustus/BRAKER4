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
    """Download an SRA run and convert it to FASTQ.

    Scratch: prefetch's .sra and fastq-dump's FASTQ files are written to a
    private directory on the node-local disk (scripts/tmp_dir.sh, [paths]
    tmp_dir); the FASTQ files are copied to output/<sample>/sra_fastq/.
    NEED 100 GB (the size is unknown before prefetch); with less free the
    job works in output/<sample>/sra_fastq/ as before.
    """
    output:
        marker="output/{sample}/sra_fastq/{sra_id}/.download_complete"
    log:
        "logs/{sample}/download_sra/{sra_id}.log"
    benchmark:
        "benchmarks/{sample}/download_sra/{sra_id}.txt"
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/sra_fastq",
        tmp_root=TMP_ROOT
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

        # .sra and FASTQ files go to the node-local disk
        source {script_dir}/tmp_dir.sh
        finalDir=$PWD/{params.outdir}
        MARKER_ABS=$PWD/{output.marker}
        scratch_dir outDir "sra_{wildcards.sample}_{wildcards.sra_id}" "{params.tmp_root}" 100 \
            "$finalDir" 2>> {log}
        trap 'rm -rf -- "$SCRATCH"' EXIT

        # Step 1: prefetch the .sra file
        prefetch \
            --max-size 35G \
            {wildcards.sra_id} \
            --output-directory "$outDir" \
            >> {log} 2>&1

        # Verify download
        if [ ! -f "$outDir/{wildcards.sra_id}/{wildcards.sra_id}.sra" ]; then
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
            "$outDir/{wildcards.sra_id}/{wildcards.sra_id}.sra" \
            --outdir "$outDir" \
            >> {log} 2>&1

        # Verify FASTQ output
        if [ -f "$outDir/{wildcards.sra_id}_1.fastq" ] && \
           [ -f "$outDir/{wildcards.sra_id}_2.fastq" ]; then
            echo "Paired-end FASTQ files created" >> {log}
            FASTQS="{wildcards.sra_id}_1.fastq {wildcards.sra_id}_2.fastq"
        elif [ -f "$outDir/{wildcards.sra_id}.fastq" ]; then
            echo "Single-end FASTQ file created" >> {log}
            FASTQS="{wildcards.sra_id}.fastq"
        else
            echo "ERROR: fastq-dump produced no FASTQ files" >> {log}
            exit 1
        fi

        # Copy the FASTQ files to the run directory (cp to .tmp, then mv)
        for fq in $FASTQS; do
            if [ -n "$SCRATCH" ]; then
                cp "$outDir/$fq" "$finalDir/$fq.tmp"
                mv "$finalDir/$fq.tmp" "$finalDir/$fq"
            fi
            if [ ! -s "$finalDir/$fq" ]; then
                echo "ERROR: $finalDir/$fq missing after the copy back" >> {log}
                exit 1
            fi
        done

        # Create marker before cleanup (marker is inside the SRA subdir)
        mkdir -p $(dirname "$MARKER_ABS")
        touch "$MARKER_ABS"

        # Clean up .sra file to save space (keep marker dir). Fallback only;
        # on scratch the EXIT trap removes it.
        if [ -z "$SCRATCH" ]; then
            rm -f "$outDir/{wildcards.sra_id}/{wildcards.sra_id}.sra"
        fi
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
