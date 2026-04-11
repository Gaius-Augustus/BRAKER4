rule run_compleasm:
    """
    Run compleasm to generate BUSCO-based hints for AUGUSTUS.

    Uses scripts/compleasm_to_hints.py which:
    1. Runs compleasm on the genome to identify BUSCO genes
    2. Filters for complete single/duplicated BUSCOs without frameshifts
    3. Converts miniprot alignments to CDSpart hints for AUGUSTUS
    4. Reports BUSCO completeness scores

    The --library_path flag points to the shared pre-downloaded lineage
    directory so compleasm doesn't re-download for every scenario.

    Input:
        genome: Genome assembly FASTA file

    Output:
        compleasm_hints: BUSCO-based CDSpart hints in GFF format
        compleasm_log: Log file
        compleasm_summary: BUSCO completeness summary
    """
    input:
        genome=lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        compleasm_hints="output/{sample}/compleasm_hints.gff",
        compleasm_log="output/{sample}/compleasm_to_hints.log",
        compleasm_summary="output/{sample}/compleasm_genome_out/summary.txt"
    log:
        "logs/{sample}/compleasm/compleasm.log"
    benchmark:
        "benchmarks/{sample}/run_compleasm/run_compleasm.txt"
    params:
        busco_lineage=lambda w: get_busco_lineage(w),
        compleasm_outdir=lambda w: f"output/{w.sample}/compleasm_genome_out",
        library_path=config['compleasm_download_path'],
        script=os.path.join(script_dir, "compleasm_to_hints.py")
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH

        echo "[INFO] Running compleasm_to_hints.py" > {log}
        echo "[INFO] Genome: {input.genome}" >> {log}
        echo "[INFO] Lineage: {params.busco_lineage}" >> {log}
        echo "[INFO] Library path: {params.library_path}" >> {log}
        echo "[INFO] Threads: {threads}" >> {log}

        # Ensure shared download directory exists
        mkdir -p {params.library_path}

        # Prevent host's ~/.local Python packages from shadowing container packages
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1

        # Run our compleasm_to_hints.py with --library_path for pre-downloaded data
        python3 {params.script} \
            -g {input.genome} \
            -d {params.busco_lineage} \
            -t {threads} \
            -o {output.compleasm_hints} \
            -s {params.compleasm_outdir} \
            -L {params.library_path} \
            > {output.compleasm_log} 2>&1 || true

        # Ensure summary exists
        # compleasm creates summary.txt inside scratch_dir or a subdirectory
        mkdir -p $(dirname {output.compleasm_summary})
        if [ ! -f {output.compleasm_summary} ]; then
            # Search for summary.txt in the scratch dir tree
            FOUND_SUMMARY=$(find {params.compleasm_outdir} -name "summary.txt" 2>/dev/null | head -1)
            if [ -n "$FOUND_SUMMARY" ] && [ -f "$FOUND_SUMMARY" ]; then
                cp "$FOUND_SUMMARY" {output.compleasm_summary}
            # Fall back to extracting from script stdout (captured in compleasm_log)
            elif grep -q "The following BUSCOs" {output.compleasm_log} 2>/dev/null; then
                grep -A20 "The following BUSCOs" {output.compleasm_log} | tail -n +2 > {output.compleasm_summary}
            elif grep -q "^S:" {output.compleasm_log} 2>/dev/null; then
                grep "^[SDFMN]:" {output.compleasm_log} > {output.compleasm_summary}
            else
                echo "Compleasm did not produce a summary. Check {log} for details." > {output.compleasm_summary}
            fi
        fi

        # Ensure hints file exists (may be empty if no BUSCOs found)
        touch {output.compleasm_hints}

        n_hints=$(wc -l < {output.compleasm_hints})
        echo "[INFO] Generated $n_hints BUSCO-based hints" >> {log}

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        COMPLEASM_VER=$(grep -oP '__version__ = "\K[^"]+' /opt/compleasm_kit/_version.py 2>/dev/null || true)
        MINIPROT_VER=$(miniprot --version 2>&1 | head -1 || true)
        ( flock 9
          printf "compleasm\t%s\n" "$COMPLEASM_VER" >> "$VERSIONS_FILE"
          printf "miniprot\t%s\n" "$MINIPROT_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite compleasm "$REPORT_DIR"
        cite miniprot "$REPORT_DIR"
        """
