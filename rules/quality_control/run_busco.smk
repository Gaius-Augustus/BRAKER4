"""
BUSCO completeness assessment of genome and predicted proteome.

Runs BUSCO twice:
1. Genome mode: assess genome assembly completeness
2. Protein mode: assess predicted protein set completeness

Then generates a combined summary text and comparison bar plot.

Container: ezlabgva/busco:v6.0.0_cv1
"""


rule busco_genome:
    """Run BUSCO on the genome assembly."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        done="output/{sample}/busco/genome/.done"
    log:
        "logs/{sample}/busco/busco_genome.log"
    benchmark:
        "benchmarks/{sample}/busco/busco_genome.txt"
    params:
        busco_lineage=lambda w: get_busco_lineage(w),
        outdir=lambda w: f"output/{w.sample}/busco",
        download_path=config['busco_download_path']
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BUSCO_CONTAINER
    shell:
        r"""
        set -euo pipefail
        OUTDIR_ABS=$(readlink -f {params.outdir})
        rm -rf "$OUTDIR_ABS/genome"

        mkdir -p {params.download_path}
        OFFLINE_FLAG=""
        # Check for dataset.cfg, not just the directory: compleasm and BUSCO
        # share the same {params.download_path}/lineages/ default, but
        # `compleasm.py download` writes a compleasm-only layout (no
        # dataset.cfg). Only switch to --offline when the directory is a
        # valid BUSCO lineage.
        if [ -f "{params.download_path}/lineages/{params.busco_lineage}/dataset.cfg" ]; then
            OFFLINE_FLAG="--offline"
            echo "[INFO] Found pre-downloaded lineage at {params.download_path}/lineages/{params.busco_lineage}; running BUSCO with --offline" >> {log}
        fi
        busco \
            -i {input.genome} \
            -o genome \
            --out_path "$OUTDIR_ABS" \
            -l {params.busco_lineage} \
            -m genome \
            -c {threads} \
            --download_path {params.download_path} \
            $OFFLINE_FLAG \
            >> {log} 2>&1

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        BUSCO_VER=$(busco --version 2>&1 | head -1 || echo "unknown")
        ( flock 9; printf "BUSCO\t%s\n" "$BUSCO_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        touch {output.done}
        """


rule get_longest_isoform:
    """Extract longest coding isoform per gene locus for BUSCO assessment.

    Uses TSEBRA's get_longest_isoform.py to select one representative
    transcript per gene (the one with the longest CDS). This avoids
    inflated BUSCO duplicate counts from alternative splicing isoforms.

    Then re-extracts proteins from the longest-isoform GTF.
    """
    input:
        gtf="output/{sample}/braker.gtf",
        genome=lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        gtf="output/{sample}/braker.longest.gtf",
        aa="output/{sample}/braker.longest.aa"
    log:
        "logs/{sample}/longest_isoform/longest_isoform.log"
    benchmark:
        "benchmarks/{sample}/get_longest_isoform/get_longest_isoform.txt"
    params:
        translation_table = config.get("translation_table", 1)
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1

        echo "Extracting longest isoform per gene locus..." > {log}

        /opt/TSEBRA/bin/get_longest_isoform.py \
            -g {input.gtf} \
            -o {output.gtf} \
            -q \
            2>> {log}

        n_genes=$(grep -cP '\tgene\t' {output.gtf} || echo 0)
        echo "Longest isoform GTF: $n_genes genes (1 transcript each)" >> {log}

        # Extract proteins from longest-isoform GTF
        STEM=$(echo {output.aa} | sed 's/\.aa$//')
        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {output.gtf} \
            -o $STEM \
            -t {params.translation_table} \
            2>> {log}

        n_prot=$(grep -c "^>" {output.aa} || echo 0)
        echo "Extracted $n_prot proteins" >> {log}
        """


rule busco_proteins:
    """Run BUSCO on the predicted proteome (longest isoforms only)."""
    input:
        proteins="output/{sample}/braker.longest.aa"
    output:
        done="output/{sample}/busco/proteins/.done"
    log:
        "logs/{sample}/busco/busco_proteins.log"
    benchmark:
        "benchmarks/{sample}/busco/busco_proteins.txt"
    params:
        busco_lineage=lambda w: get_busco_lineage(w),
        outdir=lambda w: f"output/{w.sample}/busco",
        download_path=config['busco_download_path']
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BUSCO_CONTAINER
    shell:
        r"""
        set -euo pipefail
        OUTDIR_ABS=$(readlink -f {params.outdir})
        rm -rf "$OUTDIR_ABS/proteins"

        mkdir -p {params.download_path}
        OFFLINE_FLAG=""
        # Check for dataset.cfg, not just the directory: compleasm and BUSCO
        # share the same {params.download_path}/lineages/ default, but
        # `compleasm.py download` writes a compleasm-only layout (no
        # dataset.cfg). Only switch to --offline when the directory is a
        # valid BUSCO lineage.
        if [ -f "{params.download_path}/lineages/{params.busco_lineage}/dataset.cfg" ]; then
            OFFLINE_FLAG="--offline"
            echo "[INFO] Found pre-downloaded lineage at {params.download_path}/lineages/{params.busco_lineage}; running BUSCO with --offline" >> {log}
        fi
        busco \
            -i {input.proteins} \
            -o proteins \
            --out_path "$OUTDIR_ABS" \
            -l {params.busco_lineage} \
            -m proteins \
            -c {threads} \
            --download_path {params.download_path} \
            $OFFLINE_FLAG \
            >> {log} 2>&1

        touch {output.done}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        PROT_SUMMARY=$(find "$OUTDIR_ABS/proteins" -name "short_summary*.txt" 2>/dev/null | head -1)
        BUSCO_LINE=""
        if [ -n "$PROT_SUMMARY" ] && [ -f "$PROT_SUMMARY" ]; then
            BUSCO_LINE=$(grep -oP 'C:[0-9.]+%.*' "$PROT_SUMMARY" | head -1)
        fi
        cite busco "$REPORT_DIR"
        """


rule busco_summary:
    """Generate combined BUSCO summary text.

    The BUSCO v6 container no longer ships `generate_plot.py` / R / ggplot2,
    so we do not generate BUSCO's own figure here. The pipeline's
    `completeness_plot.py` (invoked from collect_results) produces the
    combined BUSCO+compleasm figure that the HTML report consumes.
    """
    input:
        genome_done="output/{sample}/busco/genome/.done",
        proteins_done="output/{sample}/busco/proteins/.done"
    output:
        summary="output/{sample}/busco/busco_summary.txt"
    log:
        "logs/{sample}/busco/busco_summary.log"
    benchmark:
        "benchmarks/{sample}/busco_summary/busco_summary.txt"
    params:
        busco_dir=lambda w: f"output/{w.sample}/busco",
        busco_lineage=lambda w: get_busco_lineage(w)
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BUSCO_CONTAINER
    shell:
        r"""
        set -euo pipefail
        GENOME_SUMMARY=$(find {params.busco_dir}/genome -name "short_summary*.txt" | head -1)
        PROTEINS_SUMMARY=$(find {params.busco_dir}/proteins -name "short_summary*.txt" | head -1)

        # Write combined summary
        echo "===== BUSCO Completeness Summary =====" > {output.summary}
        echo "Lineage: {params.busco_lineage}" >> {output.summary}
        echo "" >> {output.summary}

        echo "--- Genome Assembly ---" >> {output.summary}
        if [ -n "$GENOME_SUMMARY" ] && [ -f "$GENOME_SUMMARY" ]; then
            grep -E "C:|S:|D:|F:|M:|n:" "$GENOME_SUMMARY" >> {output.summary}
        else
            echo "Results not found" >> {output.summary}
        fi
        echo "" >> {output.summary}

        echo "--- Predicted Proteome ---" >> {output.summary}
        if [ -n "$PROTEINS_SUMMARY" ] && [ -f "$PROTEINS_SUMMARY" ]; then
            grep -E "C:|S:|D:|F:|M:|n:" "$PROTEINS_SUMMARY" >> {output.summary}
        else
            echo "Results not found" >> {output.summary}
        fi
        """
