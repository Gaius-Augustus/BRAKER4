"""
Compute per-gene evidence support summary.

For each transcript in braker.gtf, checks how many introns and exons
are supported by hints from the hintsfile, split by evidence source
(RNA-Seq vs protein).

Output: TSV with per-transcript support statistics.
Container: teambraker/braker3:latest (contains python3 + intervaltree)

Note: Requires --singularity-args '--env PREPEND_PATH=/opt/conda/bin' to ensure
the container's conda python is used (host /opt/conda may shadow it).
"""


rule gene_support_summary:
    """Compute per-gene evidence support from hints file."""
    input:
        gtf="output/{sample}/braker.gtf",
        hints=lambda w: get_augustus_hintsfile(w)
    output:
        tsv="output/{sample}/gene_support.tsv"
    log:
        "logs/{sample}/gene_support/gene_support.log"
    benchmark:
        "benchmarks/{sample}/gene_support/gene_support.txt"
    params:
        script=os.path.join(script_dir, "gene_support_summary.py")
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
        python3 {params.script} \
            -g {input.gtf} \
            -H {input.hints} \
            -o {output.tsv} \
            2> {log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        N_SUPPORTED=$(awk -F'\\t' 'NR>1 && ($4+$5)>0' {output.tsv} | wc -l || echo 0)
        N_TOTAL=$(awk -F'\\t' 'NR>1' {output.tsv} | wc -l || echo 0)
        """
