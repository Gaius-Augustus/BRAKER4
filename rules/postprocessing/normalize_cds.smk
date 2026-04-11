"""
Normalize CDS boundaries and validate gene structures.

Fixes BRAKER issues:
- #833: Stop codon included in CDS (GeneMark convention vs AUGUSTUS convention)
- #904: Internal stop codons in protein-coding genes
- #283: Reports non-ATG start codons (warning only, does not discard)

This rule runs after filter_internal_stop_codons and produces the final braker.gtf.
Genes with broken structures after stop codon trimming are discarded entirely.

Container: teambraker/braker3:latest (has BioPython)
"""


rule normalize_cds:
    """Normalize CDS boundaries: trim stop codons, validate structures, report start codons."""
    input:
        gtf="output/{sample}/braker.filtered.gtf",
        genome=lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        gtf="output/{sample}/braker.gtf",
        log_file="output/{sample}/normalize_cds.log"
    log:
        "logs/{sample}/normalize_cds/normalize.log"
    benchmark:
        "benchmarks/{sample}/normalize_cds/normalize.txt"
    params:
        script=os.path.join(script_dir, "normalize_gtf.py")
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
            -g {input.genome} \
            -f {input.gtf} \
            -o {output.gtf} \
            -l {output.log_file} \
            2> {log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        N_GENES_OUT=$(grep -cP '\\tgene\\t' {output.gtf} || echo 0)
        """
