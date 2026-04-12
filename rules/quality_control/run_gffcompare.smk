"""
Evaluate gene predictions against a reference annotation using gffcompare.

Compares CDS features of the BRAKER predictions against a reference GTF
to compute sensitivity and specificity at gene, transcript, exon, and
intron levels.

Both the reference and prediction GTFs are filtered to CDS features only
before comparison, ensuring a fair CDS-level evaluation.

Input:
    - BRAKER predictions GTF (braker.gtf)
    - Reference annotation GTF (from samples.csv reference_gtf column)

Output:
    - gffcompare stats file with sensitivity/specificity metrics
    - gffcompare detailed output files

Container: quay.io/biocontainers/gffcompare:0.12.6
"""


rule extract_cds_reference:
    """Extract CDS features from the reference annotation."""
    input:
        ref_gtf=lambda wildcards: get_reference_gtf(wildcards.sample)
    output:
        cds_gtf=temp("output/{sample}/gffcompare/reference.CDS.gtf")
    benchmark:
        "benchmarks/{sample}/extract_cds_reference/extract_cds_reference.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.cds_gtf})
        awk '$3 == "CDS"' {input.ref_gtf} > {output.cds_gtf}
        echo "Extracted $(wc -l < {output.cds_gtf}) CDS features from reference"
        """


rule extract_cds_prediction:
    """Extract CDS features from BRAKER predictions."""
    input:
        pred_gtf="output/{sample}/braker.gtf"
    output:
        cds_gtf=temp("output/{sample}/gffcompare/prediction.CDS.gtf")
    benchmark:
        "benchmarks/{sample}/extract_cds_prediction/extract_cds_prediction.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.cds_gtf})
        awk '$3 == "CDS"' {input.pred_gtf} > {output.cds_gtf}
        echo "Extracted $(wc -l < {output.cds_gtf}) CDS features from predictions"
        """


rule run_gffcompare:
    """Run gffcompare to evaluate predictions against reference."""
    input:
        ref_cds="output/{sample}/gffcompare/reference.CDS.gtf",
        pred_cds="output/{sample}/gffcompare/prediction.CDS.gtf"
    output:
        stats="output/{sample}/gffcompare/gffcompare.stats"
    log:
        "logs/{sample}/gffcompare/gffcompare.log"
    benchmark:
        "benchmarks/{sample}/gffcompare/gffcompare.txt"
    params:
        outprefix=lambda wildcards: f"output/{wildcards.sample}/gffcompare/gffcompare"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        GFFCOMPARE_CONTAINER
    shell:
        r"""
        set -euo pipefail
        gffcompare \
            --strict-match \
            -e 3 \
            -T \
            -r {input.ref_cds} \
            -o {params.outprefix} \
            {input.pred_cds} \
            2>&1 | tee {log}

        echo "" >> {log}
        echo "=== gffcompare summary ===" >> {log}
        cat {output.stats} >> {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        GFFC_VER=$(gffcompare --version 2>&1 | head -1 || echo "unknown")
        ( flock 9; printf "gffcompare\t%s\n" "$GFFC_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite gffcompare "$REPORT_DIR"
        """
