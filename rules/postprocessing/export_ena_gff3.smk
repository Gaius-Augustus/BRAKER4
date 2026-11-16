"""
GFF3 for submission to ENA as decoupled annotation (Webin-CLI, -context
annotation, ANALYSIS_TYPE DECOUPLED_ANNOTATION;
https://ena-docs.readthedocs.io/en/latest/submit/gff3.html).

Runs only for samples with ena_locus_tag_prefix in samples.csv. Optional
columns: ena_taxid (NCBI taxonomy ID, written as ##species) and
ena_seqid_map (TSV: BRAKER sequence name, ENA sequence name or accession;
needed when the genome was submitted under other names, or when a name
contains a '.' that is not a version suffix).

Input is braker_with_ncRNA.gff3 when run_ncrna = 1, else braker.gff3.
scripts/gff3_to_ena.py keeps only INSDC qualifiers, adds locus_tag,
product, transl_table and inference, marks partial models from start/stop
codons in the genome, and leaves out what ENA rejects (see the report).

Runs on the host Python (standard library only), no container.
"""


def _ena_gff3_input(wildcards):
    if config.get('run_ncrna', False):
        return f"output/{wildcards.sample}/braker_with_ncRNA.gff3"
    return f"output/{wildcards.sample}/braker.gff3"


def _ena_inputs(wildcards):
    inputs = {
        "gff3": _ena_gff3_input(wildcards),
        "genome": get_genome(wildcards.sample),
    }
    seqid_map = get_ena_settings(wildcards.sample)["seqid_map"]
    if seqid_map:
        inputs["seqid_map"] = seqid_map
    return inputs


rule export_ena_gff3:
    """Write braker.ena.gff3 for ENA decoupled annotation (host Python)."""
    input:
        unpack(_ena_inputs)
    output:
        gff3="output/{sample}/braker.ena.gff3",
        ids="output/{sample}/braker.ena.ids.tsv",
        report="output/{sample}/braker.ena.report.txt"
    log:
        "logs/{sample}/ena/export_ena_gff3.log"
    benchmark:
        "benchmarks/{sample}/ena/export_ena_gff3.txt"
    params:
        prefix=lambda w: get_ena_settings(w.sample)["locus_tag_prefix"],
        taxid=lambda w: (f"--taxid {get_ena_settings(w.sample)['taxid']}"
                         if get_ena_settings(w.sample)["taxid"] else ""),
        seqid_map=lambda w, input: (f"--seqid-map {input.seqid_map}"
                                    if "seqid_map" in input.keys() else ""),
        table=config['translation_table'],
        versions="output/{sample}/software_versions.tsv"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        python3 {script_dir}/gff3_to_ena.py \
            --gff3 {input.gff3} \
            --genome {input.genome} \
            --locus-tag-prefix {params.prefix} \
            {params.taxid} {params.seqid_map} \
            --translation-table {params.table} \
            --software-versions {params.versions} \
            -o {output.gff3} \
            --id-map {output.ids} \
            --report {output.report} \
            > {log} 2>&1
        echo "[INFO] ENA GFF3: $(grep -c -P '\tgene\t' {output.gff3} || true) genes" >> {log}
        """
