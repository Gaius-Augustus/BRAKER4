"""
minisplice: score canonical splice sites with a small CNN.

Produces a splice-score TSV that minimap2 v2.29+ reads via --spsc to
improve junction detection during IsoSeq alignment. The model was trained
on vertebrate and insect genomes (Heng Li, Algorithms Mol Biol 2026).

Only runs when use_minisplice = 1 AND the sample has unaligned IsoSeq reads
(isoseq_fastq column). Off by default.

Container: quay.io/biocontainers/minisplice:0.4
"""


rule run_minisplice:
    """Score canonical GT/AG splice sites in the genome with minisplice."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
    output:
        scores="output/{sample}/minisplice/splice_scores.tsv",
    log:
        "logs/{sample}/minisplice/minisplice.log"
    benchmark:
        "benchmarks/{sample}/minisplice/minisplice.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        MINISPLICE_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.scores})

        echo "[INFO] Running minisplice predict on genome..." > {log}
        echo "[INFO] Threads: {threads}" >> {log}

        # Model and calibration files are bundled inside the container.
        # Default path: /usr/local/share/minisplice/vi2-7k.kan
        MODEL="${{MINISPLICE_MODEL:-/usr/local/share/minisplice/vi2-7k.kan}}"
        CALI="$MODEL.cali"
        echo "[INFO] Model: $MODEL" >> {log}

        minisplice predict \
            -t{threads} \
            -c "$CALI" \
            "$MODEL" \
            {input.genome} \
            > {output.scores} \
            2>> {log}

        n_sites=$(wc -l < {output.scores})
        echo "[INFO] Scored $n_sites splice sites" >> {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        MS_VER=$(minisplice version 2>&1 || minisplice --version 2>&1 || echo "unknown")
        ( flock 9; printf "minisplice\t%s\n" "$MS_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite minisplice "$REPORT_DIR" || true
        """
