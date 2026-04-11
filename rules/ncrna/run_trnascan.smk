"""
tRNAscan-SE: tRNA gene prediction.

Scans the genome for transfer RNA genes using tRNAscan-SE in eukaryotic mode.
Outputs GFF3 format natively via the --gff flag.

Container: quay.io/biocontainers/trnascan-se:2.0.12--pl5321h031d066_0
"""

TRNASCAN_CONTAINER = config.get(
    "trnascan_image",
    "docker://quay.io/biocontainers/trnascan-se:2.0.12--pl5321h031d066_0"
)


rule run_trnascan:
    """Run tRNAscan-SE to predict tRNA genes on the genome."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        gff="output/{sample}/ncrna/tRNAs.gff3",
        txt="output/{sample}/ncrna/tRNAs.txt"
    log:
        "logs/{sample}/ncrna/trnascan.log"
    benchmark:
        "benchmarks/{sample}/ncrna/trnascan.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        TRNASCAN_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gff})

        echo "[INFO] Running tRNAscan-SE in eukaryotic mode..." > {log}

        LC_ALL=C tRNAscan-SE \
            -E \
            --thread {threads} \
            -q \
            --forceow \
            -o {output.txt} \
            --gff {output.gff}.tmp \
            {input.genome} \
            2>> {log} || true

        # Ensure output exists even if no tRNAs found
        if [ -s {output.gff}.tmp ] && grep -qv '^#' {output.gff}.tmp; then
            # Prefix tRNA IDs with sample name for uniqueness
            awk -F'\t' -v OFS='\t' -v p="{wildcards.sample}" '
                BEGIN {{n=1}}
                /^#/ {{print; next}}
                {{
                    if ($9 ~ /ID=/) {{
                        gsub(/ID=[^;]+/, "ID=" p "-tRNA_" n, $9)
                    }} else {{
                        $9 = "ID=" p "-tRNA_" n ";" $9
                    }}
                    n++
                    print
                }}
            ' {output.gff}.tmp > {output.gff}
        else
            echo "##gff-version 3" > {output.gff}
            echo "[INFO] No tRNA genes found" >> {log}
        fi
        rm -f {output.gff}.tmp

        n_trna=$(grep -cv '^#' {output.gff} || echo 0)
        echo "[INFO] tRNAscan-SE predicted $n_trna tRNA features" >> {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        # BusyBox grep in this container has no -P; use sed/awk for portable extraction
        TRNASCAN_VER=$(LC_ALL=C tRNAscan-SE --help 2>&1 | sed -n 's/.*tRNAscan-SE \([0-9.][0-9.]*\).*/\1/p' | head -1 || echo "unknown")
        ( flock 9; printf "tRNAscan-SE\t%s\n" "$TRNASCAN_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite trnascan "$REPORT_DIR" || true
        """
