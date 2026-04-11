"""
Infernal/cmscan: ncRNA prediction using Rfam covariance models.

Scans the genome against the Rfam database to identify snoRNAs, snRNAs,
miRNAs, ribozymes, and other structured ncRNAs.

Split into two rules:
  1. run_cmscan: runs cmscan in Infernal container → tblout
  2. convert_infernal_to_gff3: converts tblout → GFF3 (no container, uses host Python)

The Rfam database (Rfam.cm + Rfam.clanin) must be pre-downloaded and
cmpress'd. The path is configured via rfam_path in config.ini [paths].

Container: quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_2
"""

INFERNAL_CONTAINER = config.get(
    "infernal_image",
    "docker://quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_2"
)


rule run_cmscan:
    """Run Infernal cmscan against Rfam to predict ncRNA genes."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
        rfam_cm=os.path.join(config.get('rfam_path', 'shared_data/rfam'), 'Rfam.cm'),
        rfam_clanin=os.path.join(config.get('rfam_path', 'shared_data/rfam'), 'Rfam.clanin')
    output:
        tblout="output/{sample}/ncrna/infernal.tblout"
    log:
        "logs/{sample}/ncrna/infernal.log"
    benchmark:
        "benchmarks/{sample}/ncrna/cmscan.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        INFERNAL_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tblout})

        RFAM_CM=$(readlink -f {input.rfam_cm})
        RFAM_CLANIN=$(readlink -f {input.rfam_clanin})
        GENOME_ABS=$(readlink -f {input.genome})

        echo "[INFO] Running Infernal cmscan against Rfam..." > {log}
        echo "[INFO] Rfam CM: $RFAM_CM" >> {log}
        echo "[INFO] Threads: {threads}" >> {log}

        cmscan \
            --cut_ga \
            --rfam \
            --nohmmonly \
            --clanin $RFAM_CLANIN \
            --fmt 2 \
            --cpu {threads} \
            --tblout {output.tblout} \
            $RFAM_CM \
            $GENOME_ABS \
            > /dev/null \
            2>> {log}

        n_hits=$(grep -cv '^#' {output.tblout} || echo 0)
        echo "[INFO] cmscan found $n_hits raw hits" >> {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        INFERNAL_VER=$(cmscan -h 2>&1 | head -2 | tail -1 | sed 's/.*INFERNAL //' | sed 's/ .*//' || echo "unknown")
        ( flock 9; printf "Infernal\t%s\n" "$INFERNAL_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"
        """


rule convert_infernal_to_gff3:
    """Convert Infernal tblout to GFF3 format (runs on host, no container)."""
    input:
        tblout="output/{sample}/ncrna/infernal.tblout"
    output:
        gff="output/{sample}/ncrna/ncRNAs_infernal.gff3"
    log:
        "logs/{sample}/ncrna/infernal_to_gff3.log"
    benchmark:
        "benchmarks/{sample}/ncrna/infernal_to_gff3.txt"
    params:
        sample="{sample}"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        r"""
        set -euo pipefail
        python3 {script_dir}/infernal_to_gff3.py \
            -i {input.tblout} \
            -o {output.gff} \
            -p {params.sample} \
            2> {log}

        n_ncrna=$(grep -cv '^#' {output.gff} || echo 0)
        echo "[INFO] Converted $n_ncrna ncRNA features to GFF3" >> {log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite infernal "$REPORT_DIR" || true
        cite rfam "$REPORT_DIR" || true
        """
