"""
Genome repeat masking with Red (REpeat Detector).

Red uses a neural-network approach to detect repeats directly from the genome
sequence, without building a repeat library first. This makes it dramatically
faster than RepeatModeler+RepeatMasker, at the cost of not classifying repeat
families.

Selected when masking_tool = red in config.ini (default: repeatmasker).

Container: quay.io/biocontainers/red:2018.09.10--h9948957_3 (~13 MB)
"""


rule run_red_masking:
    """Soft-mask repeats in the genome using Red (REpeat Detector)."""
    input:
        genome=lambda wildcards: get_genome(wildcards.sample),
    output:
        masked_genome="output/{sample}/preprocessing/genome.fa.masked",
        marker="output/{sample}/preprocessing/.masking_complete",
    log:
        "logs/{sample}/masking/masking.log"
    benchmark:
        "benchmarks/{sample}/masking/masking.txt"
    params:
        sample="{sample}",
        tmp_root=TMP_ROOT
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime']),
    container:
        RED_CONTAINER
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {log})
        LOG_ABS=$(readlink -f {log})
        GENOME_ABS=$(readlink -f {input.genome})
        OUTDIR_ABS=$(readlink -f $(dirname {output.masked_genome}))
        mkdir -p "$OUTDIR_ABS"

        # [paths] tmp_dir when set, else /tmp
        TMP_BASE="{params.tmp_root}"
        TMP_BASE=${{TMP_BASE:-/tmp}}
        WORKDIR=$(mktemp -d "$TMP_BASE/red_{params.sample}_XXXXXX")
        trap "rm -rf $WORKDIR" EXIT

        echo "[INFO] Running Red repeat detection for {params.sample}" > "$LOG_ABS"

        # Red expects a directory of FASTA files (one per sequence or one file)
        GNM_DIR="$WORKDIR/gnm"
        MSK_DIR="$WORKDIR/msk"
        RPT_DIR="$WORKDIR/rpt"
        mkdir -p "$GNM_DIR" "$MSK_DIR" "$RPT_DIR"

        cp "$GENOME_ABS" "$GNM_DIR/genome.fa"

        echo "[INFO] Running Red..." >> "$LOG_ABS"
        Red -gnm "$GNM_DIR" -msk "$MSK_DIR" -rpt "$RPT_DIR" >> "$LOG_ABS" 2>&1

        # Red produces .msk files with the same basename as input
        if [ ! -f "$MSK_DIR/genome.msk" ]; then
            echo "[ERROR] Red failed to produce masked output" >> "$LOG_ABS"
            exit 1
        fi

        # Red output is already soft-masked (lowercase repeats)
        cp "$MSK_DIR/genome.msk" "$OUTDIR_ABS/genome.fa.masked"

        # Count masked bases for the log
        TOTAL=$(grep -v '^>' "$OUTDIR_ABS/genome.fa.masked" | tr -d '\n' | wc -c)
        MASKED=$(grep -v '^>' "$OUTDIR_ABS/genome.fa.masked" | tr -d '\n' | tr -cd 'a-z' | wc -c)
        PCT=$(awk "BEGIN {{printf \"%.1f\", 100.0*$MASKED/$TOTAL}}")
        echo "[INFO] Masked $MASKED / $TOTAL bp ($PCT%)" >> "$LOG_ABS"

        # Preserve Red repeat coordinates
        cp "$RPT_DIR"/*.rpt "$OUTDIR_ABS/" 2>/dev/null || true

        touch "$OUTDIR_ABS/.masking_complete"
        echo "[INFO] Red masking complete for {params.sample}" >> "$LOG_ABS"

        # Record software version
        VERSIONS_FILE="$(dirname "$OUTDIR_ABS")/software_versions.tsv"
        RED_VER=$(Red 2>&1 | head -1 || echo "Red 2018.09.10")
        ( flock 9; printf "Red\t%s\n" "$RED_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR_ABS=$(dirname "$OUTDIR_ABS")
        mkdir -p "$REPORT_DIR_ABS"
        source {script_dir}/report_citations.sh
        cite red "$REPORT_DIR_ABS" || true
        """
