"""
tRNAscan-SE: tRNA gene prediction.

Scans the genome for transfer RNA genes using tRNAscan-SE in eukaryotic mode.
Outputs GFF3 format natively via the --gff flag.

With trnascan_high_confidence_filter = 1, tRNAscan-SE additionally writes the
score breakdown (-H --detail) and secondary structures (-f), and
EukHighConfidenceFilter (shipped with tRNAscan-SE) removes pseudogenes and
low-scoring hits such as tRNA-derived SINEs (issue #49). Only tRNAs that pass
the filter are kept in tRNAs.gff3; tRNAs.txt stays the unfiltered result and
the filter's own output goes to tRNAs.highconf.{out,ss,log}.

Container: quay.io/biocontainers/trnascan-se:2.0.12--pl5321h031d066_0
"""


rule run_trnascan:
    """Run tRNAscan-SE to predict tRNA genes on the genome."""
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        gff="output/{sample}/ncrna/tRNAs.gff3",
        txt="output/{sample}/ncrna/tRNAs.txt"
    params:
        highconf=1 if config.get('trnascan_high_confidence_filter', False) else 0,
        outdir=lambda wildcards: f"output/{wildcards.sample}/ncrna"
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

        # EukHighConfidenceFilter needs the HMM/secondary-structure score
        # columns (-H --detail) and the secondary structure file (-f).
        HC_ARGS=""
        rm -f {params.outdir}/tRNAs.ss {params.outdir}/tRNAs.highconf.out \
              {params.outdir}/tRNAs.highconf.ss {params.outdir}/tRNAs.highconf.log
        if [ "{params.highconf}" = "1" ]; then
            HC_ARGS="-H --detail -f {params.outdir}/tRNAs.ss"
        fi

        LC_ALL=C tRNAscan-SE \
            -E \
            $HC_ARGS \
            --thread {threads} \
            -q \
            --forceow \
            -o {output.txt} \
            --gff {output.gff}.tmp \
            {input.genome} \
            2>> {log} || true

        if [ "{params.highconf}" = "1" ] && [ -s {output.gff}.tmp ] && grep -qv '^#' {output.gff}.tmp; then
            echo "[INFO] Running EukHighConfidenceFilter..." >> {log}
            LC_ALL=C EukHighConfidenceFilter \
                --result {output.txt} \
                --ss {params.outdir}/tRNAs.ss \
                --output {params.outdir} \
                --prefix tRNAs.highconf \
                --remove \
                >> {log} 2>&1
            # tRNAscan-SE names its GFF features <sequence>.trna<N> after the
            # "Sequence Name" and "tRNA #" columns of the result table; keep
            # only the tRNAs (and their exons) that survived the filter.
            awk -F'\t' -v OFS='\t' '
                FNR == NR {{
                    if (FNR > 3) {{ seq = $1; gsub(/[[:space:]]/, "", seq); n = $2; gsub(/[[:space:]]/, "", n); keep[seq ".trna" n] = 1 }}
                    next
                }}
                /^#/ {{print; next}}
                {{
                    key = ""
                    if ($3 == "exon") {{ if (match($9, /Parent=[^;]+/)) key = substr($9, RSTART + 7, RLENGTH - 7) }}
                    else if (match($9, /ID=[^;]+/)) key = substr($9, RSTART + 3, RLENGTH - 3)
                    if (key in keep) print
                }}
            ' {params.outdir}/tRNAs.highconf.out {output.gff}.tmp > {output.gff}.hc
            mv {output.gff}.hc {output.gff}.tmp
            n_raw=$(awk -F'\t' 'FNR > 3 && NF > 1' {output.txt} | wc -l)
            n_hc=$(awk -F'\t' 'FNR > 3 && NF > 1' {params.outdir}/tRNAs.highconf.out | wc -l)
            echo "[INFO] EukHighConfidenceFilter kept $n_hc of $n_raw tRNAscan-SE predictions" >> {log}
        fi

        # Ensure output exists even if no tRNAs found
        if [ -s {output.gff}.tmp ] && grep -qv '^#' {output.gff}.tmp; then
            # Prefix every ID and Parent with the sample name for uniqueness, so
            # exons keep pointing to their tRNA (IDs and Parents change together).
            awk -F'\t' -v OFS='\t' -v p="{wildcards.sample}" '
                /^#/ {{print; next}}
                {{
                    $9 = ";" $9
                    gsub(/;[[:space:]]*ID=/, ";ID=" p "-", $9)
                    gsub(/;[[:space:]]*Parent=/, ";Parent=" p "-", $9)
                    $9 = substr($9, 2)
                    print
                }}
            ' {output.gff}.tmp > {output.gff}
        else
            echo "##gff-version 3" > {output.gff}
            echo "[INFO] No tRNA genes found" >> {log}
        fi
        rm -f {output.gff}.tmp

        n_trna=$(awk -F'\t' '$3 == "tRNA"' {output.gff} | wc -l)
        echo "[INFO] tRNAs.gff3 contains $n_trna tRNA genes" >> {log}

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
