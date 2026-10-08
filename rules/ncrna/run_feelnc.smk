"""
FEELnc: long non-coding RNA identification from transcriptome assembly.

FEELnc classifies transcripts from StringTie as protein-coding or lncRNA
based on coding potential, then categorizes lncRNAs by their genomic
relationship to protein-coding genes (intergenic, intronic, antisense).

Only runs when transcript evidence is available (ET, ETP, IsoSeq, dual modes).
ES and EP modes have no StringTie assembly and skip this step.

Split into two rules:
  1. run_feelnc: FEELnc_filter.pl, FEELnc_codpot.pl and FEELnc_classifier.pl
     in the FEELnc container -> lncRNAs.gtf (FEELnc's exon lines; FEELnc
     writes no transcript lines) and feelnc_classifier.txt
  2. convert_feelnc_to_gff3: lncRNAs.gtf -> lncRNAs.gff3 (lnc_RNA + exons,
     scripts/feelnc_to_gff3.py; no container, uses host Python)

FEELnc cannot train on fewer than 100 candidates or fewer than 100 BRAKER
transcripts: the outputs are then empty but for a comment line that says so
(also in the log), as when no candidate is without coding potential. Every
other FEELnc error fails the job.

Container: quay.io/biocontainers/feelnc:0.2--pl526_0
"""


def _get_feelnc_stringtie_gtf(wildcards):
    """Get StringTie GTF for FEELnc, routing by mode."""
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'dual':
        return f"output/{sample}/GeneMark-ETP-isoseq/training.gtf"
        # return f"output/{sample}/dual_etp_merged/transcripts_merged.gff"
    if mode in ('etp', 'isoseq'):
        return f"output/{sample}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff"
    elif mode == 'et':
        return f"output/{sample}/stringtie/stringtie.gtf"
    else:
        raise ValueError(f"FEELnc requires transcript evidence but sample {sample} is in {mode} mode")


rule run_feelnc:
    """Run FEELnc to identify long non-coding RNAs from StringTie assembly.

    Scratch: the FEELnc work dir (filtered candidates, codpot_out/) is a
    private directory on the node-local disk (scripts/tmp_dir.sh, [paths]
    tmp_dir); the outputs are written to their run-dir paths directly.
    NEED 10 GB (fixed: the busybox FEELnc container has no stat for
    need_gb); with less free the job works in
    output/<sample>/ncrna/feelnc_work as before. In the FEELnc image
    /var/tmp is a symlink to /tmp, so Singularity's /var/tmp bind shadows
    /tmp: with tmp_dir empty and no $TMPDIR the scratch directory (and the
    free-space check) is on the host's /var/tmp.
    """
    input:
        stringtie=_get_feelnc_stringtie_gtf,
        braker_gtf="output/{sample}/braker.gtf",
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        lncrna_gtf="output/{sample}/ncrna/lncRNAs.gtf",
        classifier="output/{sample}/ncrna/feelnc_classifier.txt"
    log:
        "logs/{sample}/ncrna/feelnc.log"
    benchmark:
        "benchmarks/{sample}/ncrna/feelnc.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        workdir=lambda wildcards: f"output/{wildcards.sample}/ncrna/feelnc_work",
        tmp_root=TMP_ROOT
    container:
        FEELNC_CONTAINER
    shell:
        r"""
        set -euo pipefail
        source {script_dir}/tmp_dir.sh
        WORKDIR=$PWD
        LOG_ABS=$WORKDIR/{log}
        CLASSIFIER_ABS=$WORKDIR/{output.classifier}
        LNCRNA_GTF_ABS=$WORKDIR/{output.lncrna_gtf}

        GENOME_ABS=$(readlink -f {input.genome})
        STRINGTIE_ABS=$(readlink -f {input.stringtie})
        BRAKER_ABS=$(readlink -f {input.braker_gtf})

        echo "[INFO] Running FEELnc lncRNA identification..." > {log}
        echo "[INFO] StringTie input: $STRINGTIE_ABS" >> {log}
        echo "[INFO] Reference mRNA: $BRAKER_ABS" >> {log}

        # FEELnc's work dir goes to the node-local disk
        scratch_dir outDir "feelnc_{wildcards.sample}" "{params.tmp_root}" 10 \
            "$WORKDIR/{params.workdir}" 2>> {log}
        trap 'rm -rf -- "$SCRATCH"' EXIT
        cd "$outDir"

        # Fix braker.gtf: AUGUSTUS transcript/gene lines have bare IDs without
        # transcript_id/gene_id attributes. FEELnc requires standard GTF format.
        awk -F'\t' -v OFS='\t' '{{
            if ($3 == "gene" && $9 !~ /gene_id/) {{
                $9 = "gene_id \"" $9 "\";"
            }} else if ($3 == "transcript" && $9 !~ /transcript_id/) {{
                tid = $9; gid = tid; sub(/\.[^.]*$/, "", gid)
                $9 = "transcript_id \"" tid "\"; gene_id \"" gid "\";"
            }}
            print
        }}' $BRAKER_ABS > braker_fixed.gtf
        BRAKER_ABS=$(readlink -f braker_fixed.gtf)
        echo "[INFO] Fixed braker.gtf attributes for FEELnc compatibility" >> "$LOG_ABS"

        # The BioContainers image sets no FEELNCPATH; FEELnc_codpot.pl dies
        # without it. Every FEELnc error fails the job (set -e); nothing is
        # swallowed with || true.
        export FEELNCPATH=${{FEELNCPATH:-/usr/local}}
        export LC_ALL=C
        # transcripts of a GTF, by the transcript_id of its exon lines
        # (FEELnc writes no transcript lines)
        count_tx() {{
            awk -F'\t' '$3 == "exon" && match($9, /transcript_id "[^"]+"/) {{
                id = substr($9, RSTART, RLENGTH); if (!(id in seen)) {{ seen[id] = 1; n++ }}
            }} END {{ print n + 0 }}' "$1"
        }}
        # why no lncRNA was called: log, and a comment line in both outputs
        note() {{
            echo "[INFO] FEELnc: $1" >> "$LOG_ABS"
            echo "# FEELnc: $1" >> "$LNCRNA_GTF_ABS"
            echo "# FEELnc: $1" >> "$CLASSIFIER_ABS"
        }}
        : > "$LNCRNA_GTF_ABS"
        : > "$CLASSIFIER_ABS"

        # Step 1: candidates, the assembled transcripts of at least 200 bp
        # with more than one exon that do not overlap a BRAKER gene
        echo "[INFO] Step 1: FEELnc_filter.pl..." >> "$LOG_ABS"
        FEELnc_filter.pl \
            -i $STRINGTIE_ABS \
            -a $BRAKER_ABS \
            --monoex=-1 \
            --size=200 \
            -p {threads} \
            > candidate_lncrna.gtf \
            2>> "$LOG_ABS"
        n_cand=$(count_tx candidate_lncrna.gtf)
        n_mrna=$(count_tx $BRAKER_ABS)
        echo "[INFO] Filter produced $n_cand candidate transcripts ($n_mrna BRAKER transcripts)" >> "$LOG_ABS"

        if [ "$n_cand" -lt 100 ] || [ "$n_mrna" -lt 100 ]; then
            note "$n_cand candidate transcripts, $n_mrna annotated transcripts; FEELnc_codpot.pl needs at least 100 of each to train, no lncRNA called"
        else
            # Step 2: coding potential; a random forest trained on the BRAKER
            # mRNAs and shuffled copies of them keeps the candidates without
            echo "[INFO] Step 2: FEELnc_codpot.pl..." >> "$LOG_ABS"
            FEELnc_codpot.pl \
                -i candidate_lncrna.gtf \
                -a $BRAKER_ABS \
                -g $GENOME_ABS \
                --mode=shuffle \
                --outdir=codpot_out \
                -p {threads} \
                >> "$LOG_ABS" 2>&1
            LNC=codpot_out/candidate_lncrna.gtf.lncRNA.gtf
            if [ ! -f "$LNC" ]; then
                echo "[ERROR] FEELnc_codpot.pl exited 0 but did not write $LNC" >> "$LOG_ABS"
                exit 1
            fi
            n_lnc=$(count_tx "$LNC")
            if [ "$n_lnc" -eq 0 ]; then
                note "none of the $n_cand candidate transcripts is without coding potential, no lncRNA called"
            else
                cp "$LNC" "$LNCRNA_GTF_ABS"
                # Step 3: the coding genes next to each lncRNA
                echo "[INFO] Step 3: FEELnc_classifier.pl..." >> "$LOG_ABS"
                FEELnc_classifier.pl \
                    -i "$LNCRNA_GTF_ABS" \
                    -a $BRAKER_ABS \
                    > "$CLASSIFIER_ABS" \
                    2>> "$LOG_ABS"
                echo "[INFO] FEELnc: $n_lnc of $n_cand candidate transcripts are lncRNAs" >> "$LOG_ABS"
            fi
        fi

        cd "$WORKDIR"

        # Clean up working directory (fallback; on scratch the EXIT trap does)
        if [ -z "$SCRATCH" ]; then
            rm -rf {params.workdir}
        fi

        # Record software version (no flock — not available in FEELnc container)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        printf "FEELnc\t0.2\n" >> "$VERSIONS_FILE"

        # Citations (no flock — not available in FEELnc container)
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite feelnc "$REPORT_DIR" || true
        """


rule convert_feelnc_to_gff3:
    """Convert the FEELnc lncRNA GTF to GFF3 (runs on host, no container)."""
    input:
        gtf="output/{sample}/ncrna/lncRNAs.gtf"
    output:
        gff="output/{sample}/ncrna/lncRNAs.gff3"
    log:
        "logs/{sample}/ncrna/feelnc_to_gff3.log"
    benchmark:
        "benchmarks/{sample}/ncrna/feelnc_to_gff3.txt"
    params:
        sample="{sample}"
    threads: 1
    resources:
        mem_mb=0 if config['slurm_args'].get('skip_mem') else 4000,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        python3 {script_dir}/feelnc_to_gff3.py \
            --stem {params.sample} \
            {input.gtf} \
            -o {output.gff} \
            2> {log}
        """
