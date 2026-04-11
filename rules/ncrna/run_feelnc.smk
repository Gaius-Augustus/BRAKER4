"""
FEELnc: long non-coding RNA identification from transcriptome assembly.

FEELnc classifies transcripts from StringTie as protein-coding or lncRNA
based on coding potential, then categorizes lncRNAs by their genomic
relationship to protein-coding genes (intergenic, intronic, antisense).

Only runs when transcript evidence is available (ET, ETP, IsoSeq, dual modes).
ES and EP modes have no StringTie assembly and skip this step.

Container: quay.io/biocontainers/feelnc:0.2.1--pl5321hdfd78af_0
"""

FEELNC_CONTAINER = config.get(
    "feelnc_image",
    "docker://quay.io/biocontainers/feelnc:0.2--pl526_0"
)


def _get_feelnc_stringtie_gtf(wildcards):
    """Get StringTie GTF for FEELnc, routing by mode."""
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'dual':
        return f"output/{sample}/dual_etp_merged/transcripts_merged.gff"
    elif mode in ('etp', 'isoseq'):
        return f"output/{sample}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff"
    elif mode == 'et':
        return f"output/{sample}/stringtie/stringtie.gtf"
    else:
        raise ValueError(f"FEELnc requires transcript evidence but sample {sample} is in {mode} mode")


rule run_feelnc:
    """Run FEELnc to identify long non-coding RNAs from StringTie assembly."""
    input:
        stringtie=_get_feelnc_stringtie_gtf,
        braker_gtf="output/{sample}/braker.gtf",
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        lncrna_gff="output/{sample}/ncrna/lncRNAs.gff3",
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
        sample="{sample}",
        workdir=lambda wildcards: f"output/{wildcards.sample}/ncrna/feelnc_work"
    container:
        FEELNC_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir}

        GENOME_ABS=$(readlink -f {input.genome})
        STRINGTIE_ABS=$(readlink -f {input.stringtie})
        BRAKER_ABS=$(readlink -f {input.braker_gtf})

        echo "[INFO] Running FEELnc lncRNA identification..." > {log}
        echo "[INFO] StringTie input: $STRINGTIE_ABS" >> {log}
        echo "[INFO] Reference mRNA: $BRAKER_ABS" >> {log}

        cd {params.workdir}

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
        echo "[INFO] Fixed braker.gtf attributes for FEELnc compatibility" >> $OLDPWD/{log}

        # Step 1: Filter — remove short transcripts and those overlapping
        # protein-coding exons
        echo "[INFO] Step 1: FEELnc_filter.pl..." >> $OLDPWD/{log}
        FEELnc_filter.pl \
            -i $STRINGTIE_ABS \
            -a $BRAKER_ABS \
            --monoex=-1 \
            --size=200 \
            -p {threads} \
            > candidate_lncrna.gtf \
            2>> $OLDPWD/{log} || true

        n_candidates=$(awk -F'\t' '$3=="transcript"' candidate_lncrna.gtf 2>/dev/null | wc -l || echo 0)
        echo "[INFO] Filter produced $n_candidates candidate transcripts" >> $OLDPWD/{log}

        FEELNC_STATUS="no_candidates"
        n_final=0

        if [ "$n_candidates" -gt 0 ]; then
            # Step 2: Coding potential — classify as coding or non-coding
            echo "[INFO] Step 2: FEELnc_codpot.pl..." >> $OLDPWD/{log}
            FEELnc_codpot.pl \
                -i candidate_lncrna.gtf \
                -a $BRAKER_ABS \
                -g $GENOME_ABS \
                --mode=shuffle \
                --outdir=codpot_out \
                -p {threads} \
                2>> $OLDPWD/{log} || true

            LNCRNA_GTF=$(find codpot_out -name "*_lncRNA.gtf" 2>/dev/null | head -1)

            if [ -n "$LNCRNA_GTF" ] && [ -s "$LNCRNA_GTF" ]; then
                n_lncrna=$(awk -F'\t' '$3=="transcript"' "$LNCRNA_GTF" 2>/dev/null | wc -l || echo 0)
                echo "[INFO] Coding potential filter: $n_lncrna lncRNA transcripts" >> $OLDPWD/{log}

                # Step 3: Classify — categorize lncRNAs by relationship to mRNAs
                echo "[INFO] Step 3: FEELnc_classifier.pl..." >> $OLDPWD/{log}
                FEELnc_classifier.pl \
                    -i "$LNCRNA_GTF" \
                    -a $BRAKER_ABS \
                    > $OLDPWD/{output.classifier} \
                    2>> $OLDPWD/{log} || true

                # Convert lncRNA GTF to GFF3 with proper IDs
                awk -F'\t' -v OFS='\t' -v p="{params.sample}" '
                    BEGIN {{
                        print "##gff-version 3"
                        n=1
                    }}
                    /^#/ {{next}}
                    $3 == "transcript" || $3 == "exon" {{
                        if ($3 == "transcript") {{
                            id = p "-lncRNA_" n
                            $9 = "ID=" id ";Name=" id ";biotype=lncRNA"
                            n++
                        }}
                        print
                    }}
                ' "$LNCRNA_GTF" > $OLDPWD/{output.lncrna_gff}

                n_final=$(grep -c 'biotype=lncRNA' $OLDPWD/{output.lncrna_gff} || echo 0)
                echo "[INFO] FEELnc identified $n_final lncRNA transcripts" >> $OLDPWD/{log}
                FEELNC_STATUS="success"
            else
                echo "[INFO] No lncRNAs identified by coding potential filter" >> $OLDPWD/{log}
                FEELNC_STATUS="no_lncrna"
            fi
        fi

        cd $OLDPWD

        # Create empty output files if FEELnc didn't produce results
        if [ ! -f {output.lncrna_gff} ]; then
            echo "##gff-version 3" > {output.lncrna_gff}
        fi
        if [ ! -f {output.classifier} ]; then
            echo "# No lncRNA candidates" > {output.classifier}
        fi

        # Clean up working directory
        rm -rf {params.workdir}

        # Record software version (no flock — not available in FEELnc container)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        printf "FEELnc\t0.2\n" >> "$VERSIONS_FILE"

        # Citations (no flock — not available in FEELnc container)
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite feelnc "$REPORT_DIR" || true
        """
