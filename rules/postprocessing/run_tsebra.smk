"""
Merge gene predictions using TSEBRA.

The semantics of TSEBRA's two input flags are critical for parity with
braker.pl:

  --keep_gtf  Forced gene sets — every transcript survives unchanged.
  --gtf       Candidate gene sets — TSEBRA filters each transcript
              against the hints file using the cfg thresholds (e.g. the
              braker3.cfg requires intron_support=1.0).

Mode-by-mode behavior in this rule mirrors braker.pl:

  ES, EP, ET    Pass everything via --keep_gtf (matches braker.pl ES/EP/ET).
                AUGUSTUS already ran with hints; the genemark.f.good set was
                already filtered upstream by filterGenemark.pl. TSEBRA only
                deduplicates.

  ETP, IsoSeq   Pass augustus.hints.fixed.gtf and GeneMark-ETP/genemark.gtf
                via --gtf so TSEBRA filters AUGUSTUS predictions against the
                hints. Pass GeneMark-ETP/training.gtf (the HC training gene
                set) via --keep_gtf to force-retain it. Matches braker.pl's
                ETP TSEBRA call (scripts/braker.pl, sub merge_transcript_sets_with_tsebra,
                called from line ~9783).

  dual          As ETP, with the GeneMark-ETP-isoseq/ outputs of the single
                combined GeneMark-ETP run (short-read and IsoSeq BAMs together).

  Config        braker3.cfg (matches braker.pl). Differs from default.cfg
                only in e_4 (0.20 vs 0.18); the load-bearing setting in both
                is intron_support 1.0.
"""


def get_tsebra_inputs(wildcards):
    """
    Mode-dependent inputs for the run_tsebra rule.

    Snakemake's string-substitution pass scans EVERY {input.X} reference in
    the shell block at parse time, even ones inside case branches that will
    never execute for the current mode. So every input key referenced
    anywhere in the shell must always exist in the dict — we alias unused
    keys to augustus_gtf (which is always present and already a real input,
    so no extra dependency is added).
    """
    sample = wildcards.sample
    mode = get_effective_mode(sample)
    augustus_gtf = f"output/{sample}/augustus.hints.fixed.gtf"

    inputs = {
        "augustus_gtf":   augustus_gtf,
        "genemark_gtf":   augustus_gtf,
        "training_gtf":   augustus_gtf,
        "hintsfile":      augustus_gtf,
        "keep_genes":     augustus_gtf,
        "genome_fai":     f"output/{sample}/genome.fa.fai",
    }

    if mode in ("etp", "isoseq"):
        inputs.update({
            "genemark_gtf": f"output/{sample}/GeneMark-ETP/genemark.gtf",
            "training_gtf": f"output/{sample}/GeneMark-ETP/training.gtf",
            "hintsfile":    get_augustus_hintsfile(sample),
        })
    elif mode == "dual":
        inputs.update({
            "genemark_gtf": f"output/{sample}/GeneMark-ETP-isoseq/genemark.gtf",
            "training_gtf": f"output/{sample}/GeneMark-ETP-isoseq/training.gtf",
            "hintsfile":    get_augustus_hintsfile(sample),
        })
    else:
        # ES / EP / ET — single keep-everything invocation against the
        # already-filtered GeneMark good set.
        inputs.update({
            "keep_genes": f"output/{sample}/genemark/genemark.f.good.gtf",
            "hintsfile":  get_augustus_hintsfile(sample),
        })

    return with_mode_decisions(wildcards, inputs)


rule run_tsebra:
    input:
        unpack(get_tsebra_inputs)
    output:
        braker_merged_gtf="output/{sample}/braker.tsebra.raw.gtf",
        tsebra_log="output/{sample}/tsebra.log"
    benchmark:
        "benchmarks/{sample}/run_tsebra/run_tsebra.txt"
    params:
        mode=lambda w: effective_mode_on_disk(w.sample),
        filter_se=config['filter_single_exon_genes']
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        OUT_DIR=$(dirname {output.braker_merged_gtf})

        echo "[INFO] ===== RUNNING TSEBRA TO MERGE PREDICTIONS =====" | tee {output.tsebra_log}
        echo "[INFO] Mode: {params.mode}" | tee -a {output.tsebra_log}

        # braker3.cfg matches braker.pl. Default has intron_support=1.0,
        # which is the load-bearing threshold that filters AUGUSTUS
        # predictions whose intron chains aren't fully supported by hints.
        TSEBRA_CFG=/opt/TSEBRA/config/braker3.cfg
        if [ ! -f "$TSEBRA_CFG" ]; then
            echo "[ERROR] braker3.cfg not found at $TSEBRA_CFG" | tee -a {output.tsebra_log}
            exit 1
        fi
        echo "[INFO] TSEBRA config: $TSEBRA_CFG" | tee -a {output.tsebra_log}

        GENOME_SIZE=$(awk '{{sum+=$2}} END {{printf "%.0f\n", sum+0}}' {input.genome_fai})
        # [PARAMS] filter_single_exon_genes: auto (> 300 Mbp), 1 (always), 0 (never)
        FILTER_SE_MODE={params.filter_se}
        FILTER_SE_ARG=""
        if [ "$FILTER_SE_MODE" = "1" ] || {{ [ "$FILTER_SE_MODE" = "auto" ] && [ "$GENOME_SIZE" -gt 300000000 ]; }}; then
            FILTER_SE_ARG="--filter_single_exon_genes"
            echo "[INFO] Genome $GENOME_SIZE bp, filter_single_exon_genes=$FILTER_SE_MODE: single-exon gene filter active" | tee -a {output.tsebra_log}
        else
            echo "[INFO] Genome $GENOME_SIZE bp, filter_single_exon_genes=$FILTER_SE_MODE: single-exon gene filter inactive" | tee -a {output.tsebra_log}
        fi

        TSEBRA_TMP={output.braker_merged_gtf}.tmp

        case "{params.mode}" in
            etp|isoseq|dual)
                # Single ETP-style merge: AUGUSTUS + raw GeneMark predictions
                # are CANDIDATES (filtered against the hints), HC training
                # genes are FORCED via --keep_gtf. Mirrors braker.pl ETP.
                AUG_GENES=$(awk '$3=="gene"{{n++}}END{{print n+0}}' {input.augustus_gtf})
                GM_GENES=$(awk -F'\t' '$3=="CDS"' {input.genemark_gtf} | grep -oE 'gene_id "[^"]+"' | sort -u | wc -l)
                KEEP_GENES=$(grep -oE 'gene_id "[^"]+"' {input.training_gtf} | sort -u | wc -l)
                echo "[INFO] Inputs: AUGUSTUS=$AUG_GENES, GeneMark=$GM_GENES (candidates), training=$KEEP_GENES (forced)" | tee -a {output.tsebra_log}

                tsebra.py \
                    --gtf {input.augustus_gtf},{input.genemark_gtf} \
                    --keep_gtf {input.training_gtf} \
                    --hintfiles {input.hintsfile} \
                    $FILTER_SE_ARG \
                    --cfg $TSEBRA_CFG \
                    --out $TSEBRA_TMP 2>&1 | tee -a {output.tsebra_log}
                ;;

            *)
                # ES / EP / ET — keep-everything (matches braker.pl). For ES,
                # genemark.f.good.gtf is the unfiltered sorted set; for EP/ET
                # it's the filterGenemark.pl-filtered "good" set.
                AUG_GENES=$(awk '$3=="gene"{{n++}}END{{print n+0}}' {input.augustus_gtf})
                KEEP_GENES=$(grep -oE 'gene_id "[^"]+"' {input.keep_genes} | sort -u | wc -l)
                echo "[INFO] Inputs: AUGUSTUS=$AUG_GENES, GeneMark good=$KEEP_GENES (both forced)" | tee -a {output.tsebra_log}

                # ES has no hints — only pass --hintfiles if non-empty.
                HINTS_ARG=""
                if [ -s {input.hintsfile} ]; then
                    HINTS_ARG="--hintfiles {input.hintsfile}"
                fi

                tsebra.py \
                    --keep_gtf {input.augustus_gtf},{input.keep_genes} \
                    $HINTS_ARG \
                    $FILTER_SE_ARG \
                    --cfg $TSEBRA_CFG \
                    --out $TSEBRA_TMP 2>&1 | tee -a {output.tsebra_log}
                ;;
        esac

        echo "[INFO] Renaming gene and transcript IDs..." | tee -a {output.tsebra_log}
        /opt/TSEBRA/bin/rename_gtf.py \
            --gtf $TSEBRA_TMP \
            --out {output.braker_merged_gtf} \
            2>&1 | tee -a {output.tsebra_log}
        rm -f $TSEBRA_TMP

        MERGED_GENES=$(awk '$3=="gene"{{n++}}END{{print n+0}}' {output.braker_merged_gtf})
        MERGED_TX=$(awk '$3=="transcript"{{n++}}END{{print n+0}}' {output.braker_merged_gtf})
        echo "[INFO] Output: $MERGED_GENES genes, $MERGED_TX transcripts" | tee -a {output.tsebra_log}
        echo "[INFO] TSEBRA merge completed successfully" | tee -a {output.tsebra_log}

        # Record software version (TSEBRA has no --version flag)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        TSEBRA_COMMIT=$(grep 'refs/remotes/origin/main' /opt/TSEBRA/.git/packed-refs 2>/dev/null | awk '{{print substr($1,1,7)}}' || true)
        ( flock 9; printf "TSEBRA\tcommit %s\n" "$TSEBRA_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite tsebra "$REPORT_DIR"
        """
