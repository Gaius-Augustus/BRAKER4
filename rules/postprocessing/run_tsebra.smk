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

  dual          Three-stage merge:
                  Stage A — ETP-style merge for the short-read run alone.
                  Stage B — ETP-style merge for the IsoSeq run alone.
                  Stage C — final merge keeping both runs' training genes
                            via --keep_gtf, with both genemark.gtf files and
                            augustus.hints.fixed.gtf via --gtf, and both
                            per-run hintsfiles via --hintfiles (passed as
                            ONE comma-separated argument; TSEBRA's argparse
                            does not accept multiple space-separated values).
                Stages A and B are produced by rule run_tsebra_etp_per_run
                and surface as side outputs braker_etp_{sr,iso}.tsebra.gtf
                for inspection. Stage C is what becomes braker.tsebra.raw.gtf.

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
    so no extra dependency is added). Same pattern as best_by_compleasm.smk.
    """
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    augustus_gtf = f"output/{sample}/augustus.hints.fixed.gtf"

    inputs = {
        "augustus_gtf":   augustus_gtf,
        "braker_etp_sr":  augustus_gtf,
        "braker_etp_iso": augustus_gtf,
        "genemark_gtf":   augustus_gtf,
        "training_gtf":   augustus_gtf,
        "genemark_sr":    augustus_gtf,
        "genemark_iso":   augustus_gtf,
        "training_sr":    augustus_gtf,
        "training_iso":   augustus_gtf,
        "hintsfile":      augustus_gtf,
        "hintsfile_sr":   augustus_gtf,
        "hintsfile_iso":  augustus_gtf,
        "keep_genes":     augustus_gtf,
    }

    if mode == "dual":
        # Stage C — depends on per-run outputs (stages A, B) plus the raw
        # augustus / genemark / training / hintsfiles needed for the final
        # merge command.
        inputs.update({
            "braker_etp_sr":  f"output/{sample}/braker_etp_sr.tsebra.gtf",
            "braker_etp_iso": f"output/{sample}/braker_etp_iso.tsebra.gtf",
            "genemark_sr":  f"output/{sample}/GeneMark-ETP/genemark.gtf",
            "genemark_iso": f"output/{sample}/GeneMark-ETP-isoseq/genemark.gtf",
            "training_sr":  f"output/{sample}/GeneMark-ETP/training.gtf",
            "training_iso": f"output/{sample}/GeneMark-ETP-isoseq/training.gtf",
            "hintsfile_sr":  f"output/{sample}/etp_hints.gff",
            "hintsfile_iso": f"output/{sample}/etp_hints_isoseq.gff",
        })
    elif mode in ("etp", "isoseq"):
        inputs.update({
            "genemark_gtf": f"output/{sample}/GeneMark-ETP/genemark.gtf",
            "training_gtf": f"output/{sample}/GeneMark-ETP/training.gtf",
            "hintsfile":    f"output/{sample}/hintsfile.gff",
        })
    else:
        # ES / EP / ET — single keep-everything invocation against the
        # already-filtered GeneMark good set.
        inputs.update({
            "keep_genes": f"output/{sample}/genemark/genemark.f.good.gtf",
            "hintsfile":  f"output/{sample}/hintsfile.gff",
        })

    return inputs


rule run_tsebra_etp_per_run:
    """
    Per-run ETP-style TSEBRA merge for dual mode (stages A and B).

    Runs once with etp_run=sr (short-read GeneMark-ETP) and once with
    etp_run=iso (IsoSeq GeneMark-ETP). Each invocation reproduces what
    braker.pl's ETP mode would do for that single run in isolation:
    augustus.hints.fixed.gtf and the per-run genemark.gtf go in via --gtf
    (TSEBRA filters them against the per-run hintsfile), and the per-run
    HC training set goes in via --keep_gtf.

    The same augustus.hints.fixed.gtf is reused for both stages because
    BRAKER4 only runs AUGUSTUS once in dual mode (with the merged
    hintsfile). The per-run filter still discriminates between the two
    by using each run's own hintsfile to score AUGUSTUS transcripts.

    This rule only fires when run_tsebra (dual mode) declares the
    braker_etp_{sr,iso}.tsebra.gtf paths as inputs.
    """
    input:
        augustus_gtf="output/{sample}/augustus.hints.fixed.gtf",
        genemark_gtf=lambda w: (
            f"output/{w.sample}/GeneMark-ETP/genemark.gtf" if w.etp_run == "sr"
            else f"output/{w.sample}/GeneMark-ETP-isoseq/genemark.gtf"
        ),
        training_gtf=lambda w: (
            f"output/{w.sample}/GeneMark-ETP/training.gtf" if w.etp_run == "sr"
            else f"output/{w.sample}/GeneMark-ETP-isoseq/training.gtf"
        ),
        hintsfile=lambda w: (
            f"output/{w.sample}/etp_hints.gff" if w.etp_run == "sr"
            else f"output/{w.sample}/etp_hints_isoseq.gff"
        ),
    output:
        braker_per_run="output/{sample}/braker_etp_{etp_run}.tsebra.gtf",
        log_file="output/{sample}/braker_etp_{etp_run}.tsebra.log",
    wildcard_constraints:
        etp_run="sr|iso",
    benchmark:
        "benchmarks/{sample}/run_tsebra_etp_per_run/{etp_run}.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== TSEBRA per-run ETP merge ({wildcards.etp_run}) =====" | tee {output.log_file}
        echo "[INFO] AUGUSTUS GTF: {input.augustus_gtf}" | tee -a {output.log_file}
        echo "[INFO] GeneMark GTF: {input.genemark_gtf}" | tee -a {output.log_file}
        echo "[INFO] Training GTF (forced): {input.training_gtf}" | tee -a {output.log_file}
        echo "[INFO] Hintsfile: {input.hintsfile}" | tee -a {output.log_file}

        # braker3.cfg matches braker.pl's TSEBRA invocation in ETP mode.
        TSEBRA_CFG=/opt/TSEBRA/config/braker3.cfg
        if [ ! -f "$TSEBRA_CFG" ]; then
            echo "[ERROR] braker3.cfg not found at $TSEBRA_CFG" | tee -a {output.log_file}
            exit 1
        fi
        echo "[INFO] TSEBRA config: $TSEBRA_CFG" | tee -a {output.log_file}

        TSEBRA_TMP={output.braker_per_run}.tmp
        tsebra.py \
            --gtf {input.augustus_gtf},{input.genemark_gtf} \
            --keep_gtf {input.training_gtf} \
            --hintfiles {input.hintsfile} \
            --cfg $TSEBRA_CFG \
            --out $TSEBRA_TMP 2>&1 | tee -a {output.log_file}

        /opt/TSEBRA/bin/rename_gtf.py \
            --gtf $TSEBRA_TMP \
            --out {output.braker_per_run} \
            2>&1 | tee -a {output.log_file}
        rm -f $TSEBRA_TMP

        N_GENES=$(grep -cP '\tgene\t' {output.braker_per_run} || echo 0)
        N_TX=$(grep -cP '\ttranscript\t' {output.braker_per_run} || echo 0)
        echo "[INFO] Per-run output ({wildcards.etp_run}): $N_GENES genes, $N_TX transcripts" | tee -a {output.log_file}
        """


rule run_tsebra:
    input:
        unpack(get_tsebra_inputs)
    output:
        braker_merged_gtf="output/{sample}/braker.tsebra.raw.gtf",
        tsebra_log="output/{sample}/tsebra.log"
    benchmark:
        "benchmarks/{sample}/run_tsebra/run_tsebra.txt"
    params:
        mode=lambda w: get_braker_mode(w.sample)
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
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

        TSEBRA_TMP={output.braker_merged_gtf}.tmp

        case "{params.mode}" in
            etp|isoseq)
                # Single ETP-style merge: AUGUSTUS + raw GeneMark predictions
                # are CANDIDATES (filtered against the hints), HC training
                # genes are FORCED via --keep_gtf. Mirrors braker.pl ETP.
                AUG_GENES=$(grep -cP '\tgene\t' {input.augustus_gtf} || echo 0)
                GM_GENES=$(awk -F'\t' '$3=="CDS"' {input.genemark_gtf} | grep -oP 'gene_id "[^"]+"' | sort -u | wc -l)
                KEEP_GENES=$(grep -oP 'gene_id "[^"]+"' {input.training_gtf} | sort -u | wc -l)
                echo "[INFO] Inputs: AUGUSTUS=$AUG_GENES, GeneMark=$GM_GENES (candidates), training=$KEEP_GENES (forced)" | tee -a {output.tsebra_log}

                tsebra.py \
                    --gtf {input.augustus_gtf},{input.genemark_gtf} \
                    --keep_gtf {input.training_gtf} \
                    --hintfiles {input.hintsfile} \
                    --cfg $TSEBRA_CFG \
                    --out $TSEBRA_TMP 2>&1 | tee -a {output.tsebra_log}
                ;;

            dual)
                # Stage C — final merge for dual mode. Stages A and B are
                # already on disk as braker_etp_{{sr,iso}}.tsebra.gtf side
                # outputs from rule run_tsebra_etp_per_run. They're listed
                # as inputs purely to enforce execution order.
                echo "[INFO] Per-run stages already complete:" | tee -a {output.tsebra_log}
                echo "[INFO]   stage A (sr):  {input.braker_etp_sr}" | tee -a {output.tsebra_log}
                echo "[INFO]   stage B (iso): {input.braker_etp_iso}" | tee -a {output.tsebra_log}

                # GeneMark and training files use MSTRG.* IDs from
                # StringTie which collide between the two ETP runs. Prefix
                # them before passing to TSEBRA so the IDs are disjoint.
                GM_SR_PFX=$OUT_DIR/_dual_gm_sr.gtf
                GM_ISO_PFX=$OUT_DIR/_dual_gm_iso.gtf
                TR_SR_PFX=$OUT_DIR/_dual_train_sr.gtf
                TR_ISO_PFX=$OUT_DIR/_dual_train_iso.gtf

                sed 's/gene_id "\([^"]*\)"/gene_id "sr_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "sr_\1"/g'  {input.genemark_sr}  > $GM_SR_PFX
                sed 's/gene_id "\([^"]*\)"/gene_id "iso_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "iso_\1"/g' {input.genemark_iso} > $GM_ISO_PFX
                sed 's/gene_id "\([^"]*\)"/gene_id "sr_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "sr_\1"/g'  {input.training_sr}  > $TR_SR_PFX
                sed 's/gene_id "\([^"]*\)"/gene_id "iso_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "iso_\1"/g' {input.training_iso} > $TR_ISO_PFX

                AUG_GENES=$(grep -cP '\tgene\t' {input.augustus_gtf} || echo 0)
                GM_SR_GENES=$(awk -F'\t' '$3=="CDS"' $GM_SR_PFX | grep -oP 'gene_id "[^"]+"' | sort -u | wc -l)
                GM_ISO_GENES=$(awk -F'\t' '$3=="CDS"' $GM_ISO_PFX | grep -oP 'gene_id "[^"]+"' | sort -u | wc -l)
                TR_SR_GENES=$(grep -oP 'gene_id "[^"]+"' $TR_SR_PFX | sort -u | wc -l)
                TR_ISO_GENES=$(grep -oP 'gene_id "[^"]+"' $TR_ISO_PFX | sort -u | wc -l)
                echo "[INFO] Stage C inputs:" | tee -a {output.tsebra_log}
                echo "[INFO]   AUGUSTUS=$AUG_GENES (candidate)" | tee -a {output.tsebra_log}
                echo "[INFO]   GeneMark sr=$GM_SR_GENES iso=$GM_ISO_GENES (candidates)" | tee -a {output.tsebra_log}
                echo "[INFO]   training sr=$TR_SR_GENES iso=$TR_ISO_GENES (forced)" | tee -a {output.tsebra_log}

                # TSEBRA's --hintfiles takes ONE comma-separated argument,
                # not multiple space-separated values (matches braker.pl
                # sub merge_transcript_sets_with_tsebra at scripts/braker.pl
                # line 8456: join(',', @hintfiles)).
                tsebra.py \
                    --gtf {input.augustus_gtf},$GM_SR_PFX,$GM_ISO_PFX \
                    --keep_gtf $TR_SR_PFX,$TR_ISO_PFX \
                    --hintfiles {input.hintsfile_sr},{input.hintsfile_iso} \
                    --cfg $TSEBRA_CFG \
                    --out $TSEBRA_TMP 2>&1 | tee -a {output.tsebra_log}

                rm -f $GM_SR_PFX $GM_ISO_PFX $TR_SR_PFX $TR_ISO_PFX
                ;;

            *)
                # ES / EP / ET — keep-everything (matches braker.pl). For ES,
                # genemark.f.good.gtf is the unfiltered sorted set; for EP/ET
                # it's the filterGenemark.pl-filtered "good" set.
                AUG_GENES=$(grep -cP '\tgene\t' {input.augustus_gtf} || echo 0)
                KEEP_GENES=$(grep -oP 'gene_id "[^"]+"' {input.keep_genes} | sort -u | wc -l)
                echo "[INFO] Inputs: AUGUSTUS=$AUG_GENES, GeneMark good=$KEEP_GENES (both forced)" | tee -a {output.tsebra_log}

                # ES has no hints — only pass --hintfiles if non-empty.
                HINTS_ARG=""
                if [ -s {input.hintsfile} ]; then
                    HINTS_ARG="--hintfiles {input.hintsfile}"
                fi

                tsebra.py \
                    --keep_gtf {input.augustus_gtf},{input.keep_genes} \
                    $HINTS_ARG \
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

        MERGED_GENES=$(grep -cP '\tgene\t' {output.braker_merged_gtf} || echo 0)
        MERGED_TX=$(grep -cP '\ttranscript\t' {output.braker_merged_gtf} || echo 0)
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
