"""
Run ProtHint iteration 2 with AUGUSTUS predictions as additional seeds.

After the first AUGUSTUS prediction in EP mode, ProtHint is re-run using
the AUGUSTUS gene models as new seeds (--geneSeeds), along with the
previous GeneMark-ES seeds (--prevGeneSeeds) and Spaln alignments
(--prevSpalnGff). This produces improved protein hints.

The updated hints replace the iteration 1 ProtHint hints in the
hintsfile, then AUGUSTUS is re-run for the final prediction.

Mirrors braker.pl's run_prothint_iter2().

Input:
    - Genome FASTA, Protein sequences
    - AUGUSTUS iteration 1 predictions (augustus.hints.gtf)
    - GeneMark-ES predictions (previous seeds)
    - ProtHint iteration 1 directory (for Spaln GFF)
    - Iteration 1 hintsfile

Output:
    - prothint_hints_iter2.gff: Improved protein hints
    - hintsfile_iter2.gff: Updated hintsfile (old ProtHint hints replaced)

Container: teambraker/braker3:latest (contains prothint.py)
"""

rule run_prothint_iter2:
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
        proteins=lambda wildcards: get_protein_fasta(wildcards.sample),
        augustus_gtf="output/{sample}/augustus.hints.gtf",
        genemark_es="output/{sample}/GeneMark-ES/genemark.gtf",
        prothint_evidence="output/{sample}/prothint/prothint.gff",
        hintsfile="output/{sample}/hintsfile.gff"
    output:
        hints="output/{sample}/prothint_hints_iter2.gff",
        hintsfile="output/{sample}/hintsfile_iter2.gff"
    log:
        "logs/{sample}/prothint_iter2/prothint_iter2.log"
    benchmark:
        "benchmarks/{sample}/prothint_iter2/prothint_iter2.txt"
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/prothint_iter2"
    threads: workflow.cores
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # ProtHint legitimately returns non-zero on success, so we wrap ONLY the
        # prothint.py invocation with set +e / +pipefail. Everything else in this
        # rule runs under strict mode so a silent failure of cp, grep, or cat
        # cannot leave a "successful" rule with missing outputs (which would
        # confuse the downstream DAG).
        set -euo pipefail
        WORKDIR=$(pwd)
        mkdir -p {params.outdir}

        GENOME_ABS=$(readlink -f {input.genome})
        PROTEINS_ABS=$(readlink -f {input.proteins})
        AUGUSTUS_GTF_ABS=$(readlink -f {input.augustus_gtf})
        GENEMARK_ES_ABS=$(readlink -f {input.genemark_es})

        # Check for Spaln GFF from iteration 1. When present, passing it via
        # --prevSpalnGff lets ProtHint reuse iter1's Spaln alignments for seeds
        # that haven't changed and only re-align the new/modified ones. Without
        # it, iter2 runs Spaln from scratch, which is much slower on large
        # genomes (potentially hours longer). Result is the same either way,
        # but the fast path is a significant speedup when available.
        #
        # braker.pl always passes --prevSpalnGff unconditionally because it
        # renames Spaln/spaln.gff to Spaln/spaln_iter1.gff in move_aug_preds
        # before iter2 runs (so the file is guaranteed to exist by that point).
        # BRAKER4 leaves iter1's output unrenamed at output/{sample}/prothint/
        # Spaln/spaln.gff and checks for it here. In normal operation the file
        # IS there because prothint.py's --cleanup defaults to False, so this
        # conditional should always take the fast path. If it doesn't, we log
        # a loud warning so the user can investigate.
        mkdir -p $(dirname {log})
        SPALN_ARG=""
        PROTHINT_DIR=$(dirname $(readlink -f {input.prothint_evidence}))
        if [ -f "$PROTHINT_DIR/Spaln/spaln.gff" ]; then
            SPALN_ARG="--prevSpalnGff $PROTHINT_DIR/Spaln/spaln.gff"
            echo "[INFO] Using previous Spaln GFF from $PROTHINT_DIR/Spaln/spaln.gff" > $WORKDIR/{log}
        else
            echo "[WARNING] iteration-1 Spaln output not found at $PROTHINT_DIR/Spaln/spaln.gff" > $WORKDIR/{log}
            echo "[WARNING] Falling back to running Spaln from scratch. This WILL be slower." >> $WORKDIR/{log}
            echo "[WARNING] If you see this in a production run, check whether iter1 ProtHint" >> $WORKDIR/{log}
            echo "[WARNING] was invoked with --cleanup, or whether something deleted the file." >> $WORKDIR/{log}
        fi

        # Run ProtHint iteration 2 (briefly disable strict mode for prothint.py only)
        cd {params.outdir}
        set +e
        set +o pipefail
        prothint.py \
            --threads={threads} \
            --geneSeeds $AUGUSTUS_GTF_ABS \
            --prevGeneSeeds $GENEMARK_ES_ABS \
            $SPALN_ARG \
            $GENOME_ABS \
            $PROTEINS_ABS \
            > prothint_iter2_run.log 2>&1
        PROTHINT_EXIT=$?
        set -euo pipefail

        cd $WORKDIR
        mkdir -p $(dirname {log})
        cp {params.outdir}/prothint_iter2_run.log {log} 2>/dev/null || true

        if [ ! -f {params.outdir}/prothint_augustus.gff ]; then
            echo "ERROR: ProtHint iter2 failed (exit=$PROTHINT_EXIT)" >> {log}
            exit 1
        fi

        cp {params.outdir}/prothint_augustus.gff {output.hints}
        echo "ProtHint iter2 generated $(wc -l < {output.hints}) hints (exit=$PROTHINT_EXIT)" >> {log}

        # Build updated hintsfile: remove old ProtHint hints, add new ones.
        # grep -v exits 1 if no non-matching lines exist, which is benign here
        # (the redirect creates the file regardless), so we explicitly tolerate
        # exit codes 0 and 1 only.
        grep -v $'\\tProtHint\\t' {input.hintsfile} > {output.hintsfile} || [ $? -eq 1 ]
        cat {output.hints} >> {output.hintsfile}

        # Final sanity check: both outputs MUST exist before the rule completes,
        # otherwise downstream rules see "successful upstream" but missing inputs.
        if [ ! -s {output.hints} ]; then
            echo "ERROR: {output.hints} is missing or empty after rule completion" >&2
            exit 1
        fi
        if [ ! -s {output.hintsfile} ]; then
            echo "ERROR: {output.hintsfile} is missing or empty after rule completion" >&2
            exit 1
        fi

        n_hints_iter2=$(wc -l < {output.hints})
        n_hints_total=$(wc -l < {output.hintsfile})
        echo "Updated hintsfile: $n_hints_total total hints" >> {log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh || true
        cite prothint "$REPORT_DIR" || true
        """
