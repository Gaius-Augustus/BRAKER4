"""
Run ProtHint to generate protein-based hints for gene prediction.

ProtHint's bundled Perl-threads Spaln dispatcher can hang forever after
the last batch is enqueued (issue #98). Before running prothint.py, this
rule builds a symlinked shadow copy of the container's ProtHint bin tree
via scripts/make_prothint_shadow.sh, with only the Spaln dispatcher
replaced by scripts/spaln_dispatcher.py, and runs prothint.py from there.
Falls back to the container's own prothint.py if the shadow cannot be
built.

Container: teambraker/braker3:latest (contains prothint.py, DIAMOND, Spaln)
"""

rule run_prothint:
    """ProtHint iteration 1.

    Scratch: prothint.py (DIAMOND db and hits, Spaln per-seed files) and
    the ProtHint shadow tree run on the node-local disk
    (scripts/tmp_dir.sh, [paths] tmp_dir). prothint.gff and
    Spaln/spaln.gff (read by run_prothint_iter2) are copied back to
    output/<sample>/prothint/, prothint_augustus.gff to the hints output,
    prothint_run.log into the rule log. NEED 10 x (genome + proteins)
    + 5 GB; with less free the job works in output/<sample>/prothint/ as
    before.
    """
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
        proteins=lambda wildcards: get_protein_fasta(wildcards.sample),
        genemark_es="output/{sample}/GeneMark-ES/genemark.gtf"
    output:
        hints="output/{sample}/prothint_hints.gff",
        evidence="output/{sample}/prothint/prothint.gff"
    log:
        "logs/{sample}/prothint/prothint.log"
    benchmark:
        "benchmarks/{sample}/prothint/prothint.txt"
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/prothint",
        tmp_root=TMP_ROOT
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Disable set -e for this rule: prothint.py returns non-zero on success
        # and various bash/Singularity/SLURM interactions make it impossible
        # to reliably capture. We check outputs explicitly instead.
        set +e
        set +o pipefail
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        WORKDIR=$(pwd)
        source {script_dir}/tmp_dir.sh
        GENOME_ABS=$(readlink -f {input.genome})
        PROTEINS_ABS=$(readlink -f {input.proteins})
        GENEMARK_GTF_ABS=$(readlink -f {input.genemark_es})
        finalDir=$(readlink -f {params.outdir})
        LOG_ABS=$WORKDIR/{log}
        HINTS_ABS=$WORKDIR/{output.hints}
        EVIDENCE_ABS=$WORKDIR/{output.evidence}

        # ProtHint's work files (DIAMOND, Spaln) go to the node-local disk.
        : > "$LOG_ABS"
        scratch_dir outDir "prothint_{wildcards.sample}" "{params.tmp_root}" \
            "$(need_gb 10 "$GENOME_ABS" "$PROTEINS_ABS")" "$finalDir" 2>> "$LOG_ABS"
        trap 'rm -rf -- "$SCRATCH"' EXIT

        # ProtHint's Perl-threads Spaln dispatcher can hang after the last
        # batch is enqueued (#98). Run ProtHint from a symlinked copy of the
        # container's bin tree that uses scripts/spaln_dispatcher.py instead.
        SHADOW_DIR=$outDir.bin
        SHADOW_MSG=$(bash {script_dir}/make_prothint_shadow.sh $SHADOW_DIR 2>&1)
        if [ -x $SHADOW_DIR/gmes/ProtHint/bin/prothint.py ]
        then
            PROTHINT=$SHADOW_DIR/gmes/ProtHint/bin/prothint.py
            SHADOW_MSG="[INFO] ProtHint run with the process-based Spaln dispatcher (#98)"
        else
            PROTHINT=prothint.py
            SHADOW_MSG="[WARNING] could not install the Spaln dispatcher fix (#98), used the container's ProtHint: $SHADOW_MSG"
        fi

        cd "$outDir"

        # Run prothint. Capture exit code without triggering set -e.
        # 'if cmd' is the ONLY set -e-safe pattern. No subshells.
        if $PROTHINT --threads={threads} --geneMarkGtf $GENEMARK_GTF_ABS $GENOME_ABS $PROTEINS_ABS > "$outDir/prothint_run.log" 2>&1
        then
            PROTHINT_EXIT=0
        else
            PROTHINT_EXIT=$?
        fi

        cd $WORKDIR

        cat "$outDir/prothint_run.log" >> {log} 2>/dev/null || true
        echo "$SHADOW_MSG" >> {log}

        if [ ! -f "$outDir/prothint_augustus.gff" ]
        then
            echo "ERROR: ProtHint failed (exit=$PROTHINT_EXIT), no prothint_augustus.gff" >> {log}
            exit 1
        fi

        cp "$outDir/prothint_augustus.gff" "$HINTS_ABS"

        if [ -f "$outDir/prothint.gff" ]
        then
            cp "$outDir/prothint.gff" "$EVIDENCE_ABS" 2>/dev/null
        else
            touch "$EVIDENCE_ABS"
        fi
        # run_prothint_iter2 reuses iteration 1's Spaln alignments.
        copy_back "$outDir" "$finalDir" Spaln/spaln.gff
        if [ ! -f "$HINTS_ABS" ] || [ ! -f "$EVIDENCE_ABS" ]
        then
            echo "ERROR: {output.hints} or {output.evidence} missing after the copy back" >> {log}
            exit 1
        fi

        n_hints=$(wc -l < {output.hints})
        echo "ProtHint generated $n_hints hints (exit=$PROTHINT_EXIT)" >> {log}

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        PH_VER=$(prothint.py --version 2>&1 | awk '{{print $NF}}' || true)
        DM_VER=$(diamond version 2>&1 | awk '{{print $NF}}' || true)
        SPALN_VER=$(/opt/ETP/bin/gmes/ProtHint/dependencies/spaln 2>&1 | grep -oP 'version \K\S+' | head -1 || true)
        ( flock 9
          printf "ProtHint\t%s\n" "$PH_VER" >> "$VERSIONS_FILE"
          printf "DIAMOND\t%s\n" "$DM_VER" >> "$VERSIONS_FILE"
          printf "Spaln\t%s\n" "$SPALN_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh || true
        cite prothint "$REPORT_DIR" || true
        cite diamond "$REPORT_DIR" || true
        cite spaln "$REPORT_DIR" || true

        # Remove ProTHint working files (DIAMOND databases, intermediate files).
        # Keep: prothint.gff (Snakemake output) and Spaln/spaln.gff (reused by
        # run_prothint_iter2 via --prevSpalnGff to avoid re-running Spaln).
        if [ -z "$SCRATCH" ]
        then
            find {params.outdir} -mindepth 1 -type f \
                ! -name 'prothint.gff' \
                ! -path '*/Spaln/spaln.gff' \
                -delete 2>/dev/null || true
            find {params.outdir} -mindepth 1 -type d -empty -delete 2>/dev/null || true
        fi
        rm -rf $SHADOW_DIR
        """
