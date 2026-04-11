"""
Run ProtHint to generate protein-based hints for gene prediction.

Container: teambraker/braker3:latest (contains prothint.py, DIAMOND, Spaln)
"""

rule run_prothint:
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
        outdir=lambda wildcards: f"output/{wildcards.sample}/prothint"
    threads: workflow.cores
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
        GENOME_ABS=$(readlink -f {input.genome})
        PROTEINS_ABS=$(readlink -f {input.proteins})
        GENEMARK_GTF_ABS=$(readlink -f {input.genemark_es})
        OUTDIR_ABS=$(readlink -f {params.outdir})

        cd $OUTDIR_ABS

        # Run prothint. Capture exit code without triggering set -e.
        # 'if cmd' is the ONLY set -e-safe pattern. No subshells.
        if prothint.py --threads={threads} --geneMarkGtf $GENEMARK_GTF_ABS $GENOME_ABS $PROTEINS_ABS > $OUTDIR_ABS/prothint_run.log 2>&1
        then
            PROTHINT_EXIT=0
        else
            PROTHINT_EXIT=$?
        fi

        cd $WORKDIR

        cp $OUTDIR_ABS/prothint_run.log {log} 2>/dev/null || true

        if [ ! -f $OUTDIR_ABS/prothint_augustus.gff ]
        then
            echo "ERROR: ProtHint failed (exit=$PROTHINT_EXIT), no prothint_augustus.gff" >> {log}
            exit 1
        fi

        cp $OUTDIR_ABS/prothint_augustus.gff {output.hints}

        if [ -f $OUTDIR_ABS/prothint.gff ]
        then
            cp $OUTDIR_ABS/prothint.gff {output.evidence}
        else
            touch {output.evidence}
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
        """
