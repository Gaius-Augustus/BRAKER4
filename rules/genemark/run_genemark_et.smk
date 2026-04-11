"""
Run GeneMark-ET for self-training gene prediction with RNA-Seq hints.

GeneMark-ET uses RNA-Seq evidence (intron hints) to train itself
on the genome, producing initial gene predictions that will be used
to train AUGUSTUS.

Input:
    - Masked genome FASTA
    - Hints file (from RNA-Seq)

Output:
    - genemark.gtf - Gene predictions
    - GeneMark model files

Container: teambraker/braker3:latest (contains gmes_petap.pl)
"""

rule run_genemark_et:
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
        hints="output/{sample}/hintsfile.gff"
    output:
        gtf="output/{sample}/genemark/genemark.gtf",
        log_file="output/{sample}/genemark/gmes.log"
    log:
        "logs/{sample}/genemark_et/genemark_et.log"
    benchmark:
        "benchmarks/{sample}/genemark_et/genemark_et.txt"
    threads: workflow.cores
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/genemark",
        min_contig=config.get("min_contig", 10000),
        fungus="--fungus" if config.get("fungus", False) else "",
        gm_max_intergenic=f"--max_intergenic {config.get('gm_max_intergenic')}" if config.get("gm_max_intergenic") else "",
        gcode=f"--gcode {config.get('translation_table')}" if config.get("translation_table", 1) != 1 else ""
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Create output directory
        mkdir -p {params.outdir}
        WORKDIR=$(pwd)

        # Resolve absolute paths before changing directory
        GENOME_ABS=$(readlink -f {input.genome})
        HINTS_ABS=$(readlink -f {input.hints})
        LOG_FILE_ABS=$(readlink -f {output.log_file})
        LOG_ABS=$(readlink -f {log})

        # Run GeneMark-ET
        cd {params.outdir}

        gmes_petap.pl \
            --verbose \
            --sequence $GENOME_ABS \
            --ET $HINTS_ABS \
            --cores {threads} \
            --min_contig {params.min_contig} \
            {params.fungus} \
            {params.gm_max_intergenic} \
            {params.gcode} \
            --soft_mask auto \
            1> $LOG_FILE_ABS \
            2> $LOG_ABS

        # Check if output was created
        if [ ! -f genemark.gtf ]; then
            echo "ERROR: GeneMark-ET failed to produce genemark.gtf" >> $LOG_ABS
            exit 1
        fi

        # Count predictions
        n_genes=$(grep -c $'\\tgene\\t' genemark.gtf || echo "0")
        echo "GeneMark-ET predicted $n_genes genes" >> $LOG_ABS

        # Record software versions
        cd $WORKDIR
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        GM_VER=$(grep -m1 '# GeneMark-ES Suite version' $(which gmes_petap.pl) 2>/dev/null | grep -oP 'version \K\S+' || true)
        GM_COMMIT=$(grep 'refs/remotes/origin/main' /opt/ETP/.git/packed-refs 2>/dev/null | awk '{{print substr($1,1,7)}}' || true)
        if [ -n "$GM_VER" ]; then
            ( flock 9; printf "GeneMark-ET\t%s (commit %s)\n" "$GM_VER" "$GM_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"
        else
            ( flock 9; printf "GeneMark-ET\tcommit %s\n" "$GM_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"
        fi

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite genemark_et "$REPORT_DIR"
        cite braker1 "$REPORT_DIR"
        """
