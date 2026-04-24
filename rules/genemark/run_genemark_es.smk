"""
Run GeneMark-ES for ab initio gene prediction.

In EP mode, GeneMark-ES provides initial gene predictions that seed
ProtHint's protein-to-genome alignment. This is necessary because
ProtHint needs an initial set of gene models to guide protein alignment.

Mirrors braker.pl's GeneMark_ES() function.

Input:
    - Masked genome FASTA

Output:
    - genemark.gtf - Ab initio gene predictions

Container: teambraker/braker3:latest (contains gmes_petap.pl)
"""

rule run_genemark_es:
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        gtf="output/{sample}/GeneMark-ES/genemark.gtf",
        log_file="output/{sample}/GeneMark-ES/gmes.log"
    log:
        "logs/{sample}/genemark_es/genemark_es.log"
    benchmark:
        "benchmarks/{sample}/genemark_es/genemark_es.txt"
    threads: workflow.cores
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/GeneMark-ES",
        min_contig=config.get("min_contig", 10000),
        fungus="--fungus" if config.get("fungus", False) else "",
        gm_max_intergenic=f"--max_intergenic {config.get('gm_max_intergenic')}" if config.get("gm_max_intergenic") else "",
        gcode=f"--gcode {config.get('translation_table')}" if config.get("translation_table", 1) != 1 else "",
        gc_donor=config.get("gc_donor", 0.001)
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        mkdir -p {params.outdir}
        WORKDIR=$(pwd)

        GENOME_ABS=$(readlink -f {input.genome})
        LOG_FILE_ABS=$(readlink -f {output.log_file})
        LOG_ABS=$(readlink -f {log})

        cd {params.outdir}

        gmes_petap.pl \
            --verbose \
            --ES \
            --sequence $GENOME_ABS \
            --cores {threads} \
            --min_contig {params.min_contig} \
            {params.fungus} \
            {params.gm_max_intergenic} \
            {params.gcode} \
            --soft_mask auto \
            --gc_donor {params.gc_donor} \
            1> $LOG_FILE_ABS \
            2> $LOG_ABS
        GM_EXIT=$?

        if [ $GM_EXIT -ne 0 ] || [ ! -f genemark.gtf ]; then
            echo "ERROR: GeneMark-ES failed (exit=$GM_EXIT), no genemark.gtf produced" >> $LOG_ABS
            echo "--- last 30 lines of gmes.log ---" >> $LOG_ABS
            tail -30 $LOG_FILE_ABS >> $LOG_ABS 2>/dev/null || true
            exit 1
        fi

        n_genes=$(grep -c $'\\tgene\\t' genemark.gtf || echo "0")
        echo "GeneMark-ES predicted $n_genes genes" >> $LOG_ABS

        if [ "$n_genes" -eq 0 ]; then
            echo "ERROR: GeneMark-ES produced genemark.gtf but predicted 0 genes." >> $LOG_ABS
            echo "       Cannot seed ProtHint with an empty gene set." >> $LOG_ABS
            echo "       Check genome quality, contig lengths (min_contig={params.min_contig})," >> $LOG_ABS
            echo "       and whether --fungus should be set." >> $LOG_ABS
            echo "--- last 30 lines of gmes.log ---" >> $LOG_ABS
            tail -30 $LOG_FILE_ABS >> $LOG_ABS 2>/dev/null || true
            exit 1
        fi

        # Record software versions
        cd $WORKDIR
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        GM_VER=$(grep -m1 '# GeneMark-ES Suite version' $(which gmes_petap.pl) 2>/dev/null | grep -oP 'version \K\S+' || true)
        GM_COMMIT=$(grep 'refs/remotes/origin/main' /opt/ETP/.git/packed-refs 2>/dev/null | awk '{{print substr($1,1,7)}}' || true)
        if [ -n "$GM_VER" ]; then
            ( flock 9; printf "GeneMark-ES\t%s (commit %s)\n" "$GM_VER" "$GM_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"
        else
            ( flock 9; printf "GeneMark-ES\tcommit %s\n" "$GM_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"
        fi

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite genemark_es "$REPORT_DIR"
        """
