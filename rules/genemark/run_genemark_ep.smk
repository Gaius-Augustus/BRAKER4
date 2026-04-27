"""
Run GeneMark-EP for gene prediction with protein hints.

GeneMark-EP uses protein-derived intron hints from ProtHint to
train itself and predict genes. This is the protein-evidence
counterpart to GeneMark-ET (RNA-Seq evidence).

Mirrors braker.pl's GeneMark_EP() function.

Input:
    - Masked genome FASTA
    - Protein intron hints (from prepare_genemark_hints_ep)
    - Evidence GFF (src=M hints, optional)

Output:
    - genemark.gtf - Gene predictions trained with protein evidence

Container: teambraker/braker3:latest (contains gmes_petap.pl)
"""

rule run_genemark_ep:
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
        hints="output/{sample}/genemark_hints_ep.gff",
        evidence="output/{sample}/genemark_evidence.gff"
    output:
        gtf="output/{sample}/GeneMark-EP/genemark.gtf",
        log_file="output/{sample}/GeneMark-EP/gmes.log"
    log:
        "logs/{sample}/genemark_ep/genemark_ep.log"
    benchmark:
        "benchmarks/{sample}/genemark_ep/genemark_ep.txt"
    threads: workflow.cores
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/GeneMark-EP",
        min_contig=config.get("min_contig", 10000),
        fungus="--fungus" if config.get("fungus", False) else "",
        gm_max_intergenic=f"--max_intergenic {config.get('gm_max_intergenic')}" if config.get("gm_max_intergenic") else "",
        gcode=f"--gcode {config.get('translation_table')}" if config.get("translation_table", 1) != 1 else "",
        gc_donor=config.get("gc_donor", 0.001)
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}
        WORKDIR=$(pwd)

        GENOME_ABS=$(readlink -f {input.genome})
        HINTS_ABS=$(readlink -f {input.hints})
        EVIDENCE_ABS=$(readlink -f {input.evidence})
        LOG_FILE_ABS=$(readlink -f {output.log_file})
        LOG_ABS=$(readlink -f {log})

        cd {params.outdir}

        echo "Running GeneMark-EP..." > $LOG_ABS

        GMES_CORES=$(python3 -c "txt=open('/proc/cpuinfo').read(); c=[l.split(':')[-1].strip() for l in txt.splitlines() if l.startswith('cpu cores')]; s=set(l.split(':')[-1].strip() for l in txt.splitlines() if l.startswith('physical id')); total=int(c[0])*max(1,len(s)) if c else 0; print(min({threads},total) if 0<total<{threads} else {threads})" 2>/dev/null || echo {threads})

        # Build command
        CMD="gmes_petap.pl --verbose"
        CMD="$CMD --seq $GENOME_ABS"
        CMD="$CMD --EP $HINTS_ABS"
        CMD="$CMD --cores $GMES_CORES"
        CMD="$CMD --min_contig {params.min_contig}"
        CMD="$CMD --gc_donor {params.gc_donor}"

        # Add evidence file if non-empty
        if [ -s "$EVIDENCE_ABS" ]; then
            CMD="$CMD --evidence $EVIDENCE_ABS"
            echo "Using evidence GFF: $EVIDENCE_ABS" >> $LOG_ABS
        fi

        # Optional flags
        if [ -n "{params.fungus}" ]; then
            CMD="$CMD {params.fungus}"
        fi

        if [ -n "{params.gm_max_intergenic}" ]; then
            CMD="$CMD {params.gm_max_intergenic}"
        fi

        if [ -n "{params.gcode}" ]; then
            CMD="$CMD {params.gcode}"
        fi

        CMD="$CMD --soft_mask auto"

        echo "Command: $CMD" >> $LOG_ABS

        eval $CMD 1>> $LOG_FILE_ABS 2>> $LOG_ABS

        if [ ! -f genemark.gtf ]; then
            echo "ERROR: GeneMark-EP failed to produce genemark.gtf" >> $LOG_ABS
            exit 1
        fi

        n_genes=$(awk '$3=="gene"{{c++}}END{{print c+0}}' genemark.gtf)
        echo "GeneMark-EP predicted $n_genes genes" >> $LOG_ABS

        # Record software versions
        cd $WORKDIR
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        GM_VER=$(grep -m1 '# GeneMark-ES Suite version' $(which gmes_petap.pl) 2>/dev/null | grep -oP 'version \K\S+' || true)
        GM_COMMIT=$(grep 'refs/remotes/origin/main' /opt/ETP/.git/packed-refs 2>/dev/null | awk '{{print substr($1,1,7)}}' || true)
        if [ -n "$GM_VER" ]; then
            ( flock 9; printf "GeneMark-EP+\t%s (commit %s)\n" "$GM_VER" "$GM_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"
        else
            ( flock 9; printf "GeneMark-EP+\tcommit %s\n" "$GM_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"
        fi

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite genemark_ep "$REPORT_DIR"
        cite braker2 "$REPORT_DIR"
        """
