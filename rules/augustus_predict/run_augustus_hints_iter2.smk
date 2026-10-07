"""
Run AUGUSTUS prediction iteration 2 for EP mode.

After ProtHint iteration 2 produces improved protein hints, AUGUSTUS
is re-run with the updated hintsfile. This reuses the genome split
and AUGUSTUS parameters from iteration 1 — only the hints change.

This rule mirrors the second augustus("off") call in braker.pl's EP mode.

Input:
    - Genome FASTA
    - Updated hintsfile from ProtHint iter2
    - AUGUSTUS training (optimize_log ensures training is done)

Output:
    - augustus.hints_iter2.gff/gtf: AUGUSTUS predictions with improved hints
"""

rule run_augustus_hints_iter2:
    """AUGUSTUS prediction with the iteration-2 hints (EP mode).

    Scratch: the genome split, the per-chunk hints, job scripts and GFFs
    are on the node-local disk (scripts/tmp_dir.sh, [paths] tmp_dir); only
    the joined predictions reach output/<sample>/. NEED 3 x (genome +
    hintsfile) + 5 GB; with less free the job works in
    output/<sample>/augustus_work_iter2/ and removes it at the end.
    """
    input:
        genome=lambda w: get_masked_genome(w.sample),
        hintsfile="output/{sample}/hintsfile_iter2.gff",
        optimize_log="output/{sample}/optimize_augustus.log"
    output:
        augustus_gff="output/{sample}/augustus.hints_iter2.gff",
        augustus_gtf="output/{sample}/augustus.hints_iter2.gtf"
    benchmark:
        "benchmarks/{sample}/run_augustus_hints_iter2/run_augustus_hints_iter2.txt"
    params:
        species_name=lambda w: get_species_name(w),
        output_dir=lambda w: get_output_dir(w),
        aug_config=augustus_config_path,
        extrinsic_cfg=lambda w: get_extrinsic_cfg(w),
        chunksize=config['augustus_chunksize'],
        overlap=config['augustus_overlap'],
        # Note: iteration 2 deliberately does NOT support /dev/shm. The
        # iter1 rule needs the speedup because it processes the entire
        # genome from scratch; iter2 only refines an already-completed
        # prediction, so the disk-vs-RAM difference is negligible.
        allow_hinted_splicesites=config.get('allow_hinted_splicesites', 'gcag,atac'),
        tmp_root=TMP_ROOT
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        export AUGUSTUS_CONFIG_PATH={params.aug_config}

        echo "[INFO] ===== AUGUSTUS PREDICTION ITERATION 2 (EP mode) ====="
        echo "[INFO] Using updated hints from ProtHint iteration 2"

        # Temporary storage: node-local scratch, removed on exit;
        # output/<sample>/augustus_work_iter2 when short of space.
        source {script_dir}/tmp_dir.sh
        WORKDIR=$PWD
        AUG_GFF_ABS=$WORKDIR/{output.augustus_gff}
        scratch_dir outDir "aughints2_{wildcards.sample}" "{params.tmp_root}" \
            "$(need_gb 3 {input.genome} {input.hintsfile})" "$WORKDIR/{params.output_dir}/augustus_work_iter2"
        trap 'rm -rf -- "$SCRATCH"' EXIT

        # Split genome (same as iter1 but into iter2 temp dir)
        GENOME_SPLIT="$outDir/genome_split"
        AUGUSTUS_TMP="$outDir/augustus_tmp"
        if [ -z "$SCRATCH" ]; then
            # Leftovers of an earlier attempt would be joined into the result
            rm -rf "$GENOME_SPLIT" "$AUGUSTUS_TMP"
        fi
        mkdir -p "$GENOME_SPLIT" "$AUGUSTUS_TMP"

        splitMfasta.pl {input.genome} --outputpath="$GENOME_SPLIT" 2>/dev/null

        # Rename split files. Glob on *.split.* because splitMfasta.pl
        # names chunks after the input file stem (genome vs genome_masked).
        cd "$GENOME_SPLIT"
        for f in *.split.*; do
            if [ -f "$f" ]; then
                NAME=$(grep ">" $f | head -1)
                mv $f ${{NAME#>}}.fa
            fi
        done
        cd -

        # Create sequence list
        GENOME_SPLIT_ABS=$(readlink -f "$GENOME_SPLIT")
        AUGUSTUS_TMP_ABS=$(readlink -f "$AUGUSTUS_TMP")
        HINTS_ABS=$(readlink -f {input.hintsfile})
        EXTRINSIC_ABS="{params.extrinsic_cfg}"

        AUG_LST="$AUGUSTUS_TMP_ABS/aug_hints_iter2.lst"
        > $AUG_LST
        for scaffold_file in $GENOME_SPLIT_ABS/*.fa; do
            if [ -f "$scaffold_file" ]; then
                SCAFFOLD_SIZE=$(grep -v "^>" "$scaffold_file" | tr -d '\n' | wc -c)
                echo -e "$scaffold_file\t$HINTS_ABS\t1\t$SCAFFOLD_SIZE" >> $AUG_LST
            fi
        done

        # Create job list
        pushd "$AUGUSTUS_TMP_ABS" > /dev/null
        createAugustusJoblist.pl \
            --sequences=$AUG_LST \
            --wrap="#!/bin/bash" \
            --overlap={params.overlap} \
            --chunksize={params.chunksize} \
            --outputdir=$AUGUSTUS_TMP_ABS \
            --joblist=$AUGUSTUS_TMP_ABS/aug_iter2.job.lst \
            --jobprefix=aug_iter2_ \
            --partitionHints \
            --command "augustus --species={params.species_name} --AUGUSTUS_CONFIG_PATH={params.aug_config} --extrinsicCfgFile=$EXTRINSIC_ABS --alternatives-from-evidence=true --UTR=off --exonnames=on --codingseq=on --allow_hinted_splicesites={params.allow_hinted_splicesites} --softmasking=1" \
            2>/dev/null
        popd > /dev/null

        NJOBS=$(wc -l < $AUGUSTUS_TMP_ABS/aug_iter2.job.lst)
        echo "[INFO] Created $NJOBS AUGUSTUS jobs in temp location"

        # Bail out loudly if createAugustusJoblist.pl produced no jobs.
        # Without this check, iter2 silently runs zero jobs and produces an
        # empty augustus.hints_iter2.gtf, which would surface as a misleading
        # "empty input" error in the next downstream rule. iter1 has the same
        # check (run_augustus_hints.smk:199-204).
        if [ $NJOBS -eq 0 ]; then
            echo "[ERROR] No AUGUSTUS jobs were created in iteration 2!"
            echo "[ERROR] createAugustusJoblist.pl produced an empty job list."
            echo "[ERROR] Inspect $AUGUSTUS_TMP_ABS/aug_hints_iter2.lst for sequence list problems."
            exit 1
        fi

        echo "[INFO] Running $NJOBS AUGUSTUS jobs (iteration 2)..."

        # Run jobs in parallel
        RUNNING=0
        while IFS= read -r job_script; do
            bash "$AUGUSTUS_TMP_ABS/$job_script" &
            RUNNING=$((RUNNING + 1))
            if [ $RUNNING -ge {threads} ]; then
                wait -n
                RUNNING=$((RUNNING - 1))
            fi
        done < $AUGUSTUS_TMP_ABS/aug_iter2.job.lst
        wait

        # Join predictions
        AUG_GFF_TMP="$outDir/augustus.hints_iter2.gff.tmp"
        rm -f "$AUG_GFF_TMP"
        for gff in "$AUGUSTUS_TMP_ABS"/*.gff; do
            [ -f "$gff" ] && cat "$gff" >> "$AUG_GFF_TMP"
        done
        cat "$AUG_GFF_TMP" | join_aug_pred.pl > "$AUG_GFF_ABS"
        rm "$AUG_GFF_TMP"

        # Convert to GTF
        grep -v "^#" {output.augustus_gff} | \
            perl -ne 'if(m/\tAUGUSTUS\t/) {{print $_;}}' | \
            grep -v "^$" > {output.augustus_gtf}

        GENES=$(awk '$3=="gene"{{n++}}END{{print n+0}}' {output.augustus_gtf})
        echo "[INFO] Iteration 2: $GENES genes predicted"

        # Cleanup (fallback; on scratch the EXIT trap removes the work dir)
        if [ -z "$SCRATCH" ]; then
            rm -rf "$outDir"
        fi

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite augustus "$REPORT_DIR"
        """
