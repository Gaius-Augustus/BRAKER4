rule run_augustus_hints:
    """
    Run AUGUSTUS with hints using data parallelization.

    This rule implements BRAKER's approach to AUGUSTUS prediction with hints:
    1. Split genome into chunks for parallel processing
    2. Create AUGUSTUS job list with createAugustusJoblist.pl
    3. Run AUGUSTUS jobs in parallel (using Snakemake's threads)
    4. Join predictions with join_aug_pred.pl
    5. Convert GFF to GTF

    The rule uses extrinsic.cfg to weight different hint sources:
    - M: Manual anchors (highest weight)
    - P: Protein alignments (high weight)
    - E: EST/RNA-Seq alignments (medium weight)
    - C: Combined protein/EST evidence

    AUGUSTUS parameters:
    - --alternatives-from-evidence=true: Generate alternative transcripts
    - --allow_hinted_splicesites=gcag,atac: Allow non-canonical splice sites
    - --UTR=off: No UTR prediction (we only rescue CDS)
    - --exonnames=on: Add exon names to output
    - --codingseq=on: Include coding sequence in output

    Resources:
        - Uses all available CPUs for parallel execution
        - Full node memory allocation
        - Submitted to SLURM cluster

    Input:
        genome: Genome FASTA file
        hintsfile: Merged hints from all sources
        optimize_log: Ensures AUGUSTUS parameters are optimized before prediction

    Output:
        augustus_gff: AUGUSTUS predictions in GFF format
        augustus_gtf: AUGUSTUS predictions in GTF format
        job_lst: Job list for parallel execution (for debugging)
    """
    input:
        genome = lambda w: get_masked_genome(w.sample),
        hintsfile = "output/{sample}/hintsfile.gff",
        optimize_log = "output/{sample}/optimize_augustus.log"
    output:
        job_lst = "output/{sample}/augustus_hints.job.lst",
        augustus_gff = "output/{sample}/augustus.hints.gff",
        augustus_gtf = "output/{sample}/augustus.hints.gtf"
    benchmark:
        "benchmarks/{sample}/run_augustus_hints/run_augustus_hints.txt"
    params:
        species_name = lambda w: get_species_name(w),
        output_dir = lambda w: get_output_dir(w),
        aug_config = augustus_config_path,
        extrinsic_cfg = lambda w: get_extrinsic_cfg(w),
        chunksize = config['augustus_chunksize'],
        overlap = config['augustus_overlap'],
        use_dev_shm = config['use_dev_shm'],
        dev_shm_path = lambda w: f"/dev/shm/{w.sample}/augustus_hints" if config['use_dev_shm'] else get_output_dir(w),
        username = config['username'],
        allow_hinted_splicesites = config.get('allow_hinted_splicesites', 'gcag,atac')
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

        echo "[INFO] ===== AUGUSTUS PREDICTION WITH HINTS ====="
        echo "[INFO] Species: {params.species_name}"
        echo "[INFO] Threads: {threads}"
        echo "[INFO] Chunksize: {params.chunksize}"
        echo "[INFO] Use /dev/shm: {params.use_dev_shm}"

        # Setup temporary storage
        if [ "{params.use_dev_shm}" = "True" ]; then
            TMP_DIR="{params.dev_shm_path}/{wildcards.sample}"
            echo "[INFO] Using /dev/shm for temporary storage: $TMP_DIR"

            # Clean up any leftover data from previous runs
            if [ -d "$TMP_DIR" ]; then
                echo "[INFO] Cleaning up leftover data from previous run in /dev/shm..."
                rm -rf "$TMP_DIR"
            fi

            mkdir -p "$TMP_DIR"

            # Setup cleanup trap that runs on EXIT (success, error, or signal)
            trap "echo '[INFO] Cleaning up /dev/shm...'; rm -rf $TMP_DIR; echo '[INFO] /dev/shm cleanup completed'" EXIT

            GENOME_SPLIT_TMP="$TMP_DIR/genome_split"
            AUGUSTUS_TMP="$TMP_DIR/augustus_tmp"
        else
            echo "[INFO] Using regular storage for temporary files"
            GENOME_SPLIT_TMP="{params.output_dir}/genome_split"
            AUGUSTUS_TMP="{params.output_dir}/augustus_tmp"
        fi

        # Step 1: Split genome into chunks for parallel processing
        echo "[INFO] Splitting genome into chunks..."
        mkdir -p "$GENOME_SPLIT_TMP"

        splitMfasta.pl {input.genome} --outputpath="$GENOME_SPLIT_TMP" 2> {params.output_dir}/splitMfasta.err

        # Rename split files according to scaffold names (BRAKER convention).
        # splitMfasta.pl names its outputs after the input file stem, so with
        # {input.genome} = genome.fa the stem is "genome" but with
        # genome_masked.fa the stem is "genome_masked". Glob on *.split.*
        # so the rename works regardless of which stem is in play.
        cd "$GENOME_SPLIT_TMP"
        for f in *.split.*; do
            if [ -f "$f" ]; then
                NAME=$(grep ">" $f | head -1)
                mv $f ${{NAME#>}}.fa
            fi
        done
        cd -

        GENOME_FILES=$(ls "$GENOME_SPLIT_TMP" | wc -l)
        echo "[INFO] Split genome into $GENOME_FILES chunks"

        # Step 2: Create AUGUSTUS sequence list file
        echo "[INFO] Creating AUGUSTUS sequence list..."
        mkdir -p "$AUGUSTUS_TMP"

        # Convert relative paths to absolute
        GENOME_SPLIT_ABS=$(readlink -f "$GENOME_SPLIT_TMP")
        AUGUSTUS_TMP_ABS=$(readlink -f "$AUGUSTUS_TMP")
        HINTS_ABS=$(readlink -f {input.hintsfile})
        OUTPUT_DIR_ABS=$(readlink -f {params.output_dir})
        # Use BRAKER's ETP extrinsic config from container
        EXTRINSIC_ABS="{params.extrinsic_cfg}"

        # aug_hints.lst also goes to temp location (must be absolute path)
        if [ "{params.use_dev_shm}" = "True" ]; then
            AUG_LST="$TMP_DIR/aug_hints.lst"
        else
            AUG_LST="$AUGUSTUS_TMP_ABS/aug_hints.lst"
        fi

        echo "[INFO] Using absolute paths:"
        echo "[INFO]   Genome split: $GENOME_SPLIT_ABS"
        echo "[INFO]   AUGUSTUS tmp: $AUGUSTUS_TMP_ABS"
        echo "[INFO]   Hints file: $HINTS_ABS"
        echo "[INFO]   Extrinsic cfg: $EXTRINSIC_ABS"
        echo "[INFO]   Aug hints list: $AUG_LST"

        # Create aug_hints.lst file in BRAKER format:
        # Format: scaffold_file<TAB>hints_file<TAB>start<TAB>end
        # Get scaffold sizes from the genome
        echo "[INFO] Creating sequence list file..."
        > $AUG_LST  # Clear file

        for scaffold_file in $GENOME_SPLIT_ABS/*.fa; do
            if [ -f "$scaffold_file" ]; then
                # Get scaffold length from FASTA file (sum of all non-header lines)
                SCAFFOLD_SIZE=$(grep -v "^>" "$scaffold_file" | tr -d '\\n' | wc -c)
                echo -e "$scaffold_file\\t$HINTS_ABS\\t1\\t$SCAFFOLD_SIZE" >> $AUG_LST
            fi
        done

        NUM_SCAFFOLDS=$(wc -l < $AUG_LST)
        echo "[INFO] Created sequence list with $NUM_SCAFFOLDS scaffolds"

        # Step 3: Create AUGUSTUS job list using createAugustusJoblist.pl
        echo "[INFO] Creating AUGUSTUS job list..."

        # Job list and scripts go to temp location if using /dev/shm
        if [ "{params.use_dev_shm}" = "True" ]; then
            JOB_LST_TMP="$TMP_DIR/augustus_hints.job.lst"
        else
            JOB_LST_TMP="{output.job_lst}"
        fi
        JOB_LST_ABS=$(readlink -f "$JOB_LST_TMP" 2>/dev/null || realpath "$JOB_LST_TMP")

        # Change to temp directory to ensure job scripts are created there
        # createAugustusJoblist.pl writes scripts relative to current directory
        pushd "$AUGUSTUS_TMP_ABS" > /dev/null

        createAugustusJoblist.pl \
            --sequences=$AUG_LST \
            --wrap="#!/bin/bash" \
            --overlap={params.overlap} \
            --chunksize={params.chunksize} \
            --outputdir=$AUGUSTUS_TMP_ABS \
            --joblist=$JOB_LST_ABS \
            --jobprefix=aug_hints_ \
            --partitionHints \
            --command "augustus --species={params.species_name} --AUGUSTUS_CONFIG_PATH={params.aug_config} --extrinsicCfgFile=$EXTRINSIC_ABS --alternatives-from-evidence=true --UTR=off --exonnames=on --codingseq=on --allow_hinted_splicesites={params.allow_hinted_splicesites}" \
            2> $OUTPUT_DIR_ABS/createAugustusJoblist.err

        # Return to original directory
        popd > /dev/null

        NJOBS=$(wc -l < $JOB_LST_ABS)
        echo "[INFO] Created $NJOBS AUGUSTUS jobs in temp location"

        # Check if jobs were created
        if [ $NJOBS -eq 0 ]; then
            echo "[ERROR] No AUGUSTUS jobs were created!"
            echo "[ERROR] Check $OUTPUT_DIR_ABS/createAugustusJoblist.err for details"
            cat $OUTPUT_DIR_ABS/createAugustusJoblist.err
            exit 1
        fi

        # Copy job list to output for Snakemake (only if using /dev/shm)
        if [ "{params.use_dev_shm}" = "True" ]; then
            cp "$JOB_LST_ABS" {output.job_lst}
        fi

        # Step 3: Run AUGUSTUS jobs in parallel
        echo "[INFO] Running AUGUSTUS jobs in parallel (threads={threads})..."

        # Read job list and run jobs in parallel
        # Job scripts are in $AUGUSTUS_TMP_ABS (temp location)
        JOB_COUNTER=0
        RUNNING_JOBS=0

        while IFS= read -r job_script; do
            JOB_COUNTER=$((JOB_COUNTER + 1))

            # Job script paths in job list are relative to outputdir
            JOB_PATH="$AUGUSTUS_TMP_ABS/$job_script"

            if [ ! -f "$JOB_PATH" ]; then
                echo "[ERROR] Cannot find job script: $JOB_PATH"
                continue
            fi

            # Run job in background
            bash "$JOB_PATH" &
            RUNNING_JOBS=$((RUNNING_JOBS + 1))

            # Wait if we've reached the thread limit
            if [ $RUNNING_JOBS -ge {threads} ]; then
                wait -n  # Wait for any job to complete
                RUNNING_JOBS=$((RUNNING_JOBS - 1))
            fi

            # Progress reporting (every 10 jobs)
            if [ $((JOB_COUNTER % 10)) -eq 0 ]; then
                echo "[INFO] Processed $JOB_COUNTER/$NJOBS jobs..."
            fi
        done < $JOB_LST_ABS

        # Wait for all remaining jobs to complete
        wait
        echo "[INFO] All AUGUSTUS jobs completed"

        # Step 4: Join AUGUSTUS predictions
        echo "[INFO] Joining AUGUSTUS predictions..."

        # Concatenate all GFF files from temp location
        for chr_dir in "$AUGUSTUS_TMP_ABS"/*.gff; do
            if [ -f "$chr_dir" ]; then
                cat "$chr_dir" >> {output.augustus_gff}.tmp
            fi
        done

        # Run join_aug_pred.pl to merge overlapping predictions
        cat {output.augustus_gff}.tmp | join_aug_pred.pl > {output.augustus_gff}
        rm {output.augustus_gff}.tmp

        GENES_PREDICTED=$(grep -c "^[^#]" {output.augustus_gff} || echo 0)
        echo "[INFO] Predicted genes in GFF: $GENES_PREDICTED lines"

        # Step 5: Convert GFF to GTF
        echo "[INFO] Converting GFF to GTF..."
        # AUGUSTUS GFF output is already in GTF-compatible format, just filter and clean
        grep -v "^#" {output.augustus_gff} | \
            perl -ne 'if(m/\tAUGUSTUS\t/) {{print $_;}}' | \
            grep -v "^$" > {output.augustus_gtf}

        GENES_GTF=$(grep -cP '\tgene\t' {output.augustus_gtf} || echo 0)
        echo "[INFO] Genes in GTF: $GENES_GTF"

        echo "[INFO] ======================================="
        echo "[INFO] AUGUSTUS prediction completed successfully"
        echo "[INFO] Output GFF: {output.augustus_gff}"
        echo "[INFO] Output GTF: {output.augustus_gtf}"

        # Cleanup empty log files
        for logfile in $OUTPUT_DIR_ABS/splitMfasta.err $OUTPUT_DIR_ABS/createAugustusJoblist.err; do
            if [ -f "$logfile" ] && [ ! -s "$logfile" ]; then
                rm "$logfile"
            fi
        done

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        AUG_VER=$(augustus --version 2>&1 | head -1 | grep -oP '\([\d.]+\)' | tr -d '()' || true)
        AUG_COMMIT=$(grep 'refs/remotes/origin/master' /opt/Augustus/.git/packed-refs 2>/dev/null | awk '{{print substr($1,1,7)}}' || true)
        ( flock 9; printf "AUGUSTUS\t%s (commit %s)\n" "$AUG_VER" "$AUG_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite augustus "$REPORT_DIR"
        """
