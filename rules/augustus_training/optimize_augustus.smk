rule optimize_augustus:
    """
    Optimize AUGUSTUS parameters using k-fold cross-validation.

    This rule implements BRAKER's parameter optimization strategy using
    optimize_augustus.pl. The process:

    1. Splits training set into train/test sets using randomSplit.pl
       - First split: creates test set for final accuracy measurement
       - Second split: creates train.train and train.test for optimization

    2. Runs optimize_augustus.pl with k-fold cross-validation
       - Default: 8-fold cross-validation (BRAKER standard)
       - Adjusts k based on available CPUs (each bucket needs >= 200 genes)
       - Uses --onlytrain for training, test set for validation
       - Optimizes exon_probs, parameters, and weightmatrix

    3. Runs final etraining on full training set with optimized parameters

    Test set sizes (following BRAKER logic):
    - < 600 genes: 200 for each test set
    - 600-1000 genes: 200 for each test set
    - > 1000 genes: 300 for each test set

    Resources:
    - Can use full node CPUs (optimize_augustus.pl parallelizes well)
    - Each fold should have >= 200 genes for robust optimization
    - Memory: standard allocation per CPU

    Input:
        gb_train: Training set after train/test split (test set already held out)
        gb_test: Test set for final accuracy measurement (held out, never used for training)
        training_log: Ensures etraining/stop codon setting is complete
        species_params_backup: Backup of parameters before modification

    Output:
        gb_train_train: Training subset for optimize_augustus.pl
        gb_train_test: Test subset for optimize_augustus.pl validation
        optimize_log: Summary of optimization process
        training_summary: Summary of training gene counts
        accuracy_after_training: Accuracy after initial training
        accuracy_after_optimize: Accuracy after optimization
    """
    input:
        gb_train = "output/{sample}/train.gb.train",
        gb_test = "output/{sample}/train.gb.test",
        training_log = "output/{sample}/augustus_training.log",
        species_params_backup = "augustus_config/species/{sample}_rescuedB3/{sample}_rescuedB3_parameters.cfg.backup"
    output:
        gb_train_train = "output/{sample}/train.gb.train.train",
        gb_train_test = "output/{sample}/train.gb.train.test",
        optimize_log = "output/{sample}/optimize_augustus.log",
        training_summary = "output/{sample}/training_gene_summary.txt",
        accuracy_after_training = "output/{sample}/accuracy_after_training.txt",
        accuracy_after_optimize = "output/{sample}/accuracy_after_optimize.txt"
    benchmark:
        "benchmarks/{sample}/optimize_augustus/optimize_augustus.txt"
    params:
        species_name = lambda w: get_species_name(w),
        output_dir = lambda w: get_output_dir(w),
        aug_config = augustus_config_path,
        skip_optimize = config['skip_optimize_augustus'],
        rounds = 5,  # Default number of optimization rounds (BRAKER default)
        use_dev_shm = config['use_dev_shm'],
        dev_shm_path = lambda w: f"/dev/shm/{w.sample}/optimize_augustus" if config['use_dev_shm'] else get_output_dir(w),
        username = config['username']
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        # Set AUGUSTUS_CONFIG_PATH
        export AUGUSTUS_CONFIG_PATH={params.aug_config}

        echo "[INFO] ========== AUGUSTUS PARAMETER OPTIMIZATION ==========" | tee {output.optimize_log}
        echo "[INFO] AUGUSTUS_CONFIG_PATH: $AUGUSTUS_CONFIG_PATH" | tee -a {output.optimize_log}
        echo "[INFO] Species: {params.species_name}" | tee -a {output.optimize_log}
        echo "[INFO] Skip optimization: {params.skip_optimize}" | tee -a {output.optimize_log}
        echo "[INFO] Training set: {input.gb_train}" | tee -a {output.optimize_log}
        echo "[INFO] Test set: {input.gb_test}" | tee -a {output.optimize_log}
        echo "[INFO] Use /dev/shm: {params.use_dev_shm}" | tee -a {output.optimize_log}

        # Setup temporary storage for optimize_augustus buckets
        if [ "{params.use_dev_shm}" = "True" ]; then
            OPTIMIZE_TMP="{params.dev_shm_path}"
            echo "[INFO] Using /dev/shm for temporary storage: $OPTIMIZE_TMP" | tee -a {output.optimize_log}

            # Clean up any leftover data from previous runs
            if [ -d "$OPTIMIZE_TMP" ]; then
                echo "[INFO] Cleaning up leftover data from previous run in /dev/shm..." | tee -a {output.optimize_log}
                rm -rf "$OPTIMIZE_TMP"
            fi

            mkdir -p "$OPTIMIZE_TMP"

            # Setup cleanup trap that runs on EXIT (success, error, or signal)
            trap "echo '[INFO] Cleaning up /dev/shm...'; rm -rf $OPTIMIZE_TMP; echo '[INFO] /dev/shm cleanup completed'" EXIT
        else
            echo "[INFO] Using regular storage for temporary files" | tee -a {output.optimize_log}
            OPTIMIZE_TMP="{params.output_dir}"
        fi

        # Count genes in training and test sets (already split by split_training_set rule)
        GENES_TEST=$(grep -c "^LOCUS" {input.gb_test} || echo 0)
        GENES_TRAIN=$(grep -c "^LOCUS" {input.gb_train} || echo 0)
        echo "[INFO] Test set genes: $GENES_TEST (held out for final accuracy)" | tee -a {output.optimize_log}
        echo "[INFO] Training set genes: $GENES_TRAIN (will be further split for optimization)" | tee -a {output.optimize_log}

        # Determine test set size for optimization split based on BRAKER logic
        if [ $GENES_TRAIN -lt 600 ]; then
            # For small datasets, use 1/3 of training genes for optimization test
            TESTSIZE2=$(($GENES_TRAIN / 3))

            # Check if we have enough genes
            REMAINING=$(($GENES_TRAIN - $TESTSIZE2))
            if [ $TESTSIZE2 -eq 0 ] || [ $REMAINING -eq 0 ]; then
                echo "[ERROR] Insufficient training genes ($GENES_TRAIN) for optimization split!" | tee -a {output.optimize_log}
                echo "[ERROR] Need at least 2 genes for optimization train/test split" | tee -a {output.optimize_log}
                exit 1
            fi
        elif [ $GENES_TRAIN -le 1000 ]; then
            TESTSIZE2=200
        else
            TESTSIZE2=300
        fi

        echo "[INFO] Optimization test set size: $TESTSIZE2 genes" | tee -a {output.optimize_log}

        # Split training set for optimization
        echo "[INFO] Creating train/test split for parameter optimization..." | tee -a {output.optimize_log}
        randomSplit.pl {input.gb_train} $TESTSIZE2 2>&1 | tee -a {output.optimize_log}
        # randomSplit.pl creates {input.gb_train}.test and {input.gb_train}.train
        # which are train.gb.train.test and train.gb.train.train

        GENES_TRAIN_TRAIN=$(grep -c "^LOCUS" {output.gb_train_train} || echo 0)
        GENES_TRAIN_TEST=$(grep -c "^LOCUS" {output.gb_train_test} || echo 0)
        echo "[INFO] Optimization training set: $GENES_TRAIN_TRAIN genes" | tee -a {output.optimize_log}
        echo "[INFO] Optimization test set: $GENES_TRAIN_TEST genes" | tee -a {output.optimize_log}

        # Step 2b: Run AUGUSTUS on test set to assess accuracy after initial training
        echo "[INFO] =======================================" | tee -a {output.optimize_log}
        echo "[INFO] Assessing AUGUSTUS accuracy after initial training (before optimization)..." | tee -a {output.optimize_log}

        augustus \
            --species={params.species_name} \
            --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
            {input.gb_test} \
            > {params.output_dir}/augustus_test_after_training.out 2>&1

        # Parse accuracy from AUGUSTUS output (values are fractions, need to multiply by 100 for percentages)
        NU_SEN=$(grep "^nucleotide level" {params.output_dir}/augustus_test_after_training.out | awk '{{printf "%.2f", $4 * 100}}')
        NU_SP=$(grep "^nucleotide level" {params.output_dir}/augustus_test_after_training.out | awk '{{printf "%.2f", $6 * 100}}')
        EX_SEN=$(grep "^exon level" {params.output_dir}/augustus_test_after_training.out | awk '{{printf "%.2f", $14 * 100}}')
        EX_SP=$(grep "^exon level" {params.output_dir}/augustus_test_after_training.out | awk '{{printf "%.2f", $16 * 100}}')
        GEN_SEN=$(grep "^gene level" {params.output_dir}/augustus_test_after_training.out | awk '{{printf "%.2f", $14 * 100}}')
        GEN_SP=$(grep "^gene level" {params.output_dir}/augustus_test_after_training.out | awk '{{printf "%.2f", $16 * 100}}')

        # Calculate target accuracy (BRAKER formula: weighted average)
        # target = (3*nu_sen + 2*nu_sp + 4*ex_sen + 3*ex_sp + 2*gen_sen + 1*gen_sp) / 15
        ACCURACY_AFTER_TRAINING=$(awk "BEGIN {{printf \"%.2f\", (3*$NU_SEN + 2*$NU_SP + 4*$EX_SEN + 3*$EX_SP + 2*$GEN_SEN + 1*$GEN_SP) / 15}}")

        echo "[INFO] Accuracy after initial training: $ACCURACY_AFTER_TRAINING%" | tee -a {output.optimize_log}
        echo "[INFO]   Nucleotide level - Sensitivity: $NU_SEN%, Specificity: $NU_SP%" | tee -a {output.optimize_log}
        echo "[INFO]   Exon level - Sensitivity: $EX_SEN%, Specificity: $EX_SP%" | tee -a {output.optimize_log}
        echo "[INFO]   Gene level - Sensitivity: $GEN_SEN%, Specificity: $GEN_SP%" | tee -a {output.optimize_log}

        # Save accuracy to file
        cat > {output.accuracy_after_training} <<ACCURACY1
AUGUSTUS Accuracy After Initial Training (before optimize_augustus.pl)

Nucleotide level:
  Sensitivity: $NU_SEN%
  Specificity: $NU_SP%

Exon level:
  Sensitivity: $EX_SEN%
  Specificity: $EX_SP%

Gene level:
  Sensitivity: $GEN_SEN%
  Specificity: $GEN_SP%

Target accuracy (weighted): $ACCURACY_AFTER_TRAINING%
Formula: (3*nu_sen + 2*nu_sp + 4*ex_sen + 3*ex_sp + 2*gen_sen + 1*gen_sp) / 15
ACCURACY1

        # Step 3: Run optimize_augustus.pl (unless skipped)
        if [ "{params.skip_optimize}" = "True" ]; then
            echo "[INFO] Skipping optimize_augustus.pl (skip_optimize_augustus=1)" | tee -a {output.optimize_log}
            echo "[INFO] This will save runtime but may reduce accuracy by 2-3%" | tee -a {output.optimize_log}

            # Create empty accuracy file indicating optimization was skipped
            cat > {output.accuracy_after_optimize} <<ACCURACY_SKIP
AUGUSTUS Accuracy After Optimization

NOTE: optimize_augustus.pl was SKIPPED (skip_optimize_augustus=1)
No optimization was performed, so there is no post-optimization accuracy to report.

The parameters from the initial etraining are being used without further optimization.
See accuracy_after_training.txt for the accuracy after initial training.
ACCURACY_SKIP
        else
            # Calculate k-fold based on available CPUs (each bucket needs >= 200 genes)
            K_FOLD=8
            if [ {threads} -gt 1 ]; then
                for ((i=1; i<={threads}; i++)); do
                    if [ $(($GENES_TRAIN_TRAIN / $i)) -ge 200 ]; then
                        K_FOLD=$i
                    fi
                done
            fi
            if [ $K_FOLD -lt 8 ]; then
                K_FOLD=8
            fi

            echo "[INFO] Using $K_FOLD-fold cross-validation with {threads} CPUs" | tee -a {output.optimize_log}
            echo "[INFO] Running optimize_augustus.pl (rounds={params.rounds})..." | tee -a {output.optimize_log}
            echo "[INFO] Working directory for buckets: $OPTIMIZE_TMP" | tee -a {output.optimize_log}

            # Convert paths to absolute before changing directory
            TRAIN_TRAIN_ABS=$(readlink -f {output.gb_train_train})
            TRAIN_TEST_ABS=$(readlink -f {output.gb_train_test})
            LOG_ABS=$(readlink -f {output.optimize_log})

            # Change to temporary directory so optimize_augustus.pl creates buckets there
            cd "$OPTIMIZE_TMP"

            # Run optimize_augustus.pl (redirect output to log)
            # Note: optimize_augustus.pl creates bucket directories in current working directory
            optimize_augustus.pl \
                --species={params.species_name} \
                --rounds={params.rounds} \
                --kfold=$K_FOLD \
                --cpus=$K_FOLD \
                --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
                --onlytrain=$TRAIN_TRAIN_ABS \
                $TRAIN_TEST_ABS \
                >> $LOG_ABS 2>&1

            # Return to original directory
            cd -

            echo "[INFO] Parameter optimization completed" | tee -a {output.optimize_log}

            # Step 4: Run final etraining on full training set with optimized parameters
            echo "[INFO] Running final etraining on full training set..." | tee -a {output.optimize_log}
            etraining \
                --species={params.species_name} \
                --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
                {input.gb_train} \
                >> {output.optimize_log} 2>&1

            echo "[INFO] Final etraining completed" | tee -a {output.optimize_log}

            # Step 5: Run AUGUSTUS on test set to assess accuracy after optimization
            echo "[INFO] =======================================" | tee -a {output.optimize_log}
            echo "[INFO] Assessing AUGUSTUS accuracy after optimization and final etraining..." | tee -a {output.optimize_log}

            augustus \
                --species={params.species_name} \
                --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
                {input.gb_test} \
                > {params.output_dir}/augustus_test_after_optimize.out 2>&1

            # Parse accuracy from AUGUSTUS output (values are fractions, need to multiply by 100 for percentages)
            NU_SEN_OPT=$(grep "^nucleotide level" {params.output_dir}/augustus_test_after_optimize.out | awk '{{printf "%.2f", $4 * 100}}')
            NU_SP_OPT=$(grep "^nucleotide level" {params.output_dir}/augustus_test_after_optimize.out | awk '{{printf "%.2f", $6 * 100}}')
            EX_SEN_OPT=$(grep "^exon level" {params.output_dir}/augustus_test_after_optimize.out | awk '{{printf "%.2f", $14 * 100}}')
            EX_SP_OPT=$(grep "^exon level" {params.output_dir}/augustus_test_after_optimize.out | awk '{{printf "%.2f", $16 * 100}}')
            GEN_SEN_OPT=$(grep "^gene level" {params.output_dir}/augustus_test_after_optimize.out | awk '{{printf "%.2f", $14 * 100}}')
            GEN_SP_OPT=$(grep "^gene level" {params.output_dir}/augustus_test_after_optimize.out | awk '{{printf "%.2f", $16 * 100}}')

            # Calculate target accuracy (BRAKER formula: weighted average)
            ACCURACY_AFTER_OPTIMIZE=$(awk "BEGIN {{printf \"%.2f\", (3*$NU_SEN_OPT + 2*$NU_SP_OPT + 4*$EX_SEN_OPT + 3*$EX_SP_OPT + 2*$GEN_SEN_OPT + 1*$GEN_SP_OPT) / 15}}")

            echo "[INFO] Accuracy after optimization: $ACCURACY_AFTER_OPTIMIZE%" | tee -a {output.optimize_log}
            echo "[INFO]   Nucleotide level - Sensitivity: $NU_SEN_OPT%, Specificity: $NU_SP_OPT%" | tee -a {output.optimize_log}
            echo "[INFO]   Exon level - Sensitivity: $EX_SEN_OPT%, Specificity: $EX_SP_OPT%" | tee -a {output.optimize_log}
            echo "[INFO]   Gene level - Sensitivity: $GEN_SEN_OPT%, Specificity: $GEN_SP_OPT%" | tee -a {output.optimize_log}

            # Calculate improvement
            IMPROVEMENT=$(awk "BEGIN {{printf \"%.2f\", $ACCURACY_AFTER_OPTIMIZE - $ACCURACY_AFTER_TRAINING}}")
            echo "[INFO] Improvement from optimization: $IMPROVEMENT%" | tee -a {output.optimize_log}

            # Save accuracy to file
            cat > {output.accuracy_after_optimize} <<ACCURACY2
AUGUSTUS Accuracy After Optimization (after optimize_augustus.pl and final etraining)

Nucleotide level:
  Sensitivity: $NU_SEN_OPT%
  Specificity: $NU_SP_OPT%

Exon level:
  Sensitivity: $EX_SEN_OPT%
  Specificity: $EX_SP_OPT%

Gene level:
  Sensitivity: $GEN_SEN_OPT%
  Specificity: $GEN_SP_OPT%

Target accuracy (weighted): $ACCURACY_AFTER_OPTIMIZE%
Formula: (3*nu_sen + 2*nu_sp + 4*ex_sen + 3*ex_sp + 2*gen_sen + 1*gen_sp) / 15

Improvement from initial training: $IMPROVEMENT%
ACCURACY2
        fi

        # Summary
        echo "[INFO] =======================================" | tee -a {output.optimize_log}
        echo "[INFO] AUGUSTUS optimization complete" | tee -a {output.optimize_log}
        echo "[INFO] Final training set: {input.gb_train} ($GENES_TRAIN genes)" | tee -a {output.optimize_log}
        echo "[INFO] Final test set: {input.gb_test} ($GENES_TEST genes)" | tee -a {output.optimize_log}
        echo "[INFO] Optimized parameters in: $AUGUSTUS_CONFIG_PATH/species/{params.species_name}/" | tee -a {output.optimize_log}

        # Create training gene summary file for easy reference
        cat > {output.training_summary} <<SUMMARY
=== AUGUSTUS Training Gene Summary ===

Test set (held out): $GENES_TEST genes
Training set (used for etraining): $GENES_TRAIN genes

Training split for optimize_augustus.pl:
  Optimization training set: $GENES_TRAIN_TRAIN genes
  Optimization test set: $GENES_TRAIN_TEST genes

Final etraining used: $GENES_TRAIN genes (from {input.gb_train})

Note: The final etraining run (after optimize_augustus.pl) uses the full
training set ({input.gb_train}) with optimized parameters.
The test set ({input.gb_test}) is NEVER used for training, only for accuracy measurement.
SUMMARY

        echo "[INFO] Training gene summary written to {output.training_summary}" | tee -a {output.optimize_log}

        # Citations
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite augustus "$REPORT_DIR"
        """
