rule split_training_set:
    """
    Split training set into train and test portions for accuracy assessment.

    This rule performs the critical first split that BRAKER does:
    - Creates train.gb.test: held-out test set for final accuracy measurement
    - Creates train.gb.train: training set used for all etraining and optimization

    The test set is NEVER used during training. It's only used to assess
    the accuracy of the trained AUGUSTUS parameters. This ensures an unbiased
    accuracy measurement.

    Test set sizes follow BRAKER logic:
    - < 600 genes: 1/3 of genes
    - 600-1000 genes: 200 genes
    - > 1000 genes: 300 genes

    This is a local rule (runs on head node) because randomSplit.pl is very fast.

    Input:
        gb_input: Clean GenBank file (after etraining filtering, optionally downsampled)

    Output:
        gb_test: Test set for accuracy measurement (never used for training)
        gb_train: Training set (used for all etraining and optimization)
        split_log: Log file documenting the split
    """
    input:
        gb_input = "output/{sample}/bonafide.f.clean.d.gb"
    output:
        gb_test = "output/{sample}/train.gb.test",
        gb_train = "output/{sample}/train.gb.train",
        split_log = "output/{sample}/training_set_split.log"
    benchmark:
        "benchmarks/{sample}/split_training_set/split_training_set.txt"
    container:
        BRAKER3_CONTAINER
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ========== TRAINING SET SPLIT ==========" | tee {output.split_log}

        # Count total genes in input
        TOTAL_GENES=$(grep -c "^LOCUS" {input.gb_input} || echo 0)
        echo "[INFO] Total genes in clean training set: $TOTAL_GENES" | tee -a {output.split_log}

        # Determine test set size based on BRAKER logic
        if [ $TOTAL_GENES -lt 600 ]; then
            # For small datasets, use 1/3 of genes for test set
            TESTSIZE=$(($TOTAL_GENES / 3))

            # Check if we have enough genes (need at least 2 genes for train and test)
            REMAINING=$(($TOTAL_GENES - $TESTSIZE))
            if [ $TESTSIZE -eq 0 ] || [ $REMAINING -eq 0 ]; then
                echo "[ERROR] Insufficient training genes ($TOTAL_GENES) for train/test split!" | tee -a {output.split_log}
                echo "[ERROR] Need at least 2 genes for train/test split" | tee -a {output.split_log}
                exit 1
            fi
        elif [ $TOTAL_GENES -le 1000 ]; then
            TESTSIZE=200
        else
            TESTSIZE=300
        fi

        echo "[INFO] Test set size: $TESTSIZE genes" | tee -a {output.split_log}
        echo "[INFO] Training set size: $(($TOTAL_GENES - $TESTSIZE)) genes" | tee -a {output.split_log}

        # Perform the split using randomSplit.pl
        echo "[INFO] Splitting training set..." | tee -a {output.split_log}
        randomSplit.pl {input.gb_input} $TESTSIZE 2>&1 | tee -a {output.split_log}

        # Move the split files to their final locations
        mv {input.gb_input}.test {output.gb_test}
        mv {input.gb_input}.train {output.gb_train}

        # Count genes in each split
        GENES_TEST=$(grep -c "^LOCUS" {output.gb_test} || echo 0)
        GENES_TRAIN=$(grep -c "^LOCUS" {output.gb_train} || echo 0)

        echo "[INFO] =======================================" | tee -a {output.split_log}
        echo "[INFO] Split complete:" | tee -a {output.split_log}
        echo "[INFO]   Test set: $GENES_TEST genes (held out for accuracy measurement)" | tee -a {output.split_log}
        echo "[INFO]   Training set: $GENES_TRAIN genes (will be used for etraining)" | tee -a {output.split_log}
        echo "[INFO] =======================================" | tee -a {output.split_log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        """
