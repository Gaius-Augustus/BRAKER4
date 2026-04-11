rule train_augustus:
    """
    Train AUGUSTUS parameters using etraining.

    This rule implements the BRAKER approach for training AUGUSTUS:
    1. Set stopCodonExcludedFromCDS to true in species parameters
    2. Run etraining on the clean training set
    3. Check error rate for "exon doesn't end in stop codon" errors
    4. If error rate >= 0.5, set stopCodonExcludedFromCDS to false and re-run
    5. Extract stop codon frequencies from etraining output
    6. Update species parameters with observed stop codon frequencies

    This ensures AUGUSTUS is properly trained with the correct stop codon
    handling and frequencies for the training data.

    Input:
        gb_train: Training set after train/test split (test set held out for accuracy)
        species_marker: Ensures AUGUSTUS species has been created
        split_log: Ensures training set has been split

    Output:
        training_stdout: etraining stdout with stop codon frequencies
        training_stderr: etraining stderr for error checking
        training_log: Log file with training statistics
        species_params_backup: Backup of original parameters before modification
    """
    input:
        gb_train = "output/{sample}/train.gb.train",
        species_marker = "augustus_config/species/{sample}_rescuedB3/.species_created",
        split_log = "output/{sample}/training_set_split.log"
    output:
        training_stdout = "output/{sample}/etraining.train.stdout",
        training_stderr = "output/{sample}/etraining.train.stderr",
        training_log = "output/{sample}/augustus_training.log",
        species_params_backup = "augustus_config/species/{sample}_rescuedB3/{sample}_rescuedB3_parameters.cfg.backup"
    benchmark:
        "benchmarks/{sample}/train_augustus/train_augustus.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        species_name = lambda w: get_species_name(w),
        aug_config = augustus_config_path,
        species_params = lambda w: f"{augustus_config_path}/species/{get_species_name(w)}/{get_species_name(w)}_parameters.cfg",
        translation_table = config.get("translation_table", 1)
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        # Set AUGUSTUS_CONFIG_PATH to our local copy
        export AUGUSTUS_CONFIG_PATH={params.aug_config}

        echo "[INFO] ========== AUGUSTUS TRAINING ==========" | tee {output.training_log}
        echo "[INFO] AUGUSTUS_CONFIG_PATH: $AUGUSTUS_CONFIG_PATH" | tee -a {output.training_log}
        echo "[INFO] Species: {params.species_name}" | tee -a {output.training_log}
        echo "[INFO] Parameters file: {params.species_params}" | tee -a {output.training_log}

        # Backup original parameters
        cp {params.species_params} {output.species_params_backup}
        echo "[INFO] Backed up original parameters" | tee -a {output.training_log}

        # Count total training genes
        TOTAL_GENES=$(grep -c "^LOCUS" {input.gb_train} || echo 0)
        echo "[INFO] Total training genes (test set held out): $TOTAL_GENES" | tee -a {output.training_log}

        # Step 1: Set stopCodonExcludedFromCDS to true
        echo "[INFO] Setting stopCodonExcludedFromCDS to true" | tee -a {output.training_log}
        sed -i 's/^\(stopCodonExcludedFromCDS\s\+\)\S\+/\1true/' {params.species_params}

        # Step 2: Run etraining with stopCodonExcludedFromCDS=true
        echo "[INFO] Running etraining (first pass, stopCodonExcludedFromCDS=true)..." | tee -a {output.training_log}
        etraining --species={params.species_name} --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
            {input.gb_train} \
            1> {output.training_stdout} \
            2> {output.training_stderr} || true

        # Step 3: Check error rate for missing stop codons
        MISSING_STOP_ERRORS=$(grep -c "exon doesn't end in stop codon" {output.training_stderr} || echo 0)

        # Avoid division by zero
        if [ "$TOTAL_GENES" -eq 0 ]; then
            ERROR_RATE="0.0000"
        else
            ERROR_RATE=$(awk "BEGIN {{printf \"%.4f\", $MISSING_STOP_ERRORS/$TOTAL_GENES}}")
        fi

        echo "[INFO] Missing stop codon errors: $MISSING_STOP_ERRORS" | tee -a {output.training_log}
        echo "[INFO] Error rate: $ERROR_RATE" | tee -a {output.training_log}

        # Step 4: If error rate >= 0.5, set stopCodonExcludedFromCDS to false and re-run
        if awk "BEGIN {{exit !($ERROR_RATE >= 0.5)}}"; then
            echo "[INFO] Error rate >= 0.5, setting stopCodonExcludedFromCDS to false" | tee -a {output.training_log}
            sed -i 's/^\(stopCodonExcludedFromCDS\s\+\)\S\+/\1false/' {params.species_params}

            echo "[INFO] Running etraining (second pass, stopCodonExcludedFromCDS=false)..." | tee -a {output.training_log}
            etraining --species={params.species_name} --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
                {input.gb_train} \
                1> {output.training_stdout} \
                2> {output.training_stderr} || true
        else
            echo "[INFO] Error rate < 0.5, keeping stopCodonExcludedFromCDS=true" | tee -a {output.training_log}
        fi

        # Step 5: Extract stop codon frequencies from etraining output
        echo "[INFO] Extracting stop codon frequencies from training output..." | tee -a {output.training_log}

        FREQ_TAG=$(grep "tag:" {output.training_stdout} | tail -1 | sed -E 's/.*\(([0-9.]+)\).*/\1/' || echo "0.33")
        FREQ_TAA=$(grep "taa:" {output.training_stdout} | tail -1 | sed -E 's/.*\(([0-9.]+)\).*/\1/' || echo "0.33")
        FREQ_TGA=$(grep "tga:" {output.training_stdout} | tail -1 | sed -E 's/.*\(([0-9.]+)\).*/\1/' || echo "0.34")

        echo "[INFO] Stop codon frequencies:" | tee -a {output.training_log}
        echo "[INFO]   TAG (amber): $FREQ_TAG" | tee -a {output.training_log}
        echo "[INFO]   TAA (ochre): $FREQ_TAA" | tee -a {output.training_log}
        echo "[INFO]   TGA (opal):  $FREQ_TGA" | tee -a {output.training_log}

        # Step 6: Update species parameters with observed frequencies
        echo "[INFO] Setting stop codon frequencies in species parameters..." | tee -a {output.training_log}

        # For alternative genetic codes, override stop codon probabilities:
        # Code 6/29 (ciliate): TGA codes for Trp, not a stop codon → opalprob=0
        if [ "{params.translation_table}" != "1" ]; then
            echo "[INFO] Alternative translation table {params.translation_table} detected" | tee -a {output.training_log}
            if [ "{params.translation_table}" = "6" ] || [ "{params.translation_table}" = "29" ]; then
                echo "[INFO] Code 6/29: TGA encodes Trp, setting opalprob to 0" | tee -a {output.training_log}
                FREQ_TGA="0"
                # Renormalize TAG and TAA to sum to 1
                FREQ_TAG=$(awk "BEGIN {{printf \"%.4f\", $FREQ_TAG / ($FREQ_TAG + $FREQ_TAA)}}")
                FREQ_TAA=$(awk "BEGIN {{printf \"%.4f\", $FREQ_TAA / ($FREQ_TAG + $FREQ_TAA)}}")
            fi
        fi

        sed -i "s#^\(/Constant/amberprob\s\+\)\S\+#\1$FREQ_TAG#" {params.species_params}
        sed -i "s#^\(/Constant/ochreprob\s\+\)\S\+#\1$FREQ_TAA#" {params.species_params}
        sed -i "s#^\(/Constant/opalprob\s\+\)\S\+#\1$FREQ_TGA#" {params.species_params}

        # Set translation_table in AUGUSTUS parameters (needed for non-standard genetic codes)
        if [ "{params.translation_table}" != "1" ]; then
            if grep -q "^translation_table" {params.species_params}; then
                sed -i "s#^translation_table\s\+\S\+#translation_table {params.translation_table}#" {params.species_params}
            else
                echo "translation_table {params.translation_table}" >> {params.species_params}
            fi
            echo "[INFO] Set translation_table to {params.translation_table} in AUGUSTUS parameters" | tee -a {output.training_log}
        fi

        echo "[INFO] Successfully updated stop codon frequencies in parameters" | tee -a {output.training_log}

        # Verify the changes
        echo "[INFO] Verifying parameter changes:" | tee -a {output.training_log}
        grep -E "^(stopCodonExcludedFromCDS|/Constant/.*prob|translation_table)" {params.species_params} | tee -a {output.training_log}

        echo "[INFO] =======================================" | tee -a {output.training_log}
        echo "[INFO] AUGUSTUS training completed successfully" | tee -a {output.training_log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        STOP_CODON_STATUS=$(grep "stopCodonExcludedFromCDS" {params.species_params} | awk '{{print $2}}')
        cite augustus "$REPORT_DIR"
        """
