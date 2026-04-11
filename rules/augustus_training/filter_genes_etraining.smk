rule filter_genes_etraining:
    """
    Filter out gene structures that cause etraining errors, then cap the
    training set to at most 8000 genes.

    Mirrors braker.pl's training_augustus flow (lines ~6320–6388):
      1. Run etraining on the hint-filtered GenBank file (bonafide.gb,
         which is now braker.pl's trainGb2 equivalent after the convert_to_genbank
         fix). etraining is expected to fail on structurally inconsistent genes.
      2. Parse stderr for bad gene IDs ("n sequence XXX:" pattern).
      3. Remove bad genes via filterGenes.pl.
      4. If the surviving set has > 8000 genes, randomly sub-sample to 8000
         using randomSplit.pl. braker.pl introduces this cap to keep
         optimize_augustus.pl tractable (memory) and to keep the subsequent
         DIAMOND self-comparison from blowing up to O(n^2).

    This rule now runs BEFORE redundancy_removal (the DIAMOND step) to match
    braker.pl's order: etraining filter → 8000-cap → DIAMOND. The previous
    order (DIAMOND → etraining) caused DIAMOND to keep representative genes
    that subsequently failed etraining, introducing non-determinism in which
    loci survived.

    Input:
        gb: hint-filtered GenBank file (post convert_to_genbank fix; this is
            braker.pl's trainGb1+trainGb2 equivalent, NOT the post-DIAMOND set)
        species_marker: Ensures AUGUSTUS species has been created

    Output:
        gb_clean: GenBank file with problematic genes removed and capped at
                  ≤ 8000 entries
        bad_genes_lst: List of gene IDs that failed etraining
        gene_count: Log showing number of genes before/after filtering
    """
    input:
        gb = "output/{sample}/bonafide.gb",
        species_marker = "augustus_config/species/{sample}_rescuedB3/.species_created"
    output:
        gb_clean = "output/{sample}/bonafide.f.clean.gb",
        bad_genes_lst = "output/{sample}/etrain.bad.lst",
        etraining_stdout = "output/{sample}/etraining.stdout",
        etraining_stderr = "output/{sample}/etraining.stderr",
        gene_count = "output/{sample}/etraining_gene_count.txt"
    benchmark:
        "benchmarks/{sample}/filter_genes_etraining/filter_genes_etraining.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        species_name = lambda w: get_species_name(w),
        aug_config = augustus_config_path
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        # Set AUGUSTUS_CONFIG_PATH to our local copy
        export AUGUSTUS_CONFIG_PATH={params.aug_config}

        echo "[INFO] Starting etraining filtering process"
        echo "[INFO] AUGUSTUS_CONFIG_PATH: $AUGUSTUS_CONFIG_PATH"
        echo "[INFO] Species: {params.species_name}"

        # Count genes before filtering
        GENES_BEFORE=$(grep -c "^LOCUS" {input.gb} || echo 0)
        echo "[INFO] ===== ETRAINING FILTERING STATISTICS ====="
        echo "[INFO] Input gene structures: $GENES_BEFORE"

        # Step 1: Run etraining to identify problematic gene structures
        # This will fail on genes with structural inconsistencies, which is expected
        echo "[INFO] Running etraining to identify problematic genes..."
        etraining --species={params.species_name} --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
            {input.gb} \
            1> {output.etraining_stdout} \
            2> {output.etraining_stderr} || true

        # Step 2: Parse stderr to extract bad gene IDs
        # Look for lines like "in sequence GENE_ID:" which indicate errors
        echo "[INFO] Extracting bad gene IDs from etraining errors..."
        grep "n sequence " {output.etraining_stderr} | \
            sed -E 's/.*n sequence ([^:]+):.*/\1/' | \
            sort -u > {output.bad_genes_lst} || touch {output.bad_genes_lst}

        BAD_GENES=$(wc -l < {output.bad_genes_lst} || echo 0)
        echo "[INFO] Problematic genes found: $BAD_GENES"

        # Step 3: Filter out bad genes using filterGenes.pl
        if [ -s {output.bad_genes_lst} ]; then
            echo "[INFO] Filtering out problematic genes with filterGenes.pl..."
            filterGenes.pl {output.bad_genes_lst} {input.gb} > {output.gb_clean}
        else
            echo "[INFO] No problematic genes found, copying input to output"
            cp {input.gb} {output.gb_clean}
        fi

        # Step 4: Cap the training set to ≤ 8000 genes via randomSplit.pl
        # (matches braker.pl lines ~6358–6388). Without this cap,
        # optimize_augustus.pl can OOM on large genomes and the subsequent
        # DIAMOND self-comparison scales O(n^2).
        N_AFTER_ETRAIN=$(grep -c "^LOCUS" {output.gb_clean} || echo 0)
        if [ "$N_AFTER_ETRAIN" -gt 8000 ]; then
            echo "[INFO] Training set has $N_AFTER_ETRAIN genes (> 8000); capping to 8000 via randomSplit.pl"
            # randomSplit.pl writes <input>.test (the 8000 sub-sample) and
            # <input>.train (the remainder). We keep .test, drop .train, and
            # move .test back to the canonical output path.
            randomSplit.pl {output.gb_clean} 8000
            mv {output.gb_clean}.test {output.gb_clean}
            rm -f {output.gb_clean}.train
        else
            echo "[INFO] Training set has $N_AFTER_ETRAIN genes (≤ 8000); no cap applied"
        fi

        # Step 5: Count genes after filtering AND cap
        GENES_AFTER=$(grep -c "^LOCUS" {output.gb_clean} || echo 0)
        GENES_RETAINED=$GENES_AFTER

        echo "[INFO] Genes removed by etraining: $BAD_GENES"
        echo "[INFO] Genes retained for training: $GENES_RETAINED"
        echo "[INFO] Retention rate: $(awk "BEGIN {{printf \"%.1f\", ($GENES_RETAINED/$GENES_BEFORE)*100}}")%"
        echo "[INFO] =========================================="

        # Write gene count statistics
        echo "$GENES_BEFORE -> $GENES_AFTER (removed: $BAD_GENES, retained: $GENES_RETAINED)" > {output.gene_count}

        if [ $GENES_AFTER -eq 0 ]; then
            echo "[ERROR] All genes were filtered out! Training cannot proceed."
            exit 1
        fi

        echo "[INFO] etraining filtering completed successfully"
        echo "[INFO] Clean training set written to: {output.gb_clean}"
        """
