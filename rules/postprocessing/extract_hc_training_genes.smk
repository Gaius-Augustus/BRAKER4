rule extract_hc_training_genes:
    """
    Extract high-confidence training genes and filter out those with in-frame stop codons.

    This rule:
    1. Uses model/hc.gff from failed BRAKER3 run (contains actual HC gene structures)
    2. Translates CDS sequences to check for in-frame stop codons
    3. Filters out genes with stop codons
    4. Creates clean HC GTF for TSEBRA

    BRAKER uses these high-confidence genes as "forced" genes in TSEBRA
    with the --keep_gtf flag, ensuring they appear in final predictions.

    Note: We use model/hc.gff instead of extracting from genome_gmst_for_HC.gtf because:
    - hc.gff contains the actual gene structures for HC genes
    - hc_regions.gtf only contains region coordinates, not gene structures
    - Gene IDs in hc_regions.gtf may not all exist in genome_gmst_for_HC.gtf

    Resources:
        - Single thread (sequential processing)
        - Minimal memory
        - Submitted to SLURM

    Input:
        hc_gtf: High-confidence gene structures from GeneMark-ETP model
        genome: Genome FASTA file for translation

    Output:
        hc_gene_ids: List of HC gene IDs from hc.gff
        hc_genes_raw: HC genes from hc.gff (before stop codon filtering)
        hc_genes_clean: HC genes without in-frame stop codons
        hc_filter_log: Log of filtering process
    """
    input:
        hc_gtf = lambda w: f"output/{w.sample}/dual_etp_merged/hc.gff" if get_braker_mode(w.sample) == 'dual' else f"output/{w.sample}/GeneMark-ETP/hc.gff",
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        hc_gene_ids = "output/{sample}/hc_gene_ids.txt",
        hc_genes_raw = "output/{sample}/hc_genes.raw.gtf",
        hc_genes_clean = "output/{sample}/hc_genes.clean.gtf",
        hc_filter_log = "output/{sample}/hc_genes_filter.log"
    benchmark:
        "benchmarks/{sample}/extract_hc_training_genes/extract_hc_training_genes.txt"
    params:
        output_dir = lambda w: get_output_dir(w),
        translation_table = config.get("translation_table", 1)
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail

        echo "[INFO] ===== EXTRACTING HIGH-CONFIDENCE TRAINING GENES =====" | tee {output.hc_filter_log}

        # Step 1: Copy HC gene structures from model/hc.gff
        echo "[INFO] Copying high-confidence genes from model/hc.gff..." | tee -a {output.hc_filter_log}
        cp {input.hc_gtf} {output.hc_genes_raw}

        # Extract gene IDs for reference
        grep -o 'gene_id "[^"]*"' {input.hc_gtf} | sed 's/gene_id "//g' | sed 's/"//g' | sort -u > {output.hc_gene_ids}

        NUM_HC_IDS=$(wc -l < {output.hc_gene_ids})
        NUM_RAW_LINES=$(grep -c "^" {output.hc_genes_raw} || echo 0)
        echo "[INFO] Found $NUM_HC_IDS high-confidence genes with $NUM_RAW_LINES GTF lines" | tee -a {output.hc_filter_log}

        if [ "$NUM_HC_IDS" -eq 0 ]; then
            echo "[ERROR] No genes found in hc.gff!" | tee -a {output.hc_filter_log}
            exit 1
        fi

        # Step 3: Extract sequences and check for stop codons
        echo "[INFO] Extracting sequences and checking for stop codons..." | tee -a {output.hc_filter_log}

        # Define temporary file paths
        HC_CODINGSEQ="{params.output_dir}/hc_genes.codingseq"
        HC_AA="{params.output_dir}/hc_genes.aa"
        STEM=$(echo $HC_CODINGSEQ | sed 's/\.codingseq$//')

        # Extract coding sequences using modified local script with -d flag
        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {output.hc_genes_raw} \
            -o {params.output_dir}/hc_genes \
            -t {params.translation_table} \
            -d {params.output_dir} \
            1>> {output.hc_filter_log} 2>&1

        # Check if bad_genes.lst was created (genes with in-frame stop codons)
        if [ ! -f {params.output_dir}/bad_genes.lst ]; then
            echo "[INFO] No genes with in-frame stop codons - all HC genes are clean!" | tee -a {output.hc_filter_log}
            cp {output.hc_genes_raw} {output.hc_genes_clean}
        else
            BAD_GENES=$(wc -l < {params.output_dir}/bad_genes.lst)
            echo "[INFO] Found $BAD_GENES genes with in-frame stop codons" | tee -a {output.hc_filter_log}

            # Extract bad gene IDs
            BAD_GENE_IDS=$(mktemp)
            grep -o 'gene_id "[^"]*"' {params.output_dir}/bad_genes.lst | sed 's/gene_id "//g' | sed 's/"//g' | sort -u > $BAD_GENE_IDS || true

            # Filter out bad genes from HC genes
            echo "[INFO] Filtering out genes with stop codons..." | tee -a {output.hc_filter_log}
            > {output.hc_genes_clean}  # Clear file

            while IFS= read -r line; do
                GENE_ID=$(echo "$line" | grep -o 'gene_id "[^"]*"' | sed 's/gene_id "//g' | sed 's/"//g' | head -1)
                if ! grep -q "^$GENE_ID$" $BAD_GENE_IDS; then
                    echo "$line" >> {output.hc_genes_clean}
                fi
            done < {output.hc_genes_raw}

            rm $BAD_GENE_IDS
            rm {params.output_dir}/bad_genes.lst
        fi

        # Count final clean genes
        NUM_CLEAN_LINES=$(grep -c "^" {output.hc_genes_clean} || echo 0)
        # HC genes are in GeneMark-ST format (no gene features), count unique gene_id values
        NUM_CLEAN_GENES=$(grep -oP 'gene_id "[^"]+"' {output.hc_genes_clean} | sort -u | wc -l || echo 0)
        echo "[INFO] Final clean HC GTF: $NUM_CLEAN_LINES lines, $NUM_CLEAN_GENES genes" | tee -a {output.hc_filter_log}

        # Re-extract sequences from clean genes using modified local script
        echo "[INFO] Extracting final sequences from clean genes..." | tee -a {output.hc_filter_log}

        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {output.hc_genes_clean} \
            -o {params.output_dir}/hc_genes \
            -t {params.translation_table} \
            -d {params.output_dir} \
            1>> {output.hc_filter_log} 2>&1

        CODING_SEQ=$(grep -c "^>" $HC_CODINGSEQ || echo 0)
        AA_SEQ=$(grep -c "^>" $HC_AA || echo 0)
        echo "[INFO] Final sequences: $CODING_SEQ CDS, $AA_SEQ proteins" | tee -a {output.hc_filter_log}

        # Clean up temporary sequence files
        echo "[INFO] Cleaning up temporary sequence files..." | tee -a {output.hc_filter_log}
        rm -f "$HC_CODINGSEQ" "$HC_AA"
        echo "[INFO] Removed temporary files: $HC_CODINGSEQ, $HC_AA" | tee -a {output.hc_filter_log}

        echo "[INFO] =======================================" | tee -a {output.hc_filter_log}
        """
