rule extract_final_sequences:
    """
    Extract coding sequences and proteins from TSEBRA-merged predictions.

    This rule extracts sequences from the final merged GTF (TSEBRA + rename_gtf.py output)
    and creates BRAKER-standard output files:
    - braker.gtf: Already created by TSEBRA rule (merged and renamed)
    - braker.codingseq: Coding sequences (nucleotides)
    - braker.aa: Protein sequences (amino acids)

    Resources:
        - Single thread
        - Submitted to SLURM

    Input:
        braker_gtf: TSEBRA-merged and renamed GTF (from merge_predictions rule)
        genome: Genome FASTA file

    Output:
        braker_codingseq: Final coding sequences
        braker_aa: Final protein sequences
    """
    input:
        braker_gtf = "output/{sample}/braker.gtf",
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        braker_codingseq = "output/{sample}/braker.codingseq",
        braker_aa = "output/{sample}/braker.aa"
    benchmark:
        "benchmarks/{sample}/extract_final_sequences/extract_final_sequences.txt"
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

        echo "[INFO] ===== EXTRACTING FINAL SEQUENCES ====="
        echo "[INFO] Using braker.gtf from TSEBRA (already merged and renamed)"

        # Extract sequences from braker.gtf using modified local script
        STEM=$(echo {output.braker_codingseq} | sed 's/\.codingseq$//')
        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {input.braker_gtf} \
            -o $STEM \
            -t {params.translation_table} \
            -d {params.output_dir} \
            1> {params.output_dir}/getAnnoFasta_final.stdout \
            2> {params.output_dir}/getAnnoFasta_final.stderr

        # Count final statistics (use grep -P for proper tab matching)
        GENES_GTF=$(grep -cP '\tgene\t' {input.braker_gtf} || echo 0)
        TRANSCRIPTS_GTF=$(grep -cP '\ttranscript\t' {input.braker_gtf} || echo 0)
        CDS_COUNT=$(grep -cP '\tCDS\t' {input.braker_gtf} || echo 0)
        CODING_SEQ=$(grep -c "^>" {output.braker_codingseq} || echo 0)
        AA_SEQ=$(grep -c "^>" {output.braker_aa} || echo 0)

        echo "[INFO] Final gene predictions:"
        echo "[INFO]   Genes: $GENES_GTF"
        echo "[INFO]   Transcripts: $TRANSCRIPTS_GTF"
        echo "[INFO]   CDS features: $CDS_COUNT"
        echo "[INFO]   Coding sequences: $CODING_SEQ"
        echo "[INFO]   Protein sequences: $AA_SEQ"

        # Verify no internal in-frame stop codons in protein sequences
        # (This is a sanity check - filtering should have removed all such genes)
        echo "[INFO] Verifying no internal stop codons in protein sequences..."

        awk '
        /^>/ {{
            if (seq != "" && seq ~ /\*.*[A-Z]/) {{
                internal_stops++
            }}
            seq=""
        }}
        !/^>/ {{seq=seq$0}}
        END {{
            if (seq != "" && seq ~ /\*.*[A-Z]/) {{
                internal_stops++
            }}
            print internal_stops + 0
        }}
        ' {output.braker_aa} > {params.output_dir}/internal_stops_verify.txt

        INTERNAL_STOPS=$(cat {params.output_dir}/internal_stops_verify.txt)

        if [ "$INTERNAL_STOPS" -gt 0 ]; then
            echo "[ERROR] VERIFICATION FAILED: Found $INTERNAL_STOPS sequences with internal stop codons!"
            echo "[ERROR] This should not happen - filtering step may have failed!"
            exit 1
        else
            echo "[INFO] ✓ Verification passed: No internal stop codons found"
            rm -f {params.output_dir}/internal_stops_verify.txt
        fi

        echo "[INFO] ======================================="

        # Cleanup empty log files
        for logfile in {params.output_dir}/getAnnoFasta_final.stdout {params.output_dir}/getAnnoFasta_final.stderr; do
            if [ -f "$logfile" ] && [ ! -s "$logfile" ]; then
                rm "$logfile"
            fi
        done

        # Cleanup bad_genes.lst if it was created (shouldn't happen after filtering)
        # Use test && to safely check files in bash strict mode
        test -f "bad_genes.lst" && {{
            echo "[WARNING] bad_genes.lst was created in root directory - moving to output"
            mv bad_genes.lst {params.output_dir}/bad_genes.lst.unexpected
        }} || true
        """


rule assess_completeness:
    """
    Assess gene set completeness using compleasm (BUSCO assessment).

    This rule runs compleasm in protein mode on the final predicted proteins
    to evaluate the completeness and quality of the gene predictions using
    single-copy ortholog benchmarks.

    Resources:
        - Uses all available CPUs for faster BUSCO searches
        - Submitted to SLURM cluster

    Input:
        braker_aa: Final protein sequences

    Output:
        compleasm_summary: Compleasm summary file with BUSCO statistics
        compleasm_log: Log file from compleasm run
    """
    input:
        braker_aa = "output/{sample}/braker.aa"
    output:
        compleasm_summary = "output/{sample}/compleasm_proteins/summary.txt",
        compleasm_log = "output/{sample}/compleasm_proteins.log"
    benchmark:
        "benchmarks/{sample}/assess_completeness/assess_completeness.txt"
    params:
        busco_lineage = lambda w: get_busco_lineage(w),
        compleasm_outdir = lambda w: f"output/{w.sample}/compleasm_proteins",
        library_path = config['compleasm_download_path']
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail

        echo "[INFO] ===== ASSESSING GENE SET COMPLETENESS WITH COMPLEASM =====" | tee {output.compleasm_log}
        echo "[INFO] BUSCO lineage: {params.busco_lineage}" | tee -a {output.compleasm_log}
        echo "[INFO] Threads: {threads}" | tee -a {output.compleasm_log}
        echo "[INFO] Input proteins: {input.braker_aa}" | tee -a {output.compleasm_log}

        # Count proteins
        PROTEIN_COUNT=$(grep -c "^>" {input.braker_aa} || echo 0)
        echo "[INFO] Total proteins: $PROTEIN_COUNT" | tee -a {output.compleasm_log}

        if [ $PROTEIN_COUNT -eq 0 ]; then
            echo "[WARNING] No proteins found - skipping compleasm" | tee -a {output.compleasm_log}
            mkdir -p {params.compleasm_outdir}
            echo "No proteins available for assessment" > {output.compleasm_summary}
            exit 0
        fi

        # Remove output directory if it exists (compleasm requires clean directory)
        if [ -d "{params.compleasm_outdir}" ]; then
            rm -rf {params.compleasm_outdir}
        fi

        mkdir -p {params.compleasm_outdir}
        mkdir -p {params.library_path}

        # compleasm only supports odb12; convert e.g. eukaryota_odb12 -> eukaryota_odb12
        COMPLEASM_LINEAGE=$(echo "{params.busco_lineage}" | sed 's/_odb[0-9]*$/_odb12/')
        echo "[INFO] Running compleasm in protein mode..." | tee -a {output.compleasm_log}
        echo "[INFO] Library path: {params.library_path}" | tee -a {output.compleasm_log}
        echo "[INFO] Lineage (odb12): $COMPLEASM_LINEAGE" | tee -a {output.compleasm_log}

        compleasm.py protein \
            -p {input.braker_aa} \
            -l $COMPLEASM_LINEAGE \
            -t {threads} \
            -o {params.compleasm_outdir} \
            -L {params.library_path} \
            2>&1 | tee -a {output.compleasm_log} || true

        # Check if summary was created (compleasm may put it in a lineage subdirectory)
        if [ ! -f {output.compleasm_summary} ]; then
            FOUND_SUMMARY=$(find {params.compleasm_outdir} -name "summary.txt" 2>/dev/null | head -1)
            if [ -n "$FOUND_SUMMARY" ] && [ -f "$FOUND_SUMMARY" ]; then
                cp "$FOUND_SUMMARY" {output.compleasm_summary}
                echo "[INFO] Copied summary from $FOUND_SUMMARY" | tee -a {output.compleasm_log}
            else
                echo "[WARNING] Compleasm did not produce summary file" | tee -a {output.compleasm_log}
                mkdir -p {params.compleasm_outdir}
                echo "Compleasm assessment failed" > {output.compleasm_summary}
            fi
        else
            echo "[INFO] Compleasm assessment completed" | tee -a {output.compleasm_log}
            echo "[INFO] Summary:" | tee -a {output.compleasm_log}
            cat {output.compleasm_summary} | tee -a {output.compleasm_log}
        fi

        echo "[INFO] =======================================" | tee -a {output.compleasm_log}
        """


rule generate_statistics:
    """
    Generate comprehensive statistics about the final gene predictions.

    This rule creates a summary report with key statistics about:
    - Number of genes, transcripts, and features
    - Gene structure statistics (exons per gene, CDS length, etc.)
    - Comparison with original failed BRAKER3 predictions
    - Training gene statistics

    Resources:
        - Local rule (no cluster submission)
        - Instant execution

    Input:
        braker_gtf: Final BRAKER GTF file
        braker_aa: Final protein sequences
        braker_codingseq: Final coding sequences
        etraining_count: Training gene statistics

    Output:
        statistics: Summary statistics file
    """
    input:
        braker_gtf = "output/{sample}/braker.gtf",
        braker_gff3 = "output/{sample}/braker.gff3",
        braker_aa = "output/{sample}/braker.aa",
        braker_codingseq = "output/{sample}/braker.codingseq",
        etraining_count = "output/{sample}/etraining_gene_count.txt",
        bonafide_gb = "output/{sample}/bonafide.f.clean.gb",
        gb_train = "output/{sample}/train.gb.train",
        gb_train_train = "output/{sample}/train.gb.train.train",
        compleasm_genome_summary = "output/{sample}/compleasm_genome_out/summary.txt",
        compleasm_proteins_summary = "output/{sample}/compleasm_proteins/summary.txt",
        hintsfile = "output/{sample}/hintsfile.gff"
    output:
        statistics = "output/{sample}/braker_stats.txt"
    benchmark:
        "benchmarks/{sample}/generate_statistics/generate_statistics.txt"
    params:
        species_name = lambda w: get_species_name(w),
        braker_dir = lambda w: get_braker_dir(w),
        output_dir = lambda w: get_output_dir(w)
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail

        echo "[INFO] ===== GENERATING FINAL STATISTICS ====="

        # Create statistics file
        cat > {output.statistics} <<'EOF'
================================================================================
           BRAKER3 RESCUE - FINAL GENE PREDICTION STATISTICS
================================================================================

RESCUE WORKFLOW SUMMARY:
------------------------
This workflow rescued a failed BRAKER3 run by:
1. Reusing GeneMark-ETP high-confidence gene predictions
2. Training AUGUSTUS parameters with these genes
3. Running AUGUSTUS with hints (RNA-Seq + protein + compleasm)
4. Fixing genes with in-frame stop codons

INPUT DATA:
-----------
EOF

        echo "Species name: {params.species_name}" >> {output.statistics}
        echo "BRAKER directory: {params.braker_dir}" >> {output.statistics}
        echo "Output directory: {params.output_dir}" >> {output.statistics}
        echo "" >> {output.statistics}

        echo "TRAINING STATISTICS:" >> {output.statistics}
        echo "--------------------" >> {output.statistics}

        # Training genes
        if [ -f {input.etraining_count} ]; then
            echo "Training genes: $(cat {input.etraining_count})" >> {output.statistics}
        fi

        TRAINING_GENES=$(grep -c "^LOCUS" {input.bonafide_gb} || echo 0)
        echo "Final training genes (after filtering): $TRAINING_GENES" >> {output.statistics}

        # Genes used for optimize_augustus.pl (if optimization was run)
        if [ -f {input.gb_train} ] && [ -f {input.gb_train_train} ]; then
            GENES_FINAL_ETRAINING=$(grep -c "^LOCUS" {input.gb_train} || echo 0)
            GENES_OPTIMIZE=$(grep -c "^LOCUS" {input.gb_train_train} || echo 0)
            echo "Genes used for final etraining: $GENES_FINAL_ETRAINING" >> {output.statistics}
            echo "Genes used for optimize_augustus.pl: $GENES_OPTIMIZE" >> {output.statistics}
        fi
        echo "" >> {output.statistics}

        echo "FINAL PREDICTIONS:" >> {output.statistics}
        echo "------------------" >> {output.statistics}

        # Gene counts (use $'\t' for actual tab character in bash)
        GENES=$(grep -cP '\tgene\t' {input.braker_gtf} || echo 0)
        TRANSCRIPTS=$(grep -cP '\ttranscript\t' {input.braker_gtf} || echo 0)
        CDS=$(grep -cP '\tCDS\t' {input.braker_gtf} || echo 0)
        EXONS=$(grep -cP '\texon\t' {input.braker_gtf} || echo 0)
        INTRONS=$(grep -cP '\tintron\t' {input.braker_gtf} || echo 0)
        START_CODONS=$(grep -cP '\tstart_codon\t' {input.braker_gtf} || echo 0)
        STOP_CODONS=$(grep -cP '\tstop_codon\t' {input.braker_gtf} || echo 0)

        echo "Total genes: $GENES" >> {output.statistics}
        echo "Total transcripts: $TRANSCRIPTS" >> {output.statistics}
        echo "CDS features: $CDS" >> {output.statistics}
        echo "Exons: $EXONS" >> {output.statistics}
        echo "Introns: $INTRONS" >> {output.statistics}
        echo "Start codons: $START_CODONS" >> {output.statistics}
        echo "Stop codons: $STOP_CODONS" >> {output.statistics}
        echo "" >> {output.statistics}

        # Calculate average exons per gene
        if [ $GENES -gt 0 ]; then
            AVG_EXONS=$(awk "BEGIN {{printf \"%.2f\", $EXONS/$GENES}}")
            echo "Average exons per gene: $AVG_EXONS" >> {output.statistics}
        fi

        # Calculate mono-exonic vs multi-exonic gene ratio
        # Count genes with exactly 1 exon (mono-exonic) vs genes with >1 exon (multi-exonic)
        # Group exons by gene_id and count
        MONO_EXONIC=$(awk -F'\t' '$3 == "exon" {{
            match($9, /gene_id "([^"]+)"/, arr);
            gene_exon_count[arr[1]]++;
        }}
        END {{
            mono = 0;
            for (gene in gene_exon_count) {{
                if (gene_exon_count[gene] == 1) mono++;
            }}
            print mono;
        }}' {input.braker_gtf})

        MULTI_EXONIC=$(awk -F'\t' '$3 == "exon" {{
            match($9, /gene_id "([^"]+)"/, arr);
            gene_exon_count[arr[1]]++;
        }}
        END {{
            multi = 0;
            for (gene in gene_exon_count) {{
                if (gene_exon_count[gene] > 1) multi++;
            }}
            print multi;
        }}' {input.braker_gtf})

        echo "Mono-exonic genes: $MONO_EXONIC" >> {output.statistics}
        echo "Multi-exonic genes: $MULTI_EXONIC" >> {output.statistics}

        if [ $MULTI_EXONIC -gt 0 ]; then
            MONO_MULTI_RATIO=$(awk "BEGIN {{printf \"%.3f\", $MONO_EXONIC/$MULTI_EXONIC}}")
            echo "Mono:Multi ratio: $MONO_MULTI_RATIO" >> {output.statistics}
        else
            echo "Mono:Multi ratio: N/A (no multi-exonic genes)" >> {output.statistics}
        fi
        echo "" >> {output.statistics}

        echo "SEQUENCE STATISTICS:" >> {output.statistics}
        echo "--------------------" >> {output.statistics}

        # Sequence counts
        CODING_SEQ=$(grep -c "^>" {input.braker_codingseq} || echo 0)
        PROTEINS=$(grep -c "^>" {input.braker_aa} || echo 0)

        echo "Coding sequences: $CODING_SEQ" >> {output.statistics}
        echo "Protein sequences: $PROTEINS" >> {output.statistics}
        echo "" >> {output.statistics}

        # Calculate total CDS length
        TOTAL_CDS_LEN=$(grep -v "^>" {input.braker_codingseq} | tr -d '\\n' | wc -c)
        if [ $CODING_SEQ -gt 0 ]; then
            AVG_CDS_LEN=$(awk "BEGIN {{printf \"%.0f\", $TOTAL_CDS_LEN/$CODING_SEQ}}")
            echo "Total CDS length: $TOTAL_CDS_LEN bp" >> {output.statistics}
            echo "Average CDS length: $AVG_CDS_LEN bp" >> {output.statistics}
        fi

        # Calculate total protein length
        TOTAL_AA_LEN=$(grep -v "^>" {input.braker_aa} | tr -d '\\n' | wc -c)
        if [ $PROTEINS -gt 0 ]; then
            AVG_AA_LEN=$(awk "BEGIN {{printf \"%.0f\", $TOTAL_AA_LEN/$PROTEINS}}")
            echo "Total protein length: $TOTAL_AA_LEN aa" >> {output.statistics}
            echo "Average protein length: $AVG_AA_LEN aa" >> {output.statistics}
        fi
        echo "" >> {output.statistics}

        # Report protein quality (should be clean after filtering step)
        echo "PROTEIN QUALITY:" >> {output.statistics}
        echo "----------------" >> {output.statistics}
        echo "Internal in-frame stop codons: None (filtered out before sequence extraction)" >> {output.statistics}
        echo "All genes with internal stop codons were removed after TSEBRA merging" >> {output.statistics}
        echo "See filter_stop_codons.log for details on how many genes were filtered" >> {output.statistics}
        echo "" >> {output.statistics}

        echo "HINT SUPPORT ANALYSIS:" >> {output.statistics}
        echo "----------------------" >> {output.statistics}
        echo "Analyzing how well genes are supported by hints (evidence)..." >> {output.statistics}
        echo "" >> {output.statistics}

        # Run hint support analysis script
        python3 scripts/analyze_hint_support.py \
            --gtf {input.braker_gtf} \
            --hints {input.hintsfile} \
            --output {params.output_dir}/hint_support_details.txt

        # Append results to main statistics file
        cat {params.output_dir}/hint_support_details.txt >> {output.statistics}
        echo "" >> {output.statistics}

        echo "BUSCO COMPLETENESS ASSESSMENT:" >> {output.statistics}
        echo "-------------------------------" >> {output.statistics}

        # Genome-level compleasm results
        echo "Genome-level assessment (compleasm on genome sequence):" >> {output.statistics}
        cat "{input.compleasm_genome_summary}" >> {output.statistics}
        echo "" >> {output.statistics}

        # Protein-level compleasm results
        echo "Protein-level assessment (compleasm on predicted proteins):" >> {output.statistics}
        cat "{input.compleasm_proteins_summary}" >> {output.statistics}
        echo "" >> {output.statistics}

        echo "OUTPUT FILES:" >> {output.statistics}
        echo "-------------" >> {output.statistics}
        echo "GTF: {input.braker_gtf}" >> {output.statistics}
        echo "GFF3: {input.braker_gff3}" >> {output.statistics}
        echo "Coding sequences: {input.braker_codingseq}" >> {output.statistics}
        echo "Protein sequences: {input.braker_aa}" >> {output.statistics}
        echo "" >> {output.statistics}

        echo "================================================================================" >> {output.statistics}
        echo "" >> {output.statistics}

        # Display statistics
        cat {output.statistics}

        echo "[INFO] Statistics written to: {output.statistics}"
        echo "[INFO] ======================================="
        """
