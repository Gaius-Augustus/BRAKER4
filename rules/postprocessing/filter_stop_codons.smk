rule filter_internal_stop_codons:
    """
    Filter out transcripts with internal in-frame stop codons from final predictions.

    This rule ensures no genes with internal stop codons make it into the final output.
    These problematic genes can come from:
    - AUGUSTUS predictions that weren't fixed by the MEA-mode re-prediction
    - GeneMark-ST (GMST) predictions merged by TSEBRA
    - Assembly errors or genuine frameshift events

    The filtering process:
    1. Extract protein sequences from GTF
    2. Identify transcripts with internal stop codons (not terminal)
    3. Remove all features (gene, transcript, exon, CDS, etc.) for problematic transcripts
    4. Generate clean GTF without any internal stop codons

    Resources:
        - Single thread
        - Minimal memory
        - Submitted to SLURM

    Input:
        braker_merged_gtf: TSEBRA-merged GTF (before stop codon filtering)
        genome: Genome FASTA file

    Output:
        braker_gtf: Final clean GTF without internal stop codons
        filter_log: Log file with filtering statistics
    """
    input:
        braker_merged_gtf = "output/{sample}/braker.rescued.gtf",
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        braker_gtf = "output/{sample}/braker.filtered.gtf",
        filter_log = "output/{sample}/filter_stop_codons.log"
    benchmark:
        "benchmarks/{sample}/filter_internal_stop_codons/filter_internal_stop_codons.txt"
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

        echo "[INFO] ===== FILTERING INTERNAL STOP CODONS =====" | tee {output.filter_log}
        echo "[INFO] Input GTF: {input.braker_merged_gtf}" | tee -a {output.filter_log}
        echo "[INFO] Output GTF: {output.braker_gtf}" | tee -a {output.filter_log}

        # Count genes before filtering
        GENES_BEFORE=$(grep -cP '\tgene\t' {input.braker_merged_gtf} || echo 0)
        TRANSCRIPTS_BEFORE=$(grep -cP '\ttranscript\t' {input.braker_merged_gtf} || echo 0)
        echo "[INFO] Genes before filtering: $GENES_BEFORE" | tee -a {output.filter_log}
        echo "[INFO] Transcripts before filtering: $TRANSCRIPTS_BEFORE" | tee -a {output.filter_log}

        # Extract protein sequences from merged GTF
        # Use modified local script that writes bad_genes.lst to specified output directory
        # This prevents parallel jobs from conflicting on bad_genes.lst in the root directory
        echo "[INFO] Extracting protein sequences to check for internal stop codons..." | tee -a {output.filter_log}

        # Use modified local script with -d flag to specify output directory for bad_genes.lst
        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {input.braker_merged_gtf} \
            -o {params.output_dir}/braker.tsebra \
            -t {params.translation_table} \
            -d {params.output_dir} \
            1> {params.output_dir}/getAnnoFasta_check.stdout \
            2> {params.output_dir}/getAnnoFasta_check.stderr || true

        STEM="{params.output_dir}/braker.tsebra"

        # Check if protein file was created
        if [ ! -f "$STEM.aa" ]; then
            echo "[ERROR] Failed to extract protein sequences!" | tee -a {output.filter_log}
            exit 1
        fi

        # Find transcript IDs with internal stop codons
        echo "[INFO] Identifying transcripts with internal stop codons..." | tee -a {output.filter_log}

        awk '
        BEGIN {{
            bad_count = 0
        }}
        /^>/ {{
            # Extract transcript ID from header (e.g., >g1.t1 -> g1.t1)
            if (seq != "" && seq ~ /\*.*[A-Z]/) {{
                bad_count++
                print tx_id
            }}
            # Parse transcript ID from FASTA header
            tx_id = substr($1, 2)  # Remove leading >
            seq = ""
        }}
        !/^>/ {{seq = seq $0}}
        END {{
            if (seq != "" && seq ~ /\*.*[A-Z]/) {{
                bad_count++
                print tx_id
            }}
            print "TOTAL_BAD:" bad_count > "/dev/stderr"
        }}
        ' "$STEM.aa" > {params.output_dir}/bad_transcript_ids.txt 2> {params.output_dir}/bad_tx_count.txt

        # Extract count and ensure it's a valid number (default to 0 if parsing fails)
        BAD_TX_COUNT=$(grep "TOTAL_BAD:" {params.output_dir}/bad_tx_count.txt | cut -d: -f2 || echo "0")
        # Remove any whitespace and ensure it's numeric
        BAD_TX_COUNT=${{BAD_TX_COUNT:-0}}
        BAD_TX_COUNT=${{BAD_TX_COUNT// /}}
        # If empty or non-numeric, set to 0
        if ! [[ "$BAD_TX_COUNT" =~ ^[0-9]+$ ]]; then
            BAD_TX_COUNT=0
        fi
        echo "[INFO] Found $BAD_TX_COUNT transcripts with internal stop codons" | tee -a {output.filter_log}

        if [ "$BAD_TX_COUNT" -gt 0 ]; then
            echo "[INFO] Filtering out problematic transcripts..." | tee -a {output.filter_log}

            # Create Python script to filter GTF by transcript ID
            cat > {params.output_dir}/filter_gtf.py <<'PYEOF'
import sys

# Load bad transcript IDs
bad_tx_ids = set()
with open(sys.argv[1], 'r') as f:
    for line in f:
        bad_tx_ids.add(line.strip())

print(f"[INFO] Loaded {{len(bad_tx_ids)}} bad transcript IDs to filter", file=sys.stderr)

# Read GTF and track which genes to keep
# A gene is kept only if ALL its transcripts are clean
gene_to_transcripts = {{}}  # gene_id -> set of transcript_ids
bad_genes = set()

# First pass: identify which genes have bad transcripts
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        # Parse attributes
        attrs = {{}}
        for attr in fields[8].split(';'):
            attr = attr.strip()
            if not attr:
                continue
            if ' ' in attr:
                key, val = attr.split(' ', 1)
                attrs[key] = val.strip('"')

        gene_id = attrs.get('gene_id')
        tx_id = attrs.get('transcript_id')

        if gene_id and tx_id:
            if gene_id not in gene_to_transcripts:
                gene_to_transcripts[gene_id] = set()
            gene_to_transcripts[gene_id].add(tx_id)

            # If this transcript is bad, mark the whole gene as bad
            if tx_id in bad_tx_ids:
                bad_genes.add(gene_id)

print(f"[INFO] Found {{len(bad_genes)}} genes with bad transcripts", file=sys.stderr)

# Second pass: write clean GTF (skip bad genes entirely)
filtered_count = 0
kept_count = 0

with open(sys.argv[2], 'r') as f_in, open(sys.argv[3], 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            f_out.write(line)
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            f_out.write(line)
            continue

        # Parse attributes
        attrs = {{}}
        for attr in fields[8].split(';'):
            attr = attr.strip()
            if not attr:
                continue
            if ' ' in attr:
                key, val = attr.split(' ', 1)
                attrs[key] = val.strip('"')

        gene_id = attrs.get('gene_id')

        # Skip entire gene if it has any bad transcripts
        if gene_id in bad_genes:
            filtered_count += 1
            continue

        f_out.write(line)
        kept_count += 1

print(f"[INFO] Filtered {{filtered_count}} lines", file=sys.stderr)
print(f"[INFO] Kept {{kept_count}} lines", file=sys.stderr)
PYEOF

            # Run the filtering script
            python3 {params.output_dir}/filter_gtf.py \
                {params.output_dir}/bad_transcript_ids.txt \
                {input.braker_merged_gtf} \
                {output.braker_gtf} \
                2>&1 | tee -a {output.filter_log}
        else
            echo "[INFO] No problematic transcripts found - copying GTF as-is" | tee -a {output.filter_log}
            cp {input.braker_merged_gtf} {output.braker_gtf}
        fi

        # Count genes after filtering
        GENES_AFTER=$(grep -cP '\tgene\t' {output.braker_gtf} || echo 0)
        TRANSCRIPTS_AFTER=$(grep -cP '\ttranscript\t' {output.braker_gtf} || echo 0)
        GENES_REMOVED=$((GENES_BEFORE - GENES_AFTER))
        TRANSCRIPTS_REMOVED=$((TRANSCRIPTS_BEFORE - TRANSCRIPTS_AFTER))

        echo "[INFO] =======================================" | tee -a {output.filter_log}
        echo "[INFO] FILTERING SUMMARY:" | tee -a {output.filter_log}
        echo "[INFO]   Genes before: $GENES_BEFORE" | tee -a {output.filter_log}
        echo "[INFO]   Genes after: $GENES_AFTER" | tee -a {output.filter_log}
        echo "[INFO]   Genes removed: $GENES_REMOVED" | tee -a {output.filter_log}
        echo "[INFO]   Transcripts before: $TRANSCRIPTS_BEFORE" | tee -a {output.filter_log}
        echo "[INFO]   Transcripts after: $TRANSCRIPTS_AFTER" | tee -a {output.filter_log}
        echo "[INFO]   Transcripts removed: $TRANSCRIPTS_REMOVED" | tee -a {output.filter_log}
        echo "[INFO] =======================================" | tee -a {output.filter_log}

        # Cleanup temporary files
        rm -f "$STEM.aa" "$STEM.codingseq"
        rm -f {params.output_dir}/getAnnoFasta_check.stdout {params.output_dir}/getAnnoFasta_check.stderr
        rm -f {params.output_dir}/bad_transcript_ids.txt {params.output_dir}/bad_tx_count.txt
        rm -f {params.output_dir}/filter_gtf.py

        # Cleanup bad_genes.lst if getAnnoFastaFromJoingenes.py created it
        # It should now be in the output directory, but check both locations for safety
        # Use test && to safely check files in bash strict mode
        test -f "{params.output_dir}/bad_genes.lst" && {{
            echo "[INFO] Renaming bad_genes.lst in output directory" | tee -a {output.filter_log}
            mv -f {params.output_dir}/bad_genes.lst {params.output_dir}/bad_genes.tsebra.lst
            echo "[INFO] This file lists genes with stop codons that were filtered out" | tee -a {output.filter_log}
        }} || true

        # Also check root directory in case it still appeared there
        # Use test && to safely check files in bash strict mode
        test -f "bad_genes.lst" && {{
            echo "[WARNING] bad_genes.lst found in root directory - moving to output" | tee -a {output.filter_log}
            mv -f bad_genes.lst {params.output_dir}/bad_genes.tsebra.lst
        }} || true

        echo "[INFO] Stop codon filtering completed successfully" | tee -a {output.filter_log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        """
