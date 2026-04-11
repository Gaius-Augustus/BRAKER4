rule redundancy_removal:
    """
    Remove redundant protein sequences from training gene set using DIAMOND.

    This rule implements a multi-step process to filter out redundant genes that
    would negatively impact AUGUSTUS training. It converts genes to protein sequences,
    uses DIAMOND (via aa2nonred.pl) to identify non-redundant sequences based on
    sequence similarity, and filters the GenBank file to retain only unique genes.

    This filtering is critical for AUGUSTUS training because:
    - Redundant genes can bias model parameters
    - Similar sequences don't provide additional information
    - Reduces training time while maintaining quality

    Process:
        1. Extract gene names from GenBank file
        2. Filter GTF for training genes
        3. Convert genes to amino acid sequences (creates temporary prot.aa and prot.codingseq)
        4. Run DIAMOND to identify non-redundant sequences (creates temporary prot.nr.aa)
        5. Map non-redundant sequences back to GenBank loci
        6. Filter GenBank file to keep only non-redundant genes
        7. Clean up temporary protein files to save disk space

    Resources:
        - Uses all available CPUs (for DIAMOND parallelization)
        - Full node memory allocation
        - Submitted to SLURM cluster

    Input:
        gb: GenBank file with all training genes
        gtf: Original GeneMark-ETP gene predictions
        genome: Genome assembly FASTA file

    Output:
        gb_filtered: Filtered GenBank file (non-redundant genes only)
        gtf_filtered: Filtered GTF file (non-redundant genes only)
        locus_count: Log file showing reduction (e.g., "1500 -> 1094")

    Note: Temporary files (protein sequences and gene lists) are automatically
          cleaned up at the end to save disk space.
    """
    input:
        # Consume the etraining-cleaned + 8000-capped set so DIAMOND only sees
        # structurally valid genes — matches braker.pl's order (etraining → cap
        # → DIAMOND). The previous order (DIAMOND → etraining) caused DIAMOND
        # to keep representative genes that subsequently failed etraining,
        # introducing non-determinism in which loci survived.
        gb = "output/{sample}/bonafide.f.clean.gb",
        gtf = get_genemark_output,
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        # Output names are unchanged. The "f" suffix now means "fully filtered
        # (etraining + DIAMOND)" rather than just "DIAMOND-filtered".
        gb_filtered = "output/{sample}/bonafide.f.gb",
        gtf_filtered = "output/{sample}/bonafide.f.gtf",
        locus_count = "output/{sample}/locus_count.txt"
    benchmark:
        "benchmarks/{sample}/redundancy_removal/redundancy_removal.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        translation_table = config.get("translation_table", 1)
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail

        # Create output directory if it doesn't exist
        OUTPUT_DIR=$(dirname {output.gb_filtered})
        mkdir -p $OUTPUT_DIR
        echo "[INFO] Starting redundancy removal process"

        # Define temporary file paths
        TRAINGENES_LST="$OUTPUT_DIR/traingenes.lst"
        NONRED_LST="$OUTPUT_DIR/nonred.lst"
        LOCI_LST="$OUTPUT_DIR/loci.lst"
        NONRED_LOCI_LST="$OUTPUT_DIR/nonred.loci.lst"

        # Step 1: Extract training gene names from GenBank file
        cat {input.gb} | perl -ne 'if(m/\/gene=\"(\S+)\"/){{print "$1\n";}}' | sort -u > $TRAINGENES_LST
        echo "[INFO] Extracted training genes to $TRAINGENES_LST"

        # Step 2: Filter GTF file for training genes
        grep -f $TRAINGENES_LST -F {input.gtf} > {output.gtf_filtered}
        echo "[INFO] Filtered GTF file to {output.gtf_filtered}"

        # Step 3: Convert GTF to amino acid sequences using modified local script
        # getAnnoFastaFromJoingenes.py creates two files: OUTPUT.codingseq and OUTPUT.aa
        # The -o parameter is the stem, so it will create OUTPUT.aa and OUTPUT.codingseq
        PROT_STEM="$(dirname {output.gtf_filtered})/prot"
        PROT_AA="$PROT_STEM.aa"
        PROT_CODINGSEQ="$PROT_STEM.codingseq"
        PROT_NR_AA="$PROT_STEM.nr.aa"

        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {output.gtf_filtered} \
            -o $PROT_STEM \
            -t {params.translation_table} \
            -d $OUTPUT_DIR
        echo "[INFO] Converted GTF to amino acids: $PROT_AA"

        # Step 4: Remove redundant sequences using DIAMOND
        # aa2nonred.pl requires --diamond flag and DIAMOND_PATH
        # DIAMOND is located at /opt/ETP/tools/diamond in the BRAKER3 container
        aa2nonred.pl $PROT_AA $PROT_NR_AA --DIAMOND_PATH=/opt/ETP/tools --diamond --cores={threads}
        echo "[INFO] Removed redundant sequences: $PROT_NR_AA"

        # Step 5: Extract non-redundant gene list
        grep '>' $PROT_NR_AA | perl -pe 's/>//' > $NONRED_LST
        echo "[INFO] Extracted non-redundant gene list: $NONRED_LST"

        # Step 6: Create loci mapping from GenBank file
        cat {input.gb} | perl -ne '
            if ($_ =~ m/LOCUS\s+(\S+)\s/) {{
                $txLocus = $1;
            }} elsif ($_ =~ m/\/gene=\"(\S+)\"/) {{
                $txInGb3{{$1}} = $txLocus
            }}
            if(eof()) {{
                foreach (keys %txInGb3) {{
                    print "$_\t$txInGb3{{$_}}\n";
                }}
            }}' > $LOCI_LST
        echo "[INFO] Created loci mapping: $LOCI_LST"

        # Step 7: Get non-redundant loci
        grep -f $NONRED_LST $LOCI_LST | cut -f2 > $NONRED_LOCI_LST
        echo "[INFO] Extracted non-redundant loci: $NONRED_LOCI_LST"

        # Step 8: Filter GenBank file to keep only non-redundant genes
        filterGenesIn.pl $NONRED_LOCI_LST {input.gb} > {output.gb_filtered}
        echo "[INFO] Filtered GenBank file: {output.gb_filtered}"

        # Step 9: Count loci before and after filtering (for logging)
        BEFORE=$(grep -c LOCUS {input.gb} || echo 0)
        AFTER=$(grep -c LOCUS {output.gb_filtered} || echo 0)
        echo "$BEFORE -> $AFTER" > {output.locus_count}
        echo "[INFO] Locus count: $BEFORE (before) -> $AFTER (after redundancy removal)"

        # Step 10: Clean up temporary files
        echo "[INFO] Cleaning up temporary files..."
        rm -f "$PROT_AA" "$PROT_CODINGSEQ" "$PROT_NR_AA" "$TRAINGENES_LST" "$NONRED_LST" "$LOCI_LST" "$NONRED_LOCI_LST"
        echo "[INFO] Removed temporary protein files: $PROT_AA, $PROT_CODINGSEQ, $PROT_NR_AA"
        echo "[INFO] Removed temporary list files: $TRAINGENES_LST, $NONRED_LST, $LOCI_LST, $NONRED_LOCI_LST"

        echo "[INFO] Redundancy removal completed successfully"

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        DIAMOND_VER=$(diamond version 2>&1 | awk '{{print $NF}}' || true)
        if [ -z "$DIAMOND_VER" ]; then
            DIAMOND_VER=$(/opt/ETP/tools/diamond version 2>&1 | awk '{{print $NF}}' || true)
        fi
        ( flock 9; printf "DIAMOND\t%s\n" "$DIAMOND_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite diamond "$REPORT_DIR"
        """
