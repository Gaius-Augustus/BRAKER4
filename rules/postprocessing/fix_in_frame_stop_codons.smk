rule get_anno_fasta:
    """
    Extract coding sequences and protein sequences from AUGUSTUS GTF.

    This rule uses getAnnoFastaFromJoingenes.py to extract:
    - Coding sequences (CDS) in FASTA format
    - Protein sequences (amino acids) in FASTA format
    - List of genes with in-frame stop codons (bad_genes.lst)

    The bad_genes.lst file is critical for the next step (fix_in_frame_stop_codons),
    which will re-predict these problematic genes using AUGUSTUS MEA mode.

    Resources:
        - Single thread (I/O bound operation)
        - Minimal memory
        - Submitted to SLURM

    Input:
        gtf: AUGUSTUS predictions in GTF format
        genome: Genome FASTA file

    Output:
        codingseq: Coding sequences in FASTA format
        aa: Protein sequences in FASTA format
        bad_genes: List of genes with in-frame stop codons
    """
    input:
        gtf = get_augustus_gtf,
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        codingseq = "output/{sample}/augustus.hints.codingseq",
        aa = "output/{sample}/augustus.hints.aa",
        bad_genes = "output/{sample}/bad_genes.lst"
    benchmark:
        "benchmarks/{sample}/get_anno_fasta/get_anno_fasta.txt"
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

        echo "[INFO] ===== EXTRACTING CODING SEQUENCES ====="
        echo "[INFO] Input GTF: {input.gtf}"

        # Use modified local script with -d flag to specify output directory for bad_genes.lst
        STEM=$(echo {output.codingseq} | sed 's/\.codingseq$//')

        python3 {script_dir}/getAnnoFastaFromJoingenes.py \
            -g {input.genome} \
            -f {input.gtf} \
            -o $STEM \
            -t {params.translation_table} \
            -d {params.output_dir} \
            1> {params.output_dir}/getAnnoFasta.stdout \
            2> {params.output_dir}/getAnnoFasta.stderr

        # Check if bad_genes.lst was created (genes with in-frame stop codons)
        if [ ! -f {output.bad_genes} ]; then
            # Create empty file if no bad genes found
            touch {output.bad_genes}
            echo "[INFO] No genes with in-frame stop codons found"
        else
            BAD_GENES=$(wc -l < {output.bad_genes})
            echo "[INFO] Found $BAD_GENES genes with in-frame stop codons"
        fi

        # Count sequences
        CODING_SEQ=$(grep -c "^>" {output.codingseq} || echo 0)
        AA_SEQ=$(grep -c "^>" {output.aa} || echo 0)

        echo "[INFO] Extracted $CODING_SEQ coding sequences"
        echo "[INFO] Extracted $AA_SEQ protein sequences"
        echo "[INFO] ======================================="
        """


rule fix_in_frame_stop_codons:
    """
    Fix AUGUSTUS genes with in-frame stop codons.

    AUGUSTUS sometimes predicts genes with in-frame stop codons when:
    - The stop codon is split across an intron (spliced stop codon)
    - There are sequencing errors in the genome
    - The gene has a non-standard genetic code

    This rule re-predicts problematic genes using AUGUSTUS in MEA (Maximum
    Expected Accuracy) mode instead of Viterbi mode. MEA mode is more robust
    to such cases.

    The process:
    1. Read list of genes with in-frame stop codons (from getAnnoFastaFromJoingenes.py)
    2. Extract genomic regions for these genes
    3. Re-predict genes with AUGUSTUS MEA mode
    4. Replace bad predictions with fixed predictions
    5. Re-extract coding sequences and proteins

    Resources:
        - Uses all available CPUs (parallel gene fixing)
        - Full node memory allocation
        - Submitted to SLURM cluster

    Input:
        gtf: AUGUSTUS predictions with some bad genes
        bad_genes: List of genes with in-frame stop codons
        genome: Genome FASTA file
        hintsfile: Merged hints file
        extrinsic_cfg: Extrinsic configuration file
        species_params: AUGUSTUS species parameters

    Output:
        gtf_fixed: AUGUSTUS predictions with fixed genes
        codingseq_fixed: Coding sequences (fixed)
        aa_fixed: Protein sequences (fixed)
        fix_log: Log file from fixing process
    """
    input:
        gtf = get_augustus_gtf,
        bad_genes = "output/{sample}/bad_genes.lst",
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa"),
        hintsfile = get_augustus_hintsfile,
        species_params = "augustus_config/species/{sample}_rescuedB3/{sample}_rescuedB3_parameters.cfg.backup"
    output:
        gtf_fixed = "output/{sample}/augustus.hints.fixed.gtf",
        codingseq_fixed = "output/{sample}/augustus.hints.fixed.codingseq",
        aa_fixed = "output/{sample}/augustus.hints.fixed.aa",
        fix_log = "output/{sample}/fix_in_frame_stop_codons.log"
    benchmark:
        "benchmarks/{sample}/fix_in_frame_stop_codons/fix_in_frame_stop_codons.txt"
    params:
        species_name = lambda w: get_species_name(w),
        output_dir = lambda w: get_output_dir(w),
        translation_table = config.get("translation_table", 1),
        aug_config = augustus_config_path,
        # Mode-aware extrinsic cfg: matches whichever cfg the original AUGUSTUS
        # prediction used for this mode (rnaseq.cfg / ep.cfg / etp.cfg). braker.pl
        # passes the same cfg to fix_ifs_genes that the upstream AUGUSTUS run
        # used. Previously hardcoded to etp.cfg, which was wrong for ET and EP.
        extrinsic_cfg = lambda w: get_extrinsic_cfg(w),
        fix_stem = lambda w: f"output/{w.sample}/augustus.hints_fix_ifs_",
        allow_hinted_splicesites = config.get('allow_hinted_splicesites', 'gcag,atac')
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

        export AUGUSTUS_CONFIG_PATH={params.aug_config}

        echo "[INFO] ===== FIXING IN-FRAME STOP CODONS ====="

        # Create local symlink to genome to allow index file creation
        LOCAL_GENOME={params.output_dir}/genome.fa
        if [ ! -e "$LOCAL_GENOME" ]; then
            ln -s {input.genome} $LOCAL_GENOME
        fi

        # Check if there are any bad genes to fix
        if [ ! -s {input.bad_genes} ]; then
            echo "[INFO] No genes with in-frame stop codons to fix"
            echo "[INFO] Copying original files..."
            cp {input.gtf} {output.gtf_fixed}

            # Extract sequences from original GTF using modified local script
            STEM=$(echo {output.codingseq_fixed} | sed 's/\.codingseq$//')
            python3 {script_dir}/getAnnoFastaFromJoingenes.py \
                -g $LOCAL_GENOME \
                -f {output.gtf_fixed} \
                -o $STEM \
                -t {params.translation_table} \
                -d {params.output_dir}

            echo "[INFO] No fixing needed" > {output.fix_log}
        else
            BAD_GENES=$(wc -l < {input.bad_genes})
            echo "[INFO] Fixing $BAD_GENES genes with in-frame stop codons..."

            # Run fix_in_frame_stop_codon_genes.py (local modified version)
            # Note: Using local script to control temp file locations
            # -m on: BRAKER4 always soft-masks, so the re-prediction must
            # honour repeat regions the same way the original AUGUSTUS run did.
            # braker.pl picks -m on/off conditionally on $soft_mask; in BRAKER4
            # the genome is always soft-masked, so the answer is always "on".
            python3 {script_dir}/fix_in_frame_stop_codon_genes.py \
                -g $LOCAL_GENOME \
                -t {input.gtf} \
                -b {input.bad_genes} \
                -o {params.fix_stem} \
                -s {params.species_name} \
                -m on \
                -u off \
                -U off \
                -a {params.aug_config} \
                -C /opt/cdbfasta \
                -A /opt/Augustus/bin \
                -S /opt/Augustus/scripts \
                -H {input.hintsfile} \
                -e {params.extrinsic_cfg} \
                -w {params.output_dir} \
                --additional_aug_args "--allow_hinted_splicesites={params.allow_hinted_splicesites} " \
                2>&1 | tee {output.fix_log}

            # Move fixed GTF to output location
            mv {params.fix_stem}.gtf {output.gtf_fixed}

            # Extract sequences from fixed GTF using modified local script
            echo "[INFO] Re-extracting sequences from fixed GTF..."
            STEM=$(echo {output.codingseq_fixed} | sed 's/\.codingseq$//')
            python3 {script_dir}/getAnnoFastaFromJoingenes.py \
                -g $LOCAL_GENOME \
                -f {output.gtf_fixed} \
                -o $STEM \
                -t {params.translation_table} \
                -d {params.output_dir}

            FIXED_GENES=$(wc -l < {output.gtf_fixed})
            echo "[INFO] Fixed GTF has $FIXED_GENES lines"
        fi

        # Count final sequences
        CODING_SEQ=$(grep -c "^>" {output.codingseq_fixed} || echo 0)
        AA_SEQ=$(grep -c "^>" {output.aa_fixed} || echo 0)

        echo "[INFO] Final coding sequences: $CODING_SEQ"
        echo "[INFO] Final protein sequences: $AA_SEQ"
        echo "[INFO] ======================================="

        # Cleanup empty log files
        for logfile in {params.output_dir}/getAnnoFasta.stdout {params.output_dir}/getAnnoFasta.stderr; do
            if [ -f "$logfile" ] && [ ! -s "$logfile" ]; then
                rm "$logfile"
            fi
        done
        """
