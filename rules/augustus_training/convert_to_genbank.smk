rule convert_to_genbank:
    """
    Convert GTF gene predictions to GenBank format for AUGUSTUS training,
    then filter the GenBank entries to keep only the hint-supported subset.

    Mirrors braker.pl's two-step training set construction (training_augustus,
    lines ~6131 and ~6280):

      1. gff2gbSmallDNA.pl on the FULL GeneMark predictions builds the
         "trainGb1" GenBank file (one entry per gene, with flanking regions).
      2. filterGenesIn_mRNAname.pl filters trainGb1 down to only those genes
         whose transcript_ids appear in the hint-supported subset
         (genemark.f.good.gtf in ET/EP/ES; training.gtf in ETP/IsoSeq/dual).
         The result is braker.pl's "trainGb2" — what AUGUSTUS will be trained on.

    The hint-based filter step was missing from earlier BRAKER4 versions, which
    silently trained AUGUSTUS on the full GeneMark output (including ab-initio
    calls with zero hint support).

    In ETP/IsoSeq/dual modes the full and filtered GTFs are the same upstream
    GeneMark-ETP training.gtf, so the filter step is a no-op there.

    Input:
        gtf: full GeneMark predictions (used to build the GenBank entries)
        filtered_gtf: hint-supported subset of GeneMark predictions (used to
                      select which entries to keep — the trainGb2 step)
        flanking_value: Flanking region size (in bp) from compute_flanking_region
        genome: Genome assembly FASTA file

    Output:
        gb: GenBank file containing only hint-supported training genes
    """
    input:
        gtf = get_genemark_output,
        filtered_gtf = get_genemark_training_gtf,
        flanking_value = "output/{sample}/flanking_dna_value.txt",
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa")
    output:
        gb = "output/{sample}/bonafide.gb"
    benchmark:
        "benchmarks/{sample}/convert_to_genbank/convert_to_genbank.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        # Read the flanking DNA value from the file
        FLANKING=$(cat {input.flanking_value})
        echo "[INFO] Using flanking DNA value: $FLANKING"

        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.gb})

        # Step 1 — gff2gbSmallDNA.pl on the FULL GeneMark output (trainGb1)
        TRAINGB1={output.gb}.full
        gff2gbSmallDNA.pl {input.gtf} {input.genome} $FLANKING $TRAINGB1
        N_GB1=$(grep -c '^LOCUS' $TRAINGB1 || echo 0)
        echo "[INFO] Step 1: gff2gbSmallDNA.pl produced $N_GB1 GenBank entries"

        # Step 2 — filter to hint-supported transcripts (trainGb2)
        # When the full and filtered GTFs are identical (ETP/IsoSeq/dual or ES),
        # this is a no-op in terms of selected genes.
        filterGenesIn_mRNAname.pl {input.filtered_gtf} $TRAINGB1 > {output.gb}
        N_GB2=$(grep -c '^LOCUS' {output.gb} || echo 0)
        echo "[INFO] Step 2: filterGenesIn_mRNAname.pl retained $N_GB2 hint-supported entries"

        # Drop the intermediate full set; downstream rules consume {output.gb} only
        rm -f $TRAINGB1

        if [ $N_GB2 -eq 0 ]; then
            echo "[ERROR] No genes survived hint-based filtering. Check that"
            echo "        {input.filtered_gtf} contains transcript_ids matching"
            echo "        the entries in the GenBank file built from {input.gtf}."
            exit 1
        fi

        echo "[INFO] GTF to GenBank conversion completed"
        echo "[INFO] Output written to {output.gb}"
        """
