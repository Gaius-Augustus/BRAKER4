
rule compute_flanking_region:
    """
    Calculate optimal flanking region size for training AUGUSTUS.

    This rule runs computeFlankingRegion.pl from BRAKER to determine the optimal
    amount of flanking DNA to include around genes for training. The script analyzes
    gene structures in the GTF file and calculates an appropriate flanking region
    based on intergenic distances (minimum of 10000bp or the calculated value).

    Input:
        gtf: GeneMark training genes of the effective mode (get_genemark_output)

    Output:
        flanking_file: Full output from computeFlankingRegion.pl (includes statistics)
        flanking_value: Single number representing the flanking DNA length in bp
                       (used by downstream rules for GenBank conversion)
    """
    input:
        gtf = get_genemark_output,
        mode_decisions = mode_decision_files
    output:
        flanking_file = "output/{sample}/genome_gmst_for_HC.flanking.gtf",
        flanking_value = "output/{sample}/flanking_dna_value.txt"
    benchmark:
        "benchmarks/{sample}/compute_flanking_region/compute_flanking_region.txt"
    container:
        BRAKER3_CONTAINER
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail

        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.flanking_file})
        echo "[INFO] Created output directory: $(dirname {output.flanking_file})"

        # Run computeFlankingRegion.pl and capture output to file
        computeFlankingRegion.pl {input.gtf} > {output.flanking_file}

        # Extract the flanking DNA value from the output file
        # Format: "The flanking_DNA value is: 1094 (the Minimum of 10 000 and 1094"
        # We want to extract the first number (1094 in this example)
        grep "The flanking_DNA value is:" {output.flanking_file} | sed -E 's/.*The flanking_DNA value is: ([0-9]+).*/\1/' > {output.flanking_value}

        echo "[INFO] computeFlankingRegion.pl completed for {input.gtf}"
        echo "[INFO] Output written to {output.flanking_file}"
        echo "[INFO] Flanking DNA value extracted: $(cat {output.flanking_value})"
        """
