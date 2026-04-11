rule downsample_training_genes:
    """
    Downsample training genes according to Poisson distribution.

    This rule uses downsample_traingenes.pl to reduce the number of training genes,
    with downsampling based on the number of exons in each gene. Genes with many
    exons are more likely to be retained, while single-exon genes are downsampled
    more aggressively. This creates a balanced training set that doesn't overweight
    single-exon genes, which are often more abundant but less informative for
    training splice site parameters.

    The Poisson distribution with lambda=2 (default in BRAKER) determines the
    probability of keeping a gene based on its exon count. This approach:
    - Prevents training sets dominated by single-exon genes
    - Retains multi-exon genes which are more informative for intron parameters
    - Creates a more balanced representation of gene structures

    Resources:
        - Uses 1 thread (single-threaded Perl script)
        - Minimal memory requirement (just reading/filtering GTF)
        - Submitted to SLURM

    Input:
        gb_diamond: GenBank file after etraining filtering AND DIAMOND
                    redundancy removal. We deliberately consume the
                    post-DIAMOND file (not bonafide.f.clean.gb) so the
                    GENES_BEFORE counter written to downsample_gene_count.txt
                    matches DIAMOND's "after" count. Otherwise the report
                    chain reads "DIAMOND: N→M loci. Downsampled to X of N"
                    with N appearing twice instead of N→M→X.

    Output:
        gb_downsampled: GenBank file with downsampled gene set
        downsample_log: Log showing downsampling statistics
    """
    input:
        gb_diamond = "output/{sample}/bonafide.f.gb",
        gtf_filtered = "output/{sample}/bonafide.f.gtf"
    output:
        gb_downsampled = "output/{sample}/bonafide.f.clean.d.gb",
        gtf_for_downsample = "output/{sample}/bonafide.f.clean.gtf",
        gtf_downsampled = "output/{sample}/bonafide.f.clean.d.gtf",
        downsample_log = "output/{sample}/downsample_traingenes.log",
        gene_count = "output/{sample}/downsample_gene_count.txt"
    benchmark:
        "benchmarks/{sample}/downsample_training_genes/downsample_training_genes.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        output_dir = lambda w: get_output_dir(w)
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== TRAINING GENE DOWNSAMPLING ====="

        # Count genes before downsampling
        GENES_BEFORE=$(grep -c "^LOCUS" {input.gb_diamond} || echo 0)
        echo "[INFO] Input genes: $GENES_BEFORE"

        # Step 1: Extract gene names from the post-DIAMOND GenBank file
        echo "[INFO] Extracting non-redundant gene names..."
        cat {input.gb_diamond} | perl -ne 'if(m/\/gene=\"(\S+)\"/){{print "$1\n";}}' | sort -u > {params.output_dir}/clean_genes.lst

        # Step 2: Filter the GTF to keep only those genes
        echo "[INFO] Filtering GTF for non-redundant genes..."
        grep -f {params.output_dir}/clean_genes.lst -F {input.gtf_filtered} > {output.gtf_for_downsample}

        # Check single-exon gene proportion
        TOTAL_GENES=$(awk -F'\t' '$3=="CDS"' {output.gtf_for_downsample} | \
            grep -oP 'transcript_id "[^"]+"' | sort -u | wc -l)
        SINGLE_EXON=$(awk -F'\t' '$3=="CDS"' {output.gtf_for_downsample} | \
            grep -oP 'transcript_id "[^"]+"' | sort | uniq -c | \
            awk '$1==1 {{count++}} END {{print count+0}}')
        if [ "$TOTAL_GENES" -gt 0 ]; then
            SINGLE_PCT=$(awk "BEGIN {{printf \"%.0f\", ($SINGLE_EXON/$TOTAL_GENES)*100}}")
        else
            SINGLE_PCT=0
        fi
        echo "[INFO] Single-exon genes: $SINGLE_EXON / $TOTAL_GENES ($SINGLE_PCT%)"

        # Skip downsampling if >= 80% single-exon genes (organism likely has mostly intronless genes)
        if [ "$SINGLE_PCT" -ge 80 ]; then
            echo "[INFO] >= 80% single-exon genes detected — skipping downsampling"
            cp {output.gtf_for_downsample} {output.gtf_downsampled}
            echo "Skipped: $SINGLE_PCT% single-exon genes" > {output.downsample_log}
        else
            # Step 3: Downsample using Poisson distribution (lambda=2, BRAKER default)
            echo "[INFO] Downsampling genes with Poisson distribution (lambda=2)..."
            downsample_traingenes.pl \
                --in_gtf={output.gtf_for_downsample} \
                --out_gtf={output.gtf_downsampled} \
                --lambda=2 \
                1> {output.downsample_log} 2>&1
        fi

        # Step 4: Extract transcript IDs from downsampled GTF
        echo "[INFO] Extracting downsampled transcript list..."
        grep -v "^#" {output.gtf_downsampled} | \
            awk -F'\t' '{{print $9}}' | \
            perl -ne 'if (/transcript_id "([^"]+)"/) {{print "$1\n"}}' | \
            sort -u > {params.output_dir}/downsampled_transcripts.lst

        # Step 5: Map transcript names to locus IDs from GenBank
        echo "[INFO] Mapping transcripts to locus IDs..."
        cat {input.gb_diamond} | perl -ne '
            if ($_ =~ m/LOCUS\s+(\S+)\s/) {{
                $txLocus = $1;
            }} elsif ($_ =~ m/\/gene=\"(\S+)\"/) {{
                $txInGb3{{$1}} = $txLocus
            }}
            if(eof()) {{
                foreach (keys %txInGb3) {{
                    print "$_\t$txInGb3{{$_}}\n";
                }}
            }}' > {params.output_dir}/gene_to_locus.lst

        grep -f {params.output_dir}/downsampled_transcripts.lst {params.output_dir}/gene_to_locus.lst | cut -f2 > {params.output_dir}/downsampled.loci.lst || true

        # Step 6: Filter GenBank file to keep only downsampled genes
        echo "[INFO] Filtering GenBank file with downsampled gene list..."
        filterGenesIn.pl {params.output_dir}/downsampled.loci.lst {input.gb_diamond} > {output.gb_downsampled}

        # Count genes after downsampling
        GENES_AFTER=$(grep -c "^LOCUS" {output.gb_downsampled} || echo 0)
        GENES_REMOVED=$((GENES_BEFORE - GENES_AFTER))

        echo "[INFO] Genes before downsampling: $GENES_BEFORE"
        echo "[INFO] Genes after downsampling: $GENES_AFTER"
        echo "[INFO] Genes removed: $GENES_REMOVED"
        echo "[INFO] Retention rate: $(awk "BEGIN {{printf \"%.1f\", ($GENES_AFTER/$GENES_BEFORE)*100}}")%"
        echo "[INFO] ======================================="

        # Write statistics
        echo "$GENES_BEFORE -> $GENES_AFTER (removed: $GENES_REMOVED)" > {output.gene_count}

        if [ $GENES_AFTER -eq 0 ]; then
            echo "[ERROR] All genes were downsampled! Training cannot proceed."
            exit 1
        fi

        echo "[INFO] Downsampling completed successfully"
        echo "[INFO] Output: {output.gb_downsampled}"
        """
