"""
Filter GeneMark predictions to extract high-quality genes for AUGUSTUS training.

Filters GeneMark predictions based on:
- Number of exons (multi-exon genes preferred)
- Intron support from RNA-Seq hints
- CDS length thresholds
- Gene completeness

Input:
    - genemark.gtf - Raw GeneMark predictions
    - hintsfile.gff - RNA-Seq hints for validation

Output:
    - genemark.f.good.gtf - High-quality genes for training
    - genemark.f.bad.gtf - Filtered out genes
    - genemark.average_gene_length.out - Statistics

Container: teambraker/braker3:latest (contains filterGenemark.pl)
"""

rule filter_genemark:
    input:
        gtf=get_genemark_output,
        hints=get_genemark_hints_for_filter
    output:
        good_gtf="output/{sample}/genemark/genemark.f.good.gtf",
        bad_gtf="output/{sample}/genemark/genemark.f.bad.gtf",
        stats="output/{sample}/genemark/genemark.average_gene_length.out"
    log:
        "logs/{sample}/filter_genemark/filter_genemark.log"
    benchmark:
        "benchmarks/{sample}/filter_genemark/filter_genemark.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        mode=lambda wildcards: get_braker_mode(wildcards.sample)
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        n_total=$(awk -F'\t' '$3=="gene"' {input.gtf} | wc -l)

        # Sort the GeneMark output by coordinates and renumber gene/transcript
        # IDs in ascending order. This makes filterGenemark.pl's tiebreaks
        # deterministic across runs. Mirrors braker.pl's sortGeneMark() call
        # before its filterGenemark.pl invocation.
        GENEMARK_DIR=$(dirname {output.good_gtf})
        mkdir -p "$GENEMARK_DIR"
        SORTED_GTF="$GENEMARK_DIR/genemark.sorted.gtf"
        sortGeneMark.py {input.gtf} > "$SORTED_GTF" 2>> {log}

        # In ES mode there is no evidence to filter by — use the sorted set
        # directly as the "good" output. (BRAKER4 does not currently apply the
        # CDS-length > 800 nt filter that braker.pl does in ES mode; see audit.)
        if [ "{params.mode}" = "es" ]; then
            echo "[INFO] ES mode: skipping evidence-based filtering, using all $n_total GeneMark genes" | tee -a {log}
            cp "$SORTED_GTF" {output.good_gtf}
            # No bad genes in ES mode
            : > {output.bad_gtf}
        else
            # filterGenemark.pl strips ".gtf" from --output= to get a prefix and
            # then writes FOUR files:
            #   <prefix>.gtf                       (the converted full set)
            #   <prefix>.f.good.gtf                (genes that pass the hint filter)
            #   <prefix>.f.bad.gtf                 (genes that fail)
            #   <prefix>.f.multi_anchored.gtf     (multi-anchored subset)
            # We need <prefix>.f.good.gtf as our declared good_gtf output, so we
            # call filterGenemark.pl with a TEMPORARY stem and then mv the right
            # files into place. Passing --output={output.good_gtf} directly is
            # WRONG: it overwrites good_gtf with the converted full set, not the
            # filtered subset.
            TEMP_STEM="$GENEMARK_DIR/_filtergenemark"

            # --randomSeed=1 makes the random tiebreaks deterministic across
            # runs (matches braker.pl).
            filterGenemark.pl \
                --genemark="$SORTED_GTF" \
                --hints={input.hints} \
                --output="$TEMP_STEM.gtf" \
                --randomSeed=1 \
                2>&1 | tee -a {log}

            # Move the actually-filtered subsets into place
            mv "$TEMP_STEM.f.good.gtf" {output.good_gtf}
            mv "$TEMP_STEM.f.bad.gtf"  {output.bad_gtf}

            # Clean up the side outputs we don't consume
            rm -f "$TEMP_STEM.gtf" "$TEMP_STEM.f.multi_anchored.gtf"

            # Top up the good set if it has fewer than 4000 genes by promoting
            # the most-supported "bad" genes back into "good". This matches
            # braker.pl's ensure_n_training_genes.py call (line ~5797). The
            # script edits both --goodGenes and --badGenes IN PLACE.
            ensure_n_training_genes.py \
                --goodGenes {output.good_gtf} \
                --badGenes  {output.bad_gtf} \
                --N 4000 \
                2>&1 | tee -a {log} || true
        fi

        # Clean up the temporary sorted file
        rm -f "$SORTED_GTF"

        # Calculate statistics from the (now correctly filtered) outputs
        n_good=$(awk -F'\t' '$3=="gene"' {output.good_gtf} | wc -l)
        n_bad=$((n_total - n_good))

        echo "Total GeneMark genes: $n_total" >> {log}
        echo "Good genes for training: $n_good" >> {log}
        echo "Filtered out: $n_bad" >> {log}

        # Calculate average gene length
        total_length=0
        count=0
        while read line; do
            if echo "$line" | grep -q $'\\tgene\\t'; then
                start=$(echo "$line" | cut -f4)
                end=$(echo "$line" | cut -f5)
                length=$((end - start + 1))
                total_length=$((total_length + length))
                count=$((count + 1))
            fi
        done < {output.good_gtf}

        if [ $count -gt 0 ]; then
            avg=$((total_length / count))
            echo "$avg" > {output.stats}
            echo "Average gene length: $avg bp" >> {log}
        else
            echo "0" > {output.stats}
            echo "WARNING: No good genes found!" >> {log}
        fi

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh || true
        """
