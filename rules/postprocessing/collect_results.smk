"""
Collect final results into a clean directory and remove intermediate files.

After the pipeline completes, the output directory contains hundreds of
intermediate files that are not useful to the end user. This rule copies
the important result files into output/{sample}/results/ and removes the
intermediate files.

The results directory contains everything a user needs:
- Gene predictions (GTF, GFF3, protein and CDS sequences)
- Quality control reports (BUSCO, compleasm, OMArk, gffcompare)
- Evidence support summary
"""


def _get_collect_inputs(wildcards):
    """Get all final output files that must exist before collection."""
    sample = wildcards.sample
    inputs = {
        "gtf": f"output/{sample}/braker.gtf",
        "aa": f"output/{sample}/braker.aa",
        "codingseq": f"output/{sample}/braker.codingseq",
        "gff3": f"output/{sample}/braker.gff3",
        "gene_support": f"output/{sample}/gene_support.tsv",
        "headers_fixed": f"output/{sample}/preprocessing/.headers_fixed"
    }

    if not config.get('skip_busco', False):
        inputs["busco_summary"] = f"output/{sample}/busco/busco_summary.txt"

    types = detect_data_types(sample)
    has_transcripts = any([types[k] for k in ['has_bam', 'has_fastq', 'has_sra', 'has_varus', 'has_isoseq']])
    if has_transcripts:
        inputs["utr_gtf"] = f"output/{sample}/braker_utr.gtf"

    # Compleasm always runs in BRAKER4. The protein-mode summary comes from
    # the best_by_compleasm rule's protein-mode invocation; the genome-mode
    # summary comes from the run_compleasm rule. We list the genome summary
    # as an unconditional input here so run_compleasm stays in the DAG even
    # when use_compleasm_hints=0 (otherwise nothing else would reference its
    # outputs and the rule would be dropped, breaking the report's QC section).
    inputs["compleasm"] = f"output/{sample}/compleasm_proteins/summary.txt"
    inputs["compleasm_genome"] = f"output/{sample}/compleasm_genome_out/summary.txt"

    if config.get('run_omark', False):
        inputs["omark"] = f"output/{sample}/omark/omark_summary.txt"

    if types.get('has_reference_gtf'):
        inputs["gffcompare"] = f"output/{sample}/gffcompare/gffcompare.stats"

    if config.get('run_fantasia', False):
        inputs["fantasia_results"]    = f"output/{sample}/fantasia/results.csv"
        inputs["fantasia_summary"]    = f"output/{sample}/fantasia/fantasia_summary.txt"
        inputs["fantasia_plot"]       = f"output/{sample}/fantasia/fantasia_go_categories.png"
        # Flat per-(transcript, GO term) table with human-readable GO names,
        # intended as the primary user-facing functional annotation file.
        inputs["fantasia_go_terms"]   = f"output/{sample}/fantasia/fantasia_go_terms.tsv"
        # GO-decorated GFF3 (Ontology_term=GO:... on mRNA + gene features).
        inputs["fantasia_go_gff3"]    = f"output/{sample}/fantasia/braker.go.gff3"
        if config.get('run_ncrna', False):
            inputs["fantasia_go_gff3_ncrna"] = (
                f"output/{sample}/fantasia/braker_with_ncRNA.go.gff3"
            )

    if config.get('run_ncrna', False):
        # All ncRNA tools (barrnap, tRNAscan-SE, Infernal, FEELnc) feed into
        # the merged braker_with_ncRNA.gff3 file produced by merge_rrna_into_gff3.
        inputs["gff3_ncrna"] = f"output/{sample}/braker_with_ncRNA.gff3"
        inputs["trna"] = f"output/{sample}/ncrna/tRNAs.gff3"
        inputs["infernal"] = f"output/{sample}/ncrna/ncRNAs_infernal.gff3"
        mode = get_braker_mode(sample)
        if mode in ('et', 'etp', 'isoseq', 'dual'):
            inputs["lncrna"] = f"output/{sample}/ncrna/lncRNAs.gff3"

    return inputs


rule collect_results:
    """Collect important output files and clean up intermediates."""
    input:
        unpack(_get_collect_inputs)
    output:
        done="output/{sample}/results/.done"
    benchmark:
        "benchmarks/{sample}/collect_results/collect_results.txt"
    params:
        sample="{sample}",
        outdir="output/{sample}",
        resultsdir="output/{sample}/results",
        no_cleanup=config.get('no_cleanup', False)
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        RESULTS="{params.resultsdir}"
        OUTDIR="{params.outdir}"

        mkdir -p "$RESULTS"
        mkdir -p "$RESULTS/quality_control"

        echo "[INFO] Collecting results for sample {params.sample}..."

        # --- Copy genome if pipeline modified it (masking or header cleaning) ---
        GENOME_PATH=$(readlink -f "$OUTDIR/genome.fa" 2>/dev/null || echo "unknown")
        echo "$GENOME_PATH" > "$RESULTS/genome_path.txt"
        COPY_GENOME=false
        if [ -f "$OUTDIR/preprocessing/.masking_complete" ]; then
            COPY_GENOME=true
        fi
        if [ -f "$OUTDIR/preprocessing/.headers_fixed" ] && grep -q "yes" "$OUTDIR/preprocessing/.headers_fixed" 2>/dev/null; then
            COPY_GENOME=true
        fi
        if [ "$COPY_GENOME" = "true" ] && [ -f "$OUTDIR/genome.fa" ]; then
            cp "$OUTDIR/genome.fa" "$RESULTS/"
            echo "[INFO] Including cleaned genome in results"
        fi

        # --- Core gene predictions (copy first, gzip AFTER report generation) ---
        for f in braker.gtf braker.gff3 braker_with_ncRNA.gff3 braker.aa braker.codingseq; do
            if [ -f "$OUTDIR/$f" ]; then
                cp "$OUTDIR/$f" "$RESULTS/"
            fi
        done

        # UTR-decorated GTF (if available)
        if [ -f "$OUTDIR/braker_utr.gtf" ]; then
            cp "$OUTDIR/braker_utr.gtf" "$RESULTS/"
        fi

        # --- Evidence support ---
        cp "$OUTDIR/gene_support.tsv" "$RESULTS/"
        cp "$OUTDIR/hintsfile.gff"    "$RESULTS/" 2>/dev/null || true

        # --- Quality control ---

        # BUSCO
        if [ -f "$OUTDIR/busco/busco_summary.txt" ]; then
            cp "$OUTDIR/busco/busco_summary.txt" "$RESULTS/quality_control/"
        fi
        # Copy full BUSCO output directories (genome + proteins short summaries).
        # Skipped entirely when BUSCO didn't run (skip_busco = 1) — otherwise the
        # `find` on a missing directory returns exit code 1, which pipefail
        # propagates through `| head -1`, killing the script under set -e.
        if [ -d "$OUTDIR/busco" ]; then
            for bmode in genome proteins; do
                if [ -d "$OUTDIR/busco/$bmode" ]; then
                    summary=$(find "$OUTDIR/busco/$bmode" -name "short_summary*.txt" 2>/dev/null | head -1 || true)
                    if [ -n "$summary" ] && [ -f "$summary" ]; then
                        cp "$summary" "$RESULTS/quality_control/busco_${{bmode}}_short_summary.txt"
                    fi
                fi
            done
        fi

        # Compleasm (proteome)
        if [ -f "$OUTDIR/compleasm_proteins/summary.txt" ]; then
            cp "$OUTDIR/compleasm_proteins/summary.txt" "$RESULTS/quality_control/compleasm_summary.txt"
        fi

        # Compleasm (genome) — preserve directory structure so report can find it
        if [ -f "$OUTDIR/compleasm_genome_out/summary.txt" ]; then
            mkdir -p "$RESULTS/quality_control/compleasm_genome_out"
            cp "$OUTDIR/compleasm_genome_out/summary.txt" "$RESULTS/quality_control/compleasm_genome_out/summary.txt"
        fi

        # OMArk
        if [ -d "$OUTDIR/omark" ]; then
            if [ -f "$OUTDIR/omark/omark_summary.txt" ]; then
                cp "$OUTDIR/omark/omark_summary.txt" "$RESULTS/quality_control/"
            fi
            # Also collect the detailed OMArk output files
            for f in "$OUTDIR"/omark/*_detailed_summary.txt "$OUTDIR"/omark/*.sum "$OUTDIR"/omark/*.tax; do
                if [ -f "$f" ]; then
                    cp "$f" "$RESULTS/quality_control/"
                fi
            done
        fi

        # gffcompare
        if [ -d "$OUTDIR/gffcompare" ]; then
            cp "$OUTDIR/gffcompare/gffcompare.stats" "$RESULTS/quality_control/" 2>/dev/null || true
            # Also collect the detailed gffcompare outputs
            for f in "$OUTDIR"/gffcompare/gffcompare.*.gtf "$OUTDIR"/gffcompare/gffcompare.tracking "$OUTDIR"/gffcompare/gffcompare.loci; do
                if [ -f "$f" ]; then
                    cp "$f" "$RESULTS/quality_control/"
                fi
            done
        fi

        # FANTASIA-Lite functional annotation (if run_fantasia=1)
        if [ -d "$OUTDIR/fantasia" ]; then
            mkdir -p "$RESULTS/quality_control/fantasia"
            for f in "$OUTDIR"/fantasia/results.csv \
                     "$OUTDIR"/fantasia/fantasia_summary.txt \
                     "$OUTDIR"/fantasia/fantasia_go_categories.png \
                     "$OUTDIR"/fantasia/failed_sequences.csv; do
                if [ -f "$f" ]; then
                    cp "$f" "$RESULTS/quality_control/fantasia/"
                fi
            done
            # topGO-compatible per-namespace exports for downstream R analysis
            if [ -d "$OUTDIR/fantasia/topgo" ]; then
                cp -r "$OUTDIR/fantasia/topgo" "$RESULTS/quality_control/fantasia/"
            fi
            # GO-decorated GFF3 files and the flat per-transcript GO term
            # table belong at top-level results/ alongside the standard
            # braker.gff3.gz, since they are primary annotation deliverables.
            for f in "$OUTDIR"/fantasia/braker.go.gff3 \
                     "$OUTDIR"/fantasia/braker_with_ncRNA.go.gff3 \
                     "$OUTDIR"/fantasia/fantasia_go_terms.tsv; do
                if [ -f "$f" ]; then
                    cp "$f" "$RESULTS/"
                fi
            done
            echo "[INFO] Collected FANTASIA-Lite functional annotation"
        fi

        # ncRNA annotations (if run_ncrna=1)
        if [ -d "$OUTDIR/ncrna" ]; then
            mkdir -p "$RESULTS/ncrna"
            for f in "$OUTDIR"/ncrna/tRNAs.gff3 "$OUTDIR"/ncrna/ncRNAs_infernal.gff3 \
                     "$OUTDIR"/ncrna/lncRNAs.gff3 "$OUTDIR"/ncrna/feelnc_classifier.txt \
                     "$OUTDIR"/ncrna/tRNAs.txt "$OUTDIR"/ncrna/infernal.tblout; do
                if [ -f "$f" ]; then
                    cp "$f" "$RESULTS/ncrna/"
                fi
            done
            echo "[INFO] Collected ncRNA annotations"
        fi

        # VARUS run list and stats (if VARUS was used)
        if [ -d "$OUTDIR/varus" ]; then
            for f in "$OUTDIR"/varus/varus_runlist.tsv "$OUTDIR"/varus/varus_stats.txt; do
                if [ -f "$f" ]; then
                    cp "$f" "$RESULTS/"
                fi
            done
        fi

        # --- Software versions ---
        if [ -f "$OUTDIR/software_versions.tsv" ]; then
            # Deduplicate (keep first occurrence of each tool) and sort
            awk -F'\t' '!seen[$1]++' "$OUTDIR/software_versions.tsv" | sort -f > "$RESULTS/software_versions.tsv"
            echo "[INFO] Collected software versions"
        fi

        # --- Training summary (must run BEFORE cleanup removes intermediate files) ---
        echo "[INFO] Generating training summary..."
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        python3 {script_dir}/training_summary.py \
            -d "$OUTDIR" \
            -o "$RESULTS/quality_control" \
            2>/dev/null || echo "[WARNING] Training summary generation failed (non-fatal)"

        # --- Gene set statistics and plots (must run BEFORE cleanup) ---
        echo "[INFO] Generating gene set statistics..."
        SUPPORT_ARG=""
        if [ -f "$OUTDIR/gene_support.tsv" ]; then
            SUPPORT_ARG="-s $OUTDIR/gene_support.tsv"
        fi
        python3 {script_dir}/gene_set_statistics.py \
            -g "$OUTDIR/braker.gtf" \
            -o "$RESULTS/quality_control" \
            $SUPPORT_ARG \
            2>/dev/null || echo "[WARNING] Gene set statistics generation failed (non-fatal)"

        # --- Completeness plot (BUSCO + compleasm combined) ---
        echo "[INFO] Generating completeness plot..."
        python3 {script_dir}/completeness_plot.py \
            -d "$RESULTS/quality_control" \
            -o "$RESULTS/quality_control/completeness.png" \
            2>/dev/null || echo "[WARNING] Completeness plot generation failed (non-fatal)"

        # --- Generate HTML report (must run BEFORE cleanup) ---
        echo "[INFO] Generating HTML report..."
        python3 {script_dir}/generate_report.py \
            -d "$OUTDIR" \
            -o "$RESULTS" \
            -s "{params.sample}" \
            2>/dev/null || echo "[WARNING] HTML report generation failed (non-fatal)"

        # Copy BibTeX to results if generated
        if [ -f "$RESULTS/braker_citations.bib" ]; then
            echo "[INFO] BibTeX citations: $RESULTS/braker_citations.bib"
        fi

        # --- Gzip large result files (AFTER report generation used them) ---
        echo "[INFO] Compressing result files..."
        for f in "$RESULTS"/braker.gtf "$RESULTS"/braker.gff3 \
                 "$RESULTS"/braker_with_ncRNA.gff3 \
                 "$RESULTS"/braker.go.gff3 "$RESULTS"/braker_with_ncRNA.go.gff3 \
                 "$RESULTS"/fantasia_go_terms.tsv \
                 "$RESULTS"/braker.aa "$RESULTS"/braker.codingseq "$RESULTS"/braker_utr.gtf \
                 "$RESULTS"/hintsfile.gff "$RESULTS"/genome.fa \
                 "$RESULTS"/quality_control/fantasia/results.csv; do
            if [ -f "$f" ]; then
                gzip -f "$f"
            fi
        done

        # --- Clean up intermediate files (unless no_cleanup is set) ---
        if [ "{params.no_cleanup}" = "True" ]; then
            echo "[INFO] Skipping cleanup (no_cleanup=1 in config.ini)"
        else
            echo "[INFO] Removing intermediate files..."
            for item in "$OUTDIR"/*; do
                basename=$(basename "$item")
                if [ "$basename" = "results" ]; then
                    continue
                fi
                rm -rf "$item"
            done

            # Clean up stray log files written to the scenario working directory
            # by AGAT, BUSCO, and other tools that ignore output path settings
            SCENARIO_WD=$(pwd)
            rm -f "$SCENARIO_WD"/*.agat.log "$SCENARIO_WD"/busco_*.log 2>/dev/null || true
        fi

        touch "$RESULTS/.done"

        # List final contents
        echo "[INFO] Results collected in $RESULTS:"
        find "$RESULTS" -type f -not -name '.done' | sort | while read f; do
            size=$(du -h "$f" | cut -f1)
            rel=$(echo "$f" | sed "s|$RESULTS/||")
            echo "  $size  $rel"
        done
        """
