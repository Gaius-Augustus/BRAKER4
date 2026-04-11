"""
Merge all ncRNA annotations into the final GFF3.

When run_ncrna=1, this rule replaces merge_rrna_into_gff3 in the DAG.
It merges protein-coding genes (braker.gff3) with:
  - rRNA (barrnap, always present)
  - tRNA (tRNAscan-SE)
  - ncRNA (Infernal/Rfam)
  - lncRNA (FEELnc, only when transcript evidence available)

ncRNA features that overlap protein-coding exons by >50% are dropped.
Output: braker_with_ncRNA.gff3
"""


def _get_ncrna_merge_inputs(wildcards):
    """Collect all ncRNA GFF3 files available for this sample."""
    sample = wildcards.sample
    inputs = {
        "gff3": f"output/{sample}/braker.gff3",
        "rrna": f"output/{sample}/barrnap/rRNA.gff3",
        "trna": f"output/{sample}/ncrna/tRNAs.gff3",
        "infernal": f"output/{sample}/ncrna/ncRNAs_infernal.gff3",
    }

    # lncRNA only available when transcript evidence exists
    mode = get_braker_mode(sample)
    if mode in ('et', 'etp', 'isoseq', 'dual'):
        inputs["lncrna"] = f"output/{sample}/ncrna/lncRNAs.gff3"

    return inputs


rule merge_ncrna_into_gff3:
    """Merge all ncRNA predictions into the final GFF3 annotation.

    Priority (kept on overlap): protein-coding > rRNA > tRNA > Infernal > lncRNA.
    ncRNA features overlapping protein-coding exons by >50% are dropped.
    """
    input:
        unpack(_get_ncrna_merge_inputs)
    output:
        merged="output/{sample}/braker_with_ncRNA.gff3"
    log:
        "logs/{sample}/ncrna/merge_ncrna.log"
    benchmark:
        "benchmarks/{sample}/ncrna/merge_ncrna.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] Merging ncRNA annotations into final GFF3..." > {log}

        # Start with protein-coding BRAKER GFF3
        cp {input.gff3} {output.merged}
        n_braker=$(grep -cv '^#' {input.gff3} || echo 0)
        echo "[INFO] Base: $n_braker protein-coding features" >> {log}

        # Extract protein-coding exon coordinates for overlap filtering
        # Format: seqid\tstart\tend (1-based)
        EXON_BED=$(mktemp)
        awk -F'\t' '$3 == "exon" || $3 == "CDS" {{print $1 "\t" $4 "\t" $5}}' {input.gff3} | \
            sort -k1,1 -k2,2n > "$EXON_BED"

        # Function: append ncRNA features, filtering those with >50% overlap
        # with protein-coding exons
        append_filtered() {{
            local ncrna_gff="$1"
            local label="$2"
            local added=0
            local dropped=0

            if [ ! -s "$ncrna_gff" ] || ! grep -qv '^#' "$ncrna_gff" 2>/dev/null; then
                echo "[INFO] $label: 0 features (empty file)" >> {log}
                return
            fi

            while IFS=$'\t' read -r seq src type start end score strand phase attrs; do
                [ -z "$seq" ] && continue
                [[ "$seq" == \#* ]] && continue

                # Calculate overlap with protein-coding exons on same sequence
                feat_len=$((end - start + 1))
                overlap=0

                while IFS=$'\t' read -r eseq estart eend; do
                    # Calculate overlap between [start,end] and [estart,eend]
                    ostart=$((start > estart ? start : estart))
                    oend=$((end < eend ? end : eend))
                    if [ "$oend" -ge "$ostart" ]; then
                        overlap=$((overlap + oend - ostart + 1))
                    fi
                done < <(awk -v s="$seq" -v a="$start" -v b="$end" \
                    '$1 == s && $3 >= a && $2 <= b' "$EXON_BED")

                # Drop if >50% of ncRNA feature overlaps coding exons
                if [ "$feat_len" -gt 0 ] && [ "$overlap" -gt 0 ]; then
                    pct=$((overlap * 100 / feat_len))
                    if [ "$pct" -gt 50 ]; then
                        dropped=$((dropped + 1))
                        continue
                    fi
                fi

                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                    "$seq" "$src" "$type" "$start" "$end" "$score" "$strand" "$phase" "$attrs" \
                    >> {output.merged}
                added=$((added + 1))
            done < "$ncrna_gff"

            echo "[INFO] $label: $added added, $dropped dropped (>50% overlap with CDS)" >> {log}
        }}

        # Append in priority order: rRNA > tRNA > Infernal > lncRNA
        append_filtered "{input.rrna}" "rRNA (barrnap)"
        append_filtered "{input.trna}" "tRNA (tRNAscan-SE)"
        append_filtered "{input.infernal}" "ncRNA (Infernal/Rfam)"

        # lncRNA is optional (only when transcript evidence available)
        LNCRNA_FILE=""
        if [ -f "output/{wildcards.sample}/ncrna/lncRNAs.gff3" ]; then
            LNCRNA_FILE="output/{wildcards.sample}/ncrna/lncRNAs.gff3"
        fi
        if [ -n "$LNCRNA_FILE" ]; then
            append_filtered "$LNCRNA_FILE" "lncRNA (FEELnc)"
        fi

        rm -f "$EXON_BED"

        n_total=$(grep -cv '^#' {output.merged} || echo 0)
        echo "[INFO] Final merged GFF3: $n_total total features" >> {log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        n_ncrna=$((n_total - n_braker))
        """
