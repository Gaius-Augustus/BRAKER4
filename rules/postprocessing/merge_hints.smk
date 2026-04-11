
def _get_merge_hints_inputs(wildcards):
    """Get mode-dependent hint sources for merging."""
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    types = detect_data_types(sample)

    inputs = {}

    if mode == 'dual':
        # Dual mode: merged hints from both ETP runs
        inputs['etp_hints'] = f"output/{sample}/dual_etp_merged/etp_hints.gff"
    elif mode in ('etp', 'isoseq'):
        # ETP/IsoSeq mode: hints come from GeneMark-ETP's get_etp_hints.py
        inputs['etp_hints'] = f"output/{sample}/etp_hints.gff"
    else:
        # ET/EP mode: separate hint sources
        has_rnaseq = any([types[k] for k in ['has_bam', 'has_fastq', 'has_sra', 'has_varus']])
        if has_rnaseq:
            inputs['bam2hints'] = f"output/{sample}/bam2hints.gff"
        if types['has_proteins']:
            inputs['prothint'] = f"output/{sample}/prothint_hints.gff"

    # Compleasm-derived CDSpart hints (optional, for all modes). Compleasm
    # ALWAYS runs because best_by_compleasm needs it for BUSCO rescue, but
    # use_compleasm_hints=0 keeps the hints out of the AUGUSTUS hintsfile.
    if config.get('use_compleasm_hints', True):
        inputs['compleasm_hints'] = f"output/{sample}/compleasm_hints.gff"

    return inputs

rule merge_hints:
    """
    Merge all hints files into a single hintsfile for AUGUSTUS, then collapse
    duplicate hints into multiplicity hints via join_mult_hints.pl.

    Mode-dependent behavior:
    - ET mode: BAM2hints + Compleasm
    - EP mode: ProtHint + Compleasm
    - ETP mode: GeneMark-ETP hints (already combined) + Compleasm

    The post-merge join_mult_hints.pl step (matches braker.pl sub join_mult_hints,
    line ~4760) is essential: without it, identical intron hints from different
    sources (RNA-seq + protein + compleasm) would be counted multiple times by
    AUGUSTUS, inflating their effective multiplicity scores.

    PARITY NOTE: braker.pl's join_mult_hints (scripts/braker.pl line 4769)
    explicitly EXCLUDES hints with src=C AND a grp= (or group=) tag from the
    join_mult_hints.pl pass. Those hints carry per-gene group linkage that
    AUGUSTUS uses to score whole-gene hint coverage; collapsing them by
    coordinate destroys the linkage and silently degrades TSEBRA's
    transcript-level intron-support scoring (which gates the
    intron_support=1.0 threshold in braker3.cfg). On A. thaliana the wrong
    treatment costs ~3 percentage points of locus-level sensitivity.

    Container: teambraker/braker3:latest (provides join_mult_hints.pl)

    Output:
        hintsfile: Merged + collapsed hints file in GFF format
        hints_stats: Statistics about hint types and sources
    """
    input:
        unpack(_get_merge_hints_inputs)
    output:
        hintsfile = "output/{sample}/hintsfile.gff",
        hints_stats = "output/{sample}/hints_statistics.txt"
    benchmark:
        "benchmarks/{sample}/merge_hints/merge_hints.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ========== MERGING HINTS FILES ==========" | tee {output.hints_stats}
        echo "[INFO] Merging hints from multiple sources..." | tee -a {output.hints_stats}

        # Concatenate all input hint files into a temporary unsorted file
        TMP_RAW={output.hintsfile}.raw
        > "$TMP_RAW"
        for hints_file in {input}; do
            if [ -f "$hints_file" ] && [ -s "$hints_file" ]; then
                cat "$hints_file" >> "$TMP_RAW"
                echo "[INFO] Added $(basename $hints_file): $(wc -l < $hints_file) hints" | tee -a {output.hints_stats}
            fi
        done

        # Split hints into two streams (mirrors braker.pl sub join_mult_hints
        # at scripts/braker.pl line 4785-4791):
        #   TMP_NOMERGE - hints with src=C AND a grp=/group= tag. These carry
        #                 per-gene group linkage that AUGUSTUS uses for
        #                 whole-gene hint scoring; coordinate-collapsing them
        #                 destroys the grp= identity and silently breaks
        #                 TSEBRA's intron_support gating. They are passed
        #                 through to the final hintsfile UNCHANGED.
        #   TMP_MERGE   - all other hints. These get sorted and run through
        #                 join_mult_hints.pl to collapse identical lines into
        #                 multiplicity-counted hints.
        # Final hintsfile is the concatenation of (joined merge stream) +
        # (passthrough nomerge stream), matching braker.pl line 4853.
        N_RAW=$(wc -l < "$TMP_RAW" || echo 0)
        TMP_MERGE={output.hintsfile}.merge
        TMP_NOMERGE={output.hintsfile}.nomerge
        awk -F'\t' '{{
            if ($9 ~ /src=C/ && ($9 ~ /grp=/ || $9 ~ /group=/)) {{
                print > "'"$TMP_NOMERGE"'"
            }} else {{
                print > "'"$TMP_MERGE"'"
            }}
        }}' "$TMP_RAW"
        touch "$TMP_MERGE" "$TMP_NOMERGE"
        N_NOMERGE=$(wc -l < "$TMP_NOMERGE")
        N_MERGE_IN=$(wc -l < "$TMP_MERGE")
        echo "[INFO] split: $N_NOMERGE src=C grp=* hints (passthrough), $N_MERGE_IN merge candidates" | tee -a {output.hints_stats}

        # Sort key order matches braker.pl sub join_mult_hints
        # (scripts/braker.pl line 4820): coordinate-based, then stable-sorted
        # by type, then by chromosome.
        TMP_JOINED={output.hintsfile}.joined
        cat "$TMP_MERGE" \
            | sort -n -k 4,4 \
            | sort -s -n -k 5,5 \
            | sort -s -n -k 3,3 \
            | sort -s -k 1,1 \
            | join_mult_hints.pl \
            > "$TMP_JOINED"
        N_JOINED_OUT=$(wc -l < "$TMP_JOINED")
        echo "[INFO] join_mult_hints.pl: $N_MERGE_IN -> $N_JOINED_OUT collapsed hints" | tee -a {output.hints_stats}

        # Concatenate joined merge stream + passthrough nomerge stream
        cat "$TMP_JOINED" "$TMP_NOMERGE" > {output.hintsfile}
        rm -f "$TMP_RAW" "$TMP_MERGE" "$TMP_NOMERGE" "$TMP_JOINED"
        N_FINAL=$(wc -l < {output.hintsfile} || echo 0)
        echo "[INFO] Final hintsfile: $N_RAW raw -> $N_FINAL hints (joined+passthrough)" | tee -a {output.hints_stats}

        # Count hints by type and source
        echo "" | tee -a {output.hints_stats}
        echo "===== HINTS STATISTICS =====" | tee -a {output.hints_stats}
        echo "" | tee -a {output.hints_stats}

        echo "Total hints: $(grep -c '^' {output.hintsfile} || echo 0)" | tee -a {output.hints_stats}
        echo "" | tee -a {output.hints_stats}

        echo "By hint type:" | tee -a {output.hints_stats}
        awk '{{print $3}}' {output.hintsfile} | sort | uniq -c | sort -rn | tee -a {output.hints_stats}
        echo "" | tee -a {output.hints_stats}

        # POSIX-portable: awk in the BRAKER container is mawk, which does
        # NOT support gawk's 3-arg match($str, regex, array). Use sub() to
        # extract the src= field instead.
        echo "By source:" | tee -a {output.hints_stats}
        awk '$9 ~ /src=/ {{ s=$9; sub(/.*src=/, "", s); sub(/;.*/, "", s); print s }}' {output.hintsfile} | sort | uniq -c | sort -rn | tee -a {output.hints_stats}
        echo "" | tee -a {output.hints_stats}

        echo "Breakdown by source and type:" | tee -a {output.hints_stats}
        awk '{{
            type=$3;
            if ($9 ~ /src=/) {{
                src=$9;
                sub(/.*src=/, "", src);
                sub(/;.*/, "", src);
            }} else {{
                src="unknown";
            }}
            print src"\t"type
        }}' {output.hintsfile} | sort | uniq -c | sort -k2,2 -k3,3 | tee -a {output.hints_stats}

        echo "" | tee -a {output.hints_stats}
        echo "[INFO] Hints merging completed: {output.hintsfile}" | tee -a {output.hints_stats}
        """
