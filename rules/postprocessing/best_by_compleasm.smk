"""
Optionally improve the TSEBRA-merged gene set using best_by_compleasm.

best_by_compleasm.py (from TSEBRA, patched for odb12 + library_path) runs
compleasm in protein mode on:
  - the TSEBRA-merged braker proteins
  - the AUGUSTUS proteins
  - the GeneMark proteins
and either keeps the merged set as-is, or rescues missing BUSCOs from the
AUGUSTUS / GeneMark sets that TSEBRA dropped.

For dual mode (short-read RNA-Seq + IsoSeq + proteins), the rule runs
best_by_compleasm THREE times:
  1. Rescue missing BUSCOs from the short-read GeneMark-ETP set
  2. Rescue missing BUSCOs from the IsoSeq GeneMark-ETP set
  3. Merge the two rescued sets via a third best_by_compleasm pass
     (using sr_rescued as "augustus" and iso_rescued as "genemark")

When run_best_by_compleasm = 0, this rule simply copies
braker.tsebra.raw.gtf -> braker.tsebra.gtf without running anything.
The default is 1 (enabled).

Container: teambraker/braker3:latest (contains compleasm, tsebra, getAnnoFastaFromJoingenes)
"""


def _bbc_genemark_inputs(wildcards):
    """Pick the GeneMark prediction GTF(s) appropriate for the sample's mode.

    Both genemark_sr and genemark_iso are ALWAYS returned. In non-dual modes,
    genemark_iso is aliased to genemark_sr — the shell block branch that
    actually reads genemark_iso only runs in dual mode, but Snakemake's string
    substitution scans every {input.X} reference at parse time and would raise
    AttributeError if the key were missing.
    """
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'es':
        sr = f"output/{sample}/GeneMark-ES/genemark.gtf"
    elif mode == 'ep':
        sr = f"output/{sample}/GeneMark-EP/genemark.gtf"
    elif mode == 'et':
        sr = f"output/{sample}/genemark/genemark.gtf"
    elif mode == 'dual':
        return {
            "genemark_sr":  f"output/{sample}/GeneMark-ETP/genemark.gtf",
            "genemark_iso": f"output/{sample}/GeneMark-ETP-isoseq/genemark.gtf",
        }
    else:  # etp / isoseq
        sr = f"output/{sample}/GeneMark-ETP/genemark.gtf"
    # Non-dual: alias genemark_iso to genemark_sr (the iso branch never runs).
    return {"genemark_sr": sr, "genemark_iso": sr}


rule best_by_compleasm:
    """Run best_by_compleasm to rescue missing BUSCOs from AUGUSTUS / GeneMark sets."""
    input:
        unpack(lambda w: {
            "braker_raw":  f"output/{w.sample}/braker.tsebra.raw.gtf",
            "augustus_gtf": f"output/{w.sample}/augustus.hints.fixed.gtf",
            "genome": os.path.join(get_braker_dir(w), "genome.fa"),
            **_bbc_genemark_inputs(w),
        })
    output:
        braker_merged_gtf="output/{sample}/braker.tsebra.gtf",
        bbc_log="output/{sample}/best_by_compleasm.log"
    log:
        "logs/{sample}/best_by_compleasm/best_by_compleasm.log"
    benchmark:
        "benchmarks/{sample}/best_by_compleasm/best_by_compleasm.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        run_bbc=config.get('run_best_by_compleasm', True),
        mode=lambda w: get_braker_mode(w.sample),
        busco_lineage=lambda w: get_busco_lineage(w),
        library_path=config['compleasm_download_path'],
        tmp_dir=lambda w: f"output/{w.sample}/best_by_compleasm_tmp",
        script=os.path.join(script_dir, "best_by_compleasm.py")
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1

        echo "[INFO] best_by_compleasm enabled: {params.run_bbc}" > {log}
        echo "[INFO] mode: {params.mode}" >> {log}

        if [ "{params.run_bbc}" != "True" ]; then
            # Disabled — pass through the raw TSEBRA output unchanged.
            echo "[INFO] run_best_by_compleasm=0; copying braker.tsebra.raw.gtf -> braker.tsebra.gtf" >> {log}
            cp {input.braker_raw} {output.braker_merged_gtf}
            echo "skipped (run_best_by_compleasm=0)" > {output.bbc_log}
            exit 0
        fi

        TMP={params.tmp_dir}
        rm -rf "$TMP"
        mkdir -p "$TMP"

        # ---- helper: run one best_by_compleasm pass ------------------------
        # $1 = stage label (e.g. "sr", "iso", "final")
        # $2 = path to braker.gtf for this pass
        # $3 = path to augustus.hints.gtf for this pass
        # $4 = path to genemark.gtf for this pass
        # writes: $TMP/<label>/better.gtf  (if produced)
        #         $TMP/<label>/bbc.stdout  (script stdout)
        # Note: augustus.hints.aa is ALWAYS regenerated from the staged
        # augustus.hints.gtf so the protein set matches the gene set the
        # pass is asked to evaluate. This matters for the dual-mode "final"
        # pass, where the staged augustus.hints.gtf is sr_rescued.gtf rather
        # than the original AUGUSTUS prediction.
        run_pass() {{
            local label="$1"
            local braker_gtf="$2"
            local aug_gtf="$3"
            local gm_gtf="$4"

            local stage="$TMP/${{label}}_stage"
            mkdir -p "$stage/GeneMark-ETP"
            cp "$braker_gtf" "$stage/braker.gtf"
            cp "$aug_gtf"    "$stage/augustus.hints.gtf"
            cp "$gm_gtf"     "$stage/GeneMark-ETP/genemark.gtf"

            # Generate braker.aa from staged braker.gtf
            python3 {script_dir}/getAnnoFastaFromJoingenes.py \
                -g {input.genome} \
                -f "$stage/braker.gtf" \
                -o "$stage/braker" \
                2>> {log} || true

            if [ ! -f "$stage/braker.aa" ]; then
                echo "[WARN] $label: could not generate braker.aa, skipping pass" >> {log}
                return 1
            fi

            # Generate augustus.hints.aa from staged augustus.hints.gtf so
            # the protein set is consistent with the GTF (critical for the
            # dual-mode "final" pass where the staged GTF is sr_rescued.gtf).
            python3 {script_dir}/getAnnoFastaFromJoingenes.py \
                -g {input.genome} \
                -f "$stage/augustus.hints.gtf" \
                -o "$stage/augustus.hints" \
                2>> {log} || true

            if [ ! -f "$stage/augustus.hints.aa" ]; then
                echo "[WARN] $label: could not generate augustus.hints.aa, skipping pass" >> {log}
                return 1
            fi

            mkdir -p "$TMP/$label"
            python3 {params.script} \
                -m "$TMP/$label" \
                -d "$stage" \
                -g {input.genome} \
                -t {threads} \
                -p {params.busco_lineage} \
                -L {params.library_path} \
                > "$TMP/$label/bbc.stdout" 2>> {log} || true

            # If the script wrote a better.gtf into $TMP/$label, copy it for clarity
            if [ -s "$TMP/$label/better.gtf" ]; then
                echo "[INFO] $label: produced better.gtf" >> {log}
            else
                echo "[INFO] $label: original was already best; no better.gtf" >> {log}
            fi
            return 0
        }}

        # ---- pick a result file from a pass --------------------------------
        # If better.gtf exists, prefer it; otherwise fall back to the input.
        pick_result() {{
            local label="$1"
            local fallback="$2"
            local target="$3"
            if [ -s "$TMP/$label/better.gtf" ]; then
                cp "$TMP/$label/better.gtf" "$target"
            else
                cp "$fallback" "$target"
            fi
        }}

        # ---------------------------------------------------------------- main
        if [ "{params.mode}" = "dual" ]; then
            echo "[INFO] dual mode: running 3 best_by_compleasm passes" >> {log}

            # Pass 1: rescue using SHORT-READ GeneMark-ETP
            run_pass "sr"  {input.braker_raw} {input.augustus_gtf} {input.genemark_sr}
            pick_result "sr"  {input.braker_raw} "$TMP/sr_rescued.gtf"

            # Pass 2: rescue using ISOSEQ GeneMark-ETP
            run_pass "iso" {input.braker_raw} {input.augustus_gtf} {input.genemark_iso}
            pick_result "iso" {input.braker_raw} "$TMP/iso_rescued.gtf"

            # Pass 3: merge the two rescued sets
            #   The script always wants three inputs (braker, augustus, genemark);
            #   we feed it the original braker as "braker", the SR-rescued set as
            #   "augustus", and the IsoSeq-rescued set as "genemark". run_pass
            #   regenerates augustus.hints.aa from sr_rescued.gtf so compleasm
            #   sees a protein set that matches the staged augustus.hints.gtf.
            run_pass "final" {input.braker_raw} "$TMP/sr_rescued.gtf" "$TMP/iso_rescued.gtf"
            pick_result "final" {input.braker_raw} {output.braker_merged_gtf}

            # Aggregate logs from all three passes into the master bbc_log
            {{
                echo "===== best_by_compleasm pass 1 (short-read) ====="
                cat "$TMP/sr/bbc.stdout" 2>/dev/null || true
                echo
                echo "===== best_by_compleasm pass 2 (IsoSeq) ====="
                cat "$TMP/iso/bbc.stdout" 2>/dev/null || true
                echo
                echo "===== best_by_compleasm pass 3 (merge) ====="
                cat "$TMP/final/bbc.stdout" 2>/dev/null || true
            }} > {output.bbc_log}

        else
            echo "[INFO] non-dual mode: running 1 best_by_compleasm pass" >> {log}
            run_pass "single" {input.braker_raw} {input.augustus_gtf} {input.genemark_sr}
            pick_result "single" {input.braker_raw} {output.braker_merged_gtf}
            cp "$TMP/single/bbc.stdout" {output.bbc_log} 2>/dev/null || true
            if [ ! -s {output.bbc_log} ]; then
                echo "best_by_compleasm produced no stdout (script may have failed)" > {output.bbc_log}
            fi
        fi

        # Final sanity check — ensure the output file exists
        if [ ! -f {output.braker_merged_gtf} ]; then
            cp {input.braker_raw} {output.braker_merged_gtf}
        fi

        # Clean up the staging area unless no_cleanup is set
        if [ "{config[no_cleanup]}" != "True" ]; then
            rm -rf "$TMP"
        fi

        # Citations
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite tsebra "$REPORT_DIR" || true
        cite compleasm "$REPORT_DIR" || true
        cite braker_book "$REPORT_DIR" || true
        """
