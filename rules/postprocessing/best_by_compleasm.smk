"""
Optionally improve the TSEBRA-merged gene set using best_by_compleasm.

best_by_compleasm.py (from TSEBRA, patched for odb12 + library_path) runs
compleasm in protein mode on:
  - the TSEBRA-merged braker proteins
  - the AUGUSTUS proteins
  - the GeneMark proteins
and either keeps the merged set as-is, or rescues missing BUSCOs from the
AUGUSTUS / GeneMark sets that TSEBRA dropped.

In dual mode (short-read RNA-Seq + IsoSeq + proteins) GeneMark-ETP runs
once with both read types (GeneMark-ETP-isoseq/), so there is one pass, as
in all other modes.

When run_best_by_compleasm = 0, this rule simply copies
braker.tsebra.raw.gtf -> braker.tsebra.gtf without running anything.
The default is 1 (enabled).

Container: teambraker/braker3:latest (contains compleasm, tsebra, getAnnoFastaFromJoingenes)
"""


def _bbc_genemark_gtf(wildcards):
    """GeneMark prediction GTF of the sample's mode."""
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'es':
        return f"output/{sample}/GeneMark-ES/genemark.gtf"
    elif mode == 'ep':
        return f"output/{sample}/GeneMark-EP/genemark.gtf"
    elif mode == 'et':
        return f"output/{sample}/genemark/genemark.gtf"
    elif mode == 'dual':
        return f"output/{sample}/GeneMark-ETP-isoseq/genemark.gtf"
    return f"output/{sample}/GeneMark-ETP/genemark.gtf"   # etp / isoseq


rule best_by_compleasm:
    """Run best_by_compleasm to rescue missing BUSCOs from AUGUSTUS / GeneMark sets.

    Scratch: the staging area (per-pass gene sets, compleasm protein-mode
    runs) is a private directory on the node-local disk
    (scripts/tmp_dir.sh, [paths] tmp_dir); only the outputs are written to
    the run directory, plus the whole staging area to
    output/<sample>/best_by_compleasm_tmp when no_cleanup = 1. NEED 10 GB;
    with less free the job works in output/<sample>/best_by_compleasm_tmp as
    before.
    """
    input:
        unpack(lambda w: {
            "braker_raw":  f"output/{w.sample}/braker.tsebra.raw.gtf",
            "augustus_gtf": f"output/{w.sample}/augustus.hints.fixed.gtf",
            "genome": os.path.join(get_braker_dir(w), "genome.fa"),
            "genemark_gtf": _bbc_genemark_gtf(w),
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
        script=os.path.join(script_dir, "best_by_compleasm.py"),
        tmp_root=TMP_ROOT
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

        # Staging area on the node-local disk
        source {script_dir}/tmp_dir.sh
        scratch_dir outDir "bbc_{wildcards.sample}" "{params.tmp_root}" 10 \
            "$PWD/{params.tmp_dir}" 2>> {log}
        trap 'rm -rf -- "$SCRATCH"' EXIT
        TMP=$outDir
        if [ -z "$SCRATCH" ]; then
            # fallback: clear an earlier attempt's staging area
            rm -rf "$TMP"
            mkdir -p "$TMP"
        fi

        # ---- helper: run one best_by_compleasm pass ------------------------
        # $1 = stage label (e.g. "sr", "iso", "final")
        # $2 = path to braker.gtf for this pass
        # $3 = path to augustus.hints.gtf for this pass
        # $4 = path to genemark.gtf for this pass
        # writes: $TMP/<label>/better.gtf  (if produced)
        #         $TMP/<label>/bbc.stdout  (script stdout)
        #         $TMP/<label>/bbc_status.txt  (one line: "OK: ...",
        #             "COMPLEASM_FAILED: ..." or "FAILED: ...")
        # On a failed pass, better.gtf is discarded so pick_result falls back
        # to the input gene set; the pipeline keeps running with a warning.
        # Note: augustus.hints.aa is ALWAYS regenerated from the staged
        # augustus.hints.gtf so the protein set matches the gene set the
        # pass is asked to evaluate.
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
            # the protein set is consistent with the GTF.
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
            local status_file="$TMP/$label/bbc_status.txt"
            local rc=0
            python3 {params.script} \
                -m "$TMP/$label" \
                -d "$stage" \
                -g {input.genome} \
                -t {threads} \
                -p {params.busco_lineage} \
                -L {params.library_path} \
                -s "$status_file" \
                > "$TMP/$label/bbc.stdout" 2>> {log} || rc=$?

            # The script writes its own status line; if it died before that
            # (killed, OOM, import error), record the exit code instead.
            if [ ! -s "$status_file" ]; then
                echo "FAILED: best_by_compleasm.py exited with code $rc without writing a status line" > "$status_file"
            fi
            local status
            status=$(head -n 1 "$status_file")
            echo "[INFO] $label: status: $status" >> {log}

            case "$status" in
                OK*)
                    if [ -s "$TMP/$label/better.gtf" ]; then
                        echo "[INFO] $label: produced better.gtf" >> {log}
                    else
                        echo "[INFO] $label: original was already best; no better.gtf" >> {log}
                    fi
                    ;;
                *)
                    # Never use a better.gtf left behind by a failed run
                    rm -f "$TMP/$label/better.gtf"
                    echo "[WARNING] $label: best_by_compleasm failed, keeping the TSEBRA gene set unchanged. $status" | tee -a {log} >&2
                    ;;
            esac
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
        echo "[INFO] running best_by_compleasm" >> {log}
        run_pass "single" {input.braker_raw} {input.augustus_gtf} {input.genemark_gtf}
        pick_result "single" {input.braker_raw} {output.braker_merged_gtf}
        # First line of the log is the run status (read by generate_report.py)
        {{
            echo "STATUS: $(head -n 1 "$TMP/single/bbc_status.txt")"
            grep -v '^STATUS: ' "$TMP/single/bbc.stdout" 2>/dev/null || true
        }} > {output.bbc_log}

        # Final sanity check — ensure the output file exists
        if [ ! -f {output.braker_merged_gtf} ]; then
            cp {input.braker_raw} {output.braker_merged_gtf}
        fi

        # Clean up the staging area unless no_cleanup is set. On scratch the
        # EXIT trap removes it; with no_cleanup it is copied to the run dir.
        if [ "{config[no_cleanup]}" = "True" ]; then
            copy_back "$TMP" "$PWD/{params.tmp_dir}"
        elif [ -z "$SCRATCH" ]; then
            rm -rf "$TMP"
        fi

        # Citations
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite tsebra "$REPORT_DIR" || true
        cite compleasm "$REPORT_DIR" || true
        cite braker_book "$REPORT_DIR" || true
        """
