"""
Optional functional annotation of BRAKER4 predicted proteins with FANTASIA-Lite.

FANTASIA-Lite assigns GO terms to each predicted protein using ProtT5
(prot_t5_xl_uniref50) protein language model embeddings and a bundled lookup
table of pre-computed reference embeddings. There is no PostgreSQL, no
RabbitMQ, and no FANTASIA repo to clone -- everything required for inference
ships inside the container.

This step is OFF BY DEFAULT and is the most fragile component of BRAKER4:
the FANTASIA-Lite container hard-requires an NVIDIA GPU with --nv. The
embedding step has only been validated on an A100 here. CPU-only execution is
not supported by the upstream container. See README.md (run_fantasia section)
for the warnings.

The container, the singularity invocation, and the FANTASIA-Lite CLI flags
mirror the validated invocation from the EukAssembly-Bin (BOUDICCA) workflow,
which has gone through extensive debugging on Hoff lab GPUs. Do not change
those flags casually.

Two rules:
    - fantasia_annotate:  GPU embedding + GO lookup, produces results.csv
    - fantasia_summarize: parses results.csv, writes summary.txt + GO bar plot
"""


FANTASIA_SIF        = config['fantasia']['sif']
FANTASIA_HF_CACHE   = config['fantasia']['hf_cache_dir']
FANTASIA_ADD_PARAMS = config['fantasia'].get('additional_params', '') or ''
FANTASIA_MIN_SCORE  = float(config['fantasia'].get('min_score', 0.5))


rule fantasia_annotate:
    """Embed proteins with ProtT5 and assign GO terms via FANTASIA-Lite (GPU)."""
    input:
        proteins="output/{sample}/braker.aa"
    output:
        results="output/{sample}/fantasia/results.csv",
        done="output/{sample}/fantasia/.fantasia_done"
    log:
        "logs/{sample}/fantasia/fantasia_annotate.log"
    benchmark:
        "benchmarks/{sample}/fantasia/fantasia_annotate.txt"
    params:
        sif=FANTASIA_SIF,
        hf_cache=FANTASIA_HF_CACHE,
        add_params=FANTASIA_ADD_PARAMS,
        outdir=lambda wc: f"output/{wc.sample}/fantasia"
    threads:
        int(config['fantasia'].get('cpus_per_task', config['slurm_args']['cpus_per_task']))
    resources:
        # GPU resource hints. These only matter when running under --executor slurm;
        # local runs ignore them. The defaults fall back to the regular SLURM_ARGS
        # if no GPU section is configured, so the rule still validates on local runs.
        mem_mb=int(config['fantasia'].get('mem_mb', config['slurm_args']['mem_of_node'])),
        runtime=int(config['fantasia'].get('max_runtime', config['slurm_args']['max_runtime'])),
        slurm_partition=config['fantasia'].get('partition', ''),
        gres="gpu:" + str(config['fantasia'].get('gpus', 1))
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        OUTDIR=$(readlink -f {params.outdir})
        PROTEINS=$(readlink -f {input.proteins})

        # FANTASIA-Lite ships an offline ProtT5 cache; force the HuggingFace
        # libraries to use it and never reach out to the network.
        export HF_HOME="{params.hf_cache}"
        export TRANSFORMERS_OFFLINE=1
        export HF_HUB_OFFLINE=1

        echo "[$(date)] Running FANTASIA-Lite on $PROTEINS" >  {log}
        nProteins=$(grep -c '^>' "$PROTEINS" || echo 0)
        echo "[$(date)] Input proteins: $nProteins" >> {log}

        # Validated invocation copied from EukAssembly-Bin/rules/fantasia.smk.
        # The flags below were debugged extensively on an A100 in the Hoff lab;
        # do not edit casually -- changes here have repeatedly broken the run.
        singularity exec --nv \
            -B "$PWD":"$PWD" \
            -B "{params.hf_cache}":"{params.hf_cache}" \
            "{params.sif}" \
            python3 /opt/fantasia-lite/src/fantasia_pipeline.py \
                --serial-models \
                --embed-models prot_t5 \
                --device cuda \
                --venv-dir /opt/venv \
                --lookup-npz /opt/fantasia-lite/data/lookup/lookup_table.npz \
                --annotations-json /opt/fantasia-lite/data/lookup/annotations.json \
                --accessions-json /opt/fantasia-lite/data/lookup/accessions.json \
                --embeddings-npz "$OUTDIR/query_embeddings.npz" \
                --config-yaml "$OUTDIR/fantasia_config.yaml" \
                --results-csv "$OUTDIR/results.csv" \
                --topgo-dir "$OUTDIR/topgo" \
                --chunk-dir "$OUTDIR/tmp/fasta_chunks" \
                --chunk-embed-dir "$OUTDIR/tmp/chunk_embeddings" \
                --chunk-results-dir "$OUTDIR/tmp/chunk_results" \
                --chunk-config-dir "$OUTDIR/tmp/chunk_configs" \
                --chunk-failure-dir "$OUTDIR/tmp/failures" \
                --failure-report "$OUTDIR/failed_sequences.csv" \
                {params.add_params} \
                "$PROTEINS" \
            >> {log} 2>&1

        echo "[$(date)] FANTASIA-Lite complete" >> {log}
        touch {output.done}

        # Citations
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite fantasia "$REPORT_DIR"
        cite fantasia_methods "$REPORT_DIR"
        """


rule fantasia_decorate_gff3:
    """Add Ontology_term=GO:... attributes to mRNA + gene features in a BRAKER GFF3.

    Wildcard `gff_base` matches `braker` and `braker_with_ncRNA` so the same
    rule handles both the protein-coding-only GFF3 and the ncRNA-merged GFF3
    (the latter is only requested when run_ncrna=1). The decorated copy is
    written alongside the FANTASIA outputs and copied to the top-level results
    directory by collect_results.
    """
    wildcard_constraints:
        gff_base="braker|braker_with_ncRNA"
    input:
        gff3="output/{sample}/{gff_base}.gff3",
        results="output/{sample}/fantasia/results.csv"
    output:
        decorated="output/{sample}/fantasia/{gff_base}.go.gff3"
    log:
        "logs/{sample}/fantasia/fantasia_decorate_{gff_base}.log"
    benchmark:
        "benchmarks/{sample}/fantasia/fantasia_decorate_{gff_base}.txt"
    params:
        min_score=FANTASIA_MIN_SCORE
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1

        python3 {script_dir}/fantasia_decorate_gff3.py \
            --gff3-in   {input.gff3} \
            --gff3-out  {output.decorated} \
            --results   {input.results} \
            --min-score {params.min_score} \
            > {log} 2>&1
        """


rule fantasia_summarize:
    """Parse FANTASIA-Lite results.csv into a summary text and a GO namespace bar plot."""
    input:
        results="output/{sample}/fantasia/results.csv"
    output:
        summary="output/{sample}/fantasia/fantasia_summary.txt",
        plot="output/{sample}/fantasia/fantasia_go_categories.png",
        go_terms="output/{sample}/fantasia/fantasia_go_terms.tsv"
    log:
        "logs/{sample}/fantasia/fantasia_summarize.log"
    benchmark:
        "benchmarks/{sample}/fantasia/fantasia_summarize.txt"
    params:
        outdir=lambda wc: f"output/{wc.sample}/fantasia",
        min_score=FANTASIA_MIN_SCORE
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1

        python3 {script_dir}/fantasia_summary.py \
            --results {input.results} \
            --out-dir {params.outdir} \
            --min-score {params.min_score} \
            > {log} 2>&1
        """
