"""
OMArk proteome quality assessment.

Runs OMArk to assess the quality and completeness of predicted proteins
using the OMA orthology database. This provides:
- Completeness (single-copy, duplicated, missing orthologs)
- Contamination detection
- Proteome quality metrics

OMArk handles isoforms natively: omamer search runs on ALL proteins,
then OMArk selects the best-mapped isoform per gene using a .splice file
that lists isoforms grouped by gene.

Requires the LUCA.h5 OMAmer database.

Three-step process:
1. omamer search: map ALL proteins to OMA hierarchical orthologous groups
2. generate_splice_file: create isoform grouping from GTF
3. omark: analyze with --isoform_file to pick best isoform per gene

Container: quay.io/biocontainers/omark:0.4.1
"""


OMAMER_DB = config.get("omamer_db", "../../test_data/LUCA.h5")


rule omamer_search:
    """Map ALL predicted proteins to OMA hierarchical orthologous groups."""
    input:
        proteins="output/{sample}/braker.aa",
        db=OMAMER_DB
    output:
        omamer="output/{sample}/omark/proteome.omamer"
    log:
        "logs/{sample}/omark/omamer_search.log"
    benchmark:
        "benchmarks/{sample}/omark/omamer_search.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        OMARK_CONTAINER
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.omamer})

        omamer search \
            --db {input.db} \
            --query {input.proteins} \
            --nthreads {threads} \
            --out {output.omamer} \
            > {log} 2>&1
        """


rule generate_splice_file:
    """Generate isoform grouping file for OMArk from BRAKER GTF.

    Creates a .splice file where each line lists the transcript IDs
    belonging to one gene, separated by semicolons. OMArk uses this
    to select the best-mapped isoform per gene for its metrics.

    Format: tx1;tx2;tx3 (one line per gene)
    """
    input:
        gtf="output/{sample}/braker.gtf"
    output:
        splice="output/{sample}/omark/isoforms.splice"
    benchmark:
        "benchmarks/{sample}/generate_splice_file/generate_splice_file.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.splice})

        # Extract gene_id -> transcript_id mapping from braker.gtf
        # Handles both AUGUSTUS bare format (col9 = "g1.t1") and
        # quoted format (gene_id "g1"; transcript_id "g1.t1";)
        # Uses Python for robust parsing (available in BRAKER3 container)
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        python3 -c "
import sys, re
from collections import defaultdict

genes = defaultdict(set)
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        if len(cols) < 9:
            continue
        feat = cols[2]
        if feat not in ('transcript', 'mRNA'):
            continue
        attr = cols[8]
        # Try quoted format first
        gm = re.search(r'gene_id \"([^\"]+)\"', attr)
        tm = re.search(r'transcript_id \"([^\"]+)\"', attr)
        if gm and tm:
            gene, tx = gm.group(1), tm.group(1)
        else:
            # Bare format: col9 is just the transcript ID
            tx = attr.strip().rstrip(';')
            # Gene = transcript minus last .tN suffix
            gene = re.sub(r'\.[^.]+$', '', tx)
        genes[gene].add(tx)

for g in sorted(genes):
    print(';'.join(sorted(genes[g])))
" {input.gtf} > {output.splice}

        n_genes=$(wc -l < {output.splice})
        echo "Generated splice file: $n_genes genes"
        """


rule run_omark:
    """Run OMArk quality assessment with isoform-aware analysis."""
    input:
        omamer="output/{sample}/omark/proteome.omamer",
        splice="output/{sample}/omark/isoforms.splice",
        db=OMAMER_DB
    output:
        summary="output/{sample}/omark/omark_summary.txt",
        done="output/{sample}/omark/.done"
    log:
        "logs/{sample}/omark/omark.log"
    benchmark:
        "benchmarks/{sample}/omark/omark.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        OMARK_CONTAINER
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(readlink -f output/{wildcards.sample}/omark)

        omark \
            -f {input.omamer} \
            -d {input.db} \
            -i {input.splice} \
            -o "$OUTDIR" \
            > {log} 2>&1

        # Copy the detailed summary to a predictable location
        DETAILED=$(find "$OUTDIR" -name "*_detailed_summary.txt" | head -1)
        if [ -n "$DETAILED" ] && [ -f "$DETAILED" ]; then
            cp "$DETAILED" {output.summary}
        else
            SUM=$(find "$OUTDIR" -name "*.sum" | head -1)
            if [ -n "$SUM" ] && [ -f "$SUM" ]; then
                cp "$SUM" {output.summary}
            else
                echo "OMArk completed but no summary found" > {output.summary}
            fi
        fi

        touch {output.done}

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        OMARK_VER=$(pip show omark 2>/dev/null | grep -oP '^Version: \K.*' || true)
        OMAMER_VER=$(omamer --version 2>&1 | head -1 || true)
        ( flock 9
          printf "OMArk\t%s\n" "$OMARK_VER" >> "$VERSIONS_FILE"
          printf "OMAmer\t%s\n" "$OMAMER_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite omark "$REPORT_DIR"
        cite omamer "$REPORT_DIR"
        """
