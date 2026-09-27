"""
Fix and convert GTF to GFF3 format using AGAT.

First fixes the GTF (normalizes attributes, adds missing Parent relationships),
then converts to GFF3 format.

Container: quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0
"""


def _get_gtf_for_gff3(wildcards):
    """Get the best GTF for GFF3 conversion: UTR-decorated if available."""
    sample = wildcards.sample
    types = detect_data_types(sample)
    has_transcripts = any([types[k] for k in ['has_bam', 'has_fastq', 'has_sra', 'has_varus', 'has_isoseq']])
    if has_transcripts:
        return f"output/{sample}/braker_utr.gtf"
    return f"output/{sample}/braker.gtf"


rule fix_gtf:
    """Fix GTF format issues using AGAT (normalize attributes, add Parent relationships)."""
    input:
        gtf=_get_gtf_for_gff3
    output:
        gtf="output/{sample}/braker.fixed.gtf"
    log:
        "logs/{sample}/agat/fix_gtf.log"
    benchmark:
        "benchmarks/{sample}/agat/fix_gtf.txt"
    threads: 1
    resources:
        mem_mb=max(int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']), 16000),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        AGAT_CONTAINER
    shell:
        r"""
        set -euo pipefail
        agat_convert_sp_gxf2gxf.pl \
            -g {input.gtf} \
            -o {output.gtf} \
            > {log} 2>&1
        """


rule convert_gtf_to_gff3:
    """Convert fixed GTF to GFF3 format using AGAT."""
    input:
        gtf="output/{sample}/braker.fixed.gtf"
    output:
        gff3="output/{sample}/braker.gff3"
    log:
        "logs/{sample}/agat/convert_to_gff3.log"
    benchmark:
        "benchmarks/{sample}/agat/convert_to_gff3.txt"
    threads: 1
    resources:
        mem_mb=max(int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']), 16000),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        AGAT_CONTAINER
    shell:
        r"""
        set -euo pipefail
        agat_convert_sp_gxf2gxf.pl \
            -g {input.gtf} \
            -o {output.gff3} \
            > {log} 2>&1

        # Protein-coding genes as in NCBI/Ensembl GFF3 (Annotrieve): transcripts
        # with CDS children become mRNA (AGAT promotes some but not all), and
        # their genes get gene_biotype=protein_coding unless already set.
        # ncRNA transcripts (no CDS children) are left unchanged.
        # Two passes: pass 1 collects CDS parents and transcript->gene links.
        _TMP="{output.gff3}.tmp"
        awk '
        function attr(s, key,   n, a, i) {{
          n = split(s, a, ";")
          for (i = 1; i <= n; i++) {{
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", a[i])
            if (index(a[i], key "=") == 1) return substr(a[i], length(key) + 2)
          }}
          return ""
        }}
        BEGIN {{ FS = "\t"; OFS = "\t" }}
        NR == FNR {{
          if ($3 == "CDS") {{
            m = split(attr($9, "Parent"), p, ",")
            for (j = 1; j <= m; j++) coding_tx[p[j]] = 1
          }} else if ($3 == "transcript" || $3 == "mRNA") {{
            tx_gene[attr($9, "ID")] = attr($9, "Parent")
          }}
          next
        }}
        !done {{
          for (t in tx_gene) if (t in coding_tx) coding_gene[tx_gene[t]] = 1
          done = 1
        }}
        /^#/ {{ print; next }}
        {{
          if ($3 == "transcript" && (attr($9, "ID") in coding_tx)) $3 = "mRNA"
          else if ($3 == "gene" && (attr($9, "ID") in coding_gene) && attr($9, "gene_biotype") == "")
            {{ sub(/;$/, "", $9); $9 = $9 ";gene_biotype=protein_coding" }}
          print
        }}
        ' "{output.gff3}" "{output.gff3}" > "$_TMP" && mv "$_TMP" "{output.gff3}"

        # Record software version (LC_ALL=C avoids locale warnings in biocontainer)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        AGAT_VER=$(LC_ALL=C agat --version 2>&1 | head -1 || true)
        ( flock 9; printf "AGAT\t%s\n" "$AGAT_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite agat "$REPORT_DIR"
        """
