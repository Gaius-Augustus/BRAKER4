rule prepare_genome:
    """
    Copy genome to output directory and clean FASTA headers.

    FASTA headers are simplified to the first word (accession only).
    This is critical because GeneMark uses the full header as the
    sequence name in GTF output, which must match the FASTA key for
    all downstream tools (getAnnoFastaFromJoingenes.py, AUGUSTUS, etc.).

    Lines starting with ';' (and blank lines) before the first '>' record
    are dropped. These are Pearson FASTA comments, used e.g. by the
    FAIR-bioHeaders metadata block (FHR, ';~' YAML lines,
    https://doi.org/10.1093/bib/bbae122). samtools faidx, GeneMark
    (probuild), ProtHint and splitMfasta.pl reject or mis-parse them
    (Gaius-Augustus/BRAKER#982). Only the leading block is removed, so
    sequences, headers and coordinates are unchanged.

    If headers were modified, a marker file is created so the cleaned
    genome can be included in the final results.
    """
    input:
        genome = lambda wildcards: samples_df[samples_df["sample_name"] == wildcards.sample].iloc[0]["genome"]
    output:
        genome = "output/{sample}/genome.fa",
        genome_fai = "output/{sample}/genome.fa.fai",
        headers_fixed = "output/{sample}/preprocessing/.headers_fixed"
    benchmark:
        "benchmarks/{sample}/prepare_genome/prepare_genome.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        pipeline_version=config['pipeline_version']
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        mkdir -p output/{wildcards.sample}/preprocessing

        # Drop leading ';' comment lines (FHR metadata, BRAKER#982) and
        # clean FASTA headers: keep only first word (accession)
        N_COMMENT=$(awk '/^>/ {{exit}} /^[[:space:]]*;/ {{n++}} END {{print n+0}}' {input.genome})
        if [ "$N_COMMENT" -gt 0 ]; then
            echo "[INFO] Dropped $N_COMMENT leading ';' comment line(s) from genome FASTA (e.g. FAIR-bioHeaders metadata)"
        fi
        awk '!seen && /^[[:space:]]*(;|$)/ {{next}} {{seen=1; if (/^>/) sub(/ .*/, ""); print}}' \
            {input.genome} > {output.genome}
        # Duplicated sequence names pass samtools faidx with only a warning
        # and fail hours later (randomSplit.pl: "LOCUS names in genbank file
        # are not unique", Gaius-Augustus/BRAKER#977). Stop here instead.
        DUPS=$(awk '/^>/ {{ n = substr($1, 2); if (c[n]++ == 1 && k++ < 5) printf "%s ", n }}' {output.genome})
        if [ -n "$DUPS" ]; then
            echo "ERROR: sequence names occur more than once in the genome FASTA (after cutting headers at the first space), e.g.: $DUPS" >&2
            echo "  Rename or remove the duplicated sequences and start again." >&2
            exit 1
        fi
        samtools faidx {output.genome}

        # Check if headers were actually modified
        if diff -q <(grep '^>' {input.genome}) <(grep '^>' {output.genome}) > /dev/null 2>&1; then
            echo "no" > {output.headers_fixed}
        else
            echo "yes" > {output.headers_fixed}
            echo "[INFO] Cleaned FASTA headers in genome (stripped descriptions)"
        fi

        # Record pipeline and Snakemake versions (runs on host, no container)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        mkdir -p $(dirname $VERSIONS_FILE)
        SM_VER=$(snakemake --version 2>/dev/null || python3 -c "import snakemake; print(snakemake.__version__)" 2>/dev/null || echo "unknown")
        ( flock 9
          printf "BRAKER4 pipeline\t%s\n" "{params.pipeline_version}" >> "$VERSIONS_FILE"
          printf "Snakemake\t%s\n" "$SM_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Cite pipeline and workflow manager
        REPORT_DIR=output/{wildcards.sample}
        mkdir -p $REPORT_DIR
        source {script_dir}/report_citations.sh
        cite braker4 "$REPORT_DIR"
        cite snakemake "$REPORT_DIR"
        """


rule prepare_masked_genome:
    """
    Copy pre-masked genome to output directory and clean FASTA headers.

    Only runs when user provides a genome_masked column in samples.csv.
    Leading ';' comment lines are dropped and headers are cleaned
    identically to prepare_genome.
    """
    input:
        genome_masked = lambda wildcards: samples_df[samples_df["sample_name"] == wildcards.sample].iloc[0]["genome_masked"]
    output:
        genome_masked = "output/{sample}/genome_masked.fa"
    benchmark:
        "benchmarks/{sample}/prepare_masked_genome/prepare_masked_genome.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        # Drop leading ';' comment lines and clean FASTA headers
        # (same as prepare_genome, BRAKER#982)
        awk '!seen && /^[[:space:]]*(;|$)/ {{next}} {{seen=1; if (/^>/) sub(/ .*/, ""); print}}' \
            {input.genome_masked} > {output.genome_masked}
        # Duplicated sequence names pass samtools faidx with only a warning
        # and fail hours later (randomSplit.pl: "LOCUS names in genbank file
        # are not unique", Gaius-Augustus/BRAKER#977). Stop here instead.
        DUPS=$(awk '/^>/ {{ n = substr($1, 2); if (c[n]++ == 1 && k++ < 5) printf "%s ", n }}' {output.genome_masked})
        if [ -n "$DUPS" ]; then
            echo "ERROR: sequence names occur more than once in the genome FASTA (after cutting headers at the first space), e.g.: $DUPS" >&2
            echo "  Rename or remove the duplicated sequences and start again." >&2
            exit 1
        fi
        """


rule create_empty_hints:
    """Create empty hints file for ES mode (no evidence)."""
    output:
        "output/{sample}/empty_hints.gff"
    benchmark:
        "benchmarks/{sample}/create_empty_hints/create_empty_hints.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        "touch {output}"
