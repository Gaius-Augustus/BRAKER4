"""
Merge and decompress protein FASTA files.

When multiple protein files are specified (colon-separated in samples.csv),
they must be concatenated before being passed to ProtHint or GeneMark-ETP.
A single gzipped file (.gz, e.g. an OrthoDB partition as downloaded) also
goes through this rule, because ProtHint/DIAMOND and GeneMark-ETP read only
plain FASTA (DIAMOND reports "Inflate error", Gaius-Augustus/BRAKER#789).

gzip -cdf decompresses .gz files and passes plain files through unchanged.
This is a simple concatenation, no deduplication.
"""


rule merge_proteins:
    """Concatenate protein FASTA files into one, decompressing .gz files."""
    input:
        proteins=lambda wildcards: get_protein_fasta_files(wildcards.sample)
    output:
        merged="output/{sample}/preprocessing/proteins_merged.fa"
    log:
        "logs/{sample}/merge_proteins/merge.log"
    benchmark:
        "benchmarks/{sample}/merge_proteins/merge_proteins.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged})
        for f in {input.proteins}; do gzip -cdf -- "$f"; done > {output.merged}
        n_seqs=$(grep -c '^>' {output.merged})
        echo "Merged/decompressed $(echo {input.proteins} | wc -w) protein file(s): $n_seqs sequences" > {log}
        """
