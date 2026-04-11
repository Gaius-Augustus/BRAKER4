"""
Merge multiple protein FASTA files into one.

When multiple protein files are specified (colon-separated in samples.csv),
they must be concatenated before being passed to ProtHint or GeneMark-ETP.

This is a simple concatenation — no deduplication.
"""


rule merge_proteins:
    """Concatenate multiple protein FASTA files into one."""
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
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged})
        cat {input.proteins} > {output.merged}
        n_seqs=$(grep -c '^>' {output.merged})
        echo "Merged $(echo {input.proteins} | wc -w) protein files: $n_seqs sequences" > {log}
        """
