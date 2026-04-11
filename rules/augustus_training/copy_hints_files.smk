
rule copy_bam2hints:
    """
    Copy BAM2hints file from failed BRAKER3 run.

    This file contains hints derived from RNA-Seq BAM files via bam2hints.
    We can reuse this directly from the GeneMark-ETP output directory.

    This is a local rule as it's just copying an existing file.
    """
    input:
        bam2hints = lambda w: os.path.join(get_braker_dir(w), "GeneMark-ETP/rnaseq/hints/hintsfile_merged.gff")
    output:
        bam2hints = "output/{sample}/bam2hints.gff"
    benchmark:
        "benchmarks/{sample}/copy_bam2hints/copy_bam2hints.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        cp {input.bam2hints} {output.bam2hints}
        echo "[INFO] Copied BAM2hints file: {output.bam2hints}"
        """

rule copy_prothint_hints:
    """
    Copy ProtHint hints file from failed BRAKER3 run.

    This file contains hints P (protein alignment) and C (chained alignment)
    derived from protein evidence via ProtHint.
    We can reuse this directly from the GeneMark-ETP output directory.

    This is a local rule as it's just copying an existing file.
    """
    input:
        prothint = lambda w: os.path.join(get_braker_dir(w), "GeneMark-ETP/rnaseq/hints/proteins.fa/prothint/prothint_augustus.gff")
    output:
        prothint = "output/{sample}/prothint_hints.gff"
    benchmark:
        "benchmarks/{sample}/copy_prothint_hints/copy_prothint_hints.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        cp {input.prothint} {output.prothint}
        echo "[INFO] Copied ProtHint file: {output.prothint}"
        """
