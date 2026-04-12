"""
Rescue multi-exon genes where AUGUSTUS and GeneMark agree but TSEBRA dropped.

When both gene finders independently predict the same CDS intron chain and
TSEBRA still filtered it out, the agreement is a strong confidence signal.
Only genes whose every intron is supported by at least one hint are rescued.

Benchmarked on A. thaliana ETP: +29 TP, +4 FP at gene level.  Only adds
new genes (never modifies existing predictions).
"""


def _both_agree_genemark_input(wildcards):
    """Return the GeneMark GTF path for the sample's pipeline mode."""
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'es':
        return f"output/{sample}/GeneMark-ES/genemark.gtf"
    elif mode == 'ep':
        return f"output/{sample}/GeneMark-EP/genemark.gtf"
    elif mode == 'et':
        return f"output/{sample}/genemark/genemark.gtf"
    else:  # etp, isoseq, dual
        return f"output/{sample}/GeneMark-ETP/genemark.gtf"


def _both_agree_hints_input(wildcards):
    """Return the hints file path.  ES mode has no hints, so we return
    an empty placeholder that the script handles gracefully (no rescues)."""
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'es':
        # ES has no hints — rule will find zero intron hints and rescue nothing.
        # Return the GeneMark GTF as a dummy dep so the rule still runs.
        return f"output/{sample}/GeneMark-ES/genemark.gtf"
    return f"output/{sample}/hintsfile.gff"


rule both_agree_rescue:
    input:
        braker_gtf="output/{sample}/braker.tsebra.gtf",
        augustus_gtf="output/{sample}/augustus.hints.fixed.gtf",
        genemark_gtf=_both_agree_genemark_input,
        hints=_both_agree_hints_input,
    output:
        rescued_gtf="output/{sample}/braker.rescued.gtf",
        rescue_log="output/{sample}/both_agree_rescue.log",
    benchmark:
        "benchmarks/{sample}/both_agree_rescue/both_agree_rescue.txt"
    params:
        script=os.path.join(script_dir, "both_agree_rescue.py"),
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        python3 {params.script} \
            --augustus {input.augustus_gtf} \
            --genemark {input.genemark_gtf} \
            --braker {input.braker_gtf} \
            --hints {input.hints} \
            --output {output.rescued_gtf} \
            2>&1 | tee {output.rescue_log}
        """
