"""
Score the final gene set on the held-out AUGUSTUS test set (#54).

split_training_set keeps train.gb.test away from AUGUSTUS training, and
optimize_augustus scores AUGUSTUS ab initio on it. This rule scores the final
braker.gtf on the same loci, so the report can show both side by side.

QC only: the test genes are GeneMark training genes and the final gene set
contains GeneMark transcripts selected by TSEBRA, so the values are biased
upwards and are not comparable to an evaluation against a reference
annotation (see run_gffcompare for that).

Output: accuracy_final_gene_set.txt in the format of
accuracy_after_optimize.txt, parsed by training_summary.py.
"""


rule eval_test_set_accuracy:
    """Nucleotide, exon and gene level accuracy of braker.gtf on train.gb.test."""
    input:
        gtf="output/{sample}/braker.gtf",
        gb_test="output/{sample}/train.gb.test"
    output:
        acc="output/{sample}/accuracy_final_gene_set.txt"
    log:
        "logs/{sample}/eval_test_set_accuracy/eval_test_set_accuracy.log"
    benchmark:
        "benchmarks/{sample}/eval_test_set_accuracy/eval_test_set_accuracy.txt"
    params:
        script=os.path.join(script_dir, "eval_test_set_accuracy.py")
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        python3 {params.script} \
            -t {input.gb_test} \
            -g {input.gtf} \
            -o {output.acc} \
            > {log} 2>&1
        """
