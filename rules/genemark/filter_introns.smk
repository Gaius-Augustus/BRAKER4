"""
Filter intron hints and determine strand from splice site sequences.

Uses AUGUSTUS filterIntronsFindStrand.pl to:
1. Determine strand from splice site dinucleotides (GT-AG, GC-AG, AT-AC)
2. Set score column to multiplicity value
3. Filter out introns without canonical splice sites

This is a critical step that braker.pl performs between bam2hints and
GeneMark training. Without it, intron hints lack strand information
and contain spurious introns that can cause GeneMark to fail.

Input:
    - Raw joined hints from bam2hints (no strand, score=0)
    - Genome FASTA (for splice site lookup)

Output:
    - Filtered hints with strand and multiplicity score

Container: teambraker/braker3:latest (contains filterIntronsFindStrand.pl)
"""

rule filter_introns_find_strand:
    input:
        hints="output/{sample}/bam2hints.raw.gff",
        genome=lambda wildcards: get_masked_genome(wildcards.sample)
    output:
        hints="output/{sample}/bam2hints.gff"
    log:
        "logs/{sample}/filter_introns/filter_introns.log"
    benchmark:
        "benchmarks/{sample}/filter_introns/filter_introns.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        filterIntronsFindStrand.pl \
            {input.genome} \
            {input.hints} \
            --score \
            > {output.hints} \
            2> {log}

        n_before=$(wc -l < {input.hints})
        n_after=$(wc -l < {output.hints})
        echo "Filtered introns: $n_before -> $n_after" >> {log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        """
