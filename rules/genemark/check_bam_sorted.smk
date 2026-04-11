"""
Check if BAM file is coordinate-sorted and sort if necessary.

This rule fixes a common BRAKER issue where unsorted BAM files cause
cryptic failures hours into the pipeline. It also parallelizes sorting
(original BRAKER uses single-threaded sorting).

Input:
    - BAM file (may or may not be sorted)

Output:
    - Coordinate-sorted BAM file
    - Index file (.bai)

Container: teambraker/braker3:latest (contains samtools)
"""

def get_input_bam(wildcards):
    """Get the input BAM file path for a given bam_id."""
    bam_files = get_bam_files(wildcards.sample)
    bam_ids = get_bam_ids(wildcards.sample)
    # Find the BAM file corresponding to this bam_id
    for bam_path, bam_id in zip(bam_files, bam_ids):
        if bam_id == wildcards.bam_id:
            return bam_path
    raise ValueError(f"No BAM file found for bam_id {wildcards.bam_id}")

rule check_bam_sorted:
    input:
        bam=get_input_bam
    output:
        bam=temp("output/{sample}/bam_sorted/{bam_id}.sorted.bam"),
        bai=temp("output/{sample}/bam_sorted/{bam_id}.sorted.bam.bai")
    log:
        "logs/{sample}/check_bam_sorted/{bam_id}.log"
    benchmark:
        "benchmarks/{sample}/check_bam_sorted/{bam_id}.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Check if BAM is already coordinate-sorted and indexable
        # Even if header says sorted, unmapped reads might be in wrong position
        IS_SORTED=false

        if samtools view -H {input.bam} | grep -q '@HD.*SO:coordinate'; then
            echo "BAM file {input.bam} header indicates coordinate-sorted" > {log}

            # Try to create a symlink and index it
            ln -sf $(readlink -f {input.bam}) {output.bam} 2>> {log}

            # Test if we can index it
            if samtools index -@ {threads} {output.bam} 2>> {log}; then
                echo "BAM file is properly sorted and indexable" >> {log}
                IS_SORTED=true
            else
                echo "BAM file claims to be sorted but cannot be indexed (unmapped reads in wrong position)" >> {log}
                echo "Will re-sort to fix the issue" >> {log}
                rm -f {output.bam} {output.bai}
                IS_SORTED=false
            fi
        else
            echo "BAM file {input.bam} is not coordinate-sorted" > {log}
            IS_SORTED=false
        fi

        # If not properly sorted, re-sort
        if [ "$IS_SORTED" = "false" ]; then
            echo "Sorting BAM file..." >> {log}

            # Sort with parallel threads (fixes BRAKER issue: was -@ 0)
            samtools sort \
                -@ {threads} \
                -o {output.bam} \
                {input.bam} \
                2>> {log}

            echo "Sorting complete" >> {log}

            # Create index
            samtools index -@ {threads} {output.bam} 2>> {log}
            echo "Indexing complete" >> {log}
        fi

        echo "Final BAM file: {output.bam}" >> {log}
        echo "Final index file: {output.bai}" >> {log}
        """
