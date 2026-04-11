"""
Join and merge multiple hint files into a single hints file.

Combines hints from multiple BAM files, filters and sorts them,
then joins identical hints into multiplicity hints using
AUGUSTUS join_mult_hints.pl script.

Input:
    - Multiple hint GFF files (from bam2hints)

Output:
    - Single merged hintsfile.gff with multiplicities

Container: teambraker/braker3:latest (contains join_mult_hints.pl)
"""

rule join_hints:
    input:
        hints=lambda wildcards: expand(
            "output/{sample}/hints/individual/{bam_id}.hints.gff",
            sample=wildcards.sample,
            bam_id=get_all_rnaseq_ids(wildcards.sample)
        )
    output:
        hintsfile="output/{sample}/bam2hints.raw.gff"
    log:
        "logs/{sample}/join_hints/join_hints.log"
    benchmark:
        "benchmarks/{sample}/join_hints/join_hints.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Concatenate all individual hint files
        cat {input.hints} > output/{wildcards.sample}/hints/hintsfile.temp.gff

        # Separate hints: src=C with grp tags should not be merged
        grep 'src=C' output/{wildcards.sample}/hints/hintsfile.temp.gff | \
            grep -E 'grp=|group=' > output/{wildcards.sample}/hints/no_merge.gff || true

        # Everything else can be merged
        grep -v 'src=C' output/{wildcards.sample}/hints/hintsfile.temp.gff > output/{wildcards.sample}/hints/to_merge.gff || true
        grep 'src=C' output/{wildcards.sample}/hints/hintsfile.temp.gff | \
            grep -vE 'grp=|group=' >> output/{wildcards.sample}/hints/to_merge.gff || true

        # Sort hints for join_mult_hints.pl
        # Sort by: scaffold, type, start, end
        if [ -s output/{wildcards.sample}/hints/to_merge.gff ]; then
            cat output/{wildcards.sample}/hints/to_merge.gff | \
                sort -k1,1 -k3,3 -k4,4n -k5,5n \
                > output/{wildcards.sample}/hints/sorted.gff

            # Join identical hints into multiplicity hints
            join_mult_hints.pl \
                < output/{wildcards.sample}/hints/sorted.gff \
                > output/{wildcards.sample}/hints/joined.gff \
                2>> {log}

            # Combine joined hints with non-merged hints
            cat output/{wildcards.sample}/hints/joined.gff \
                output/{wildcards.sample}/hints/no_merge.gff \
                > {output.hintsfile}
        else
            # If nothing to merge, just use non-merged hints
            cp output/{wildcards.sample}/hints/no_merge.gff {output.hintsfile}
        fi

        # Cleanup temp files
        rm -f output/{wildcards.sample}/hints/hintsfile.temp.gff \
              output/{wildcards.sample}/hints/to_merge.gff \
              output/{wildcards.sample}/hints/no_merge.gff \
              output/{wildcards.sample}/hints/sorted.gff \
              output/{wildcards.sample}/hints/joined.gff

        echo "Final hints file contains $(wc -l < {output.hintsfile}) hints" >> {log}
        """
