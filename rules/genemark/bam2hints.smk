"""
Convert BAM file to intron hints in GFF format.

Uses AUGUSTUS bam2hints tool to extract splice junction evidence
from RNA-Seq alignments. This creates intron hints that guide
GeneMark-ET and AUGUSTUS gene prediction.

Input:
    - Coordinate-sorted BAM file
    - BAM index (.bai)

Output:
    - GFF file with intron hints

Container: teambraker/braker3:latest (contains augustus/bin/bam2hints)
"""

rule bam2hints:
    input:
        unpack(get_rnaseq_bam)
    output:
        hints=temp("output/{sample}/hints/individual/{bam_id}.hints.gff")
    log:
        "logs/{sample}/bam2hints/{bam_id}.log"
    benchmark:
        "benchmarks/{sample}/bam2hints/{bam_id}.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Run bam2hints from AUGUSTUS
        bam2hints \
            --intronsonly \
            --in={input.bam} \
            --out={output.hints} \
            2>&1 | tee {log}

        # Check output was created
        if [ ! -s {output.hints} ]; then
            echo "ERROR: bam2hints produced empty output" >> {log}
            exit 1
        fi

        n_hints=$(grep -c 'intron' {output.hints} || echo "0")
        echo "Successfully extracted $n_hints intron hints" >> {log}

        # Record software version (bam2hints is part of AUGUSTUS)
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        SAM_VER=$(samtools --version 2>&1 | head -1 | awk '{{print $2}}' || true)
        ( flock 9; printf "SAMtools\t%s\n" "$SAM_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite samtools "$REPORT_DIR"
        """
