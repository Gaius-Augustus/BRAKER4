"""
Convert BAM file to intron hints in GFF format.

Uses AUGUSTUS bam2hints tool to extract splice junction evidence
from RNA-Seq alignments. This creates intron hints that guide
GeneMark-ET and AUGUSTUS gene prediction.

Input:
    - Coordinate-sorted BAM file
    - BAM index (.csi)

Output:
    - GFF file with intron hints

Container: BRAKER3_CONTAINER (teambraker/braker3, tag pinned in the Snakefile; contains augustus/bin/bam2hints)
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
    params:
        # 1: a library without spliced reads is not an error, because the
        # RNA-Seq check of the automatic mode fallback (check_evidence.smk)
        # then sends an ET sample to ES and an ETP sample to EP mode
        allow_empty=lambda wildcards: 1 if (config.get("auto_mode_fallback", True)
                                            and get_braker_mode(wildcards.sample) in ("et", "etp")) else 0
    threads: 1
    resources:
        mem_mb=mem_single_thread,
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
            if [ "{params.allow_empty}" = "1" ]; then
                echo "WARNING: bam2hints found no spliced reads in {input.bam} (empty output); the RNA-Seq check of the automatic mode fallback decides how the sample continues" >> {log}
                : > {output.hints}
            else
                echo "ERROR: bam2hints produced empty output" >> {log}
                exit 1
            fi
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
