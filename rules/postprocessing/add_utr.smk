"""
Add UTR features to BRAKER gene predictions using StringTie evidence.

Two-step process:
1. StringTie transcript assembly:
   - ET mode: run StringTie on RNA-Seq BAMs (separate rule)
   - ETP/IsoSeq mode: use StringTie assembly from GeneMark-ETP (already exists)
   - dual mode: additionally a long-read-only assembly of the IsoSeq BAM
     (stringtie -L); its UTRs take priority (#29)
2. stringtie2utr.py: decorate BRAKER CDS predictions with UTRs from StringTie;
   UTRs stop at the nearest same-strand neighbour gene that has a StringTie
   match of its own (#29). check_gtf_loci.py then fails the rule if a model
   spans several loci (#62)

Container: teambraker/braker3:latest (contains stringtie, python3 + intervaltree)
"""


def _get_stringtie_bam(wildcards):
    """Get BAM file(s) for StringTie assembly (ET mode only).

    In ETP/IsoSeq mode, StringTie is run internally by GeneMark-ETP,
    so this rule only handles ET mode (short-read RNA-Seq without proteins).
    Collects BAMs from all short-read RNA-Seq sources: BAM, SRA, FASTQ, VARUS.
    """
    sample = wildcards.sample
    bams = []

    for bid in get_bam_ids(sample):
        bams.append(f"output/{sample}/bam_sorted/{bid}.sorted.bam")
    for sid in get_sra_ids(sample):
        bams.append(f"output/{sample}/hisat2_aligned/{sid}.sorted.bam")
    for fid in get_fastq_ids(sample):
        bams.append(f"output/{sample}/hisat2_aligned/{fid}.sorted.bam")
    for vid in get_varus_ids(sample):
        bams.append(f"output/{sample}/varus/{vid}.sorted.bam")

    return bams


rule run_stringtie:
    """Run StringTie to assemble transcripts from RNA-Seq BAMs (ET mode only)."""
    input:
        bams=_get_stringtie_bam
    output:
        assembly="output/{sample}/stringtie/stringtie.gtf"
    log:
        "logs/{sample}/stringtie/stringtie.log"
    benchmark:
        "benchmarks/{sample}/stringtie/stringtie.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        mkdir -p $(dirname {output.assembly})

        if [ $(echo {input.bams} | wc -w) -gt 1 ]; then
            echo "Merging $(echo {input.bams} | wc -w) BAM files..." > {log}
            samtools merge -@ {threads} output/{wildcards.sample}/stringtie/merged.bam {input.bams} 2>> {log}
            samtools sort -@ {threads} -T {resources.tmpdir}/{wildcards.sample}_addutr -o output/{wildcards.sample}/stringtie/merged.s.bam output/{wildcards.sample}/stringtie/merged.bam 2>> {log}
            samtools index -c -@ {threads} output/{wildcards.sample}/stringtie/merged.s.bam 2>> {log}
            INPUT_BAM=output/{wildcards.sample}/stringtie/merged.s.bam
            rm -f output/{wildcards.sample}/stringtie/merged.bam
        else
            INPUT_BAM={input.bams}
            echo "Using single BAM: $INPUT_BAM" > {log}
        fi

        echo "Running StringTie..." >> {log}
        stringtie "$INPUT_BAM" \
            -o {output.assembly} \
            -p {threads} \
            2>> {log}

        n_tx=$(awk '$3=="transcript"{{n++}}END{{print n+0}}' {output.assembly})
        echo "StringTie assembled $n_tx transcripts" >> {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        ST_VER=$(stringtie --version 2>&1 | head -1 || echo "unknown")
        ( flock 9; printf "StringTie\t%s\n" "$ST_VER" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        rm -f output/{wildcards.sample}/stringtie/merged.s.bam output/{wildcards.sample}/stringtie/merged.s.bam.csi
        """


rule run_stringtie_isoseq:
    """Assemble IsoSeq reads alone (dual mode only).

    GeneMark-ETP assembles short and long reads together (stringtie --mix),
    which keeps short-read read-through joins of neighbouring genes. The
    long-read-only assembly gives add_utr a second source that wins where
    both match a gene (#29).
    """
    input:
        bam=lambda w: get_isoseq_bam_for_etp(w.sample)
    output:
        assembly="output/{sample}/stringtie/isoseq.gtf"
    log:
        "logs/{sample}/stringtie/stringtie_isoseq.log"
    benchmark:
        "benchmarks/{sample}/stringtie/stringtie_isoseq.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        mkdir -p $(dirname {output.assembly})

        echo "Running StringTie in long-read mode on {input.bam}..." > {log}
        stringtie -L "{input.bam}" \
            -o {output.assembly} \
            -p {threads} \
            2>> {log}

        n_tx=$(awk '$3=="transcript"{{n++}}END{{print n+0}}' {output.assembly})
        echo "StringTie assembled $n_tx transcripts from IsoSeq reads" >> {log}
        """


def _get_isoseq_stringtie_gtf(wildcards):
    """Long-read-only assembly for add_utr, dual mode only (#29)."""
    if get_braker_mode(wildcards.sample) == 'dual':
        return f"output/{wildcards.sample}/stringtie/isoseq.gtf"
    return []


def _get_stringtie_gtf(wildcards):
    """Get StringTie assembly GTF, routing by mode.

    - ET mode: our run_stringtie rule output
    - ETP/IsoSeq mode: GeneMark-ETP's internal StringTie output
    - EP/ES mode: should not be called (no transcript evidence)
    """
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'dual':
        return f"output/{sample}/GeneMark-ETP-isoseq/rnaseq/stringtie/transcripts_merged.gff"
    elif mode in ('etp', 'isoseq'):
        return f"output/{sample}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff"
    else:  # et
        return f"output/{sample}/stringtie/stringtie.gtf"



rule add_utr:
    """Decorate BRAKER CDS predictions with UTRs from StringTie assembly.

    Local rule (no container) — uses host Python3 with intervaltree package.
    Install with: pip install intervaltree
    """
    input:
        genes="output/{sample}/braker.gtf",
        stringtie=_get_stringtie_gtf,
        isoseq_stringtie=_get_isoseq_stringtie_gtf
    output:
        utr_gtf="output/{sample}/braker_utr.gtf"
    log:
        "logs/{sample}/add_utr/add_utr.log"
    benchmark:
        "benchmarks/{sample}/add_utr/add_utr.txt"
    params:
        script=os.path.join(script_dir, "stringtie2utr.py"),
        check=os.path.join(script_dir, "check_gtf_loci.py")
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        echo "Decorating gene predictions with UTRs from StringTie..." > {log}

        # Dual mode: UTRs from the IsoSeq-only assembly take priority (#29)
        LR_OPT=""
        if [ -n "{input.isoseq_stringtie}" ]; then
            LR_OPT="-l {input.isoseq_stringtie}"
        fi

        python3 {params.script} \
            -g {input.genes} \
            -s {input.stringtie} \
            $LR_OPT \
            -o {output.utr_gtf} \
            2>> {log}

        # Every gene and transcript on one sequence and strand, spans consistent (#62)
        python3 {params.check} {output.utr_gtf} 2>> {log}

        n_genes=$(awk '$3=="gene"{{n++}}END{{print n+0}}' {output.utr_gtf})
        n_utr=$(awk '$3=="five_prime_UTR"||$3=="three_prime_UTR"{{n++}}END{{print n+0}}' {output.utr_gtf})
        echo "Output: $n_genes genes, $n_utr UTR features" >> {log}

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite stringtie "$REPORT_DIR"
        """
