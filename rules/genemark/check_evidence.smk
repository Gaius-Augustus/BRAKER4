"""
Automatic mode fallback: is the evidence enough for the planned mode?

A sample whose evidence is too thin for the planned GeneMark mode continues
in the next simpler mode instead of failing ([PARAMS] auto_mode_fallback,
default on):

    ETP -> EP -> ES        ET -> ES

Five checkpoints write one decision file each to
output/<sample>/mode_selection/ (kept by collect_results, so that a finished
run stays finished). get_effective_mode() in rules/common.smk reads them:

    check_rnaseq_evidence   planned ET and ETP. Fewer than
                            fallback_min_introns (150) RNA-Seq intron hints
                            with multiplicity >= 10 in bam2hints.gff:
                            ET -> ES, ETP -> EP (GeneMark-ETP is not started).
    check_etp_result        planned ETP that passed the RNA-Seq check.
                            GeneMark-ETP ran, but stopped for lack of
                            training data (run_genemark_etp wrote
                            mode_selection/etp_insufficient.txt): ETP -> EP.
                            The intron count cannot predict this case: the
                            number of training genes is only known after
                            GeneMarkS-T, ProtHint and the filter have run.
    check_et_result         planned ET that passed the RNA-Seq check.
                            GeneMark-ET stopped with an error of its own in
                            training (run_genemark_et wrote
                            mode_selection/et_insufficient.txt, see
                            scripts/genemark_fallback.sh): ET -> ES.
    check_protein_evidence  EP, planned or after a fallback from ETP. Fewer
                            than fallback_min_introns protein intron hints
                            with multiplicity >= 4 in genemark_hints_ep.gff
                            (ProtHint, seeded by GeneMark-ES): EP -> ES; the
                            GeneMark-ES prediction is already there. When
                            DIAMOND finds no input protein in the genome,
                            run_prothint writes empty hints instead of
                            failing, and this check then counts 0 introns.
    check_ep_result         EP that passed the protein check. GeneMark-EP
                            stopped with an error of its own in training
                            (mode_selection/ep_insufficient.txt): EP -> ES.

After a fallback the evidence that was too thin is not used any more (no
RNA-Seq hints, UTRs or lncRNA after ETP -> EP; no protein hints after
EP -> ES), with one exception: hintsfile.gff of an ET sample is the input of
GeneMark-ET and cannot depend on its result, so after ET -> ES the RNA-Seq
intron hints stay in the AUGUSTUS hints file (GeneMark-ES trains AUGUSTUS).

The limits (150 introns, multiplicity 10 for RNA-Seq and 4 for proteins) are
those of check_genemark_hints() in braker.pl, which warns that GeneMark-ET/EP
may crash below them. IsoSeq and dual mode have no fallback.

A library without spliced reads gives an empty bam2hints file; for ET and ETP
samples the bam2hints rule accepts that (the RNA-Seq check then counts what
the other libraries have).

Calling Snakemake again in the same directory (helpers in rules/common.smk):

    - Every rule whose inputs or params depend on the effective mode has the
      decision files as inputs (mode_decision_files, with_mode_decisions).
      Only then does Snakemake evaluate the rule again when a check runs
      again and decides differently. tests/test_mode_fallback_helpers.py
      guards this.
    - Decision files are only valid for the evidence they were made from.
      reset_stale_mode_decisions() compares the sample's row of samples.csv,
      size and modification time of the files named in it and a few settings
      with mode_selection/evidence.fingerprint when the workflow starts; if
      they differ it removes the decision files, the markers and
      results/.done, and the run proceeds like a new one.
    - A check that will run again while its old decision file exists does
      not make the input functions wait (IncompleteCheckpointException):
      they plan with the old decision. And once the decisions of a sample
      were complete in the running process, a decision that changes the
      path does not make them wait either; the next decision file is just
      an input. Reason: Snakemake 8 treats a job whose inputs wait for a
      checkpoint as up to date when its old output exists, and it keeps
      such jobs in its job cache and runs them without inputs when another
      job asks for their output later.

Container: teambraker/braker3 (python3 only)
"""


checkpoint check_rnaseq_evidence:
    """Enough RNA-Seq intron hints for GeneMark-ET / GeneMark-ETP?"""
    input:
        hints="output/{sample}/bam2hints.gff"
    output:
        decision="output/{sample}/mode_selection/rnaseq_evidence.tsv"
    log:
        "logs/{sample}/mode_selection/rnaseq_evidence.log"
    benchmark:
        "benchmarks/{sample}/mode_selection/rnaseq_evidence.txt"
    params:
        script=os.path.join(script_dir, "check_evidence.py"),
        min_introns=config.get("fallback_min_introns", 150)
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        python3 {params.script} introns \
            --hints {input.hints} \
            --evidence RNA-Seq \
            --min-mult 10 \
            --min-introns {params.min_introns} \
            --out {output.decision} > {log} 2>&1
        """


checkpoint check_etp_result:
    """Did GeneMark-ETP find enough training data?"""
    input:
        gtf="output/{sample}/GeneMark-ETP/genemark.gtf"
    output:
        decision="output/{sample}/mode_selection/etp_result.tsv"
    log:
        "logs/{sample}/mode_selection/etp_result.log"
    benchmark:
        "benchmarks/{sample}/mode_selection/etp_result.txt"
    params:
        script=os.path.join(script_dir, "check_evidence.py"),
        marker=lambda wildcards: f"output/{wildcards.sample}/mode_selection/etp_insufficient.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        python3 {params.script} result \
            --tool GeneMark-ETP \
            --gtf {input.gtf} \
            --marker {params.marker} \
            --out {output.decision} > {log} 2>&1
        """


checkpoint check_et_result:
    """Could GeneMark-ET train on the RNA-Seq hints?"""
    input:
        gtf="output/{sample}/genemark/genemark.gtf"
    output:
        decision="output/{sample}/mode_selection/et_result.tsv"
    log:
        "logs/{sample}/mode_selection/et_result.log"
    benchmark:
        "benchmarks/{sample}/mode_selection/et_result.txt"
    params:
        script=os.path.join(script_dir, "check_evidence.py"),
        marker=lambda wildcards: f"output/{wildcards.sample}/mode_selection/et_insufficient.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        python3 {params.script} result \
            --tool GeneMark-ET \
            --gtf {input.gtf} \
            --marker {params.marker} \
            --out {output.decision} > {log} 2>&1
        """


checkpoint check_ep_result:
    """Could GeneMark-EP train on the protein hints?"""
    input:
        gtf="output/{sample}/GeneMark-EP/genemark.gtf"
    output:
        decision="output/{sample}/mode_selection/ep_result.tsv"
    log:
        "logs/{sample}/mode_selection/ep_result.log"
    benchmark:
        "benchmarks/{sample}/mode_selection/ep_result.txt"
    params:
        script=os.path.join(script_dir, "check_evidence.py"),
        marker=lambda wildcards: f"output/{wildcards.sample}/mode_selection/ep_insufficient.txt"
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        python3 {params.script} result \
            --tool GeneMark-EP \
            --gtf {input.gtf} \
            --marker {params.marker} \
            --out {output.decision} > {log} 2>&1
        """


checkpoint check_protein_evidence:
    """Enough protein intron hints for GeneMark-EP?"""
    input:
        hints="output/{sample}/genemark_hints_ep.gff"
    output:
        decision="output/{sample}/mode_selection/protein_evidence.tsv"
    log:
        "logs/{sample}/mode_selection/protein_evidence.log"
    benchmark:
        "benchmarks/{sample}/mode_selection/protein_evidence.txt"
    params:
        script=os.path.join(script_dir, "check_evidence.py"),
        min_introns=config.get("fallback_min_introns", 150)
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        python3 {params.script} introns \
            --hints {input.hints} \
            --evidence protein \
            --min-mult 4 \
            --min-introns {params.min_introns} \
            --out {output.decision} > {log} 2>&1
        """
