"""
Merge outputs from dual GeneMark-ETP runs (short-read + IsoSeq).

In dual mode, GeneMark-ETP is run twice:
1. Short-read RNA-Seq + proteins (braker3:latest) → GeneMark-ETP/
2. IsoSeq + proteins (braker3:isoseq) → GeneMark-ETP-isoseq/

This rule:
- Runs TSEBRA to remove redundancy between the two training gene sets
  (--keep_gtf forces both, so only true redundancies are removed)
- Merges HC genes from both runs
- Merges hints from both runs (with multiplicity joining)
- Merges StringTie assemblies from both runs
"""


rule merge_dual_etp:
    """Merge and de-duplicate outputs from both ETP runs using TSEBRA."""
    input:
        training_sr="output/{sample}/GeneMark-ETP/training.gtf",
        training_iso="output/{sample}/GeneMark-ETP-isoseq/training.gtf",
        hc_sr="output/{sample}/GeneMark-ETP/hc.gff",
        hc_iso="output/{sample}/GeneMark-ETP-isoseq/hc.gff",
        hints_sr="output/{sample}/etp_hints.gff",
        hints_iso="output/{sample}/etp_hints_isoseq.gff",
        stringtie_sr="output/{sample}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff",
        stringtie_iso="output/{sample}/GeneMark-ETP-isoseq/rnaseq/stringtie/transcripts_merged.gff"
    output:
        training="output/{sample}/dual_etp_merged/training.gtf",
        hc_gff="output/{sample}/dual_etp_merged/hc.gff",
        etp_hints="output/{sample}/dual_etp_merged/etp_hints.gff",
        stringtie="output/{sample}/dual_etp_merged/transcripts_merged.gff"
    log:
        "logs/{sample}/merge_dual_etp/merge.log"
    benchmark:
        "benchmarks/{sample}/merge_dual_etp/merge.txt"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        mkdir -p $(dirname {output.training})

        echo "Merging dual ETP outputs..." > {log}

        # Step 1: Remove redundancy from training genes using TSEBRA
        # --keep_gtf forces both gene sets, so TSEBRA only removes true redundancies
        # (overlapping genes at the same locus from both runs)
        echo "  Running TSEBRA to remove redundant training genes..." >> {log}

        # Merge hints first (needed by TSEBRA for conflict resolution).
        # NOTE: do NOT call join_mult_hints.pl here. The downstream
        # merge_hints rule does the join correctly with the src=C grp=
        # split that braker.pl requires (see run_genemark_etp.smk for
        # the rationale). Premature collapsing destroys the per-gene
        # grp= linkage and silently degrades AUGUSTUS hint scoring.
        cat {input.hints_sr} {input.hints_iso} > {output.etp_hints}

        TSEBRA_CFG=""
        if [ -f "/opt/TSEBRA/config/default.cfg" ]; then
            TSEBRA_CFG="--cfg /opt/TSEBRA/config/default.cfg"
        fi

        # Prefix IDs to avoid MSTRG.* collisions between the two ETP runs
        TRAINING_SR_PREFIXED=$(dirname {output.training})/training_sr_prefixed.gtf
        TRAINING_ISO_PREFIXED=$(dirname {output.training})/training_iso_prefixed.gtf
        sed 's/gene_id "\([^"]*\)"/gene_id "sr_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "sr_\1"/g' {input.training_sr} > $TRAINING_SR_PREFIXED
        sed 's/gene_id "\([^"]*\)"/gene_id "iso_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "iso_\1"/g' {input.training_iso} > $TRAINING_ISO_PREFIXED

        tsebra.py \
            --keep_gtf $TRAINING_SR_PREFIXED,$TRAINING_ISO_PREFIXED \
            --hintfiles {output.etp_hints} \
            $TSEBRA_CFG \
            --out {output.training} \
            >> {log} 2>&1

        rm -f $TRAINING_SR_PREFIXED $TRAINING_ISO_PREFIXED

        n_training=$(awk '$3=="gene"{{n++}}END{{print n+0}}' {output.training})
        echo "  Training genes after redundancy removal: $n_training" >> {log}

        # Step 2: Merge HC genes (prefix IDs to avoid collisions between runs)
        # Both runs use MSTRG.* transcript IDs from StringTie, so we must
        # make them unique before concatenation.
        sed 's/gene_id "\([^"]*\)"/gene_id "sr_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "sr_\1"/g' {input.hc_sr} > {output.hc_gff}
        sed 's/gene_id "\([^"]*\)"/gene_id "iso_\1"/g; s/transcript_id "\([^"]*\)"/transcript_id "iso_\1"/g' {input.hc_iso} >> {output.hc_gff}
        n_hc=$(wc -l < {output.hc_gff})
        echo "  HC gene features: $n_hc (merged with prefixed IDs)" >> {log}

        # Step 3: Hints already merged above
        n_hints=$(wc -l < {output.etp_hints})
        echo "  Hints: $n_hints (merged + joined)" >> {log}

        # Step 4: Merge StringTie assemblies
        cat {input.stringtie_sr} {input.stringtie_iso} > {output.stringtie}
        n_tx=$(awk '$3=="transcript"{{n++}}END{{print n+0}}' {output.stringtie})
        echo "  StringTie transcripts: $n_tx (merged)" >> {log}

        echo "Dual ETP merge complete" >> {log}
        """
