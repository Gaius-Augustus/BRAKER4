"""
Run GeneMark-ETP for IsoSeq evidence in dual mode.

In dual mode (short-read RNA-Seq + IsoSeq + proteins), GeneMark-ETP is run
twice: once with short-read BAMs (braker3:latest) and once with IsoSeq BAMs
(braker3:isoseq). This rule handles the IsoSeq run.

Output goes to GeneMark-ETP-isoseq/ to avoid conflict with the short-read run.

Container: teambraker/braker3:isoseq (GeneMark-ETP build for long-read evidence)
"""


def _get_etp_isoseq_bam_files(wildcards):
    """Get sorted IsoSeq BAM file for the IsoSeq ETP run."""
    sample = wildcards.sample
    isoseq_bam = get_isoseq_bam_for_etp(sample)
    if isoseq_bam:
        return [isoseq_bam]
    return []

def _get_etp_bam_files(wildcards):
    """Get sorted short-read BAM files for GeneMark-ETP.
    """
    sample = wildcards.sample
    mode = get_braker_mode(sample)
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


rule run_genemark_etp_isoseq:
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
        proteins=lambda wildcards: get_protein_fasta(wildcards.sample),
        bams=_get_etp_isoseq_bam_files,
        sr_bams=_get_etp_bam_files
    output:
        gtf="output/{sample}/GeneMark-ETP-isoseq/genemark.gtf",
        training="output/{sample}/GeneMark-ETP-isoseq/training.gtf",
        hc_gff="output/{sample}/GeneMark-ETP-isoseq/hc.gff",
        etp_hints="output/{sample}/etp_hints_isoseq.gff",
        stringtie_gff="output/{sample}/GeneMark-ETP-isoseq/rnaseq/stringtie/transcripts_merged.gff"
    log:
        "logs/{sample}/genemark_etp_isoseq/genemark_etp_isoseq.log"
    benchmark:
        "benchmarks/{sample}/genemark_etp_isoseq/genemark_etp_isoseq.txt"
    threads: workflow.cores
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        outdir="output/{sample}/GeneMark-ETP-isoseq",
        species_name=lambda wildcards: get_species_name(wildcards),
        fungus="--fungus" if config.get("fungus", False) else ""
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Disable set -e and pipefail: gmetp.pl, get_etp_hints.py, and
        # join_mult_hints.pl may return non-zero. See README.md Developer Notes.
        set +e
        set +o pipefail
        WORKDIR=$(pwd)
        mkdir -p {params.outdir}/etp_lr_data # isoseq reads
        mkdir -p {params.outdir}/etp_sr_data # short reads

        OUTDIR_ABS=$(readlink -f {params.outdir})
        GENOME_ABS=$(readlink -f {input.genome})
        PROTEINS_ABS=$(readlink -f {input.proteins})

        # Step 1: Copy IsoSeq BAM into etp_lr_data/
        echo "Preparing IsoSeq BAM for GeneMark-ETP (isoseq)..." > {log}
        BAM_IDS=""
        for bam in {input.bams}; do
            BAM_ABS=$(readlink -f $bam)
            LIB=$(basename $bam .sorted.bam)
            cp $BAM_ABS $OUTDIR_ABS/etp_lr_data/${{LIB}}.bam
            if [ -z "$BAM_IDS" ]; then
                BAM_IDS="$LIB"
            else
                BAM_IDS="$BAM_IDS,$LIB"
            fi
            echo "  Prepared BAM: $LIB" >> {log}
        done

        echo "Preparing SR BAM for GeneMark-ETP (isoseq)..." > {log}
        BAM_IDS=""
        for bam in {input.sr_bams}; do
            BAM_ABS=$(readlink -f $bam)
            LIB=$(basename $bam .sorted.bam)
            cp $BAM_ABS $OUTDIR_ABS/etp_sr_data/${{LIB}}.bam
            if [ -z "$BAM_IDS" ]; then
                BAM_IDS="$LIB"
            else
                BAM_IDS="$BAM_IDS,$LIB"
            fi
            echo "  Prepared BAM: $LIB" >> {log}
        done

        # Step 2: Prepare protein file
        PROT_FILE=$WORKDIR/output/{wildcards.sample}/proteins_isoseq.fa
        sed '/^>/!s/\\.$//' $PROTEINS_ABS > $PROT_FILE

        # Step 3: Create YAML config
        cat > $OUTDIR_ABS/etp_config.yaml << YAMLEOF
---
RepeatMasker_path: ''
annot_path: ''
genome_path: $GENOME_ABS
protdb_path: $(readlink -f $PROT_FILE)
rnaseq_sets: [$BAM_IDS]
species: {params.species_name}_isoseq
YAMLEOF

        echo "YAML config created with rnaseq_sets: [$BAM_IDS]" >> {log}

        GMES_CORES={threads}
        # Step 4: Run GeneMark-ETP with isoseq container
        cd $OUTDIR_ABS

        if gmetp.pl \
            --cfg $OUTDIR_ABS/etp_config.yaml \
            --workdir $OUTDIR_ABS \
            --long_bam $OUTDIR_ABS/etp_lr_data/RNAseq.bam \
            --bam $OUTDIR_ABS/etp_sr_data/ \
            --cores $GMES_CORES \
            --softmask \
            {params.fungus} \
            >> $WORKDIR/{log} 2>&1
        then
            ETP_EXIT=0
        else
            ETP_EXIT=$?
        fi

        cd $WORKDIR

        if [ ! -f $OUTDIR_ABS/genemark.gtf ]; then
            echo "ERROR: GeneMark-ETP (isoseq) failed (exit=$ETP_EXIT)" >> {log}
            if grep -q "Illegal division by zero" $WORKDIR/{log} 2>/dev/null; then
                N_TRAIN=$(grep "genes found for training:" $WORKDIR/{log} | tail -1 | awk '{{print $NF}}' 2>/dev/null || echo "unknown")
                echo "HINT: GeneMark-ETP crashed in model training (parse_set.pl / train_super.pl division by zero)." >> {log}
                echo "  Training genes found: $N_TRAIN (minimum needed: ~100 multi-exon genes)." >> {log}
                echo "  Likely causes:" >> {log}
                echo "    1. Fungal organism: add 'fungus: true' to config.yaml" >> {log}
                echo "    2. Low IsoSeq coverage: too few reads to build reliable gene models" >> {log}
                echo "    3. Protein database too distant: ProtHint yields too few HC introns" >> {log}
                echo "  See https://github.com/gatech-genemark/GeneMark-ETP/issues" >> {log}
            else
                TSEQ=$OUTDIR_ABS/rnaseq/stringtie/transcripts_merged.fasta
                if [ -f "$TSEQ" ]; then
                    TSEQ_SIZE=$(wc -c < "$TSEQ")
                    echo "DIAGNOSTIC: transcripts_merged.fasta size: $TSEQ_SIZE bytes" >> {log}
                    if [ "$TSEQ_SIZE" -eq 0 ]; then
                        echo "HINT: transcripts_merged.fasta is empty -- StringTie produced no transcripts from the IsoSeq BAM. Check alignment quality and coverage." >> {log}
                    else
                        echo "HINT: exit 139 = segfault in gmhmmp (internal GeneMark binary). This is a known issue when the transcript set is very large. Try subsampling the IsoSeq BAM and rerunning." >> {log}
                    fi
                else
                    echo "DIAGNOSTIC: transcripts_merged.fasta not found -- GeneMark-ETP likely crashed before StringTie completed." >> {log}
                fi
            fi
            GMS_LOG=$(find $OUTDIR_ABS -name "gms.log" | head -1)
            if [ -n "$GMS_LOG" ]; then
                echo "DIAGNOSTIC: last lines of gms.log ($GMS_LOG):" >> {log}
                tail -10 "$GMS_LOG" >> {log}
            fi
            exit 1
        fi

        n_genes=$(awk '$3=="gene"{{c++}}END{{print c+0}}' $OUTDIR_ABS/genemark.gtf)
        echo "GeneMark-ETP (isoseq) predicted $n_genes genes (exit=$ETP_EXIT)" >> {log}

        # Step 5: Find and copy training genes and HC genes
        ETP_MODEL=$(find $OUTDIR_ABS -path "*/model/training.gtf" -not -path "*/etr/*" | head -1 | xargs dirname 2>/dev/null || echo "")

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/training.gtf" ]; then
            cp "$ETP_MODEL/training.gtf" {output.training}
        else
            echo "WARNING: No model/training.gtf found, using genemark.gtf" >> {log}
            cp $OUTDIR_ABS/genemark.gtf {output.training}
        fi

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/hc.gff" ]; then
            cp "$ETP_MODEL/hc.gff" {output.hc_gff}
        else
            touch {output.hc_gff}
        fi

        # Step 6: Extract hints
        # CRITICAL: --genemark_scripts must point at /opt/ETP/bin (where
        # format_back.pl lives), NOT /opt/ETP/bin/gmes/. See
        # run_genemark_etp.smk for the full explanation. Also note that
        # get_etp_hints.py uses >> (append) for output, so the file
        # must be truncated first to make the rule re-runnable.
        # get_etp_hints.py probes for proteins.fa in --etp_wdir to detect a
        # valid GeneMark-ETP run. Our isoseq variant writes proteins_isoseq.fa,
        # so symlink the expected name. -f makes the rule re-runnable.
        ln -sf $OUTDIR_ABS/proteins_isoseq.fa $OUTDIR_ABS/proteins.fa
        ln -sf $OUTDIR_ABS/rnaseq/hints/proteins_isoseq.fa $OUTDIR_ABS/rnaseq/hints/proteins.fa
        rm -f {output.etp_hints}
        if get_etp_hints.py \
            --genemark_scripts /opt/ETP/bin \
            --out {output.etp_hints} \
            --etp_wdir $OUTDIR_ABS \
            >> {log} 2>&1
        then
            HINTS_EXIT=0
        else
            HINTS_EXIT=$?
        fi

        # Hard-fail if get_etp_hints.py didn't produce a file. The manual
        # fallback that lived here historically was structurally wrong
        # (raw nonhc coordinates and only one hintsfile_merged copy
        # instead of two).
        if [ ! -s {output.etp_hints} ]; then
            echo "[ERROR] get_etp_hints.py failed (exit=$HINTS_EXIT) and produced no hints." >> {log}
            echo "[ERROR] Inspect $OUTDIR_ABS for the GeneMark-ETP run state." >> {log}
            exit 1
        fi

        # NOTE: do NOT call join_mult_hints.pl here. See run_genemark_etp.smk
        # for the rationale. The downstream merge_hints rule does the join
        # correctly with the src=C grp= split that braker.pl requires.

        echo "GeneMark-ETP (isoseq) completed" >> {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        GMETP_VER=$(grep -oP 'my \$version = "\K[^"]+' $(which gmetp.pl) 2>/dev/null || true)
        GM_COMMIT=$(grep 'refs/remotes/origin/main' /opt/ETP/.git/packed-refs 2>/dev/null | awk '{{print substr($1,1,7)}}' || true)
        ( flock 9; printf "GeneMark-ETP (IsoSeq)\t%s (commit %s)\n" "$GMETP_VER" "$GM_COMMIT" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh || true
        cite genemark_etp "$REPORT_DIR" || true
        cite genemarks_t "$REPORT_DIR" || true
        cite braker3 "$REPORT_DIR" || true
        cite braker_book "$REPORT_DIR" || true
        """
