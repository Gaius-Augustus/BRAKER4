"""
Run GeneMark-ETP for dual mode (short-read + IsoSeq combined).

In dual mode (short-read RNA-Seq + IsoSeq + proteins), GeneMark-ETP is run
ONCE with both short-read BAMs (via --bam) and the IsoSeq BAM (via --long_bam).
This single combined run outputs to GeneMark-ETP-isoseq/.

As in run_genemark_etp, gmetp.pl is run from a shadow copy of the
container's ProtHint bin tree (scripts/make_prothint_shadow.sh) with
the Spaln dispatcher replaced by scripts/spaln_dispatcher.py, to avoid
the Perl-threads dispatcher hang (issue #98).

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
    """GeneMark-ETP on IsoSeq plus short-read BAMs.

    Scratch: the whole gmetp.pl work dir (BAM copies, StringTie per
    library, DIAMOND, Spaln, model dirs) and the ProtHint shadow tree are
    on the node-local disk (scripts/tmp_dir.sh, [paths] tmp_dir).
    get_etp_hints.py reads the work dir there; genemark.gtf and
    rnaseq/stringtie/transcripts_merged.gff are copied back to
    output/<sample>/GeneMark-ETP-isoseq/, training.gtf and hc.gff are copied to
    their output paths. On failure gms.log, etp_config.yaml and the
    StringTie GFFs go to output/<sample>/GeneMark-ETP-isoseq/failed_run_debug/.
    NEED 1 x BAMs + 10 x genome + 5 x proteins (+ 5 GB each); with less
    free the job works in output/<sample>/GeneMark-ETP-isoseq/ as before.
    """
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
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        outdir="output/{sample}/GeneMark-ETP-isoseq",
        species_name=lambda wildcards: get_species_name(wildcards),
        fungus="--fungus" if config.get("fungus", False) else "",
        translation_table=config.get("translation_table", 1),
        tmp_root=TMP_ROOT
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        # Disable set -e and pipefail: gmetp.pl, get_etp_hints.py, and
        # join_mult_hints.pl may return non-zero. See README.md Developer Notes.
        set +e
        set +o pipefail
        WORKDIR=$(pwd)
        source {script_dir}/tmp_dir.sh
        mkdir -p {params.outdir}

        finalDir=$(readlink -f {params.outdir})
        rm -rf "$finalDir/failed_run_debug"
        GENOME_ABS=$(readlink -f {input.genome})
        PROTEINS_ABS=$(readlink -f {input.proteins})
        LOG_ABS=$WORKDIR/{log}
        TRAINING_ABS=$WORKDIR/{output.training}
        HC_GFF_ABS=$WORKDIR/{output.hc_gff}
        ETP_HINTS_ABS=$WORKDIR/{output.etp_hints}

        # Step 1: Copy IsoSeq BAM into etp_lr_data/
        echo "Preparing IsoSeq BAM for GeneMark-ETP (isoseq)..." > {log}

        # The gmetp.pl work dir runs on the node-local disk; only the outputs
        # are copied back to $finalDir.
        scratch_dir outDir "gmetp_isoseq_{wildcards.sample}" "{params.tmp_root}" \
            "$(( $(need_gb 1 {input.bams} {input.sr_bams}) + $(need_gb 10 "$GENOME_ABS") + $(need_gb 5 "$PROTEINS_ABS") ))" \
            "$finalDir" 2>> "$LOG_ABS"
        # also the ProtHint shadow tree next to the scratch dir ($outDir.bin)
        trap 'rm -rf -- "$SCRATCH" ${{SCRATCH:+"$SCRATCH.bin"}}' EXIT
        mkdir -p "$outDir/etp_lr_data" # isoseq reads
        mkdir -p "$outDir/etp_sr_data" # short reads

        # The scratch dir is gone when the job ends: keep the logs needed to
        # debug a failed GeneMark-ETP run.
        save_debug() {{
            local _dbg=(etp_config.yaml) f
            for f in $(find "$outDir" -name gms.log -not -path "*/failed_run_debug/*") "$outDir"/rnaseq/stringtie/*.gff; do
                [ -e "$f" ] && _dbg+=("${{f#"$outDir"/}}")
            done
            mkdir -p "$finalDir/failed_run_debug"
            copy_back "$outDir" "$finalDir/failed_run_debug" "${{_dbg[@]}}"
            echo "Debug files (gms.log, etp_config.yaml, StringTie GFFs) copied to $finalDir/failed_run_debug/" >> "$LOG_ABS"
        }}

        BAM_IDS=""
        for bam in {input.bams}; do
            BAM_ABS=$(readlink -f $bam)
            LR_LIB=$(basename $bam .sorted.bam)
            cp $BAM_ABS $outDir/etp_lr_data/${{LR_LIB}}.bam
            if [ -z "$BAM_IDS" ]; then
                BAM_IDS="$LR_LIB"
            else
                BAM_IDS="$BAM_IDS,$LR_LIB"
            fi
            echo "  Prepared BAM: $LR_LIB" >> {log}
        done

        echo "Preparing SR BAM for GeneMark-ETP (isoseq)..." >> {log}
        BAM_IDS=""
        for bam in {input.sr_bams}; do
            BAM_ABS=$(readlink -f $bam)
            LIB=$(basename $bam .sorted.bam)
            cp $BAM_ABS $outDir/etp_sr_data/${{LIB}}.bam
            if [ -z "$BAM_IDS" ]; then
                BAM_IDS="$LIB"
            else
                BAM_IDS="$BAM_IDS,$LIB"
            fi
            echo "  Prepared BAM: $LIB" >> {log}
        done

        # Step 2: Prepare protein file. The copy is as large as the protein
        # database: write it to the node-local work dir, not to the run dir;
        # protdb/ keeps it apart from the proteins_isoseq.fa/ dir gmetp.pl creates.
        mkdir -p "$outDir/protdb"
        PROT_FILE=$outDir/protdb/proteins_isoseq.fa
        sed '/^>/!s/\\.$//' $PROTEINS_ABS > $PROT_FILE

        # GeneMark-ETP only says "error in protein file parsing" on duplicated
        # IDs or unexpected characters; report the offending records (#99).
        if ! python3 {script_dir}/check_protein_fasta.py $PROT_FILE >> {log} 2>&1; then
            echo "ERROR: fix the protein file {input.proteins} and rerun." >> {log}
            exit 1
        fi

        # Step 3: Create YAML config
        cat > $outDir/etp_config.yaml << YAMLEOF
---
RepeatMasker_path: ''
annot_path: ''
genome_path: $GENOME_ABS
protdb_path: $(readlink -f $PROT_FILE)
rnaseq_sets: [$BAM_IDS]
species: {params.species_name}_isoseq
translation_table: {params.translation_table}
gcode: {params.translation_table}
YAMLEOF

        echo "YAML config created with rnaseq_sets: [$BAM_IDS]" >> {log}

        # Spaln dispatcher fix for ProtHint inside gmetp.pl (#98), see
        # run_genemark_etp.smk.
        SHADOW_DIR=$outDir.bin
        if bash {script_dir}/make_prothint_shadow.sh $SHADOW_DIR >> {log} 2>&1 && [ -x $SHADOW_DIR/gmetp.pl ]; then
            GMETP=$SHADOW_DIR/gmetp.pl
            echo "[INFO] ProtHint will use the process-based Spaln dispatcher (#98)" >> {log}
        else
            GMETP=gmetp.pl
            echo "[WARNING] could not install the Spaln dispatcher fix (#98), using the container's ProtHint" >> {log}
        fi

        GMES_CORES={threads}
        # Step 4: Run GeneMark-ETP with isoseq container
        cd $outDir

        if $GMETP \
            --cfg $outDir/etp_config.yaml \
            --workdir $outDir \
            --long_bam $outDir/etp_lr_data/${{LR_LIB}}.bam \
            --bam $outDir/etp_sr_data/ \
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

        if [ ! -f $outDir/genemark.gtf ]; then
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
                TSEQ=$outDir/rnaseq/stringtie/transcripts_merged.fasta
                if [ -f "$TSEQ" ]; then
                    TSEQ_SIZE=$(wc -c < "$TSEQ")
                    echo "DIAGNOSTIC: transcripts_merged.fasta size: $TSEQ_SIZE bytes" >> {log}
                    if [ "$TSEQ_SIZE" -eq 0 ]; then
                        echo "HINT: transcripts_merged.fasta is empty -- StringTie produced no transcripts from the IsoSeq BAM. Check alignment quality and coverage." >> {log}
                    else
                        if [ "$ETP_EXIT" -eq 139 ]; then
                            echo "HINT: exit 139 = segmentation fault inside GeneMark-ETP." >> {log}
                        else
                            echo "HINT: GeneMark-ETP exited with code $ETP_EXIT without writing genemark.gtf." >> {log}
                        fi
                        echo "  Common causes (check these before subsampling):" >> {log}
                        echo "    1. Thread count: GeneMark-ETP ran with --cores {threads}. Very high counts (e.g. 256) are known to crash it; lower slurm_args: cpus_per_task in config.yaml." >> {log}
                        echo "    2. Sequence name mismatch: @SQ names in the BAM header differ from the genome FASTA headers." >> {log}
                        echo "       Check: samtools view -H <bam> | grep ^@SQ  vs  grep '>' <genome> | head" >> {log}
                        echo "    3. Stranded library: GeneMark-ETP assembles transcripts unstranded, which can produce conflicting models from stranded data." >> {log}
                        echo "    4. Very large transcript set: only if 1-3 are ruled out, subsample the IsoSeq BAM and rerun." >> {log}
                        echo "  The failing step is usually visible in gms.log (tail printed below)." >> {log}
                    fi
                else
                    echo "DIAGNOSTIC: transcripts_merged.fasta not found -- GeneMark-ETP likely crashed before StringTie completed." >> {log}
                fi
            fi
            GMS_LOG=$(find $outDir -name "gms.log" | head -1)
            if [ -n "$GMS_LOG" ]; then
                echo "DIAGNOSTIC: last lines of gms.log ($GMS_LOG):" >> {log}
                tail -10 "$GMS_LOG" >> {log}
            fi
            save_debug
            exit 1
        fi

        # genemark.gtf from gmetp.pl has no gene lines: count distinct gene_id
        n_genes=$(awk 'match($0, /gene_id "[^"]+"/) {{ id = substr($0, RSTART, RLENGTH); if (!(id in seen)) {{ seen[id] = 1; n++ }} }} END {{ print n + 0 }}' $outDir/genemark.gtf)
        echo "GeneMark-ETP (isoseq) predicted $n_genes genes (exit=$ETP_EXIT)" >> {log}

        # Step 5: Find and copy training genes and HC genes
        ETP_MODEL=$(find $outDir -path "*/model/training.gtf" -not -path "*/etr/*" | head -1 | xargs dirname 2>/dev/null || echo "")

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/training.gtf" ]; then
            cp "$ETP_MODEL/training.gtf" "$TRAINING_ABS"
        else
            echo "WARNING: No model/training.gtf found, using genemark.gtf" >> {log}
            cp $outDir/genemark.gtf "$TRAINING_ABS"
        fi

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/hc.gff" ]; then
            cp "$ETP_MODEL/hc.gff" "$HC_GFF_ABS"
        else
            touch "$HC_GFF_ABS"
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
        ln -sf $outDir/proteins_isoseq.fa $outDir/proteins.fa
        ln -sf $outDir/rnaseq/hints/proteins_isoseq.fa $outDir/rnaseq/hints/proteins.fa
        rm -f {output.etp_hints}
        if get_etp_hints.py \
            --genemark_scripts /opt/ETP/bin \
            --out "$ETP_HINTS_ABS" \
            --etp_wdir $outDir \
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
            save_debug
            exit 1
        fi

        copy_back "$outDir" "$finalDir" genemark.gtf rnaseq/stringtie/transcripts_merged.gff
        if [ ! -s "$finalDir/genemark.gtf" ] || [ ! -f "$finalDir/rnaseq/stringtie/transcripts_merged.gff" ]; then
            echo "ERROR: genemark.gtf or rnaseq/stringtie/transcripts_merged.gff missing in $finalDir" >> {log}
            exit 1
        fi
        # the protein copy is not needed any more (fallback; on scratch the
        # EXIT trap removes it)
        if [ -z "$SCRATCH" ]; then
            rm -f "$PROT_FILE"
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

        # Remove GeneMark-ETP-isoseq internal working files and etp_data/ BAM copies.
        # Tracked outputs kept: genemark.gtf, training.gtf, hc.gff,
        # rnaseq/stringtie/transcripts_merged.gff (etp_hints_isoseq.gff is outside outdir).
        if [ -z "$SCRATCH" ]; then
            find "$finalDir" -mindepth 1 -type f \
                ! -name 'genemark.gtf' \
                ! -name 'training.gtf' \
                ! -name 'hc.gff' \
                ! -path '*/rnaseq/stringtie/transcripts_merged.gff' \
                -delete 2>/dev/null || true
            find "$finalDir" -mindepth 1 -type d -empty -delete 2>/dev/null || true
        fi
        rm -rf $SHADOW_DIR
        """
