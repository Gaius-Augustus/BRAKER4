"""
Run GeneMark-ETP for combined RNA-Seq + protein evidence training.

GeneMark-ETP integrates both RNA-Seq and protein evidence for gene
prediction training. It runs internally:
- StringTie for transcriptome assembly from BAM files
- DIAMOND/Spaln for protein-to-genome alignment
- Combined evidence-based gene model training

The pipeline:
1. Sort BAM files into etp_data/ directory
2. Create YAML config for gmetp.pl
3. Run gmetp.pl with RNA-Seq BAMs + protein FASTA
4. Extract hints for AUGUSTUS via get_etp_hints.py

Key outputs:
- genemark.gtf: Gene predictions
- proteins.fa/model/training.gtf: HC training genes for AUGUSTUS
- proteins.fa/model/hc.gff: HC gene structures (used by TSEBRA)
- hintsfile: Combined hints for AUGUSTUS (RNA-Seq + protein)

The ProtHint step run internally by gmetp.pl uses the same Spaln
dispatcher shadow tree as run_prothint (scripts/make_prothint_shadow.sh,
scripts/spaln_dispatcher.py) to avoid the Perl-threads dispatcher hang
(issue #98); gmetp.pl is invoked from the shadow tree instead of the
container's own copy. The shadow tree also holds a patched super.pl that
writes exon counts into set.out (#110).

Dual mode (short-read RNA-Seq + IsoSeq + proteins) is not run by this rule:
it is one combined GeneMark-ETP run in run_genemark_etp_isoseq, writing to
GeneMark-ETP-isoseq/.

Mirrors braker.pl's GeneMark_ETP() + get_etp_hints_for_Augustus().

Input:
    - Masked genome FASTA
    - Sorted RNA-Seq BAM file(s)
    - Protein sequences FASTA

Output:
    - genemark.gtf: Gene predictions
    - training.gtf: HC training genes
    - hc.gff: HC gene structures
    - etp_hints.gff: Combined hints for AUGUSTUS

Container: GENEMARK_ETP_CONTAINER (same image as BRAKER3_CONTAINER: teambraker/braker3, tag pinned in the Snakefile; contains gmetp.pl, get_etp_hints.py)
"""


def _get_etp_bam_files(wildcards):
    """Get sorted BAM files for GeneMark-ETP.

    In etp and dual mode, only short-read BAMs are returned (in dual mode the
    IsoSeq BAM is added by run_genemark_etp_isoseq, which runs one combined
    GeneMark-ETP).
    In isoseq mode, only the IsoSeq BAM is returned.
    """
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    bams = []

    # Short-read BAMs (for etp and dual modes)
    if mode in ('etp', 'dual'):
        for bid in get_bam_ids(sample):
            bams.append(f"output/{sample}/bam_sorted/{bid}.sorted.bam")
        for sid in get_sra_ids(sample):
            bams.append(f"output/{sample}/hisat2_aligned/{sid}.sorted.bam")
        for fid in get_fastq_ids(sample):
            bams.append(f"output/{sample}/hisat2_aligned/{fid}.sorted.bam")
        for vid in get_varus_ids(sample):
            bams.append(f"output/{sample}/varus/{vid}.sorted.bam")

    # IsoSeq BAM (isoseq mode only; dual mode adds it in run_genemark_etp_isoseq)
    if mode == 'isoseq':
        isoseq_bam = get_isoseq_bam_for_etp(sample)
        if isoseq_bam:
            bams.append(isoseq_bam)

    return bams


rule run_genemark_etp:
    """GeneMark-ETP.

    Scratch: the whole gmetp.pl work dir (BAM copies, StringTie per
    library, DIAMOND, Spaln, model dirs) and the ProtHint shadow tree are
    on the node-local disk (scripts/tmp_dir.sh, [paths] tmp_dir).
    get_etp_hints.py reads the work dir there; genemark.gtf and
    rnaseq/stringtie/transcripts_merged.gff are copied back to
    output/<sample>/GeneMark-ETP/, training.gtf and hc.gff are copied to
    their output paths. On failure gms.log, etp_config.yaml and the
    StringTie GFFs go to output/<sample>/GeneMark-ETP/failed_run_debug/.
    NEED 1 x BAMs + 10 x genome + 5 x proteins (+ 5 GB each); with less
    free the job works in output/<sample>/GeneMark-ETP/ as before.

    Automatic mode fallback (auto_mode_fallback, planned ETP mode only):
    when GeneMark-ETP stops for lack of training data (no spliced
    transcripts, no gene seeds, or fewer than 100 training genes), the job
    writes the reason to
    output/<sample>/mode_selection/etp_insufficient.txt, leaves empty
    placeholders for its outputs and ends without error; the checkpoint
    check_etp_result (check_evidence.smk) then sends the sample to EP mode.
    Every other failure, and every failure in IsoSeq mode, still fails.
    """
    input:
        genome=lambda wildcards: get_masked_genome(wildcards.sample),
        proteins=lambda wildcards: get_protein_fasta(wildcards.sample),
        bams=_get_etp_bam_files
    output:
        gtf="output/{sample}/GeneMark-ETP/genemark.gtf",
        training="output/{sample}/GeneMark-ETP/training.gtf",
        hc_gff="output/{sample}/GeneMark-ETP/hc.gff",
        etp_hints="output/{sample}/etp_hints.gff",
        stringtie_gff="output/{sample}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff"
    log:
        "logs/{sample}/genemark_etp/genemark_etp.log"
    benchmark:
        "benchmarks/{sample}/genemark_etp/genemark_etp.txt"
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/GeneMark-ETP",
        species_name=lambda wildcards: get_species_name(wildcards),
        fungus="--fungus" if config.get("fungus", False) else "",
        translation_table=config.get("translation_table", 1),
        tmp_root=TMP_ROOT,
        # 1: when GeneMark-ETP stops for lack of training data, write the
        # marker for check_etp_result and end without error, so that the
        # sample continues in EP mode (auto_mode_fallback, planned ETP only)
        fallback=lambda wildcards: 1 if (config.get("auto_mode_fallback", True)
                                         and get_braker_mode(wildcards.sample) == "etp") else 0,
        thin_marker=lambda wildcards: f"output/{wildcards.sample}/mode_selection/etp_insufficient.txt"
    container:
        GENEMARK_ETP_CONTAINER
    shell:
        r"""
        # Disable strict mode: gmetp.pl, get_etp_hints.py, and join_mult_hints.pl
        # may return non-zero on success. Same issue as ProtHint — see README.md.
        # Must disable both -e and pipefail; Snakemake 8 SLURM executor re-enables them.
        set +e
        set +o pipefail
        WORKDIR=$(pwd)
        source {script_dir}/tmp_dir.sh
        mkdir -p {params.outdir}

        finalDir=$(readlink -f {params.outdir})
        rm -rf "$finalDir/failed_run_debug"
        mkdir -p "$(dirname {params.thin_marker})"
        THIN_MARKER=$WORKDIR/{params.thin_marker}
        rm -f "$THIN_MARKER"
        GENOME_ABS=$(readlink -f {input.genome})
        PROTEINS_ABS=$(readlink -f {input.proteins})
        LOG_ABS=$WORKDIR/{log}
        TRAINING_ABS=$WORKDIR/{output.training}
        HC_GFF_ABS=$WORKDIR/{output.hc_gff}
        ETP_HINTS_ABS=$WORKDIR/{output.etp_hints}

        # Step 1: Copy/sort BAM files into etp_data/
        echo "Preparing BAM files for GeneMark-ETP..." > {log}

        # The gmetp.pl work dir runs on the node-local disk; only the outputs
        # are copied back to $finalDir.
        scratch_dir outDir "gmetp_{wildcards.sample}" "{params.tmp_root}" \
            "$(( $(need_gb 1 {input.bams}) + $(need_gb 10 "$GENOME_ABS") + $(need_gb 5 "$PROTEINS_ABS") ))" \
            "$finalDir" 2>> "$LOG_ABS"
        # also the ProtHint shadow tree next to the scratch dir ($outDir.bin)
        trap 'rm -rf -- "$SCRATCH" ${{SCRATCH:+"$SCRATCH.bin"}}' EXIT
        mkdir -p "$outDir/etp_data"

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
            LIB=$(basename $bam .sorted.bam)
            cp $BAM_ABS $outDir/etp_data/${{LIB}}.bam
            if [ -z "$BAM_IDS" ]; then
                BAM_IDS="$LIB"
            else
                BAM_IDS="$BAM_IDS,$LIB"
            fi
            echo "  Prepared BAM: $LIB" >> {log}
        done

        # Step 2: Prepare protein file (remove trailing dots from sequences)
        # Note: gmetp.pl creates a directory named after the protein file basename
        # inside the workdir (e.g. workdir/proteins.fa/). The file must be named
        # proteins.fa for get_etp_hints.py to find the output directory.
        # Place it in a subdirectory to avoid conflict with the created directory.
        # The copy is as large as the protein database: write it to the
        # node-local work dir, not to the run dir.
        mkdir -p "$outDir/protdb"
        PROT_FILE=$outDir/protdb/proteins.fa
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
species: {params.species_name}
translation_table: {params.translation_table}
gcode: {params.translation_table}
YAMLEOF

        echo "YAML config created with rnaseq_sets: [$BAM_IDS]" >> {log}

        # ProtHint (run inside gmetp.pl) has a Perl-threads Spaln dispatcher
        # that can hang after the last batch is enqueued (#98). Run gmetp.pl
        # from a symlinked copy of the container's ETP bin tree in which only
        # ProtHint's run_spliced_alignment.pl is replaced.
        SHADOW_DIR=$outDir.bin
        if bash {script_dir}/make_prothint_shadow.sh $SHADOW_DIR >> {log} 2>&1 && [ -x $SHADOW_DIR/gmetp.pl ]; then
            GMETP=$SHADOW_DIR/gmetp.pl
            echo "[INFO] ProtHint will use the process-based Spaln dispatcher (#98)" >> {log}
        else
            GMETP=gmetp.pl
            echo "[WARNING] could not install the Spaln dispatcher fix (#98), using the container's ProtHint" >> {log}
        fi

        # Step 4: Run GeneMark-ETP
        cd $outDir

        if $GMETP \
            --cfg $outDir/etp_config.yaml \
            --workdir $outDir \
            --bam $outDir/etp_data/ \
            --cores {threads} \
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
            echo "ERROR: GeneMark-ETP failed (exit=$ETP_EXIT), no genemark.gtf" >> {log}
            # Stranded reads assembled unstranded (or poor alignments) can leave
            # StringTie with single-exon transcripts only; GeneMark-ETP then
            # fails much later, in training (#67).
            # THIN_REASON is set when GeneMark-ETP stopped for lack of
            # training data; with auto_mode_fallback the sample then
            # continues in EP mode (see the end of this block).
            THIN_REASON=""
            N_TX=0
            N_MULTI=0
            TGFF=$outDir/rnaseq/stringtie/transcripts_merged.gff
            if [ -s "$TGFF" ]; then
                read -r N_TX N_MULTI < <(awk -F'\t' '$3 == "exon" && match($9, /transcript_id "[^"]+"/) {{ n[substr($9, RSTART, RLENGTH)]++ }} END {{ for (t in n) {{ tx++; if (n[t] > 1) m++ }} print tx + 0, m + 0 }}' "$TGFF")
                echo "DIAGNOSTIC: StringTie assembled $N_TX transcripts, $N_MULTI of them multi-exon." >> {log}
                if [ "$N_TX" -gt 0 ] && [ "$N_MULTI" -eq 0 ]; then
                    echo "HINT: StringTie assembled 0 multi-exon transcripts, so GeneMark-ETP has no spliced genes to train on." >> {log}
                    echo "  Check library strandedness and the alignments: GeneMark-ETP runs StringTie without a strand option (no --rf/--fr), and BAM files of stranded reads need the XS tag (HISAT2 --rna-strandness, STAR --outSAMstrandField intronMotif)." >> {log}
                    echo "  Fallback: ET mode (RNA-Seq only) or EP mode (proteins only) do not use StringTie." >> {log}
                    THIN_REASON="StringTie assembled $N_TX transcripts from the RNA-Seq data, none of them spliced, so GeneMark-ETP had no spliced genes to train on"
                fi
            fi
            # GeneMarkS-T found no usable genes in the assembled transcripts:
            # gms2hints.pl writes an empty seed file, ProtHint stops at once and
            # gmetp.pl reports only "evidence.gff not found".
            GMST_DIR=$outDir/rnaseq/gmst
            NO_SEEDS=0
            if [ -f "$GMST_DIR/genome_gmst.gtf" ] && [ ! -s "$GMST_DIR/genome_gmst.gtf" ]; then
                NO_SEEDS=1
                N_GMST=$(awk -F'\t' '$3 == "CDS"' "$GMST_DIR/transcripts_merged.fasta.gff" 2>/dev/null | wc -l)
                echo "DIAGNOSTIC: GeneMarkS-T predicted $N_GMST genes in the assembled transcripts, but none of them maps to the genome as a usable gene (rnaseq/gmst/genome_gmst.gtf is empty)." >> {log}
                if grep -q "no start sites in training" $WORKDIR/{log} "$GMST_DIR/gms.log" 2>/dev/null; then
                    echo "DIAGNOSTIC: GeneMarkS-T found no start codons to train its model on (\"no start sites in training\")." >> {log}
                fi
                echo "HINT: The RNA-Seq data contain too few complete transcripts of this genome's genes for GeneMark-ETP. ProtHint then gets no gene seeds, so \"evidence.gff not found\" above is a consequence, not the cause." >> {log}
                echo "  Typical reasons: few reads from this species (SRA runs from host tissue, a mixed sample or a related species), low coverage, or degraded RNA." >> {log}
                echo "  Check the fraction of uniquely mapped reads (with VARUS: avgUmrPercent in output/{wildcards.sample}/varus/RunStatistics.csv)." >> {log}
                echo "  Without the automatic mode fallback: run EP mode (proteins only) by leaving the RNA-Seq columns of samples.csv empty." >> {log}
                THIN_REASON="GeneMarkS-T found no usable genes in the $N_TX transcripts assembled from the RNA-Seq data"
            fi
            # Number of training genes GeneMark-ETP selected (gmetp.pl prints
            # "# from file training.list parsed IDs: N"). With a handful of
            # genes training fails in different places: division by zero in
            # parse_set.pl, or no model and then no pred_m/genemark.gtf.
            N_TRAIN=$(grep "from file training.list parsed IDs:" $WORKDIR/{log} | tail -1 | awk '{{print $NF}}')
            case "$N_TRAIN" in
                ''|*[!0-9]*) N_TRAIN=unknown ;;
                *) if [ "$N_TRAIN" -lt 100 ]; then
                       echo "DIAGNOSTIC: GeneMark-ETP selected only $N_TRAIN training genes; about 100 spliced genes are needed to train." >> {log}
                       if [ -z "$THIN_REASON" ]; then
                           THIN_REASON="GeneMark-ETP found only $N_TRAIN training genes ($N_TX transcripts assembled from the RNA-Seq data, $N_MULTI of them spliced)"
                       fi
                   fi ;;
            esac
            if [ "$NO_SEEDS" -eq 1 ]; then
                :
            elif grep -q "Illegal division by zero" $WORKDIR/{log} 2>/dev/null; then
                echo "HINT: GeneMark-ETP crashed in model training (parse_set.pl / train_super.pl division by zero)." >> {log}
                echo "  Training genes found: $N_TRAIN (minimum needed: ~100 multi-exon genes)." >> {log}
                echo "  Likely causes:" >> {log}
                echo "    1. Fungal organism: set fungus = 1 in the [PARAMS] section of config.ini" >> {log}
                echo "    2. Low RNA-seq coverage: too few reads to build reliable gene models" >> {log}
                echo "    3. Protein database too distant: ProtHint yields too few HC introns" >> {log}
                echo "  See https://github.com/gatech-genemark/GeneMark-ETP/issues" >> {log}
            else
                TSEQ=$outDir/rnaseq/stringtie/transcripts_merged.fasta
                if [ -f "$TSEQ" ]; then
                    TSEQ_SIZE=$(wc -c < "$TSEQ")
                    echo "DIAGNOSTIC: transcripts_merged.fasta size: $TSEQ_SIZE bytes" >> {log}
                    if [ "$TSEQ_SIZE" -eq 0 ]; then
                        echo "HINT: transcripts_merged.fasta is empty. Individual per-library GFFs were produced by StringTie, but gffread extracted no sequences." >> {log}
                        echo "  Most likely causes:" >> {log}
                        echo "    1. Chromosome name mismatch: contig names in the BAM/GFF differ from those in the genome FASTA." >> {log}
                        echo "       Check: samtools view -H <bam> | grep ^@SQ  vs  grep '>' <genome> | head" >> {log}
                        echo "    2. StringTie merge filtered all transcripts (low coverage/FPKM across libraries)." >> {log}
                        echo "  File an issue at https://github.com/gatech-genemark/GeneMark-ETP/issues with this log." >> {log}
                        echo "  Fallback: use ET mode (RNA-Seq only) or EP mode (protein only) instead of ETP." >> {log}
                    else
                        if [ "$ETP_EXIT" -eq 139 ]; then
                            echo "HINT: exit 139 = segmentation fault inside GeneMark-ETP." >> {log}
                        else
                            echo "HINT: GeneMark-ETP exited with code $ETP_EXIT without writing genemark.gtf." >> {log}
                        fi
                        echo "  Common causes (check these before subsampling):" >> {log}
                        echo "    1. Thread count: GeneMark-ETP ran with --cores {threads}. Very high counts (e.g. 256) are known to crash it; lower cpus_per_task in the [SLURM_ARGS] section of config.ini." >> {log}
                        echo "    2. Sequence name mismatch: @SQ names in the BAM header differ from the genome FASTA headers." >> {log}
                        echo "       Check: samtools view -H <bam> | grep ^@SQ  vs  grep '>' <genome> | head" >> {log}
                        echo "    3. Stranded library: GeneMark-ETP assembles transcripts unstranded, which can produce conflicting models from stranded data." >> {log}
                        echo "    4. Very large transcript set: only if 1-3 are ruled out, subsample the RNA-Seq BAM and rerun." >> {log}
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
            # A technical failure is not thin evidence (as in
            # scripts/genemark_fallback.sh for ET and EP): killed by a
            # signal, a crash of a program that GeneMark-ETP called (e.g.
            # "Error on last system call, error code 11" of gmst.pl), full
            # disk or exhausted memory.
            if [ -n "$THIN_REASON" ]; then
                if [ "$ETP_EXIT" -ge 128 ] || grep -qiE 'No space left on device|Disk quota exceeded|Cannot allocate memory|Out of memory|Segmentation fault|core dumped|error code (9|11|134|137|139)([^0-9]|$)' $WORKDIR/{log} 2>/dev/null; then
                    echo "DIAGNOSTIC: GeneMark-ETP or a program it called was killed, crashed or ran out of disk space or memory (exit=$ETP_EXIT). This is a technical failure, so there is no automatic mode fallback." >> {log}
                    THIN_REASON=""
                fi
            fi
            if [ "{params.fallback}" = "1" ] && [ -n "$THIN_REASON" ]; then
                # Not enough RNA-Seq evidence for GeneMark-ETP: leave the
                # reason for check_etp_result and empty placeholders for
                # the outputs of this rule (not read in EP mode).
                echo "$THIN_REASON" > "$THIN_MARKER"
                echo "FALLBACK: $THIN_REASON." >> {log}
                echo "FALLBACK: auto_mode_fallback is on, the sample continues in EP mode (proteins only). Set auto_mode_fallback = 0 in config.ini to fail here instead." >> {log}
                mkdir -p "$finalDir/rnaseq/stringtie"
                for f in "$finalDir/genemark.gtf" "$TRAINING_ABS" "$HC_GFF_ABS" "$ETP_HINTS_ABS" "$finalDir/rnaseq/stringtie/transcripts_merged.gff"; do
                    [ -e "$f" ] || : > "$f"
                done
                rm -rf $SHADOW_DIR
                exit 0
            fi
            exit 1
        fi

        # genemark.gtf from gmetp.pl has no gene lines: count distinct gene_id
        n_genes=$(awk 'match($0, /gene_id "[^"]+"/) {{ id = substr($0, RSTART, RLENGTH); if (!(id in seen)) {{ seen[id] = 1; n++ }} }} END {{ print n + 0 }}' $outDir/genemark.gtf)
        echo "GeneMark-ETP predicted $n_genes genes (exit=$ETP_EXIT)" >> {log}

        # Step 5: Find and copy training genes and HC genes
        # GeneMark-ETP creates a directory named after the protein file
        # (e.g. proteins.fa/) containing model/training.gtf and model/hc.gff
        ETP_MODEL=$(find $outDir -path "*/model/training.gtf" -not -path "*/etr/*" | head -1 | xargs dirname 2>/dev/null || echo "")

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/training.gtf" ]; then
            cp "$ETP_MODEL/training.gtf" "$TRAINING_ABS"
            echo "Found training.gtf at $ETP_MODEL" >> {log}
        else
            echo "WARNING: No model/training.gtf found, using genemark.gtf" >> {log}
            cp $outDir/genemark.gtf "$TRAINING_ABS"
        fi

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/hc.gff" ]; then
            cp "$ETP_MODEL/hc.gff" "$HC_GFF_ABS"
            echo "Found hc.gff at $ETP_MODEL" >> {log}
        else
            echo "WARNING: No model/hc.gff found" >> {log}
            touch "$HC_GFF_ABS"
        fi

        # Step 6: Extract hints for AUGUSTUS using get_etp_hints.py
        echo "Extracting hints for AUGUSTUS..." >> {log}

        # CRITICAL: --genemark_scripts must point at /opt/ETP/bin (where
        # format_back.pl lives), NOT /opt/ETP/bin/gmes/. The wrong path
        # silently triggers ScriptNotFound inside get_etp_hints.py and
        # falls through to the manual fallback below, which produces a
        # subtly wrong hintsfile (raw nonhc coordinates that should have
        # been transformed via nonhc.trace by format_back.pl, plus only
        # one copy of hintsfile_merged.gff instead of two). On the
        # athaliana ETP benchmark this silent fallback cost ~3pp of
        # locus-level sensitivity vs native braker.pl. Match braker.pl's
        # GeneMark_PATH which resolves to /opt/ETP/bin in the container.
        #
        # Note: get_etp_hints.py uses >> (append) for ALL its output
        # writes, so any leftover file from a previous attempt would
        # be appended to. Truncate first to make the rule re-runnable.
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
        # (see comment above) and silently degraded ETP accuracy. If
        # get_etp_hints.py fails, the right action is to surface the
        # error and let the user investigate, not paper over it.
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

        # NOTE: do NOT call join_mult_hints.pl here. braker.pl's
        # get_etp_hints_for_Augustus (scripts/braker.pl line 5638) calls
        # get_etp_hints.py and then defers all multiplicity joining to its
        # join_mult_hints sub (line 4769), which has the critical
        # src=C-with-grp= split. Running join_mult_hints.pl on the raw
        # etp_hints stream destroys the per-gene grp= linkage that AUGUSTUS
        # uses to score whole-gene hint coverage and silently costs ~3
        # percentage points of locus-level sensitivity in ETP mode. The
        # downstream merge_hints rule does the join correctly with the
        # src=C grp= split — leave the file uncollapsed here.
        n_hints=$(wc -l < {output.etp_hints})
        echo "Extracted $n_hints hints for AUGUSTUS" >> {log}
        echo "GeneMark-ETP completed successfully" >> {log}

        # Record software versions
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        GMETP_VER=$(grep -oP 'my \$version = "\K[^"]+' $(which gmetp.pl) 2>/dev/null || true)
        GM_COMMIT=$(grep 'refs/remotes/origin/main' /opt/ETP/.git/packed-refs 2>/dev/null | awk '{{print substr($1,1,7)}}' || true)
        ST_VER=$(stringtie --version 2>&1 || true)
        DM_VER=$(diamond version 2>&1 | awk '{{print $NF}}' || true)
        ( flock 9
          printf "GeneMark-ETP\t%s (commit %s)\n" "$GMETP_VER" "$GM_COMMIT" >> "$VERSIONS_FILE"
          printf "StringTie\t%s\n" "$ST_VER" >> "$VERSIONS_FILE"
          printf "DIAMOND\t%s\n" "$DM_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh || true
        cite genemark_etp "$REPORT_DIR" || true
        cite genemarks_t "$REPORT_DIR" || true
        cite braker3 "$REPORT_DIR" || true

        # Remove GeneMark-ETP internal working files and etp_data/ BAM copies.
        # Tracked outputs kept: genemark.gtf, training.gtf, hc.gff,
        # rnaseq/stringtie/transcripts_merged.gff (etp_hints.gff is outside outdir).
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
