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

Container: teambraker/braker3:latest (contains gmetp.pl, get_etp_hints.py)
"""


def _get_etp_bam_files(wildcards):
    """Get sorted BAM files for GeneMark-ETP (short-read + IsoSeq for non-dual, short-read only for dual).

    In dual mode, IsoSeq goes through a separate ETP run (run_genemark_etp_isoseq).
    In isoseq mode, only IsoSeq BAMs are used.
    In etp mode, only short-read BAMs are used.
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

    # IsoSeq BAMs (for isoseq mode only — dual mode has separate rule)
    if mode == 'isoseq':
        isoseq_bam = get_isoseq_bam_for_etp(sample)
        if isoseq_bam:
            bams.append(isoseq_bam)

    return bams


rule run_genemark_etp:
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
    threads: workflow.cores
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    params:
        outdir=lambda wildcards: f"output/{wildcards.sample}/GeneMark-ETP",
        species_name=lambda wildcards: get_species_name(wildcards),
        fungus="--fungus" if config.get("fungus", False) else ""
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
        mkdir -p {params.outdir}/etp_data

        OUTDIR_ABS=$(readlink -f {params.outdir})
        GENOME_ABS=$(readlink -f {input.genome})
        PROTEINS_ABS=$(readlink -f {input.proteins})

        # Step 1: Copy/sort BAM files into etp_data/
        echo "Preparing BAM files for GeneMark-ETP..." > {log}
        BAM_IDS=""
        for bam in {input.bams}; do
            BAM_ABS=$(readlink -f $bam)
            LIB=$(basename $bam .sorted.bam)
            cp $BAM_ABS $OUTDIR_ABS/etp_data/${{LIB}}.bam
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
        # Place it one level up to avoid conflict with the created directory.
        PROT_FILE=$WORKDIR/output/{wildcards.sample}/proteins.fa
        sed '/^>/!s/\\.$//' $PROTEINS_ABS > $PROT_FILE

        # Step 3: Create YAML config
        cat > $OUTDIR_ABS/etp_config.yaml << YAMLEOF
---
RepeatMasker_path: ''
annot_path: ''
genome_path: $GENOME_ABS
protdb_path: $(readlink -f $PROT_FILE)
rnaseq_sets: [$BAM_IDS]
species: {params.species_name}
YAMLEOF

        echo "YAML config created with rnaseq_sets: [$BAM_IDS]" >> {log}

        # Step 4: Run GeneMark-ETP
        cd $OUTDIR_ABS

        if gmetp.pl \
            --cfg $OUTDIR_ABS/etp_config.yaml \
            --workdir $OUTDIR_ABS \
            --bam $OUTDIR_ABS/etp_data/ \
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

        if [ ! -f $OUTDIR_ABS/genemark.gtf ]; then
            echo "ERROR: GeneMark-ETP failed (exit=$ETP_EXIT), no genemark.gtf" >> {log}
            exit 1
        fi

        n_genes=$(grep -c $'\\tgene\\t' $OUTDIR_ABS/genemark.gtf || echo "0")
        echo "GeneMark-ETP predicted $n_genes genes (exit=$ETP_EXIT)" >> {log}

        # Step 5: Find and copy training genes and HC genes
        # GeneMark-ETP creates a directory named after the protein file
        # (e.g. proteins.fa/) containing model/training.gtf and model/hc.gff
        ETP_MODEL=$(find $OUTDIR_ABS -path "*/model/training.gtf" -not -path "*/etr/*" | head -1 | xargs dirname 2>/dev/null || echo "")

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/training.gtf" ]; then
            cp "$ETP_MODEL/training.gtf" {output.training}
            echo "Found training.gtf at $ETP_MODEL" >> {log}
        else
            echo "WARNING: No model/training.gtf found, using genemark.gtf" >> {log}
            cp $OUTDIR_ABS/genemark.gtf {output.training}
        fi

        if [ -n "$ETP_MODEL" ] && [ -f "$ETP_MODEL/hc.gff" ]; then
            cp "$ETP_MODEL/hc.gff" {output.hc_gff}
            echo "Found hc.gff at $ETP_MODEL" >> {log}
        else
            echo "WARNING: No model/hc.gff found" >> {log}
            touch {output.hc_gff}
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
        # (see comment above) and silently degraded ETP accuracy. If
        # get_etp_hints.py fails, the right action is to surface the
        # error and let the user investigate, not paper over it.
        if [ ! -s {output.etp_hints} ]; then
            echo "[ERROR] get_etp_hints.py failed (exit=$HINTS_EXIT) and produced no hints." >> {log}
            echo "[ERROR] Inspect $OUTDIR_ABS for the GeneMark-ETP run state." >> {log}
            exit 1
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
        """
