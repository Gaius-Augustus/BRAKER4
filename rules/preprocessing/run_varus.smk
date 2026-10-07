"""
Run pyVARUS to automatically select and download RNA-Seq data from SRA.

pyVARUS:
1. Searches NCBI SRA for RNA-Seq data matching the species
2. Iteratively downloads and aligns complementary RNA-Seq reads
3. Produces a coordinate-sorted BAM file for use in gene prediction

Container: gaiusaugustus/pyvarus:v2.0.0a0 (built from
https://github.com/Gaius-Augustus/pyVARUS/tree/main/docker; BRAKER4 pins a
version tag, bump it in Snakefile, rules/common.smk, config.ini.example and
README.md together)
"""


def get_varus_genus(sample):
    """Get VARUS genus for a sample."""
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    return row["varus_genus"]


def get_varus_species(sample):
    """Get VARUS species for a sample."""
    row = samples_df[samples_df["sample_name"] == sample].iloc[0]
    return row["varus_species"]


rule run_varus:
    """Run VARUS to auto-select, download, and align RNA-Seq from SRA.

    Scratch: the pyVARUS outdir (HISAT2 index, batches/, VARUS.bam, the
    sorted BAM) is a private directory on the node-local disk
    (scripts/tmp_dir.sh, [paths] tmp_dir). The sorted BAM and its .csi go
    to their output paths; Coverage.csv, RunStatistics.csv and introns.gff
    are copied back to output/<sample>/varus/. NEED 10 x genome + 55 GB;
    with less free the job works in output/<sample>/varus/ as before.
    """
    input:
        genome=lambda wildcards: get_genome(wildcards.sample)
    output:
        bam="output/{sample}/varus/varus.sorted.bam",
        csi="output/{sample}/varus/varus.sorted.bam.csi"
    log:
        "logs/{sample}/varus/varus.log"
    benchmark:
        "benchmarks/{sample}/varus/varus.txt"
    params:
        genus=lambda wildcards: get_varus_genus(wildcards.sample),
        species=lambda wildcards: get_varus_species(wildcards.sample),
        varus_dir=lambda wildcards: f"output/{wildcards.sample}/varus",
        use_logan=1 if config.get('varus_logan', False) else 0,
        wrapper=os.path.join(script_dir, "run_varus_wrapper.sh"),
        tmp_root=TMP_ROOT
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        VARUS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        source {script_dir}/tmp_dir.sh
        WORKDIR=$PWD
        finalDir=$WORKDIR/{params.varus_dir}
        BAM_ABS=$WORKDIR/{output.bam}
        mkdir -p "$finalDir" "$(dirname $WORKDIR/{log})"
        GENOME_ABS=$(readlink -f {input.genome})

        # pyVARUS works on the node-local disk; downloads are not known in
        # advance, hence the fixed 50 GB on top of the index.
        : > {log}
        scratch_dir outDir "varus_{wildcards.sample}" "{params.tmp_root}" \
            "$(( $(need_gb 10 "$GENOME_ABS") + 50 ))" "$finalDir" 2>> {log}
        trap 'rm -rf -- "$SCRATCH"' EXIT

        bash {params.wrapper} \
            "$outDir" \
            {input.genome} \
            {params.genus} \
            {params.species} \
            {threads} \
            "$BAM_ABS" \
            {log} \
            {params.use_logan}

        copy_back "$outDir" "$finalDir" Coverage.csv RunStatistics.csv introns.gff
        if [ ! -s "$BAM_ABS" ] || [ ! -s "$BAM_ABS.csi" ]; then
            echo "[ERROR] {output.bam} or its .csi missing after the copy back" >> {log}
            exit 1
        fi

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        ( flock 9; printf "pyVARUS\tv2.0.0a0\n" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite varus "$REPORT_DIR"

        # Remove pyVARUS working files: per-batch tree, HISAT2 index, unsorted BAM.
        #
        # batches/<acc>/N<n>X<x>/ holds one directory per downloaded batch. pyVARUS
        # unlinks the FASTAs after alignment and the BAMs after the final merge, but
        # never removes the directories or the per-batch Log.final.out, so a default
        # run (--max-batches 1000) leaves ~2000 inodes behind. Under the old C++ VARUS
        # these lived inside <Genus>_<species>/ and were swallowed by that rm -rf; the
        # pyVARUS switch moved them to the top level and they lost their coverage.
        #
        # Keep: varus.sorted.bam, .csi, varus_stats.txt, varus_runlist.tsv,
        #       Coverage.csv, RunStatistics.csv, introns.gff (small, diagnostic).
        # Fallback only; on scratch the EXIT trap removes the work dir.
        if [ -z "$SCRATCH" ]; then
            VARUS_DIR_ABS=$(readlink -f output/{wildcards.sample}/varus)
            rm -rf "$VARUS_DIR_ABS/batches"    2>/dev/null || true
            rm -rf "$VARUS_DIR_ABS/genome"     2>/dev/null || true
            rm -f  "$VARUS_DIR_ABS/VARUS.bam"  2>/dev/null || true
            rm -f  "$VARUS_DIR_ABS/Runlist.tsv" 2>/dev/null || true
            rm -f  "$VARUS_DIR_ABS/intronDB.splice_sites" \
                   "$VARUS_DIR_ABS/intronDB.junc.bed" 2>/dev/null || true
        fi
        """
