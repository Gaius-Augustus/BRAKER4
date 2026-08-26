"""
Run pyVARUS to automatically select and download RNA-Seq data from SRA.

pyVARUS:
1. Searches NCBI SRA for RNA-Seq data matching the species
2. Iteratively downloads and aligns complementary RNA-Seq reads
3. Produces a coordinate-sorted BAM file for use in gene prediction

Container: katharinahoff/pyvarus:v1.0.0
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
    """Run VARUS to auto-select, download, and align RNA-Seq from SRA."""
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
        wrapper=os.path.join(script_dir, "run_varus_wrapper.sh")
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        VARUS_CONTAINER
    shell:
        r"""
        set -euo pipefail
        bash {params.wrapper} \
            {params.varus_dir} \
            {input.genome} \
            {params.genus} \
            {params.species} \
            {threads} \
            {output.bam} \
            {log}

        # Record software version
        VERSIONS_FILE=output/{wildcards.sample}/software_versions.tsv
        ( flock 9; printf "pyVARUS\tv1.0.0\n" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

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
        VARUS_DIR_ABS=$(readlink -f output/{wildcards.sample}/varus)
        rm -rf "$VARUS_DIR_ABS/batches"    2>/dev/null || true
        rm -rf "$VARUS_DIR_ABS/genome"     2>/dev/null || true
        rm -f  "$VARUS_DIR_ABS/VARUS.bam"  2>/dev/null || true
        rm -f  "$VARUS_DIR_ABS/Runlist.tsv" 2>/dev/null || true
        rm -f  "$VARUS_DIR_ABS/intronDB.splice_sites" \
               "$VARUS_DIR_ABS/intronDB.junc.bed" 2>/dev/null || true
        """
