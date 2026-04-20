"""
Run VARUS to automatically select and download RNA-Seq data from SRA.

VARUS automatically:
1. Searches SRA for RNA-Seq data matching the species
2. Downloads and aligns complementary RNA-Seq reads
3. Produces BAM file for use in gene prediction

Important: VARUS internally changes directories (chdir), which conflicts
with Snakemake's working directory management. Therefore, VARUS must be
invoked from a wrapper script that resolves all paths to absolute paths
before cd-ing into the VARUS working directory.

Input:
    - Genome FASTA file

Output:
    - Coordinate-sorted BAM file with index

Container: katharinahoff/varus-notebook:v0.0.5
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
        bai="output/{sample}/varus/varus.sorted.bam.bai"
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
        ( flock 9; printf "VARUS\tcontainer v0.0.6\n" >> "$VERSIONS_FILE" ) 9>"$VERSIONS_FILE.lock"

        # Report
        REPORT_DIR=output/{wildcards.sample}
        source {script_dir}/report_citations.sh
        cite varus "$REPORT_DIR"
        """
