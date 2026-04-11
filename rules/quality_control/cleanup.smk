
rule cleanup_tmp_opt:
    """
    Clean up temporary optimization files created by optimize_augustus.pl.

    The optimize_augustus.pl script creates a tmp_opt_<species> directory
    in the AUGUSTUS config species directory during parameter optimization.
    This directory contains temporary files used during cross-validation
    and is no longer needed after optimization completes.

    This rule deletes the tmp_opt directory to save disk space.

    BRAKER performs this cleanup automatically after successful completion
    (see braker.pl:10133-10136).

    Resources:
        - Local rule (instant execution)
        - No SLURM submission needed

    Input:
        optimize_log: Ensures optimization is complete before cleanup

    Output:
        cleanup_marker: Marker file indicating cleanup was performed
    """
    input:
        optimize_log = "output/{sample}/optimize_augustus.log"
    output:
        cleanup_marker = "output/{sample}/.tmp_opt_cleaned"
    benchmark:
        "benchmarks/{sample}/cleanup_tmp_opt/cleanup_tmp_opt.txt"
    params:
        species_name = lambda w: get_species_name(w),
        tmp_opt_dir = lambda w: f"augustus_config/species/{get_species_name(w)}/tmp_opt_{get_species_name(w)}"
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] ===== CLEANING UP TEMPORARY OPTIMIZATION FILES ====="
        echo "[INFO] Species: {params.species_name}"
        echo "[INFO] tmp_opt directory: {params.tmp_opt_dir}"

        if [ -d "{params.tmp_opt_dir}" ]; then
            echo "[INFO] Removing tmp_opt directory..."
            rm -rf "{params.tmp_opt_dir}"
            echo "[INFO] Successfully removed {params.tmp_opt_dir}"
        else
            echo "[INFO] tmp_opt directory does not exist (already cleaned or optimization was skipped)"
        fi

        # Create marker file
        touch {output.cleanup_marker}
        echo "[INFO] Cleanup completed"
        echo "[INFO] ======================================="
        """
