
rule create_new_species:
    """
    Create new AUGUSTUS species parameter set.

    This rule initializes a new species in the AUGUSTUS config directory by running
    new_species.pl. Two names are involved, on purpose (#108):
      - the AUGUSTUS species is get_species_name(sample), i.e. the sample name
        with "-" and "." replaced by "_"; its parameter files are written to
        <augustus_config_path>/species/<species>/
      - augustus_config/species/{sample}_rescuedB3/ in the run directory only
        holds Snakemake bookkeeping files (.species_created, the parameter
        backup); it is not an AUGUSTUS species.
    new_species.pl is skipped if the species already has a parameters file.

    Input:
        augustus_config: Local copy of AUGUSTUS config directory (dependency)

    Output:
        species_dir: Bookkeeping directory ({sample}_rescuedB3), not the species itself
        marker: Marker file indicating successful species creation

    Side effects:
        Creates new species parameter files in augustus_config/species/{species_name}/
    """
    input:
        augustus_config = "augustus_config/.copy_complete"
    output:
        species_dir = directory("augustus_config/species/{sample}_rescuedB3"),
        marker = "augustus_config/species/{sample}_rescuedB3/.species_created"
    log:
        "logs/{sample}/create_new_species/create_new_species.log"
    benchmark:
        "benchmarks/{sample}/create_new_species/create_new_species.txt"
    params:
        aug_config = augustus_config_path,
        species_name = lambda w: get_species_name(w)
    container:
        BRAKER3_CONTAINER
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        # Set AUGUSTUS_CONFIG_PATH to our local copy (absolute path)
        export AUGUSTUS_CONFIG_PATH={params.aug_config}

        # Use file-based locking to prevent race conditions when creating species
        LOCKFILE="$AUGUSTUS_CONFIG_PATH/species/.species_creation.lock"
        exec 200>"$LOCKFILE"
        flock -x 200

        echo "[INFO] Lock acquired for species creation"
        echo "[INFO] AUGUSTUS_CONFIG_PATH set to: $AUGUSTUS_CONFIG_PATH"
        echo "[INFO] Creating new AUGUSTUS species: {params.species_name}"
        echo "[INFO] Current working directory: $(pwd)"
        echo "[INFO] Checking if config directory is writable..."
        test -w "$AUGUSTUS_CONFIG_PATH" && echo "[INFO] Config directory is writable" || echo "[ERROR] Config directory is NOT writable"

        # Check if species parameter files already exist
        if [ -f "$AUGUSTUS_CONFIG_PATH/species/{params.species_name}/{params.species_name}_parameters.cfg" ]; then
            echo "[INFO] Species {params.species_name} already exists with parameter files. Skipping creation."
        else
            echo "[INFO] Running new_species.pl..."
            # Create new species using new_species.pl
            if ! new_species.pl --species={params.species_name}; then
                echo "[ERROR] new_species.pl failed for species {params.species_name} in $AUGUSTUS_CONFIG_PATH/species" >&2
                echo "  Check that $AUGUSTUS_CONFIG_PATH is writable and contains species/generic, and that no" >&2
                echo "  incomplete $AUGUSTUS_CONFIG_PATH/species/{params.species_name} is left over from an earlier run." >&2
                exit 1
            fi
            if [ ! -f "$AUGUSTUS_CONFIG_PATH/species/{params.species_name}/{params.species_name}_parameters.cfg" ]; then
                echo "[ERROR] new_species.pl did not write {params.species_name}_parameters.cfg in $AUGUSTUS_CONFIG_PATH/species/{params.species_name}" >&2
                exit 1
            fi

            echo "[INFO] Successfully created species {params.species_name}"
        fi

        # Create marker file
        touch {output.marker}

        echo "[INFO] Species directory: $AUGUSTUS_CONFIG_PATH/species/{params.species_name}"
        echo "[INFO] Contents:"
        ls -la "$AUGUSTUS_CONFIG_PATH/species/{params.species_name}" | head -10 || true

        # Release lock
        flock -u 200
        rm -f "$LOCKFILE"
        """
