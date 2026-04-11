
rule create_new_species:
    """
    Create new AUGUSTUS species parameter set.

    This rule initializes a new species in the AUGUSTUS config directory by running
    new_species.pl. The species name is derived from the output directory basename
    plus "_rescuedB3" suffix. The rule only executes if the species directory doesn't
    already exist, preventing accidental overwrites of existing parameter sets.

    This is a local rule because it's a quick configuration step that doesn't require
    cluster resources.

    Input:
        augustus_config: Local copy of AUGUSTUS config directory (dependency)

    Output:
        species_dir: Directory containing species-specific AUGUSTUS parameters
        marker: Marker file indicating successful species creation

    Side effects:
        Creates new species parameter files in augustus_config/species/{species_name}/
    """
    input:
        augustus_config = "augustus_config/.copy_complete"
    output:
        species_dir = directory("augustus_config/species/{sample}_rescuedB3"),
        marker = "augustus_config/species/{sample}_rescuedB3/.species_created"
    benchmark:
        "benchmarks/{sample}/create_new_species/create_new_species.txt"
    params:
        aug_config = augustus_config_path,
        species_name = lambda w: get_species_name(w)
    container:
        BRAKER3_CONTAINER
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail

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
            new_species.pl --species={params.species_name}

            echo "[INFO] Successfully created species {params.species_name}"
        fi

        # Create marker file
        touch {output.marker}

        echo "[INFO] Species directory: $AUGUSTUS_CONFIG_PATH/species/{params.species_name}"
        echo "[INFO] Contents:"
        ls -la "$AUGUSTUS_CONFIG_PATH/species/{params.species_name}" | head -10

        # Release lock
        flock -u 200
        rm -f "$LOCKFILE"
        """
