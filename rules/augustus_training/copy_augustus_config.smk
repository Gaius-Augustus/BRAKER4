
rule copy_augustus_config:
    """
    Copy AUGUSTUS configuration from container to local directory.

    This rule extracts the AUGUSTUS config directory from the BRAKER3 container
    to the local filesystem, making it writable for subsequent modifications.
    Since containers are read-only, we need a local copy to create new species
    and train AUGUSTUS models.

    This is a local rule (runs on head node, not submitted to SLURM) because:
    - It's a simple file copy operation
    - Needs to run only once at workflow start
    - Doesn't require significant compute resources

    Output:
        config_dir: Directory containing AUGUSTUS configuration files
        marker: Marker file indicating successful copy (prevents re-copying)
    """
    output:
        config_dir = directory("augustus_config"),
        marker = "augustus_config/.copy_complete"
    benchmark:
        "benchmarks/copy_augustus_config.txt"
    container:
        BRAKER3_CONTAINER
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    shell:
        r"""
        set -euo pipefail

        # Use file-based locking to prevent race conditions when multiple jobs run in parallel
        LOCKFILE="{output.config_dir}.lock"
        LOCKFD=200

        # Acquire exclusive lock (waits if another process has the lock)
        exec 200>"$LOCKFILE"
        flock -x 200

        echo "[INFO] Lock acquired for AUGUSTUS config copy"

        # Check if the directory already exists and is complete
        if [ -d "{output.config_dir}" ] && [ -f "{output.marker}" ]; then
            echo "[INFO] AUGUSTUS config directory already exists. Validating integrity..."

            # Validate critical subdirectories exist
            CRITICAL_DIRS="species model extrinsic parameters"
            VALID=true
            for dir in $CRITICAL_DIRS; do
                if [ ! -d "{output.config_dir}/$dir" ]; then
                    echo "[WARN] Critical directory missing: $dir"
                    VALID=false
                    break
                fi
            done

            # Validate minimum number of species (container has ~145 species)
            if [ "$VALID" = true ]; then
                SPECIES_COUNT=$(ls -1 {output.config_dir}/species 2>/dev/null | wc -l)
                if [ "$SPECIES_COUNT" -lt 100 ]; then
                    echo "[WARN] Species directory incomplete (found $SPECIES_COUNT, expected >100)"
                    VALID=false
                fi
            fi

            # If validation passed, skip copy
            if [ "$VALID" = true ]; then
                echo "[INFO] Config directory validation passed. Skipping copy."
                flock -u 200
                rm -f "$LOCKFILE"
                exit 0
            else
                echo "[WARN] Config directory validation FAILED. Will re-copy."
                rm -rf "{output.config_dir}"
            fi
        fi

        echo "[INFO] Copying AUGUSTUS config from container to local directory"

        # Remove any incomplete previous attempts (check for actual content, not just existence)
        # Note: Snakemake creates empty directory before running shell, so we check for files
        if [ -d "{output.config_dir}" ] && [ "$(ls -A {output.config_dir} 2>/dev/null)" ]; then
            echo "[WARN] Found incomplete config directory with files. Removing..."
            rm -rf "{output.config_dir}"
        fi

        # Copy the AUGUSTUS config directory from container
        # If directory exists but is empty (created by Snakemake), remove and recreate
        if [ -d "{output.config_dir}" ]; then
            rmdir "{output.config_dir}"
        fi
        cp -r /opt/Augustus/config "{output.config_dir}"

        # Make it writable so we can modify it later
        chmod -R u+w "{output.config_dir}"

        # Validate the copy was successful before creating marker
        echo "[INFO] Validating copied config directory..."
        COPY_VALID=true

        # Check critical directories exist
        CRITICAL_DIRS="species model extrinsic parameters"
        for dir in $CRITICAL_DIRS; do
            if [ ! -d "{output.config_dir}/$dir" ]; then
                echo "[ERROR] Critical directory missing after copy: $dir"
                COPY_VALID=false
                break
            fi
        done

        # Check minimum number of species
        if [ "$COPY_VALID" = true ]; then
            SPECIES_COUNT=$(ls -1 {output.config_dir}/species 2>/dev/null | wc -l)
            if [ "$SPECIES_COUNT" -lt 100 ]; then
                echo "[ERROR] Species directory incomplete after copy (found $SPECIES_COUNT, expected >100)"
                COPY_VALID=false
            fi
        fi

        # Check total file count is reasonable (container has ~850 files)
        if [ "$COPY_VALID" = true ]; then
            FILE_COUNT=$(find {output.config_dir} -type f | wc -l)
            if [ "$FILE_COUNT" -lt 500 ]; then
                echo "[ERROR] Too few files copied (found $FILE_COUNT, expected >500)"
                COPY_VALID=false
            fi
        fi

        # Only create marker if validation passed
        if [ "$COPY_VALID" = true ]; then
            touch "{output.marker}"
            echo "[INFO] AUGUSTUS config successfully copied and validated"
            echo "[INFO] Total files copied: $(find {output.config_dir} -type f | wc -l)"
        else
            echo "[ERROR] Config copy validation FAILED. Removing incomplete copy."
            rm -rf "{output.config_dir}"
            flock -u 200
            rm -f "$LOCKFILE"
            exit 1
        fi

        # Release lock and remove lockfile
        flock -u 200
        rm -f "$LOCKFILE"
        """
