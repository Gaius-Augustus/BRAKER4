rule run_masking:
    """
    Run RepeatModeler2 + RepeatMasker on genome.

    Uses /dev/shm or /tmp for HPC-friendly I/O.
    The working directory is cleaned up automatically on exit.

    Input:
        genome: Input genome FASTA file

    Output:
        masked_genome: Soft-masked genome (lowercase for repeats)
        marker: Completion marker file

    Container:
        dfam/tetools:latest (contains RepeatModeler and RepeatMasker)
    """
    input:
        genome = lambda wildcards: get_genome(wildcards.sample)
    output:
        masked_genome = "output/{sample}/preprocessing/genome.fa.masked",
        marker = "output/{sample}/preprocessing/.masking_complete"
    log:
        "logs/{sample}/masking/masking.log"
    benchmark:
        "benchmarks/{sample}/masking/masking.txt"
    params:
        sample = "{sample}",
        use_dev_shm = 1 if config['use_dev_shm'] else 0
    threads: int(config['slurm_args']['cpus_per_task'])
    resources:
        mem_mb = int(config['slurm_args']['mem_of_node']),
        runtime = int(config['slurm_args']['max_runtime'])
    container:
        "docker://dfam/tetools:latest"
    shell:
        r"""
        set -euo pipefail

        # Resolve all paths to absolute before cd-ing into workdir
        mkdir -p $(dirname {log})
        LOG_ABS=$(readlink -f {log})
        GENOME_ABS=$(readlink -f {input.genome})
        OUTDIR_ABS=$(readlink -f $(dirname {output.masked_genome}))
        mkdir -p "$OUTDIR_ABS"

        # Clean up stale masking workdirs from killed jobs (older than 3 days)
        find /dev/shm -maxdepth 1 -name "masking_*" -type d -mtime +3 -exec rm -rf {{}} \; 2>/dev/null || true
        find /tmp -maxdepth 1 -name "masking_*" -type d -mtime +3 -exec rm -rf {{}} \; 2>/dev/null || true

        # Choose working directory: /dev/shm for fast I/O or /tmp
        if [ "{params.use_dev_shm}" = "1" ] && [ -d /dev/shm ]; then
            WORKDIR=$(mktemp -d /dev/shm/masking_{params.sample}_XXXXXX)
            echo "[INFO] Using /dev/shm: $WORKDIR" > "$LOG_ABS"
        else
            WORKDIR=$(mktemp -d /tmp/masking_{params.sample}_XXXXXX)
            echo "[INFO] Using /tmp: $WORKDIR" > "$LOG_ABS"
        fi

        # Cleanup on exit
        trap "rm -rf $WORKDIR" EXIT

        cp "$GENOME_ABS" "$WORKDIR/genome.fa"
        cd "$WORKDIR"

        echo "[INFO] Running BuildDatabase..." >> "$LOG_ABS"
        BuildDatabase -name {params.sample} genome.fa >> "$LOG_ABS" 2>&1

        echo "[INFO] Running RepeatModeler2..." >> "$LOG_ABS"
        RepeatModeler -database {params.sample} -threads {threads} -LTRStruct >> "$LOG_ABS" 2>&1

        if [ ! -f "{params.sample}-families.fa" ]; then
            echo "[ERROR] RepeatModeler failed to produce {params.sample}-families.fa" >> "$LOG_ABS"
            exit 1
        fi
        echo "[INFO] RepeatModeler found $(grep -c '^>' {params.sample}-families.fa) repeat families" >> "$LOG_ABS"

        echo "[INFO] Running RepeatMasker..." >> "$LOG_ABS"
        RepeatMasker -pa {threads} -xsmall -lib {params.sample}-families.fa genome.fa >> "$LOG_ABS" 2>&1

        if [ ! -f "genome.fa.masked" ]; then
            echo "[ERROR] RepeatMasker failed to produce genome.fa.masked" >> "$LOG_ABS"
            exit 1
        fi

        # --- TRF tandem repeat masking (additional layer on top of RepeatMasker) ---
        echo "[INFO] Running TRF tandem repeat masking..." >> "$LOG_ABS"

        # Split masked genome into chunks (min 25 Mb each)
        perl {script_dir}/splitMfasta.pl --minsize=25000000 genome.fa.masked >> "$LOG_ABS" 2>&1

        # If no split files were produced (genome smaller than minsize), use whole genome as single chunk
        if ! ls genome.split.*.fa 1>/dev/null 2>&1; then
            echo "[INFO] No split files produced (genome smaller than 25 Mb), using whole genome" >> "$LOG_ABS"
            cp genome.fa.masked genome.split.0.fa
        fi

        # Run TRF on each chunk (sequential — trf is single-threaded per file)
        for chunk in genome.split.*.fa; do
            if [ -f "$chunk" ]; then
                echo "[INFO] TRF: $chunk" >> "$LOG_ABS"
                trf "$chunk" 2 7 7 80 10 50 500 -d -m -h >> "$LOG_ABS" 2>&1 || true
            fi
        done

        # Parse TRF output
        for dat in genome.split.*.fa.2.7.7.80.10.50.500.dat; do
            if [ -f "$dat" ]; then
                python3 {script_dir}/parseTrfOutput.py "$dat" --minCopies 1 \
                    --statistics "$dat.STATS" > "$dat.raw.gff" 2>> "$LOG_ABS" || true
            fi
        done

        # Sort parsed output
        for gff in genome.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff; do
            if [ -f "$gff" ]; then
                sort -k1,1 -k4,4n -k5,5n "$gff" > "$gff.sorted" 2>> "$LOG_ABS"
            fi
        done

        # Merge overlapping intervals and convert to GFF
        # Using perl instead of bedtools merge (not in tetools container)
        for sorted in genome.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff.sorted; do
            if [ -f "$sorted" ]; then
                perl -e '
                    my ($chr, $s, $e);
                    while (<>) {{
                        chomp;
                        my @f = split /\t/;
                        if (!defined($chr) || $f[0] ne $chr || $f[3] > $e) {{
                            print "$chr\ttrf\trepeat\t$s\t$e\t.\t.\t.\t.\n" if defined($chr);
                            ($chr, $s, $e) = ($f[0], $f[3], $f[4]);
                        }} else {{
                            $e = $f[4] if $f[4] > $e;
                        }}
                    }}
                    print "$chr\ttrf\trepeat\t$s\t$e\t.\t.\t.\t.\n" if defined($chr);
                ' "$sorted" > "$sorted.merged.gff" 2>> "$LOG_ABS"
            fi
        done

        # Soft-mask FASTA chunks with TRF regions
        for chunk in genome.split.*.fa; do
            if [ -f "$chunk" ] && [ -f "$chunk.2.7.7.80.10.50.500.dat.raw.gff.sorted.merged.gff" ]; then
                # Soft-mask using perl (replace matching regions with lowercase)
                perl -e '
                    use strict;
                    my %regions;
                    open(my $gff, "<", $ARGV[1]) or die "Cannot open GFF: $!";
                    while (<$gff>) {{
                        chomp;
                        my @f = split /\t/;
                        push @{{$regions{{$f[0]}}}}, [$f[3]-1, $f[4]]; # 0-based start
                    }}
                    close($gff);

                    open(my $fa, "<", $ARGV[0]) or die "Cannot open FASTA: $!";
                    my ($hdr, $seq);
                    while (<$fa>) {{
                        chomp;
                        if (/^>(.*)/) {{
                            if (defined $hdr) {{
                                foreach my $r (@{{$regions{{$hdr}} // []}}) {{
                                    substr($seq, $r->[0], $r->[1]-$r->[0]) =
                                        lc(substr($seq, $r->[0], $r->[1]-$r->[0]));
                                }}
                                print ">$hdr\n";
                                print substr($seq, $_, 80) . "\n" for grep {{ $_ < length($seq) }} map {{ $_ * 80 }} 0 .. int(length($seq)/80);
                            }}
                            $hdr = $1; $seq = "";
                        }} else {{ $seq .= $_; }}
                    }}
                    if (defined $hdr) {{
                        foreach my $r (@{{$regions{{$hdr}} // []}}) {{
                            substr($seq, $r->[0], $r->[1]-$r->[0]) =
                                lc(substr($seq, $r->[0], $r->[1]-$r->[0]));
                        }}
                        print ">$hdr\n";
                        print substr($seq, $_, 80) . "\n" for grep {{ $_ < length($seq) }} map {{ $_ * 80 }} 0 .. int(length($seq)/80);
                    }}
                    close($fa);
                ' "$chunk" "$chunk.2.7.7.80.10.50.500.dat.raw.gff.sorted.merged.gff" \
                    > "$chunk.combined.masked" 2>> "$LOG_ABS"
            else
                # No TRF regions for this chunk, keep as-is
                cp "$chunk" "$chunk.combined.masked" 2>/dev/null || true
            fi
        done

        # Concatenate all masked chunks into final output
        cat genome.split.*.fa.combined.masked > genome.fa.combined.masked 2>> "$LOG_ABS"

        N_TRF=$(cat genome.split.*.fa.2.7.7.80.10.50.500.dat.raw.gff.sorted.merged.gff 2>/dev/null | wc -l || echo 0)
        echo "[INFO] TRF found $N_TRF tandem repeat regions" >> "$LOG_ABS"

        # Use TRF+RepeatMasker combined output as final masked genome
        cp genome.fa.combined.masked "$OUTDIR_ABS/genome.fa.masked"
        cp {params.sample}-families.fa "$OUTDIR_ABS/{params.sample}-families.fa" 2>/dev/null || true

        # Preserve RepeatMasker report files (before workdir cleanup)
        cp genome.fa.out "$OUTDIR_ABS/repeatmasker.out" 2>/dev/null || true
        cp genome.fa.tbl "$OUTDIR_ABS/repeatmasker.tbl" 2>/dev/null || true

        touch "$OUTDIR_ABS/.masking_complete"
        echo "[INFO] Repeat masking complete for {params.sample}" >> "$LOG_ABS"

        # Record software versions
        VERSIONS_FILE="$(dirname "$OUTDIR_ABS")/software_versions.tsv"
        RM_VER=$(RepeatMasker -v 2>&1 | head -1 || true)
        RMOD_VER=$(RepeatModeler --version 2>&1 | head -1 || true)
        TRF_VER=$(trf 2>&1 | grep -oP 'Version \K\S+' || echo "4.09.1")
        ( flock 9
          printf "RepeatMasker\t%s\n" "$RM_VER" >> "$VERSIONS_FILE"
          printf "RepeatModeler\t%s\n" "$RMOD_VER" >> "$VERSIONS_FILE"
          printf "TRF\t%s\n" "$TRF_VER" >> "$VERSIONS_FILE"
        ) 9>"$VERSIONS_FILE.lock"

        # Report (still cd'd into $WORKDIR, use absolute paths)
        REPORT_DIR_ABS=$(dirname "$OUTDIR_ABS")
        mkdir -p "$REPORT_DIR_ABS"
        source {script_dir}/report_citations.sh
        N_FAMILIES=$(grep -c '^>' {params.sample}-families.fa 2>/dev/null || echo 0)
        cite repeatmodeler "$REPORT_DIR_ABS"
        cite repeatmasker "$REPORT_DIR_ABS"
        cite trf "$REPORT_DIR_ABS"
        """
