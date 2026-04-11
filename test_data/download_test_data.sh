#!/bin/bash
# Download test data files needed for the BRAKER4 Snakemake pipeline.
#
# Usage:
#   bash test_data/download_test_data.sh             # download everything
#   bash test_data/download_test_data.sh --local-only # only local-scenario data
#
# The --local-only flag (or BRAKER4_LOCAL_ONLY=1 when sourcing) skips HPC-only
# datasets: the O. tauri genome, Viridiplantae proteins, OMAmer LUCA.h5 and
# Rfam. It still downloads RNAseq.bam, extracts paired FASTQ reads, and creates
# the isoseq_lib2.bam symlink -- everything the test_scenarios_local/ suite
# needs on top of the files already shipped in the repo.
#
# Downloads (full mode):
#   - O. tauri genome (~12.6 MB) for HPC test scenarios (08-11)
#   - Viridiplantae proteins for HPC protein scenarios (09, 11)
#   - RNAseq.bam (164 MB) for local A. thaliana scenarios (03, 04)
#   - Extracts FASTQ files from BAM for local FASTQ scenario (03)
#   - OMAmer LUCA.h5 (~15 GB) for HPC scenario 08
#   - Rfam covariance models (~30 MB) for run_ncrna=1

set -euo pipefail

TESTDATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

LOCAL_ONLY="${BRAKER4_LOCAL_ONLY:-0}"
for arg in "$@"; do
    case "$arg" in
        --local-only) LOCAL_ONLY=1 ;;
        -h|--help)
            sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            return 0 2>/dev/null || exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            return 1 2>/dev/null || exit 1
            ;;
    esac
done

if [ "$LOCAL_ONLY" = "1" ]; then
    echo "Running in --local-only mode: skipping HPC-only datasets."
fi

# ============================================================================
# Ostreococcus tauri genome (~12.6 MB, for HPC test scenarios 08-11)
# ============================================================================
if [ "$LOCAL_ONLY" != "1" ]; then
    OT_GENOME="$TESTDATA_DIR/Ostreococcus_tauri.fa"
    OT_GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz"

    if [ ! -f "$OT_GENOME" ]; then
        echo "Downloading O. tauri genome (GCF_000214015.3 from NCBI)..."
        echo "  URL: $OT_GENOME_URL"
        wget -q "$OT_GENOME_URL" -O "$OT_GENOME.gz"
        gunzip "$OT_GENOME.gz"
        echo "  Downloaded: $OT_GENOME ($(du -h "$OT_GENOME" | cut -f1))"
    else
        echo "O. tauri genome already exists: $OT_GENOME"
    fi
fi

# ============================================================================
# Viridiplantae proteins (for HPC scenarios with protein evidence)
# ============================================================================
if [ "$LOCAL_ONLY" != "1" ]; then
    FULL_PROTEINS="$TESTDATA_DIR/Viridiplantae.fa"
    FULL_PROTEINS_URL="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Viridiplantae.fa.gz"

    if [ ! -f "$FULL_PROTEINS" ]; then
        echo "Downloading Viridiplantae protein set (ODB12)..."
        echo "  URL: $FULL_PROTEINS_URL"
        wget -q "$FULL_PROTEINS_URL" -O "$FULL_PROTEINS.gz"
        gunzip "$FULL_PROTEINS.gz"
        echo "  Downloaded: $FULL_PROTEINS ($(du -h "$FULL_PROTEINS" | cut -f1))"
    else
        echo "Viridiplantae proteins already exist: $FULL_PROTEINS"
    fi
fi

# ============================================================================
# RNA-Seq BAM (for local BAM/FASTQ scenarios)
# ============================================================================
BAM_FILE="$TESTDATA_DIR/RNAseq.bam"
BAM_URL="http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam"

if [ ! -f "$BAM_FILE" ]; then
    echo "Downloading test BAM file (164 MB)..."
    echo "  URL: $BAM_URL"
    wget -q "$BAM_URL" -O "$BAM_FILE"
    echo "  Download complete."
else
    echo "Test BAM file already exists: $BAM_FILE"
fi

# ============================================================================
# Extract paired FASTQ files from BAM (for local FASTQ scenario)
# ============================================================================
FASTQ_R1="$TESTDATA_DIR/reads_R1.fastq.gz"
FASTQ_R2="$TESTDATA_DIR/reads_R2.fastq.gz"

if [ ! -f "$FASTQ_R1" ] || [ ! -f "$FASTQ_R2" ]; then
    echo "Extracting paired FASTQ files from BAM..."
    samtools sort -n -@ 4 "$BAM_FILE" | \
        samtools fastq -@ 4 -1 "$FASTQ_R1" -2 "$FASTQ_R2" -0 /dev/null -s /dev/null -
    echo "  Created: $FASTQ_R1 ($(du -h "$FASTQ_R1" | cut -f1))"
    echo "  Created: $FASTQ_R2 ($(du -h "$FASTQ_R2" | cut -f1))"
else
    echo "Test FASTQ files already exist: $FASTQ_R1, $FASTQ_R2"
fi

# ============================================================================
# OMAmer database for OMArk (large, ~15 GB; only needed by scenario 08)
# ============================================================================
if [ "$LOCAL_ONLY" != "1" ]; then
    OMAMER_DB="$TESTDATA_DIR/LUCA.h5"

    if [ ! -f "$OMAMER_DB" ]; then
        echo "Downloading OMAmer LUCA.h5 database (~15 GB)..."
        echo "  URL: https://omabrowser.org/All/LUCA.h5"
        wget -q "https://omabrowser.org/All/LUCA.h5" -O "$OMAMER_DB"
        echo "  Downloaded: $OMAMER_DB ($(du -h "$OMAMER_DB" | cut -f1))"
    else
        echo "OMAmer database already exists: $OMAMER_DB"
    fi
fi

# ============================================================================
# Rfam database for Infernal ncRNA annotation (run_ncrna=1)
# ============================================================================
if [ "$LOCAL_ONLY" != "1" ]; then
    RFAM_DIR="$(cd "$TESTDATA_DIR/.." && pwd)/shared_data/rfam"
    RFAM_CM="$RFAM_DIR/Rfam.cm"
    RFAM_CLANIN="$RFAM_DIR/Rfam.clanin"

    if [ ! -f "$RFAM_CM" ]; then
        mkdir -p "$RFAM_DIR"

        if [ ! -f "$RFAM_CM" ]; then
            echo "Downloading Rfam covariance models (~30 MB)..."
            wget -q "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz" -O "$RFAM_CM.gz"
            gunzip "$RFAM_CM.gz"
            echo "  Downloaded: $RFAM_CM ($(du -h "$RFAM_CM" | cut -f1))"
        fi

        if [ ! -f "$RFAM_CLANIN" ]; then
            echo "Downloading Rfam clan info..."
            wget -q "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin" -O "$RFAM_CLANIN"
            echo "  Downloaded: $RFAM_CLANIN"
        fi

        # Note: cmpress indexing is handled by the Snakemake Infernal rule
        # (inside the container where cmpress is available)
        echo "  Rfam database downloaded. cmpress will run inside the Infernal container at pipeline time."
    else
        echo "Rfam database already exists: $RFAM_DIR/"
    fi
fi

# ============================================================================
# FANTASIA-Lite container + ProtT5 HuggingFace cache (optional, GPU-only)
# ============================================================================
# This block is intentionally OPT-IN -- the SIF is ~6 GB and the ProtT5 weights
# are another ~5 GB, and FANTASIA-Lite is the most fragile, GPU-only component
# of BRAKER4 (validated only on an A100 in the Hoff lab). Set
# BRAKER4_DOWNLOAD_FANTASIA=1 to stage both, otherwise this section is skipped.
#
# After staging, point the [fantasia] section of config.ini at:
#     sif          = <repo>/shared_data/fantasia/fantasia_lite.sif
#     hf_cache_dir = <repo>/shared_data/fantasia/hf_cache
# (these are the BRAKER4 defaults if both keys are left unset).
#
# The Snakefile checks the SIF exists at parse time when fantasia.enable=1, so
# you will get an immediate error -- not a wasted GPU job -- if the staging
# step has not run.
if [ "${BRAKER4_DOWNLOAD_FANTASIA:-0}" = "1" ]; then
    FANTASIA_DIR="$(cd "$TESTDATA_DIR/.." && pwd)/shared_data/fantasia"
    FANTASIA_SIF="$FANTASIA_DIR/fantasia_lite.sif"
    FANTASIA_HF_CACHE="$FANTASIA_DIR/hf_cache"
    FANTASIA_IMAGE="docker://katharinahoff/fantasia_for_brain:lite.v0.0.2"

    mkdir -p "$FANTASIA_DIR" "$FANTASIA_HF_CACHE"

    # Most HPC sites need `module load singularity`, and on systems without
    # environment modules the command silently no-ops. Mirror the pattern from
    # EukAssembly-Bin/rules/fantasia.smk so this script behaves the same on
    # both kinds of host.
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load singularity 2>/dev/null || true

    if ! command -v singularity >/dev/null 2>&1; then
        echo "ERROR: singularity not on PATH after attempting 'module load singularity'." >&2
        echo "       FANTASIA-Lite requires Singularity/Apptainer with --nv GPU support." >&2
        echo "       Install singularity (or load the appropriate module) and re-run with" >&2
        echo "       BRAKER4_DOWNLOAD_FANTASIA=1." >&2
        exit 1
    fi

    if [ ! -f "$FANTASIA_SIF" ]; then
        echo "Pulling FANTASIA-Lite container (~6 GB)..."
        echo "  Image: $FANTASIA_IMAGE"
        echo "  Target: $FANTASIA_SIF"
        singularity pull "$FANTASIA_SIF" "$FANTASIA_IMAGE"
        echo "  Pulled: $FANTASIA_SIF ($(du -h "$FANTASIA_SIF" | cut -f1))"
    else
        echo "FANTASIA-Lite container already exists: $FANTASIA_SIF"
    fi

    # Pre-cache the ProtT5 weights inside the SIF so HF_HUB_OFFLINE=1 works at
    # pipeline runtime. The download is ~5 GB and only happens once per host.
    if [ ! -d "$FANTASIA_HF_CACHE/hub" ] && [ -z "$(ls -A "$FANTASIA_HF_CACHE" 2>/dev/null)" ]; then
        echo "Pre-caching ProtT5 weights (Rostlab/prot_t5_xl_uniref50, ~5 GB)..."
        echo "  Target: $FANTASIA_HF_CACHE"
        HF_HOME="$FANTASIA_HF_CACHE" singularity exec \
            -B "$FANTASIA_HF_CACHE":"$FANTASIA_HF_CACHE" \
            "$FANTASIA_SIF" \
            python3 -c "
from transformers import T5Tokenizer, T5EncoderModel
T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50')
T5EncoderModel.from_pretrained('Rostlab/prot_t5_xl_uniref50')
print('ProtT5 weights cached.')
"
        echo "  Cache populated: $(du -sh "$FANTASIA_HF_CACHE" | cut -f1)"
    else
        echo "ProtT5 HuggingFace cache already populated: $FANTASIA_HF_CACHE"
    fi

    echo ""
    echo "FANTASIA-Lite staging complete. Add to your config.ini:"
    echo "    [fantasia]"
    echo "    enable       = 1"
    echo "    sif          = $FANTASIA_SIF"
    echo "    hf_cache_dir = $FANTASIA_HF_CACHE"
else
    echo "Skipping FANTASIA-Lite staging (set BRAKER4_DOWNLOAD_FANTASIA=1 to enable)."
fi

# ============================================================================
# Create symlink for multi-IsoSeq BAM testing (scenario 05)
# ============================================================================
if [ ! -e "$TESTDATA_DIR/isoseq_lib2.bam" ] && [ -f "$TESTDATA_DIR/isoseq.bam" ]; then
    ln -sf isoseq.bam "$TESTDATA_DIR/isoseq_lib2.bam"
    echo "Created symlink: isoseq_lib2.bam -> isoseq.bam"
fi

echo ""
echo "Test data ready."
