#!/bin/bash
# Download a pre-partitioned OrthoDB v12 protein FASTA for BRAKER4.
#
# Usage:
#   bash scripts/download_orthodb.sh <clade> [output_dir]
#
# Arguments:
#   clade       Clade name (case-insensitive). One of:
#               Metazoa, Vertebrata, Viridiplantae, Arthropoda, Fungi,
#               Alveolata, Stramenopiles, Amoebozoa, Euglenozoa, Eukaryota
#   output_dir  Directory to save the file (default: current directory)
#
# Examples:
#   bash scripts/download_orthodb.sh Viridiplantae
#   bash scripts/download_orthodb.sh Metazoa /data/orthodb
#
# Notes:
#   - The file is downloaded from the partitioned_odb12 mirror at
#     https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/
#   - Partitions are large (100s of MB to a few GB compressed).
#   - Stramenopiles, Amoebozoa, and Euglenozoa are smaller clades; for
#     these, BRAKER's authors recommend either combining with Eukaryota
#     or using Eukaryota directly.

set -euo pipefail

VALID_CLADES=(
    Metazoa Vertebrata Viridiplantae Arthropoda Fungi
    Alveolata Stramenopiles Amoebozoa Euglenozoa Eukaryota
)

usage() {
    echo "Usage: $0 <clade> [output_dir]"
    echo ""
    echo "Available clades:"
    printf "  %s\n" "${VALID_CLADES[@]}"
    exit 1
}

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    usage
fi

# Normalize input: first letter upper, rest lower (matches mirror filenames)
INPUT_CLADE="$1"
CLADE="$(echo "$INPUT_CLADE" | tr '[:upper:]' '[:lower:]')"
CLADE="$(tr '[:lower:]' '[:upper:]' <<< "${CLADE:0:1}")${CLADE:1}"

# Validate
VALID=false
for c in "${VALID_CLADES[@]}"; do
    if [ "$c" = "$CLADE" ]; then
        VALID=true
        break
    fi
done

if [ "$VALID" = false ]; then
    echo "ERROR: '$INPUT_CLADE' is not a recognized clade." >&2
    echo ""
    usage
fi

OUTDIR="${2:-.}"
mkdir -p "$OUTDIR"
OUT_FA="$OUTDIR/$CLADE.fa"
OUT_GZ="$OUT_FA.gz"
URL="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/$CLADE.fa.gz"

if [ -f "$OUT_FA" ]; then
    echo "OrthoDB file already exists: $OUT_FA"
    echo "(Delete it to force re-download.)"
    exit 0
fi

echo "Downloading OrthoDB v12 partition: $CLADE"
echo "  From: $URL"
echo "  To:   $OUT_GZ"

if command -v wget &> /dev/null; then
    wget -q --show-progress "$URL" -O "$OUT_GZ"
elif command -v curl &> /dev/null; then
    curl -# -L "$URL" -o "$OUT_GZ"
else
    echo "ERROR: neither wget nor curl is available." >&2
    exit 1
fi

echo "Decompressing..."
gunzip "$OUT_GZ"

echo ""
echo "Done. OrthoDB partition ready at:"
echo "  $OUT_FA"
echo ""
echo "Use this path in the 'protein_fasta' column of your samples.csv."
