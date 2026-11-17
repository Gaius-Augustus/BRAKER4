#!/bin/bash
# Data for scenario 21 (reproduces BRAKER4 issue #110: GeneMark-ETP
# fails in the masking-penalty step on S. cerevisiae with fungus = 1).
#   - S. cerevisiae S288C genome, RefSeq GCF_000146045.2 (R64), unmasked
#   - OrthoDB v12 Fungi partition (as used in the issue)
# RNA-Seq (SRR5819434, SRR5819435) is fetched by the pipeline via sra_ids.
set -euo pipefail
TESTDATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../test_data" && pwd)"
GENOME="$TESTDATA_DIR/Saccharomyces_cerevisiae_R64.fa"
PROT="$TESTDATA_DIR/Fungi.fa"
if [ ! -s "$GENOME" ]; then
    echo "Downloading S. cerevisiae R64 genome..."
    curl -fsSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz \
        | gunzip -c > "$GENOME.tmp"
    mv "$GENOME.tmp" "$GENOME"
fi
if [ ! -s "$PROT" ]; then
    echo "Downloading OrthoDB v12 Fungi partition..."
    curl -fsSL https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Fungi.fa.gz \
        | gunzip -c > "$PROT.tmp"
    mv "$PROT.tmp" "$PROT"
fi
