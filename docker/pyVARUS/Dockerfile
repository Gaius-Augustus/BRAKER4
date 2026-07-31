# pyVARUS: Python reimplementation of VARUS for automatic RNA-Seq selection
# from SRA and alignment to a genome for annotation support.
#
# Includes: pyVARUS, HISAT2, minimap2, samtools, sra-tools (fasterq-dump)
#
# Build:
#   docker build -t katharinahoff/pyvarus:v1.0.0 .
#
# Convert to Singularity:
#   singularity build pyvarus_v1.0.0.sif docker-daemon://katharinahoff/pyvarus:v1.0.0
#
# Usage (short-read mode):
#   varus runlist "Genus species" --outdir /out --email user@host
#   varus index genome.fa --outdir /out --threads 8
#   varus run "Genus species" genome.fa --runlist /out/runlist.tsv \
#       --index /out/genome_index --outdir /out --threads 8

FROM mambaorg/micromamba:1.5.8

USER root

# Install bioconda tools: hisat2, minimap2, samtools, sra-tools
# Pin samtools to >=1.17 for CSI index support (-c flag)
RUN micromamba install -y -n base \
        -c conda-forge -c bioconda \
        python=3.11 \
        hisat2 \
        minimap2 \
        "samtools>=1.17" \
        sra-tools \
        pip \
        git && \
    micromamba clean -afy

ENV PATH="/opt/conda/bin:$PATH"

# Clone and install pyVARUS with the [align] extra (provides pysam for
# intron extraction from BAM; requires htslib, Linux/macOS only)
RUN git clone https://github.com/Gaius-Augustus/pyVARUS.git /opt/pyvarus && \
    pip install --no-cache-dir -e "/opt/pyvarus[align]"

# Smoke-test: confirm the varus entry-point is on PATH
RUN varus --help

ENTRYPOINT ["varus"]
CMD ["--help"]
