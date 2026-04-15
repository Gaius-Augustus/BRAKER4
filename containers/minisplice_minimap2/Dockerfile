FROM ubuntu:22.04

LABEL maintainer="Katharina Hoff <katharina.hoff@uni-greifswald.de>"
LABEL description="minimap2 2.29+ with minisplice splice-site scoring for BRAKER4"

ARG MINIMAP2_VERSION=2.29
ARG MINISPLICE_VERSION=0.4

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        wget \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install minimap2
RUN wget -q https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2 \
    && tar xjf minimap2-${MINIMAP2_VERSION}.tar.bz2 \
    && cd minimap2-${MINIMAP2_VERSION} \
    && make -j$(nproc) \
    && cp minimap2 /usr/local/bin/ \
    && cd .. \
    && rm -rf minimap2-${MINIMAP2_VERSION} minimap2-${MINIMAP2_VERSION}.tar.bz2

# Install minisplice binary from git
RUN git clone --depth 1 --branch v${MINISPLICE_VERSION} https://github.com/lh3/minisplice.git \
    && cd minisplice \
    && make -j$(nproc) \
    && cp minisplice /usr/local/bin/ \
    && cd .. \
    && rm -rf minisplice

# Download pre-trained model and calibration from Zenodo (not in the git repo)
RUN mkdir -p /usr/local/share/minisplice \
    && wget -O- https://zenodo.org/records/15814006/files/vi2-7k.tgz | tar zxf - -C /usr/local/share/minisplice

ENV MINISPLICE_MODEL=/usr/local/share/minisplice/vi2-7k.kan

# Install samtools (for BAM sorting in the sort_isoseq_sam rule)
ARG SAMTOOLS_VERSION=1.21
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --without-curses \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

# Verify installs
RUN minimap2 --version \
    && minisplice version 2>/dev/null || minisplice 2>&1 | head -1 \
    && ls /usr/local/share/minisplice/

WORKDIR /data
