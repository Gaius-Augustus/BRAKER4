# FANTASIA-Lite: lightweight protein functional annotation via embedding lookup.
# No PostgreSQL, no RabbitMQ -- just PyTorch + transformer models.
#
# Build:
#   docker build -t fantasia-lite .
#
# Convert to Singularity:
#   singularity build fantasia-lite.sif docker-daemon://fantasia-lite:latest

FROM nvidia/cuda:12.4.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

# System packages + Python with pip inside a venv (avoids PEP 668 issues)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 python3-pip python3-venv git && \
    rm -rf /var/lib/apt/lists/*

# Use a venv so pip works without --break-system-packages
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Clone FANTASIA-Lite
RUN git clone https://github.com/CBBIO/FANTASIA-Lite.git /opt/fantasia-lite

WORKDIR /opt/fantasia-lite

# Install runtime dependencies into the venv.
# Do NOT pip-install the fantasia-lite package itself -- its code uses
# FANTASIA_ROOT = Path(__file__).parent.parent which must resolve to
# /opt/fantasia-lite/, not site-packages/.
#
# PyTorch must be pinned to an explicit version that has +cu124 wheels.
# Without pinning, pip resolves to the latest PyTorch (2.11+) which only
# ships +cu130 wheels.  Installing a cu130 binary inside a cuda:12.4 base
# image causes SIGBUS (bus error) at the first CUDA kernel call that uses
# a CUDA 13.0 API not present in the 12.4 runtime.
# To update PyTorch: first check https://download.pytorch.org/whl/cu124/
# for the highest torch version with a cu124 wheel, then update the pin below.
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir \
        numpy requests biopython pandas scipy \
        sentencepiece protobuf h5py && \
    pip install --no-cache-dir \
        "torch==2.6.0+cu124" "torchvision==0.21.0+cu124" "torchaudio==2.6.0+cu124" \
        --extra-index-url https://download.pytorch.org/whl/cu124 && \
    pip install --no-cache-dir \
        transformers safetensors huggingface_hub requests

# Wrap pip so that "pip install --upgrade pip" is a silent no-op at runtime.
# All packages are pre-installed above; FANTASIA's install_packages() will see
# every dependency already satisfied and exit 0 without touching the filesystem.
# This makes the container work on read-only Singularity mounts (no --writable-tmpfs needed).
RUN mv /opt/venv/bin/pip /opt/venv/bin/_pip_real && \
    echo '#!/bin/bash' > /opt/venv/bin/pip && \
    echo 'if printf "%s\n" "$*" | grep -qE "install.*(--upgrade pip|pip --upgrade)"; then' >> /opt/venv/bin/pip && \
    echo '    echo "Requirement already satisfied: pip"; exit 0' >> /opt/venv/bin/pip && \
    echo 'fi' >> /opt/venv/bin/pip && \
    echo 'exec /opt/venv/bin/_pip_real "$@"' >> /opt/venv/bin/pip && \
    chmod +x /opt/venv/bin/pip

# Make the source scripts available on PATH (fantasia_pipeline, generate_embeddings, etc.)
ENV PATH="/opt/fantasia-lite/src:$PATH"
ENV PYTHONPATH="/opt/fantasia-lite/src:$PYTHONPATH"
RUN chmod +x /opt/fantasia-lite/src/*.py

# NOTE: lookup data (lookup_table.npz, annotations.json, accessions.json) is NOT
# baked into this image. Download it separately from Zenodo record 17720428 and
# provide the directory via bind-mount at runtime (configured in config.ini [FANTASIA] lookup_dir).

# Include the test FASTA so we can smoke-test without bind-mounts
COPY test.fa /opt/fantasia-lite/test.fa

ENTRYPOINT ["python3", "/opt/fantasia-lite/src/fantasia_pipeline.py"]
CMD ["--help"]
