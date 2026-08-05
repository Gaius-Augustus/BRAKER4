#!/bin/bash
# Build and push the fixed FANTASIA-Lite image.
# Run on greifserv3 (or any machine with Docker and docker login done).
#
# Fix: pin PyTorch to 2.6.0+cu124 so the binary matches the CUDA 12.4
# runtime in the base image (nvidia/cuda:12.4.0-runtime-ubuntu22.04).
# The old image had torch 2.11.0+cu130 installed into a CUDA 12.4 base,
# causing SIGBUS on T5 forward passes (CUDA 13.0 kernel APIs missing in 12.4).
#
# After pushing, on brain:
#   singularity pull fantasia_lite_v1.0.1.sif \
#       docker://katharinahoff/fantasia_for_brain:lite.v1.0.1
#   cp fantasia_lite_v1.0.1.sif /projects/BOUDICCA/db/images/fantasia_lite.sif

set -euo pipefail

IMAGE=katharinahoff/fantasia_for_brain
TAG=lite.v1.0.1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "[$(date)] Building $IMAGE:$TAG"
docker build -t "$IMAGE:$TAG" "$SCRIPT_DIR"

echo "[$(date)] Pushing $IMAGE:$TAG"
docker push "$IMAGE:$TAG"

echo "[$(date)] Done."
echo ""
echo "Next steps on brain:"
echo "  singularity pull /projects/BOUDICCA/db/images/fantasia_lite_v1.0.1.sif \\"
echo "      docker://katharinahoff/fantasia_for_brain:lite.v1.0.1"
echo "  # then either replace the existing SIF or update run_fantasia.smk"
