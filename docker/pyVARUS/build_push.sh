#!/bin/bash
# Build and push the pyVARUS image.
# Run on a machine with Docker and docker login done.
#
# After pushing, pull to brain with:
#   singularity pull pyvarus_v1.0.0.sif docker://katharinahoff/pyvarus:v1.0.0
#
# Then update varus_image in Snakefile _container_defaults to:
#   'varus_image': 'docker://katharinahoff/pyvarus:v1.0.0'

set -euo pipefail

IMAGE=katharinahoff/pyvarus
TAG=v1.0.0

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "[$(date)] Building $IMAGE:$TAG"
docker build -t "$IMAGE:$TAG" "$SCRIPT_DIR"

echo "[$(date)] Pushing $IMAGE:$TAG"
docker push "$IMAGE:$TAG"

echo "[$(date)] Done."
echo ""
echo "Next steps on brain:"
echo "  singularity pull pyvarus_v1.0.0.sif docker://$IMAGE:$TAG"
