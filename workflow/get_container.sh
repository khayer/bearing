#!/usr/bin/env bash
# Pull the prebuilt BEARING container from GHCR into the working directory.
#
#   bash workflow/get_container.sh            # pulls :latest
#   bash workflow/get_container.sh v1.0.0     # pulls a specific tag
#
# Set BEARING_OWNER if your fork is under a different account.
set -euo pipefail

OWNER="${BEARING_OWNER:-khayer}"
TAG="${1:-latest}"
OUT="${BEARING_SIF:-bearing.sif}"

echo "Pulling oras://ghcr.io/${OWNER}/bearing:${TAG} -> ${OUT}"
apptainer pull --force "${OUT}" "oras://ghcr.io/${OWNER}/bearing:${TAG}"
echo "Done. Use it with: snakemake --use-apptainer ..."
