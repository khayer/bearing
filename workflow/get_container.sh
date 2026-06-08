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

# Detect container runtime: apptainer or singularity (the renamed predecessor).
if command -v apptainer >/dev/null 2>&1; then
  RT=apptainer
elif command -v singularity >/dev/null 2>&1; then
  RT=singularity
else
  echo "ERROR: neither 'apptainer' nor 'singularity' is on PATH." >&2
  echo "Install one of them, or build locally with: apptainer build bearing.sif Apptainer.def" >&2
  exit 1
fi

echo "Pulling oras://ghcr.io/${OWNER}/bearing:${TAG} -> ${OUT}  (using ${RT})"
"$RT" pull --force "${OUT}" "oras://ghcr.io/${OWNER}/bearing:${TAG}"
echo "Done. Use it with: snakemake --use-apptainer ..."
