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

if command -v apptainer >/dev/null 2>&1; then
	RUNTIME="apptainer"
elif command -v singularity >/dev/null 2>&1; then
	RUNTIME="singularity"
else
	cat >&2 <<'EOF'
ERROR: neither apptainer nor singularity is installed on PATH.

To pull the prebuilt image, install Apptainer/Singularity first, or build a
local SIF with:
	apptainer build bearing.sif Apptainer.def
EOF
	exit 127
fi

echo "Pulling oras://ghcr.io/${OWNER}/bearing:${TAG} -> ${OUT}"
"${RUNTIME}" pull --force "${OUT}" "oras://ghcr.io/${OWNER}/bearing:${TAG}"
echo "Done. Use it with: snakemake --use-apptainer ..."
