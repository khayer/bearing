#!/usr/bin/env bash
# Preflight: verify the scripts in THIS checkout match sync_manifest.sha256 (the
# versions that shipped together). Catches the recurring "synced the Snakefile
# but not the .py" partial copy, and any stale cluster file.
#
# Run on the cluster from the repo root, before launching anything:
#     bash sync_check.sh
#
# Exit 0 = all in sync; 1 = stale/missing files (listed with copy commands);
# 2 = manifest absent. ASCII only.
set -uo pipefail
cd "$(dirname "$0")"
M=sync_manifest.sha256
if [ ! -f "$M" ]; then
  echo "[sync_check] no $M here -- copy it from the source repo (it ships in the zip)."
  exit 2
fi
n=$(wc -l < "$M")
out=$(sha256sum -c "$M" 2>/dev/null) || true
bad=$(printf '%s\n' "$out" | grep -vE ': OK$' || true)
if [ -z "$bad" ]; then
  echo "[sync_check] OK -- all $n tracked files match the shipped version."
  exit 0
fi
echo "[sync_check] STALE or MISSING (cluster copy differs from shipped version):"
printf '%s\n' "$bad" | sed -E 's/: (FAILED.*)?$/   <-- \1/' | sed 's/^/  /'
echo
echo "[sync_check] copy these from the source repo, then re-run sync_check.sh:"
printf '%s\n' "$bad" | sed -E 's/:.*$//' | sed 's/^/  cp /'
exit 1
