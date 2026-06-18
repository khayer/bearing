#!/usr/bin/env bash
# Regenerate sync_manifest.sha256 over the tracked scripts (root *.py, the
# Snakefile, launchers, and workflow wrappers). Run from the source repo AFTER
# any code change and BEFORE shipping, so the manifest matches what gets synced.
# ASCII only.
set -euo pipefail
cd "$(dirname "$0")"
files=$( {
  ls *.py 2>/dev/null || true
  ls *.R 2>/dev/null || true
  ls *.sh 2>/dev/null || true
  ls workflow/Snakefile 2>/dev/null || true
  ls workflow/*.sh 2>/dev/null || true
  ls workflow/*.py 2>/dev/null || true
  ls slurm/*.sh 2>/dev/null || true
  ls slurm/*.py 2>/dev/null || true
} | sort -u )
[ -n "$files" ] || { echo "no tracked files found"; exit 1; }
sha256sum $files > sync_manifest.sha256
echo "wrote sync_manifest.sha256 ($(wc -l < sync_manifest.sha256) files)"
