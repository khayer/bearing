#!/usr/bin/env bash
#
# run_cluster.sh -- preflight the inputs, then launch the BEARING workflow on
# SLURM (Snakemake submits each rule as its own sbatch job via the slurm
# profile). Fails fast if any input file is missing.
#
#   bash workflow/run_cluster.sh                 # full pipeline
#   bash workflow/run_cluster.sh --core-only     # skip Hi-C inputs/phases
#   bash workflow/run_cluster.sh --dry           # dry run only (no submission)
#
# Run from the directory your sheet's relative paths resolve against (the same
# place you would launch snakemake). Requires the `bearing` conda env active.
# ASCII-only.

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG="$REPO/workflow/config/config.yaml"
PROFILE="$REPO/workflow/profiles/slurm"

CORE_ONLY=0
DRY=0
for a in "$@"; do
  case "$a" in
    --core-only) CORE_ONLY=1 ;;
    --dry)       DRY=1 ;;
    *) echo "unknown arg: $a" >&2; exit 2 ;;
  esac
done

echo "== Preflight: checking every input file exists =="
PF_ARGS=(--configfile "$CONFIG")
[ "$CORE_ONLY" = "1" ] && PF_ARGS+=(--core-only)
python3 "$REPO/workflow/preflight.py" "${PF_ARGS[@]}"
echo

echo "== Snakemake dry run (validate DAG) =="
SM_TARGETS=()
if [ "$CORE_ONLY" = "1" ]; then
  # core path only: stop at the p-value sentinel, skip Hi-C targets
  OUT=$(python3 -c "import yaml;print(yaml.safe_load(open('$CONFIG'))['outdir'])")
  SM_TARGETS=("$REPO/workflow/$OUT/pvalue.done")
fi
snakemake -s "$REPO/workflow/Snakefile" --configfile "$CONFIG" -n "${SM_TARGETS[@]}"

if [ "$DRY" = "1" ]; then
  echo; echo "Dry run only (--dry). Not submitting."
  exit 0
fi

echo
echo "== Submitting to SLURM (one sbatch job per rule, via the slurm profile) =="
snakemake -s "$REPO/workflow/Snakefile" --configfile "$CONFIG" \
  --profile "$PROFILE" "${SM_TARGETS[@]}"
