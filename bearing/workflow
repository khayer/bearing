#!/usr/bin/env bash
#
# run_cluster_v1.sh -- preflight + launch the BEARING workflow for the V1-region
# use case (Allyn 2025 mutants: DN/R1KO, V1P, RCTKO), 5-track panel, on SLURM.
#
# This is a SEPARATE run from the core six-track pipeline: it uses its own config
# (config_v1.yaml -> samples_v1.tsv, v1_5track_panel.yaml) and writes to its own
# outdir (results_v1/), so it never collides with the core "DN" results. The V1
# conditions have no genome-wide Hi-C (their 3D side is capture, handled
# separately), so this is a core-only run: per-condition p-values + region QC.
#
#   bash workflow/run_cluster_v1.sh         # core path -> results_v1/
#   bash workflow/run_cluster_v1.sh --dry   # dry run only (no submission)
#
# Run from the directory your sheet's relative paths resolve against (the same
# place you would launch snakemake). Requires the `bearing` conda env active.
# ASCII-only.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG="$REPO/workflow/config/config_v1.yaml"
PROFILE="$REPO/workflow/profiles/slurm"
DRY=0
for a in "$@"; do
  case "$a" in
    --dry) DRY=1 ;;
    *) echo "unknown arg: $a" >&2; exit 2 ;;
  esac
done

OUT=$(python3 -c "import yaml;print(yaml.safe_load(open('$CONFIG'))['outdir'])")
echo "== V1-region run (5-track, core-only) -> $OUT =="
echo

echo "== Preflight: checking every input file exists =="
python3 "$REPO/workflow/preflight.py" --configfile "$CONFIG" --core-only
echo

echo "== Snakemake dry run (validate DAG) =="
# core path only: p-values + per-region QC PNGs (no Hi-C targets for V1)
SM_TARGETS=("$REPO/workflow/$OUT/pvalue.done"
            "$REPO/workflow/$OUT/region_qc.done")
snakemake -s "$REPO/workflow/Snakefile" --configfile "$CONFIG" -n "${SM_TARGETS[@]}"
if [ "$DRY" = "1" ]; then
  echo; echo "Dry run only (--dry). Not submitting."
  exit 0
fi

echo
echo "== Submitting to SLURM (one sbatch job per rule, via the slurm profile) =="
snakemake -s "$REPO/workflow/Snakefile" --configfile "$CONFIG" \
  --profile "$PROFILE" "${SM_TARGETS[@]}"

echo; echo "Done -> workflow/$OUT/. Differentials are DN-vs-V1P and DN-vs-RCTKO."
