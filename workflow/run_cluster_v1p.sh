#!/usr/bin/env bash
#
# run_cluster_v1p.sh -- preflight + launch the BEARING workflow for the
# V1P-vs-DN targeted validation run (6-track panel: RNA +/-, CTCF, Cohesin,
# NIPBL, H3K27ac), on SLURM.
#
# This is a SEPARATE run from both the core six-track pipeline and the 5-track
# V1-region run. It uses its own config (config_v1p.yaml -> samples_v1p.tsv,
# v1_6track_panel.yaml) and writes to its own outdir (results_v1p/), so it never
# collides with results_v1/ or the core DN results. RCTKO is dropped so H3K27ac
# can be included; scoring uses normalize=nonzero-quantile and min_signal=0.01
# (set in config_v1p.yaml) so the Trbv1 on->off bins survive and RNA/H3K27ac
# count comparably. Core-only run (no genome-wide Hi-C): p-values + region QC.
#
#   bash workflow/run_cluster_v1p.sh         # core path -> results_v1p/
#   bash workflow/run_cluster_v1p.sh --dry   # dry run only (no submission)
#
# Run from the directory your sheet's relative paths resolve against (the same
# place you would launch snakemake). Requires the `bearing` conda env active.
#
# ONE-TIME SETUP before first run:
#   1) cp workflow/config/categories/v1_5track_panel.yaml \
#         workflow/config/categories/v1_6track_panel.yaml
#      then add a 6th category, H3K27ac (e.g. color #008000), as the LAST entry
#      so it matches the bw order in samples_v1p.tsv.
#   2) Fill the H3K27ac bigwig paths in samples_v1p.tsv (placeholders for now).
# ASCII-only.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG="$REPO/workflow/config/config_v1p.yaml"
PROFILE="$REPO/workflow/profiles/slurm"
DRY=0
for a in "$@"; do
  case "$a" in
    --dry) DRY=1 ;;
    *) echo "unknown arg: $a" >&2; exit 2 ;;
  esac
done
OUT=$(python3 -c "import yaml;print(yaml.safe_load(open('$CONFIG'))['outdir'])")
echo "== V1P-vs-DN validation run (6-track incl H3K27ac) -> $OUT =="
echo
echo "== Preflight: checking every input file exists =="
python3 "$REPO/workflow/preflight.py" --configfile "$CONFIG" --core-only
echo
echo "== Snakemake dry run (validate DAG) =="
# core path only: p-values + per-region QC PNGs (no Hi-C targets)
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
echo; echo "Done -> workflow/$OUT/. Only differential is DN-vs-V1P."
echo "Check: awk -F'\t' '\$1==\"chr6\" && \$2>40891000 && \$2<40894000' $OUT/pvalue/diff_DN_vs_V1P.stats.tsv"
