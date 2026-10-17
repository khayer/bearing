#!/usr/bin/env bash
# End-to-end workflow test on simulated data (no real data needed).
#
#   bash workflow/test/run_sim_test.sh [WORKDIR] [CORES]
#
# 1. builds the simulated inputs (make_sim_inputs.sh),
# 2. runs the full BEARING Snakemake workflow on them (rule all; the Hi-C
#    and manuscript-figure phases are off because the sheet has no Hi-C and
#    paper_figures/paper_tables are false),
# 3. checks that the key outputs exist and are non-empty.
# Needs the conda env (environment.yml): snakemake, bedtools, pyGenomeTracks.
# Takes ~10 min on 4 cores. ASCII only.
set -euo pipefail

WORK="${1:-sim_workflow_test}"
CORES="${2:-4}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

for tool in snakemake bedtools pyGenomeTracks python; do
  command -v "$tool" > /dev/null || { echo "ERROR: '$tool' not on PATH (activate the bearing conda env)" >&2; exit 1; }
done

bash "$REPO/workflow/test/make_sim_inputs.sh" "$WORK"
WORK="$(cd "$WORK" && pwd)"

echo "== running workflow (log: $WORK/workflow.log)"
start=$(date +%s)
if ! snakemake -s "$REPO/workflow/Snakefile" --cores "$CORES" --keep-going \
     --configfile "$REPO/workflow/config/config.yaml" "$WORK/config.sim.yaml" \
     -d "$WORK" > "$WORK/workflow.log" 2>&1; then
  echo "FAILED. Failing rules:" >&2
  grep "Error in rule" "$WORK/workflow.log" | sort | uniq -c >&2 || true
  exit 1
fi
echo "   workflow finished in $(( $(date +%s) - start )) s"

R="$WORK/results"
fail=0
for f in \
  score.done compare.done pvalue.done region_qc.done calibration.done decomposition.done \
  A_rep1.qcat.bgz B_rep2.qcat.bgz \
  compare/diff_A_vs_B.qcat.bgz compare/pca.pdf compare/jsd_heatmap.pdf \
  pvalue/diff_A_vs_B.stats.tsv \
  regional/consolidated_enrichment_sim.tsv regional/significant_bins_summary.tsv \
  benchmark/recovery.done benchmarks_summary.tsv
do
  if [[ -s "$R/$f" || ( -e "$R/$f" && "$f" == *.done ) ]]; then
    echo "   ok      $f"
  else
    echo "   MISSING $f"; fail=1
  fi
done

nsig=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="pval_adj_bh") c=i; next} c && $c<0.05 {n++} END{print n+0}' \
       "$R/pvalue/diff_A_vs_B.stats.tsv")
echo "   A vs B: $nsig bins with FDR < 0.05 (expected > 0: B's planted blocks are shifted)"
[[ "$nsig" -gt 0 ]] || fail=1

if [[ $fail -eq 0 ]]; then echo "PASS: simulated end-to-end workflow"; else echo "FAIL" >&2; exit 1; fi
