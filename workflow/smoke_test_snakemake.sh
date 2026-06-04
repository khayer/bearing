#!/usr/bin/env bash
#
# smoke_test_snakemake.sh -- runs the REAL Snakemake workflow end to end on a
# tiny synthetic dataset, through the core (no-Hi-C) path. Unlike smoke_test.sh
# (which calls the scripts directly), this exercises the Snakefile itself: rule
# wiring, .done sentinel ordering, the BEARING-sheet derivation, wildcard
# expansion, and config plumbing. It is hermetic -- it synthesizes its own
# reference files and needs no real data, no Hi-C, no cluster.
#
#   bash workflow/smoke_test_snakemake.sh
#
# Exit 0 = `snakemake` drove score -> compare -> perm -> p-values to completion
# and the expected sentinels/outputs exist.
# ASCII-only.

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SANDBOX="$(mktemp -d)"
echo "repo    : $REPO"
echo "sandbox : $SANDBOX"

CHROM=chrSim
LEN=2000000
CONF="$SANDBOX/config"
RES="$SANDBOX/resources"
mkdir -p "$CONF" "$RES" "$SANDBOX/annotations"

echo "== 1. synthesize 4 samples + minimal reference files =="
cd "$SANDBOX"
for s in A_rep1 A_rep2 B_rep1 B_rep2; do
  python3 "$REPO/benchmark/simulate_bearing_tracks.py" \
    --prefix "$s" --chrom "$CHROM" --length "$LEN" --bin-size 200 \
    --seed "${#s}" >/dev/null
done
printf "%s\t%d\n" "$CHROM" "$LEN" > "$RES/sim.chrom.sizes"
# minimal stand-ins for the reference inputs the core rules name
: > "$RES/empty_blacklist.bed"
: > "$RES/genes.bed"
: > "$RES/agr.bed"
printf "name\tregion\tresolution\tlabel\nnew_wide\t%s:0-%d\t2000\tsim\n" "$CHROM" "$LEN" > "$CONF/regions.tsv"

echo "== 2. write sample sheet (absolute bw paths) =="
{
  printf "sample\tcondition\treplicate\tout\tbw\n"
  for s in A_rep1 A_rep2 B_rep1 B_rep2; do
    cond="${s%%_*}"; rep="${s##*rep}"
    bws="$SANDBOX/${s}_atac.bw,$SANDBOX/${s}_rnaplus.bw,$SANDBOX/${s}_rnaminus.bw,$SANDBOX/${s}_ctcf.bw,$SANDBOX/${s}_rad21.bw,$SANDBOX/${s}_h3k27ac.bw"
    printf "%s\t%s\t%s\t%s.qcat.bgz\t%s\n" "$s" "$cond" "$rep" "$s" "$bws"
  done
} > "$CONF/samples.tsv"

echo "== 3. write a minimal config (no Hi-C, synthetic refs, abs paths) =="
cat > "$CONF/config.yaml" << EOF
samples_sheet: "$CONF/samples.tsv"
categories:    "$REPO/categories/mm10_6track_panel.yaml"
chrom_sizes:   "$RES/sim.chrom.sizes"
regions_file:  "$CONF/regions.tsv"
genes_bed:     "$RES/genes.bed"
agr_genes_bed: "$RES/agr.bed"
gtf:           ""
blacklist_external: "$RES/empty_blacklist.bed"
min_signal: 0.1
normalize: ""
n_perms: 2
min_shift: 500000
seed: 1
blacklist_merge_distance: 3003
reference_condition: A
hic:
  resolutions: [10000]
  contact_resolution: 10000
  tad_dir: ""
  tad_anchor: ""
  features_bed: ""
  compartment_region: ""
  isolation_target: ""
  isolation_window: 0
  isolation_max_distance: 0
  isolation_perm_n: 0
  crosslocus_targets: ""
outdir: "$SANDBOX/results"
bearing_dir: "$REPO"
container: "bearing.sif"
EOF

echo "== 4. snakemake DRY RUN (validates the real DAG) =="
snakemake -s "$REPO/workflow/Snakefile" --configfile "$CONF/config.yaml" -n >/dev/null
echo "  dry run OK"

echo "== 5. snakemake EXECUTE core targets (real rules, 2 cores) =="
# Target the core sentinels explicitly: this drives bearing_sheet -> blacklist
# -> score -> compare -> perm nulls -> diff p-values through the actual rules.
snakemake -s "$REPO/workflow/Snakefile" --configfile "$CONF/config.yaml" \
  --cores 2 \
  "$SANDBOX/results/compare.done" \
  "$SANDBOX/results/pvalue.done" \
  >/dev/null 2>"$SANDBOX/snakemake.err" \
  || { echo "  FAIL: snakemake execution errored"; tail -20 "$SANDBOX/snakemake.err"; exit 1; }

echo "== 6. verify outputs produced by the workflow =="
for f in \
  "$SANDBOX/results/samples.bearing.tsv" \
  "$SANDBOX/results/blacklist/merged_blacklist.bed" \
  "$SANDBOX/results/score.done" \
  "$SANDBOX/results/compare.done" \
  "$SANDBOX/results/pvalue.done"; do
  test -e "$f" || { echo "  FAIL: missing $f"; exit 1; }
  echo "  ok: ${f#$SANDBOX/}"
done

echo
echo "=========================================================="
echo "PASS: the Snakemake workflow ran the core path end to end."
echo "  bearing_sheet -> blacklist -> score -> compare -> perm -> p-values"
echo "  (driven through the real Snakefile, not the scripts directly)"
echo "results: $SANDBOX/results"
echo "=========================================================="
