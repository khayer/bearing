#!/usr/bin/env bash
# ----------------------------------------------------------------------
# File     : paper/run_edger_validation.sh
# Version  : 1.1.0
# Date     : 2026-10-12 (1.1.0: optional DN vs S3T3 arm)
# Authors  : Katharina E. Hayer (katharinaehayer@gmail.com) and Claude
#            (Anthropic), co-created
# ----------------------------------------------------------------------
# One command for the independent count-based (edgeR) validation of BEARING's
# strand-resolved RNA differentials: Table S8 and the Results sentences that
# quote it. Runs, against one finished BEARING results directory:
#
#   1. DN vs DP, genome-wide, edgeR filterByExpr min.count = 10   (Table S8, headline)
#   2. DN vs DP, genome-wide, min.count = 3                         (Table S8, sensitivity)
#   3. DN vs Pro-B, concordance restricted to the Igh locus, 10 and 3
#   4. DN vs EbKO, recombination-center characterisation (region_match_ebko.py)
#   5. OPTIONAL, only when the manifest has S3T3 rows: DN vs 3T3, genome-wide,
#      min.count = 10 and 3 (Col1a1 and recombination-center calls)
#
# For each, featureCounts counts reads per strand on BEARING's own 200 bp grid
# (blacklist from the same run), edgeR quasi-likelihood tests each bin, and the
# BEARING RNA+ / RNA- differential bins of that comparison are scored for
# recovery on the TESTABLE subset, direction agreement and rank enrichment
# (rna_concordance_stranded.R -> dev/extract_track_diff_bed.py ->
# dev/bearing_diff_concordance.py). featureCounts/edgeR results are cached per
# comparison and strand, so the min.count = 3 and region runs reuse them.
#
# The two Table S8 summaries are copied to <results>/table_sources/ under the
# names paper/build_tables.py expects, so rebuilding the workbook picks them up.
#
# Requirements: R with Rsubread, edgeR, GenomicRanges, rtracklayer; python3 with
# pandas and numpy (the BEARING environment); the RNA-seq BAMs.
#
# Usage (from the repository root):
#   bash paper/run_edger_validation.sh \
#       --bams paper/edger_bam_manifest.tsv \
#       [--results workflow/results_prior1_n100] [--library reverse] \
#       [--single-end] [--no-dup-flag] [--threads 4] [--outdir DIR] [--dry]
#
# The BAM manifest is a TSV with a header and columns: sample  condition  bam
# (condition one of DN, DP, EbKO, ProB, and optionally S3T3). Template:
# paper/edger_bam_manifest.template.tsv. Library type (--library, --single-end,
# --no-dup-flag) applies to every BAM; check that the S3T3 BAMs share it.
#
# ASCII-only.

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS="$REPO/workflow/results_prior1_n100"
MANIFEST=""
LIBRARY="reverse"
PAIRED="TRUE"
DUPS="TRUE"
THREADS=4
OUTDIR=""
DRY=0

while [ $# -gt 0 ]; do
  case "$1" in
    --bams)        MANIFEST="$2"; shift 2 ;;
    --results)     RESULTS="$2"; shift 2 ;;
    --library)     LIBRARY="$2"; shift 2 ;;
    --single-end)  PAIRED="FALSE"; shift ;;
    --no-dup-flag) DUPS="FALSE"; shift ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --outdir)      OUTDIR="$2"; shift 2 ;;
    --dry)         DRY=1; shift ;;
    -h|--help)     sed -n 2,40p "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

RESULTS="$(cd "$RESULTS" && pwd)"
[ -z "$OUTDIR" ] && OUTDIR="$RESULTS/edger_validation"
TS="$RESULTS/table_sources"

die() { echo "ERROR: $*" >&2; exit 1; }

# ---- preflight --------------------------------------------------------------
echo "== preflight"
[ -n "$MANIFEST" ] || die "--bams <manifest.tsv> is required (template: paper/edger_bam_manifest.template.tsv)"
[ -f "$MANIFEST" ] || die "BAM manifest not found: $MANIFEST"
case "$LIBRARY" in forward|reverse) ;; *) die "--library must be forward or reverse (unstranded cannot separate RNA+/RNA-)";; esac
command -v Rscript >/dev/null || die "Rscript not on PATH"
command -v python3 >/dev/null || die "python3 not on PATH"
Rscript -e 'for (p in c("Rsubread","edgeR","GenomicRanges","rtracklayer")) if (!requireNamespace(p, quietly=TRUE)) stop("missing R package: ", p)' \
  || die "install the R packages above (Bioconductor)"
for f in "$REPO/rna_concordance_stranded.R" "$REPO/dev/extract_track_diff_bed.py" \
         "$REPO/dev/bearing_diff_concordance.py" "$REPO/dev/region_match_ebko.py" \
         "$REPO/dev/make_bearing_bins_bed.py" "$REPO/workflow/config/regions.tsv"; do
  [ -f "$f" ] || die "missing repository file: $f"
done
CATS="$RESULTS/DN_rep1_cats.json"
[ -f "$CATS" ] || die "missing $CATS (finished BEARING run expected in --results)"
HAS_3T3=0
if awk -F'\t' 'NR>1 && $2=="S3T3"{f=1} END{exit !f}' "$MANIFEST"; then HAS_3T3=1; fi
COND_REQ="DP ProB EbKO"; [ "$HAS_3T3" = 1 ] && COND_REQ="$COND_REQ S3T3"
for c in $COND_REQ; do
  [ -f "$RESULTS/pvalue/diff_DN_vs_${c}.stats.tsv" ] || die "missing $RESULTS/pvalue/diff_DN_vs_${c}.stats.tsv"
done

# BAM vectors per condition, from the manifest (header: sample condition bam)
bams_for() {
  awk -F'\t' -v c="$1" 'NR>1 && $2==c {printf "%s\"%s\"", (n++?", ":""), $3}' "$MANIFEST"
}
for c in DN $COND_REQ; do
  n=$(awk -F'\t' -v c="$c" 'NR>1 && $2==c' "$MANIFEST" | wc -l)
  [ "$n" -ge 2 ] || die "manifest has $n BAM(s) for $c; need at least 2 replicates"
  while IFS= read -r b; do [ -f "$b" ] || die "BAM not found: $b"; done \
    < <(awk -F'\t' -v c="$c" 'NR>1 && $2==c {print $3}' "$MANIFEST")
done
echo "   results : $RESULTS"
echo "   outdir  : $OUTDIR"
echo "   library : $LIBRARY  paired=$PAIRED  dups_marked=$DUPS  threads=$THREADS"
echo "   S3T3    : $([ "$HAS_3T3" = 1 ] && echo "yes (DN vs S3T3 arm runs)" || echo "no S3T3 rows; arm skipped")"
mkdir -p "$OUTDIR/configs" "$TS"

# ---- the BEARING bin grid (same chrom sizes and blacklist as the run) --------
BINS="$OUTDIR/bearing_bins_200bp.bed"
BL="$RESULTS/blacklist/merged_blacklist.bed"
[ -f "$BL" ] || BL="$REPO/workflow/resources/mm10-blacklist.v2.bed"
if [ ! -s "$BINS" ]; then
  echo "== bin grid -> $BINS (blacklist: $BL)"
  [ "$DRY" = 1 ] || python3 "$REPO/dev/make_bearing_bins_bed.py" \
      --chrom-sizes "$REPO/workflow/resources/mm10.chrom.sizes" \
      --blacklist "$BL" --out "$BINS"
fi

# ---- one rna_concordance_stranded.R run per config ---------------------------
run_one() {  # comp condA condB min_count region tag
  local comp="$1" A="$2" B="$3" mc="$4" region="$5" tag="$6"
  local cfg="$OUTDIR/configs/${tag}.R"
  cat > "$cfg" <<EOF
# generated by paper/run_edger_validation.sh $(date +%Y-%m-%d)
comp   <- "$comp"; condA <- "$A"; condB <- "$B"
bam_A  <- c($(bams_for "$A"))
bam_B  <- c($(bams_for "$B"))
paired_end <- $PAIRED; dups_marked <- $DUPS
library_strand <- "$LIBRARY"
bearing_dir <- "$REPO"
helper_dir  <- "$REPO/dev"
regions_tsv <- "$REPO/workflow/config/regions.tsv"
bins_bed    <- "$BINS"
bearing_diff_tsv <- "$RESULTS/pvalue/diff_${comp}.stats.tsv"
categories_json  <- "$CATS"
bearing_select <- "significant"
min_count <- $mc
concord_region <- "$region"
outdir <- "$OUTDIR"
prefix <- "rna_concordance_${comp}"
EOF
  echo "== $tag  ($comp, min.count=$mc${region:+, region=$region})"
  if [ "$DRY" = 1 ]; then echo "   [dry] Rscript rna_concordance_stranded.R $cfg"; return; fi
  ( cd "$OUTDIR" && PYTHON=python3 Rscript "$REPO/rna_concordance_stranded.R" "$cfg" ) \
      > "$OUTDIR/${tag}.log" 2>&1 || die "$tag failed; see $OUTDIR/${tag}.log"
  grep -E "^\[summary\]|recovery|direction" "$OUTDIR/${tag}.log" | sed 's/^/   /' | head -8 || true
}

run_one DN_vs_DP   DN DP   10 ""    S8_DN_vs_DP_mc10
run_one DN_vs_DP   DN DP   3  ""    S8_DN_vs_DP_mc3
run_one DN_vs_ProB DN ProB 10 "igh" igh_DN_vs_ProB_mc10
run_one DN_vs_ProB DN ProB 3  "igh" igh_DN_vs_ProB_mc3
run_one DN_vs_EbKO DN EbKO 10 ""    ebko_counts_mc10
if [ "$HAS_3T3" = 1 ]; then
  run_one DN_vs_S3T3 DN S3T3 10 ""  S3T3_DN_vs_S3T3_mc10
  run_one DN_vs_S3T3 DN S3T3 3  ""  S3T3_DN_vs_S3T3_mc3
fi

# ---- EbKO recombination center ----------------------------------------------
echo "== EbKO recombination center (region_match_ebko.py)"
if [ "$DRY" != 1 ]; then
  python3 "$REPO/dev/region_match_ebko.py" \
      --diff "$RESULTS/pvalue/diff_DN_vs_EbKO.stats.tsv" \
      --edger-pos "$OUTDIR/rna_concordance_DN_vs_EbKO_pos_edgeR_allbins.csv" \
      --edger-neg "$OUTDIR/rna_concordance_DN_vs_EbKO_neg_edgeR_allbins.csv" \
      --out "$OUTDIR/ebko_rc_characterisation.tsv" > "$OUTDIR/ebko_rc.log" 2>&1 \
    || die "region_match_ebko.py failed; see $OUTDIR/ebko_rc.log"
fi

# ---- hand the Table S8 sources to the workbook ------------------------------
if [ "$DRY" != 1 ]; then
  cp "$OUTDIR/rna_concordance_DN_vs_DP_summary.tsv"     "$TS/rna_concordance_DN_vs_DP_summary.tsv"
  cp "$OUTDIR/rna_concordance_DN_vs_DP_mc3_summary.tsv" "$TS/rna_concordance_DN_vs_DP_mc3_summary.tsv"
fi

cat <<EOF

Done. Outputs in $OUTDIR:
  rna_concordance_DN_vs_DP_summary.tsv, ..._mc3_summary.tsv   (Table S8; copied to table_sources/)
  rna_concordance_DN_vs_ProB_igh_summary.tsv, ..._mc3_igh_summary.tsv   (Igh sentence)
  ebko_rc_characterisation.tsv                                 (EbKO recombination-center paragraph)
  rna_concordance_DN_vs_S3T3_summary.tsv, ..._mc3_summary.tsv (only with S3T3 rows; DN vs 3T3)
  *_missed.tsv per strand (testable BEARING bins edgeR did not call, with logFC/FDR)
Rebuild the workbook:  cd workflow && snakemake ... paper_tables  (or rerun all)
Report recovery as a fraction of the TESTABLE BEARING bins, not the raw count.
EOF
