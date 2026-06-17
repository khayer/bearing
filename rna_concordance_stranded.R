#!/usr/bin/env Rscript
# rna_concordance_stranded.R
#
# Strand-matched RNA concordance between BEARING's per-track RNA differential and
# an independent edgeR test, for ONE comparison, BOTH strands, in one call.
#
# Why this exists: BEARING scores RNAseq+ and RNAseq- as SEPARATE strand-specific
# tracks. An unstranded featureCounts run (Strand="+", strandSpecific unset)
# counts reads from both strands, so it cannot validate a strand-specific track.
# This wrapper runs featureCounts with the correct strandSpecific setting and a
# per-strand feature Strand, so the plus-strand count validates RNAseq+ and the
# minus-strand count validates RNAseq-. It then hands each edgeR table to
# bearing_diff_concordance.py, which reports recovery on the TESTABLE subset
# (not the raw BEARING count), direction agreement, and rank enrichment with the
# saturated-universe handling -- i.e. the corrected statistics, not the inflated
# hypergeometric.
#
# Usage:
#   Rscript rna_concordance_stranded.R [config.R]
# where config.R (optional) overrides any variable in the CONFIG block below.
# Without a config file, edit the CONFIG block and run.
#
# Outputs (in outdir):
#   <prefix>_pos_edgeR_allbins.csv / <prefix>_neg_edgeR_allbins.csv
#   <prefix>_pos_bearing.bed       / <prefix>_neg_bearing.bed   (directional)
#   <prefix>_pos_concordance.tsv   / <prefix>_neg_concordance.tsv
#   <prefix>_summary.tsv           (both strands, one row each)
#
# ASCII only.

suppressPackageStartupMessages({
  library(Rsubread)
  library(edgeR)
  library(GenomicRanges)
  library(rtracklayer)
})

## ===========================================================================
## CONFIG  (override any of these in an optional config.R passed as arg 1)
## ===========================================================================

# Comparison: condA is the POSITIVE-logFC / A-enriched direction (matches
# BEARING's A-B sign for diff_<A>_vs_<B>). For diff_DN_vs_DP -> condA="DN".
comp   <- "DN_vs_DP"
condA  <- "DN"
condB  <- "DP"

# BAMs per condition (replicates).
bam_A <- c("../bam/DN_rep1.RNA.bam", "../bam/DN_rep2.RNA.bam")
bam_B <- c("../bam/DP_rep1.RNA.bam", "../bam/DP_rep2.RNA.bam")

paired_end  <- TRUE
dups_marked <- TRUE

# Library strandedness. CRITICAL -- getting this wrong swaps + and -.
#   "reverse"    : dUTP / most modern RNA-seq (read1 is antisense)   -> strandSpecific 2
#   "forward"    : read1 is sense                                    -> strandSpecific 1
#   "unstranded" : cannot separate strands (script will abort)       -> strandSpecific 0
library_strand <- "reverse"

bearing_dir <- "/mnt/isilon/bassing_lab/integration_paper/bearing"
bins_bed    <- file.path(bearing_dir, "bearing_bins_200bp.bed")

# BEARING differential to validate against (the diff stats for THIS comparison)
# and the cats JSON, used to extract the directional per-strand BEDs.
bearing_diff_tsv <- file.path(bearing_dir,
  sprintf("workflow/results/pvalue/diff_%s.stats.tsv", comp))
categories_json  <- file.path(bearing_dir, "DN_rep1_cats.json")

# How to pick BEARING's significant per-strand bins for the BED:
#   bearing_select="significant" -> bins flagged FDR-significant by the run
#   bearing_select="fdr"         -> bins with pval_adj_bh < bearing_fdr
#   bearing_select="top"         -> top bearing_top_percent of |kl|, |kl|>=bearing_min_abs
bearing_select     <- "significant"
bearing_fdr        <- 0.05
bearing_top_percent <- 0.01
bearing_min_abs    <- 0.0

fdr_cut  <- 0.05     # edgeR FDR for "independently significant"
rank_pct <- 10.0     # top-% of edgeR universe for rank enrichment
outdir   <- "."
prefix   <- sprintf("rna_concordance_%s", comp)

## ---- apply optional config file -------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1 && nzchar(.args[1])) {
  if (!file.exists(.args[1])) stop("config file not found: ", .args[1])
  message("[config] sourcing ", .args[1])
  source(.args[1], local = FALSE)
}

## ===========================================================================
## Setup + validation
## ===========================================================================

strand_specific <- switch(library_strand,
  "forward"    = 1L,
  "reverse"    = 2L,
  "unstranded" = 0L,
  stop("library_strand must be 'forward', 'reverse', or 'unstranded'"))
if (strand_specific == 0L) {
  stop("library_strand='unstranded' cannot separate RNAseq+/-; ",
       "use a stranded library or test only total RNA.")
}

bam_files <- c(bam_A, bam_B)
condition <- factor(c(rep(condA, length(bam_A)), rep(condB, length(bam_B))),
                    levels = c(condB, condA))   # coef 2 = condA-enriched positive

for (f in c(bam_files, bins_bed, bearing_diff_tsv, categories_json)) {
  if (!file.exists(f)) stop("missing input: ", f)
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
extract_py <- file.path(bearing_dir, "extract_track_diff_bed.py")
concord_py <- file.path(bearing_dir, "bearing_diff_concordance.py")
for (f in c(extract_py, concord_py)) if (!file.exists(f)) stop("missing helper: ", f)
python_bin <- Sys.getenv("PYTHON", unset = "python")

message(sprintf("[setup] %s : condA=%s (positive) vs condB=%s | library=%s -> strandSpecific=%d",
                comp, condA, condB, library_strand, strand_specific))
message("[setup] coef 2 sign: positive logFC = enriched in ", condA,
        " (matches BEARING A-B). Verify on a known stranded gene if unsure.")

bins <- import(bins_bed)
saf_base <- data.frame(
  GeneID = paste0("bin_", seq_along(bins)),
  Chr    = as.character(seqnames(bins)),
  Start  = start(bins),
  End    = end(bins),
  stringsAsFactors = FALSE)

## ===========================================================================
## Per-strand worker
## ===========================================================================

run_strand <- function(strand_label, strand_sym, bearing_track) {
  message(sprintf("\n==== %s strand (%s) : BEARING track '%s' ====",
                  strand_label, strand_sym, bearing_track))
  saf <- saf_base
  saf$Strand <- strand_sym   # count this transcription strand

  fc <- featureCounts(
    files                  = bam_files,
    annot.ext              = saf,
    isPairedEnd            = paired_end,
    strandSpecific         = strand_specific,
    countMultiMappingReads = FALSE,
    ignoreDup              = dups_marked,
    minMQS                 = 20,
    nthreads               = 4)
  counts <- fc$counts
  colnames(counts) <- sub("\\.bam$", "", basename(bam_files))

  y <- DGEList(counts = counts, group = condition)
  keep <- filterByExpr(y, group = condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~condition)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf, sort.by = "none")$table
  kept_idx  <- as.integer(sub("bin_", "", rownames(res)))
  res$chr   <- as.character(seqnames(bins))[kept_idx]
  res$start <- start(bins)[kept_idx]
  res$end   <- end(bins)[kept_idx]

  edger_csv <- file.path(outdir, sprintf("%s_%s_edgeR_allbins.csv", prefix, strand_label))
  write.csv(res[, c("logFC", "logCPM", "F", "PValue", "FDR", "chr", "start", "end")],
            edger_csv, row.names = FALSE)
  message("  wrote ", edger_csv, " (", nrow(res), " tested bins)")

  # Directional BEARING BED for this strand's track.
  bearing_bed <- file.path(outdir, sprintf("%s_%s_bearing.bed", prefix, strand_label))
  ex_args <- c(extract_py, "--diff-tsv", bearing_diff_tsv,
               "--categories", categories_json, "--track", bearing_track,
               "--out", bearing_bed)
  if (bearing_select == "significant") {
    ex_args <- c(ex_args, "--significant")          # use the run's FDR-sig flag
  } else if (bearing_select == "fdr") {
    ex_args <- c(ex_args, "--fdr", format(bearing_fdr))
  } else {
    ex_args <- c(ex_args, "--top-percent", format(bearing_top_percent),
                 "--min-abs", format(bearing_min_abs))
  }
  st <- system2(python_bin, ex_args, stdout = TRUE, stderr = TRUE)
  cat(paste(st, collapse = "\n"), "\n")

  # Concordance (recovery-on-testable + direction + rank enrichment).
  concord_tsv <- file.path(outdir, sprintf("%s_%s_concordance.tsv", prefix, strand_label))
  co_args <- c(concord_py, "--bearing-bed", bearing_bed, "--edger-csv", edger_csv,
               "--fdr", format(fdr_cut), "--rank-pct", format(rank_pct),
               "--out", concord_tsv)
  co <- system2(python_bin, co_args, stdout = TRUE, stderr = TRUE)
  cat(paste(co, collapse = "\n"), "\n")

  list(strand = strand_label, track = bearing_track, edger_csv = edger_csv,
       bearing_bed = bearing_bed, concordance_tsv = concord_tsv,
       n_tested = nrow(res))
}

## ===========================================================================
## Run both strands + combine
## ===========================================================================

res_pos <- run_strand("pos", "+", "RNAseq +")
res_neg <- run_strand("neg", "-", "RNAseq -")

read_concord <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.delim(path, stringsAsFactors = FALSE), error = function(e) NULL)
}
sum_pos <- read_concord(res_pos$concordance_tsv)
sum_neg <- read_concord(res_neg$concordance_tsv)
combined <- do.call(rbind, Filter(Negate(is.null), list(
  if (!is.null(sum_pos)) cbind(strand = "pos", track = "RNAseq+", sum_pos),
  if (!is.null(sum_neg)) cbind(strand = "neg", track = "RNAseq-", sum_neg))))

summary_tsv <- file.path(outdir, sprintf("%s_summary.tsv", prefix))
if (!is.null(combined)) {
  write.table(combined, summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  message("\n[summary] wrote ", summary_tsv)
  print(combined)
} else {
  message("\n[summary] concordance tools produced no parseable summary; ",
          "see the per-strand *_concordance.tsv and stdout above.")
}

message("\nDone. Report recovery as a fraction of the TESTABLE BEARING bins ",
        "(not the raw BEARING count), and quote the rank enrichment rather than ",
        "the hypergeometric p when the edgeR universe is saturated.")
