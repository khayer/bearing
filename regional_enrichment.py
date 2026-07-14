#!/usr/bin/env python3
"""
regional_enrichment.py -- Regional enrichment testing for BEARING differential analyses.

Given a set of pre-specified genomic regions and a pairwise differential analysis,
this script reports whether each region is significantly enriched for low-p bins
(spatial concentration) and whether the differential signs in that region are
concordant (directional consistency). This complements the per-bin BH-corrected
q-values from bearing_pvalue.py, which are conservative for sub-locus claims.

The test combines two component statistics via Fisher's method:
- Spatial concentration: binomial test for enrichment of low-p bins within the region
  relative to the locus-wide proportion (conditional on total low-p bins in locus).
- Directional concordance: binomial test for deviation from 50/50 A>B vs B>A signs
  among low-p bins in the region.
- Combined: Fisher's chi-squared combination, then BH correction across all tests.
"""

import argparse
import csv
import json
import logging
import math
import sys
from pathlib import Path

import numpy as np
import scipy.stats

logger = logging.getLogger(__name__)


def bh_fdr(pvals, alpha=0.05):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    pvals_sorted = pvals[order]
    ranks = np.arange(1, n + 1)
    pvals_adj = np.minimum(1.0, pvals_sorted * n / ranks)
    pvals_adj_sorted = pvals_adj.copy()
    for i in range(n - 2, -1, -1):
        pvals_adj_sorted[i] = min(pvals_adj_sorted[i], pvals_adj_sorted[i + 1])
    pvals_adj[order] = pvals_adj_sorted
    rejected = pvals_adj <= alpha
    return rejected, pvals_adj


def parse_ucsc_region(region_str):
    """Parse UCSC region string (chr:start-end) to tuple."""
    parts = region_str.split(':')
    if len(parts) != 2:
        raise ValueError("Region must be in format chrom:start-end")
    chrom = parts[0]
    coords = parts[1].split('-')
    if len(coords) != 2:
        raise ValueError("Region must be in format chrom:start-end")
    start, end = int(coords[0]), int(coords[1])
    return chrom, start, end


def parse_diff_qcat(path):
    """
    Parse differential qcat file into dict keyed by (chrom, start, end).
    Supports both JSON format (qcat files) and simple bearing_score column (stats files).
    Handles files with or without headers.
    """
    bins = {}
    
    with open(path, 'r') as f:
        first_line = f.readline().strip()
        f.seek(0)  # Reset to start
        
        # Check if first line is a header
        has_header = first_line.startswith('chrom') or '\t' in first_line and first_line.split('\t')[0] == 'chrom'
        
        if has_header:
            # Parse as DictReader (header-based)
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                chrom = row['chrom']
                start = int(row['start'])
                end = int(row['end'])
                
                # Try JSON format first (qcat files)
                if 'qcat' in row:
                    try:
                        qcat_data = json.loads(row['qcat'])
                        score = qcat_data.get('qcat', [[0, 0]])[0][0]
                    except (json.JSONDecodeError, IndexError, KeyError, TypeError):
                        score = 0.0
                # Fall back to bearing_score column (stats files)
                elif 'bearing_score' in row:
                    try:
                        score = float(row['bearing_score'])
                    except (ValueError, KeyError):
                        score = 0.0
                else:
                    score = 0.0
                
                bins[(chrom, start, end)] = score
        else:
            # Parse as raw lines (no header)
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4:
                    continue
                try:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    qcat_json = parts[3]
                    qcat_data = json.loads(qcat_json)
                    score = qcat_data.get('qcat', [[0, 0]])[0][0]
                except (json.JSONDecodeError, IndexError, KeyError, ValueError):
                    score = 0.0
                
                bins[(chrom, start, end)] = score
    
    return bins


def parse_diff_pvals(path):
    """
    Parse differential p-values file into dict keyed by (chrom, start, end).
    Returns dict of bin -> (p_value, signed_score).
    """
    bins = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = row['chrom']
            start = int(row['start'])
            end = int(row['end'])
            pval = float(row['pval'])
            score = float(row['bearing_score'])
            bins[(chrom, start, end)] = (pval, score)
    return bins


def parse_regions_bed(path):
    """Parse BED file of regions."""
    regions = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            regions.append((chrom, start, end, name))
    return regions


def compute_regional_enrichment(diff_qcat_bins, diff_pvals_bins, regions,
                                locus_chrom, locus_start, locus_end, p_thresh,
                                verbose=False, region_assign='start'):
    """
    Compute regional enrichment for each region.
    Returns list of result dicts, one per region.
    """
    # Filter bins to locus
    locus_bins = [
        coord for coord in diff_qcat_bins.keys()
        if coord[0] == locus_chrom and locus_start <= coord[1] < locus_end
    ]
    locus_pvals = {
        coord: data for coord, data in diff_pvals_bins.items()
        if coord[0] == locus_chrom and locus_start <= coord[1] < locus_end
    }

    if not locus_bins:
        raise ValueError("No bins found in locus {}:{}->{}".format(
            locus_chrom, locus_start, locus_end))

    # Count total locus bins and low-p bins
    n_locus = sum(1 for p_val, _ in locus_pvals.values() if p_val < p_thresh)
    L_locus_bins = len(locus_bins)

    # Genome-wide directional baseline g for the directional-concordance null.
    # Computed over ALL significant bins in the diff table (genome-wide when the
    # full diff stats are supplied), NOT within the locus, which would be
    # circular. g = P(signed_score > 0 | bin significant). The directional test
    # below uses g rather than 0.5 so that comparisons with a global sign
    # imbalance (e.g. a condition with systematically higher architectural
    # signal genome-wide) are not flagged as regionally concordant merely for
    # reproducing that genome-wide skew.
    gw_signs = [s for (p_val, s) in diff_pvals_bins.values() if p_val < p_thresh]
    n_gw_sig = len(gw_signs)
    if n_gw_sig > 0:
        g_dir = sum(1 for s in gw_signs if s > 0) / n_gw_sig
    else:
        g_dir = 0.5
    # Guard against degenerate 0/1 baselines that would make the binomial null
    # ill-defined.
    g_dir = min(max(g_dir, 1e-6), 1.0 - 1e-6)

    results = []
    for chrom, start, end, name in regions:
        # Validate region is within locus
        if chrom != locus_chrom or start < locus_start or end > locus_end:
            raise ValueError(
                "Region {} ({}:{}-{}) extends beyond locus {}:{}-{}".format(
                    name, chrom, start, end, locus_chrom, locus_start, locus_end))

        # Count region bins. Two assignment modes:
        #   'start'   (default): a bin belongs to the region if its start falls
        #             inside [start, end). Correct and fast for regions larger
        #             than the bin grid (the standard locus sub-regions); this
        #             preserves the original, validated behaviour.
        #   'overlap': a bin belongs if it overlaps the region by >=1 bp. Needed
        #             for features SMALLER than the bin size (e.g. 18 bp CTCF-
        #             binding elements), which under 'start' capture zero bins
        #             because no 200 bp grid start lands inside them.
        if region_assign == 'overlap':
            region_bins = [
                coord for coord in locus_bins
                if coord[0] == chrom and coord[1] < end and coord[2] > start
            ]
        else:
            region_bins = [
                coord for coord in locus_bins
                if coord[0] == chrom and start <= coord[1] < end
            ]
        L_region_bins = len(region_bins)
        if L_region_bins == 0:
            logger.warning(
                "region %s (%s:%d-%d) contains zero bins under assign='%s'; "
                "its q_combined will be a trivial 1.0. If this feature is "
                "smaller than the bin size, re-run with --region-assign overlap.",
                name, chrom, start, end, region_assign)

        # Count low-p bins in region
        k = 0
        k_pos = 0
        k_neg = 0
        for coord in region_bins:
            if coord in locus_pvals:
                p_val, signed_score = locus_pvals[coord]
                if p_val < p_thresh:
                    k += 1
                    if signed_score > 0:
                        k_pos += 1
                    else:
                        k_neg += 1

        # Compute proportions
        pi = L_region_bins / L_locus_bins if L_locus_bins > 0 else 0.0

        # Spatial concentration test
        if n_locus == 0:
            p_spatial = 1.0
        else:
            p_spatial = scipy.stats.binom.sf(k - 1, n=n_locus, p=pi)
            p_spatial = min(p_spatial, 1.0)

        # Directional concordance test. One-sided exact binomial test for EXCESS
        # concordance in the region's own majority direction, relative to the
        # genome-wide directional baseline g_dir (not 0.5). A region is called
        # concordant only if its dominant-direction bin count exceeds what the
        # genome-wide sign balance predicts, in whichever direction the region
        # leans (so genuine against-the-skew regional effects are still
        # detected); a region that is merely more balanced than a skewed
        # background is not flagged.
        if k == 0:
            p_directional = 1.0
        else:
            if k_pos >= k_neg:
                p_directional = float(scipy.stats.binom.sf(k_pos - 1, n=k, p=g_dir))
            else:
                p_directional = float(scipy.stats.binom.sf(k_neg - 1, n=k, p=1.0 - g_dir))
            p_directional = min(p_directional, 1.0)

        # Fisher's combined
        if p_spatial == 0 or p_directional == 0:
            p_combined = 0.0
        else:
            T = -2 * (math.log(p_spatial) + math.log(p_directional))
            p_combined = scipy.stats.chi2.sf(T, df=4)

        result = {
            'region_name': name,
            'chrom': chrom,
            'start': start,
            'end': end,
            'L_region_bins': L_region_bins,
            'L_locus_bins': L_locus_bins,
            'pi': pi,
            'k': k,
            'n_locus': n_locus,
            'k_pos': k_pos,
            'k_neg': k_neg,
            'g_dir': g_dir,
            'p_spatial': p_spatial,
            'p_directional': p_directional,
            'p_combined': p_combined,
        }

        if verbose:
            logger.info("[regional_enrichment] {}".format(name))
            logger.info("  region: {}:{}-{}  ({} valid bins, {:.4f} of locus)".format(
                chrom, start, end, L_region_bins, pi))
            logger.info("  k = {} bins below {} in region".format(k, p_thresh))
            logger.info("  n_locus = {} bins below {} in locus".format(n_locus, p_thresh))
            logger.info("  pi = {:.4f}".format(pi))
            logger.info("  directional: {} A>B, {} B>A".format(k_pos, k_neg))
            logger.info("  p_spatial = {:.3e}, p_directional = {:.3e}, p_combined = {:.3e}".format(
                p_spatial, p_directional, p_combined))

        results.append(result)

    return results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest='mode', help='Mode of operation')

    # Single comparison mode
    single_parser = subparsers.add_parser('single', help='Single comparison analysis')
    single_parser.add_argument('--diff-qcat', required=True,
                               help='Differential qcat file from compare_qcat.py')
    single_parser.add_argument('--diff-pvals', required=True,
                               help='Per-bin p-values TSV from bearing_pvalue.py --diff')
    single_parser.add_argument('--regions', required=True,
                               help='BED file of regions to test')
    single_parser.add_argument('--locus', required=True,
                               help='Analysis locus (chrom:start-end)')
    single_parser.add_argument('--p-thresh', type=float, default=0.05,
                               help='P-value threshold for "significant" bins (default: 0.05)')
    single_parser.add_argument('--region-assign', choices=['start', 'overlap'], default='start',
                               help="Bin-to-region assignment. 'start' (default): a bin belongs "
                                    "if its start falls inside the region. 'overlap': a bin "
                                    "belongs if it overlaps the region by >=1 bp; use for features "
                                    "smaller than the bin size, e.g. CTCF-binding elements.")
    single_parser.add_argument('--out', required=True,
                               help='Output TSV path')
    single_parser.add_argument('--bh-by', choices=['none', 'comparison', 'all'], default='all',
                               help='BH correction scope (default: all)')
    single_parser.add_argument('--verbose', '-v', action='store_true',
                               help='Print diagnostic info per region')

    # Batch mode
    batch_parser = subparsers.add_parser('batch', help='Batch analysis of multiple comparisons')
    batch_parser.add_argument('--diff-table', required=True,
                              help='TSV with columns: comparison_label, diff_qcat_path, diff_pvals_path')
    batch_parser.add_argument('--regions', required=True,
                              help='BED file of regions to test')
    batch_parser.add_argument('--locus', required=True,
                              help='Analysis locus (chrom:start-end)')
    batch_parser.add_argument('--p-thresh', type=float, default=0.05,
                              help='P-value threshold for "significant" bins (default: 0.05)')
    batch_parser.add_argument('--region-assign', choices=['start', 'overlap'], default='start',
                              help="Bin-to-region assignment. 'start' (default): a bin belongs "
                                   "if its start falls inside the region. 'overlap': a bin "
                                   "belongs if it overlaps the region by >=1 bp; use for features "
                                   "smaller than the bin size, e.g. CTCF-binding elements.")
    batch_parser.add_argument('--out', required=True,
                              help='Output TSV path')
    batch_parser.add_argument('--bh-by', choices=['none', 'comparison', 'all'], default='all',
                              help='BH correction scope (default: all)')
    batch_parser.add_argument('--verbose', '-v', action='store_true',
                              help='Print diagnostic info per region')

    args = parser.parse_args()

    if not args.mode:
        parser.error("Must specify mode: single or batch")

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format='[regional_enrichment] %(message)s'
    )

    # Parse locus
    locus_chrom, locus_start, locus_end = parse_ucsc_region(args.locus)

    # Parse regions
    regions = parse_regions_bed(args.regions)

    # Collect all results for BH correction
    all_results = []

    if args.mode == 'single':
        comparison_label = Path(args.diff_qcat).stem.replace('.qcat', '').replace('compare_', '')

        # Load data
        diff_qcat_bins = parse_diff_qcat(args.diff_qcat)
        diff_pvals_bins = parse_diff_pvals(args.diff_pvals)

        # Compute
        results = compute_regional_enrichment(
            diff_qcat_bins, diff_pvals_bins, regions,
            locus_chrom, locus_start, locus_end, args.p_thresh, args.verbose, region_assign=args.region_assign
        )

        # Add comparison label
        for r in results:
            r['comparison'] = comparison_label

        all_results.extend(results)

    elif args.mode == 'batch':
        # Check if diff_table is a metadata file (with comparison_label column) or a stats file
        with open(args.diff_table, 'r') as f:
            first_line = f.readline().strip()
            headers = first_line.split('\t')
            has_comparison_label = 'comparison_label' in headers
        
        if has_comparison_label:
            # Old format: metadata file with paths to multiple comparisons
            comparisons = []
            with open(args.diff_table, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    comparisons.append(row)

            for row in comparisons:
                comparison_label = row['comparison_label']
                diff_qcat_path = row['diff_qcat_path']
                diff_pvals_path = row['diff_pvals_path']

                logger.info("Processing comparison: {}".format(comparison_label))

                # Load data
                diff_qcat_bins = parse_diff_qcat(diff_qcat_path)
                diff_pvals_bins = parse_diff_pvals(diff_pvals_path)

                # Compute
                results = compute_regional_enrichment(
                    diff_qcat_bins, diff_pvals_bins, regions,
                    locus_chrom, locus_start, locus_end, args.p_thresh, args.verbose, region_assign=args.region_assign
                )

                # Add comparison label
                for r in results:
                    r['comparison'] = comparison_label

                all_results.extend(results)
        else:
            # New format: direct stats file with pval and bearing_score columns
            comparison_label = Path(args.diff_table).stem
            logger.info("Processing stats file: {}".format(comparison_label))
            
            # Use same file for both qcat and pvals (it has both columns)
            diff_qcat_bins = parse_diff_qcat(args.diff_table)
            diff_pvals_bins = parse_diff_pvals(args.diff_table)
            
            # Compute
            results = compute_regional_enrichment(
                diff_qcat_bins, diff_pvals_bins, regions,
                locus_chrom, locus_start, locus_end, args.p_thresh, args.verbose, region_assign=args.region_assign
            )
            
            # Add comparison label
            for r in results:
                r['comparison'] = comparison_label
            
            all_results.extend(results)

    # BH correction
    if all_results:
        if args.bh_by == 'all':
            p_combined_all = [r['p_combined'] for r in all_results]
            _, q_combined_all = bh_fdr(p_combined_all, alpha=0.05)
            for r, q in zip(all_results, q_combined_all):
                r['q_combined'] = q
        else:
            for r in all_results:
                r['q_combined'] = r['p_combined']

        n_sig = sum(1 for r in all_results if r['q_combined'] < 0.05)
        logger.info("BH correction across {} tests (scope: {})".format(
            len(all_results), args.bh_by))
        logger.info("Significant at q < 0.05: {}/{}".format(n_sig, len(all_results)))

    # Write output
    with open(args.out, 'w', newline='') as f:
        fieldnames = ['comparison', 'region_name', 'chrom', 'start', 'end',
                      'L_region_bins', 'L_locus_bins', 'pi', 'k', 'n_locus',
                      'k_pos', 'k_neg', 'p_spatial', 'p_directional',
                      'p_combined', 'q_combined', 'g_dir']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for r in all_results:
            formatted = {}
            for k, v in r.items():
                if isinstance(v, float):
                    formatted[k] = "{:.6g}".format(v)
                else:
                    formatted[k] = str(v)
            writer.writerow(formatted)

    logger.info("Output written to {}".format(args.out))


if __name__ == '__main__':
    main()
