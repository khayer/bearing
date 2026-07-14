#!/usr/bin/env python3
"""
run_bes_hic_matrix.py -- run bes_hic_correlation.py across a matrix of
(locus x condition-pair) comparisons and tabulate the headline result into one
manuscript-ready table.

The headline is the shift-null CONTACT co-localization: BES (p95) vs
delta_contact, reporting spearman rho and the circular-shift empirical p
(p_emp). The parametric Spearman p and the insulation columns are written to a
separate secondary file, since they are inflated by spatial autocorrelation and
do not belong in the headline.

INPUT
  --comparisons TSV with a header and columns: locus  region  condA  condB
      e.g.
        locus  region                      condA  condB
        Tcrb   chr6:40400000-42400000      DN     DP
        Igh    chr12:112934363-116108354   DN     ProB
        EBKO   chr6:40400000-42400000      DN     EBKO

  Path templates (use {cond} for single-condition files, {a}/{b} for the pair):
  --cool-template   e.g. ../hic_files/merged_corrected_KR_{cond}_bs_10000.cool
  --insul-template  e.g. .../07tad/merged_corrected_KR_{cond}_bs_10000_tad_score.bm
  --bes-template    e.g. workflow/results_adaptive/pvalue/diff_{a}_vs_{b}.stats.tsv

OUTPUT (into --out-dir)
  bes_hic_main_summary.tsv     one row per (comparison, locus): contact headline
  bes_hic_pertrack_summary.tsv one row per (comparison, locus, track): contact strata
  bes_hic_insulation.tsv       secondary: insulation rho/p_emp (not headline)
  <prefix>.stats.json/.tsv/.pdf per comparison (from bes_hic_correlation.py)

ASCII only.
"""
import argparse
import json
import os
import subprocess
import sys


def get(d, *keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d and d[k] is not None:
            return d[k]
    return default


def fmt(x, nd=4):
    if x is None:
        return "NA"
    try:
        return ("%." + str(nd) + "g") % float(x)
    except (TypeError, ValueError):
        return str(x)


def read_comparisons(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            if header is None:
                header = [c.lower() for c in f]
                continue
            rec = dict(zip(header, f))
            for need in ("locus", "region", "conda", "condb"):
                if need not in rec:
                    sys.exit("[ERROR] comparisons TSV needs columns "
                             "locus/region/condA/condB; missing %s" % need)
            rows.append(rec)
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--comparisons", required=True)
    ap.add_argument("--script", default="bes_hic_correlation.py")
    ap.add_argument("--cool-template", required=True)
    ap.add_argument("--insul-template", required=True)
    ap.add_argument("--bes-template", required=True)
    ap.add_argument("--bes-score-col", default="bearing_score")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--hic-bin", type=int, default=10000)
    ap.add_argument("--min-distance", type=int, default=50000)
    ap.add_argument("--max-distance", type=int, default=500000)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--reuse", action="store_true",
                    help="skip the run if <prefix>.stats.json already exists")
    ap.add_argument("--dry-run", action="store_true",
                    help="print commands, do not run")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    comps = read_comparisons(args.comparisons)

    main_rows, track_rows, insul_rows = [], [], []
    for rec in comps:
        locus, region = rec["locus"], rec["region"]
        a, b = rec["conda"], rec["condb"]
        prefix = os.path.join(args.out_dir, "%s_vs_%s_%s" % (a, b, locus))
        statp = prefix + ".stats.json"

        cool_a = args.cool_template.format(cond=a)
        cool_b = args.cool_template.format(cond=b)
        insul_a = args.insul_template.format(cond=a)
        insul_b = args.insul_template.format(cond=b)
        bes = args.bes_template.format(a=a, b=b)

        cmd = [sys.executable, args.script,
               "--bes", bes, "--bes-score-col", args.bes_score_col,
               "--insul-A", insul_a, "--insul-B", insul_b,
               "--cool-A", cool_a, "--cool-B", cool_b,
               "--region", region, "--hic-bin", str(args.hic_bin),
               "--min-distance", str(args.min_distance),
               "--max-distance", str(args.max_distance),
               "--n-perm", str(args.n_perm), "--prefix", prefix]

        if args.dry_run:
            print(" ".join(cmd))
            continue
        if not (args.reuse and os.path.exists(statp)):
            sys.stderr.write("[run] %s vs %s @ %s\n" % (a, b, locus))
            r = subprocess.run(cmd)
            if r.returncode != 0 or not os.path.exists(statp):
                sys.stderr.write("[WARN] failed: %s vs %s @ %s\n" % (a, b, locus))
                continue
        else:
            sys.stderr.write("[reuse] %s\n" % statp)

        with open(statp) as fh:
            S = json.load(fh)
        corr = S.get("correlations", {})
        nbin = S.get("n_hic_bins")

        c = corr.get("p95_bes__delta_contact", {})
        enr = c.get("top_decile_enrichment", {}) or {}
        main_rows.append([
            "%s_vs_%s" % (a, b), locus, region, nbin,
            fmt(c.get("spearman_rho")), fmt(c.get("empirical_p")),
            fmt(c.get("spearman_p")),
            fmt(get(enr, "odds_ratio", "OR", "oddsratio")),
            fmt(get(enr, "fisher_p", "p")),
        ])

        ins = corr.get("p95_bes__delta_insulation", {})
        einr = ins.get("top_decile_enrichment", {}) or {}
        insul_rows.append([
            "%s_vs_%s" % (a, b), locus,
            fmt(ins.get("spearman_rho")), fmt(ins.get("empirical_p")),
            fmt(ins.get("spearman_p")),
            fmt(get(einr, "odds_ratio", "OR", "oddsratio")),
            fmt(get(einr, "fisher_p", "p")),
        ])

        for tc, row in (S.get("track_stratified", {}) or {}).items():
            dc = row.get("delta_contact", {}) or {}
            track_rows.append([
                "%s_vs_%s" % (a, b), locus, tc,
                row.get("n_bins_dominated"),
                fmt(dc.get("spearman_rho")), fmt(dc.get("spearman_p")),
                fmt(dc.get("empirical_p")),
            ])

    def write(path, header, rows):
        with open(path, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join("NA" if x is None else str(x) for x in r) + "\n")
        sys.stderr.write("[table] %s (%d rows)\n" % (path, len(rows)))

    write(os.path.join(args.out_dir, "bes_hic_main_summary.tsv"),
          ["comparison", "locus", "region", "n_hic_bins",
           "rho_contact", "p_emp_contact", "p_param_contact",
           "topdecile_OR", "fisher_p"], main_rows)
    write(os.path.join(args.out_dir, "bes_hic_pertrack_summary.tsv"),
          ["comparison", "locus", "track", "n_bins_dominated",
           "rho_contact", "p_param_contact", "p_emp_contact"], track_rows)
    write(os.path.join(args.out_dir, "bes_hic_insulation.tsv"),
          ["comparison", "locus", "rho_insul", "p_emp_insul",
           "p_param_insul", "topdecile_OR", "fisher_p"], insul_rows)


if __name__ == "__main__":
    main()
