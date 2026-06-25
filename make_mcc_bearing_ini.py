#!/usr/bin/env python3
"""
make_mcc_bearing_ini.py -- write a pyGenomeTracks .ini that stacks, at one locus:

  [MCC matrices]  one capture cool per condition (triangle)
  [BES]           the BEARING differential track (diff_*.neglog10p.bw)
  [panel]         the five V1 assay tracks (RNAseq +/-, CTCF, Cohesin, NIPBL),
                  colored from cats_v1.json, optionally overlaying DN (faint)
                  under the mutant (solid) so the assay-level change is visible
  [genes/features] Tcrb annotation + the CBE anchors

The point of the figure: BEARING on a DIFFERENT sample set (Allyn 2025 V1 mutants)
and a DIFFERENT assay mix (NIPBL in, ATAC/H3K27ac out, capture Hi-C as the 3D
side) -- a portability supplement.

  python make_mcc_bearing_ini.py --sheet workflow/config/samples_v1.tsv \
    --cats categories/v1_5track_panel.yaml \
    --bes workflow/results_v1/pvalue/diff_DN_vs_RCTKO.neglog10p.bw \
    --cools DN:capture_cools/cool/merged_capture_DN_bs_2000.cool,\
RCTKO:capture_cools/cool/merged_capture_RCTKO_bs_2000.cool \
    --panel-condition RCTKO --overlay-condition DN \
    --genes annotations/tcrb_genes.bed --features annotations/cbe_mm10.bed \
    --region chr6:40800000-41650000 --out-ini tracks_mcc_rctko.ini
  pyGenomeTracks --tracks tracks_mcc_rctko.ini --region chr6:40800000-41650000 \
    -o mcc_bearing_RCTKO.pdf

cats accepts the repo YAML panel OR a cats JSON. ASCII only.
"""
import argparse
import json
import os
import sys


def read_sheet_conditions(path):
    """Return {condition: [bw paths]} using the first replicate seen per condition.
    Skips comment/blank lines before and within the sheet (matches the BEARING
    convention)."""
    cond = {}
    with open(path) as fh:
        header = None
        for line in fh:
            s = line.rstrip("\n")
            if not s.strip() or s.startswith("#"):
                continue
            header = s.split("\t")
            break
        idx = {h: i for i, h in enumerate(header)}
        ci, bi = idx.get("condition", 1), idx.get("bw", 4)
        for line in fh:
            s = line.rstrip("\n")
            if not s.strip() or s.startswith("#"):
                continue
            f = s.split("\t")
            if len(f) <= bi:
                continue
            c = f[ci]
            if c not in cond:
                cond[c] = [p.strip() for p in f[bi].split(",") if p.strip()]
    return cond


def read_cats(path):
    """Return ordered [(label, color, negative_strand)] from repo YAML or cats JSON."""
    if path.endswith(".json"):
        d = json.load(open(path))["categories"]
        out = []
        for k in sorted(d, key=lambda x: int(x)):
            name, color = d[k][0], d[k][1]
            out.append((name, color, name.strip().endswith("-")))
        return out
    # minimal YAML reader for the repo panel format (no pyyaml dependency)
    out, name, color, neg = [], None, None, False
    for line in open(path):
        t = line.strip()
        if t.startswith("- name:"):
            if name is not None:
                out.append((name, color or "#000000", neg))
            name = t.split(":", 1)[1].strip().strip('"')
            color, neg = None, False
        elif t.startswith("color:"):
            color = t.split(":", 1)[1].strip().strip('"')
        elif t.startswith("negative_strand:"):
            neg = t.split(":", 1)[1].strip().lower() == "true"
    if name is not None:
        out.append((name, color or "#000000", neg))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--cats", required=True)
    ap.add_argument("--bes", required=True)
    ap.add_argument("--cools", required=True,
                    help="comma list of LABEL:path.cool (matrices, top to bottom)")
    ap.add_argument("--panel-condition", required=True,
                    help="which condition's 5 assay tracks to draw")
    ap.add_argument("--overlay-condition", default=None,
                    help="optional: overlay this condition faint under each track")
    ap.add_argument("--genes", default=None)
    ap.add_argument("--features", default=None)
    ap.add_argument("--region", required=True, help="chrom:start-end (for the title)")
    ap.add_argument("--depth", type=int, default=400000, help="matrix triangle depth bp")
    ap.add_argument("--bes-label", default=None)
    ap.add_argument("--out-ini", required=True)
    args = ap.parse_args()

    cond_bw = read_sheet_conditions(args.sheet)
    cats = read_cats(args.cats)
    if args.panel_condition not in cond_bw:
        sys.exit("[ERROR] condition %s not in sheet (have: %s)"
                 % (args.panel_condition, ", ".join(cond_bw)))
    panel = cond_bw[args.panel_condition]
    overlay = cond_bw.get(args.overlay_condition) if args.overlay_condition else None
    if len(panel) != len(cats):
        sys.exit("[ERROR] %d bw but %d categories" % (len(panel), len(cats)))

    L = []
    cmaps = ["Reds", "Blues", "Purples", "Greens"]
    for i, item in enumerate(args.cools.split(",")):
        label, path = item.split(":", 1)
        L += ["[mcc_%s]" % label,
              "file = %s" % path,
              "title = MCC %s" % label,
              "depth = %d" % args.depth,
              "transform = log1p",
              "colormap = %s" % cmaps[i % len(cmaps)],
              "height = 4",
              "file_type = hic_matrix",
              "show_masked_bins = false", ""]

    L += ["[spacer]", "height = 0.3", ""]
    bes_lab = args.bes_label or "BEARING BES -log10p"
    L += ["[bes]",
          "file = %s" % args.bes,
          "title = %s" % bes_lab,
          "color = #222222",
          "min_value = 0",
          "height = 2.5",
          "file_type = bigwig", ""]
    L += ["[spacer]", "height = 0.3", ""]

    for (label, color, _neg), bw in zip(cats, panel):
        tag = label.replace(" ", "_").replace("+", "pos").replace("-", "neg")
        if overlay is not None:
            ov = overlay[panel.index(bw)]
            L += ["[ov_%s]" % tag,
                  "file = %s" % ov,
                  "title = %s" % label,
                  "color = #cfcfcf",
                  "alpha = 0.6",
                  "height = 2",
                  "file_type = bigwig", ""]
            L += ["[%s]" % tag,
                  "file = %s" % bw,
                  "color = %s" % color,
                  "overlay_previous = share-y",
                  "file_type = bigwig", ""]
        else:
            L += ["[%s]" % tag,
                  "file = %s" % bw,
                  "title = %s" % label,
                  "color = %s" % color,
                  "height = 2",
                  "file_type = bigwig", ""]

    if args.genes:
        L += ["[genes]", "file = %s" % args.genes, "title = Tcrb",
              "file_type = bed", "height = 1.5", "fontsize = 9",
              "gene_rows = 2", ""]
    if args.features:
        L += ["[cbe]", "file = %s" % args.features, "title = CBE",
              "color = #d00000", "file_type = bed", "display = collapsed",
              "height = 0.6", "labels = true", "fontsize = 8", ""]
    L += ["[x-axis]", "fontsize = 10", ""]

    with open(args.out_ini, "w") as fh:
        fh.write("\n".join(L) + "\n")
    sys.stderr.write("[ini] %s  (%d matrices, %d panel tracks%s)\n"
                     % (args.out_ini, len(args.cools.split(",")), len(cats),
                        ", DN overlay" if overlay is not None else ""))
    sys.stderr.write("[run] pyGenomeTracks --tracks %s --region %s -o figure.pdf\n"
                     % (args.out_ini, args.region))


if __name__ == "__main__":
    main()
