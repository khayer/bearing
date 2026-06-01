#!/usr/bin/env python3
"""
montage_panels.py

Stitch already-rendered PNG panels into a labeled grid composite. General
fallback for layout gaps where the per-panel figures exist but no single script
emits the combined arrangement (e.g. Supp S9 per-replicate stack, or any
multi-panel composite assembled outside a deck tool).

For Supp S9 specifically: the recommended source of the 10-replicate stack is
the compare_qcat region plot at the new_wide region (it already renders all
samples). Use this montage only if you batch per-replicate panels separately and
need them stacked with row labels.

Usage:
  python3 montage_panels.py \\
      --panels DN_rep1.png DN_rep2.png DP_rep1.png ... S3T3_rep2.png \\
      --labels DN_rep1 DN_rep2 DP_rep1 ... S3T3_rep2 \\
      --ncols 1 --out suppS9_perreplicate_domain.png

ASCII-only. Requires Pillow.
"""

import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="+", required=True)
    ap.add_argument("--labels", nargs="*", default=None,
                    help="optional per-panel labels (same count as --panels)")
    ap.add_argument("--ncols", type=int, default=1)
    ap.add_argument("--pad", type=int, default=8, help="pixels between panels")
    ap.add_argument("--label-height", type=int, default=22)
    ap.add_argument("--bg", default="white")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        sys.exit("ERROR: Pillow required (pip install pillow)")

    if args.labels and len(args.labels) != len(args.panels):
        sys.exit("ERROR: --labels count must match --panels count")

    imgs = [Image.open(p).convert("RGB") for p in args.panels]
    cell_w = max(im.width for im in imgs)
    cell_h = max(im.height for im in imgs)
    lh = args.label_height if args.labels else 0
    ncols = max(1, args.ncols)
    nrows = (len(imgs) + ncols - 1) // ncols
    pad = args.pad

    W = ncols * cell_w + (ncols + 1) * pad
    H = nrows * (cell_h + lh) + (nrows + 1) * pad
    canvas = Image.new("RGB", (W, H), args.bg)
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 14)
    except Exception:
        font = ImageFont.load_default()

    for i, im in enumerate(imgs):
        r, c = divmod(i, ncols)
        x = pad + c * (cell_w + pad)
        y = pad + r * (cell_h + lh + pad)
        if args.labels:
            draw.text((x + 4, y + 2), args.labels[i], fill="black", font=font)
            y += lh
        # center the panel in its cell
        ox = x + (cell_w - im.width) // 2
        oy = y + (cell_h - im.height) // 2
        canvas.paste(im, (ox, oy))

    canvas.save(args.out)
    sys.stderr.write("wrote %s (%dx%d, %d panels)\n" %
                     (args.out, W, H, len(imgs)))


if __name__ == "__main__":
    main()
