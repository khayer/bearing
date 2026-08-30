#!/usr/bin/env python3
"""
bearing_extract_attribution.py  (v3 - adds --per-bin normalized concentration)

Return the per-track attribution numbers for the V1P/4C result WITHOUT uploading
the big file. Runs locally against the DN-vs-dV1P DIFFERENTIAL qcat
(compare_qcat.py --diff output) that INCLUDES the observed/expected virtual-4C
contact channel(s) as tracks.

v3 adds --per-bin: report each track's |diff| PER BIN (summed |contribution| /
bin count) so regions of very different sizes are comparable, plus a
concentration summary for each viewpoint/contact channel: its per-bin peak
region vs the out-of-TAD floor, and the ratio. Use it to watch a reciprocal
viewpoint (e.g. 4C Trbv13 at DJ_RC) pull away from its genome-wide floor as you
fix normalization. A clean loop-specific channel (e.g. 4C RC at Trbv1) shows a
large ratio; a noisy one shows a ratio near 1.

v2 (kept): if the input is bgzip+tabix indexed (a .tbi/.csi sits next to it),
only the requested windows are fetched via tabix (CLI), then pysam, then a plain
gzip stream (slow) as a last resort.

ASCII only. Standard library + (optional) tabix/pysam. No repo import required.

USAGE
  python3 bearing_extract_attribution.py diff_DN_vs_V1P.qcat.bgz --describe
  python3 bearing_extract_attribution.py diff_DN_vs_V1P.qcat.bgz \
      --track-names "RNAseq +,RNAseq -,CTCF,Cohesin,NIPBL,H3K27ac,4C RC,4C Trbv13" \
      --contact-name "4C" --regions v1p_controls.bed --per-bin
"""
import argparse, gzip, io, json, os, re, shutil, subprocess, sys

DEFAULT_REGIONS = [
    ("chr6", 40880000, 40905000, "Trbv1"),
    ("chr6", 41040000, 41290000, "Vbeta_cluster"),
    ("chr6", 41313087, 41486888, "NEGCTL"),
    ("chr6", 41500000, 41551500, "DJ_RC"),
    ("chr6", 41552000, 41562000, "Ebeta_Trbv31_CBE3"),
]


# ---------------- fast region-restricted line iterator ----------------------
def is_indexed(path):
    return (path.endswith(".bgz") or path.endswith(".gz")) and (
        os.path.exists(path + ".tbi") or os.path.exists(path + ".csi"))


def iter_region_tabix(path, chrom, start, end):
    """Yield lines for chrom:start-end using the `tabix` CLI (1-based inclusive)."""
    reg = f"{chrom}:{start+1}-{end}"
    proc = subprocess.Popen(["tabix", path, reg], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        yield line.rstrip("\n")
    proc.stdout.close()
    proc.wait()


def iter_region_pysam(tbx, chrom, start, end):
    for row in tbx.fetch(chrom, start, end):
        yield row


def iter_region_stream(path, chrom, start, end):
    """Fallback: gzip stream, keep only rows in the region. Slow (full scan)."""
    op = io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace") \
        if path.endswith((".gz", ".bgz")) else open(path, encoding="utf-8", errors="replace")
    with op as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            p = line.split("\t", 3)
            if len(p) < 2:
                continue
            if p[0] != chrom:
                continue
            try:
                s = int(p[1])
            except ValueError:
                continue
            if start <= s < end:
                yield line.rstrip("\n")


def make_region_iter(path):
    """Return a function(chrom,start,end)->iterable of lines, using the fastest path."""
    if is_indexed(path) and shutil.which("tabix"):
        return lambda c, s, e: iter_region_tabix(path, c, s, e), "tabix"
    try:
        import pysam  # noqa
        if is_indexed(path):
            tbx = pysam.TabixFile(path)
            return lambda c, s, e: iter_region_pysam(tbx, c, s, e), "pysam"
    except Exception:
        pass
    return lambda c, s, e: iter_region_stream(path, c, s, e), "gzip-stream(SLOW)"


# ---------------- qcat payload parsing --------------------------------------
def parse_qcat_payload(payload):
    """From 'id:N,qcat:[[s,t],...],raw:[...]' return list of (track_id:int, score:float)."""
    m = payload.find("qcat:")
    if m < 0:
        return None
    rest = payload[m + len("qcat:"):]
    cut = rest.find(",raw:")
    arr = rest[:cut] if cut >= 0 else rest
    arr = arr.strip().rstrip(",")
    try:
        pairs = json.loads(arr)
    except Exception:
        return None
    out = []
    for p in pairs:
        try:
            out.append((int(p[1]), float(p[0])))
        except Exception:
            continue
    return out


def row_to_pairs(line, layout, ncol_names=None):
    parts = line.split("\t")
    if len(parts) < 2:
        return None, None, None
    chrom = parts[0]
    try:
        start = int(parts[1])
    except ValueError:
        return None, None, None
    if layout == "qcat":
        pairs = parse_qcat_payload(parts[3] if len(parts) > 3 else "")
    else:
        pairs = []
        for j in range(3, len(parts)):
            try:
                pairs.append((j - 3, float(parts[j])))
            except ValueError:
                pass
    return chrom, start, pairs


def sniff_layout(region_iter):
    """Grab one row from the Trbv1 window to decide qcat vs tsv and show a sample."""
    for line in region_iter("chr6", 40880000, 40905000):
        return ("qcat" if "qcat:" in line else "tsv"), line
    for line in region_iter("chr6", 40790000, 41690000):
        return ("qcat" if "qcat:" in line else "tsv"), line
    return None, None


def main():
    ap = argparse.ArgumentParser(description="Per-track (contact-channel) attribution for V1P/4C.")
    ap.add_argument("diffqcat")
    ap.add_argument("--track-names", help="comma list in index order")
    ap.add_argument("--contact-name", help="name (or substring) of the contact channel track")
    ap.add_argument("--regions", help="BED override: chrom start end name")
    ap.add_argument("--describe", action="store_true", help="print access path, layout, a sample row, track ids/means in Trbv1, then exit")
    ap.add_argument("--per-bin", action="store_true",
                    help="also print |diff|-per-bin matrix and a per-viewpoint concentration summary (peak region vs out-of-TAD floor)")
    ap.add_argument("--floor-substr", default="out",
                    help="regions whose NAME contains this substring form the out-of-TAD floor set for --per-bin (default: 'out', matches outTAD/out4C)")
    ap.add_argument("--viewpoint-names",
                    help="comma list of track names (or substrings) to treat as viewpoint/contact channels in the --per-bin summary; default: auto-detect (4c/contact/oe/hic)")
    ap.add_argument("--anchor-map",
                    help="pin the concentration anchor per channel: 'CHAN:REGION,CHAN:REGION' "
                         "(e.g. '4C RC:Trbv1_anchor,4C Trbv13:DJ_RC_anchor'). CHAN and REGION "
                         "match by substring. For mapped channels the ratio is anchor/floor at the "
                         "named region instead of the auto per-bin peak; unmapped channels still use auto-peak.")
    args = ap.parse_args()

    region_iter, access = make_region_iter(args.diffqcat)
    layout, sample = sniff_layout(region_iter)
    if layout is None:
        sys.exit("ERROR: no rows returned for chr6 Tcrb window. Check the chrom name (chr6 vs 6) and that the file is indexed.")

    names = args.track_names.split(",") if args.track_names else None

    # ---- detect 0-based vs 1-based track ids by pre-scanning the Trbv1 window ----
    prescan_ids = set()
    for line in region_iter("chr6", 40790000, 41690000):
        if not line or line[0] == "#":
            continue
        _c, _s, _pairs = row_to_pairs(line, layout)
        if _pairs:
            for tid, _ in _pairs:
                prescan_ids.add(tid)
        if len(prescan_ids) >= 8:
            break
    base = 0 if (0 in prescan_ids or not prescan_ids) else 1

    regions = DEFAULT_REGIONS
    if args.regions:
        regions = []
        with open(args.regions) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    p = line.split()
                    regions.append((p[0], int(p[1]), int(p[2]), p[3] if len(p) > 3 else f"{p[0]}:{p[1]}"))

    acc = {r[3]: {} for r in regions}
    nbins = {r[3]: 0 for r in regions}
    top_contact = {r[3]: 0 for r in regions}
    seen_tracks = set()

    def track_name(tid):
        j = tid - base
        if names and 0 <= j < len(names):
            return names[j].strip()
        return f"track{tid}"

    for (rc, rs, re_, name) in regions:
        for line in region_iter(rc, rs, re_):
            if not line or line[0] == "#":
                continue
            chrom, start, pairs = row_to_pairs(line, layout)
            if chrom is None or pairs is None:
                continue
            if not (chrom == rc and rs <= start < re_):
                continue
            nbins[name] += 1
            best_tid, best_abs = None, -1.0
            for tid, sc in pairs:
                seen_tracks.add(tid)
                acc[name][tid] = acc[name].get(tid, 0.0) + abs(sc)
                if abs(sc) > best_abs:
                    best_abs, best_tid = abs(sc), tid
            cid = contact_id(names, args.contact_name, seen_tracks, track_name)
            if best_tid is not None and cid is not None and best_tid == cid:
                top_contact[name] += 1

    if args.describe:
        print(f"# access path: {access}")
        print(f"# layout: {layout}")
        print(f"# sample row (first Tcrb row):\n   {sample[:240]}")
        print(f"# track ids observed in Tcrb: {sorted(seen_tracks)}")
        print("# id -> name (as resolved):")
        for tid in sorted(seen_tracks):
            print(f"   {tid} -> {track_name(tid)}")
        print("\n# Trbv1 per-track summed |contribution|:")
        for tid in sorted(acc.get('Trbv1', {}), key=lambda t: -acc['Trbv1'][t]):
            print(f"   {track_name(tid):14s} {acc['Trbv1'][tid]:.4f}")
        print("\nNow re-run with --track-names ... --contact-name NAME (the 4C/contact track).")
        return

    cid = contact_id(names, args.contact_name, seen_tracks, track_name)
    print(f"# bearing_extract_attribution.py  file={args.diffqcat}  access={access}  layout={layout}")
    print(f"# contact channel resolved to: "
          f"{track_name(cid) if cid is not None else 'UNRESOLVED (pass --contact-name)'}")
    print()
    for (rc, rs, re_, name) in regions:
        tot = sum(acc[name].values())
        print(f"== {name} ==  bins={nbins[name]}  total|diff|={tot:.4f}")
        if tot == 0:
            print("   (no scored signal in region)\n")
            continue
        for tid in sorted(acc[name], key=lambda t: -acc[name][t]):
            frac = 100.0 * acc[name][tid] / tot
            star = "  <== CONTACT" if (cid is not None and tid == cid) else ""
            print(f"   {track_name(tid):14s} {acc[name][tid]:9.4f}  {frac:5.1f}%{star}")
        if cid is not None and cid in acc[name]:
            cfrac = 100.0 * acc[name][cid] / tot
            print(f"   -> contact channel = {cfrac:.1f}% of |diff|; "
                  f"top contributor in {top_contact[name]}/{nbins[name]} bins")
        print()

    if args.per_bin:
        print_per_bin(regions, acc, nbins, seen_tracks, track_name,
                      names, args.viewpoint_names, args.contact_name, args.floor_substr,
                      args.anchor_map)

    tr = acc.get("Trbv1", {})
    if cid is not None and tr and cid in tr:
        cfrac = 100.0 * tr[cid] / sum(tr.values())
        print("# READY-TO-PASTE for R3 (verify against numbers above):")
        print(f'   "At Trbv1 the differential is carried predominantly by the contact channel '
              f'({cfrac:.0f}% of the summed per-track differential; top contributor in '
              f'{top_contact["Trbv1"]}/{nbins["Trbv1"]} bins), with the linear chromatin assays low -- '
              f'the loop loss is largely invisible on linear tracks and surfaces once the contact '
              f'channel is added to the panel."')


def per_bin_value(acc_region, nbins_region, tid):
    if nbins_region <= 0:
        return 0.0
    return acc_region.get(tid, 0.0) / nbins_region


def viewpoint_ids(names, viewpoint_names, contact_name, seen_tracks, track_name):
    """Return sorted track ids to summarize as viewpoint/contact channels."""
    subs = []
    if viewpoint_names:
        subs = [s.strip().lower() for s in viewpoint_names.split(",") if s.strip()]
    elif contact_name:
        subs = [contact_name.strip().lower()]
    ids = []
    if subs:
        for tid in sorted(seen_tracks):
            nm = track_name(tid).lower()
            if any(s in nm for s in subs):
                ids.append(tid)
    if not ids:  # auto-detect
        for tid in sorted(seen_tracks):
            if any(k in track_name(tid).lower() for k in ("4c", "contact", "oe", "hic")):
                ids.append(tid)
    return ids


def parse_anchor_map(anchor_map, order, track_name, seen_tracks):
    """Return {tid: region_name} from 'CHAN:REGION,...' (substring match on both)."""
    out = {}
    if not anchor_map:
        return out
    for item in anchor_map.split(","):
        if ":" not in item:
            continue
        chan_key, reg_key = item.split(":", 1)
        chan_key = chan_key.strip().lower()
        reg_key = reg_key.strip().lower()
        if not chan_key or not reg_key:
            continue
        reg = None
        for n in order:                       # exact (case-insensitive) first
            if n.lower() == reg_key:
                reg = n
                break
        if reg is None:                        # then substring
            for n in order:
                if reg_key in n.lower():
                    reg = n
                    break
        if reg is None:
            sys.stderr.write(f"# --anchor-map: no region matches '{reg_key}'; known: {order}\n")
            continue
        for tid in sorted(seen_tracks):        # map every channel matching the substring
            if chan_key in track_name(tid).lower():
                out[tid] = reg
    return out


def print_per_bin(regions, acc, nbins, seen_tracks, track_name,
                  names, viewpoint_names, contact_name, floor_substr, anchor_map=None):
    order = [r[3] for r in regions]
    # rank tracks by total |diff| across all regions for row order
    tot_by_tid = {}
    for tid in seen_tracks:
        tot_by_tid[tid] = sum(acc[n].get(tid, 0.0) for n in order)
    tracks_sorted = sorted(seen_tracks, key=lambda t: -tot_by_tid[t])

    print("# ---- per-bin normalized attribution: |diff| per bin (summed |contrib| / bins) ----")
    colw = 11
    header = "  ".join(f"{n[:colw]:>{colw}s}" for n in order)
    print(f"{'track':14s}  {header}")
    print(f"{'bins ->':14s}  " + "  ".join(f"{nbins[n]:>{colw}d}" for n in order))
    for tid in tracks_sorted:
        row = "  ".join(f"{per_bin_value(acc[n], nbins[n], tid):>{colw}.3f}" for n in order)
        print(f"{track_name(tid):14s}  {row}")
    print()

    floor_key = (floor_substr or "").lower()
    floor_regions = [n for n in order if floor_key and floor_key in n.lower()]
    vids = viewpoint_ids(names, viewpoint_names, contact_name, seen_tracks, track_name)
    pinned = parse_anchor_map(anchor_map, order, track_name, seen_tracks)

    print("# ---- viewpoint-channel concentration (per-bin anchor vs out-of-TAD floor) ----")
    if floor_regions:
        print(f"# floor set (name contains '{floor_substr}'): {', '.join(floor_regions)}")
    else:
        print(f"# floor set: NONE matched '{floor_substr}'. Pass --floor-substr or name your out-of-TAD "
              f"regions with that substring; ratio uses the max non-anchor region instead.")
    if pinned:
        print("# anchors pinned via --anchor-map: "
              + ", ".join(f"{track_name(t)}->{pinned[t]}" for t in sorted(pinned)))
    print("# target for a clean loop-specific channel: ratio >> 5x (RC/Trbv1 is the reference).")
    print()
    for tid in vids:
        pb = {n: per_bin_value(acc[n], nbins[n], tid) for n in order}
        auto_region = max(order, key=lambda n: pb[n])
        if tid in pinned:
            anchor_region = pinned[tid]
            pinned_note = "" if anchor_region == auto_region else \
                f"   (auto per-bin peak is elsewhere: {auto_region} = {pb[auto_region]:.3f})"
        else:
            anchor_region = auto_region
            pinned_note = ""
        anchor_val = pb[anchor_region]
        # floor excludes the anchor region so a pinned anchor inside the floor set is not its own floor
        if floor_regions:
            floor_pairs = [(n, pb[n]) for n in floor_regions if n != anchor_region]
        else:
            floor_pairs = [(n, pb[n]) for n in order if n != anchor_region]
        floor_max_region, floor_max = (max(floor_pairs, key=lambda x: x[1])
                                       if floor_pairs else ("-", 0.0))
        floor_mean = (sum(v for _, v in floor_pairs) / len(floor_pairs)) if floor_pairs else 0.0
        ratio = (anchor_val / floor_max) if floor_max > 0 else float("inf")
        ratio_s = f"{ratio:.2f}x" if ratio != float("inf") else "inf (floor=0)"
        label = "anchor (pinned)" if tid in pinned else "anchor (auto-peak)"
        print(f"{track_name(tid)}")
        print("   per-bin by region: " + "  ".join(f"{n}={pb[n]:.3f}" for n in order))
        print(f"   {label}: {anchor_region} = {anchor_val:.3f}")
        if pinned_note:
            print(pinned_note)
        print(f"   out-of-TAD floor:  max={floor_max:.3f} ({floor_max_region}), mean={floor_mean:.3f}")
        print(f"   concentration:     anchor/floor_max = {ratio_s}")
        print()


def contact_id(names, contact_name, seen_tracks, track_name):
    key = (contact_name or "").lower()
    if key:
        for tid in sorted(seen_tracks):
            if key in track_name(tid).lower():
                return tid
        return None
    for tid in sorted(seen_tracks):
        if any(k in track_name(tid).lower() for k in ("4c", "contact", "oe", "hic")):
            return tid
    return max(seen_tracks) if seen_tracks else None


if __name__ == "__main__":
    main()
