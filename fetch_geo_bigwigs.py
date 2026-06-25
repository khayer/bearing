#!/usr/bin/env python3
"""
fetch_geo_bigwigs.py -- pull processed bigWig tracks from a GEO series and
organize them into a BEARING-ready layout.

Reads the GEO family SOFT file (one record per GSM, listing Sample_title and
Sample_supplementary_file_* URLs), keeps the bigWig supplementary files,
classifies each by condition / assay / strand from the sample title, and writes:
  - <out>/geo_bigwig_manifest.tsv   gsm, title, condition, assay, strand, url
  - <out>/download_bigwigs.sh       wget commands (resumable) into <out>/bw/

Run it once to inspect the manifest (titles -> classification), correct the
classifier rules or the manifest if needed, then run the download script (or
pass --download to fetch here). No download happens without --download or
running the emitted script.

This runs where there is internet (your cluster), NOT in a sandbox. Uses only
the Python standard library.

USAGE
  python3 fetch_geo_bigwigs.py --gse GSE296315 --out geo_v1 \
      --conditions WT,Rag1,R1KO,513TKO,RCTKO,V1PKO,V1TxS \
      --assays CTCF,RAD21,NIPBL,RNA,GRO,H3K27ac,H3K4me3,Pol2
  # inspect geo_v1/geo_bigwig_manifest.tsv, then:
  bash geo_v1/download_bigwigs.sh

ASCII only.
"""
import argparse
import gzip
import io
import os
import re
import sys
import urllib.request


def soft_url(gse):
    stub = gse[:-3] + "nnn"
    return ("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/soft/%s_family.soft.gz"
            % (stub, gse, gse))


def fetch_soft(gse):
    url = soft_url(gse)
    sys.stderr.write("[soft] %s\n" % url)
    with urllib.request.urlopen(url, timeout=120) as r:
        raw = r.read()
    return gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode("utf-8", "replace")


def parse_samples(text):
    """Yield dicts: gsm, title, strategy, [supp urls]."""
    cur = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                yield cur
            cur = {"gsm": line.split("=", 1)[1].strip(), "title": "",
                   "strategy": "", "urls": []}
        elif cur is None:
            continue
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_library_strategy"):
            cur["strategy"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_supplementary_file"):
            u = line.split("=", 1)[1].strip()
            if re.search(r"\.(bw|bigwig|bigWig)$", u, re.I):
                cur["urls"].append(u)
    if cur:
        yield cur


# classification rules (title is messy and lab-specific; tweak as needed)
ASSAY_PATTERNS = [
    ("CTCF", r"ctcf"),
    ("RAD21", r"rad21|cohesin"),
    ("NIPBL", r"nipbl"),
    ("H3K27ac", r"h3k27ac|k27ac"),
    ("H3K4me3", r"h3k4me3|k4me3"),
    ("Pol2Ser2", r"pol\s*2.*ser2|polii.*ser2|ser2|rnapol"),
    ("GRO", r"gro[-_ ]?seq|gro\b"),
    ("RNA", r"rna[-_ ]?seq|rnaseq|\brna\b|total\s*rna"),
]
STRAND_PATTERNS = [
    ("plus", r"\bplus\b|\bsense\b|\bfwd\b|forward|\+|pos\b"),
    ("minus", r"\bminus\b|antisense|\brev\b|reverse|\-|neg\b"),
]


def classify(title, conditions):
    t = title.lower()
    assay = next((a for a, pat in ASSAY_PATTERNS if re.search(pat, t)), "NA")
    strand = "NA"
    if assay in ("RNA", "GRO"):
        strand = next((s for s, pat in STRAND_PATTERNS if re.search(pat, t)), "NA")
    cond = "NA"
    for c in conditions:
        if re.search(re.escape(c.lower()), t):
            cond = c
            break
    rep = "NA"
    m = re.search(r"rep[\s_]?(\d+)|[_\s]r(\d)\b", t)
    if m:
        rep = m.group(1) or m.group(2)
    return cond, assay, strand, rep


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gse", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--conditions", default="",
                    help="comma list of condition tokens to match in titles")
    ap.add_argument("--assays", default="",
                    help="comma list of assay names to KEEP (default: keep all classified)")
    ap.add_argument("--soft-file", default=None,
                    help="use a local family.soft(.gz) instead of downloading")
    ap.add_argument("--download", action="store_true",
                    help="download here via urllib (else just emit download_bigwigs.sh)")
    args = ap.parse_args()

    os.makedirs(os.path.join(args.out, "bw"), exist_ok=True)
    conditions = [c.strip() for c in args.conditions.split(",") if c.strip()]
    keep_assays = set(a.strip() for a in args.assays.split(",") if a.strip())

    if args.soft_file:
        op = gzip.open if args.soft_file.endswith(".gz") else open
        text = op(args.soft_file, "rt", encoding="utf-8", errors="replace").read()
    else:
        text = fetch_soft(args.gse)

    rows = []
    for s in parse_samples(text):
        if not s["urls"]:
            continue
        cond, assay, strand, rep = classify(s["title"], conditions)
        if keep_assays and assay not in keep_assays:
            continue
        for u in s["urls"]:
            rows.append([s["gsm"], s["title"], cond, assay, strand, rep,
                         s["strategy"], u])

    man = os.path.join(args.out, "geo_bigwig_manifest.tsv")
    with open(man, "w") as fh:
        fh.write("\t".join(["gsm", "title", "condition", "assay", "strand",
                            "rep", "strategy", "url"]) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    sys.stderr.write("[manifest] %s (%d bigwigs)\n" % (man, len(rows)))

    # emit a resumable wget script with informative local names
    dl = os.path.join(args.out, "download_bigwigs.sh")
    with open(dl, "w") as fh:
        fh.write("#!/usr/bin/env bash\nset -euo pipefail\ncd \"$(dirname \"$0\")/bw\"\n")
        for gsm, title, cond, assay, strand, rep, strat, url in rows:
            base = "_".join([x for x in [cond, assay, strand, "rep" + rep if rep != "NA" else "", gsm]
                             if x and x != "NA"])
            ext = ".bw" if url.lower().endswith(".bw") else ".bigwig"
            fh.write('wget -c -O "%s%s" "%s"\n' % (base, ext, url))
    os.chmod(dl, 0o755)
    sys.stderr.write("[script] %s\n" % dl)

    n_na = sum(1 for r in rows if r[2] == "NA" or r[3] == "NA")
    if n_na:
        sys.stderr.write("[warn] %d/%d rows have NA condition or assay - "
                         "check titles and tweak classifier before downloading\n"
                         % (n_na, len(rows)))

    if args.download:
        for gsm, title, cond, assay, strand, rep, strat, url in rows:
            base = "_".join([x for x in [cond, assay, strand,
                             "rep" + rep if rep != "NA" else "", gsm] if x and x != "NA"])
            ext = ".bw" if url.lower().endswith(".bw") else ".bigwig"
            dst = os.path.join(args.out, "bw", base + ext)
            if os.path.exists(dst):
                continue
            sys.stderr.write("[get] %s\n" % base)
            urllib.request.urlretrieve(url.replace("ftp://", "https://"), dst)


if __name__ == "__main__":
    main()
