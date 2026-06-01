#!/usr/bin/env python3
"""
stage_inputs.py

Download remote inputs referenced in a BEARING sample sheet to a local cache
using wget, ONLY if not already present (idempotent), and emit a "localized"
sample sheet whose paths point at the cache. Use when the cluster has no
aws/mount-s3: host the .cool / .bw / pca1 files behind HTTPS (public bucket URLs
or presigned S3 URLs) and let wget pull them.

  python3 stage_inputs.py \\
      --sheet workflow/config/samples.tsv \\
      --resolutions 10000 25000 100000 250000 500000 \\
      --cache-dir resources/staged \\
      --out-sheet workflow/config/samples.local.tsv

Then run the workflow against the localized sheet:
  snakemake ... --config samples_sheet=workflow/config/samples.local.tsv

Notes
- Only http:// and https:// are downloaded. s3:// is NOT supported by wget;
  convert to an https bucket URL or a presigned URL first (the script errors
  with guidance if it sees s3://).
- Local paths in the sheet are passed through unchanged.
- {res} placeholders in cool/insul/pca1 are expanded over --resolutions; each
  expansion is downloaded, and the localized cell keeps {res} pointing at the
  cache so the workflow substitutes it per resolution as before.
- Idempotent: a file already present and non-empty in the cache is skipped.
  Downloads are atomic (.part then rename) so an interrupted run re-downloads
  cleanly rather than leaving a truncated file.

ASCII-only.
"""

import argparse
import hashlib
import os
import subprocess
import sys

import pandas as pd

URL_COLS_SINGLE = ["bw"]          # comma-separated list of files
URL_COLS_RES = ["cool", "insul", "pca1"]  # may contain {res}


def is_remote(tok):
    return tok.startswith("http://") or tok.startswith("https://")


def cache_path(url_or_pattern, cache_dir):
    """Deterministic cache location: cache_dir/<sha1(dir)>/<basename>.
    Keeps basenames (incl. any {res}) readable while avoiding collisions
    between different remote directories that share a basename."""
    base = url_or_pattern.rsplit("/", 1)[-1]
    parent = url_or_pattern[: -len(base)] if base else url_or_pattern
    h = hashlib.sha1(parent.encode("utf-8")).hexdigest()[:10]
    return os.path.join(cache_dir, h, base)


def wget_one(url, dest, dry_run=False):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        sys.stderr.write("  skip (present): %s\n" % dest)
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    sys.stderr.write("  get : %s\n        -> %s\n" % (url, dest))
    if dry_run:
        return True
    tmp = dest + ".part"
    try:
        subprocess.run(["wget", "-q", "-O", tmp, url], check=True)
        os.replace(tmp, dest)
        return True
    except subprocess.CalledProcessError as e:
        if os.path.exists(tmp):
            os.remove(tmp)
        sys.stderr.write("  FAILED (%s): %s\n" % (e.returncode, url))
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out-sheet", required=True)
    ap.add_argument("--cache-dir", default="resources/staged")
    ap.add_argument("--resolutions", nargs="*", type=int,
                    default=[10000, 25000, 100000, 250000, 500000])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    df = pd.read_csv(args.sheet, sep="\t", dtype=str, comment="#").fillna("")
    n_get, n_fail = 0, 0

    def localize_token(tok, expand_res):
        nonlocal n_get, n_fail
        tok = tok.strip()
        if not tok:
            return tok
        if tok.startswith("s3://"):
            sys.exit("ERROR: s3:// not supported by wget. Convert to an https "
                     "bucket URL or a presigned URL: %s" % tok)
        if not is_remote(tok):
            return tok  # local path, leave alone
        local = cache_path(tok, args.cache_dir)
        if expand_res and "{res}" in tok:
            for res in args.resolutions:
                ok = wget_one(tok.replace("{res}", str(res)),
                              local.replace("{res}", str(res)), args.dry_run)
                n_get += 1
                n_fail += 0 if ok else 1
        else:
            ok = wget_one(tok, local, args.dry_run)
            n_get += 1
            n_fail += 0 if ok else 1
        return local

    for col in URL_COLS_SINGLE:
        if col in df.columns:
            df[col] = df[col].map(
                lambda cell: ",".join(localize_token(t, False)
                                      for t in cell.split(",")) if cell else cell)
    for col in URL_COLS_RES:
        if col in df.columns:
            df[col] = df[col].map(lambda cell: localize_token(cell, True))

    df.to_csv(args.out_sheet, sep="\t", index=False)
    sys.stderr.write("\nlocalized sheet: %s\n" % args.out_sheet)
    sys.stderr.write("downloads attempted: %d, failed: %d\n" % (n_get, n_fail))
    if n_fail:
        sys.exit("ERROR: %d downloads failed (see above)" % n_fail)


if __name__ == "__main__":
    main()
