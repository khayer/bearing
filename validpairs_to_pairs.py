#!/usr/bin/env python3
"""
validpairs_to_pairs.py -- turn a name-sorted Hi-C valid-pairs BAM into a
sortable, upper-triangle pairs stream for `cooler cload pairs`.

Reads `samtools view` output on stdin (SAM records, name-sorted so the two mates
of a template are consecutive). For each template it takes the two PRIMARY
alignments (skipping secondary/supplementary, FLAG & 0x900) and emits one
contact, using each mate's own RNAME/POS (does NOT rely on RNEXT/PNEXT, which
HiC-Pro valid-pairs BAMs do not always set).

Output columns (tab-separated), designed so a numeric sort gives cooler the
block-by-chrom, then-by-position order it wants:

    idx1  idx2  readid  chrom1  pos1  chrom2  pos2

idx1/idx2 are the chromosome ranks from the chrom.sizes file (sort keys only;
cooler ignores them via -c/-p column indices). Each contact is oriented so
(idx1,pos1) <= (idx2,pos2), i.e. upper triangle.

USAGE
    samtools view sample_valid_pairs.bam \\
      | python3 validpairs_to_pairs.py mm10.chrom.sizes \\
      | sort -k1,1n -k2,2n -k5,5n -k7,7n \\
      | bgzip > sample.pairs.gz
    cooler cload pairs -c1 4 -p1 5 -c2 6 -p2 7 \\
      mm10.chrom.sizes:2000 sample.pairs.gz sample.2000.cool

ASCII only.
"""
import sys


def load_order(path):
    order = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line:
                order[line.split()[0]] = i
    return order


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: samtools view BAM | %s chrom.sizes\n" % sys.argv[0])
        sys.exit(2)
    order = load_order(sys.argv[1])
    out = sys.stdout.write

    cur = None
    grp = []
    n_emit = 0
    n_skip = 0

    def flush(grp):
        nonlocal n_emit, n_skip
        if len(grp) != 2:
            if grp:
                n_skip += 1
            return
        (c1, p1), (c2, p2) = grp
        i1 = order.get(c1)
        i2 = order.get(c2)
        if i1 is None or i2 is None:
            n_skip += 1
            return
        if (i1, p1) > (i2, p2):
            i1, p1, c1, i2, p2, c2 = i2, p2, c2, i1, p1, c1
        out("%d\t%d\t.\t%s\t%d\t%s\t%d\n" % (i1, i2, c1, p1, c2, p2))
        n_emit += 1

    for line in sys.stdin:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & 0x900:          # secondary / supplementary
            continue
        qname, rname, pos = f[0], f[2], f[3]
        if qname != cur:
            flush(grp)
            cur, grp = qname, []
        if rname != "*":
            grp.append((rname, int(pos)))
    flush(grp)

    sys.stderr.write("[pairs] emitted %d contacts, skipped %d incomplete templates\n"
                     % (n_emit, n_skip))


if __name__ == "__main__":
    main()
