#!/usr/bin/env python3
"""
hic_probe.py
Quick diagnostic: print available resolutions, chromosomes,
and which normalizations are present in a .hic file.

Usage:
    python hic_probe.py myfile.hic
    python hic_probe.py myfile.hic chr6   # also probe one chrom
"""
import sys

def probe(hic_path, test_chrom=None):
    try:
        import hicstraw
    except ImportError:
        sys.exit("pip install hic-straw")

    hic = hicstraw.HiCFile(hic_path)

    print("=" * 60)
    print("File:", hic_path)
    print("=" * 60)

    resolutions = hic.getResolutions()
    print("\nAvailable resolutions (bp):")
    for r in resolutions:
        print(f"  {r:>10,}")

    chroms = hic.getChromosomes()
    print("\nChromosomes in file:")
    for c in chroms:
        print(f"  {c.name:15s}  length={c.length:,}")

    print(dir(hic))
    #norms = hic.normalizations()
    #print("\nNormalization in file:")
    #for n in norms:
    #    print(f"  {n}")

    mzd = hic.getMatrixZoomData(chroms[1].name, chroms[1].name, "observed", "NONE", "BP", resolutions[0])

    # 3. Use the mzd object to see available normalizations
    # In many versions, the available norms are found here:
    print(f"Available Normalizations: {hic.getNormalizationOptions() if hasattr(hic, 'getNormalizationOptions') else 'Check documentation for version compatibility'}")
    
    if test_chrom is None:
        # Pick the first non-ALL chromosome
        candidates = [c.name for c in chroms
                      if c.name.lower() not in ("all", "assembly")]
        test_chrom = candidates[0] if candidates else None

    if test_chrom and resolutions:
        print(f"\nProbing normalization x resolution for {test_chrom}:")
        norms = ["KR", "SCALE", "VC", "VC_SQRT", "NONE"]
        for res in resolutions:
            for norm in norms:
                try:
                    chrom_str = test_chrom.replace("chr", "")
                    mzd = hic.getMatrixZoomData(
                        chrom_str, chrom_str, "observed", norm, "BP", res
                    )
                    # Probe a tiny region at the start of the chromosome
                    length = next(
                        c.length for c in chroms
                        if c.name == test_chrom or c.name == chrom_str
                    )
                    probe_end = min(length, res * 10)
                    recs = mzd.getRecords(0, probe_end, 0, probe_end)
                    n = len(recs)
                    status = "OK (n=" + str(n) + " records)"
                except Exception as e:
                    status = "FAIL: " + str(e)[:60]
                print(f"  {res:>10,} bp  {norm:<8s}  {status}")

    print("\nDone.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python hic_probe.py file.hic [chrom]")
    probe(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
