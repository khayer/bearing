import pyBigWig, os

BWROOT = ".."   # so BWROOT/bigwigs/... matches the sheet; adjust if needed

CTCF = {"DN_1":"bigwigs/CTCF/R1KO_CTCF0821CR_DN_1.bw",
        "DN_2":"bigwigs/CTCF/R1KO_CTCF0905CR_DN_2.bw",
        "TKO_1":"bigwigs/CTCF/HTS_SO85_01_RCTKO_CTCF1030CR_DN_1_S1.bw",
        "TKO_2":"bigwigs/CTCF/HTS_SO85_09_RCTKO_CTCF1206CR_DN_2_S9.bw"}
RAD21 = {"DN_1":"bigwigs/Rad21/R1KO_Rad210712CR_DN_1.bw",
         "DN_2":"bigwigs/Rad21/R1KO_Rad210712CR_DN_2.bw",
         "TKO_1":"bigwigs/Rad21/HTS_SO85_02_RCTKO_Rad211030CR_DN_1_S2.bw",
         "TKO_2":"bigwigs/Rad21/HTS_SO85_06_RCTKO_Rad211206CR_DN_2_S6.bw"}
NIPBL = {"DN_1":"bigwigs/NIPBL/R1KO_Nipbl0816CR_DN_1.bw",
         "DN_2":"bigwigs/NIPBL/R1KO_Nipbl0816CR_DN_2.bw",
         "TKO_1":"bigwigs/NIPBL/HTS_SO85_03_RCTKO_Nipbl1030CR_DN_1_S3.bw",
         "TKO_2":"bigwigs/NIPBL/HTS_SO85_07_RCTKO_Nipbl1206CR_DN_2_S7.bw"}

# 200 bp covering bin for each element (bin start = floor(pos/200)*200)
elts = [("5'PC  (DEL)", 41505538), ("CBE_1 (DEL)", 41530176),
        ("CBE_3 (DEL)", 41556571), ("CBE_2 (kept)", 41532266)]

def binmean(path, s, e):
    p = os.path.join(BWROOT, path)
    if not os.path.exists(p): return "MISSING"
    bw = pyBigWig.open(p)
    v = bw.stats("chr6", s, e, type="mean")[0]
    bw.close()
    return 0.0 if v is None else round(v, 3)

for label, mark in [("CTCF", CTCF), ("RAD21", RAD21), ("NIPBL", NIPBL)]:
    print("\n=== %s (mean over covering 200 bp bin) ===" % label)
    print("%-14s %8s %8s | %8s %8s" % ("element","DN_1","DN_2","TKO_1","TKO_2"))
    for name, pos in elts:
        b = (pos // 200) * 200
        row = [binmean(mark[s], b, b+200) for s in ("DN_1","DN_2","TKO_1","TKO_2")]
        print("%-14s %8s %8s | %8s %8s" % (name, row[0], row[1], row[2], row[3]))