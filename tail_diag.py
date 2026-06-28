import numpy as np

def load(path, lo=40800000, hi=41650000):
    rows=[]
    with open(path) as f:
        next(f)
        for ln in f:
            c=ln.split("\t")
            if c[0]=="chr6" and int(c[1])>=lo and int(c[2])<=hi:
                rows.append((int(c[1]), abs(float(c[4]))))
    return rows

txs=load("mcc_results_v1txs/binwise_DN_vs_V1TxS.tsv")
v1p=load("mcc_results_v1_4track/binwise_DN_vs_V1P.tsv")
rct=load("mcc_results_v1_4track/binwise_DN_vs_RCTKO.tsv")

null=np.array([v for _,v in txs])
p95=np.percentile(null,95); p99=np.percentile(null,99)
print("V1TxS null:  p50=%.4f  p95=%.4f  p99=%.4f  (n=%d)"%(np.percentile(null,50),p95,p99,len(null)))
print()
for tag,rows in [("V1P",v1p),("RCTKO",rct),("V1TxS",txs)]:
    a=np.array([v for _,v in rows])
    print("%-6s  %%bins>null_p95=%4.1f%%   %%bins>null_p99=%4.1f%%   max=%.4f"%(
        tag,100*np.mean(a>p95),100*np.mean(a>p99),a.max()))
print("(chance is 5%/1%; V1TxS-vs-itself is the by-definition reference)")
print()
print("Top-10 V1P core bins by |dContact| (are they at the Trbv1 stripe ~40.9M / RC ~41.5M?):")
for start,v in sorted(v1p,key=lambda x:-x[1])[:10]:
    print("  chr6:%d  |dContact|=%.4f"%(start,v))
