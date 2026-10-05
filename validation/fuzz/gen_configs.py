"""Randomized comparison of benchsiDE differential expression with limma/edgeR: configuration generator.
Usage: python3 validation/fuzz/gen_configs.py [seed]   (needs tests/engine/data from fetch_data.py)
Writes validation/fuzz/work/<load>.matrix.tsv, <load>.design.tsv and configs.json."""
import json, os, sys, random
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); E = os.path.join(ROOT, "tests", "engine")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20261004
rng = random.Random(SEED); nrng = np.random.default_rng(SEED)
W = os.path.join(HERE, "work"); os.makedirs(W, exist_ok=True)

def eng_design(n):
    D = pd.read_csv(os.path.join(E, "designs", n + "_design.csv"), header=None, index_col=0, dtype=str, keep_default_na=False).T
    D.columns = ["sample"] + list(D.columns[1:]); return D.reset_index(drop=True)
def eng_matrix(n): return pd.read_csv(os.path.join(E, "data", n + ".tsv"), sep="\t", index_col=0)

LOADS = []   # (name, species, matrix df, design df, factor, real covariates, normsel, fmode)
demo = pd.read_csv(os.path.join(ROOT, "demo", "GSE63310_counts.tsv"), sep="\t", index_col=0)
ddemo = pd.read_csv(os.path.join(ROOT, "demo", "GSE63310_design.tsv"), sep="\t", dtype=str, keep_default_na=False)
ddemo["lane"] = ddemo["sample"].map(dict(GSM1545535="L004",GSM1545536="L004",GSM1545538="L004",GSM1545539="L006",GSM1545540="L006",GSM1545541="L006",GSM1545542="L006",GSM1545544="L008",GSM1545545="L008"))
LOADS.append(("GSE63310", "mouse", demo, ddemo[["sample","celltype","lane"]], "celltype", ["lane"], "tmm", "fbe"))
LOADS.append(("GSE63310_cpm", "mouse", demo, ddemo[["sample","celltype","lane"]], "celltype", ["lane"], "cpm", "cpm"))
LOADS.append(("GSE186063", "human", eng_matrix("GSE186063"), eng_design("GSE186063"), "Type", ["Diag","Sex","Age","Pair"], "tmm", "fbe"))
LOADS.append(("GSE121212", "human", eng_matrix("GSE121212"), eng_design("GSE121212"), "Group", ["Disease","Skin","Patient"], "tmm", "fbe"))
LOADS.append(("GSE143688_all", "mouse", eng_matrix("GSE143688_all"), eng_design("GSE143688_all"), "Cell", ["Treatment","Genotype","Day"], "tmm", "fbe"))
LOADS.append(("GSE83645", "human", eng_matrix("GSE83645"), eng_design("GSE83645"), "Cond", ["Patient","Site"], "tmm", "fbe"))
LOADS.append(("GSE41745", "human", eng_matrix("GSE41745"), eng_design("GSE41745"), "Cond", ["Patient"], "tmm", "fbe"))
d54 = eng_design("GSE54456"); pick = sorted(rng.sample(list(d54.index[d54.Cond=="normal_skin"]), 20) + rng.sample(list(d54.index[d54.Cond!="normal_skin"]), 20))
d54 = d54.loc[pick].reset_index(drop=True); m54 = eng_matrix("GSE54456")[d54["sample"].tolist()]
LOADS.append(("GSE54456_sub40", "human", m54, d54, "Cond", [], "tmm", "fbe"))
d171 = eng_design("GSE171012"); d171 = d171[d171.CellType=="bulk"].reset_index(drop=True); m171 = eng_matrix("GSE171012")[d171["sample"].tolist()]
LOADS.append(("GSE171012_bulk", "human", m171, d171[["sample","Status","Subject","Batch"]], "Status", ["Subject","Batch"], "tmm", "fbe"))
LOADS.append(("GSE54456_sub40_cpmfilter", "human", m54, d54, "Cond", [], "tmm", "cpm"))
# already-normalized and log-scale inputs (the counts converted outside the application)
cpm186 = eng_matrix("GSE186063"); cpm186 = cpm186 / cpm186.sum(axis=0) * 1e6
LOADS.append(("GSE186063_cpmvalues", "human", cpm186.round(4), eng_design("GSE186063"), "Type", ["Diag","Sex","Age","Pair"], "tmm", "cpm", "norm"))
LOADS.append(("GSE186063_log2p1", "human", np.log2(cpm186 + 1).round(6), eng_design("GSE186063"), "Type", ["Diag","Sex","Age","Pair"], "tmm", "cpm", "log"))
cpm143 = eng_matrix("GSE143688_all"); cpm143 = np.log2((cpm143 + 0.5) / (cpm143.sum(axis=0) + 1) * 1e6)
LOADS.append(("GSE143688_log2cpm", "mouse", cpm143.round(6), eng_design("GSE143688_all"), "Cell", ["Treatment","Genotype","Day"], "tmm", "cpm", "log"))
LOADS = [l if len(l) == 9 else l + ("counts",) for l in LOADS]

SYN = ["zCont", "zCat3", "zAlias", "zMiss", "zDonor", "zSingle", "zPartial"]
def add_synthetic(D, fac):
    n = len(D); g = D[fac].tolist()
    D = D.copy()
    D["zCont"] = [f"{v:.2f}" for v in nrng.normal(50, 10, n)]
    D["zCat3"] = [rng.choice(["x", "y", "z"]) for _ in range(n)]
    D["zAlias"] = ["k_" + v for v in g]                                   # identical to the groups
    D["zMiss"] = [("NA" if rng.random() < 0.15 else f"{v:.1f}") for v in nrng.normal(7, 1, n)]
    D["zDonor"] = [f"D{(i // 2) + 1:02d}" for i in rng.sample(range(n), n)]  # named as an identifier: categorical
    one = rng.randrange(n); D["zSingle"] = [("u" if i == one else rng.choice(["v", "w"])) for i in range(n)]
    levs = sorted(set(g)); first = levs[0]
    D["zPartial"] = [("b1" if v == first else rng.choice(["b1", "b2"])) for v in g]  # confounded with one group only
    return D

configs = []
for (name, sp, M, D, fac, real, norm, fmode, units) in LOADS:
    D = add_synthetic(D, fac)
    M = M.loc[:, D["sample"].tolist()]
    M.to_csv(os.path.join(W, name + ".matrix.tsv"), sep="\t")
    D.to_csv(os.path.join(W, name + ".design.tsv"), sep="\t", index=False)
    groups = D[fac].value_counts(); ok = [k for k, v in groups.items() if v >= 2 and k != "NA"]
    pool = real + SYN
    nconf = 30 if len(D) <= 80 else 22
    for i in range(nconf):
        A, B = rng.sample(ok, 2)
        m = rng.choices(["mod", "voom", "welch"], [0.42, 0.42, 0.16])[0]
        k = rng.choices([0, 1, 2], [0.3, 0.45, 0.25])[0]
        cov = rng.sample(pool, k)
        treat = m != "welch" and rng.random() < 0.35
        thr = rng.choice([0.5, 1, 1.5]) if treat else rng.choice([0, 0.58, 1, 2])
        nex = rng.choices([0, 1, 2, 3], [0.6, 0.15, 0.15, 0.1])[0]
        excl = rng.sample(D["sample"].tolist(), nex)
        configs.append(dict(id=f"{name}#{i}", load=name, A=A, B=B, m=m, fit=rng.choice(["all", "pair"]), cov=cov, treat=treat, thr=thr,
                            fdr=rng.choice([0.01, 0.05, 0.1]), excl=excl))
loads = [dict(name=n, species=sp, factor=fac, norm=norm, fmode=fmode, units=u) for (n, sp, M, D, fac, real, norm, fmode, u) in LOADS]
json.dump(dict(seed=SEED, loads=loads, configs=configs), open(os.path.join(W, "configs.json"), "w"), indent=0)
print(len(loads), "loads,", len(configs), "configurations")
