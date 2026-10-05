"""Randomized checks of the Gene Explorer tests (against scipy/statsmodels) and of the per-gene meta-analysis
(against metafor). Usage: python3 validation/fuzz/stat_tests.py   (RSCRIPT=path to Rscript with metafor)
Writes work/gene_tests.json, work/meta_fuzz.json, work/meta_fuzz_ref.json and prints a summary."""
import json, os, sys, math, subprocess
import numpy as np
from scipy import stats as SS
from statsmodels.stats.multitest import multipletests
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); W = os.path.join(HERE, "work")
exec(open(os.path.join(ROOT, "tests", "engine", "harness.py")).read())
APP = app_of(open(os.path.join(ROOT, "index.html"), encoding="utf-8").read())
def lw(o, pre): return [json.loads(l[len(pre):]) for l in o.split("\n") if l.startswith(pre)]
def rel(a, b): return abs(a - b) / max(abs(b), 1e-300)

# ---- Gene Explorer: one-way ANOVA, Kruskal-Wallis, Spearman trend, pairwise Welch with Holm ----
rng = np.random.default_rng(5); CASES = []
for i in range(400):
    k = int(rng.integers(2, 6)); gs = []
    for j in range(k):
        n = int(rng.integers(2, 9)); u = rng.random(); v = rng.normal(rng.normal(0, 1), rng.uniform(0.1, 2), n)
        if u < 0.25: v = np.round(v, 0)
        if u > 0.95: v = np.full(n, 3.0)
        gs.append([float(x) for x in v])
    CASES.append(gs)
o = run(APP, r"""
  const C=""" + json.dumps(CASES) + r""";
  print("G " + JSON.stringify(C.map(gs=>{ const a=anova(gs), k=kruskal(gs); const xs=[], ys=[]; gs.forEach((g,i)=>g.forEach(v=>{ xs.push(i); ys.push(v); }));
    const rho=spearman(xs,ys); const tt=Math.abs(rho)<1? rho*Math.sqrt((xs.length-2)/(1-rho*rho)) : null;
    const pw=[]; for(let i=0;i<gs.length;i++) for(let j=i+1;j<gs.length;j++){ const w=welch(gs[i],gs[j]); pw.push(w&&w.se>0? w.p : null); }
    const ti=pw.map((p,i)=>i).filter(i=>pw[i]!==null), h=holm(ti.map(i=>pw[i]));
    return {pA:a&&a.p, pK:k&&k.p, rho, tp:tt===null? 0 : tPval(tt, xs.length-2), pw, holm:h}; })));
""", timeout=600)
G = lw(o, "G ")[0]; bad = {"anova": 0, "kruskal": 0, "spearman": 0, "welch": 0, "holm": 0}; worst = dict.fromkeys(bad, 0.0); n_ut = 0
for gs, r in zip(CASES, G):
    arr = [np.array(g) for g in gs]; x = np.concatenate([np.full(len(g), i) for i, g in enumerate(gs)]); y = np.concatenate(arr)
    with np.errstate(all="ignore"):
        refs = {"anova": SS.f_oneway(*arr).pvalue, "spearman": SS.spearmanr(x, y).pvalue}
        try: refs["kruskal"] = SS.kruskal(*arr).pvalue
        except ValueError: refs["kruskal"] = float("nan")
    for key, av in (("anova", r["pA"]), ("kruskal", r["pK"]), ("spearman", r["tp"])):
        rv = refs[key]
        if rv is None or np.isnan(rv): continue
        if av is None: bad[key] += 1; continue
        d = rel(av, rv); worst[key] = max(worst[key], d); bad[key] += int(d > 1e-6)
    tested = []
    for (i, j), av in zip([(i, j) for i in range(len(gs)) for j in range(i + 1, len(gs))], r["pw"]):
        constant = np.ptp(arr[i]) == 0 and np.ptp(arr[j]) == 0
        if constant: n_ut += 1; bad["welch"] += int(av is not None); continue
        with np.errstate(all="ignore"): rv = SS.ttest_ind(arr[j], arr[i], equal_var=False).pvalue
        d = rel(av, rv); worst["welch"] = max(worst["welch"], d); bad["welch"] += int(d > 1e-6); tested.append(rv)
    if tested:
        ho = multipletests(tested, method="holm")[1]; d = max(abs(a - b) for a, b in zip(ho, r["holm"])); worst["holm"] = max(worst["holm"], d); bad["holm"] += int(d > 1e-9)
json.dump(dict(cases=len(CASES), bad=bad, worst={k: float(v) for k, v in worst.items()}, untestable_pairs=n_ut), open(os.path.join(W, "gene_tests.json"), "w"))
print("Gene Explorer tests:", len(CASES), "random designs; disagreements", bad, "; worst", {k: f"{v:.1e}" for k, v in worst.items()}, "; pairs without variance (not testable):", n_ut)

# ---- meta-analysis: 3000 random cases x 6 model/test variants against metafor ----
rng = np.random.default_rng(99); MC = []
for i in range(3000):
    k = int(rng.choice([2, 2, 3, 3, 4, 5, 6, 8, 12])); u = rng.random()
    v = np.exp(rng.normal(-2, 1.5 if u < 0.3 else 0.6, k))
    if u < 0.1: y = np.full(k, rng.normal())
    elif u < 0.25: y = rng.normal(1, 0.01, k)
    elif u < 0.4: y = rng.normal(0, 3, k)
    else: y = rng.normal(rng.normal(0, 1), np.sqrt(v + rng.exponential(0.2)))
    if u > 0.95: v[0] *= 1e4
    MC.append(dict(y=[float(x) for x in y], v=[float(x) for x in v]))
methods = [["FEM", False], ["REM", False], ["REML", False], ["REM", "hk"], ["REML", "hk"], ["REML", "adhoc"]]
o = run(APP, r"""const C=""" + json.dumps(MC) + r"""; const M=""" + json.dumps(methods) + r""";
  print("M " + JSON.stringify(C.map(c=>M.map(([m,hk])=>{ const r=metaGene(c.y, c.v, m, hk); return [r.mu, r.seTest, r.p, r.tau2, r.I2, r.ciLo, r.ciHi, r.piLo, r.piHi, r.Q]; }))));""", timeout=600)
json.dump(dict(cases=MC, methods=methods, app=lw(o, "M ")[0]), open(os.path.join(W, "meta_fuzz.json"), "w"))
subprocess.run([os.environ.get("RSCRIPT", "Rscript"), os.path.join(HERE, "meta_reference.R")], check=True)
J = json.load(open(os.path.join(W, "meta_fuzz.json"))); REF = json.load(open(os.path.join(W, "meta_fuzz_ref.json")))
fields = ["mu", "se", "p", "tau2", "I2", "lo", "hi", "plo", "phi", "Q"]; names = ["FE z", "DL z", "REML z", "DL HK", "REML HK", "REML HK truncated"]
ident = {i for i, c in enumerate(J["cases"]) if max(c["y"]) - min(c["y"]) < 1e-12}
out = {}
for mi, nm in enumerate(names):
    n = b = 0; wv = 0.0
    for ci, (a_, r_) in enumerate(zip(J["app"], REF)):
        a, r = a_[mi], r_[mi]
        if r.get("err") or ci in ident: continue
        n += 1
        for fi, f in enumerate(fields):
            rv = r[f]
            if rv is None: continue
            d = abs(a[fi] - rv) if f in ("mu", "tau2", "I2", "lo", "hi", "plo", "phi") else rel(a[fi], rv)
            wv = max(wv, d); b += d > (1e-6 if f != "tau2" else 1e-7 * max(1, abs(rv)))
    out[nm] = dict(cases=n, disagreements=b, worst=wv)
json.dump(dict(variants=out, identical_estimate_cases=len(ident)), open(os.path.join(W, "meta_summary.json"), "w"), indent=1)
print("meta-analysis vs metafor:", {k: (v["cases"], v["disagreements"], f"{v['worst']:.0e}") for k, v in out.items()}, "; cases with identical estimates (excluded; Q and the Hartung-Knapp SE are rounding-level, not zero):", len(ident))
