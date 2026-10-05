"""Compare FRY and ORA results (work/*.sets.jsonl) with the reference (work/*.setsref.json)."""
import json, os, math
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(HERE, "work")
spec = json.load(open(os.path.join(W, "configs.json")))
tot = dict(fry=0, ora=0); worst = dict(fry_p=0.0, fry_pm=0.0, ora_p=0.0); findings = []
def rel(a, b): return abs(a - b) / max(abs(b), 1e-300)
for L in spec["loads"]:
    fa, fr = os.path.join(W, L["name"] + ".sets.jsonl"), os.path.join(W, L["name"] + ".setsref.json")
    if not (os.path.exists(fa) and os.path.exists(fr)): continue
    ref = json.load(open(fr))
    for l in open(fa):
        if not l.startswith("C "): continue
        a = json.loads(l[2:]); r = ref.get(a["id"], {})
        if a.get("err") or r.get("err"):
            if bool(a.get("err")) != bool(r.get("err")): findings.append((a["id"], "DE refusal differs", f"{a.get('err')} | {r.get('err')}"))
            continue
        if a.get("fryErr"): findings.append((a["id"], "FRY refused, limma fitted", a["fryErr"])); continue
        rk = r.get("kept") or []; rk = rk if isinstance(rk, list) else [rk]
        if sorted(a.get("fryCov") or []) != sorted(rk): findings.append((a["id"], "FRY covariates differ", f"app {a.get('fryCov')}, R {rk}"))
        if a.get("fryN") != r.get("n"): findings.append((a["id"], "FRY samples differ", f"app {a.get('fryN')}, R {r.get('n')}"))
        R_ = {x[0]: x for x in r["fry"]}
        for nm, m, p, pm, d in a["fry"]:
            x = R_.get(nm); tot["fry"] += 1
            if x is None: findings.append((a["id"], "FRY set missing in R", nm)); continue
            dp, dpm = rel(p, x[2]), rel(pm, x[3]); worst["fry_p"] = max(worst["fry_p"], dp); worst["fry_pm"] = max(worst["fry_pm"], dpm)
            if m != x[1] or d != x[4] or dp > 1e-6 or dpm > 1e-6: findings.append((a["id"], "FRY differs", f"{nm}: n {m}/{x[1]} dir {d}/{x[4]} p {p:.4g}/{x[2]:.4g} pm {pm:.4g}/{x[3]:.4g}"))
        ro = r.get("ora"); ao = a.get("ora")
        if (ro is None) != (ao is None): findings.append((a["id"], "ORA run differs", f"app {'ran' if ao else 'not run'}, R {'ran' if ro else 'not run'}: {a.get('oraNote','')[:120]}")); continue
        if ro is None: continue
        RO = {x[0]: x for x in ro}
        for nm, m, k, p in ao:
            x = RO.get(nm); tot["ora"] += 1
            if x is None: findings.append((a["id"], "ORA set missing in R", nm)); continue
            dp = rel(p, x[3]); worst["ora_p"] = max(worst["ora_p"], dp)
            if m != x[1] or k != x[2] or dp > 1e-6: findings.append((a["id"], "ORA differs", f"{nm}: m {m}/{x[1]} k {k}/{x[2]} p {p:.4g}/{x[3]:.4g}"))
# FRY: separate configurations that agree with limma exactly from those whose design makes FRY depend on the
# order of the samples (limma itself changes when the samples are reordered: r[5], r[6])
from statsmodels.stats.multitest import multipletests
fry_rows = []
for L in spec["loads"]:
    fa, fr = os.path.join(W, L["name"] + ".sets.jsonl"), os.path.join(W, L["name"] + ".setsref.json")
    if not (os.path.exists(fa) and os.path.exists(fr)): continue
    ref = json.load(open(fr))
    for l in open(fa):
        if not l.startswith("C "): continue
        a = json.loads(l[2:]); r = ref.get(a["id"], {})
        if a.get("err") or r.get("err") or a.get("fryErr") or not r.get("fry"): continue
        R_ = {x[0]: x for x in r["fry"]}; nm = [x[0] for x in a["fry"]]
        ap = [x[2] for x in a["fry"]]; rp = [R_[n][2] for n in nm]; r2 = [R_[n][5] for n in nm]
        dA = max(rel(x, y) for x, y in zip(ap, rp)); dL = max(rel(x, y) for x, y in zip(r2, rp))
        sa, sr, s2 = [multipletests(v, method="fdr_bh")[1] <= 0.05 for v in (ap, rp, r2)]
        fry_rows.append((a["id"], dA, dL, int((sa != sr).sum()), int((s2 != sr).sum()), int(sr.sum())))
ex = [x for x in fry_rows if x[1] <= 1e-6]; df_ = [x for x in fry_rows if x[1] > 1e-6]
print(f"FRY: {len(fry_rows)} configurations; {len(ex)} equal to limma (rel. diff <= 1e-6); {len(df_)} differ, of which {sum(1 for x in df_ if x[2] > 1e-6)} are designs where limma itself changes when the samples are reordered")
if df_:
    print(f"  in those: app vs limma max rel. diff {max(x[1] for x in df_):.2f}, limma vs reordered limma {max(x[2] for x in df_):.2f}; sets called differently at FDR 0.05: app {sum(x[3] for x in df_)}, reordered limma {sum(x[4] for x in df_)} (of {sum(x[5] for x in df_)} limma calls)")
json.dump(fry_rows, open(os.path.join(W, "fry_order.json"), "w"))
findings = [f for f in findings if f[1] != "FRY differs" or f[0] not in {x[0] for x in df_ if x[2] > 1e-6}]
print("set tests compared:", tot, {k: f"{v:.1e}" for k, v in worst.items()})
from collections import Counter
print(Counter(f[1] for f in findings))
for f in findings[:12]: print(" ", f)
json.dump(findings, open(os.path.join(W, "sets_findings.json"), "w"), indent=1)
