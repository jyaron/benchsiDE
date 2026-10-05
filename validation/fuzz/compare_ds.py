"""Compare comparison-dataset results (work/*.ds.jsonl) with the limma reference (work/*.ref.json)."""
import json, os, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(HERE, "work")
spec = json.load(open(os.path.join(W, "configs.json"))); cfg = {c["id"]: c for c in spec["configs"]}
n_fit = n_eq = 0; worst_t = worst_p = 0.0; findings = []
for L in spec["loads"]:
    fa, fr = os.path.join(W, L["name"] + ".ds.jsonl"), os.path.join(W, L["name"] + ".ref.json")
    if not (os.path.exists(fa) and os.path.exists(fr)): continue
    ref = json.load(open(fr)); load = None; rows = []
    for l in open(fa):
        if l.startswith("L "): load = json.loads(l[2:])
        elif l.startswith("C "): rows.append(json.loads(l[2:]))
    if not rows: continue
    if load["ids"] != ref["_load"]["ids"]: findings.append((L["name"], "filtered genes differ", f"{load['nG']} vs {ref['_load']['nG']}")); continue
    n = load["nG"]; step = max(1, n // 400); samp = list(range(0, n, step))
    for a in rows:
        r = ref[a["id"]]; c = cfg[a["id"]]
        if a.get("err") or r.get("err"):
            if not (a.get("err") and r.get("err")): findings.append((a["id"], "refusal differs", f"app {a.get('err')} | R {r.get('err')}"))
            continue
        n_fit += 1
        ca, cr = set(a["calls"]), set(r["calls"] if isinstance(r["calls"], list) else [r["calls"]])
        n_eq += ca == cr
        t = np.array(a["t"])[samp]; p = np.array(a["p"])[samp]
        dt = float(np.max(np.abs(t - np.array(r["t"], float)))); dp = float(np.max(np.abs(p - np.array(r["p"], float)) / np.maximum(np.array(r["p"], float), 1e-300)))
        worst_t, worst_p = max(worst_t, dt), max(worst_p, dp)
        if ca != cr or dt > 1e-6 or dp > 1e-6: findings.append((a["id"], "statistics differ", f"calls {len(ca)} vs {len(cr)}, t {dt:.1e}, p {dp:.1e}"))
        kept = r.get("kept") or []; kept = kept if isinstance(kept, list) else [kept]
        stated = [x.split(" (")[0].strip() for x in (a["cov"] or "").split(" + ") if x]
        if sorted(stated) != sorted(kept): findings.append((a["id"], "stated covariates differ from the fitted model", f"stated {stated}, fitted {kept}"))
        dropped = r.get("dropped") or []; dropped = dropped if isinstance(dropped, list) else [dropped]
        if sorted(a.get("notEst") or []) != sorted(dropped): findings.append((a["id"], "not-estimable covariates differ", f"app {a.get('notEst')}, R {dropped}"))
        m = re.search(r"adjusted within each dataset for (.*?)\); raw expression", a.get("meth", ""))
        seg = m.group(1) if m else ""
        for cv in dropped:
            if re.search(r"(:|\+) " + re.escape(cv) + r" \(", seg): findings.append((a["id"], "methods state an unfitted covariate as adjusted", seg[:150]))
            if cv not in seg: findings.append((a["id"], "methods omit the not-estimable covariate", seg[:150]))
        mi = a.get("mi") or {}
        if mi.get("subj") and mi["subj"]["col"] not in kept: findings.append((a["id"], "subjects counted for a covariate not in the model", str(mi["subj"])))
print("comparison datasets:", n_fit, "fitted,", n_eq, "with identical calls; max |dt|", f"{worst_t:.1e}", "max rel dp", f"{worst_p:.1e}")
for f in findings: print(" ", f)
json.dump(findings, open(os.path.join(W, "ds_findings.json"), "w"), indent=1)
