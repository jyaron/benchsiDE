"""Randomized comparison: compare application results and texts with the limma/edgeR reference.
Usage: python3 validation/fuzz/compare.py   -> work/fuzz_results.csv, work/fuzz_findings.json"""
import json, os, re, math
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(HERE, "work")
spec = json.load(open(os.path.join(W, "configs.json"))); cfg = {c["id"]: c for c in spec["configs"]}
rows, findings = [], []
def find(cid, kind, detail): findings.append(dict(id=cid, kind=kind, detail=detail))
def num(s): return float(s.replace(",", ""))
for L in spec["loads"]:
    n = L["name"]; fa, fr = os.path.join(W, n + ".app.jsonl"), os.path.join(W, n + ".ref.json")
    if not (os.path.exists(fa) and os.path.exists(fr)): find(n, "missing output", ""); continue
    app = {}; load = None
    for l in open(fa):
        if l.startswith("L "): load = json.loads(l[2:])
        elif l.startswith("C "): o = json.loads(l[2:]); app[o["id"]] = o
    ref = json.load(open(fr)); rl = ref["_load"]
    aslist = lambda v: [] if v is None else (v if isinstance(v, list) else [v])
    for k_, v_ in ref.items():
        if k_ != "_load":
            for f_ in ("kept", "dropped", "groups", "calls", "lfc", "t", "p"): v_[f_] = aslist(v_.get(f_))
    if load["nG"] != rl["nG"] or load["ids"] != rl["ids"]: find(n, "filtered genes differ", f"app {load['nG']}, R {rl['nG']}")
    for cid, c in cfg.items():
        if c["load"] != n: continue
        a, r = app.get(cid), ref.get(cid)
        row = dict(id=cid, method=c["m"], fit=c["fit"], cov="+".join(c["cov"]), treat=c["treat"], thr=c["thr"], fdr=c["fdr"], n_excluded=len(c["excl"]))
        if a is None: find(cid, "no app output", ""); continue
        if "exc" in a: find(cid, "app exception", a["exc"]); row["status"] = "app exception"; rows.append(row); continue
        ae, re_ = a.get("err"), r.get("err")
        if ae or re_:
            row["status"] = "both refused" if (ae and re_) else ("app refused, R fitted" if ae else "app fitted, R refused")
            row["app_err"], row["ref_err"] = ae, re_
            if not (ae and re_): find(cid, row["status"], f"app: {ae} | R: {re_}")
            rows.append(row); continue
        row["status"] = "both fitted"
        ca, cr = set(a["calls"]), set(r["calls"])
        row.update(calls_app=len(ca), calls_R=len(cr), calls_equal=ca == cr)
        la, lr = np.array(a["lfc"], float), np.array(r["lfc"], float); ta, tr = np.array(a["t"], float), np.array(r["t"], float)
        pa, pr_ = np.array(a["p"], float), np.array(r["p"], float)
        row.update(d_lfc=float(np.nanmax(np.abs(la - lr))), d_t=float(np.nanmax(np.abs(ta - tr))), d_p_rel=float(np.nanmax(np.abs(pa - pr_) / np.maximum(pr_, 1e-300))))
        if ca != cr: find(cid, "significant genes differ", f"app {len(ca)}, R {len(cr)}, symmetric difference {len(ca ^ cr)}")
        if row["d_t"] > 1e-6 or row["d_p_rel"] > 1e-6 or row["d_lfc"] > 1e-6: find(cid, "statistics differ", f"lfc {row['d_lfc']:.2e} t {row['d_t']:.2e} p {row['d_p_rel']:.2e}")
        P = a["prior"]; S, MT, H = a["summary"], a["methods"], " ".join(a["csvhead"])
        if c["m"] != "welch":
            if P["dfr"] != r["dfr"]: find(cid, "residual df differs", f"app {P['dfr']}, R {r['dfr']}")
            d0a, d0r = P["df0"], r["df0"]   # an infinite prior df is serialized as null by both
            if (d0a is None) != (d0r is None) or (d0a is not None and abs(d0a - d0r) > 1e-6 * max(1, abs(d0r))): find(cid, "prior df differs", f"app {d0a}, R {d0r}")
            if P["nFit"] != r["n"]: find(cid, "fitted samples differ", f"app {P['nFit']}, R {r['n']}")
            # covariate typing and the text describing the model
            kinds = dict(re.findall(r"([^+()]+?) \((continuous|categorical|2-level)\)", P["cov"] or ""))
            kinds = {k.strip(): ("categorical" if v != "continuous" else v) for k, v in kinds.items()}
            for cv, ty in (r.get("types") or {}).items():
                if kinds.get(cv) != ty: find(cid, "covariate type differs", f"{cv}: app {kinds.get(cv)}, rule {ty}")
            kept, dropped = r.get("kept") or [], r.get("dropped") or []
            m = re.search(r"Adjusted for (.+?) \((\d+) covariate terms? in the model; residual df = (\d+)\)", S)
            stated = [x.split(" (")[0].strip() for x in m.group(1).split(" + ")] if m else []
            if sorted(stated) != sorted(kept): find(cid, "summary: covariates stated", f"stated {stated}, fitted {kept}")
            for cv in dropped:
                if not re.search(re.escape(cv) + r"[^.]*could not be estimated", S): find(cid, "summary: dropped covariate not stated", cv)
                if not re.search(re.escape(cv) + r"[^.]*could not be estimated", MT): find(cid, "methods: dropped covariate not stated", cv)
            mh = re.search(r"covariates: ([^;]+);", H)
            hs = [x.split(" (")[0].strip() for x in mh.group(1).split(" + ")] if mh else []
            if sorted(hs) != sorted(kept): find(cid, "CSV header: covariates stated", f"stated {hs}, fitted {kept}")
            mm = re.search(r"(?:plus|with) ([^.(]+?(?:\([a-z0-9-]+\)[^.(]*?)*?)(?: as additional terms|\),| \()", MT)
            for cv in kept:
                if cv not in MT: find(cid, "methods: fitted covariate not named", cv)
            for cv in dropped:
                if re.search(r"(plus|with) [^.]*\b" + re.escape(cv) + r" \((continuous|categorical|2-level)\)", MT): find(cid, "methods: dropped covariate listed as fitted", cv)
            if len(r["groups"]) > 2:
                if f"all {r['n']} samples of the {len(r['groups'])} groups" not in MT: find(cid, "methods: sample statement", f"expected all {r['n']} samples of {len(r['groups'])} groups")
            else:
                if f"{r['nA']} {c['A']} and {r['nB']} {c['B']} samples" not in MT: find(cid, "methods: sample statement", f"expected {r['nA']} {c['A']} and {r['nB']} {c['B']}")
            if f"residual df = {r['dfr']}" not in S and kept: find(cid, "summary: residual df", f"expected {r['dfr']}")
        else:
            if c["cov"] and "UNADJUSTED" not in S: find(cid, "Welch with covariates: no warning", S[:120])
            if f"{r['nA']} {c['A']} and {r['nB']} {c['B']} samples" not in MT: find(cid, "methods: sample statement", f"expected {r['nA']} {c['A']} and {r['nB']} {c['B']}")
            if "fitted " in H.split("call:")[0] and "coefficient" in H: find(cid, "CSV header: Welch log2FC described as a fitted coefficient", "")
        # call rule in every text
        if c["treat"]:
            if f"TREAT test of |log2FC| > {c['thr']:g}" not in H.replace(".0 ", " ") and f"TREAT test of |log2FC| > {c['thr']}" not in H: find(cid, "CSV header: call rule", H[-200:])
            if "TREAT" not in MT: find(cid, "methods: TREAT not stated", "")
        else:
            exp = f"FDR <= {c['fdr']:g}" + (f" and |log2FC| >= {c['thr']:g}" if c["thr"] > 0 else "")
            if exp not in H: find(cid, "CSV header: call rule", f"expected '{exp}' in {H[-160:]}")
        if a["csvCalls"] != len(cr): find(cid, "CSV calls differ", f"CSV {a['csvCalls']}, R {len(cr)}")
        mc = re.findall(r"([\d,]+) up in .*?([\d,]+) down", a["counts"])
        if mc:
            up, dn = num(mc[0][0]), num(mc[0][1]); rup = sum(1 for g in cr if g >= 0); rdn = len(cr) - rup
            if (up, dn) != (rup, rdn): find(cid, "counts shown differ", f"shown {up:g}/{dn:g}, R {rup}/{rdn}")
        else: find(cid, "counts not parsed", a["counts"][:100])
        GM = a.get("genMethods", "")
        if c["excl"]:
            nex = sum(1 for x in c["excl"])
            if not re.search(r"Gene filtering and normalization used all \d+ loaded samples; the " + str(nex) + r" excluded sample", GM): find(cid, "methods: excluded samples not stated", f"{nex} excluded")
            if any(x not in GM for x in c["excl"]): find(cid, "methods: excluded sample names missing", "")
        rows.append(row)
R = pd.DataFrame(rows); R.to_csv(os.path.join(W, "fuzz_results.csv"), index=False)
json.dump(findings, open(os.path.join(W, "fuzz_findings.json"), "w"), indent=1)
print(len(R), "configurations;", R.status.value_counts().to_dict())
if "calls_equal" in R: print("fitted:", int((R.status == "both fitted").sum()), "calls equal:", int((R.calls_equal == True).sum()), "max d_t", R.d_t.max(), "max d_p_rel", R.d_p_rel.max())
F = pd.DataFrame(findings); print(F.kind.value_counts().to_string() if len(F) else "no findings")
