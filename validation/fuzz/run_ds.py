"""Randomized comparison: comparison-dataset differential expression (modTwoGroupB) for the moderated-t
configurations without TREAT or sample exclusion. Usage: python3 validation/fuzz/run_ds.py <load> -> work/<load>.ds.jsonl"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
exec(open(os.path.join(ROOT, "tests", "engine", "harness.py")).read())
APP = app_of(open(os.path.join(ROOT, "index.html"), encoding="utf-8").read())
W = os.path.join(HERE, "work"); spec = json.load(open(os.path.join(W, "configs.json")))
name = sys.argv[1]; L = next(l for l in spec["loads"] if l["name"] == name)
C = [c for c in spec["configs"] if c["load"] == name and c["m"] == "mod" and not c["treat"] and not c["excl"]]
if L.get("units", "counts") != "counts" or L["norm"] != "tmm" or L["fmode"] != "fbe": C = []
M = open(os.path.join(W, name + ".matrix.tsv")).read(); D = open(os.path.join(W, name + ".design.tsv")).read()
sess = open(os.path.join(ROOT, "demo", "GSE63310_counts.tsv")).read(); sdes = open(os.path.join(ROOT, "demo", "GSE63310_design.tsv")).read()
js = r"""
  el("species").value=""" + json.dumps(L["species"]) + r"""; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe"; el("unitsel").value="counts";
  loadFromText(""" + json.dumps(sess) + ", " + json.dumps(sdes) + r"""); GROUPTYPE="categorical"; analyzeNow(); __drain();
  el("cmpspecies").value="same"; CMPB_TXT=""" + json.dumps(M) + r"""; CMPB_DES=""" + json.dumps(D) + r"""; CMPB_NAME="fz"; onAddDataset(); __drain();
  const d=DSETS[DSETS.length-1]; if(d.gfactor!==""" + json.dumps(L["factor"]) + r""") setDsFactor(d, """ + json.dumps(L["factor"]) + r""");
  print("L " + JSON.stringify({nG:d.nG, ids:d.ids, gfactor:d.gfactor}));
  const CF=""" + json.dumps(C) + r""";
  for(const c of CF){
    d.selA=c.A; d.selB=c.B; d.covSel=c.cov.slice(); DEFIT=c.fit; d.resKey="";
    const r=dsResult(DSETS.length); const out={id:c.id};
    if(r.err){ out.err=r.err; print("C " + JSON.stringify(out)); continue; }
    const calls=[]; for(let g=0;g<d.nG;g++) if(r.q[g]<=c.fdr && Math.abs(r.lfc[g])>=c.thr) calls.push(r.lfc[g]>0? g : -g-1);
    Object.assign(out, {calls, lfc:Array.from(r.lfc), t:Array.from(r.tv), p:Array.from(r.p), cov:d.res.cov, dfr:d.res.dfr, df0:d.res.df0, nFit:d.res.nFit, notEst:d.res.covNotEstimable||[], meth:(genMethods(), el("methodstext").textContent), mi:metaSampleInfo(DSETS.length)});
    print("C " + JSON.stringify(out));
  }
"""
o = run(APP, js, timeout=20000)
with open(os.path.join(W, name + ".ds.jsonl"), "w") as f:
    for l in o.split("\n"):
        if l.startswith(("L ", "C ")): f.write(l + "\n")
print(name, len(C), "configs", [l for l in o.split("\n") if l.startswith(("Exception", "RUNTIME"))][:3])
