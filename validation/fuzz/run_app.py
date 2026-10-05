"""Randomized comparison: run every configuration in the application (JavaScriptCore engine).
Usage: python3 validation/fuzz/run_app.py <load name>   -> work/<load>.app.jsonl"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
exec(open(os.path.join(ROOT, "tests", "engine", "harness.py")).read())
APP = app_of(open(os.path.join(ROOT, "index.html"), encoding="utf-8").read())
W = os.path.join(HERE, "work"); spec = json.load(open(os.path.join(W, "configs.json")))
name = sys.argv[1]; L = next(l for l in spec["loads"] if l["name"] == name); C = [c for c in spec["configs"] if c["load"] == name]
M = open(os.path.join(W, name + ".matrix.tsv")).read(); D = open(os.path.join(W, name + ".design.tsv")).read()
js = r"""
  el("species").value=""" + json.dumps(L["species"]) + r"""; useBuiltinAnno(); el("normsel").value=""" + json.dumps(L["norm"]) + r"""; el("fmode").value=""" + json.dumps(L["fmode"]) + r"""; el("unitsel").value=""" + json.dumps(L.get("units","counts")) + r""";
  loadFromText(""" + json.dumps(M) + ", " + json.dumps(D) + r"""); applyFactor(""" + json.dumps(L["factor"]) + r"""); GROUPTYPE="categorical"; analyzeNow(); __drain();
  switchTab("de"); __drain();
  const inc0=included.slice();
  const step=Math.max(1, Math.floor(nG/400)), SAMP=[]; for(let g=0;g<nG;g+=step) SAMP.push(g);
  print("L " + JSON.stringify({nG, ids:IDS, samples:S.map(x=>x.id), inc0, samp:SAMP, normNote:NORMNOTE}));
  const CF=""" + json.dumps(C) + r""";
  for(const c of CF){
    for(let q=0;q<nS;q++) included[q]=inc0[q];
    for(const sn of c.excl){ const q=S.findIndex(x=>x.id===sn); if(q>=0) included[q]=false; }
    buildSampleBar && buildSampleBar();
    el("selA").value=c.A; el("selB").value=c.B; el("selFDR").value=String(c.fdr); el("fcthr").value=String(c.thr);
    deMethod=c.m; DEFIT=c.fit; DECOVSEL=c.cov.slice(); DETREAT=c.treat; deKey=""; syncTreatUI();
    let out={id:c.id};
    try{
      drawDE(); const r=computeDE();
      if(r.err){ out.err=String(r.err); out.summary=el("deprior").textContent; print("C " + JSON.stringify(out)); continue; }
      const fdr=c.fdr, thr=c.thr, calls=[];
      for(let g=0;g<nG;g++){ const ok = r.treat? r.q[g]<=fdr : (r.q[g]<=fdr && Math.abs(r.lfc[g])>=thr); if(ok) calls.push(r.lfc[g]>0? g : -g-1); }
      let csv=null; const _d=downloadCSV; downloadCSV=(t,f)=>{csv=t;}; el("deexport")._h["click"](); downloadCSV=_d;
      const pr=r.prior||{};
      Object.assign(out, {treat:!!r.treat, calls, lfc:SAMP.map(g=>r.lfc[g]), t:SAMP.map(g=>r.tv[g]), p:SAMP.map(g=>r.p[g]),
        prior:{dfr:pr.dfr, df0:pr.df0, nFit:pr.nFit, nA:pr.nA, nB:pr.nB, cov:pr.cov||"", dropNA:pr.dropNA||0, groupsFit:pr.groupsFit||[], dropped:pr.covDropped||[], terms:pr.covTerms||[]},
        summary:el("deprior").textContent, counts:el("decounts").textContent, methods:deMethodsText(), genMethods:(genMethods(), el("methodstext").textContent), csvhead:csv.split("\n").slice(0,3), rule:deRuleText(fdr,thr),
        csvCalls: csv.split("\n").filter(l=>l&&!l.startsWith("#")).slice(1).filter(l=>/,(up|down)$/.test(l)).length});
    }catch(e){ out.exc=String(e); }
    print("C " + JSON.stringify(out));
  }
"""
o = run(APP, js, timeout=20000)
with open(os.path.join(W, name + ".app.jsonl"), "w") as f:
    for l in o.split("\n"):
        if l.startswith(("L ", "C ")): f.write(l + "\n")
errs = [l for l in o.split("\n") if l.startswith(("Exception", "RUNTIME", "SyntaxError"))]
print(name, sum(1 for l in o.split("\n") if l.startswith("C ")), "configs", errs[:3])
