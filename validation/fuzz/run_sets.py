"""Randomized comparison: FRY and over-representation (ORA) for the fitted moderated-t and voom configurations.
Usage: python3 validation/fuzz/run_sets.py <load> -> work/<load>.sets.jsonl (Hallmark sets of the load's species)"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
exec(open(os.path.join(ROOT, "tests", "engine", "harness.py")).read())
APP = app_of(open(os.path.join(ROOT, "index.html"), encoding="utf-8").read())
W = os.path.join(HERE, "work"); spec = json.load(open(os.path.join(W, "configs.json")))
name = sys.argv[1]; L = next(l for l in spec["loads"] if l["name"] == name)
C = [c for c in spec["configs"] if c["load"] == name and c["m"] in ("mod", "voom")]
M = open(os.path.join(W, name + ".matrix.tsv")).read(); D = open(os.path.join(W, name + ".design.tsv")).read()
lib = "hallmark_hs" if L["species"] == "human" else "hallmark_mm"
js = r"""
  el("species").value=""" + json.dumps(L["species"]) + r"""; useBuiltinAnno(); el("normsel").value=""" + json.dumps(L["norm"]) + r"""; el("fmode").value=""" + json.dumps(L["fmode"]) + r"""; el("unitsel").value=""" + json.dumps(L.get("units","counts")) + r""";
  loadFromText(""" + json.dumps(M) + ", " + json.dumps(D) + r"""); applyFactor(""" + json.dumps(L["factor"]) + r"""); GROUPTYPE="categorical"; analyzeNow(); __drain();
  useBuiltin(""" + json.dumps(lib) + r"""); __drain();
  el("enrmin").value="5"; el("enrmax").value="2000";
  const SETS=[]; for(const st of GMT){ const ix=new Set(); for(const gn of st.genes){ const i=nameToIdx[gn.toLowerCase()]??idToIdx[gn.toLowerCase()]; if(i!==undefined) ix.add(i); } if(ix.size>=5&&ix.size<=2000) SETS.push({name:st.name, rows:[...ix].sort((a,b)=>a-b)}); }
  const U=[...oraUniverse()].sort((a,b)=>a-b);
  print("S " + JSON.stringify({sets:SETS, universe:U, nG, groupOrder:GROUPS}));
  const inc0=included.slice(); const CF=""" + json.dumps(C) + r""";
  for(const c of CF){
    for(let q=0;q<nS;q++) included[q]=inc0[q];
    for(const sn of c.excl){ const q=S.findIndex(x=>x.id===sn); if(q>=0) included[q]=false; }
    el("selA").value=c.A; el("selB").value=c.B; el("selFDR").value=String(c.fdr); el("fcthr").value=String(c.thr);
    deMethod=c.m; DEFIT=c.fit; DECOVSEL=c.cov.slice(); DETREAT=c.treat; deKey=""; syncTreatUI();
    const r=computeDE(); const out={id:c.id};
    if(r.err){ out.err=r.err; print("C " + JSON.stringify(out)); continue; }
    const fr=fryTest(SETS);
    if(fr.err) out.fryErr=fr.err; else { out.fry=fr.rows.map(x=>[x.name, x.m, x.p, x.pm, x.dir]); out.fryCov=fr.covUsed; out.fryOut=fr.covOut; out.fryN=fr.nFit; }
    el("enrsrc").value="de_up"; runEnrichment();
    out.ora = enrResults&&enrResults.rows? enrResults.rows.map(x=>[x.name, x.m, x.k, x.p]) : null;
    out.oraNote = el("enrinfo")? el("enrinfo").textContent.slice(0,300) : "";
    print("C " + JSON.stringify(out));
  }
"""
o = run(APP, js, timeout=20000)
with open(os.path.join(W, name + ".sets.jsonl"), "w") as f:
    for l in o.split("\n"):
        if l.startswith(("S ", "C ")): f.write(l + "\n")
print(name, len(C), "configs", [l[:300] for l in o.split("\n") if l.startswith(("Exception", "RUNTIME"))][:3])
