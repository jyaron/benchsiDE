"""Snapshot of analysis results on public data, for comparing builds.
Usage: python3 tests/engine/snapshot.py OUT.json   (requires tests/engine/data from fetch_data.py)"""
import json, os, sys, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
exec(open(os.path.join(HERE, "harness.py")).read())
APP = app_of(open(os.path.join(ROOT, "index.html"), encoding="utf-8").read())
def mat(a): return open(os.path.join(HERE, "data", a + ".tsv")).read()
def des(a): return open(os.path.join(HERE, "designs", a + "_design.csv")).read()
J = json.dumps
SETUP = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE54456")) + ", " + J(des("GSE54456")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="normal_skin"; el("selB").value="Psoriasis_skin"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DECOVSEL=[]; deKey="";
  el("cmpspecies").value="human";
  CMPB_TXT=""" + J(mat("GSE121212")) + r"""; CMPB_DES=""" + J(des("GSE121212")) + r"""; CMPB_NAME="GSE121212"; onAddDataset(); __drain();
  CMPB_TXT=""" + J(mat("GSE186063")) + r"""; CMPB_DES=""" + J(des("GSE186063")) + r"""; CMPB_NAME="GSE186063"; onAddDataset(); __drain();
  CMPB_TXT=""" + J(mat("GSE83645")) + r"""; CMPB_DES=""" + J(des("GSE83645")) + r"""; CMPB_NAME="GSE83645"; onAddDataset(); __drain();
  el("cmpspecies").value="mouse";
  CMPB_TXT=""" + J(mat("GSE143688")) + r"""; CMPB_DES=""" + J(des("GSE143688")) + r"""; CMPB_NAME="GSE143688"; onAddDataset(); __drain();
  { const d=DSETS[0]; if(d.gfactor!=="Group") setDsFactor(d,"Group"); d.selA="PSO_non_lesional"; d.selB="PSO_lesional"; d.covSel=["Patient"]; d.resKey=""; }
  { const d=DSETS[1]; setDsFactor(d,"Type"); d.selA="non-lesion"; d.selB="lesion"; d.covSel=["Pair"]; d.resKey=""; }
  { const d=DSETS[2]; d.selA="uninvolved"; d.selB="psoriasis"; d.covSel=["Patient"]; d.resKey=""; }
  { const d=DSETS[3]; if(d.gfactor!=="Treatment") setDsFactor(d,"Treatment"); d.selA="Control"; d.selB="Aldara"; d.covSel=["Line","Day"]; d.resKey=""; }
"""
BODY = SETUP + r"""
  const R=[dsResult(0)].concat(DSETS.map((d,i)=>dsResult(i+1))), names=[GENES].concat(DSETS.map(d=>d.genes||d.geneNames||[]));
  const out={app:APP_VERSION, de:{}};
  R.forEach((r,i)=>{ const lab=i===0?"GSE54456":DSETS[i-1].label; let u=0,dn=0; const n=r.lfc.length;
    for(let g=0;g<n;g++) if(r.q[g]<=0.05&&Math.abs(r.lfc[g])>=1){ if(r.lfc[g]>0)u++; else dn++; }
    out.de[lab]={n, up:u, down:dn, lfc:Array.from(r.lfc), q:Array.from(r.q)}; });
  document.querySelectorAll = (sel)=> sel==="#metads input.metaon"? [0,1,2,3].map(i=>({checked:true, dataset:{i:String(i)}})) : [];
  // every meta-analysis control is set explicitly (the 0.23.0 baseline used DerSimonian-Laird with the z test)
  el("metamethod").value="REM"; el("metahk").value=""; el("metamink").value="3"; el("metafdr").value="0.05"; el("metalfc").value="1"; el("metapi").checked=false; runMeta();
  out.meta={method:META.method, tested:META.rows.length, sig:META.sig.length, up:META.sig.filter(r=>r.mu>0).length,
    rows:META.rows.map(r=>[GENES[r.g], r.mu, r.se, r.p, r.q, r.k]), sigGenes:META.sig.map(r=>GENES[r.g])};
  const sides=famAllSides().filter(Z=>Z.ok); out.fam=[]; sides.forEach(Z=>{ const s2=famAllScore(Z); ORTHO_FAM.forEach((f,k)=>{ const r=s2.res[k]; if(r) out.fam.push([k, Z.label, r.d, r.p, r.q, r.k]); }); });
  useBuiltin("hallmark_hs"); el("metaenrset").value="sig"; el("metaenrdir").value="up"; el("metaenrgo")._h["click"]();
  out.enr=enrResults.rows.map(r=>[r.name, r.m, r.k, r.p, r.q]);
  print("SNAP " + JSON.stringify(out));
"""
if __name__ == "__main__":
    o = run(APP, BODY, timeout=3000); S = lines_with(o, "SNAP ")
    if not S: print(o[-2000:]); sys.exit(1)
    json.dump(S[0], open(sys.argv[1], "w")); print("wrote", sys.argv[1], S[0]["app"], {k: (v["up"], v["down"]) for k, v in S[0]["de"].items()}, S[0]["meta"]["sig"])
