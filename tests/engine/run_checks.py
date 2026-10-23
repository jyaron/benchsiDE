"""Engine checks: run the application script in JavaScriptCore (macOS) with a stub page and stub Plotly.
This is not a browser test (see tests/*.spec.js for Playwright); it tests computation and plot specifications.
Usage: python3 tests/engine/run_checks.py [--only NAME]"""
import json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
exec(open(os.path.join(HERE, "harness.py")).read())
SRC = open(os.path.join(ROOT, "index.html"), encoding="utf-8").read(); APP = app_of(SRC)
CMP = open(os.path.join(ROOT, "validation", "compare.js")).read()
def mat(acc): return open(os.path.join(HERE, "data", acc + ".tsv")).read()
def des(acc): return open(os.path.join(HERE, "designs", acc + "_design.csv")).read()
CHECKS = {}
def check(fn): CHECKS[fn.__name__] = fn; return fn
def J(x): return json.dumps(x)

@check
def core():
    m = open(os.path.join(HERE, "..", "..", "demo", "GSE63310_counts.tsv")).read(); d = open(os.path.join(ROOT, "demo", "GSE63310_design.tsv")).read()
    ref = open(os.path.join(ROOT, "validation", "gse63310_reference.json")).read()
    o = run(APP, r"""
  const res = runExternalValidation(""" + J(m) + ", " + J(d) + ", " + ref + r""");
  print("RES " + JSON.stringify({ext:[res.rows.filter(x=>x.pass).length, res.rows.length]}));
  const st=runSelfTest(); __drain(); print("RES " + JSON.stringify({self:[st.rows.filter(x=>x.pass).length, st.rows.length]}));
  let errs=[]; for(const t of ["overview","gene","de","coexp","heat","enrich","compare","patterns","disco","traits"]){ try{ switchTab(t); __drain(); }catch(e){ errs.push(t+": "+e); } }
  print("RES " + JSON.stringify({tabErrors:errs}));""", CMP, timeout=900)
    r = {}; [r.update(x) for x in lines_with(o, "RES ")]
    ok = r.get("ext", [0, 1])[0] == r.get("ext", [0, 1])[1] and r.get("self", [0, 1])[0] == r.get("self", [0, 1])[1] and r.get("tabErrors") == []
    return ok, r

@check
def interaction():
    import gzip
    ref = json.loads(gzip.decompress(open(os.path.join(ROOT, "validation", "freeze", "interaction_reference_GSE143688.json.gz"), "rb").read()))
    o = run(APP, r"""
  el("species").value="mouse"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE143688_all")) + ", " + J(des("GSE143688_all")) + r"""); applyFactor("Cell"); GROUPTYPE="auto"; analyzeNow(); __drain();
  const gi=groupsIdx(), out={ids:Array.from(IDS)};
  const cv=(cells)=>{ const v=new Float64Array(cells.length); v[cells.indexOf("Aldara_KO")]=1; v[cells.indexOf("Control_KO")]=-1; v[cells.indexOf("Aldara_WT")]=-1; v[cells.indexOf("Control_WT")]=1; return v; };
  const four=["Control_WT","Aldara_WT","Control_KO","Aldara_KO"];
  for(const [nm,cells,cov,meth] of [["mod_all_day",GROUPS.slice(),["Day"],"mod"],["voom_four_day",four,["Day"],"voom"]]){
    const r=lmContrastDE({cells, gi, cvecs:[cv(cells)], covNames:cov, method:meth}); out[nm]={logFC:Array.from(r.lfc), t:Array.from(r.tv), p:Array.from(r.p), df0:r.prior.df0}; }
  print("RES " + JSON.stringify(out));""", timeout=900)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-400:]}
    A = R[0]; info = {"genes_identical": A["ids"] == ref["genes"]}
    ok = info["genes_identical"]
    for nm in ["mod_all_day", "voom_four_day"]:
        a, b = A[nm], ref[nm]
        dl = max(abs(x - y) for x, y in zip(a["logFC"], b["logFC"])); dt = max(abs(x - y) for x, y in zip(a["t"], b["t"]))
        dp = max(abs(x - y) / max(abs(y), 1e-300) for x, y in zip(a["p"], b["p"])); dd = abs(a["df0"] - b["df0"][0] if isinstance(b["df0"], list) else a["df0"] - b["df0"])
        info[nm] = {"logFC": dl, "t": dt, "p_rel": dp, "df0": dd}
        ok = ok and dl < 1e-8 and dt < 1e-8 and dp < 1e-6 and dd < 1e-8
    return ok, info

PSO_SETUP = r"""
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
  document.querySelectorAll = (sel)=> sel==="#metads input.metaon"? [0,1,2,3].map(i=>({checked:true, dataset:{i:String(i)}})) : [];
  el("metamethod").value="REML"; el("metahk").value="hk"; el("metamink").value="3"; el("metafdr").value="0.05"; el("metalfc").value="1"; el("metapi").checked=false; runMeta();
"""

@check
def cross_dataset():
    # meta-analysis default, model-based intervals, homology-family dot plot over all datasets, meta-analysis ORA vs hypergeometric
    from scipy.stats import hypergeom
    CL = ["defb4a","serpinb4","pi3","s100a9","il36a","s100a8","sprr2c","s100a7","s100a12","serpinb3","cxcl8","lce3d","il36g","lcn2","ccl20","krt6a","krt16","il17a","il23a","krt77"]
    o = run(APP, PSO_SETUP + r"""
  const sg=new Set(META.sig.map(r=>r.g)), CL=""" + J(CL) + r""";
  const SEX=["xist","rps4y1","ddx3y","uty","kdm5d","eif1ay","usp9y","zfy","nlgn4y","txlngy"], sexIn=SEX.filter(g=>sg.has(nameToIdx[g])), sexTested=SEX.filter(g=>META.rows.some(r=>r.g===nameToIdx[g])).length;
  const mito=META.sig.map(r=>GENES[r.g]).filter(g=>/^MT-|^RP[LS]\d+[A-Z]?\d*$/.test(g));
  const ciOK=META.rows.every(r=>Math.abs((r.ciHi-r.mu)-(r.mu-r.ciLo))<1e-9 && r.ciHi>=r.mu && r.dfTest===r.k-1);
  const fam=[]; el("fdfams").value="SERPINB3, S100A7, S100A9, LCE3, KRT6A, DEFB4A"; drawFamDot(); const c=plotlyCalls.filter(x=>x.id==="p_famdot").pop();
  useBuiltin("hallmark_hs"); el("metaenrset").value="sig"; el("metaenrdir").value="up"; el("metaenrgo")._h["click"]();
  const E=enrResults;
  print("RES " + JSON.stringify({method:META.method, hk:META.hk, sig:META.sig.length, cl:CL.filter(g=>sg.has(nameToIdx[g])).length, ciOK, sexIn, sexTested, mito,
    famTraces:c? c.traces.filter(t=>t.showlegend).map(t=>t.name) : [], famRows:c? c.layout.yaxis.ticktext.length : 0,
    enr:{q:E.q, N:E.N, rows:E.rows.map(r=>[r.m, r.k, r.p])}}));""", timeout=2400)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]; dev = max(abs(hypergeom.sf(k - 1, r["enr"]["N"], m, r["enr"]["q"]) - p) / max(hypergeom.sf(k - 1, r["enr"]["N"], m, r["enr"]["q"]), 1e-300) for m, k, p in r["enr"]["rows"])
    info = {k: r[k] for k in ["method", "hk", "sig", "cl", "ciOK", "famRows", "sexIn", "sexTested", "mito"]}; info["famDatasets"] = len(r["famTraces"]); info["ora_rel_dev"] = dev
    ok = r["method"] == "REML" and r["hk"] == "hk" and r["cl"] == 20 and r["sexIn"] == [] and r["sexTested"] == 10 and r["mito"] == [] and r["ciOK"] and len(r["famTraces"]) == 5 and r["famRows"] >= 6 and dev < 1e-6
    return ok, info

@check
def figures_and_qc():
    # angled-label margins, very-low-depth exclusion, module heatmap layout, volcano selection, Prism export, session round trip
    o = run(APP, r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE171012")) + ", " + J(des("GSE171012")) + r"""); applyFactor("Status"); GROUPTYPE="auto"; analyzeNow(); __drain();
  const out={};
  out.vlow=QCVLOW.map(i=>S[i].id);
  S.forEach((x,i)=>{ included[i]=(pending.factors.CellType[x.id]==="bulk"); });
  const want=["Healthy","Psoriasis_PreTreatment","Psoriasis_SecukinumabTreatmentWeek2","Psoriasis_SecukinumabTreatmentWeek4","Psoriasis_SecukinumabTreatmentWeek12"];
  GROUPS.splice(0, GROUPS.length, ...want); deKey="";
  curGene=nameToIdx["serpinb4"]; out.over=[];
  for(const ang of ["45","-45","auto"]) for(const w of [620,1000,1460]){
    el("p_gene").clientWidth=w; LABELCFG={angle:ang,size:"auto",long:"full"}; plotlyCalls=[]; GPLOT="bar"; drawGene();
    const c=plotlyCalls.filter(x=>x.id==="p_gene").pop(), xa=c.layout.xaxis, m=c.layout.margin, tt=xa.ticktext||want, n=tt.length, pw=w-m.l-m.r, a=Math.abs(xa.tickangle||0)*Math.PI/180;
    let worst=-1e9; tt.forEach((t,i)=>{ const xi=m.l+(i+0.5)*pw/n, wl=textW(String(t).split("<br>")[0],12); const o2= (xa.tickangle||0)>0? xi+wl*Math.cos(a)-w : (xa.tickangle||0)<0? wl*Math.cos(a)-xi : 0; worst=Math.max(worst,o2); });
    out.over.push(Math.round(worst)); }
  if(!Plotly.toImage) Plotly.toImage=async()=>"data:image/png;base64,AAAA";
  figPreview=async(o2)=>false;
  const genes=["xrcc6","xrcc5","xrcc2","xrcc4","xrcc1","xrcc3"].map(g=>nameToIdx[g]);
  (async()=>{ await exportModuleHeatmap(genes, "t", "Figure"); const c=plotlyCalls.filter(x=>x.id==="p_panel_tmp").pop();
    out.mh={zmin:c.traces[0].zmin, zmax:c.traces[0].zmax, szmin:c.traces[1].zmin, szmax:c.traces[1].zmax, cbLen:c.traces[0].colorbar.len, nX:c.traces[0].x.length, nLab:(c.layout.xaxis.ticktext||[]).length, hdr:c.layout.annotations.map(a=>a.text)}; })();
  __drain();
  el("selA").value="Healthy"; el("selB").value="Psoriasis_PreTreatment"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DECOVSEL=[]; deKey="";
  el("dedrag").value="select"; switchTab("de"); __drain(); drawDE();
  const v=plotlyCalls.filter(x=>x.id==="p_volcano").pop(), tr=v.traces[0], pts=[]; for(let k=0;k<tr.x.length;k++) if(tr.x[k]>=3&&tr.y[k]>=6) pts.push({curveNumber:0, pointIndex:k, text:tr.text[k]});
  el("p_volcano")._on["plotly_selected"]({points:pts}); out.sel=[pts.length, DESEL?DESEL.size:0];
  const ps=plotToSheets(v.traces, v.layout); out.prism=[ps.sheets.length, ps.sheets[0].rows.length];
  const ss=JSON.parse(JSON.stringify(sessionState())); DESEL=null; applySession(ss); __drain(); out.sessSel=DESEL?DESEL.size:0; out.sessMeta=!!ss.meta; out.sessIx=!!ss.ix;
  print("RES " + JSON.stringify(out));""", timeout=2400)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]; mh = r.get("mh", {})
    ok = (r["vlow"] == ["GSM5216179"] and max(r["over"]) <= 1 and mh.get("zmin") == mh.get("szmin") and mh.get("zmax") == mh.get("szmax") and mh.get("zmin") == -mh.get("zmax")
          and mh.get("cbLen", 0) >= 90 and r["sel"][0] > 0 and r["sel"][0] == r["sel"][1] and r["prism"][1] > 0 and r["sessSel"] == r["sel"][1] and r["sessMeta"] and r["sessIx"])
    return ok, r

@check
def discovery_text():
    # captions of Discovery figure exports name the test; the covariate note appears when DE covariates are selected
    o = run(APP, r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(des("GSE121212")) + r"""); applyFactor("Group"); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="PSO_non_lesional"; el("selB").value="PSO_lesional"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DECOVSEL=["Patient"]; deKey="";
  S.forEach((x,i)=>{ included[i]=(x.group==="PSO_lesional"||x.group==="PSO_non_lesional"); });
  el("dscfam").checked=true; el("dsccur").checked=true; el("dscgmt").checked=true; el("dscmin").value="3"; el("dscmax").value="100"; el("dscperm").value="200";
  useBuiltin("hallmark_hs"); const H={};
  const _qs=document.querySelectorAll;
  document.querySelectorAll=function(sel){ const b=makeEl("_btn"); b.dataset={fig:"0",hm:"0",more:"0",g:"0",r:"0"}; b.addEventListener=(ev,fn)=>{ (H[sel]=H[sel]||[]).push(fn); }; b.closest=()=>makeEl("_box"); return [b]; };
  const CAPS=[]; composeGenePanels=function(rows,name,cap){ CAPS.push(cap); }; exportModuleHeatmap=function(rows,name,cap){ CAPS.push(cap); };
  discoveryRun(); __drain();
  (H['#dscout [data-fig]']||[]).forEach(f=>f()); (H['#dscout [data-hm]']||[]).forEach(f=>f());
  print("RES " + JSON.stringify({camera:!!(DISCO.rows[0]&&DISCO.rows[0].camera), caps:CAPS.map(c=>String(c).slice(0,400)), info:el("dscinfo").textContent}));""", timeout=900)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]; caps = r["caps"]
    info = {"camera": r["camera"], "n_captions": len(caps), "captions_name_test": [("CAMERA" in c) for c in caps], "covariate_note": "does not adjust for the covariates" in r["info"]}
    return (r["camera"] and len(caps) >= 2 and all(info["captions_name_test"]) and info["covariate_note"]), info

@check
def rc_fixes():
    # defects found in the 1.0.0-rc.1 check, kept fixed
    o = run(APP, r"""
  el("species").value="mouse"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE143688_all")) + ", " + J(des("GSE143688_all")) + r"""); applyFactor("Cell"); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="Aldara_WT"; el("selB").value="Aldara_KO"; el("fcthr").value="0"; deMethod="welch"; DECOVSEL=[]; deKey="";
  const r=computeDE(), g=nameToIdx["slc28a2b"]; const out={welchConstP:r.p[g]};
  drawBiotypes(); const b=plotlyCalls.filter(x=>x.id==="p_bio").pop(); out.msgNoBar= !!(b&&b.traces.length===0? (b.cfg&&b.cfg.displayModeBar===false) : true);
  print("RES1 " + JSON.stringify(out));""", timeout=900)
    o2 = run(APP, PSO_SETUP + r"""
  const bad=META.rows.filter(r=>r.tau2===0&&r.I2>0).length;
  el("mfgenes").value="DEFB4A, SERPINB4, KRT16"; drawMetaForest(); const f=plotlyCalls.filter(x=>x.id==="p_metaforest").pop(), xr=f.layout.xaxis.range;
  renderMetaLOO(); genMethods(); const mt=el("methodstext").textContent;
  print("RES2 " + JSON.stringify({i2bad:bad, forestSym:Math.abs(xr[0]+xr[1])<1e-9, looDir:META.loo.every(o=>isFinite(o.dirFrac)), methodsHK:/Hartung-Knapp t test on k − 1 degrees of freedom and Benjamini/.test(mt), noWald:!/a Wald z test and Benjamini/.test(mt)}));""", timeout=2400)
    R1, R2 = lines_with(o, "RES1 "), lines_with(o2, "RES2 ")
    if not R1 or not R2: return False, {"error": (o + o2)[-600:]}
    titles = "click a cell for the shared genes" not in APP and "click a ridge for its barcode)" not in APP
    info = dict(R1[0], **R2[0], titlesClean=titles)
    ok = info["welchConstP"] == 1 and info["msgNoBar"] and info["i2bad"] == 0 and info["forestSym"] and info["looDir"] and info["methodsHK"] and info["noWald"] and titles
    return ok, info


@check
def covariates_treat():
    # 1.0.0-rc.2: voom with covariates and TREAT against limma (contrasts.fit workflow), plus the user-facing text
    import gzip
    REFD = json.load(gzip.open(os.path.join(HERE, "de_rc2_reference.json.gz"), "rt"))
    sets = {"GSE63310": (open(os.path.join(HERE, "..", "..", "demo", "GSE63310_counts.tsv")).read(), open(os.path.join(HERE, "designs", "GSE63310_lane_rin_design.tsv")).read(), "mouse", "celltype"),
            "GSE186063": (mat("GSE186063"), des("GSE186063"), "human", "Type")}
    info = {}; ok = True
    for ds, (M, D, sp, fac) in sets.items():
        refs = REFD[ds]
        o = run(APP, r"""
  el("species").value=""" + J(sp) + r"""; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe"; loadFromText(""" + J(M) + ", " + J(D) + r"""); applyFactor(""" + J(fac) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  const CF=""" + J([r["cfg"] for r in refs]) + r""";
  print("IDS " + JSON.stringify(IDS));
  CF.forEach((c,i)=>{ el("selA").value=c.A; el("selB").value=c.B; el("selFDR").value="0.05"; el("fcthr").value=String(c.tr==null?1:c.tr); deMethod=c.m; DEFIT=c.fit; DECOVSEL=c.cov; DETREAT=c.tr!=null; deKey="";
    const r=computeDE(); print("R " + JSON.stringify(r.err? {i, err:r.err} : {i, t:Array.from(r.tv), p:Array.from(r.p)})); });
""", timeout=2400)
        ids = lines_with(o, "IDS ")[0]; res = lines_with(o, "R ")
        dmax = 0.0; pmax = 0.0
        for r in res:
            ref = refs[r["i"]]
            if "err" in r or ref["ids"] != ids: ok = False; info[ds + "_err"] = r.get("err", "gene order differs"); continue
            dmax = max(dmax, max(abs(a - b) for a, b in zip(r["t"], ref["t"])))
            pmax = max(pmax, max(abs(a - b) / max(abs(b), 1e-300) for a, b in zip(r["p"], ref["p"])))
        info[ds] = {"configs": len(res), "t": dmax, "p_rel": pmax}
        ok = ok and len(res) == len(refs) and dmax <= 1e-6 and pmax <= 1e-6
    # user-facing behaviour on the demo with lane
    o = run(APP, r"""
  el("species").value="mouse"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe"; loadFromText(""" + J(sets["GSE63310"][0]) + ", " + J(sets["GSE63310"][1]) + r"""); applyFactor("celltype"); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="Basal"; el("selB").value="LP"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="voom"; DEFIT="all"; DECOVSEL=["lane"]; DETREAT=true; deKey="";
  switchTab("de"); __drain(); drawDE(); const r=computeDE();
  let csv=null; const _d=downloadCSV; downloadCSV=(txt,fn)=>{ csv=txt; }; el("deexport")._h["click"](); downloadCSV=_d;
  const hdr=csv? csv.split("\n").find(l=>l.startsWith("gene,")) : "";
  const vc=plotlyCalls.filter(x=>x.id==="p_volcano").pop();
  const mt=deMethodsText(); const ss=JSON.parse(JSON.stringify(sessionState())); DETREAT=false; applySession(ss); __drain();
  // ranking-based and cross-dataset functions use the ordinary moderated t
  DETREAT=true; deKey=""; const rT=computeDE(); DETREAT=false; deKey=""; const r0=computeDE(); DETREAT=true; deKey="";
  let same=true; for(let g=0;g<nG;g++) if(rT.ebT[g]!==r0.tv[g]||rT.ebP[g]!==r0.p[g]) { same=false; break; }
  const ds0=dsResult(0); let dsStd=true; for(let g=0;g<nG;g++) if(ds0.tv[g]!==r0.tv[g]) { dsStd=false; break; }
  let below=0; for(let g=0;g<nG;g++) if(rT.q[g]<=0.05&&Math.abs(rT.lfc[g])<=1) below++;
  print("UI " + JSON.stringify({treat:!!r.treat, csvTreatCols:/t_treat,p_treat,FDR_treat/.test(hdr), csvRule:/TREAT test of \|log2FC\| > 1/.test(csv||""),
    yTitle:/TREAT/.test(vc.layout.yaxis.title), methodsTreat:/McCarthy DJ & Smyth GK/.test(mt), methodsVoomCov:/voom weights were estimated from the fit of this design/.test(mt),
    sessionTreat:DETREAT===true, ebStandard:same, dsStandard:dsStd, belowThreshold:below, rule:deRuleText(0.05,1)}));
""", timeout=900)
    U = lines_with(o, "UI ")
    if not U: return False, dict(info, error=o[-500:])
    u = U[0]; info["ui"] = u
    ok = ok and u["treat"] and u["csvTreatCols"] and u["csvRule"] and u["yTitle"] and u["methodsTreat"] and u["methodsVoomCov"] and u["sessionTreat"] and u["ebStandard"] and u["dsStandard"] and u["belowThreshold"] == 0
    return ok, info
@check
def venn_labels():
    # long contrast names are wrapped and anchored away from each other, so that the set labels cannot overlap
    o = run(APP, r"""
  const labs=["Psoriasis_T vs Healthy","Psoriasis_Secukinumab_TreatmentWeek12 vs Healthy","AVeryLongNameWithoutAnyBreakCharacters vs B"];
  const out={};
  for(const n of [2,3]){
    vennSets=function(){ return {sets2:labs.slice(0,n).map((_,k)=>new Set([1,2,3,k+10])), fdr:0.05, thr:1}; };
    VENNSTATE={picks:labs.slice(0,n).map(l=>({label:l})), dir2:"any"}; drawVenn();
    const c=plotlyCalls.filter(x=>x.id==="p_venn").pop(); const a=c.layout.annotations;
    out["n"+n]={anchors:a.map(x=>x.xanchor), maxLine:Math.max(...a.map(x=>Math.max(...x.text.replace(/<b>|<\/b>/g,"").split("<br>").map(l=>l.length)))),
      yTop:c.layout.yaxis.range[1]};
  }
  out.short=vennLabel("LP vs Basal");
  print("RES " + JSON.stringify(out));""", timeout=300)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = (r["n2"]["anchors"] == ["right", "left"] and r["n3"]["anchors"] == ["center", "right", "left"]
          and r["n2"]["maxLine"] <= 22 and r["n3"]["maxLine"] <= 22 and r["short"] == "LP vs Basal")
    return ok, r


@check
def angled_labels():
    # slanted group labels lean towards the side with room, so a long last group name is not hidden at the edge
    def one(groups, cfg, width):
        samp = [f"S{i:02d}" for i in range(10*len(groups))]
        m = "gene\t" + "\t".join(samp) + "\n" + "\n".join(f"G{g}\t" + "\t".join(str(50 + (g*7 + i*13) % 400) for i in range(len(samp))) for g in range(200))
        d = "sample\tgroup\n" + "\n".join(f"{x}\t{groups[i//10]}" for i, x in enumerate(samp))
        o = run(APP, r"""
  el("species").value="human"; el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(m) + ", " + J(d) + r"""); GROUPTYPE="categorical"; analyzeNow(); __drain();
  Object.assign(LABELCFG, """ + J(cfg) + r"""); el("p_gene").clientWidth=""" + str(width) + r"""; curGene=0; switchTab("gene"); __drain(); drawGene(); __drain();
  const c=plotlyCalls.filter(x=>x.id==="p_gene").pop(); print("RES " + JSON.stringify({ang:c.layout.xaxis.tickangle, r:c.layout.margin.r, l:c.layout.margin.l}));""", timeout=300)
        R = lines_with(o, "RES "); return R[0] if R else {"error": o[-300:]}
    late = ["Healthy", "Psoriasis_PreTreatment", "Psoriasis_SecukinumabTreatmentWeek2", "Psoriasis_SecukinumabTreatmentWeek4", "Psoriasis_SecukinumabTreatmentWeekTwelve"]
    early = ["A_Psoriasis_SecukinumabTreatmentWeekTwelve", "B_Healthy", "C_Pre", "D_Wk2", "E_Wk4"]   # groups are shown in sorted order
    a = one(late, {"angle": "45", "size": "auto", "long": "wrap"}, 680)
    b = one(early, {"angle": "45", "size": "auto", "long": "wrap"}, 680)
    ok = a.get("ang") == -45 and a.get("r", 999) <= 20 and b.get("ang") == 45
    return ok, {"long_last": a, "long_first": b}


MOUSE_SIG_SETUP = r"""
  el("species").value="mouse"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE143688")) + ", " + J(des("GSE143688")) + r"""); applyFactor("Treatment"); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="Control"; el("selB").value="Aldara"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=["Line","Day"]; deKey="";
  deCache=computeDE(); el("cmpspecies").value="human";
  CMPB_TXT=""" + J(mat("GSE121212")) + r"""; CMPB_DES=""" + J(des("GSE121212")) + r"""; CMPB_NAME="GSE121212"; onAddDataset(); __drain();
  CMPB_TXT=""" + J(mat("GSE83645")) + r"""; CMPB_DES=""" + J(des("GSE83645")) + r"""; CMPB_NAME="GSE83645"; onAddDataset(); __drain();
  { const d=DSETS[0]; if(d.gfactor!=="Group") setDsFactor(d,"Group"); d.selA="PSO_non_lesional"; d.selB="PSO_lesional"; d.resKey=""; }
  { const d=DSETS[1]; d.selA="uninvolved"; d.selB="psoriasis"; d.resKey=""; }
  el("metamethod").value="REML"; el("metahk").value="hk"; refreshCmpUI(); buildSigSelector();
"""

@check
def signature_transfer():
    # mouse imiquimod signature scored in two human psoriasis cohorts (two datasets: the Hartung-Knapp pooled test on 1 df is
    # expected not to reach p <= 0.05, which the check also confirms): scores recomputed independently, Hedges' g against the
    # formula used by metafor escalc("SMD"), pooling through the metafor-validated meta-analysis, family rescue of genes
    # without a one-to-one ortholog, verdicts, and determinism of the seeded random signatures
    import numpy as np
    from math import lgamma, exp, log
    o = run(APP, MOUSE_SIG_SETUP + r"""
  el("sigsel").value="de_up"; el("sigfam").checked=true; runSigTransfer(); const R=SIGRES, d=DSETS[1], r1=R.rows[1];
  const rowsU=[...new Set(r1.M.units.flatMap(u=>u.rows))];
  const sc=sigScoreUnits(d, r1.M.units, dsZ(d));
  const out={units:r1.M.units.map(u=>({rows:u.rows, sign:u.sign})), L:Object.fromEntries(rowsU.map(r=>[r, Array.from(d.L.slice(r*d.nS,(r+1)*d.nS))])),
    sc:Array.from(sc), A:d.gi[d.selA], B:d.gi[d.selB], g:r1.g, v:r1.v,
    eff:R.rows.map(r=>[r.g, r.v]), pooled:{mu:R.pooled.mu, check:metaGene(R.rows.map(r=>r.g), R.rows.map(r=>r.v), "REML", "hk").mu},
    verdicts:R.rows.map(r=>r.verdict), pEmp:R.rows.map(r=>r.pEmp), fams:R.rows.map(r=>r.M.nFamUnits), rep:R.pooled.replicates};
  { const c=plotlyCalls.filter(x=>x.id==="p_sigforest").pop(); out.leg=c.layout.legend; out.mt=c.layout.margin.t; }
  runSigTransfer(); out.det=SIGRES.rows.map(r=>r.pEmp);
  const fidx=sigFamIndex("m"); const g=[]; for(let i=0;i<nG;i++) if(deCache.q[i]<=0.05&&deCache.lfc[i]>=1&&DSETS[0].matchAtoB[i]<0&&fidx[GENES[i].split(" (")[0].toLowerCase()]!==undefined) g.push(i);
  sigSpec=()=>({genes:g, sign:g.map(()=>1), expect:1}); SIGLABELS.x="x"; el("sigsel").value="x";
  runSigTransfer(); out.famOn=SIGRES.rows.map(r=>r.verdict); el("sigfam").checked=false; runSigTransfer(); out.famOff=SIGRES.rows.length;
  print("RES " + JSON.stringify(out));""", timeout=1800)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]; nS_ = len(r["sc"]); sc = np.zeros(nS_)
    for u in r["units"]:
        acc = np.zeros(nS_)
        for row in u["rows"]:
            x = np.array(r["L"][str(row)]); acc += (x - x.mean()) / x.std()
        sc += u["sign"] * acc / len(u["rows"])
    sc /= len(r["units"]); dsc = float(np.abs(sc - np.array(r["sc"])).max())
    a = sc[r["A"]]; b = sc[r["B"]]; n1, n2 = len(a), len(b); m = n1 + n2 - 2
    sp = np.sqrt(((n1-1)*a.var(ddof=1) + (n2-1)*b.var(ddof=1)) / m); Jc = exp(lgamma(m/2) - 0.5*log(m/2) - lgamma((m-1)/2))
    g = Jc*(b.mean()-a.mean())/sp; v = 1/n1 + 1/n2 + g*g/(2*(n1+n2))
    dg, dv = abs(g - r["g"])/abs(g), abs(v - r["v"])/v
    ok = (dsc < 1e-12 and dg < 1e-12 and dv < 1e-12 and abs(r["pooled"]["mu"] - r["pooled"]["check"]) < 1e-12
          and all(x == "replicates" for x in r["verdicts"]) and not r["rep"] and r["det"] == r["pEmp"] and min(r["fams"]) > 0
          and all(x in ("replicates", "partial replication") for x in r["famOn"][:1]) and r["famOff"] == 0
          and r["leg"].get("y", 0) >= 1 and r["leg"].get("yanchor") == "bottom" and r["mt"] >= 56)
    return ok, {"score_max_diff": dsc, "g_rel": dg, "v_rel": dv, "verdicts": r["verdicts"], "pooled_replicates": r["rep"], "families": r["fams"],
                "deterministic": r["det"] == r["pEmp"], "family_only_with": r["famOn"], "family_only_without_scored": r["famOff"], "legend_above_plot": r["leg"].get("y", 0) >= 1}


@check
def export_full_names():
    # figure data exports give full pathway and group names, not the shortened or wrapped labels drawn on the axis
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE54456")) + ", " + J(des("GSE54456")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="normal_skin"; el("selB").value="Psoriasis_skin"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey=""; deCache=computeDE();
  const out={}; const chk=(id, names)=>{ const c=plotlyCalls.filter(x=>x.id===id).pop(); if(!c){ out[id]=null; return; }
    const sh=plotToSheets(c.traces, c.layout).sheets; const cells=[]; sh.forEach(S=>S.rows.forEach(r=>r.forEach(v=>{ if(typeof v==="string") cells.push(v); })));
    out[id]={ell:cells.filter(v=>v.indexOf("\u2026")>=0).length, full:names.filter(n=>cells.some(v=>v.indexOf(n)>=0)).length, of:names.length}; };
  useBuiltin("gobp_hs"); el("enrsrc").value="de_up"; el("enrmin").value="5"; el("enrmax").value="2000"; switchTab("enrich"); __drain(); runEnrichment(); __drain();
  const nm=enrResults.rows.slice(0,15).map(r=>r.name);
  chk("p_enrlolli", nm.filter(n=>n.length>44)); drawSetJaccard(); chk("p_setjaccard", nm.filter(n=>n.length>34));
  drawRidgeline(); chk("p_ridge", enrResults.rows.slice(0,12).map(r=>r.name).filter(n=>n.length>30));
  useBuiltin("hallmark_hs"); enrMethod="gsea"; runSync(runGSEAGen()); __drain(); drawLeadingEdge();
  chk("p_leadedge", enrResults.rows.slice(0,10).filter(r=>r.leIdx&&r.leIdx.length).map(r=>r.name).filter(n=>n.length>30));
  print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=1800)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = all(v is not None and v["ell"] == 0 and v["of"] > 0 and v["full"] == v["of"] for v in r.values())
    return ok, r


@check
def relabel_reorder():
    # renaming and reordering groups (session and comparison datasets), hiding the sample list, reordering datasets
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(des("GSE121212")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  const out={}; const cnt=()=>{ deKey=""; const r=computeDE(); let u=0,dn=0; for(let g=0;g<nG;g++) if(r.q[g]<=0.05&&Math.abs(r.lfc[g])>=1){ if(r.lfc[g]>0)u++; else dn++; } return [u,dn]; };
  el("selA").value="CTRL_healthy"; el("selB").value="PSO_lesional"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DECOVSEL=[];
  const de0=cnt(); moveGroup("PSO_lesional",-1); out.ren=renameSessionGroup("CTRL_healthy","Healthy skin");
  out.sameDE=JSON.stringify(cnt())===JSON.stringify(de0); out.A=el("selA").value; out.groups=GROUPS.slice();
  out.merge=!!renameSessionGroup("Healthy skin","PSO_lesional");
  const ss=JSON.parse(JSON.stringify(sessionState())); out.sess=ss.design.filter(x=>x.group==="Healthy skin").length;
  SBHIDE=true; buildSampleBar(); out.hidden=el("samplebar").style.display; SBHIDE=false;
  el("cmpspecies").value="human";
  CMPB_TXT=""" + J(mat("GSE186063")) + r"""; CMPB_DES=""" + J(des("GSE186063")) + r"""; CMPB_NAME="GSE186063"; onAddDataset(); __drain();
  CMPB_TXT=""" + J(mat("GSE83645")) + r"""; CMPB_DES=""" + J(des("GSE83645")) + r"""; CMPB_NAME="GSE83645"; onAddDataset(); __drain();
  const d=DSETS[0]; setDsFactor(d,"Type"); d.selA="non-lesion"; d.selB="lesion";
  const dr=()=>{ d.resKey=""; const r=dsResult(1); let u=0; for(let g=0;g<d.nG;g++) if(r.q[g]<=0.05&&Math.abs(r.lfc[g])>=1) u++; return u; };
  const ds0=dr(); out.dsRen=renameDsGroup(d,"lesion","Lesional skin"); moveDsGroup(d,"normal_skin",-1);
  out.dsSame=dr()===ds0; out.dsKeys=Object.keys(d.gi); out.dsB=d.selB; out.dsMerge=!!renameDsGroup(d,"Lesional skin","non-lesion");
  { const d2=DSETS[1]; d2.selA="uninvolved"; d2.selB="psoriasis"; }
  document.querySelectorAll = (sel)=> sel==="#metads input.metaon"? [0,1,2].map(i=>({checked:true, dataset:{i:String(i)}})) : [];
  el("metamethod").value="REML"; el("metahk").value="hk"; el("metamink").value="2"; el("metafdr").value="0.05"; el("metalfc").value="1"; el("metapi").checked=false; runMeta();
  const lab0=META.idxs.map(i=>dsLabel(i)); moveDataset(0,1);
  out.moveOrder=DSETS.map(x=>x.label); out.metaLabelsKept=JSON.stringify(META.idxs.map(i=>dsLabel(i)))===JSON.stringify(lab0);
  print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=1800)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = (r["ren"] is None and r["sameDE"] and r["A"] == "Healthy skin" and r["groups"].index("PSO_lesional") < r["groups"].index("Healthy skin")
          and r["merge"] and r["sess"] == 37 and r["hidden"] == "none" and r["dsRen"] is None and r["dsSame"]
          and r["dsKeys"] == ["Lesional skin", "normal_skin", "non-lesion"] and r["dsB"] == "Lesional skin" and r["dsMerge"]
          and r["moveOrder"] == ["GSE83645", "GSE186063"] and r["metaLabelsKept"])
    return ok, r


@check
def compare_stats():
    # multi-gene comparison: brackets from the DE model equal the DE table; Welch/Holm brackets; Statistics sheet;
    # up to 30 genes; volcano selection survives redraws
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(des("GSE121212")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="CTRL_healthy"; el("selB").value="PSO_lesional"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey="";
  const res=computeDE(); const out={};
  const order=Array.from({length:nG},(_,g)=>g).sort((a,b)=>res.p[a]-res.p[b]); compareList=order.slice(0,25); CMPCHART="bar"; CMPARR="cond_in_gene";
  el("cmpstat").value="de"; el("cmptest").value="de"; el("cmpstatfmt").value="p"; switchTab("gene"); __drain(); drawCompare();
  let c=plotlyCalls.filter(x=>x.id==="p_compare").pop(), sr=c.layout.meta.statsRows;
  out.n=compareList.length; out.br=c.layout.annotations.length; out.qOK=compareList.every((g,i)=>sr[i+1][5]===res.q[g]&&sr[i+1][3]===res.lfc[g]);
  out.sheet=plotToSheets(c.traces, c.layout).sheets.some(x=>x.name==="Statistics");
  el("cmpstat").value="all"; el("cmptest").value="welch"; drawCompare(); c=plotlyCalls.filter(x=>x.id==="p_compare").pop(); sr=c.layout.meta.statsRows;
  const gi=groupsIdx(), g0=compareList[0];
  const rows0=sr.slice(1).filter(x=>x[0]===GENES[g0].split(" (")[0]);
  out.welch=rows0.map(x=>[x[1],x[2],x[4],x[5]]); out.vals=Object.fromEntries(Object.entries(gi).map(([k,ix])=>[k, ix.map(j=>L[g0*nS+j])]));
  // selection persistence
  el("dedrag").value="select"; switchTab("de"); __drain(); drawDE();
  const v=plotlyCalls.filter(x=>x.id==="p_volcano").pop(), tr=v.traces[0], pts=[]; for(let k=0;k<tr.x.length;k++) if(tr.x[k]>=3&&tr.y[k]>=15) pts.push({curveNumber:0, pointIndex:k});
  const V=el("p_volcano"); V._lastDraw=0; V._on["plotly_selected"]({points:pts}); const n0=DESEL?DESEL.size:0;
  switchTab("gene"); __drain(); switchTab("de"); __drain(); V._on["plotly_deselect"](); V._on["plotly_selected"]({points:[]});
  out.sel=[n0, DESEL?DESEL.size:0, el("deselbar").style.display];
  print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=1800)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    from scipy import stats as SS
    from statsmodels.stats.multitest import multipletests
    ps = [SS.ttest_ind(r["vals"][b], r["vals"][a], equal_var=False).pvalue for a, b, _, _ in r["welch"]]
    adj = multipletests(ps, method="holm")[1]
    dp = max(abs(x - w[2]) / x for x, w in zip(ps, r["welch"])); da = max(abs(x - w[3]) for x, w in zip(adj, r["welch"]))
    ok = r["n"] == 25 and r["br"] == 25 and r["qOK"] and r["sheet"] and dp < 1e-9 and da < 1e-9 and r["sel"][0] > 0 and r["sel"][1] == r["sel"][0] and r["sel"][2] == "flex"
    return ok, {"n": r["n"], "brackets": r["br"], "de_q_identical": r["qOK"], "stats_sheet": r["sheet"], "welch_p_rel": dp, "holm_abs": da, "selection": r["sel"]}


def _synth_counts(cols, n=2000, seed=1, eff=()):
    import random
    rnd = random.Random(seed); rows = ["Geneid\t" + "\t".join(cols)]
    for i in range(n):
        mu = rnd.gammavariate(0.6, 300)
        vals = []
        for j in range(len(cols)):
            m = mu * (4 if (j in eff and i < 100) else 1)
            vals.append(str(max(0, int(round(rnd.gauss(m, (m + m * m / 10) ** 0.5))))))
        rows.append(f"g{i+1}\t" + "\t".join(vals))
    return "\n".join(rows) + "\n"

@check
def sample_names_design():
    # repeated column names, group guessing from names, time factor chosen as the grouping, sample toggles recompute
    dup = ["m5", "m5", "m5", "m5+baf", "m5+baf", "m5+baf"]
    M1 = _synth_counts(dup, eff=(3, 4, 5))
    D1 = "sample\ttreatment\n" + "".join(f"{a}\t{b}\n" for a, b in zip(dup, ["ctrl"] * 3 + ["baf"] * 3))
    tc = [f"S{i}" for i in range(1, 13)]
    M2 = _synth_counts(tc, seed=2)
    D2 = "\t" + "\t".join(tc) + "\nStrain\t" + "\t".join(["A"] * 6 + ["B"] * 6) + "\nDay\t" + "\t".join(["0", "0", "4", "4", "10", "10"] * 2) + "\n"
    cases = ["WT_1,WT_2,KO_1,KO_2", "WT1,WT2,KO1,KO2", "S1_WT_rep1,S2_WT_rep2,S3_KO_rep1,S4_KO_rep2", "ctrl_rep1,ctrl_rep2,IMQ_rep1,IMQ_rep2", "GSM1,GSM2,GSM3,GSM4"]
    body = r"""
  el("species").value="mouse"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  const out={};
  loadFromText(""" + J(M1) + r""", null); GROUPTYPE="auto"; analyzeNow(); __drain();
  out.ids=S.map(x=>x.id); out.groups=GROUPS.slice(); out.warn=/occur more than once/.test(el("parseinfo").innerHTML);
  const lib=plotlyCalls.filter(x=>x.id==="p_lib").pop(); out.libDistinct=lib? new Set(lib.traces.map(t=>t.x).flat()).size : 0;
  loadFromText(""" + J(M1) + ", " + J(D1) + r"""); out.dGroups=pending.groups.slice();
  loadFromText(""" + J(M2) + ", " + J(D2) + r"""); out.dflt=FACTOR; applyFactor("Day"); GROUPTYPE="auto"; analyzeNow(); __drain();
  out.factor=FACTOR; out.numeric=NUMERIC; out.tGroups=GROUPS.slice(); out.cov=Object.keys(COVARS);
  el("selA").value="0"; el("selB").value="10"; deMethod="mod"; DEFIT="all"; DECOVSEL=[]; deKey=""; switchTab("de"); __drain();
  const p0=deCache.p[0], nA0=S.filter((x,i)=>included[i]&&x.group==="0").length;
  el("samplebar")._h["change"]({target:{dataset:{i:"0"}, checked:false, classList:{contains:()=>false}}}); __drain();
  out.toggle={off:!included[0], recomputed:deCache.p[0]!==p0, nA:[nA0, S.filter((x,i)=>included[i]&&x.group==="0").length]};
  applyFactor("Strain × Day"); out.combo=[...new Set(pending.groups)].sort(natCmp);
  out.infer={}; for(const c of """ + J(cases) + r"""){ const nm=c.split(","); out.infer[c]=[...new Set(inferGroups(uniqueSampleNames(nm), nm))]; }
  print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=900)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    exp_inf = {"WT_1,WT_2,KO_1,KO_2": ["WT", "KO"], "WT1,WT2,KO1,KO2": ["WT", "KO"], "S1_WT_rep1,S2_WT_rep2,S3_KO_rep1,S4_KO_rep2": ["WT", "KO"],
               "ctrl_rep1,ctrl_rep2,IMQ_rep1,IMQ_rep2": ["ctrl", "IMQ"], "GSM1,GSM2,GSM3,GSM4": ["GSM1", "GSM2", "GSM3", "GSM4"]}
    ok = (r["ids"] == ["m5_1", "m5_2", "m5_3", "m5+baf_1", "m5+baf_2", "m5+baf_3"] and r["groups"] == ["m5", "m5+baf"] and r["warn"] and r["libDistinct"] == 6
          and r["dGroups"] == ["ctrl"] * 3 + ["baf"] * 3 and r["dflt"] == "Strain × Day" and r["factor"] == "Day" and r["numeric"] and r["tGroups"] == ["0", "4", "10"] and r["cov"] == ["Strain"]
          and r["toggle"]["off"] and r["toggle"]["recomputed"] and r["toggle"]["nA"] == [4, 3]
          and r["combo"] == ["A · 0", "A · 4", "A · 10", "B · 0", "B · 4", "B · 10"] and all(sorted(r["infer"][k]) == sorted(v) for k, v in exp_inf.items()))
    return ok, r


@check
def two_factor():
    # combined grouping (Treatment x Day): per-level contrasts equal the DE tab's contrast; the interaction F with one
    # contrast equals t^2 of the Interaction card (checked against limma); a crossed time design defaults to the combination
    tc = [f"S{i}" for i in range(1, 13)]
    M2 = _synth_counts(tc, seed=3)
    D2 = "\t" + "\t".join(tc) + "\nStrain\t" + "\t".join(["A"] * 6 + ["B"] * 6) + "\nDay\t" + "\t".join(["0", "0", "4", "4", "10", "10"] * 2) + "\n"
    body = r"""
  el("species").value="mouse"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe"; const out={};
  loadFromText(""" + J(M2) + ", " + J(D2) + r"""); out.defaultFactor=FACTOR;
  loadFromText(""" + J(mat("GSE143688")) + ", " + J(des("GSE143688")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  regroupBy("Treatment × Day"); __drain(); const T2=twoFactor(); out.T=[T2.xName, T2.sName, T2.xLv, T2.sLv];
  el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey="";
  switchTab("de"); __drain(); el("tfref").value="Control"; el("tfcmp").value="Aldara"; runTwoFactor(); __drain();
  el("selA").value="Control · d7"; el("selB").value="Aldara · d7"; deKey=""; const de=computeDE(); const r7=TFRES.per[1];
  let dq=0; for(let g=0;g<nG;g++) dq=Math.max(dq, Math.abs(de.q[g]-r7.q[g]), Math.abs(de.lfc[g]-r7.lfc[g])); out.perLevelMaxDiff=dq;
  ixRefreshUI(); el("ixA1").value="Control · d3"; el("ixB1").value="Aldara · d3"; el("ixA2").value="Control · d7"; el("ixB2").value="Aldara · d7"; el("ixmeth").value="mod"; el("ixfit").value="all"; el("ixfdr").value="0.05"; el("ixlfc").value="0"; runInteraction();
  let mt=0; for(let g=0;g<nG;g++){ const t=(IXRES.t||IXRES.tv)[g]; mt=Math.max(mt, Math.abs(t*t-TFRES.IX.F[g])/Math.max(TFRES.IX.F[g],1e-12)); } out.tF=mt;
  curGene=0; switchTab("gene"); __drain(); drawGene(); const c=plotlyCalls.filter(x=>x.id==="p_gene").pop();
  out.lines=c.traces.filter(t=>t.meta==="summary").map(t=>t.name); out.sheets=plotToSheets(c.traces,c.layout).sheets.map(s2=>s2.name); out.labels=c.layout.annotations.length;
  print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=900)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = (r["defaultFactor"] == "Strain × Day" and r["T"] == ["Day", "Treatment", ["d3", "d7"], ["Aldara", "Control"]] and r["perLevelMaxDiff"] < 1e-12
          and r["tF"] < 1e-9 and r["lines"] == ["Aldara", "Control"] and "Statistics" in r["sheets"] and "Prism" in r["sheets"] and r["labels"] == 2)
    return ok, r


@check
def export_series_names():
    # every drawn figure's data export names its series; ridgeline Prism sheet = one column per gene set with all its genes
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(des("GSE121212")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="CTRL_healthy"; el("selB").value="PSO_lesional"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey=""; deCache=computeDE();
  for(const t of ["overview","gene","de","coexp","heat","patterns","disco","traits"]){ try{ switchTab(t); __drain(); }catch(e){} }
  try{ curGene=nameToIdx["pi3"]; switchTab("gene"); __drain(); drawGene(); for(const gp of ["box","violin","bar"]){ GPLOT=gp; drawGene(); } }catch(e){ print("ERR gene "+e); }
  try{ compareList=["pi3","s100a9","defb4a"].map(n=>nameToIdx[n]); for(const ch of ["bar","box","violin"]){ CMPCHART=ch; drawCompare(); } }catch(e){ print("ERR cmp "+e); }
  try{ useBuiltin("hallmark_hs"); el("enrsrc").value="de_up"; el("enrmin").value="5"; el("enrmax").value="2000"; switchTab("enrich"); __drain(); runEnrichment(); __drain(); drawSetJaccard(); drawRidgeline(); }catch(e){ print("ERR enr "+e); }
  try{ runAllPairs&&runAllPairs(); __drain(); }catch(e){ print("ERR pairs "+e); }
  try{ computeModF(); drawModF&&drawModF(); }catch(e){ print("ERR modf "+e); }
  const ids=[...new Set(plotlyCalls.map(x=>x.id))]; const rep=[];
  for(const id of ids){ const c=plotlyCalls.filter(x=>x.id===id).pop(); if(!c.traces||!c.traces.length) continue;
    let r; try{ r=plotToSheets(c.traces,c.layout); }catch(e){ rep.push([id,"THROW "+e]); continue; }
    const P=r.sheets.find(x=>x.name==="Prism"), A=r.sheets.find(x=>x.name==="All values");
    const hdr=P? P.rows[0] : []; const generic=(c.layout.meta&&c.layout.meta.prismRows)? 0 : hdr.filter((h,j)=>/^series \d+$/.test(String(h))||(j>0&&(h===""||h===null))).length;
    const blankSeries=A.rows.slice(1).filter(x=>x[0]===""||x[0]===null||/^series \d+$/.test(x[0])).length, nLong=A.rows.length-1;
    const types=[...new Set(c.traces.map(t=>t.type||"scatter"))];
    rep.push([id, r.kind||"none", generic, hdr.length, blankSeries, nLong, types.join("/")]); }
  
  const c=plotlyCalls.filter(x=>x.id==="p_ridge").pop(), P=plotToSheets(c.traces,c.layout).sheets[0];
  const ridge={hdr:P.rows[0], sets:RIDGEROWS.map(r=>r.name), counts:P.rows[0].map((_,j)=>P.rows.slice(1).filter(r=>r[j]!==null).length), nGenes:RIDGEROWS.map(r=>(r.allIdx||r.rows||r.ovIdx).length)};
  print("RES " + JSON.stringify({rep, ridge}));
"""
    o = run(APP, body, timeout=1200)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]; bad = [x for x in r["rep"] if isinstance(x[1], str) and x[1].startswith("THROW") or (len(x) > 4 and (x[2] > 0 or x[4] > 0))]
    rd = r["ridge"]
    ok = len(r["rep"]) >= 20 and not bad and rd["hdr"] == rd["sets"] and rd["counts"] == rd["nGenes"]
    return ok, {"figures": len(r["rep"]), "bad": bad, "ridge_sets": len(rd["sets"]), "ridge_counts_match": rd["counts"] == rd["nGenes"]}


@check
def gene_lists_detable():
    # Venn region -> heatmap / enrichment (identical to the same genes pasted), DE table sorted over all significant genes
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(des("GSE121212")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="CTRL_healthy"; el("selB").value="PSO_lesional"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey="";
  useBuiltin("hallmark_hs"); enrMethod="ora"; switchTab("enrich"); __drain();
  VENNSTATE={picks:[{A:"CTRL_healthy",B:"PSO_lesional",label:"P"},{A:"CTRL_healthy",B:"AD_lesional",label:"A"}], dir2:"any"}; drawVenn(); vennShowRegion(3);
  const genes=VENNSEL.genes2; setGeneList("P ∩ A", genes); refreshEnrSrc(); el("enrsrc").value="list"; el("enrmin").value="5"; el("enrmax").value="2000"; runEnrichment(); __drain();
  const a=enrResults.rows.map(r=>[r.name, r.p, r.q]);
  el("enrsrc").value="paste"; el("enrpaste").value=genes.map(g=>IDS[g]).join("\n"); runEnrichment(); __drain();
  const b=enrResults.rows.map(r=>[r.name, r.p, r.q]);
  const same=a.length===b.length && a.every((r,i)=>r[0]===b[i][0]&&r[1]===b[i][1]&&r[2]===b[i][2]);
  el("heatsrc").value="list"; switchTab("heat"); __drain(); const hc=plotlyCalls.filter(x=>x.id==="p_heat").pop();
  const hset=new Set(hc.traces[0].y), gset=new Set(genes.map(g=>GENES[g]));
  const out={n:genes.length, nSets:a.length, same, heatRows:hc.traces[0].y.length, heatSubset:[...hset].every(x=>gset.has(x))};
  // DE table sorting
  switchTab("de"); __drain(); drawDE();
  const rowsOf=()=>{ const h=el("detable").innerHTML; return [...h.matchAll(/data-g="(\d+)"/g)].map(m=>+m[1]); };
  const r0=rowsOf(); out.nShown=r0.length; out.title=el("detabtitle").textContent;
  const chk=(col,dir)=>{ DETAB.col=col; DETAB.dir=dir; renderDETable(); const r=rowsOf(); const v=r.map(g=>({gene:GENES[g].toLowerCase(), lfc:DETAB.lfc[g], abs:Math.abs(DETAB.lfc[g]), p:DETAB.p[g], q:DETAB.q[g]}[col]));
    let ok=true; for(let i=1;i<v.length;i++) if((v[i]<v[i-1]&&dir>0)||(v[i]>v[i-1]&&dir<0)) ok=false;
    // first shown must be the extreme over ALL significant genes
    const all=DETAB.rows.map(g=>({gene:GENES[g].toLowerCase(), lfc:DETAB.lfc[g], abs:Math.abs(DETAB.lfc[g]), p:DETAB.p[g], q:DETAB.q[g]}[col])); const ext=dir>0? all.reduce((a,b)=>a<b?a:b) : all.reduce((a,b)=>a>b?a:b);
    return ok && v[0]===ext; };
  out.sorts={}; for(const [c,d] of [["gene",1],["gene",-1],["lfc",-1],["lfc",1],["abs",-1],["p",1],["q",1],["q",-1]]) out.sorts[c+d]=chk(c,d);
  DETAB.shown=50; DETAB.col="q"; DETAB.dir=1; renderDETable(); const nb=DETAB.shown; DETAB.shown+=200; renderDETable(); out.more=[nb, rowsOf().length];
    print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=1200)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = r["same"] and r["heatRows"] == min(500, r["n"]) and r["heatSubset"] and all(r["sorts"].values()) and r["more"] == [50, 250]
    return ok, r


@check
def enrichment_dot_export():
    # ORA dot plot: no zero-overlap sets drawn; export in drawn order with overlap, p, FDR and genes; GSEA leading edges
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE54456")) + ", " + J(des("GSE54456")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="normal_skin"; el("selB").value="Psoriasis_skin"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey=""; deCache=computeDE();
  useBuiltin("hallmark_hs"); enrMethod="ora"; el("enrmin").value="5"; el("enrmax").value="2000"; switchTab("enrich"); __drain();
  const out={};
  for(const [lab,src,paste] of [["tiny","paste","TP53 CDKN1A MDM2 BAX"],["de","de_up",""]]){
    el("enrsrc").value=src; el("enrpaste").value=paste; runEnrichment(); __drain();
    const c=plotlyCalls.filter(x=>x.id==="p_enrlolli").pop(), sh=plotToSheets(c.traces,c.layout).sheets;
    const P=sh.find(x=>x.name==="Prism"), St=sh.find(x=>x.name==="Statistics"), G=sh.find(x=>x.name==="Genes");
    const exp=enrResults.rows.filter(r=>r.k>0).slice(0,15);
    out[lab]={zeros:c.traces[1].x.filter(v=>v===0).length, order:JSON.stringify(P.rows.slice(1).map(r=>r[0]))===JSON.stringify(exp.map(r=>r.name)),
      stats:St.rows.slice(1).every((r,i)=>r[2]===exp[i].k&&r[6]===exp[i].p&&r[7]===exp[i].q&&r[8].split("; ").length===exp[i].k), genes:G.rows.length-1===exp.reduce((a,r)=>a+r.k,0),
      first:St.rows[1].slice(1,8)};
  }
  el("gseaperm").value="gene"; enrMethod="gsea"; runSync(runGSEAGen()); __drain();
  const c=plotlyCalls.filter(x=>x.id==="p_enrlolli").pop(), sh=plotToSheets(c.traces,c.layout).sheets, St=sh.find(x=>x.name==="Statistics");
  out.gsea={nes:St.rows.slice(1).every((r,i)=>r[4]===enrResults.rows[i].NES), fdrText:St.rows[1][6], colour:c.traces[1].marker.colorbar.title.text};
  print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=1800)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = all(r[k]["zeros"] == 0 and r[k]["order"] and r[k]["stats"] and r[k]["genes"] for k in ("tiny", "de")) and r["gsea"]["nes"] and r["gsea"]["fdrText"].startswith("not reported") and r["gsea"]["colour"] == "−log₁₀ p"
    return ok, r


@check
def app_sweep():
    # every tab and analysis on two public datasets: no exceptions, no NaN/undefined in the page or in plotted
    # values, statistics drawn on a figure are in its data export, and every figure exports data
    sweep = open(os.path.join(HERE, "audit_sweep.js")).read()
    out = {}
    for ds, sp, lib in [("GSE121212", "human", "hallmark_hs"), ("GSE143688", "mouse", "hallmark_hs")]:
        body = r"""
  el("species").value=""" + J(sp) + r"""; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat(ds)) + ", " + J(_design_with_score(ds) if ds == "GSE121212" else des(ds)) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  const LIB=""" + J(lib) + ";\n" + sweep
        o = run(APP, body, timeout=3000)
        R = lines_with(o, "AUD ")
        if not R: return False, {"error": ds + " " + o[-600:]}
        r = R[0]
        # with two groups the sweep asks for a third group on purpose; the resulting messages are expected
        err = [e for e in r["ERR"] if not (ds == "GSE143688" and e.startswith("venn"))]
        txt = [t for t in r["textHits"] if not (ds == "GSE143688" and t[0] in ("venngenes",))]
        out[ds] = {"plots": r["nPlots"], "ERR": err, "text": txt, "plot": r["plotHits"], "statsMissing": r["statsMissing"], "export": r["exportHits"], "banners": r["banners"]}
    ok = all(not (v["ERR"] or v["text"] or v["plot"] or v["statsMissing"] or v["export"]) and v["plots"] >= 25 for v in out.values())
    return ok, out

@check
def stale_results():
    # after a sample is excluded: results that need Run are marked out of date until rerun; inexpensive ones recompute
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(des("GSE121212")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="CTRL_healthy"; el("selB").value="PSO_lesional"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey="";
  useBuiltin("hallmark_hs"); enrMethod="ora"; el("enrmin").value="5"; el("enrmax").value="2000"; switchTab("enrich"); __drain(); el("enrsrc").value="de_up"; runEnrichment(); __drain();
  el("hubM").value="600"; el("hubrun")._h.click(); __drain(); switchTab("de"); __drain(); runAllPairs(); computeModF(); drawModF(); __drain();
  const ban=id=>/^<div class="stalebar"/.test(el(id).innerHTML);
  const out={fresh:[ban("enrinfo"), ban("hubinfo")]}; const ap0=ALLPAIRS, mf0=MODF, nEnr0=enrResults.q;
  const j=S.findIndex(x=>x.group==="PSO_lesional");
  el("samplebar")._h["change"]({target:{dataset:{i:String(j)}, checked:false, classList:{contains:()=>false}}}); __drain();
  switchTab("de"); __drain(); switchTab("enrich"); __drain();
  out.after=[ban("enrinfo"), ban("hubinfo")]; out.autoPairs=ALLPAIRS!==ap0&&ALLPAIRS._sig===sampleSig(); out.autoModF=MODF!==mf0&&MODF._sig===sampleSig();
  runEnrichment(); __drain(); out.rerun=ban("enrinfo"); out.nEnr=[nEnr0, enrResults.q];
  el("selB").value="AD_lesional"; deKey=""; switchTab("de"); __drain(); switchTab("enrich"); __drain(); out.deChange=ban("enrinfo");
  print("RES " + JSON.stringify(out));"""
    o = run(APP, body, timeout=1500)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = r["fresh"] == [False, False] and r["after"] == [True, True] and r["autoPairs"] and r["autoModF"] and r["rerun"] is False and r["deChange"] is True
    return ok, r


@check
def enrichment_table_export():
    # the enrichment results table downloads (CSV and Excel) hold every tested set, not only the 60 shown, for ORA,
    # GSEA and FRY, with no NaN; ORA fold enrichment = overlap / expected; the workbook opens as a valid xlsx
    body = r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE54456")) + ", " + J(des("GSE54456")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="normal_skin"; el("selB").value="Psoriasis_skin"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey=""; deCache=computeDE();
  useBuiltin("gobp_hs"); el("enrmin").value="5"; el("enrmax").value="2000"; switchTab("enrich"); __drain();
  const out={}; let csv=null; const _d=downloadCSV; downloadCSV=(t,f)=>{ csv={t,f}; }; let xl=null;
  const grab=lab=>{ csv=null; let err=null; try{ el("enrtcsv")._h.click(); }catch(e){ err=String(e); }
    const T=enrTableSheets(); const b=xlsxBytes(T.sheets); if(lab==="ora") xl=Array.from(b);
    const body=csv? csv.t.split("\n").filter(l=>l&&l[0]!=="#") : [];
    out[lab]={err, dataLines:body.length-1, nRows:enrResults.rows.length, nan:csv? (csv.t.match(/NaN|undefined/g)||[]).length : -1, sheets:T.sheets.map(s=>s.name), bar:el("enrtbar").style.display}; };
  enrMethod="ora"; el("enrsrc").value="de_up"; runEnrichment(); __drain(); grab("ora");
  out.ratioErr=Math.max(...enrResults.rows.map(r=>Math.abs(r.ratio - r.k/(enrResults.q*r.m/enrResults.N))));
  el("gseaperm").value="gene"; enrMethod="gsea"; runSync(runGSEAGen()); __drain(); grab("gsea");
  enrMethod="fry"; runSync(runFRYGen()); __drain(); grab("fry");
  downloadCSV=_d;
  print("RES " + JSON.stringify(out)); print("XL " + JSON.stringify(xl));
"""
    o = run(APP, body, timeout=2400)
    R = [json.loads(l[4:]) for l in o.split("\n") if l.startswith("RES ")]
    X = [json.loads(l[3:]) for l in o.split("\n") if l.startswith("XL ")]
    if not R or not X: return False, {"error": o[-600:]}
    r = R[0]
    import io, zipfile as _zf
    z = _zf.ZipFile(io.BytesIO(bytes(X[0]))); names = z.namelist()
    ok = all(r[m]["err"] is None and r[m]["nRows"] > 60 and r[m]["dataLines"] == r[m]["nRows"] and r[m]["nan"] == 0 and r[m]["bar"] == "flex" for m in ("ora", "gsea", "fry"))
    ok = ok and r["ratioErr"] < 1e-9 and r["ora"]["sheets"] == ["Results", "Genes", "Notes"] and r["fry"]["sheets"] == ["Results", "Notes"] and "xl/workbook.xml" in names and z.testzip() is None
    return ok, {m: {k: r[m][k] for k in ("nRows", "dataLines", "nan", "sheets")} for m in ("ora", "gsea", "fry")} | {"ratioErr": r["ratioErr"], "xlsx_parts": len(names)}


def _design_with_score(ds):
    # adds a numeric per-sample measure (deterministic) so that the sample-trait analyses run
    lines = des(ds).rstrip("\n").split("\n"); n = len(lines[0].split(",")) - 1
    return "\n".join(lines + ["score," + ",".join("%.3f" % (10 + ((i * 7919) % 97) / 9.7) for i in range(n))]) + "\n"

def _adversarial_design(ds):
    # group names with a comma and with quotes, written as a spreadsheet program writes them
    out = []
    for l in _design_with_score(ds).rstrip("\n").split("\n"):
        if l.startswith("Group,"):
            v = l.split(",")
            l = ",".join([v[0]] + ['"AD, lesional"' if x == "AD_lesional" else '"CTRL ""healthy"""' if x == "CTRL_healthy" else x for x in v[1:]])
        out.append(l)
    return "\n".join(out) + "\n"

def _csv_problems(t):
    import csv as _csv, io as _io
    body = "\n".join(l for l in t.split("\n") if not l.startswith("#"))
    rows = [r for r in _csv.reader(_io.StringIO(body)) if r]
    if not rows: return {"empty": True}
    h = rows[0]; bad = re.compile(r"^(NaN|-?Infinity|undefined|null|\[object Object\])$|<span|<b>|&amp;|&lt;")
    return {"rows": len(rows) - 1, "ragged": sum(1 for r in rows[1:] if len(r) != len(h)),
            "bad": sum(1 for r in rows[1:] for v in r if bad.search(v)), "empty_header": sum(1 for x in h if not x.strip())}

@check
def export_audit():
    # every export button, after every analysis has been run (group and gene-set names containing commas and quotes):
    # each file is written without error, parses as CSV with one value per column, holds no NaN/undefined/markup;
    # every figure's Excel export is a valid workbook; after a sample is excluded, exports of results that need a
    # Run click are refused and inexpensive results are recomputed before export
    import io as _io, zipfile as _zf
    sweep = open(os.path.join(HERE, "audit_sweep.js")).read()
    rd = lambda f: open(os.path.join(HERE, f)).read()
    pre = r"""
  const _ub=useBuiltin; useBuiltin=function(k){ _ub(k); if(GMT) GMT=GMT.map((x,i)=>Object.assign({}, x, {name: i%2? x.name.replace("_",", ") : x.name+' "q"'})); };
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(_adversarial_design("GSE121212")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  const LIB="hallmark_hs";
""" + sweep.split("  // ---------- scans ----------")[0] + rd("export_extra.js")
    scen = {"GSE121212": pre, "cohorts": PSO_SETUP + rd("export_cohorts.js")}
    out = {"files": 0, "problems": [], "clickErrors": [], "stale_not_refused": [], "groups_parsed": None}
    for nm, body in scen.items():
        o = run(APP, body + rd("export_capture.js"), timeout=3600)
        F = [json.loads(l[6:]) for l in o.split("\n") if l.startswith("FILES ")]
        C = [json.loads(l[6:]) for l in o.split("\n") if l.startswith("CLICK ")]
        if not F or not C: return False, {"error": nm + " " + o[-600:]}
        out["clickErrors"] += C[0]["CLICKERR"] + C[0]["ERR"]
        for x in F[0]:
            out["files"] += 1
            if x["kind"] == "csv" and not x["f"].endswith(".gmt"):
                p = _csv_problems(x["t"])
                if x["trunc"]: p["ragged"] = max(0, p.get("ragged", 0) - 1)
                if p.get("ragged") or p.get("bad") or p.get("empty_header"): out["problems"].append([nm, x["f"], p])
            if x["kind"] == "bin" and x.get("bytes"):
                z = _zf.ZipFile(_io.BytesIO(bytes(x["bytes"])))
                if z.testzip() is not None or "xl/workbook.xml" not in z.namelist(): out["problems"].append([nm, x["f"], "invalid xlsx"])
        o = run(APP, body + rd("export_stale.js"), timeout=3600)
        S_ = [json.loads(l[6:]) for l in o.split("\n") if l.startswith("STALE ")]
        if not S_: return False, {"error": nm + " stale " + o[-600:]}
        RUNKEYS = ("enr", "emat", "hub", "net", "nb", "dsc", "tr", "cmpexp", "famexp", "cons", "meta", "sig", "ix", "tf")
        for b, nA, nB, same, after, thr, hidA in S_[0]["rep"]:
            if thr: out["clickErrors"].append(b + " " + thr)
            if nA and not hidA and b.startswith(RUNKEYS) and after not in ("blocked", "hidden") and same:
                out["stale_not_refused"].append([nm, b])
    # every figure's data workbook opens (adversarial names)
    o = run(APP, pre + r"""
  const ids=[...new Set(plotlyCalls.map(x=>x.id))]; let n=0; const bad=[];
  for(const id of ids){ const c=plotlyCalls.filter(x=>x.id===id).pop(); try{ const b=xlsxBytes(plotToSheets(c.traces, c.layout).sheets); if(!(b&&b.length>100)) bad.push(id); n++; }catch(e){ bad.push(id+": "+e); } }
  print("XLS " + JSON.stringify({n, bad, groups:GROUPS}));
""", timeout=3000)
    X = [json.loads(l[4:]) for l in o.split("\n") if l.startswith("XLS ")]
    if not X: return False, {"error": "xlsx " + o[-600:]}
    out["figures"] = X[0]["n"]; out["figure_errors"] = X[0]["bad"]; out["groups_parsed"] = X[0]["groups"]
    ok = (out["files"] >= 40 and not out["problems"] and not out["clickErrors"] and not out["stale_not_refused"]
          and out["figures"] >= 30 and not out["figure_errors"] and "AD, lesional" in out["groups_parsed"] and 'CTRL "healthy"' in out["groups_parsed"])
    return ok, out


@check
def kmeans_k():
    # heatmap and expression patterns: a k chosen by the user is used; auto k is unchanged (mean silhouette, k = 2-10);
    # the elbow plot's within-cluster sum of squares and mean silhouette equal independent computations from the
    # clusters drawn; the elbow plot shows 0 on both axes and its data export names both series
    import numpy as np
    o = run(APP, r"""
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(""" + J(mat("GSE121212")) + ", " + J(des("GSE121212")) + r"""); GROUPTYPE="auto"; analyzeNow(); __drain();
  const last=id=>plotlyCalls.filter(x=>x.id===id).pop(), out={};
  switchTab("heat"); __drain(); el("heatN").value="100"; el("heatord").value="kmeans"; el("heatK").value="5"; drawHeat();
  const he=last("p_heatelbow"); out.heat={nC:new Set(heatClusters.map(r=>r.cluster)).size, wrap:el("heatkwrap").style.display, wcss:he.traces[0].y, xr:he.layout.xaxis.range, yr:he.layout.yaxis.range, vline:he.layout.shapes[0].x0};
  switchTab("patterns"); __drain(); el("patK").value="auto"; drawPatterns(); out.autoK=patState.order.length; out.autoSil=autoK(patState.profiles).k;
  el("patK").value="8"; drawPatterns(); const pe=last("p_patelbow");
  out.pat={k:patState.order.length, asg:Array.from(patState.km.asg), P:patState.profiles, wcss:pe.traces[0].y, sil:pe.traces[1].y, silk:pe.traces[1].x, names:pe.traces.map(t=>t.name), sheet:plotToSheets(pe.traces, pe.layout).sheets[0].rows[0]};
  el("heatord").value="hier"; drawHeat(); out.hier=[el("heatkwrap").style.display, el("p_heatelbow").style.display];
  print("RES " + JSON.stringify(out));
""", timeout=1800)
    R = lines_with(o, "RES ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]; P = np.array(r["pat"]["P"]); a = np.array(r["pat"]["asg"])
    w = float(sum(((P[a == c] - P[a == c].mean(0)) ** 2).sum() for c in np.unique(a)))
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)); sil = []
    for i in range(len(P)):
        same = a == a[i]
        if same.sum() <= 1: sil.append(0.0); continue
        ai = D[i, same].sum() / (same.sum() - 1); bi = min(D[i, a == c].mean() for c in np.unique(a) if c != a[i])
        sil.append((bi - ai) / max(ai, bi))
    s_ind = float(np.mean(sil)); s_app = r["pat"]["sil"][r["pat"]["silk"].index(8)]
    res = {"heat_k": r["heat"]["nC"], "pat_k": r["pat"]["k"], "wcss_rel": abs(r["pat"]["wcss"][7] - w) / w, "sil_diff": abs(s_app - s_ind),
           "autoK": [r["autoK"], r["autoSil"]], "origin": [r["heat"]["xr"][0], r["heat"]["yr"][0]], "hier": r["hier"], "sheet": r["pat"]["sheet"]}
    ok = (r["heat"]["nC"] == 5 and r["heat"]["vline"] == 5 and r["heat"]["wrap"] == "inline" and r["pat"]["k"] == 8 and res["wcss_rel"] < 1e-9
          and res["sil_diff"] < 1e-6 and r["autoK"] == r["autoSil"] and res["origin"] == [0, 0] and r["hier"] == ["none", "none"]
          and all(np.diff(r["pat"]["wcss"]) <= 1e-9) and "mean silhouette" in res["sheet"] and "within-cluster sum of squares" in res["sheet"])
    return ok, res


GROUP_ORDER_JS = r'''
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(@MAT121212@, @DES121212@); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="CTRL_healthy"; el("selB").value="PSO_lesional"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey=""; const r0=computeDE();
  const out={g0:GROUPS.slice()};
  const tgt=orderAfterDrop(GROUPS, "PSO_lesional", GROUPS[0], false); out.tgt=tgt;
  setGroupOrder(tgt); out.g1=GROUPS.slice(); out.sel=[el("selA").value, el("selB").value];
  deKey=""; const r1=computeDE(); let md=0; for(let g=0;g<nG;g++) md=Math.max(md, Math.abs(r1.tv[g]-r0.tv[g]), Math.abs(r1.q[g]-r0.q[g])); out.deDiff=md;
  GROUPORDER=null; deriveGroups(S.map(x=>x.group)); while(GROUPS.indexOf("PSO_lesional")>0) moveGroup("PSO_lesional",-1); out.arrows=GROUPS.slice();
  setGroupOrder(orderAfterDrop(GROUPS, "AD_lesional", GROUPS[GROUPS.length-1], true)); out.g2=GROUPS.slice();
  out.bad=setGroupOrder(["x"]); out.sess=sessionState().groupOrder;
  el("cmpspecies").value="human"; CMPB_TXT=@MAT83645@; CMPB_DES=@DES83645@; CMPB_NAME="GSE83645"; onAddDataset(); __drain();
  const d=DSETS[0]; out.d0=(d.gorder||Object.keys(d.gi)).slice(); setDsGroupOrder(d, orderAfterDrop(out.d0, out.d0[1], out.d0[0], false)); out.d1=d.gorder.slice(); out.dsel=[d.selA,d.selB];
  el("species").value="mouse"; useBuiltinAnno();
  loadFromText(@MAT143688@, @DESNUM@); applyFactor("Day"); GROUPTYPE="auto"; analyzeNow(); __drain();
  out.num={NUMERIC, groups:GROUPS.slice(), refused:setGroupOrder(GROUPS.slice().reverse())===false, after:GROUPS.slice()};
  print("R " + JSON.stringify(out));
'''

@check
def group_order():
    # group order set by dragging (orderAfterDrop + setGroupOrder) equals the same order reached with the arrows; no
    # statistic changes; the contrast is kept; the order is saved in the session; comparison datasets likewise;
    # numeric (time) groups are not reordered
    lines = des("GSE143688").rstrip("\n").split("\n")
    dn = "\n".join((l.replace(",d3", ",3").replace(",d7", ",7") if l.split(",")[0] == "Day" else l) for l in lines) + "\n"
    o = run(APP, GROUP_ORDER_JS.replace("@MAT121212@", J(mat("GSE121212"))).replace("@DES121212@", J(des("GSE121212"))).replace("@MAT83645@", J(mat("GSE83645"))).replace("@DES83645@", J(des("GSE83645"))).replace("@MAT143688@", J(mat("GSE143688"))).replace("@DESNUM@", J(dn)), timeout=1800)
    R = lines_with(o, "R ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = (r["g1"] == r["tgt"] == r["arrows"] and r["g1"][0] == "PSO_lesional" and r["deDiff"] == 0 and r["sel"] == ["CTRL_healthy", "PSO_lesional"]
          and r["g2"][-1] == "AD_lesional" and r["bad"] is False and r["sess"] == r["g2"] and r["d1"] == r["d0"][::-1] and r["dsel"] == ["psoriasis", "uninvolved"]
          and r["num"]["NUMERIC"] is True and r["num"]["refused"] is True and r["num"]["after"] == r["num"]["groups"])
    return ok, {k: r[k] for k in ("g1", "g2", "deDiff", "d1", "num")}


CAT_AXES_JS = r'''
  el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
  loadFromText(@MAT@, @DES@); GROUPTYPE="auto"; analyzeNow(); __drain();
  el("selA").value="normal_skin"; el("selB").value="Psoriasis_skin"; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey=""; deCache=computeDE();
  useBuiltin("gobp_hs"); el("enrmin").value="5"; el("enrmax").value="2000";
  const a={};
  switchTab("heat"); __drain(); el("heatsrc").value="de"; el("heatN").value="100"; el("heatord").value="kmeans"; el("heatK").value="auto"; drawHeat(); __drain();
  const byC={}; heatClusters.forEach(r=>(byC[r.cluster]=byC[r.cluster]||[]).push(r.gidx)); a.nClusters=Object.keys(byC).length;
  HEATENRBG="all"; renderHeatEnrichment(); const big=Object.keys(byC).sort((x,y)=>byC[y].length-byC[x].length)[0];
  const mine=heatEnr.filter(r=>r.cluster===big).slice(0,5), ref=oraRows(byC[big]).filter(r=>r.k>=2).slice(0,5);
  a.same=mine.length===5&&mine.every((r,i)=>r.name===ref[i].name&&r.p===ref[i].p&&r.q===ref[i].q);
  HEATENRBG="clustered"; renderHeatEnrichment(); a.untestable=byC[big].length/heatClusters.length<0.9 || /no set can be enriched/.test(el("heatenr")._c||""); HEATENRBG="all";
  switchTab("enrich"); __drain(); el("enrsrc").value="de_up"; runEnrichment(); __drain();
  const c=plotlyCalls.filter(x=>x.id==="p_enrlolli").pop(); a.lolli={xr:c.layout.xaxis.rangemode, xz:c.layout.xaxis.zeroline, yr:c.layout.yaxis.range, yz:c.layout.yaxis.zeroline};
  a.heatXz=plotlyCalls.filter(x=>x.id==="p_heat").pop().layout.xaxis.zeroline;
  curGene=nameToIdx["pi3"]; switchTab("gene"); __drain(); drawGene(); a.geneXz=plotlyCalls.filter(x=>x.id==="p_gene").pop().layout.xaxis.zeroline;
  print("R " + JSON.stringify(a));
'''

@check
def axes_cluster_enrichment():
    # category axes have no zero line through the first category; the enrichment dot plot shows 0 and pads its rows;
    # heatmap cluster enrichment against all filtered genes equals the Enrichment-tab test of the same gene list, and
    # against the clustered genes a cluster holding >= 90% of them is reported as untestable
    o = run(APP, CAT_AXES_JS.replace("@MAT@", J(mat("GSE54456"))).replace("@DES@", J(des("GSE54456"))), timeout=1500)
    R = lines_with(o, "R ")
    if not R: return False, {"error": o[-600:]}
    r = R[0]
    ok = (r["lolli"] == {"xr": "tozero", "xz": True, "yr": [-0.6, 14.6], "yz": False} and r["heatXz"] is False and r["geneXz"] is False
          and r["same"] and r["untestable"] and r["nClusters"] == 2)
    return ok, r


if __name__ == "__main__":
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    fails = 0
    for name, fn in CHECKS.items():
        if only and name != only: continue
        try: ok, info = fn()
        except Exception as e: ok, info = False, {"error": str(e)[:500]}
        fails += not ok; print(("PASS " if ok else "FAIL ") + name + " " + json.dumps(info)[:600])
    sys.exit(1 if fails else 0)
