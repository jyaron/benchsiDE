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


if __name__ == "__main__":
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    fails = 0
    for name, fn in CHECKS.items():
        if only and name != only: continue
        try: ok, info = fn()
        except Exception as e: ok, info = False, {"error": str(e)[:500]}
        fails += not ok; print(("PASS " if ok else "FAIL ") + name + " " + json.dumps(info)[:600])
    sys.exit(1 if fails else 0)
