  // ---------- drive every analysis ----------
  const ERR=[]; const tryit=(lab,f)=>{ try{ f(); __drain(); }catch(e){ ERR.push(lab+": "+String(e).slice(0,160)); } };
  const click=id=>{ const h=el(id)._h&&el(id)._h.click; if(!h) throw "no click handler on "+id; h({preventDefault(){}, target:el(id)}); };
  el("selA").value=GROUPS[0]; el("selB").value=GROUPS[1]; el("selFDR").value="0.05"; el("fcthr").value="1"; deMethod="mod"; DEFIT="all"; DETREAT=false; DECOVSEL=[]; deKey="";
  for(const t of ["overview","gene","de","coexp","heat","enrich","compare","patterns","disco","traits"]) tryit("tab "+t, ()=>switchTab(t));
  tryit("gene", ()=>{ curGene=0; switchTab("gene"); drawGene(); for(const gp of ["box","violin","bar"]){ GPLOT=gp; drawGene(); } el("pwsel").value="all"; drawGene(); });
  tryit("compare", ()=>{ compareList=[0,1,2,3]; for(const ch of ["bar","box","violin"]){ CMPCHART=ch; drawCompare(); } el("cmpstat").value="all"; el("cmptest").value="welch"; CMPARR="cond_in_gene"; drawCompare(); });
  tryit("de", ()=>{ switchTab("de"); drawDE(); });
  tryit("allpairs", ()=>{ runAllPairs(); drawVolcanoGrid(); drawUpset(); });
  tryit("venn", ()=>{ VENNSTATE={picks:[{A:GROUPS[0],B:GROUPS[1],label:GROUPS[1]+" vs "+GROUPS[0]},{A:GROUPS[0],B:GROUPS[2],label:GROUPS[2]+" vs "+GROUPS[0]}], dir2:"any"}; drawVenn(); vennShowRegion(3); });
  tryit("modf", ()=>click("modfrun"));
  tryit("ix", ()=>{ ixRefreshUI(); if(GROUPS.length>=4){ el("ixA1").value=GROUPS[0]; el("ixB1").value=GROUPS[1]; el("ixA2").value=GROUPS[2]; el("ixB2").value=GROUPS[3]; runInteraction(); if(!IXRES) throw "interaction did not run"; } });
  tryit("cc", ()=>{ el("ccA1").value=GROUPS[0]; el("ccB1").value=GROUPS[1]; el("ccA2").value=GROUPS[0]; el("ccB2").value=GROUPS[2]; runCC(); });
  tryit("coexp", ()=>{ coGene=0; switchTab("coexp"); drawCoexp(); });
  tryit("heat", ()=>{ for(const v of ["var","de"]){ el("heatsrc").value=v; drawHeat(); } el("heatord").value="kmeans"; drawHeat(); el("heatord").value="hier"; });
  tryit("patterns", ()=>{ switchTab("patterns"); drawPatterns(); });
  useBuiltin(LIB); el("enrmin").value="5"; el("enrmax").value="2000";
  tryit("ora", ()=>{ enrMethod="ora"; switchTab("enrich"); el("enrsrc").value="de_up"; runEnrichment(); drawSetJaccard(); drawRidgeline(); });
  tryit("gsea", ()=>{ el("gseaperm").value="gene"; enrMethod="gsea"; runSync(runGSEAGen()); });
  tryit("fry", ()=>{ enrMethod="fry"; runSync(runFRYGen()); });
  tryit("enrmat", ()=>{ enrMethod="ora"; runSync(runEnrMatrixGen()); });
  tryit("disco", ()=>{ switchTab("disco"); for(const c of ["dscfam","dsccur","dscgmt"]) el(c).checked=true; click("dscrun"); for(let i=0;i<40;i++){ if(typeof drainMicrotasks==="function") drainMicrotasks(); __drain(); } if(!DISCO) throw "discovery screen did not run"; });
  tryit("hubs", ()=>{ el("hubM").value="800"; click("hubrun"); });
  tryit("traits", ()=>{ switchTab("traits"); const tc=traitColumns(); if(tc.length){ el("trsel").value=tc[0]; click("trrun"); for(let i=0;i<40;i++){ if(typeof drainMicrotasks==="function") drainMicrotasks(); __drain(); } if(!TRAITRES) throw "trait analysis did not run"; } });
  tryit("power", ()=>{ switchTab("overview"); drawPower(); });
  tryit("loadings", ()=>drawLoadings());

  // ---------- scans ----------
  const BAD=/\bNaN\b|\bundefined\b|\bInfinity\b|\[object Object\]|>null<|: null\b/;
  const textHits=[];
  for(const id of Object.keys(elCache)){ const c=elCache[id]._c||""; const m=BAD.exec(c); if(m) textHits.push([id, c.replace(/<[^>]+>/g," ").replace(/\s+/g," ").slice(Math.max(0,m.index-140>0? 0 : 0), 0)+(c.replace(/<[^>]+>/g," ").replace(/\s+/g," ").match(/.{0,70}(NaN|undefined|Infinity|\[object Object\]|null).{0,50}/)||[""])[0]]); }
  const last={}; plotlyCalls.forEach(x=>{ last[x.id]=x; });
  const plotHits=[], statsMissing=[], originHits=[], exportHits=[];
  const nonfin=v=>typeof v==="number"&&!isFinite(v);
  for(const id of Object.keys(last)){ const c=last[id]; const T=c.traces||[], Lay=c.layout||{}; if(!T.length) continue;
    // non-finite values that are drawn or shown on hover
    T.forEach((t,k)=>{ for(const f of ["x","y","z"]){ const a=t[f]; if(!Array.isArray(a)) continue; const flat=f==="z"? a.flat() : a; const n=flat.filter(nonfin).length; if(n) plotHits.push([id,k,t.type||"scatter",t.name||"",f+" non-finite",n,flat.length]); }
      const mc=t.marker&&t.marker.color; if(Array.isArray(mc)){ const n=mc.filter(v=>typeof v==="number"&&!isFinite(v)).length; if(n) plotHits.push([id,k,t.type||"scatter",t.name||"","marker.color non-finite",n,mc.length]); }
      const ms=t.marker&&t.marker.size; if(Array.isArray(ms)){ const n=ms.filter(v=>typeof v==="number"&&!isFinite(v)).length; if(n) plotHits.push([id,k,t.type||"scatter",t.name||"","marker.size non-finite",n,ms.length]); }
      if(Array.isArray(t.customdata)&&typeof t.hovertemplate==="string"&&/customdata/.test(t.hovertemplate)){ const n=t.customdata.flat().filter(v=>nonfin(v)||v==="NaN").length; if(n) plotHits.push([id,k,t.type||"scatter",t.name||"","hover customdata non-finite",n,t.customdata.length]); }
      if(Array.isArray(t.text)){ const n=t.text.filter(v=>typeof v==="string"&&/\bNaN\b|undefined/.test(v)).length; if(n) plotHits.push([id,k,t.type||"scatter",t.name||"","text NaN/undefined",n,t.text.length]); } });
    // statistics drawn on the figure but absent from its data export
    let sh=null; try{ sh=plotToSheets(T, Lay).sheets; }catch(e){ exportHits.push([id,"plotToSheets threw "+e]); continue; }
    const annTxt=(Lay.annotations||[]).map(a=>String(a.text||"")).join(" | ");
    const tickTxt=[].concat(...["xaxis","yaxis"].map(a=>(Lay[a]&&Lay[a].ticktext)||[])).map(String).join(" | ");
    const hasStatTxt=/(^|[\s(|])(\*{1,4}|ns)(\s|$|\|)|FDR\s*[\d<=≤]|\bp\s*[=<≤]\s*[\d.]|\bq\s*[=<≤]\s*[\d.]|\bρ\s*=|\br\s*=\s*-?[\d.]/.test(annTxt+" | "+tickTxt);
    if(hasStatTxt && !sh.some(x=>x.name==="Statistics")) statsMissing.push([id, (annTxt+" | "+tickTxt).replace(/<[^>]+>/g,"").slice(0,160)]);
    const P=sh.find(x=>x.name==="Prism"); if(!P||P.rows.length<2) exportHits.push([id,"no Prism rows", P? P.rows.length : 0]);
    // origin: numeric axis whose data are all >= 0 (or all <= 0) but 0 is not inside the axis range
    for(const ax of ["x","y"]){ const A=Lay[ax+"axis"]||{}; if(A.type==="log"||A.type==="category"||A.visible===false) continue;
      const vals=[].concat(...T.filter(t=>t.type!=="heatmap"&&!(t.xaxis&&t.xaxis!=="x")&&!(t.yaxis&&t.yaxis!=="y")).map(t=>(Array.isArray(t[ax])? t[ax] : []))).filter(v=>typeof v==="number"&&isFinite(v));
      if(vals.length<2) continue; const mn=Math.min(...vals), mx=Math.max(...vals);
      if(Array.isArray(A.tickvals)&&Array.isArray(A.ticktext)&&A.ticktext.some(x=>isNaN(Number(String(x).replace(/<[^>]+>/g,""))))) continue;   // categorical positions
      const r=A.range; let shows0;
      if(Array.isArray(r)) shows0= Math.min(r[0],r[1])<=0 && Math.max(r[0],r[1])>=0;
      else shows0= mn<=0&&mx>=0 || A.rangemode==="tozero";
      if(!shows0) originHits.push([id, ax, +mn.toPrecision(3), +mx.toPrecision(3), Array.isArray(r)? r.map(v=>+(+v).toPrecision(3)) : (A.rangemode||"auto"), String((A.title&&(A.title.text||A.title))||"").replace(/<[^>]+>/g,"").slice(0,40)]); }
  }
  const snap=()=>{ const o={}; for(const id of Object.keys(elCache)) o[id]=elCache[id]._c||""; const p={}; plotlyCalls.forEach(x=>{ p[x.id]=JSON.stringify(x.traces).length+":"+JSON.stringify(x.traces).slice(0,4000); }); return {o,p}; };
  const before=snap();
  // ---------- staleness: exclude one sample, revisit every tab ----------
  const j0=S.findIndex((x,i)=>included[i]&&x.group===GROUPS[1]);
  el("samplebar")._h["change"]({target:{dataset:{i:String(j0)}, checked:false, classList:{contains:()=>false}}}); __drain();
  for(const t of ["overview","gene","de","coexp","heat","enrich","compare","patterns","disco","traits"]) tryit("tab2 "+t, ()=>switchTab(t));
  const after=snap();
  const WATCH=["enrtable","enrinfo","pairsout","allpairs","venngenes","modfinfo","modftable","ixinfo","ixtable","ccinfo","dscinfo","dscout","dsctable","hubinfo","hubtable","trinfo","trtable","pwtable","genetable","detable","copos","heatnote","patlist","tfout","powerout","pcaloadout"];
  const PLOTS=["p_enrlolli","p_ridge","p_setjaccard","p_venn","p_upset","p_vgrid","p_modf","p_ixvolcano","p_cc","p_hubs","p_trvolc","p_heat","p_pat","p_gene","p_volcano","p_compare","p_coexp","p_pca","p_dsc","p_enrmat","p_es"];
  const stale=[]; WATCH.forEach(id=>{ if(before.o[id]!==undefined&&before.o[id]!==""&&before.o[id]===after.o[id]) stale.push(["el",id, before.o[id].replace(/<[^>]+>/g," ").replace(/\s+/g," ").slice(0,90)]); });
  PLOTS.forEach(id=>{ if(before.p[id]!==undefined&&before.p[id]===after.p[id]) stale.push(["plot",id]); });
  const banners=[]; for(const id of Object.keys(elCache)){ if(/^<div class="stalebar"/.test(elCache[id]._c||"")) banners.push(id); }
  const existing=Object.keys(last).sort();
  print("AUD " + JSON.stringify({ERR, textHits, plotHits, statsMissing, originHits, exportHits, stale, banners, plots:existing, nPlots:existing.length}));
