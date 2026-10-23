// export, exclude one sample, export again: results that need a Run click must be refused

  const SNAPF=(tag)=>{ const files={}; const _dc=downloadCSV, _dt=downloadText, _sb=saveBlob; let CUR="";
    const put=(t)=>{ (files[CUR]=files[CUR]||[]).push(t); };
    downloadCSV=(t,f)=>put(String(t)); downloadText=(t,f)=>put(String(t)); saveBlob=(b,f)=>{ const p=(b&&b.p)||[]; put(p.every(x=>typeof x==="string")? p.join("") : "[bin "+(p[0]&&p[0].length)+"]"); };
    for(const id of Object.keys(elCache)){ const fn=(elCache[id]._h&&elCache[id]._h.click)||elCache[id].onclick; if(typeof fn!=="function") continue; const src=String(fn);
      if(!/download|saveBlob|export|Export|xlsx|CSV|GMT|cyto/i.test(src)) continue;
      if(/^(enrrun|hubrun|dscrun|trrun|modfrun|savesess|reportbtn|sumfigbtn|selftestbtn|panelfig|netdraw|famrun|trcamrun|ixenr|metaenrgo)$/.test(id)) continue;
      if(elCache[id].style&&elCache[id].style.display==="none") { (files[id]=files[id]||[]).push("HIDDEN"); continue; }
      CUR=id; const _al=alert; alert=(m)=>put("BLOCKED "+m); try{ fn({preventDefault(){}, target:elCache[id], stopPropagation(){}}); __drain(); }catch(e){ put("THROW "+e); } alert=_al; }
    downloadCSV=_dc; downloadText=_dt; saveBlob=_sb; return files; };
  const hash=t=>{ let h=0; for(let i=0;i<t.length;i++){ h=(h*31+t.charCodeAt(i))|0; } return h; };
  const A=SNAPF("A");
  el("samplebar")._h["change"]({target:{dataset:{i:"0"}, checked:false, classList:{contains:()=>false}}}); __drain();
  for(const t of ["overview","gene","de","coexp","heat","enrich","compare","patterns","disco","traits"]){ try{ switchTab(t); __drain(); }catch(e){} }
  const B=SNAPF("B");
  const boxes={}; for(const id of Object.keys(elCache)){ const c=elCache[id]._c||""; if(/class="stalebar"/.test(c)) boxes[id]=1; }
  const rep=[]; for(const id of new Set(Object.keys(A).concat(Object.keys(B)))){ const a=(A[id]||[]).map(hash).join(","), b=(B[id]||[]).map(hash).join(","); const bt=(B[id]||[]).join("\n");
    rep.push([id, (A[id]||[]).length, (B[id]||[]).length, a===b, /^BLOCKED/.test(bt)? "blocked" : /^HIDDEN/.test(bt)? "hidden" : "", /THROW/.test(bt)? bt.slice(0,200) : "", /^HIDDEN/.test((A[id]||[""])[0])]); }
  print("STALE " + JSON.stringify({rep, boxes:Object.keys(boxes), incl0:included[0]}));
