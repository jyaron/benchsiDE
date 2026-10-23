// click every visible export button and capture the file it writes

  const FILES=[]; const _dc=downloadCSV, _dt=downloadText, _sb=saveBlob;
  let CUR="";
  downloadCSV=(t,f)=>{ FILES.push({src:CUR, f, kind:"csv", t:String(t)}); };
  downloadText=(t,f)=>{ FILES.push({src:CUR, f, kind:"text", t:String(t)}); };
  saveBlob=(b,f)=>{ const p=(b&&b.p)||[]; const txt=p.every(x=>typeof x==="string"); FILES.push({src:CUR, f, kind:txt? "text" : "bin", t: txt? p.join("") : "", bytes: !txt&&p[0]&&p[0].length<4e6? Array.from(p[0]) : null, nb: !txt&&p[0]? p[0].length : 0}); };
  const CLICKERR=[], NODL=[];
  const ids=Object.keys(elCache).filter(id=>(elCache[id]._h&&elCache[id]._h.click)||typeof elCache[id].onclick==="function");
  for(const id of ids){ const fn=(elCache[id]._h&&elCache[id]._h.click)||elCache[id].onclick; const src=String(fn);
    if(!/download|saveBlob|export|Export|xlsx|CSV|GMT|cyto/i.test(src)) continue;
    if(/^(enrrun|hubrun|dscrun|trrun|modfrun|savesess|reportbtn|sumfigbtn)$/.test(id)) continue;
    if(elCache[id].style&&elCache[id].style.display==="none") continue;
    CUR=id; const n0=FILES.length;
    try{ fn({preventDefault(){}, target:elCache[id], stopPropagation(){}}); __drain(); }catch(e){ CLICKERR.push(id+": "+String(e).slice(0,200)); }
    if(FILES.length===n0) NODL.push(id); }
  downloadCSV=_dc; downloadText=_dt; saveBlob=_sb;
  print("FILES " + JSON.stringify(FILES.map(x=>({src:x.src, f:x.f, kind:x.kind, n:x.t.length, t:x.t.length>1500000? x.t.slice(0,x.t.lastIndexOf("\n",1500000)) : x.t, trunc:x.t.length>1500000, bytes:x.bytes, nb:x.nb}))));
  print("CLICK " + JSON.stringify({ERR, CLICKERR, NODL, ids:ids.length}));
