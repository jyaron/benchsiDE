// analyses not run by audit_sweep.js, run before the exports are collected

  try{ el("heatord").value="kmeans"; switchTab("heat"); drawHeat(); }catch(e){ print("E heat "+e); }
  try{ for(const c of ["dscfam","dsccur","dscgmt"]) el(c).checked=true; switchTab("disco"); el("dscrun")._h.click(); __drain(); }catch(e){ print("E disco "+e); }
  try{ el("ixA1").value=GROUPS[0]; el("ixB1").value=GROUPS[1]; el("ixA2").value=GROUPS[2]; el("ixB2").value=GROUPS[3]; runInteraction(); }catch(e){ print("E ix "+e); }
  try{ switchTab("de"); DESEL=new Set([0,1,2,3,4,5,6,7,8,9]); }catch(e){ print("E desel "+e); }
  try{ switchTab("patterns"); drawPatterns(); showPatCluster(0); }catch(e){ print("E pat "+e); }
  try{ openNeighborhood(HUBS.rows[0].g); }catch(e){ print("E nb "+e); }
  print("ST " + JSON.stringify({heat:!!heatClusters, disco:!!DISCO, ix:!!IXRES, hubnb:!!HUBNB}));

  const ERRX=[]; const T_=(lab,f)=>{ try{ f(); for(let i=0;i<40;i++){ if(typeof drainMicrotasks==="function") drainMicrotasks(); __drain(); } }catch(e){ ERRX.push(lab+": "+String(e).slice(0,200)); } };
  T_("trait", ()=>{ switchTab("traits"); __drain(); el("trsel").value="score"; el("trrun")._h.click(); });
  T_("tmod", ()=>{ el("trmodrun")._h.click(); }); T_("tcam", ()=>{ el("trcamrun")._h.click(); });
  //print("TR " + JSON.stringify(TRAITRES? {n:TRAITRES.n, df:TRAITRES.df, keys:Object.keys(TRAITRES).slice(0,20)} : null)); print("ST " + JSON.stringify({tc:traitColumns(), trait:!!TRAITRES, tmod:!!TRAITMOD, tcam:!!TRAITCAM, ERRX, trinfo:(el("trinfo")._c||"").replace(/<[^>]+>/g,"").slice(0,200), tm:(el("trmodinfo")._c||"").replace(/<[^>]+>/g,"").slice(0,200), tcm:(el("trcaminfo")._c||"").replace(/<[^>]+>/g,"").slice(0,200)}));
