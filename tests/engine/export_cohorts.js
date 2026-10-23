// cross-dataset analyses on the four psoriasis cohorts and the mouse model
const ERR=[];

  const ERRX=[]; const T_=(lab,f)=>{ try{ f(); __drain(); }catch(e){ ERRX.push(lab+": "+String(e).slice(0,200)); } };
  T_("fam", ()=>{ el("famx").value="s"; el("famds").value="3"; runFamilies(); });
  T_("cmp", ()=>runCompare()); T_("multi", ()=>runMulti());
  T_("sig", ()=>{ buildSigSelector(); el("sigsel").value="de_up"; el("sigfam").checked=true; runSigTransfer(); });
  T_("hubs", ()=>{ el("hubM").value="800"; el("hubrun")._h.click(); }); T_("nb", ()=>openNeighborhood(HUBS.rows[0].g));
  T_("hubx", ()=>{ runSync(hubxGen()); });
  print("ST " + JSON.stringify({fam:!!FAMROWS, famInfo:(el("faminfo")._c||"").replace(/<[^>]+>/g,"").slice(0,200), meta:!!META, sig:!!SIGRES, cmp:!!CMPLAST, multi:!!CONSENSUS, hubx:!!HUBX, ERRX}));
