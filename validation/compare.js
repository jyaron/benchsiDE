// benchsiDE external validation: compares the application's computations on the public
// demo dataset (GEO GSE63310) with the R reference values in gse63310_reference.json.
// Runs inside the application page. Used by tests/validation.spec.js and by VALIDATION.md.
//
// runExternalValidation(matrixText, designText, REF) -> {rows, fixture}
//   rows: [{name, n, delta, tol, pass}]  delta = maximum deviation over the n compared values
//   (absolute unless the name says "relative"; exact counts use tolerance 0).
function runExternalValidation(M, D, REF) {
  const rows = [];
  const ok = (name, n, delta, tol) => rows.push({ name, n, delta, tol, pass: delta <= tol });
  const absMax = (a, b) => a.reduce((m, v, i) => Math.max(m, Math.abs(v - b[i])), 0);
  const relMax = (a, b) => a.reduce((m, v, i) => Math.max(m, Math.abs(v - b[i]) / Math.max(Math.abs(b[i]), 1e-300)), 0);
  const G = REF.compared_genes;
  const at = (arr) => G.map(g => arr[g]);
  const hitsOf = (r) => { let u = 0, d = 0;
    for (let g = 0; g < nG; g++) if (r.q[g] <= 0.05 && Math.abs(r.lfc[g]) >= 1) { if (r.lfc[g] > 0) u++; else d++; }
    return [u, d]; };

  // design with the two test covariates, in matrix column order
  const cols = M.split("\n")[0].replace(/\r$/, "").split("\t").slice(1);
  const cell = {};
  D.trim().split("\n").slice(1).forEach(l => { const f = l.replace(/\r$/, "").split("\t"); cell[f[0]] = f[1]; });
  const Dcov = "sample\tcelltype\tlcov\tblk\tlane\n" +
    cols.map((s, q) => [s, cell[s], REF.covariates.lcov[q], REF.covariates.blk[q], REF.covariates.lane[q]].join("\t")).join("\n");

  el("normsel").value = "tmm"; el("fmode").value = "fbe";
  loadFromText(M, Dcov); GROUPTYPE = "categorical"; analyzeNow();
  const fixture = { factor: FACTOR, groups: GROUPS.join(","), nS, nG, covariates: Object.keys(COVARS).join(",") };

  // filtering and normalization
  const allIds = M.trim().split("\n").slice(1).map(l => l.split("\t")[0]);
  const pos = {}; allIds.forEach((id, i) => { pos[id] = i + 1; });
  const checksum = IDS.reduce((s, id) => s + pos[id], 0);
  ok("filterByExpr: kept-gene count and index checksum (exact)", 2,
    Math.abs(nG - REF.n_kept) + Math.abs(checksum - REF.kept_checksum), 0);
  ok("TMM normalization factors", nS, absMax(TMMF, REF.tmm), 1e-9);
  ok("log-CPM (edgeR, prior.count 2), gene 1, all samples", nS, absMax(Array.from({ length: nS }, (_, q) => L[q]), REF.logcpm_1), 1e-9);

  // differential expression, LP vs Basal
  el("selA").value = "Basal"; el("selB").value = "LP"; el("selFDR").value = "0.05"; el("fcthr").value = "1";
  const de = (method, cov) => { deMethod = method; DECOVSEL = cov; deKey = ""; return computeDE(); };
  const savFit = DEFIT; DEFIT = "pair";   // checks below up to the full-design block use the two-group fit option
  let r = de("mod", []);
  const tLP = Array.from(r.tv);
  ok("moderated t", G.length, absMax(at(r.tv), REF.mod.t), 1e-6);
  ok("moderated t p-values (relative)", G.length, relMax(at(r.p), REF.mod.p), 1e-6);
  ok("BH-adjusted p-values", G.length, absMax(at(r.q), REF.mod.q), 1e-8);
  ok("empirical-Bayes prior, trend (d0, median s0^2)", 2, Math.max(Math.abs(r.prior.df0 - REF.mod.df0), Math.abs(r.prior.s20 - REF.mod.s20)), 1e-6);
  const hm = hitsOf(r);
  ok("moderated t: significant genes up/down (exact)", 2, absMax(hm, REF.mod.hits), 0);
  // over-representation of the up hits in SET_TOP_LFC
  const upSet = new Set(); for (let g = 0; g < nG; g++) if (r.q[g] <= 0.05 && r.lfc[g] >= 1) upSet.add(g);
  const kO = REF.sets.SET_TOP_LFC.filter(g => upSet.has(g)).length;
  ok("hypergeometric ORA: overlap and hit counts (exact)", 2, Math.abs(kO - REF.ora.k) + Math.abs(upSet.size - REF.ora.q), 0);
  ok("hypergeometric ORA p (relative)", 1, relMax([hyperUpperP(kO, upSet.size, REF.ora.m, nG)], [REF.ora.p]), 1e-6);

  r = de("voom", []);
  ok("voom: moderated t", G.length, absMax(at(r.tv), REF.voom.t), 1e-6);
  ok("voom: prior d0", 1, Math.abs(r.prior.df0 - REF.voom.df0), 1e-6);
  ok("voom: significant genes up/down (exact)", 2, absMax(hitsOf(r), REF.voom.hits), 0);

  r = de("mod", ["lcov"]);
  ok("continuous covariate: moderated t", G.length, absMax(at(r.tv), REF.cov_continuous.t), 1e-6);
  ok("continuous covariate: prior d0", 1, Math.abs(r.prior.df0 - REF.cov_continuous.df0), 1e-6);
  ok("continuous covariate: significant genes (exact)", 2, absMax(hitsOf(r), REF.cov_continuous.hits), 0);
  r = de("mod", ["blk"]);
  ok("blocking factor: moderated t", G.length, absMax(at(r.tv), REF.cov_block.t), 1e-6);
  ok("blocking factor: prior d0", 1, Math.abs(r.prior.df0 - REF.cov_block.df0), 1e-6);
  ok("blocking factor: significant genes (exact)", 2, absMax(hitsOf(r), REF.cov_block.hits), 0);

  // moderated F over all three groups
  DECOVSEL = []; deMethod = "mod"; deKey = "";
  computeModF();
  ok("moderated F, three groups (relative)", G.length, relMax(at(MODF.F), REF.modF.F), 1e-7);
  ok("moderated F: prior d0", 1, Math.abs(MODF.df0 - REF.modF.df0), 1e-6);
  MODF = null;

  // display-only batch removal
  BATCHDISP = "blk"; LADJKEY = "";
  const LM = dispM(), mine = [], refv = [];
  for (let q = 0; q < nS; q++) for (let k = 0; k < 10; k++) { mine.push(LM[G[k] * nS + q]); refv.push(REF.rbe_rows[k + 10 * q]); }
  ok("removeBatchEffect (display values)", mine.length, absMax(mine, refv), 1e-9);
  BATCHDISP = ""; LADJKEY = "";

  // MDS (the sign of each dimension is arbitrary)
  const md = computeMDS();
  const sx = Math.sign(md.x[0]) === Math.sign(REF.mds.x[0]) ? 1 : -1, sy = Math.sign(md.y[0]) === Math.sign(REF.mds.y[0]) ? 1 : -1;
  ok("MDS coordinates and variance explained", 2 * nS + 2, Math.max(
    absMax(md.x.map(v => sx * v), REF.mds.x), absMax(md.y.map(v => sy * v), REF.mds.y), absMax(md.ve.slice(0, 2), REF.mds.ve)), 1e-8);

  // FRY, without and with the continuous covariate
  const setList = Object.keys(REF.sets).map(k => ({ name: k, rows: REF.sets[k] }));
  const fryCheck = (label, cov, ref) => {
    DECOVSEL = cov;
    const fr = fryTest(setList);
    if (fr.err) { ok(`FRY ${label}: runs (error: ${fr.err})`, 1, 1, 0); return; }
    ok(`FRY ${label}: p and mixed p (relative)`, 2 * fr.rows.length,
      Math.max(relMax(fr.rows.map(x => x.p), ref.p), relMax(fr.rows.map(x => x.pm), ref.pm)), 1e-5);
    ok(`FRY ${label}: direction (exact)`, fr.rows.length, fr.rows.filter((x, i) => x.dir !== ref.dir[i]).length, 0);
  };
  fryCheck("without covariate", [], REF.fry);
  fryCheck("with continuous covariate", ["lcov"], REF.fry_cov);
  DECOVSEL = [];

  // regression engine (co-expression tab)
  const xv = [], yv = [];
  for (let q = 0; q < nS; q++) { xv.push(L[q]); yv.push(L[nS + q]); }
  const f = olsCI(xv, yv);
  ok("OLS regression: intercept, slope, SE, p", 4, absMax([f.a, f.b, f.seB, f.pB], REF.ols), 1e-8);

  // co-expression table: r, slope and p of every compared gene against gene 1
  if (REF.coexp) {
    coGene = 0; coKey = ""; drawCoexp();
    const cg = REF.coexp.genes;
    ok("co-expression: Pearson r", cg.length, absMax(cg.map(g => coR[g]), REF.coexp.r), 1e-10);
    ok("co-expression: slope", cg.length, absMax(cg.map(g => coB[g]), REF.coexp.slope), 1e-10);
    ok("co-expression: p for r = 0 (relative)", cg.length, relMax(cg.map(g => coP[g]), REF.coexp.p), 1e-6);
    coGene = null; coKey = "";
  }
  // hub neighbourhood, within-group mode: partial correlation given cell type and its t-test p
  if (REF.pcor) {
    const inc = S.map((_, j) => j).filter(j => included[j]);
    const gl = [...new Set(inc.map(j => S[j].group))], gc = inc.map(j => gl.indexOf(S[j].group));
    const u1 = residNorm(inc.map(j => L[0 * nS + j]), gc, gl.length), df = inc.length - gl.length - 1;
    const rr = REF.pcor.genes.map(g => { const u = residNorm(inc.map(j => L[g * nS + j]), gc, gl.length); let s2 = 0; for (let q = 0; q < u.length; q++) s2 += u[q] * u1[q]; return s2; });
    ok("hub neighbourhood: partial correlation given group", rr.length, absMax(rr, REF.pcor.r), 1e-10);
    ok("hub neighbourhood: p for partial r = 0 (relative)", rr.length, relMax(rr.map(r => corrP(r, df)), REF.pcor.p), 1e-6);
  }
  // hub across datasets: soft connectivity, WGCNA module-preservation statistics, robust rank aggregation
  if (REF.kconn) {
    const inc = S.map((_, j) => j).filter(j => included[j]);
    const u = unitRowsGeneric((g, c) => L[g * nS + c], Array.from({ length: 200 }, (_, g) => g), inc, inc.map(() => 0), 1);
    ok("connectivity: sum of |r|^6", 200, relMax(Array.from(connectivityOver(u, 200, 6)), REF.kconn), 1e-9);
  }
  if (REF.pres) {
    const P = REF.pres, uni = (gs, cols) => unitRowsGeneric((g, c) => L[g * nS + c], gs, cols, cols.map(() => 0), 1);
    const ur = uni(P.mod, P.refS), Rr = presRefPart(corFromUnit(ur.U, ur.n, P.mod.map((_, i) => i)), P.mod.length);
    const names = ["propVarExplained", "meanSignAwareKME", "meanSignAwareCorDat", "meanAdj", "corkIM", "corkME", "corcor"];
    const st = gs => { const ut = uni(gs, P.tstS); const o = presStats(Rr, corFromUnit(ut.U, ut.n, gs.map((_, i) => i)), gs.length); return names.map(k => o[k]); };
    const js = st(P.mod).concat(...P.sets.map(st)), rf = P.obs.concat(...P.perm);
    ok("module preservation: 7 statistics, observed and 5 random sets (WGCNA " + P.wgcna + ")", js.length, absMax(js, rf), 1e-10);
  }
  if (REF.trait) {
    const T = REF.trait, saveC = COVARS;
    COVARS = Object.assign({}, COVARS, { TraitX: {}, BlockX: {} });
    S.forEach((_, j) => { COVARS.TraitX[j] = String(T.x[j]); COVARS.BlockX[j] = "b" + T.block[j]; });
    const R1 = runSync(traitDEGen({ trait: "TraitX", transform: "none", groups: [], adjustGroup: false, covs: [], block: null }));
    ok("trait association: moderated t (limma-trend, ~ trait)", 50, absMax(Array.from(R1.t.slice(0, 50)), T.t_all), 1e-8);
    ok("trait association: p (relative)", 50, relMax(Array.from(R1.p.slice(0, 50)), T.p_all), 1e-6);
    const R2 = runSync(traitDEGen({ trait: "TraitX", transform: "none", groups: [], adjustGroup: true, covs: [], block: "BlockX", blockMode: "dupcor" }));
    ok("trait association with blocks: duplicateCorrelation consensus", 1, Math.abs(R2.rho - T.rho), 1e-8);
    ok("trait association with blocks: moderated t (lmFit block, correlation)", 50, absMax(Array.from(R2.t.slice(0, 50)), T.t_blk), 1e-6);
    const sets = [Array.from({ length: 30 }, (_, i) => i), Array.from({ length: 31 }, (_, i) => 499 + 50 * i), Array.from({ length: 101 }, (_, i) => 2999 + i)];
    const cam = cameraOnTrait(R1, sets.map((rows, i) => ({ name: "s" + i, rows })));
    ok("trait gene sets: CAMERA p (relative)", 3, relMax(cam.map(c => c.p), T.cam_p), 1e-6);
    if (T.me) {
      const e = moduleEigengene(Array.from({ length: 40 }, (_, i) => i), R1.Y, R1.n), sc = T.me.reduce((a, v, i) => a + v * e.me[i], 0) / e.me.reduce((a, v) => a + v * v, 0);
      ok("module eigengene and its trait t (WGCNA moduleEigengenes, lm)", 10, Math.max(absMax(Array.from(e.me, v => v * sc), T.me), Math.abs(moduleTraitTest(R1, e.me).t - T.me_t)), 1e-8);
    }
    COVARS = saveC;
  }
  if (REF.rra) {
    const n = 300, ranks = REF.rra.lists.map(l => { const r = new Int32Array(n); l.forEach((it, p) => r[it] = p + 1); return r; });
    const sc = REF.rra.item.map(it => rraRho(ranks.map(r => r[it] / n)));
    ok("robust rank aggregation: RRA score (relative)", n, relMax(sc, REF.rra.score), 1e-9);
  }

  // density
  const v1 = Array.from({ length: nG }, (_, g) => L[g * nS]);
  const kd = kdeExact(v1, 512);
  ok("density: nrd0 bandwidth", 1, Math.abs(kd.bw - REF.kde.bw), 1e-12);
  ok("density: curve values at 9 grid points", 9, absMax(REF.kde.probes.map(j => kd.ys[j]), REF.kde.y), 1e-8);

  // RRHO: LP vs Basal against ML vs Basal
  el("selB").value = "ML"; deKey = "";
  const tML = Array.from(de("mod", []).tv);
  CMPLAST = { tA: tLP, tB: tML, N: nG };
  drawRRHO();
  const pr = el("p_rrho");
  const z = (pr && pr.data && pr.data[0] && pr.data[0].z) ||
    (typeof plotlyCalls !== "undefined" ? plotlyCalls.filter(c => c.id === "p_rrho").pop().traces[0].z : null);
  const cells = [], refc = [];
  const push = (mv, rv) => { if (mv >= 0 && isFinite(rv)) { cells.push(mv); refc.push(rv); } };
  z.forEach((row, i) => push(row[i], REF.rrho.diag[i]));
  z[0].forEach((mv, j) => push(mv, REF.rrho.row1[j]));
  ok("RRHO grid dimensions (exact)", 2, Math.abs(z.length - REF.rrho.dim[0]) + Math.abs(z[0].length - REF.rrho.dim[1]), 0);
  ok("RRHO -log10 p, diagonal and first row", cells.length, absMax(cells, refc), 1e-8);
  el("selB").value = "LP"; deKey = "";
 // full design (application default): every group fitted, LP vs Basal by contrast
  DEFIT = "all";
  r = de("mod", []);
  ok("full design: moderated t", G.length, absMax(at(r.tv), REF.full_mod.t), 1e-6);
  ok("full design: moderated t p-values (relative)", G.length, relMax(at(r.p), REF.full_mod.p), 1e-6);
  ok("full design: prior d0 and residual df", 2, Math.abs(r.prior.df0 - REF.full_mod.df0) + Math.abs(r.prior.dfr - REF.full_mod.dfr), 1e-6);
  ok("full design: significant genes up/down (exact)", 2, absMax(hitsOf(r), REF.full_mod.hits), 0);
  r = de("mod", ["lcov"]);
  ok("full design, continuous covariate: moderated t", G.length, absMax(at(r.tv), REF.full_covc.t), 1e-6);
  ok("full design, continuous covariate: significant genes (exact)", 2, absMax(hitsOf(r), REF.full_covc.hits), 0);
  r = de("mod", ["blk"]);
  ok("full design, blocking factor: moderated t", G.length, absMax(at(r.tv), REF.full_covb.t), 1e-6);
  ok("full design, blocking factor: significant genes (exact)", 2, absMax(hitsOf(r), REF.full_covb.hits), 0);
  r = de("voom", []);
  ok("full design, voom: moderated t", G.length, absMax(at(r.tv), REF.full_voom.t), 1e-6);
  ok("full design, voom: significant genes up/down (exact)", 2, absMax(hitsOf(r), REF.full_voom.hits), 0);
  r = de("voom", ["blk"]);
  ok("full design, voom with blocking factor (contrasts.fit): moderated t", G.length, absMax(at(r.tv), REF.full_voom_blk.t), 1e-6);
  ok("full design, voom with blocking factor: significant genes up/down (exact)", 2, absMax(hitsOf(r), REF.full_voom_blk.hits), 0);
  // TREAT (McCarthy & Smyth 2009), limma-trend, all groups
  const savTreat = DETREAT; DETREAT = true;
  r = de("mod", []);
  ok("TREAT, |log2FC| > 1, moderated t: t", G.length, absMax(at(r.tv), REF.treat_mod.t), 1e-6);
  ok("TREAT, |log2FC| > 1, moderated t: p (relative)", G.length, relMax(at(r.p), REF.treat_mod.p), 1e-6);
  { let u = 0, d = 0; for (let g = 0; g < nG; g++) if (r.q[g] <= 0.05) { if (r.lfc[g] > 0) u++; else d++; }
    ok("TREAT, moderated t: significant genes up/down (exact)", 2, absMax([u, d], REF.treat_mod.hits), 0); }
  // the analysis of the limma/Glimma/edgeR workflow article (Law et al. 2016): ~0 + group + lane, voom, TREAT lfc 1
  DETREAT = false; r = de("voom", ["lane"]);
  ok("workflow article, voom with lane: moderated t", G.length, absMax(at(r.tv), REF.article.t), 1e-6);
  ok("workflow article, voom with lane: prior d0", 1, Math.abs(r.prior.df0 - REF.article.df0), 1e-6);
  ok("workflow article, voom with lane: significant genes up/down (exact)", 2, absMax(hitsOf(r), REF.article.hits), 0);
  DETREAT = true; r = de("voom", ["lane"]);
  ok("workflow article, TREAT: t", G.length, absMax(at(r.tv), REF.article.treat_t), 1e-6);
  ok("workflow article, TREAT: p (relative)", G.length, relMax(at(r.p), REF.article.treat_p), 1e-6);
  { const sig = (A, B) => { el("selA").value = A; el("selB").value = B; deKey = ""; const rr = computeDE(); const s = new Set(); for (let g = 0; g < nG; g++) if (rr.q[g] <= 0.05) s.add(g); return s; };
    const sLP = sig("Basal", "LP"), sML = sig("Basal", "ML"); let both = 0; sLP.forEach(g => { if (sML.has(g)) both++; });
    ok("workflow article, TREAT: LP vs Basal, ML vs Basal and shared genes (exact)", 3, absMax([sLP.size, sML.size, both], REF.article.venn), 0);
    el("selA").value = "Basal"; el("selB").value = "LP"; }
  DETREAT = savTreat;
  deMethod = "mod"; DECOVSEL = []; deKey = "";
  fryCheck("full design", [], REF.full_fry);
  DEFIT = savFit;

  // experiment-planning power model
  ok("power model: n per group and power (relative)", 2, relMax(
    [rnapowerJS(20, null, null, 0.4, 0.4, 2, 0.05, 0.8), rnapowerJS(20, 3, 3, 0.4, 0.4, 2, 0.05, null)], REF.rnapower), 1e-8);

  return { rows, fixture };
}

 