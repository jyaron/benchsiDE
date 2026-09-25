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
  const Dcov = "sample\tcelltype\tlcov\tblk\n" +
    cols.map((s, q) => [s, cell[s], REF.covariates.lcov[q], REF.covariates.blk[q]].join("\t")).join("\n");

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
  ok("log2(CPM + 1), gene 1, all samples", nS, absMax(Array.from({ length: nS }, (_, q) => L[q]), REF.logcpm_1), 1e-9);

  // differential expression, LP vs Basal
  el("selA").value = "Basal"; el("selB").value = "LP"; el("selFDR").value = "0.05"; el("fcthr").value = "1";
  const de = (method, cov) => { deMethod = method; DECOVSEL = cov; deKey = ""; return computeDE(); };
  let r = de("mod", []);
  const tLP = Array.from(r.tv);
  ok("moderated t", G.length, absMax(at(r.tv), REF.mod.t), 1e-6);
  ok("moderated t p-values (relative)", G.length, relMax(at(r.p), REF.mod.p), 1e-6);
  ok("BH-adjusted p-values", G.length, absMax(at(r.q), REF.mod.q), 1e-8);
  ok("empirical-Bayes prior (d0, s0^2)", 2, Math.max(Math.abs(r.prior.df0 - REF.mod.df0), Math.abs(r.prior.s20 - REF.mod.s20)), 1e-6);
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

  // experiment-planning power model
  ok("power model: n per group and power (relative)", 2, relMax(
    [rnapowerJS(20, null, null, 0.4, 0.4, 2, 0.05, 0.8), rnapowerJS(20, 3, 3, 0.4, 0.4, 2, 0.05, null)], REF.rnapower), 1e-8);

  return { rows, fixture };
}
