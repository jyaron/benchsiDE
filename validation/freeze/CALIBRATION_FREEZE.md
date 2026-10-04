# benchsiDE 1.0.0-rc.1 (2026-10-03): null calibration and binomial-thinning simulation (CALIBRATION_FREEZE)

Build under test: `benchside/index.html`, APP_VERSION 1.0.0-rc.1 (2026-10-03); sha256 of the last `<script>` block `cadbda80bd67db0e3cb15fcbcf138c23fb75cb24ff17cf6c28a79c68886c7d5b`. The file was not modified. All application results were produced by the application's own functions executed in macOS JavaScriptCore with the stub DOM and stub Plotly of `tests/engine/harness.py`; this is not a browser, so rendering and browser-specific behaviour are not covered. References: R 4.5.3 with edgeR 4.8.2, limma 3.66.0, metafor 5.0.1 and seqgendiff 1.2.4 (with sva 3.58.0, clue 0.3.68, irlba 2.3.7, matchingR 2.0.0, pdist 1.2.1); Python with numpy, pandas and scipy. Count matrices GSE54456, GSE121212 and GSE186063 matched the sha256 values in `tests/engine/datasets.json` (GSE54456 match, GSE121212 match, GSE186063 match).

Summary. 499 checks were run; 430 passed. Every count and p-value the app produced on the null splits agreed with the R references (32 parity checks in total), and every meta-analysis p-value agreed with metafor after resolving 2 genes in which metafor stopped at a local maximum of the restricted likelihood or did not converge. The 69 calibration checks that failed concern methods, not the implementation: gene-level tests at 3 per group (6 cells), ORA of a top-200 list (30 cells), FRY (1 cell), fixed-effect and z-test meta-analyses with 3 per group (8 cells), and z-based or fixed-effect meta-analyses when effects are heterogeneous (24 of 24 settings). The default meta-analysis (REML + Hartung-Knapp) stayed within the binomial bound in 12 of 12 null configurations and its observed FDR was not significantly above 0.05 in 20 of 20 simulated settings, but its power at meta FDR 0.05 was at most 0.0011 with 2 or 3 datasets.

## 1. Design

**Null pools.** Skin (as in VAL_CAL): GSE54456 normal skin (81 samples), GSE121212 healthy controls (37), GSE121212 psoriasis non-lesional (27), GSE186063 non-lesional (28). Non-skin (new, Section 4): GSE57945 non-IBD ileal biopsies (42). For each split, 2n samples were drawn without replacement from one pool (`numpy.random.default_rng([sum(map(ord, pool)), n, seed, 77])`, as in VAL_CAL, so seeds 0-39 reproduce the VAL_CAL splits), the first n assigned to group A, and only those samples were loaded as a session (filterByExpr, TMM and log-CPM computed on the split). n = 3, 5, 10; seeds 0-99 (R = 100 per pool and size; 1,500 sessions).

**Outputs per session (application defaults unless stated).** DE: moderated t (limma-trend), voom, Welch; genes at BH FDR <= 0.05 (no fold-change threshold), fraction p < 0.05 and p < 0.001. Gene sets, built-in Hallmark (50 sets) and Reactome (1,736 sets; about 1,630-1,680 testable at size 5-2000): ORA of a pasted list of the 200 genes with the smallest moderated-t p; ORA of the DE hits (FDR 0.05, |log2FC| 1; runs only when at least 3 hits); FRY (directional p); CAMERA through the Discovery screen with the library as the only source and size 5-2000; GSEA with sample permutation, run at 10 v 10 only (at 3 and 5 per group fewer than 1,000 distinct labelings exist and this build switches to gene permutation, which reports no FDR). Discovery screen with its defaults (gene families + HGNC groups, size 3-100).

**Meta-analysis null.** K = 2, 3, 4 datasets were drawn from K different skin pools chosen at random per replicate; K = 5 used the four skin pools plus a second, disjoint GSE54456 split (only GSE54456 is large enough to supply two datasets). Each dataset had n v n samples (n = 3, 5, 10), was analysed with the default moderated t, and was added as a comparison dataset. `runMeta` was run with minimum datasets K (genes measured in all), no fold-change threshold, prediction-interval criterion off, for six variants: REML + Hartung-Knapp (default), REML + truncated Hartung-Knapp (`adhoc`), REML + z, DerSimonian-Laird + z (previous default), fixed effect + z, and DL + Hartung-Knapp. R = 120 per K and n (1,440 replicates); seeds `default_rng([K, n, seed, 4242])`. No dataset was flagged as a duplicate.

**Criteria.** A cell fails when the number of replicates with at least one discovery at FDR <= 0.05 exceeds the 99th percentile of Binomial(R, 0.05) (11 for R = 100, 12 for R = 120). Mean fraction p < 0.05 fails when mean - 3 SE > 0.05. Simulation FDR fails when mean FDP - 2 SE > 0.05. Parity: tolerance 0 for counts, 1e-8 absolute for fractions, 1e-6 relative for p-values. Across the 232 binomial cells, 1.20 exceedances are expected by chance if every method were exactly calibrated at 5%.

**Binomial thinning (Gerard 2020).** `seqgendiff::thin_diff(mat, design_fixed = group indicator (B = 1), coef_fixed = log2 effects, relative = TRUE)` on the raw counts of real null splits; genes with zero effect are left unchanged. Gene level: in each split, 10% of the genes that pass filterByExpr on that split received log2FC ~ N(0, 0.8^2). Gene sets: per replicate, 2 Hallmark and 8 Reactome sets (20-200 testable members, from the 10 v 10 seed-0 split of the same pool) were chosen at random with a random direction; half of each set's expressed members received log2FC = direction x |N(0, 0.8^2)|; all other genes were null. Meta-analysis: K datasets from different skin pools as above; 10% of the genes passing filterByExpr in all K splits received one shared log2FC ~ N(0, 0.8^2) applied to every dataset (homogeneous effects). A second meta-analysis scenario added heterogeneity: each non-null gene received its shared effect plus a dataset-specific N(0, 0.4^2) deviation, and a further 10% of genes received dataset-specific N(0, 0.4^2) effects with mean zero (null for the pooled effect). Pools: GSE54456 normal (skin) and GSE57945 non-IBD (non-skin) for the gene and set level; skin pools for the meta-analysis. n = 3, 5, 10 and R = 30 per setting (heterogeneous scenario: n = 5, 10 and R = 20). FDR = mean false discovery proportion V/max(R,1); power = true positives / non-null genes tested. The inserted effects are on the log2 scale: the slope of the observed change in log2 mean ratio on the inserted coefficient was 0.9889-1.0090 over all 2180 thinned datasets.

## 2. Parity with R on the same splits

| Check | Reference | n compared | Max deviation | Tolerance | Pass |
|---|---|---|---|---|---|
| PAR-01: Null splits: genes kept by filterByExpr | edgeR::filterByExpr | 1500 | 0 | 0 | yes |
| PAR-02: Null splits: genes at BH FDR<=0.05, mod | limma lmFit + eBayes(trend=TRUE) | 1500 | 0 | 0 | yes |
| PAR-03: Null splits: fraction of genes with p<0.05, mod | limma lmFit + eBayes(trend=TRUE) | 1500 | 5e-16 | 1e-08 | yes |
| PAR-04: Null splits: genes at BH FDR<=0.05, voom | limma voom + eBayes | 1500 | 0 | 0 | yes |
| PAR-05: Null splits: fraction of genes with p<0.05, voom | limma voom + eBayes | 1500 | 5e-16 | 1e-08 | yes |
| PAR-06: Null splits: genes at BH FDR<=0.05, welch | stats::t.test (Welch) | 1500 | 0 | 0 | yes |
| PAR-07: Null splits: fraction of genes with p<0.05, welch | stats::t.test (Welch) | 1500 | 5e-16 | 1e-08 | yes |
| PAR-08: Null splits: per-gene p-values (relative), mod | limma-trend | 75 | 1.64e-09 | 1e-06 | yes |
| PAR-09: Null splits: per-gene p-values (relative), voom | limma voom | 75 | 2.12e-10 | 1e-06 | yes |
| PAR-10: Null splits: per-gene p-values (relative), welch | stats::t.test | 75 | 2.92e-10 | 1e-06 | yes |
| PAR-11: Null splits: FRY Hallmark sets at FDR<=0.05 and at p<=0.05 | limma::fry | 3750 | 0 | 0 | yes |
| PAR-12: Null splits: FRY Hallmark set p-values (relative) | limma::fry | 3750 | 3.68e-11 | 1e-06 | yes |
| PAR-13: Null splits: FRY Reactome sets at FDR<=0.05 and at p<=0.05 | limma::fry | 124030 | 0 | 0 | yes |
| PAR-14: Null splits: FRY Reactome set p-values (relative) | limma::fry | 124030 | 1.77e-10 | 1e-06 | yes |
| PAR-15: Null splits: CAMERA Hallmark sets at FDR<=0.05 and at p<=0.05 | limma::camera(inter.gene.cor=NA) | 3750 | 0 | 0 | yes |
| PAR-16: Null splits: CAMERA Hallmark set p-values (relative) | limma::camera(inter.gene.cor=NA) | 3750 | 6.13e-12 | 1e-06 | yes |
| PAR-17: Null splits: CAMERA Reactome sets at FDR<=0.05 and at p<=0.05 | limma::camera(inter.gene.cor=NA) | 124030 | 0 | 0 | yes |
| PAR-18: Null splits: CAMERA Reactome set p-values (relative) | limma::camera(inter.gene.cor=NA) | 124030 | 7.84e-11 | 1e-06 | yes |
| PAR-19: Regression: splits with FDR hits and total hits, seeds 0-39, same split definition as VAL_CAL 4.1 | VAL_CAL.md table 4.1 (build 0.22.0-beta) | 36 | 0 | 0 | yes |
| PAR-20: Meta null dumps: per-gene p-value (relative), REML + Hartung-Knapp (default) | metafor::rma (REML_hk) | 18547 | 3.15e-09 | 1e-06 | yes |
| PAR-21: Meta null dumps: per-gene p-value (relative), REML + truncated HK | metafor::rma (REML_adhoc) | 18550 | 3.15e-09 | 1e-06 | yes |
| PAR-22: Meta null dumps: per-gene p-value (relative), REML + z | metafor::rma (REML_none) | 18552 | 1.01e-08 | 1e-06 | yes |
| PAR-23: Meta null dumps: per-gene p-value (relative), DL + z | metafor::rma (REM_none) | 18550 | 2.26e-13 | 1e-06 | yes |
| PAR-24: Meta null dumps: per-gene p-value (relative), Fixed effect + z | metafor::rma (FEM_none) | 18550 | 2.01e-13 | 1e-06 | yes |
| PAR-25: Meta null dumps: per-gene p-value (relative), DL + Hartung-Knapp | metafor::rma (REM_hk) | 18547 | 9.12e-13 | 1e-06 | yes |
| PAR-26: Four psoriasis cohorts (all and leave-one-out): genes at meta FDR<=0.05, shared-signature size and canonical genes recovered, REML + Hartung-Knapp (default) | metafor::rma (REML_hk) + BH, same yi/vi | 15 | 0 | 0 | yes |
| PAR-27: Four psoriasis cohorts (all and leave-one-out): genes at meta FDR<=0.05, shared-signature size and canonical genes recovered, REML + truncated HK | metafor::rma (REML_adhoc) + BH, same yi/vi | 15 | 0 | 0 | yes |
| PAR-28: Four psoriasis cohorts (all and leave-one-out): genes at meta FDR<=0.05, shared-signature size and canonical genes recovered, REML + z | metafor::rma (REML_none) + BH, same yi/vi | 15 | 0 | 0 | yes |
| PAR-29: Four psoriasis cohorts (all and leave-one-out): genes at meta FDR<=0.05, shared-signature size and canonical genes recovered, DL + z | metafor::rma (REM_none) + BH, same yi/vi | 15 | 0 | 0 | yes |
| PAR-30: Four psoriasis cohorts (all and leave-one-out): genes at meta FDR<=0.05, shared-signature size and canonical genes recovered, Fixed effect + z | metafor::rma (FEM_none) + BH, same yi/vi | 15 | 0 | 0 | yes |
| PAR-31: Four psoriasis cohorts (all and leave-one-out): genes at meta FDR<=0.05, shared-signature size and canonical genes recovered, DL + Hartung-Knapp | metafor::rma (REM_hk) + BH, same yi/vi | 15 | 0 | 0 | yes |
| PAR-32: Four cohorts, REML: genes whose p differs from metafor by >1e-6 are refitted at the global maximum of the restricted likelihood (profiled in R); app p vs p at that tau2 | metafor::rma with tau2 fixed at the profiled global REML maximum | 4 | 4.51e-08 | 1e-06 | yes |
| THN-01: Thinned matrices: counts of genes with zero effect unchanged | original GEO counts | 2180 | 0 | 0 | yes |
| THN-02: Thinned matrices: slope of observed change in log2 mean ratio on inserted log2FC (genes with mean count>=50) | inserted coefficients (seqgendiff coef_fixed) | 2180 | 0.0111 | 0.02 | yes |
| WRN-01: DE small-sample warning shown when the smaller group has <=3 samples (moderated t, voom), not otherwise | CHANGELOG 0.23.0 specification | 12 | 0 | 0 | yes |
| WRN-02: Warning text range (10-25% of null 3v3 splits of healthy skin with an FDR hit) covers the frozen-build moderated-t fractions | observed fractions, R=100 per pool | 4 | 0 | 0 | yes |

The REML parity on the four psoriasis cohorts flagged two genes. For 100507412 (cohorts without GSE186063) the restricted likelihood is bimodal; metafor converged to tau2 = 0.3220 (log-likelihood -1.354356), the app to the global maximum tau2 = 0.02798 (log-likelihood -1.350050); at that tau2 metafor gives the app's p-values (relative difference 3.0e-08). For TMEM126B (all four cohorts) metafor's default Fisher scoring did not converge; at the profiled maximum the p-values agree to 4.5e-08. All other genes agree within 1e-6, and every count of genes at FDR and of signature genes is identical.

## 3. Null calibration on the frozen build
### 3.1 Gene level (R = 100 per cell)

Entries: splits with at least one gene at FDR <= 0.05 / R (bold above the bound of 11); total genes in parentheses; mean fraction of genes with p < 0.05 (nominal 0.05).

| Pool | n | moderated t | voom | Welch | frac p<0.05 mod / voom / Welch | frac p<0.001 mod |
|---|---|---|---|---|---|---|
| GSE54456_normal | 3 | **21/100** (4408) | **18/100** (4162) | 2/100 (4) | 0.0680 / 0.0672 / 0.0442 | 0.00222 |
| GSE54456_normal | 5 | 3/100 (3185) | 3/100 (3001) | 4/100 (5) | 0.0551 / 0.0546 / 0.0470 | 0.00143 |
| GSE54456_normal | 10 | 2/100 (2274) | 2/100 (2333) | 1/100 (1) | 0.0544 / 0.0542 / 0.0517 | 0.00096 |
| GSE121212_CTRL | 3 | **17/100** (116) | 9/100 (45) | 1/100 (1) | 0.0473 / 0.0471 / 0.0330 | 0.00094 |
| GSE121212_CTRL | 5 | 4/100 (94) | 3/100 (8) | 1/100 (1) | 0.0453 / 0.0451 / 0.0394 | 0.00077 |
| GSE121212_CTRL | 10 | 4/100 (6) | 3/100 (5) | 2/100 (2) | 0.0517 / 0.0524 / 0.0496 | 0.00085 |
| GSE121212_PSOnl | 3 | **21/100** (279) | **12/100** (214) | 1/100 (1) | 0.0477 / 0.0488 / 0.0342 | 0.00118 |
| GSE121212_PSOnl | 5 | 5/100 (29) | 5/100 (15) | 4/100 (4) | 0.0489 / 0.0500 / 0.0428 | 0.00109 |
| GSE121212_PSOnl | 10 | 7/100 (743) | 7/100 (715) | 7/100 (428) | 0.0547 / 0.0567 / 0.0529 | 0.00143 |
| GSE186063_nonles | 3 | **14/100** (757) | 10/100 (69) | 2/100 (2) | 0.0491 / 0.0496 / 0.0333 | 0.00148 |
| GSE186063_nonles | 5 | 5/100 (43) | 3/100 (4) | 2/100 (2) | 0.0466 / 0.0464 / 0.0399 | 0.00086 |
| GSE186063_nonles | 10 | 4/100 (5) | 4/100 (4) | 1/100 (1) | 0.0546 / 0.0548 / 0.0524 | 0.00100 |
| GSE57945_nonIBD | 3 | 8/100 (2869) | 8/100 (2727) | 1/100 (1) | 0.0502 / 0.0497 / 0.0299 | 0.00142 |
| GSE57945_nonIBD | 5 | 4/100 (5168) | 4/100 (4709) | 2/100 (83) | 0.0478 / 0.0477 / 0.0392 | 0.00137 |
| GSE57945_nonIBD | 10 | 3/100 (6242) | 4/100 (6089) | 3/100 (4573) | 0.0494 / 0.0481 / 0.0466 | 0.00178 |

Seeds 0-39 reproduce VAL_CAL table 4.1 exactly (PAR-19). At 3 per group the moderated t exceeded the bound in 4 of 4 skin pools (21, 17, 21, 14 of 100 splits) and voom in 2 (18, 9, 12, 10); the non-skin pool stayed within it (moderated t 8/100). From 5 per group upward every gene-level cell is within the bound (maximum 7/100). No cell failed the mean-fraction criterion (largest mean fraction p < 0.05: 0.0680). The DE tab warns when the smaller group has 3 or fewer samples (WRN-01); its quoted range, 10-25% of null 3 v 3 skin splits, covers the frozen-build moderated-t values (WRN-02). Some splits gave thousands of genes (for example GSE57945 at 10 v 10: 3 of 100 splits with hits but 6242 genes in total); these are splits whose two groups differ systematically for reasons not recorded in GEO, and they affect all three tests.

### 3.2 Gene sets (R = 100 per cell)

Entries: splits with at least one set at FDR <= 0.05 / R (bold above the bound); H = Hallmark, R = Reactome. GSEA uses sample permutation (10 v 10 only).

| Pool | n | ORA top-200 H | ORA top-200 R | FRY H | FRY R | CAMERA H | CAMERA R | GSEA H | GSEA R |
|---|---|---|---|---|---|---|---|---|---|
| GSE54456_normal | 3 | **42/100** | **62/100** | 10/100 | **12/100** | 1/100 | 0/100 | - | - |
| GSE54456_normal | 5 | **40/100** | **56/100** | 9/100 | 5/100 | 2/100 | 0/100 | - | - |
| GSE54456_normal | 10 | **36/100** | **45/100** | 3/100 | 3/100 | 2/100 | 0/100 | 9/100 | 7/100 |
| GSE121212_CTRL | 3 | **45/100** | **73/100** | 4/100 | 4/100 | 0/100 | 0/100 | - | - |
| GSE121212_CTRL | 5 | **31/100** | **58/100** | 3/100 | 3/100 | 0/100 | 0/100 | - | - |
| GSE121212_CTRL | 10 | **25/100** | **40/100** | 3/100 | 5/100 | 0/100 | 0/100 | 8/100 | 9/100 |
| GSE121212_PSOnl | 3 | **52/100** | **64/100** | 7/100 | 8/100 | 2/100 | 0/100 | - | - |
| GSE121212_PSOnl | 5 | **40/100** | **54/100** | 6/100 | 3/100 | 1/100 | 1/100 | - | - |
| GSE121212_PSOnl | 10 | **33/100** | **39/100** | 9/100 | 10/100 | 2/100 | 0/100 | 7/100 | 8/100 |
| GSE186063_nonles | 3 | **44/100** | **73/100** | 4/100 | 3/100 | 0/100 | 0/100 | - | - |
| GSE186063_nonles | 5 | **36/100** | **55/100** | 3/100 | 3/100 | 1/100 | 0/100 | - | - |
| GSE186063_nonles | 10 | **27/100** | **39/100** | 7/100 | 6/100 | 2/100 | 0/100 | 11/100 | 10/100 |
| GSE57945_nonIBD | 3 | **26/100** | **58/100** | 6/100 | 7/100 | 0/100 | 0/100 | - | - |
| GSE57945_nonIBD | 5 | **39/100** | **55/100** | 2/100 | 7/100 | 1/100 | 0/100 | - | - |
| GSE57945_nonIBD | 10 | **25/100** | **34/100** | 3/100 | 2/100 | 0/100 | 0/100 | 4/100 | 6/100 |

FRY: 29 of 30 cells within the bound (maximum 12/100; mean fraction p < 0.05 0.041-0.076). CAMERA: maximum 2/100 (mean fraction p < 0.05 0.013-0.052). Sample-permutation GSEA at 10 v 10: 4-11/100 splits with a set at FDR <= 0.05, all within the bound; at the FDR <= 0.25 threshold the app shades, 23-36/100 splits had at least one set (not a 5% criterion). ORA of a top-200 list exceeded the bound in all 30 cells (25-73/100), as described in VAL_CAL; the app's ORA note now states that the test is anti-conservative for co-regulated lists. ORA of the DE hits could be run only in the 29 of 1,500 sessions with at least 3 DE hits, which are themselves false positives; in those runs 43 of 85 (Hallmark and Reactome together) gave sets at FDR <= 0.05 (descriptive only; the runs are conditional on false DE hits).

Discovery screen (families + HGNC groups): modules at FDR <= 0.05 in 0 of 1,500 null sessions; the mean number of modules at p <= 0.01 was 0.39-9.27 against 27.0-28.7 expected by chance (the screen is conservative, as in VAL_CAL).

### 3.3 Meta-analysis null (R = 120 per cell; bound 12)

| K | n | REML + Hartung-Knapp (default) | REML + truncated HK | REML + z | DL + z | Fixed effect + z | DL + Hartung-Knapp |
|---|---|---|---|---|---|---|---|
| 2 | 3 | 5/120 (0.049) | 0/120 (0.000) | 9/120 (0.035) | 9/120 (0.035) | **18/120** (0.047) | 5/120 (0.049) |
| 2 | 5 | 1/120 (0.051) | 0/120 (0.000) | 5/120 (0.041) | 5/120 (0.041) | 10/120 (0.054) | 1/120 (0.051) |
| 2 | 10 | 9/120 (0.051) | 0/120 (0.000) | 7/120 (0.042) | 7/120 (0.042) | 9/120 (0.056) | 9/120 (0.051) |
| 3 | 3 | 5/120 (0.051) | 0/120 (0.000) | 11/120 (0.035) | 10/120 (0.035) | **24/120** (0.047) | 5/120 (0.051) |
| 3 | 5 | 5/120 (0.051) | 0/120 (0.000) | 5/120 (0.036) | 5/120 (0.036) | 10/120 (0.048) | 5/120 (0.051) |
| 3 | 10 | 6/120 (0.053) | 0/120 (0.000) | 2/120 (0.040) | 1/120 (0.040) | 3/120 (0.054) | 6/120 (0.053) |
| 4 | 3 | 11/120 (0.053) | 0/120 (0.001) | **13/120** (0.038) | **13/120** (0.037) | **25/120** (0.050) | 11/120 (0.053) |
| 4 | 5 | 7/120 (0.052) | 0/120 (0.001) | 6/120 (0.040) | 5/120 (0.040) | 9/120 (0.054) | 7/120 (0.052) |
| 4 | 10 | 5/120 (0.051) | 0/120 (0.001) | 1/120 (0.038) | 1/120 (0.037) | 4/120 (0.051) | 5/120 (0.051) |
| 5 | 3 | 5/120 (0.052) | 0/120 (0.004) | **21/120** (0.039) | **22/120** (0.038) | **34/120** (0.053) | 5/120 (0.052) |
| 5 | 5 | 11/120 (0.055) | 0/120 (0.004) | 8/120 (0.040) | 8/120 (0.040) | 10/120 (0.053) | 11/120 (0.055) |
| 5 | 10 | 10/120 (0.047) | 0/120 (0.003) | 2/120 (0.035) | 3/120 (0.035) | 4/120 (0.047) | 10/120 (0.047) |

Entries: replicates with at least one gene at meta FDR <= 0.05 / R, and mean fraction of genes with p < 0.05 in parentheses. REML + Hartung-Knapp stayed within the bound in all 12 configurations (maximum 11/120; mean fraction p < 0.05 0.047-0.055). REML + z and DL + z exceeded it with 4 and 5 datasets of 3 v 3 (REML + z 13 and 21/120); the fixed-effect model exceeded it in every configuration with 3 per group (18-34/120). The truncated Hartung-Knapp test gave no discovery in any of the 1,440 replicates (mean fraction p < 0.05 at most 0.0044).

## 4. Non-skin null pool: GSE57945

GSE57945 ("Core Ileal Transcriptome in Pediatric Crohn Disease"; RISK cohort) has 322 samples in its series matrix, all `tissue: Ileal biopsy`, all Illumina HiSeq 2000, from treatment-naive children and adolescents undergoing diagnostic colonoscopy; diagnosis: cd 218, not ibd 42, uc 62. The series design states that non-IBD controls had suspected IBD but no microscopic or macroscopic inflammation and normal radiographic, endoscopic and histologic findings. The 42 samples annotated `diagnosis: Not IBD` (or `not IBD`) form the pool: one tissue, one condition, one study, one platform; sex Male 26, Female 16; age at diagnosis 2.3-16.9 years; histopathology field NA 42. They are symptomatic non-IBD controls rather than healthy volunteers, which is the usual control group of a diagnostic biopsy study. NCBI-generated counts: `GSE57945_raw_counts_GRCh38.p13_NCBI.tsv.gz` (sha256 of the gzip file `a88a3aece0c9435ac93d6e4e130729b3356a4af78cf85dd719973eb9b2271126`, of the decompressed table `48b89ba9e8b6a6fa948bec306f415cfee77a498ab51539de789fe336957a3e9d`; series matrix sha256 `8fac343c1554d4f49c7fbbc329b38158ff76671bfffa5d1846520487c7214804`); all 42 control samples are in the count table.

Gene-level and gene-set results for this pool are in the tables of Section 3 (last three rows of each). All gene-level cells were within the bound, including 3 v 3 (moderated t 8/100, voom 8/100, Welch 1/100); FRY, CAMERA and sample-permutation GSEA were within the bound; ORA of the top-200 list exceeded it as in skin.

## 5. Binomial-thinning simulation
### 5.1 Gene level (R = 30)

| Pool | n | Method | Observed FDR (SE) | Power | Power, abs(true log2FC) >= 1 | Mean calls | Non-null genes tested |
|---|---|---|---|---|---|---|---|
| GSE54456_normal | 3 | mod | 0.0791 (0.0218) | 0.070 | 0.281 | 149 | 1905 |
| GSE54456_normal | 3 | voom | 0.0704 (0.0208) | 0.065 | 0.263 | 136 | 1905 |
| GSE54456_normal | 3 | welch | 0.0333 (0.0333) | 0.000 | 0.001 | 0 | 1905 |
| GSE54456_normal | 5 | mod | 0.0479 (0.0252) | 0.243 | 0.690 | 564 | 1885 |
| GSE54456_normal | 5 | voom | 0.0472 (0.0250) | 0.242 | 0.687 | 557 | 1885 |
| GSE54456_normal | 5 | welch | 0.0367 (0.0212) | 0.150 | 0.467 | 322 | 1885 |
| GSE54456_normal | 10 | mod | 0.0488 (0.0242) | 0.439 | 0.914 | 930 | 1869 |
| GSE54456_normal | 10 | voom | 0.0479 (0.0239) | 0.439 | 0.912 | 927 | 1869 |
| GSE54456_normal | 10 | welch | 0.0458 (0.0224) | 0.418 | 0.892 | 868 | 1869 |
| GSE57945_nonIBD | 3 | mod | 0.0440 (0.0209) | 0.020 | 0.077 | 40 | 1724 |
| GSE57945_nonIBD | 3 | voom | 0.0217 (0.0129) | 0.018 | 0.068 | 33 | 1724 |
| GSE57945_nonIBD | 3 | welch | 0.0000 (0.0000) | 0.000 | 0.000 | 0 | 1724 |
| GSE57945_nonIBD | 5 | mod | 0.0520 (0.0272) | 0.126 | 0.431 | 290 | 1713 |
| GSE57945_nonIBD | 5 | voom | 0.0619 (0.0285) | 0.132 | 0.449 | 306 | 1713 |
| GSE57945_nonIBD | 5 | welch | 0.0367 (0.0198) | 0.062 | 0.216 | 126 | 1713 |
| GSE57945_nonIBD | 10 | mod | 0.0436 (0.0281) | 0.306 | 0.803 | 623 | 1698 |
| GSE57945_nonIBD | 10 | voom | 0.0456 (0.0288) | 0.321 | 0.825 | 660 | 1698 |
| GSE57945_nonIBD | 10 | welch | 0.0396 (0.0261) | 0.284 | 0.750 | 566 | 1698 |

All gene-level settings controlled the FDR (largest mean FDP 0.0791, GSE54456_normal mod 3v3). The pooled ratio of false to all calls over replicates is larger (up to 0.261) because a few splits with systematic group differences contribute many false calls; the FDR is the mean of the per-replicate proportion. Welch had essentially no power at 3 per group.

### 5.2 Gene sets (R = 30; 2 Hallmark and 8 Reactome sets with inserted effects per replicate)

| Pool | n | Library | Method | Touched sets called | Clean null sets called | Untouched sets overlapping non-null genes called | Called sets that are clean null (SE) | Replicates |
|---|---|---|---|---|---|---|---|---|
| GSE54456_normal | 3 | Hallmark | cam | 0.067 | 0.0000 | 0.003 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 3 | Hallmark | fry | 0.183 | 0.0000 | 0.007 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 3 | Hallmark | ora200 | 0.917 | 0.0000 | 0.028 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 3 | Hallmark | oraDE | 0.393 | 0.0000 | 0.006 | 0.0000 (0.0000) | 14 |
| GSE54456_normal | 3 | Reactome | cam | 0.000 | 0.0000 | 0.000 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 3 | Reactome | fry | 0.071 | 0.0009 | 0.010 | 0.0079 (0.0037) | 30 |
| GSE54456_normal | 3 | Reactome | ora200 | 0.783 | 0.0026 | 0.166 | 0.0121 (0.0034) | 30 |
| GSE54456_normal | 3 | Reactome | oraDE | 0.289 | 0.0001 | 0.064 | 0.0003 (0.0003) | 16 |
| GSE54456_normal | 5 | Hallmark | cam | 0.483 | 0.0000 | 0.002 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 5 | Hallmark | fry | 0.683 | 0.0111 | 0.013 | 0.0111 (0.0111) | 30 |
| GSE54456_normal | 5 | Hallmark | ora200 | 1.000 | 0.0000 | 0.029 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 5 | Hallmark | oraDE | 0.933 | 0.0000 | 0.014 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 5 | Reactome | cam | 0.142 | 0.0001 | 0.013 | 0.0002 (0.0002) | 30 |
| GSE54456_normal | 5 | Reactome | fry | 0.571 | 0.0049 | 0.085 | 0.0241 (0.0050) | 30 |
| GSE54456_normal | 5 | Reactome | ora200 | 0.929 | 0.0007 | 0.226 | 0.0029 (0.0012) | 30 |
| GSE54456_normal | 5 | Reactome | oraDE | 0.667 | 0.0002 | 0.119 | 0.0011 (0.0011) | 30 |
| GSE54456_normal | 10 | Hallmark | cam | 0.750 | 0.0000 | 0.004 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 10 | Hallmark | fry | 0.950 | 0.0000 | 0.029 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 10 | Hallmark | ora200 | 1.000 | 0.0000 | 0.047 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 10 | Hallmark | oraDE | 0.883 | 0.0000 | 0.022 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 10 | Reactome | cam | 0.675 | 0.0002 | 0.058 | 0.0005 (0.0004) | 30 |
| GSE54456_normal | 10 | Reactome | fry | 0.892 | 0.0054 | 0.205 | 0.0145 (0.0051) | 30 |
| GSE54456_normal | 10 | Reactome | ora200 | 0.971 | 0.0000 | 0.279 | 0.0000 (0.0000) | 30 |
| GSE54456_normal | 10 | Reactome | oraDE | 0.779 | 0.0000 | 0.137 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 3 | Hallmark | cam | 0.033 | 0.0000 | 0.000 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 3 | Hallmark | fry | 0.083 | 0.0000 | 0.008 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 3 | Hallmark | ora200 | 0.817 | 0.0000 | 0.024 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 3 | Hallmark | oraDE | 0.500 | 0.0000 | 0.156 | 0.0000 (0.0000) | 1 |
| GSE57945_nonIBD | 3 | Reactome | cam | 0.000 | 0.0000 | 0.000 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 3 | Reactome | fry | 0.013 | 0.0043 | 0.005 | 0.0436 (0.0345) | 30 |
| GSE57945_nonIBD | 3 | Reactome | ora200 | 0.746 | 0.0020 | 0.138 | 0.0130 (0.0053) | 30 |
| GSE57945_nonIBD | 3 | Reactome | oraDE | 0.188 | 0.0084 | 0.032 | 0.0617 (0.0617) | 2 |
| GSE57945_nonIBD | 5 | Hallmark | cam | 0.133 | 0.0000 | 0.000 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 5 | Hallmark | fry | 0.250 | 0.0000 | 0.010 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 5 | Hallmark | ora200 | 0.933 | 0.0000 | 0.023 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 5 | Hallmark | oraDE | 0.540 | 0.0000 | 0.006 | 0.0000 (0.0000) | 25 |
| GSE57945_nonIBD | 5 | Reactome | cam | 0.071 | 0.0000 | 0.003 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 5 | Reactome | fry | 0.283 | 0.0048 | 0.037 | 0.0160 (0.0068) | 30 |
| GSE57945_nonIBD | 5 | Reactome | ora200 | 0.875 | 0.0009 | 0.198 | 0.0031 (0.0017) | 30 |
| GSE57945_nonIBD | 5 | Reactome | oraDE | 0.404 | 0.0001 | 0.059 | 0.0003 (0.0003) | 26 |
| GSE57945_nonIBD | 10 | Hallmark | cam | 0.617 | 0.0000 | 0.002 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 10 | Hallmark | fry | 0.767 | 0.0000 | 0.021 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 10 | Hallmark | ora200 | 0.983 | 0.0000 | 0.051 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 10 | Hallmark | oraDE | 0.900 | 0.0000 | 0.014 | 0.0000 (0.0000) | 30 |
| GSE57945_nonIBD | 10 | Reactome | cam | 0.350 | 0.0001 | 0.032 | 0.0004 (0.0004) | 30 |
| GSE57945_nonIBD | 10 | Reactome | fry | 0.637 | 0.0074 | 0.104 | 0.0227 (0.0079) | 30 |
| GSE57945_nonIBD | 10 | Reactome | ora200 | 0.967 | 0.0003 | 0.226 | 0.0009 (0.0008) | 30 |
| GSE57945_nonIBD | 10 | Reactome | oraDE | 0.742 | 0.0000 | 0.109 | 0.0000 (0.0000) | 30 |

Clean null sets are sets without any member with an inserted effect (on average 2.1 Hallmark and 608 Reactome sets per replicate; Hallmark sets overlap heavily, so few are clean and the Hallmark clean-set rates rest on few sets). Rows for ORA of DE hits include only replicates with at least 3 DE hits. FRY and ORA called few clean null sets (largest proportion of called sets that were clean null: 0.0617); untouched sets that share non-null genes with a touched set are called at a higher rate, which is expected for both tests.

### 5.3 Meta-analysis, homogeneous shared effects (R = 30): observed FDR | power

| K | n | REML + Hartung-Knapp (default) | REML + truncated HK | REML + z | DL + z | Fixed effect + z | DL + Hartung-Knapp |
|---|---|---|---|---|---|---|---|
| 2 | 3 | 0.033 / 0.000 | 0.000 / 0.000 | 0.032 / 0.294 | 0.032 / 0.294 | 0.048 / 0.319 | 0.033 / 0.000 |
| 2 | 5 | 0.000 / 0.000 | 0.000 / 0.000 | 0.035 / 0.434 | 0.035 / 0.434 | 0.047 / 0.481 | 0.000 / 0.000 |
| 2 | 10 | 0.067 / 0.000 | 0.000 / 0.000 | 0.033 / 0.572 | 0.033 / 0.572 | 0.045 / 0.617 | 0.067 / 0.000 |
| 3 | 3 | 0.083 / 0.000 | 0.000 / 0.000 | 0.037 / 0.433 | 0.036 / 0.433 | 0.051 / 0.465 | 0.083 / 0.000 |
| 3 | 5 | 0.010 / 0.001 | 0.000 / 0.000 | 0.027 / 0.536 | 0.027 / 0.536 | 0.039 / 0.574 | 0.010 / 0.001 |
| 3 | 10 | 0.030 / 0.001 | 0.000 / 0.000 | 0.029 / 0.659 | 0.028 / 0.659 | 0.046 / 0.696 | 0.030 / 0.001 |
| 4 | 3 | 0.047 / 0.050 | 0.000 / 0.000 | 0.037 / 0.483 | 0.036 / 0.482 | 0.056 / 0.517 | 0.047 / 0.050 |
| 4 | 5 | 0.048 / 0.170 | 0.000 / 0.000 | 0.037 / 0.604 | 0.035 / 0.603 | 0.049 / 0.636 | 0.048 / 0.170 |
| 4 | 10 | 0.046 / 0.369 | 0.000 / 0.000 | 0.040 / 0.712 | 0.039 / 0.712 | 0.058 / 0.744 | 0.046 / 0.368 |
| 5 | 3 | 0.044 / 0.256 | 0.000 / 0.000 | 0.042 / 0.535 | 0.041 / 0.534 | 0.064 / 0.570 | 0.044 / 0.256 |
| 5 | 5 | 0.048 / 0.384 | 0.000 / 0.000 | 0.037 / 0.634 | 0.036 / 0.633 | 0.058 / 0.669 | 0.048 / 0.384 |
| 5 | 10 | 0.045 / 0.571 | 0.000 / 0.423 | 0.024 / 0.745 | 0.024 / 0.745 | 0.035 / 0.766 | 0.045 / 0.570 |

### 5.4 Meta-analysis, heterogeneous effects (tau = 0.4; R = 20): observed FDR | power

| K | n | REML + Hartung-Knapp (default) | REML + truncated HK | REML + z | DL + z | Fixed effect + z | DL + Hartung-Knapp |
|---|---|---|---|---|---|---|---|
| 2 | 5 | 0.025 / 0.000 | 0.000 / 0.000 | **0.164** / 0.316 | **0.164** / 0.316 | **0.204** / 0.455 | 0.025 / 0.000 |
| 2 | 10 | 0.050 / 0.000 | 0.000 / 0.000 | **0.222** / 0.396 | **0.222** / 0.396 | **0.298** / 0.620 | 0.050 / 0.000 |
| 3 | 5 | 0.050 / 0.000 | 0.000 / 0.000 | **0.114** / 0.362 | **0.113** / 0.362 | **0.187** / 0.541 | 0.050 / 0.000 |
| 3 | 10 | 0.100 / 0.000 | 0.000 / 0.000 | **0.149** / 0.422 | **0.148** / 0.423 | **0.282** / 0.672 | 0.100 / 0.000 |
| 4 | 5 | 0.050 / 0.000 | 0.000 / 0.000 | **0.092** / 0.409 | **0.092** / 0.409 | **0.180** / 0.593 | 0.050 / 0.000 |
| 4 | 10 | 0.132 / 0.000 | 0.000 / 0.000 | **0.106** / 0.459 | **0.105** / 0.465 | **0.267** / 0.710 | 0.132 / 0.000 |
| 5 | 5 | 0.038 / 0.008 | 0.000 / 0.000 | **0.081** / 0.446 | **0.080** / 0.445 | **0.176** / 0.628 | 0.038 / 0.008 |
| 5 | 10 | 0.055 / 0.030 | 0.000 / 0.000 | **0.099** / 0.495 | **0.098** / 0.498 | **0.266** / 0.733 | 0.055 / 0.029 |

Bold: mean FDP - 2 SE > 0.05. With homogeneous effects, 72 of 72 variant-settings passed the FDR criterion. The default REML + Hartung-Knapp had power 0.0000-0.0000 with K = 2 and 0.0002-0.0011 with K = 3 at any group size, against 0.294-0.572 and 0.433-0.659 for REML + z; with K = 4 its power was 0.050-0.369 and with K = 5 0.256-0.571. When effects differ between datasets, every z-based model and the fixed-effect model had observed FDR 0.080-0.298, while the observed FDR of the Hartung-Knapp variants was not significantly above 0.05 (largest mean FDP 0.132, from very few calls) but their power was at most 0.030. The truncated Hartung-Knapp test had zero power except K = 5 with 10 per group (homogeneous: 0.423).

### 5.5 Four psoriasis cohorts: shared signature (meta FDR <= 0.05, |pooled log2FC| >= 1, same direction) and canonical genes (of 20)

Setup of `PSO_SETUP` in `tests/engine/run_checks.py` (GSE54456 session; GSE121212 and GSE186063 paired by patient; GSE83645 paired by patient), minimum 3 datasets, prediction-interval criterion off. Entries: signature genes (canonical genes recovered); values identical to metafor (PAR-26 to PAR-31).

| Cohorts | REML + Hartung-Knapp (default) | REML + truncated HK | REML + z | DL + z | Fixed effect + z | DL + Hartung-Knapp |
|---|---|---|---|---|---|---|
| all four | 1773 (20) | 1721 (20) | 2157 (20) | 2151 (20) | 2212 (20) | 1778 (20) |
| without GSE54456 | 1 (0) | 0 (0) | 1769 (20) | 1788 (20) | 1730 (20) | 1 (0) |
| without GSE121212 | 1 (0) | 0 (0) | 2119 (20) | 2126 (20) | 2296 (20) | 1 (0) |
| without GSE186063 | 909 (4) | 0 (0) | 2306 (20) | 2301 (20) | 2346 (20) | 853 (4) |
| without GSE83645 | 0 (0) | 0 (0) | 2205 (20) | 2206 (20) | 2234 (20) | 0 (0) |

With all four cohorts the default gives 1773 signature genes and all 20 canonical genes; with any one cohort left out (three datasets, Hartung-Knapp t on 2 df) it keeps 1, 1, 909, 0 genes, whereas REML + z keeps 1769, 2119, 2306, 2205 and recovers all 20 canonical genes in every subset. This matches the simulation: with three datasets the Hartung-Knapp test rarely reaches meta FDR 0.05 over about 19,000 genes.

## 7. Not assessed

- Behaviour in real browsers: all runs used JavaScriptCore with the stub DOM and stub Plotly.
- GSEA at 3 and 5 per group: sample permutation is not available (fewer than 1,000 labelings) and the gene-permutation fallback reports no FDR in this build, so there is nothing to calibrate; the fallback itself was not re-run.
- Paired nulls, three-group moderated F and Discovery random-set path, co-expression FDR, RRHO, consensus and signature transfer: not repeated on the frozen build (outside this track; VAL_CAL results apply to code paths not changed by the freeze).
- Meta-analysis on the non-skin pool or mixing tissues: not run; the meta-analysis null and simulation use skin pools only.
- Heterogeneous-effect meta-analysis simulation at 3 per group and with other values of tau: only tau = 0.4 at 5 and 10 per group was run.
- ORA of DE hits under the null: only descriptive (it could run only in sessions with at least 3 false DE hits).
- Covariate-adjusted or paired designs in the thinning simulation: not run.

## 8. Deviations from the requested design

- Replicates: gene level and gene sets R = 100 per pool and size (requested >= 40); meta-analysis null R = 120 (>= 60); thinning R = 30 per setting (>= 20); heterogeneous meta-analysis scenario (added) R = 20 at 5 and 10 per group only.
- Meta-analysis K = 5 (added on request) used the four skin pools plus a second disjoint GSE54456 split; K = 2-4 used different pools chosen at random per replicate.
- GSEA with sample permutation was run only at 10 per group (n allows it only there).
- The non-skin pool, GSE57945 non-IBD ileum (42 samples), consists of symptomatic children without intestinal inflammation rather than healthy volunteers.
- R parity of per-gene p-values, FRY and CAMERA used seeds 0-4 of every cell (75 of 1,500 splits); counts, kept genes and fractions were compared on all 1,500. metafor parity on meta-analysis null dumps used 1,550 genes per replicate (1,500 random plus the 50 smallest p) in 12 replicates.
- seqgendiff is not available as a conda package; it was installed from the CRAN source (1.2.4) into a session library, with sva 3.58.0 from Bioconductor (conda) and clue, irlba, matchingR and pdist from CRAN.

## 9. Files

- `CALIBRATION_FREEZE_checks.csv`: all 499 checks (check_id, name, dataset, reference, n_compared, max_deviation, tolerance, pass, note).
- `CALIBRATION_FREEZE_null_replicates.csv`: one row per null session (pool, n, seed, samples, every output).
- `CALIBRATION_FREEZE_meta_null_replicates.csv`: one row per meta-analysis null replicate and variant.
- `CALIBRATION_FREEZE_sim_gene_replicates.csv`, `..._sim_set_replicates.csv`, `..._sim_meta_replicates.csv`, `..._sim_meta_het_replicates.csv`: per-replicate simulation results; `CALIBRATION_FREEZE_power_summary.csv`, `CALIBRATION_FREEZE_set_sim_summary.csv`: aggregated FDR and power; `CALIBRATION_FREEZE_real_meta.csv`: four-cohort results.
- `calibration_freeze_code.zip`: drivers (JavaScript and Python), R references (`ref_null.R`, `ref_meta.R`, `ref_real_meta.R`, `thin_sims.R`) and analysis scripts.

