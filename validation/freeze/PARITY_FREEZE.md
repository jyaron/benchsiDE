# R parity of benchsiDE 1.0.0-rc.1 (frozen build)

Build under test: `benchside/index.html`, APP_VERSION 1.0.0-rc.1 (2026-10-03); SHA-256 of the last `<script>` block cadbda80bd67db0e3cb15fcbcf138c23fb75cb24ff17cf6c28a79c68886c7d5b (computed by the run and asserted equal to the frozen value). The build was not modified. The application script was executed in macOS JavaScriptCore (`jsc`) with the repository's stub page and stub Plotly (`tests/engine/harness.py`); this is not a browser, so rendering and browser-specific behaviour are not tested here.

References: R version 4.5.3 (2026-03-11); limma 3.66.0; edgeR 4.8.2; metafor 5.0.1; fgsea 1.36.2; WGCNA 1.74 functions sourced from the CRAN source tarball (SHA-256 4f61e117b82287f01b61bf936d5c5ea3814c9c4b0d0b90d5016dda77f1f955d6), not installed. Every reference value was computed in R by package functions on the same input files (count matrices, designs) or, for the gene-set, meta-analysis and preservation tests, on the inputs the application itself used (set membership, per-dataset log2FC and variances, module expression values). All checks were produced in one end-to-end run (`run_all.py`, 320 s wall time, at most 4 parallel workers); each row of `PARITY_FREEZE_checks.csv` compares application output with the reference computed in that run.

Tolerances (VAL_DE.md): log2FC, t, SE, group means, log-CPM, TMM factors, enrichment scores and preservation statistics 1e-8 absolute; p-values, FDR, F, Q and d0 1e-6 relative; kept-gene sets and significant-gene counts exact.

## 1. Summary

- 1096 checks were run; 1085 passed and 11 did not.
  - `GSE143688_all.c1_welch.t`: t statistic: Aldara_KO vs Aldara_WT welch — max deviation 3.35 (tolerance 1e-08). 
  - `GSE143688_all.c1_welch.p`: p-value: Aldara_KO vs Aldara_WT welch — max deviation 0.992 (tolerance 1e-06). 
  - `GSE143688_all.c1_welch.q`: BH FDR: Aldara_KO vs Aldara_WT welch — max deviation 0.954 (tolerance 1e-06). 
  - `GSE143688_all.c1_welch.sigcounts`: significant up/down counts at FDR 0.01/0.05/0.10 x |log2FC| >= 0/1: Aldara_KO vs Aldara_WT welch — max deviation 1.0 (tolerance 0.0). FDR<=0.05,|lfc|>=1 up/down app 268/520, R 268/520; t.test errors set to p=1: 1
  - `GSE121212.fry_all_cov.p`: p-value (directional): FRY, fit=all +Patient, Hallmark (human, 50 sets) — max deviation 0.183 (tolerance 1e-06). limma fry itself changes by up to 0.165 (directional p) and 0.32 (mixed p), relative, when the covariate columns of the same design are reversed
  - `GSE121212.fry_all_cov.fdr`: FDR (directional): FRY, fit=all +Patient, Hallmark (human, 50 sets) — max deviation 0.183 (tolerance 1e-06). limma fry itself changes by up to 0.165 (directional p) and 0.32 (mixed p), relative, when the covariate columns of the same design are reversed
  - `GSE121212.fry_all_cov.pm`: p-value (mixed): FRY, fit=all +Patient, Hallmark (human, 50 sets) — max deviation 0.125 (tolerance 1e-06). limma fry itself changes by up to 0.165 (directional p) and 0.32 (mixed p), relative, when the covariate columns of the same design are reversed
  - `GSE121212.fry_all_cov.fdrm`: FDR (mixed): FRY, fit=all +Patient, Hallmark (human, 50 sets) — max deviation 0.125 (tolerance 1e-06). limma fry itself changes by up to 0.165 (directional p) and 0.32 (mixed p), relative, when the covariate columns of the same design are reversed
  - `meta.REML_knha.I2`: I^2 (%): REML, Hartung-Knapp — max deviation 66.6 (tolerance 1e-06). app I2 = (Q - df)/Q for every model; metafor I2 = 100 tau2/(tau2 + s2) from the model's tau2
  - `meta.REML_z.I2`: I^2 (%): REML, z test — max deviation 66.6 (tolerance 1e-06). app I2 = (Q - df)/Q for every model; metafor I2 = 100 tau2/(tau2 + s2) from the model's tau2
  - `meta.REML_adhoc.I2`: I^2 (%): REML, truncated Hartung-Knapp — max deviation 66.6 (tolerance 1e-06). app I2 = (Q - df)/Q for every model; metafor I2 = 100 tau2/(tau2 + s2) from the model's tau2
- Gene-level DE (moderated t, voom, Welch; two-group and multi-group with the all-groups and two-groups-only fits; categorical, continuous and patient/pair covariates; CPM-filter runs): 612 checks over 31 distinct run tags in 8 datasets. Maximum deviations: log2FC 2.0e-13, t 3.4e-12 (Welch 3.4e+00), SE 6.0e-13, p 6.7e-11 relative, FDR 6.7e-11 relative, group means 2.0e-13, d0 4.0e-13 relative.
- Significant-gene counts (up and down at FDR 0.01/0.05/0.10 × |log2FC| ≥ 0/1, 12 counts per run): exact in 77 of 78 DE runs, 12 of 12 interaction fits and 7 of 7 meta-analysis configurations.
- Filtering and normalization: kept-gene sets identical in 16 of 16 sessions (filterByExpr and CPM rule); TMM factors 2.2e-16; log-CPM matrices 6.2e-15 over 20,657,427 values.
- Moderated F (4 multi-group datasets): F, p and FDR 1.4e-11 relative; d0 5.7e-14; counts exact in 4 of 4.
- Interaction contrasts (GSE143688, 74 samples; genotype × treatment with and without Day, all eight groups and four groups; treatment × time on a derived Cell×Day factor): log2FC 1.8e-13, t 4.6e-13, p 9.3e-12 relative, simple effects 1.8e-13.
- Gene-set tests (7 datasets, Hallmark): ORA p 6.9e-11; CAMERA p 1.5e-11, inter-gene correlation 3.6e-14; FRY p (excluding the GSE121212 all-groups + Patient fit) 2.1e-10; GSEA enrichment scores 4.6e-14; gene-permutation p within Monte-Carlo error of fgseaSimple (max |z| 3.0e+00).
- Meta-analysis (four psoriasis cohorts, 19,401 genes in ≥ 3 cohorts; REML, DerSimonian-Laird and fixed effect; z, Hartung-Knapp and truncated Hartung-Knapp): pooled log2FC 4.9e-10, SE 5.6e-10, p 4.5e-08 relative, CI 1.9e-09, prediction interval 6.1e-09, Q 2.1e-15 relative, tau² 1.1e-09. I² agrees for DerSimonian-Laird and fixed effect (5.7e-14) but not for REML (6.7e+01 percentage points).
- Hub module preservation (15 hub × mode × test-dataset combinations, 4 hubs, 3 test datasets): the seven statistics 7.4e-14, Z-scores and Zsummary 5.8e-13.

## 2. Methods

**Application side.** For each dataset the session was loaded with species, built-in annotation, TMM and `fmode = fbe` (re-run with `fmode = cpm`, CPM ≥ 1 in ≥ 2 samples); every control was set explicitly. DE results were read from `computeDE()` (globals `deMethod`, `DEFIT` = `all` or `pair`, `DECOVSEL`), the moderated F from `computeModF()`, interaction contrasts from `lmContrastDE()`, ORA from `runEnrichment()` (query DE up or down at FDR ≤ 0.05, |log2FC| ≥ 1), FRY from `runFRY()`, GSEA from `runGSEA()` (gene permutation forced, then automatic mode), CAMERA from `discoveryRun()` with only the two contrasted groups included, the meta-analysis from `runMeta()` with the `PSO_SETUP` of `tests/engine/run_checks.py` (minimum 3 datasets), and preservation from `preservationInDataset()` after `computeHubs()` (2,000 most variable genes, β = 6, unsigned) and `hubNeighborhood(hub, 15, mode)`, 50 random sets each.

**R side.** `DGEList` → `filterByExpr(group = factor)` or `rowSums(cpm >= 1) >= 2` → `y[keep, , keep.lib.sizes = FALSE]` → `calcNormFactors(TMM)` → `cpm(log = TRUE, prior.count = 2)`. Moderated t: `lmFit(~ group [+ covariates])` on the fitted samples (all groups, or the two groups), `eBayes(trend = TRUE)`, `topTable(coef = group B)`. voom: `voom(y[, samples], design)`, `lmFit`, `eBayes`; group means as weighted coefficients of `~ 0 + group`. Welch: `t.test(var.equal = FALSE)` per gene; genes that `t.test()` rejects as constant were given t = 0, p = 1 (the application's convention for zero variance). Moderated F: `topTable(coef = 2:k)`. Interactions: `contrasts.fit` for the moderated t and for voom without covariates; for voom with Day the interaction was fitted as a coefficient of a reparametrised `~ 0 + cell + day` design (exact per-gene weighted fit; validation/freeze/INTERACTION.md), cross-checked against `~ trt * gen + day`. ORA: `phyper(k - 1, m, N - m, q, lower.tail = FALSE)` with R's own DE calls and the application's set mapping and universe, `p.adjust(BH)`. FRY: `fry()` on the application's sample order and design-column order (FRY's robust standardization depends on both). CAMERA: `camera(inter.gene.cor = NA)` on the same candidate modules. GSEA: `fgsea::calcGseaStat(gseaParam = 1)` on the R limma-trend t ranking; gene-permutation p against `fgseaSimple(nperm = 10000)`. Meta-analysis: `metafor::rma(yi, vi, method, test)` per gene on the application's yi, vi (REML with `threshold = 1e-12`), prediction interval `predict(predtype = 'Riley')`; the variances themselves were recomputed in R from each dataset's log2FC, t and df with `pt`/`qnorm`, and each cohort's log2FC, t, SE and p were checked against the R limma fit of the same contrast. Preservation: WGCNA `.coreCalcForExpr` with the reference module data and the test data of the observed module and of each random set exported by the application (within-group mode: group means removed first), Z = (observed − mean)/sd over the 50 sets.

## 3. Filtering, TMM and log-CPM

| Dataset | filter | kept genes (app) | kept-set difference | TMM max dev | log-CPM max dev | log-CPM values |
|---|---|---|---|---|---|---|
| GSE63310 | fbe | 16624 | 0.0 | 0.0 | 5.33e-15 | 149616 |
| GSE63310 | cpm | 14490 | 0.0 | 1.11e-16 | 5.77e-15 | 130410 |
| GSE54456 | fbe | 19518 | 0.0 | 2.22e-16 | 5.33e-15 | 3337578 |
| GSE54456 | cpm | 19902 | 0.0 | 2.22e-16 | 5.77e-15 | 3403242 |
| GSE121212 | fbe | 23168 | 0.0 | 2.22e-16 | 6.22e-15 | 3336192 |
| GSE121212 | cpm | 21755 | 0.0 | 1.11e-16 | 5.33e-15 | 3132720 |
| GSE186063 | fbe | 21112 | 0.0 | 2.22e-16 | 5.33e-15 | 1372280 |
| GSE186063 | cpm | 20505 | 0.0 | 2.22e-16 | 5.33e-15 | 1332825 |
| GSE41745 | fbe | 18036 | 0.0 | 2.22e-16 | 5.33e-15 | 90180 |
| GSE41745 | cpm | 17409 | 0.0 | 1.11e-16 | 3.55e-15 | 87045 |
| GSE83645 | fbe | 20956 | 0.0 | 1.11e-16 | 5.33e-15 | 523900 |
| GSE83645 | cpm | 18591 | 0.0 | 0.0 | 5.33e-15 | 464775 |
| GSE143688 | fbe | 14885 | 0.0 | 0.0 | 5.77e-15 | 535860 |
| GSE143688 | cpm | 14603 | 0.0 | 0.0 | 5.33e-15 | 525708 |
| GSE143688_all | fbe | 15492 | 0.0 | 2.22e-16 | 5.33e-15 | 1146408 |
| GSE143688_all | cpm | 14712 | 0.0 | 2.22e-16 | 5.33e-15 | 1088688 |

## 4. Gene-level differential expression

Maximum deviation per run (absolute for log2FC, t, SE and group means; relative for p and FDR). Counts: FDR ≤ 0.05 and |log2FC| ≥ 1, up/down, application vs R; all 12 counts per run are compared in the CSV.

| Dataset | run | genes | log2FC | t | SE | p | FDR | means | d0 | up/down app (R) | 12 counts exact |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE63310 | LP vs Basal mod fit=all | 16624 | 8.88e-15 | 1.5e-12 | 3.48e-14 | 5.26e-12 | 5.26e-12 | 5.33e-15 | 5.72e-14 | 2763/3305 (2763/3305) | yes |
| GSE63310 | LP vs Basal voom fit=all | 16624 | 5.37e-14 | 7.85e-13 | 3.02e-14 | 1.45e-12 | 1.45e-12 | 4.62e-14 | 6.09e-14 | 2781/3301 (2781/3301) | yes |
| GSE63310 | LP vs Basal mod fit=pair | 16624 | 6.22e-15 | 1.09e-12 | 4.81e-14 | 1.59e-12 | 1.59e-12 | 5.33e-15 | 6.01e-14 | 2780/3357 (2780/3357) | yes |
| GSE63310 | LP vs Basal voom fit=pair | 16624 | 4.93e-14 | 1.56e-12 | 2.46e-14 | 1.4e-12 | 1.4e-12 | 4.44e-14 | 6.47e-14 | 2810/3338 (2810/3338) | yes |
| GSE63310 | LP vs Basal welch | 16624 | 6.22e-15 | 1.96e-12 | 2.19e-15 | 9.22e-13 | 9.22e-13 | 5.33e-15 | – | 1968/2418 (1968/2418) | yes |
| GSE63310 | ML vs Basal mod fit=all | 16624 | 7.11e-15 | 1.38e-12 | 3.47e-14 | 1.45e-12 | 1.45e-12 | 5.33e-15 | 5.74e-14 | 2962/3307 (2962/3307) | yes |
| GSE63310 | ML vs Basal voom fit=all | 16624 | 6.84e-14 | 8.38e-13 | 3.24e-14 | 1.47e-12 | 1.47e-12 | 6.39e-14 | 5.97e-14 | 2973/3329 (2973/3329) | yes |
| GSE63310 | ML vs Basal mod fit=pair | 16624 | 5.33e-15 | 1.1e-12 | 4.34e-14 | 1.24e-12 | 1.24e-12 | 5.33e-15 | 6.08e-14 | 2969/3315 (2969/3315) | yes |
| GSE63310 | ML vs Basal voom fit=pair | 16624 | 4.53e-14 | 9.95e-13 | 5.44e-14 | 1.66e-12 | 1.66e-12 | 4.62e-14 | 6.11e-14 | 2966/3301 (2966/3301) | yes |
| GSE63310 | ML vs Basal welch | 16624 | 5.33e-15 | 2.84e-12 | 2.41e-15 | 9.47e-13 | 9.47e-13 | 5.33e-15 | – | 2282/2422 (2282/2422) | yes |
| GSE63310 | ML vs LP mod fit=all | 16624 | 9.1e-15 | 2.01e-13 | 3.49e-14 | 5.66e-12 | 1.49e-12 | 5.33e-15 | 5.72e-14 | 1418/1387 (1418/1387) | yes |
| GSE63310 | ML vs LP voom fit=all | 16624 | 8.55e-14 | 4.07e-13 | 3.11e-14 | 4.66e-12 | 4.66e-12 | 6.39e-14 | 6.25e-14 | 1331/1290 (1331/1290) | yes |
| GSE63310 | ML vs LP mod fit=pair | 16624 | 5.33e-15 | 1.35e-13 | 3.11e-14 | 2.76e-12 | 1.14e-12 | 5.33e-15 | 2.8e-14 | 1299/1189 (1299/1189) | yes |
| GSE63310 | ML vs LP voom fit=pair | 16624 | 5.33e-14 | 2.49e-13 | 2.55e-14 | 1.29e-12 | 1.29e-12 | 5.15e-14 | 3.56e-14 | 1287/1212 (1287/1212) | yes |
| GSE63310 | ML vs LP welch | 16624 | 7.11e-15 | 1.19e-12 | 2.39e-15 | 2.24e-12 | 9.31e-13 | 5.33e-15 | – | 0/0 (0/0) | yes |
| GSE63310 | LP vs Basal mod fit=all (CPM filter) | 14490 | 7.55e-15 | 2.06e-12 | 5.31e-14 | 1.58e-12 | 1.57e-12 | 5.33e-15 | 3.94e-14 | 2392/2804 (2392/2804) | yes |
| GSE54456 | Psoriasis_skin vs normal_skin mod fit=all | 19518 | 1.42e-14 | 2.61e-13 | 4.16e-16 | 6.76e-12 | 5.43e-12 | 7.11e-15 | 2.25e-14 | 1293/1513 (1293/1513) | yes |
| GSE54456 | Psoriasis_skin vs normal_skin voom fit=all | 19518 | 3.02e-14 | 5.83e-13 | 8.77e-15 | 7.52e-12 | 7.52e-12 | 2.31e-14 | 2.48e-14 | 1340/1653 (1340/1653) | yes |
| GSE54456 | Psoriasis_skin vs normal_skin welch | 19518 | 8.88e-15 | 2.82e-13 | 2.22e-16 | 5.43e-12 | 5.43e-12 | 7.11e-15 | – | 1293/1513 (1293/1513) | yes |
| GSE54456 | Psoriasis_skin vs normal_skin mod fit=all (CPM filter) | 19902 | 1.24e-14 | 2.88e-13 | 6.11e-16 | 3.99e-12 | 3.99e-12 | 5.33e-15 | 2.5e-14 | 1344/1696 (1344/1696) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional mod fit=all | 23168 | 9.33e-15 | 1.1e-13 | 1.14e-15 | 1.16e-11 | 1.16e-11 | 4.44e-15 | 2.38e-14 | 1401/1845 (1401/1845) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional voom fit=all | 23168 | 3.61e-14 | 3.22e-13 | 7.22e-15 | 7.3e-12 | 7.3e-12 | 3.02e-14 | 2.35e-14 | 1657/2481 (1657/2481) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional mod fit=pair | 23168 | 9.77e-15 | 1.41e-13 | 3.58e-15 | 8.47e-12 | 8.47e-12 | 4.44e-15 | 2.32e-14 | 1401/1846 (1401/1846) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional voom fit=pair | 23168 | 3.94e-14 | 2.45e-13 | 6.99e-15 | 1.04e-11 | 8.51e-12 | 3.2e-14 | 2.38e-14 | 1664/2499 (1664/2499) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional welch | 23168 | 7.11e-15 | 1.55e-13 | 2.78e-16 | 8.44e-12 | 8.44e-12 | 4.44e-15 | – | 1401/1846 (1401/1846) | yes |
| GSE121212 | PSO_lesional vs CTRL_healthy mod fit=all | 23168 | 6.11e-15 | 1.34e-13 | 9.16e-16 | 4.53e-12 | 4.53e-12 | 5.33e-15 | 2.31e-14 | 1750/2640 (1750/2640) | yes |
| GSE121212 | PSO_lesional vs CTRL_healthy voom fit=all | 23168 | 6.57e-14 | 3.53e-13 | 1.2e-14 | 9.11e-12 | 9.11e-12 | 6.93e-14 | 2.4e-14 | 1966/3419 (1966/3419) | yes |
| GSE121212 | PSO_lesional vs CTRL_healthy mod fit=pair | 23168 | 6.55e-15 | 1.3e-13 | 1.9e-15 | 8.49e-12 | 8.49e-12 | 5.33e-15 | 1.37e-14 | 1751/2640 (1751/2640) | yes |
| GSE121212 | PSO_lesional vs CTRL_healthy voom fit=pair | 23168 | 5.4e-14 | 5.54e-13 | 1.68e-14 | 8.26e-12 | 8.26e-12 | 4.44e-14 | 2.13e-14 | 1967/3424 (1967/3424) | yes |
| GSE121212 | PSO_lesional vs CTRL_healthy welch | 23168 | 7.11e-15 | 1.33e-13 | 2.5e-16 | 8.6e-12 | 8.6e-12 | 5.33e-15 | – | 1750/2640 (1750/2640) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional mod fit=all +Patient | 23168 | 3.13e-14 | 8.37e-13 | 5.58e-15 | 8.67e-12 | 8.67e-12 | 4.44e-15 | 1.97e-14 | 1420/1825 (1420/1825) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional mod fit=pair +Patient | 23168 | 1.07e-14 | 2.62e-13 | 9.16e-15 | 4.32e-12 | 4.32e-12 | 4.44e-15 | 2.69e-14 | 1420/1825 (1420/1825) | yes |
| GSE121212 | PSO_lesional vs PSO_non_lesional mod fit=all (CPM filter) | 21755 | 7.99e-15 | 1.07e-13 | 1.44e-15 | 5.51e-12 | 5.51e-12 | 5.33e-15 | 2.04e-14 | 1293/1890 (1293/1890) | yes |
| GSE186063 | lesion vs non-lesion mod fit=all | 21112 | 7.31e-15 | 8.93e-14 | 1.55e-15 | 2.15e-11 | 2.15e-11 | 5.33e-15 | 2.57e-14 | 1016/1399 (1016/1399) | yes |
| GSE186063 | lesion vs non-lesion voom fit=all | 21112 | 6.66e-14 | 2.54e-13 | 1.68e-14 | 8.55e-12 | 8.55e-12 | 6.39e-14 | 2.83e-14 | 1115/1703 (1115/1703) | yes |
| GSE186063 | lesion vs non-lesion mod fit=pair | 21112 | 8.31e-15 | 8.39e-14 | 1.94e-15 | 8.68e-12 | 8.68e-12 | 5.33e-15 | 2.32e-14 | 1015/1399 (1015/1399) | yes |
| GSE186063 | lesion vs non-lesion voom fit=pair | 21112 | 5.15e-14 | 2.08e-13 | 1.95e-14 | 6.72e-11 | 6.72e-11 | 6.39e-14 | 2.8e-14 | 1114/1708 (1114/1708) | yes |
| GSE186063 | lesion vs non-lesion welch | 21112 | 8.88e-15 | 1.01e-13 | 3.33e-16 | 8.41e-12 | 8.41e-12 | 5.33e-15 | – | 1015/1398 (1015/1398) | yes |
| GSE186063 | lesion vs non-lesion mod fit=all +Diag | 21112 | 8.27e-15 | 1.06e-13 | 1.75e-15 | 8.47e-12 | 8.47e-12 | 5.33e-15 | 2.42e-14 | 1016/1406 (1016/1406) | yes |
| GSE186063 | lesion vs non-lesion mod fit=pair +Diag | 21112 | 8.44e-15 | 1.21e-13 | 1.92e-15 | 8.66e-12 | 8.66e-12 | 5.33e-15 | 2.48e-14 | 1015/1407 (1015/1407) | yes |
| GSE186063 | lesion vs non-lesion mod fit=all +Sex+Age | 21112 | 9.88e-15 | 9.89e-14 | 1.78e-15 | 2.33e-11 | 2.33e-11 | 5.33e-15 | 2.52e-14 | 1019/1400 (1019/1400) | yes |
| GSE186063 | lesion vs non-lesion mod fit=pair +Sex+Age | 21112 | 1.07e-14 | 1.17e-13 | 1.86e-15 | 8.59e-12 | 8.59e-12 | 5.33e-15 | 2.86e-14 | 1015/1404 (1015/1404) | yes |
| GSE186063 | lesion vs non-lesion mod fit=all +Pair | 21112 | 6.55e-15 | 9.24e-14 | 2.44e-15 | 4.35e-12 | 4.35e-12 | 5.33e-15 | 2.7e-14 | 997/1429 (997/1429) | yes |
| GSE186063 | lesion vs non-lesion mod fit=pair +Pair | 21112 | 6.55e-15 | 9.24e-14 | 2.44e-15 | 4.35e-12 | 4.35e-12 | 5.33e-15 | 2.7e-14 | 997/1429 (997/1429) | yes |
| GSE186063 | lesion vs non-lesion mod fit=all (CPM filter) | 20505 | 7.91e-15 | 8.97e-14 | 3.52e-15 | 8.42e-12 | 8.42e-12 | 7.11e-15 | 3.3e-14 | 1027/1385 (1027/1385) | yes |
| GSE41745 | lesional vs non_lesional mod fit=all | 18036 | 4.95e-15 | 2.95e-13 | 1.73e-13 | 2.86e-12 | 2.49e-12 | 3.55e-15 | 8.43e-14 | 634/301 (634/301) | yes |
| GSE41745 | lesional vs non_lesional voom fit=all | 18036 | 7.39e-14 | 2.38e-13 | 7.59e-14 | 3.51e-11 | 3.51e-11 | 6.26e-14 | 5.88e-14 | 627/182 (627/182) | yes |
| GSE41745 | lesional vs non_lesional welch | 18036 | 5.33e-15 | 1.64e-12 | 2.86e-15 | 7.5e-13 | 7.5e-13 | 3.55e-15 | – | 0/0 (0/0) | yes |
| GSE41745 | lesional vs non_lesional mod fit=all +Patient | 18036 | 6.77e-15 | 3.39e-12 | 6e-13 | 2.04e-11 | 2.04e-11 | 3.55e-15 | 3.96e-13 | 867/690 (867/690) | yes |
| GSE41745 | lesional vs non_lesional mod fit=all (CPM filter) | 17409 | 6.22e-15 | 2.63e-13 | 6.95e-14 | 1.9e-12 | 1.9e-12 | 5.33e-15 | 6.7e-14 | 636/306 (636/306) | yes |
| GSE83645 | psoriasis vs uninvolved mod fit=all | 20956 | 5.33e-15 | 6.06e-14 | 3.39e-15 | 7.53e-12 | 7.53e-12 | 5.33e-15 | 1.65e-14 | 932/919 (932/919) | yes |
| GSE83645 | psoriasis vs uninvolved voom fit=all | 20956 | 5.77e-14 | 4.09e-13 | 9.53e-14 | 3.22e-12 | 3.19e-12 | 6.75e-14 | 1.38e-14 | 873/1041 (873/1041) | yes |
| GSE83645 | psoriasis vs uninvolved welch | 20956 | 5.33e-15 | 9.95e-14 | 1.11e-15 | 1.78e-12 | 1.78e-12 | 5.33e-15 | – | 818/855 (818/855) | yes |
| GSE83645 | psoriasis vs uninvolved mod fit=all +Patient | 20956 | 8.88e-15 | 1.19e-13 | 3.94e-15 | 3.58e-12 | 3.58e-12 | 5.33e-15 | 2.19e-14 | 1026/1120 (1026/1120) | yes |
| GSE83645 | psoriasis vs uninvolved mod fit=all +Site | 20956 | 3.21e-14 | 2.42e-13 | 5.77e-15 | 8.42e-12 | 8.42e-12 | 5.33e-15 | 2.93e-14 | 600/280 (600/280) | yes |
| GSE83645 | psoriasis vs uninvolved mod fit=all +Patient+Site | 20956 | 3.51e-14 | 3.37e-13 | 7.55e-15 | 2.61e-12 | 2.61e-12 | 5.33e-15 | 3.15e-14 | 767/602 (767/602) | yes |
| GSE83645 | psoriasis vs uninvolved mod fit=all (CPM filter) | 18591 | 5.05e-15 | 8.88e-14 | 7.47e-15 | 3e-12 | 3e-12 | 5.33e-15 | 1.03e-14 | 859/748 (859/748) | yes |
| GSE143688 | Aldara vs Control mod fit=all | 14885 | 8.44e-15 | 1.33e-13 | 1.11e-15 | 6.51e-12 | 6.51e-12 | 5.33e-15 | 8.62e-15 | 2138/2414 (2138/2414) | yes |
| GSE143688 | Aldara vs Control voom fit=all | 14885 | 1.39e-13 | 1.44e-12 | 5.68e-14 | 6.58e-12 | 6.58e-12 | 1.39e-13 | 1.03e-14 | 2134/2394 (2134/2394) | yes |
| GSE143688 | Aldara vs Control welch | 14885 | 7.11e-15 | 1.81e-13 | 4.44e-16 | 1.06e-11 | 1.06e-11 | 5.33e-15 | – | 2137/2413 (2137/2413) | yes |
| GSE143688 | Aldara vs Control mod fit=all +Day | 14885 | 5.83e-15 | 1.49e-13 | 1.11e-15 | 6.43e-12 | 6.43e-12 | 5.33e-15 | 1.1e-14 | 2138/2415 (2138/2415) | yes |
| GSE143688 | Aldara vs Control mod fit=all +Line+Day | 14885 | 5.83e-15 | 1.58e-13 | 1.55e-15 | 1.07e-11 | 1.07e-11 | 5.33e-15 | 1.59e-14 | 2139/2415 (2139/2415) | yes |
| GSE143688 | Aldara vs Control mod fit=all (CPM filter) | 14603 | 7.16e-15 | 1.63e-13 | 6.66e-16 | 6.48e-12 | 6.48e-12 | 5.33e-15 | 7.53e-15 | 2148/2363 (2148/2363) | yes |
| GSE143688_all | Aldara_WT vs Control_WT mod fit=all | 15492 | 5.33e-15 | 1.1e-13 | 8.88e-16 | 8.36e-12 | 8.36e-12 | 5.33e-15 | 9.9e-15 | 2280/2589 (2280/2589) | yes |
| GSE143688_all | Aldara_WT vs Control_WT voom fit=all | 15492 | 2e-13 | 6.41e-13 | 1.57e-14 | 5.25e-11 | 5.25e-11 | 2.03e-13 | 9.47e-15 | 2298/2618 (2298/2618) | yes |
| GSE143688_all | Aldara_WT vs Control_WT mod fit=pair | 15492 | 5.77e-15 | 2.77e-13 | 4.88e-15 | 1.77e-12 | 1.77e-12 | 5.33e-15 | 2.93e-14 | 2278/2587 (2278/2587) | yes |
| GSE143688_all | Aldara_WT vs Control_WT voom fit=pair | 15492 | 1.08e-13 | 5.15e-13 | 1.62e-14 | 2.27e-11 | 2.27e-11 | 1.07e-13 | 3.38e-14 | 2298/2618 (2298/2618) | yes |
| GSE143688_all | Aldara_WT vs Control_WT welch | 15492 | 6.22e-15 | 1.92e-13 | 7.22e-16 | 1.85e-12 | 1.85e-12 | 5.33e-15 | – | 2275/2581 (2275/2581) | yes |
| GSE143688_all | Aldara_KO vs Aldara_WT mod fit=all | 15492 | 4.77e-15 | 8.17e-14 | 1.33e-15 | 8.37e-12 | 8.37e-12 | 3.55e-15 | 9.78e-15 | 401/634 (401/634) | yes |
| GSE143688_all | Aldara_KO vs Aldara_WT voom fit=all | 15492 | 1.24e-13 | 3.61e-13 | 3.06e-14 | 1.02e-11 | 1.02e-11 | 1.14e-13 | 8.26e-15 | 413/613 (413/613) | yes |
| GSE143688_all | Aldara_KO vs Aldara_WT mod fit=pair | 15492 | 7.33e-15 | 9.33e-14 | 4.44e-15 | 2.29e-12 | 2.29e-12 | 3.55e-15 | 2.21e-14 | 324/564 (324/564) | yes |
| GSE143688_all | Aldara_KO vs Aldara_WT voom fit=pair | 15492 | 5.85e-14 | 2.55e-13 | 1.73e-14 | 2.1e-12 | 2.1e-12 | 6.75e-14 | 2.24e-14 | 334/539 (334/539) | yes |
| GSE143688_all | Aldara_KO vs Aldara_WT welch | 15492 | 5.33e-15 | 3.35 | 7.77e-16 | 0.992 | 0.954 | 3.55e-15 | – | 268/520 (268/520) | NO |
| GSE143688_all | Aldara_WT vs Control_WT mod fit=all +Day | 15492 | 1.69e-14 | 2.41e-13 | 9.99e-16 | 8.34e-12 | 8.34e-12 | 5.33e-15 | 9.66e-15 | 2283/2589 (2283/2589) | yes |
| GSE143688_all | Aldara_WT vs Control_WT mod fit=pair +Day | 15492 | 5.33e-15 | 1.74e-13 | 3.55e-15 | 1.83e-12 | 1.83e-12 | 5.33e-15 | 2.02e-14 | 2282/2588 (2282/2588) | yes |
| GSE143688_all | Aldara_KO vs Aldara_WT mod fit=all +Day | 15492 | 1.13e-14 | 1.94e-13 | 9.58e-16 | 1.55e-11 | 1.55e-11 | 3.55e-15 | 9.54e-15 | 379/632 (379/632) | yes |
| GSE143688_all | Aldara_KO vs Aldara_WT mod fit=pair +Day | 15492 | 8.44e-15 | 1.59e-13 | 9.33e-15 | 1.04e-11 | 2.68e-12 | 3.55e-15 | 3.92e-14 | 341/591 (341/591) | yes |
| GSE143688_all | Aldara_WT vs Control_WT mod fit=all (CPM filter) | 14712 | 7.11e-15 | 7.99e-14 | 1.33e-15 | 8.51e-12 | 8.51e-12 | 4.44e-15 | 6.57e-15 | 2177/2419 (2177/2419) | yes |

Moderated F:

| Dataset | groups | F | p | FDR | d0 | genes at FDR 0.01/0.05/0.10, app (R) |
|---|---|---|---|---|---|---|
| GSE63310 | 3 | 1.94e-12 | 5.41e-13 | 5.41e-13 | 5.72e-14 | [8789, 11335, 12560] ([8789, 11335, 12560]) |
| GSE121212 | 6 | 8.85e-14 | 7.72e-12 | 7.71e-12 | 2.38e-14 | [18408, 20107, 20889] ([18408, 20107, 20889]) |
| GSE186063 | 3 | 8.12e-13 | 9.21e-12 | 9.2e-12 | 2.57e-14 | [9337, 11824, 13144] ([9337, 11824, 13144]) |
| GSE143688_all | 8 | 6.27e-14 | 1.41e-11 | 1.41e-11 | 9.9e-15 | [12097, 12940, 13336] ([12097, 12940, 13336]) |

## 5. Interaction contrasts (lmContrastDE), GSE143688

| Fit | log2FC | t | SE | p | FDR | simple effects | d0 | up/down app (R) | counts exact |
|---|---|---|---|---|---|---|---|---|---|
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), mod, all groups + Day | 1.6e-14 | 1.42e-13 | 2.66e-15 | 8.51e-12 | 8.51e-12 | 1.07e-14 | 9.66e-15 | 349/522 (349/522) | yes |
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), mod, all groups | 1.07e-14 | 9.68e-14 | 1.78e-15 | 8.4e-12 | 8.4e-12 | 7.11e-15 | 9.42e-15 | 287/470 (287/470) | yes |
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), mod, four groups + Day | 1.42e-14 | 1.47e-13 | 2.83e-15 | 6.5e-12 | 6.5e-12 | 9.77e-15 | 1.09e-14 | 340/508 (340/508) | yes |
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), mod, four groups | 8.88e-15 | 7.73e-14 | 2.28e-15 | 6.77e-12 | 6.77e-12 | 6.22e-15 | 1.05e-14 | 239/462 (239/462) | yes |
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), voom, all groups + Day | 6.4e-14 | 2.74e-13 | 1.46e-14 | 8.31e-12 | 8.31e-12 | – | 1.06e-14 | 312/449 (312/449) | yes |
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), voom, all groups | 9.15e-14 | 1.9e-13 | 1.38e-14 | 8.58e-12 | 8.58e-12 | 5.51e-14 | 8.75e-15 | 273/371 (273/371) | yes |
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), voom, four groups + Day | 6.03e-14 | 3.31e-13 | 1.73e-14 | 6.58e-12 | 6.58e-12 | – | 8.33e-15 | 307/453 (307/453) | yes |
| interaction (Aldara_KO - Control_KO) - (Aldara_WT - Control_WT), voom, four groups | 9.41e-14 | 2.74e-13 | 1.6e-14 | 7.8e-12 | 7.8e-12 | 8.7e-14 | 8.9e-15 | 250/368 (250/368) | yes |
| interaction (Aldara_WT_d7 - Control_WT_d7) - (Aldara_WT_d3 - Control_WT_d3), mod, all groups | 8.88e-15 | 9.26e-14 | 1.78e-15 | 8.65e-12 | 8.61e-12 | 7.11e-15 | 1.09e-14 | 585/806 (585/806) | yes |
| interaction (Aldara_WT_d7 - Control_WT_d7) - (Aldara_WT_d3 - Control_WT_d3), mod, four groups | 8.88e-15 | 9.77e-14 | 1.51e-14 | 7.29e-12 | 7.29e-12 | 6.22e-15 | 3.2e-14 | 396/521 (396/521) | yes |
| interaction (Aldara_WT_d7 - Control_WT_d7) - (Aldara_WT_d3 - Control_WT_d3), voom, all groups | 1.79e-13 | 4.57e-13 | 6.13e-14 | 8.59e-12 | 8.59e-12 | 1.76e-13 | 9.62e-15 | 526/743 (526/743) | yes |
| interaction (Aldara_WT_d7 - Control_WT_d7) - (Aldara_WT_d3 - Control_WT_d3), voom, four groups | 1.38e-13 | 3.14e-13 | 2.35e-14 | 9.31e-12 | 9.31e-12 | 1.14e-13 | 3.46e-14 | 371/569 (371/569) | yes |

voom with Day, four groups: factorial ~trt*gen+day agrees with the reparametrised fit to 2.9e-13 in t; contrasts.fit differs by up to 0.310 in t (765 vs 760 significant). The application reproduces the exact coefficient fit, not limma's approximate `contrasts.fit`.

GSE171012 treatment × time: not assessed. Whole-skin (CellType = bulk) samples have the Status levels Healthy, Psoriasis_PreTreatment, Psoriasis_SecukinumabTreatmentWeek12, Psoriasis_SecukinumabTreatmentWeek2, Psoriasis_SecukinumabTreatmentWeek4 (Healthy 11; Psoriasis_PreTreatment 15; Psoriasis_SecukinumabTreatmentWeek12 13; Psoriasis_SecukinumabTreatmentWeek2 14; Psoriasis_SecukinumabTreatmentWeek4 13 samples). Every psoriasis sample belongs to the secukinumab-treated series, healthy skin has a single Status level, and no healthy subject appears in the psoriasis series (overlap 0 subjects), so there is no untreated or second-arm time course. (Week12 − PreTreatment) − (Week2 − PreTreatment) reduces algebraically to Week12 − Week2, a simple contrast, not a difference of differences; no valid 2 × 2 can be formed from these groups.

## 6. Gene-set tests

| Dataset | ORA p (up, down) | ORA counts exact | FRY p / p mixed (all fits) | CAMERA p / FDR / cor | modules | GSEA ES (preranked, auto) | auto mode | gene-perm p, max abs z |
|---|---|---|---|---|---|---|---|---|
| GSE54456 | 5.06e-11, 2.98e-11 | yes | 1.5e-12 / 7.0e-12 | 1.47e-11 / 1.47e-11 / 2.11e-15 | 2905 | 8.33e-16, 4.6e-14 | samples | 2.41 |
| GSE121212 | 6.03e-11, 3.87e-11 | yes | 1.8e-01 / 1.2e-01 | 8.43e-12 / 4.16e-12 / 5.11e-15 | 2983 | 1.67e-15, 1.88e-14 | samples | 1.68 |
| GSE186063 | 5.29e-11, 1.89e-11 | yes | 3.3e-11 / 4.1e-11 | 8.46e-12 / 5.33e-12 / 2.66e-15 | 2968 | 8.88e-16, 9.21e-15 | samples | 1.68 |
| GSE83645 | 3.85e-11, 1.72e-11 | yes | 1.7e-11 / 5.0e-11 | 2.02e-12 / 5.82e-13 / 6.83e-15 | 2994 | 1.78e-15, 1.08e-14 | samples | 3.0 |
| GSE41745 | 3.67e-11, 1.81e-11 | yes | 2.1e-10 / 1.3e-10 | 7.4e-13 / 4.04e-14 / 3.59e-14 | 2576 | 2.39e-15, 2.39e-15 | genes | 2.02 |
| GSE63310 | 5.63e-11, 6.3e-11 | yes | 5.5e-11 / 5.6e-11 | 1.01e-12 / 1.01e-12 / 2.33e-14 | 2661 | 2.89e-15, 2.89e-15 | genes | 2.12 |
| GSE143688 | 6.9e-11, 5.65e-11 | yes | 4.7e-11 / 6.6e-11 | 6.15e-12 / 6.15e-12 / 4e-15 | 2680 | 8.33e-16, 2.79e-14 | samples | 1.77 |

Gene-permutation p-values: z = (p_app − p_fgsea)/SE with SE from the pooled p and half of each method's permutations as same-sign draws; criterion max |z| ≤ 4 over the 50 sets of each dataset. Sample-permutation p-values were not compared (section 10).

FRY, GSE121212 all six groups + Patient (100 design columns, 144 samples): p up to 0.183 relative, mixed p up to 0.125. sets 50 vs 50; R design 100 columns, 144 samples.

## 7. Meta-analysis

| Model, test | pooled log2FC | SE | p (rel) | 95% CI | 95% PI | Q (rel) | Q p (rel) | tau² | I² (points) | FDR (rel) | genes FDR ≤ 0.05 / and abs(log2FC) ≥ 1, app (R) | shared signature |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| REML, Hartung-Knapp | 4.86e-10 | 4.03e-10 | 1.27e-09 | 1.8e-09 | 5.1e-09 | 2.12e-15 | 2.49e-14 | 1.06e-09 | 66.6 | 1.27e-09 | 7,016 / 1,773 (7,016 / 1,773) | 1,773 |
| REML, z test | 4.86e-10 | 5.61e-10 | 4.49e-08 | 1.4e-09 | 6.1e-09 | 2.12e-15 | 2.49e-14 | 1.06e-09 | 66.6 | 4.49e-08 | 12,846 / 2,166 (12,846 / 2,166) | 2,157 |
| REML, truncated Hartung-Knapp | 4.86e-10 | 5.61e-10 | 5.42e-09 | 1.9e-09 | 6.1e-09 | 2.12e-15 | 2.49e-14 | 1.06e-09 | 66.6 | 5.42e-09 | 6,254 / 1,721 (6,254 / 1,721) | 1,721 |
| DerSimonian-Laird, Hartung-Knapp | 5.33e-15 | 4.44e-16 | 7.44e-13 | 1.1e-13 | 1.9e-12 | 2.12e-15 | 2.49e-14 | 3.55e-15 | 5.68e-14 | 7.44e-13 | 7,013 / 1,778 (7,013 / 1,778) | 1,778 |
| DerSimonian-Laird, z test | 5.33e-15 | 4.44e-16 | 3.35e-13 | 5.3e-15 | 1.8e-12 | 2.12e-15 | 2.49e-14 | 3.55e-15 | 5.68e-14 | 3.35e-13 | 13,175 / 2,163 (13,175 / 2,163) | 2,151 |
| DerSimonian-Laird, truncated Hartung-Knapp | 5.33e-15 | 4.44e-16 | 7.43e-13 | 1.1e-13 | 1.9e-12 | 2.12e-15 | 2.49e-14 | 3.55e-15 | 5.68e-14 | 7.43e-13 | 6,145 / 1,691 (6,145 / 1,691) | 1,691 |
| fixed effect, z test | 7.11e-15 | 1.11e-16 | 5.15e-13 | 7.1e-15 | – | 2.12e-15 | 2.49e-14 | 0.0 | 5.68e-14 | 5.14e-13 | 16,636 / 2,232 (16,636 / 2,232) | 2,212 |

Inputs: each cohort's log2FC, t, SE and p agree with the R limma fit of the same contrast (log2FC 1.4e-14, t 2.8e-13, SE 5.6e-15, p 8.2e-12 relative), and the variances carried into the meta-analysis agree with R's `pt`/`qnorm` recomputation to 7.87e-08 relative over 77,193 gene × cohort values. The REML differences (order 1e-9) reflect the convergence tolerances of the two iterative solutions.

## 8. Hub module preservation

| Hub | correlation | test dataset | matched members | statistics max dev | Z max dev | Zsummary (app) |
|---|---|---|---|---|---|---|
| IL36A | all samples | GSE121212 | 16 | 1.67e-14 | 1.73e-13 | 8.247 |
| IL36A | all samples | GSE186063 | 16 | 2.99e-14 | 3.81e-13 | 8.835 |
| IL36A | all samples | GSE83645 | 16 | 7.38e-14 | 5.84e-13 | 6.928 |
| IL36A | within groups | GSE121212 | 16 | 5.88e-15 | 8.53e-14 | 10.579 |
| IL36A | within groups | GSE186063 | 16 | 4e-15 | 4.69e-13 | 7.879 |
| IL36A | within groups | GSE83645 | 16 | 6.22e-15 | 6.57e-14 | 5.745 |
| IL36G | all samples | GSE121212 | 16 | 3.14e-14 | 8.44e-14 | 9.244 |
| IL36G | all samples | GSE186063 | 16 | 1.83e-14 | 7.11e-14 | 8.322 |
| IL36G | all samples | GSE83645 | 16 | 2.65e-14 | 9.24e-14 | 6.375 |
| GJB2 | all samples | GSE121212 | 16 | 6.49e-15 | 2.49e-14 | 7.464 |
| GJB2 | all samples | GSE186063 | 16 | 3.25e-14 | 1.98e-13 | 9.049 |
| GJB2 | all samples | GSE83645 | 16 | 1.98e-14 | 1.55e-13 | 7.622 |
| SPRR2A | all samples | GSE121212 | 16 | 2.22e-14 | 3.46e-13 | 8.962 |
| SPRR2A | all samples | GSE186063 | 16 | 2.42e-14 | 1.23e-13 | 10.013 |
| SPRR2A | all samples | GSE83645 | 16 | 6.11e-15 | 4.16e-13 | 7.343 |

## 10. Not assessed

- Sample-permutation GSEA p-values, NES and FDR: no package implements label permutation with the application's ranking statistic (moderated t with the trend prior held fixed); fgsea::fgseaLabel ranks genes by correlation with the label, so its p-values are not a parity reference. Enrichment scores of the sample-permutation runs were compared exactly (5 datasets) and gene-permutation p-values statistically against fgseaSimple.
- Treatment x time interaction on GSE171012 whole-skin samples: no valid 2 x 2 exists. Bulk samples have the levels Healthy, Psoriasis_PreTreatment, Psoriasis_SecukinumabTreatmentWeek12, Psoriasis_SecukinumabTreatmentWeek2, Psoriasis_SecukinumabTreatmentWeek4; all psoriasis samples belong to the secukinumab-treated series, healthy skin has no time course, and (Week12 - PreTreatment) - (Week2 - PreTreatment) reduces to the simple contrast Week12 - Week2. A treatment x time interaction was assessed instead on GSE143688 (Aldara vs control, day 7 vs day 3, wild type) with a Cell x Day grouping derived from the GEO annotation.
- Browser execution, rendering and interface behaviour: all application runs used JavaScriptCore with a stub DOM and stub Plotly, not a browser.

## 11. Scope

- Gene-set tests used one library (built-in Hallmark, human or mouse) and the first contrast of seven datasets (GSE143688_all not included); CAMERA was run on the Discovery candidate modules (gene families, HGNC groups and Hallmark sets) with only the two contrasted groups included, as the Discovery CAMERA path requires.
- Meta-analysis used the PSO_SETUP of run_checks.py with a minimum of 3 datasets per gene (19,401 genes); genes measured in only 2 cohorts were not tested.
- Hub preservation: 4 hubs (the 4 most connected genes of the screen: IL36A, IL36G, GJB2, SPRR2A), correlation across all samples for all hubs and within groups for IL36A, 3 comparison datasets, 50 random gene sets per test.
- Welch reference: genes that stats::t.test rejects as constant were assigned t = 0, p = 1 (the application's zero-variance convention); this affected one gene in one run.
- The GSE143688 treatment x time interaction used a derived design row (Cell_Day) built from the GEO Cell and Day annotations of the same samples.

## 12. Files

- `PARITY_FREEZE_checks.csv`: one row per check (check_id, name, dataset, reference, n_compared, max_deviation, tolerance, pass, note).
- `PARITY_FREEZE_scripts.zip`: the pipeline (`run_all.py`, application drivers `app_*.py`, R references `*_ref.R`, comparison modules `cmp_*.py`, report builder).
