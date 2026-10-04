# VAL_DE: end-to-end concordance of the benchsiDE DE tab with R

Build under test: `benchside/index.html`, APP_VERSION 0.22.0-beta (2026-09-30), SHA-256 17511676c3096bda2f4ba964dea8ce5b78faa300a29ad0d1007c210da7c7ae65 (matches MANIFEST.json). The build was not modified. Count matrices were checked against the SHA-256 values in MANIFEST.json (all six match).
References: R version 4.5.3 (2026-03-11); edgeR 4.8.2; limma 3.66.0; orthogonal methods: R version 4.5.3 (2026-03-11); edgeR 4.8.2; DESeq2 1.50.2. The application was run headlessly in macOS JavaScriptCore with the kit's DOM/Plotly stub (`harness.py`); no browser was used.

## 1. Summary

- 1,481 checks were run (1,259 gene-level and display parity checks, 210 methods-text checks, 4 selection-table checks, 8 diagnostic and unit checks). 151 did not pass: 121 methods-text checks, 14 parity checks in one voom contrast, 11 Welch checks whose R reference excludes genes that `t.test()` cannot test (a convention difference, section 4.4; all pass when the app's convention is applied in R), and 5 diagnostic or unit checks.
- Every gene-level statistic of the moderated t (limma-trend), covariate-adjusted and blocked fits, voom, Welch t and the moderated F agrees with R on all kept genes in all six datasets: maximum deviations are ≤ 3.4×10⁻¹² for log2FC and t (absolute) and ≤ 6.7×10⁻¹¹ for p and FDR (relative). Kept gene sets, TMM factors, the full log-CPM matrices, MDS and PCA coordinates agree to ≤ 3.1×10⁻¹². Significant up/down counts agree exactly in all 42 runs at all six thresholds (FDR 0.05 and 0.01 × |log2FC| ≥ 0, 1, 2), except three Welch counts explained by the convention in 4.4. No p-value below 10⁻³⁰⁰ occurred in any run.
- One numerical defect changes results on the supplied data. One numerical defect is reachable only when the prior degrees of freedom are infinite. The generated methods text misreports what was computed in several configurations, most seriously by stating a covariate adjustment that was not applied when voom or Welch is selected.
- Methodological concerns: pairwise contrasts in multi-group data are fitted on the two groups only, and numeric-looking identifiers are fitted as continuous slopes.
- Orthogonal agreement: under the pre-stated rule (section 6) the app's default method does not disagree markedly with both edgeR QL and DESeq2 in any contrast. The lowest agreement for the default method is GSE41745 paired by Patient (n = 5, one residual df): Jaccard 0.688 with QL and 0.488 with DESeq2; QL and DESeq2 agree with each other at 0.671.

## 2. What was compared

For each dataset (factor in parentheses) the app was loaded with `species`, built-in annotation, `normsel = tmm`, and `fmode = fbe`; the primary contrast was re-run with `fmode = cpm` (CPM ≥ 1 in ≥ 2 samples). R: `DGEList` → `filterByExpr(y, group = factor)` (or `rowSums(cpm(y) >= 1) >= 2`) → `y[keep, , keep.lib.sizes = FALSE]` → `calcNormFactors(TMM)` → `cpm(log = TRUE, prior.count = 2)` on all samples. DE fits in R mirror the app: the samples of the two contrasted groups (A = reference, B = test), minus samples with a missing value in a selected covariate.

- moderated t: `lmFit(logCPM[, A+B], ~ group [+ covariates])`, `eBayes(trend = TRUE)`, `topTable(coef = 2, sort.by = "none")`; numeric covariates with more than two distinct values as numeric (Age), others as factors (Diag, Sex, Pair, Patient, Site).
- voom: `voom(y[, A+B], ~ group)` (limma default `adaptive.span = TRUE`), `lmFit`, `eBayes()`; group means as the weighted group coefficients of `~ 0 + group`.
- Welch: `stats::t.test(var.equal = FALSE)` per gene on log-CPM.
- moderated F: all samples of groups with n ≥ 2, `~ group`, `eBayes(trend = TRUE)`, `topTable(coef = 2:k)`.
- MDS: `plotMDS(logCPM, top = 500, plot = FALSE)`. PCA: `prcomp(t(logCPM[top2000, ]))`, top 2,000 genes by variance. p-value histogram: 20 bins, R counts both as `floor(20p)` (the app's [a, b) rule) and with `hist()` (right-closed).

Contrasts (42 runs): GSE63310 (celltype) LP vs Basal, ML vs Basal, ML vs LP × {moderated t, voom, Welch} + moderated F over 3 groups; GSE54456 (Cond) Psoriasis_skin vs normal_skin × 3 methods; GSE186063 (Type) lesion vs non-lesion × 3 methods, + Diag, + Sex + Age (Age continuous), + Pair (2 samples with Pair NA dropped, 52 samples), voom with Pair selected, and moderated F over Type; GSE41745 (Cond) lesional vs non_lesional (2 v 3) × 3 methods and + Patient; GSE83645 (Cond) psoriasis vs uninvolved × 3 methods, + Patient, + Site, + Patient + Site; GSE121212 (Group) PSO_lesional vs PSO_non_lesional × 3 methods and + Patient, PSO_lesional vs CTRL_healthy × 3 methods, and moderated F over the 6 Group levels; plus the `fmode = cpm` run of each dataset's first contrast (moderated t).

Compared per run, on all kept genes: log2FC, group means A and B, t, p, BH FDR, SE of log2FC, per-gene posterior variance s2.post (moderated t runs: derived from the app's SE and R's unscaled variance), prior d0, median s0², residual and total df, up/down counts at six thresholds (read from the rendered count pills), the 20 p-value histogram bins as plotted, the volcano trace values, the top-50 DE table as rendered (gene, log2FC, p, FDR strings), the CSV export (columns, rows, values), the eBayes prior line, and the methods text.

Tolerances: 1×10⁻⁸ absolute for log-scale statistics, 1×10⁻⁶ relative for p, FDR, variances and d0, 0 for sets and counts; rendered or exported values at their printed precision (3 decimals: 5×10⁻⁴; 4 significant digits: 5×10⁻⁴ relative).

## 3. Session-level outputs

| Dataset (filter) | Kept genes app / R | Kept-set differences | TMM max dev | log-CPM max dev | log-CPM values | MDS coord max dev | MDS var. expl. max dev | PCA top-2000 set differences | PCA PC1-2 max dev | PCA axis % dev |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE63310 (fbe) | 16624 / 16624 | 0 | 2.2e-16 | 5.3e-15 | 149616 | 1.1e-14 | 4.0e-16 | 0 | 1.8e-13 | 0 |
| GSE63310 (cpm) | 14490 / 14490 | 0 | 2.2e-16 | 5.8e-15 | 130410 | – | – | – | – | – |
| GSE54456 (fbe) | 19518 / 19518 | 0 | 4.4e-16 | 5.3e-15 | 3337578 | 8.3e-14 | 1.9e-15 | 0 | 3.1e-12 | 0 |
| GSE54456 (cpm) | 19902 / 19902 | 0 | 2.2e-16 | 5.8e-15 | 3403242 | – | – | – | – | – |
| GSE186063 (fbe) | 21112 / 21112 | 0 | 2.2e-16 | 5.3e-15 | 1372280 | 3.6e-14 | 5.1e-16 | 0 | 1.1e-12 | 0 |
| GSE186063 (cpm) | 20505 / 20505 | 0 | 4.4e-16 | 5.3e-15 | 1332825 | – | – | – | – | – |
| GSE41745 (fbe) | 18036 / 18036 | 0 | 2.2e-16 | 5.3e-15 | 90180 | 1.8e-15 | 8.3e-16 | 0 | 1.5e-13 | 0 |
| GSE41745 (cpm) | 17409 / 17409 | 0 | 1.1e-16 | 3.6e-15 | 87045 | – | – | – | – | – |
| GSE83645 (fbe) | 20956 / 20956 | 0 | 2.2e-16 | 5.3e-15 | 523900 | 1.2e-14 | 1.9e-16 | 0 | 4.0e-13 | 0 |
| GSE83645 (cpm) | 18591 / 18591 | 0 | 2.2e-16 | 5.3e-15 | 464775 | – | – | – | – | – |
| GSE121212 (fbe) | 23168 / 23168 | 0 | 3.3e-16 | 6.2e-15 | 3336192 | 6.0e-14 | 2.4e-15 | 0 | 2.6e-12 | 0 |
| GSE121212 (cpm) | 21755 / 21755 | 0 | 3.3e-16 | 5.3e-15 | 3132720 | – | – | – | – | – |

All kept-gene sets are identical in count and identity. The TMM factors, the full log-CPM matrices (all genes × all samples), the MDS coordinates (sign-aligned) and PCA scores agree to machine precision. The plotted MDS points equal the computed coordinates (deviation 0), and the PCA axis titles give the prcomp variance explained to the printed 1 decimal.

## 4. Differential expression, gene by gene

### 4.1 Two-group runs

Maximum deviations over all kept genes (absolute unless marked relative):

| Run (contrast · method · covariates, filter) | Genes | log2FC | t | p (rel) | FDR (rel) | Group means | d0 (app) | d0 dev (rel) | median s0² dev (rel) | Up/down FDR≤0.05, abs(log2FC)≥1: app vs R | Failed checks | Checks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE63310 LP_vs_Basal · mod (fbe) | 16624 | 8.9e-15 | 1.1e-12 | 1.6e-12 | 1.6e-12 | 7.1e-15 | 5.128372 | 5.9e-14 | 3.4e-14 | app 2780/3357, R 2780/3357 | 0 | 28 |
| GSE63310 LP_vs_Basal · voom (fbe) | 16624 | 6.1e-14 | 9.3e-13 | 1.4e-12 | 1.4e-12 | 6.6e-14 | 5.480555 | 6.8e-14 | 2.8e-14 | app 2810/3338, R 2810/3338 | 0 | 27 |
| GSE63310 LP_vs_Basal · welch (fbe) | 16624 | 8.0e-15 | 2.0e-12 | 9.2e-13 | 9.2e-13 | 7.1e-15 | – | – | – | app 1968/2418, R 1968/2418 | 0 | 23 |
| GSE63310 ML_vs_Basal · mod (fbe) | 16624 | 8.0e-15 | 1.0e-12 | 1.2e-12 | 1.2e-12 | 7.1e-15 | 4.411452 | 6.0e-14 | 2.4e-14 | app 2969/3315, R 2969/3315 | 0 | 28 |
| GSE63310 ML_vs_Basal · voom (fbe) | 16624 | 7.6e-14 | 1.7e-12 | 1.7e-12 | 1.7e-12 | 7.8e-14 | 4.727329 | 6.0e-14 | 2.4e-14 | app 2966/3301, R 2966/3301 | 0 | 27 |
| GSE63310 ML_vs_Basal · welch (fbe) | 16624 | 7.1e-15 | 2.8e-12 | 9.7e-13 | 6.0e-05 | 7.1e-15 | – | – | – | app 2282/2422, R 2282/2422 | 7 | 32 |
| GSE63310 ML_vs_LP · mod (fbe) | 16624 | 5.6e-15 | 1.4e-13 | 2.8e-12 | 1.1e-12 | 7.1e-15 | 3.046528 | 2.8e-14 | 2.1e-15 | app 1299/1189, R 1299/1189 | 0 | 28 |
| GSE63310 ML_vs_LP · voom (fbe) | 16624 | 2.8e-05 | 1.4e-04 | 1.5e-04 | 1.2e-04 | 3.6e-05 | 3.578226 | 4.6e-07 | 2.0e-06 | app 1287/1212, R 1287/1212 | 14 | 27 |
| GSE63310 ML_vs_LP · welch (fbe) | 16624 | 7.3e-15 | 1.2e-12 | 2.2e-12 | 9.0e-04 | 7.1e-15 | – | – | – | app 0/0, R 0/0 | 4 | 32 |
| GSE63310 LP_vs_Basal · mod (cpm) | 14490 | 8.0e-15 | 1.2e-12 | 1.6e-12 | 1.2e-12 | 5.3e-15 | 4.914532 | 6.5e-14 | 4.2e-14 | app 2404/2835, R 2404/2835 | 0 | 28 |
| GSE54456 Psoriasis_vs_normal · mod (fbe) | 19518 | 1.0e-14 | 2.7e-13 | 6.8e-12 | 5.4e-12 | 7.1e-15 | 4.848483 | 2.3e-14 | 3.8e-15 | app 1293/1513, R 1293/1513 | 0 | 28 |
| GSE54456 Psoriasis_vs_normal · voom (fbe) | 19518 | 2.7e-14 | 8.0e-13 | 7.5e-12 | 7.5e-12 | 5.0e-14 | 4.651368 | 2.4e-14 | 9.8e-15 | app 1340/1653, R 1340/1653 | 0 | 27 |
| GSE54456 Psoriasis_vs_normal · welch (fbe) | 19518 | 7.1e-15 | 2.8e-13 | 5.4e-12 | 5.4e-12 | 7.1e-15 | – | – | – | app 1293/1513, R 1293/1513 | 0 | 23 |
| GSE54456 Psoriasis_vs_normal · mod (cpm) | 19902 | 1.0e-14 | 2.9e-13 | 4.0e-12 | 4.0e-12 | 6.2e-15 | 4.860823 | 2.5e-14 | 5.6e-15 | app 1344/1696, R 1344/1696 | 0 | 28 |
| GSE186063 lesion_vs_nonlesion · mod (fbe) | 21112 | 8.3e-15 | 9.2e-14 | 8.7e-12 | 8.7e-12 | 7.1e-15 | 4.361219 | 2.3e-14 | 1.9e-14 | app 1015/1399, R 1015/1399 | 0 | 28 |
| GSE186063 lesion_vs_nonlesion · mod · Diag (fbe) | 21112 | 8.2e-15 | 1.0e-13 | 8.7e-12 | 8.7e-12 | 7.1e-15 | 4.373029 | 2.3e-14 | 1.2e-14 | app 1015/1407, R 1015/1407 | 0 | 28 |
| GSE186063 lesion_vs_nonlesion · mod · Pair (fbe) | 21112 | 6.3e-15 | 9.1e-14 | 4.3e-12 | 4.3e-12 | 7.1e-15 | 4.413064 | 2.7e-14 | 1.3e-14 | app 997/1429, R 997/1429 | 0 | 28 |
| GSE186063 lesion_vs_nonlesion · mod · Sex+Age (fbe) | 21112 | 1.0e-14 | 1.2e-13 | 8.6e-12 | 8.6e-12 | 7.1e-15 | 4.444961 | 2.7e-14 | 1.4e-14 | app 1015/1404, R 1015/1404 | 0 | 28 |
| GSE186063 lesion_vs_nonlesion · voom (fbe) | 21112 | 5.2e-14 | 2.1e-13 | 6.7e-11 | 6.7e-11 | 4.8e-14 | 4.194359 | 2.9e-14 | 1.3e-14 | app 1114/1708, R 1114/1708 | 0 | 27 |
| GSE186063 lesion_vs_nonlesion · voom · Pair (fbe) | 21112 | 5.2e-14 | 2.1e-13 | 6.7e-11 | 6.7e-11 | 4.8e-14 | 4.194359 | 2.9e-14 | 1.3e-14 | app 1114/1708, R 1114/1708 | 0 | 27 |
| GSE186063 lesion_vs_nonlesion · welch (fbe) | 21112 | 7.2e-15 | 1.0e-13 | 8.4e-12 | 8.4e-12 | 7.1e-15 | – | – | – | app 1015/1398, R 1015/1398 | 0 | 23 |
| GSE186063 lesion_vs_nonlesion · mod (cpm) | 20505 | 8.4e-15 | 9.1e-14 | 8.5e-12 | 8.5e-12 | 7.1e-15 | 4.326684 | 3.2e-14 | 9.1e-15 | app 1026/1386, R 1026/1386 | 0 | 28 |
| GSE41745 lesional_vs_nonlesional · mod (fbe) | 18036 | 5.0e-15 | 3.0e-13 | 2.9e-12 | 2.5e-12 | 3.6e-15 | 7.339579 | 8.5e-14 | 1.4e-14 | app 634/301, R 634/301 | 0 | 28 |
| GSE41745 lesional_vs_nonlesional · mod · Patient (fbe) | 18036 | 6.7e-15 | 3.4e-12 | 2.0e-11 | 2.0e-11 | 3.6e-15 | 16.376311 | 4.0e-13 | 8.8e-14 | app 867/690, R 867/690 | 0 | 28 |
| GSE41745 lesional_vs_nonlesional · voom (fbe) | 18036 | 1.1e-13 | 2.7e-13 | 3.5e-11 | 3.5e-11 | 9.1e-14 | 7.519045 | 6.2e-14 | 1.8e-14 | app 627/182, R 627/182 | 0 | 27 |
| GSE41745 lesional_vs_nonlesional · welch (fbe) | 18036 | 5.4e-15 | 1.6e-12 | 1.6e-12 | 7.5e-13 | 3.6e-15 | – | – | – | app 0/0, R 0/0 | 0 | 23 |
| GSE41745 lesional_vs_nonlesional · mod (cpm) | 17409 | 6.2e-15 | 2.6e-13 | 1.9e-12 | 1.9e-12 | 5.3e-15 | 7.318816 | 6.7e-14 | 2.5e-14 | app 636/306, R 636/306 | 0 | 28 |
| GSE83645 psoriasis_vs_uninvolved · mod (fbe) | 20956 | 9.8e-15 | 7.0e-14 | 7.5e-12 | 7.5e-12 | 7.1e-15 | 2.930815 | 1.9e-14 | 9.9e-16 | app 932/919, R 932/919 | 0 | 28 |
| GSE83645 psoriasis_vs_uninvolved · mod · Patient (fbe) | 20956 | 9.3e-15 | 9.6e-14 | 3.6e-12 | 3.6e-12 | 7.1e-15 | 3.018758 | 2.2e-14 | 5.6e-15 | app 1026/1120, R 1026/1120 | 0 | 28 |
| GSE83645 psoriasis_vs_uninvolved · mod · Patient+Site (fbe) | 20956 | 3.4e-14 | 3.5e-13 | 2.6e-12 | 2.6e-12 | 7.1e-15 | 3.056074 | 3.3e-14 | 2.4e-14 | app 767/602, R 767/602 | 0 | 28 |
| GSE83645 psoriasis_vs_uninvolved · mod · Site (fbe) | 20956 | 3.5e-14 | 2.7e-13 | 8.4e-12 | 8.4e-12 | 7.1e-15 | 2.93962 | 2.9e-14 | 1.2e-14 | app 600/280, R 600/280 | 0 | 28 |
| GSE83645 psoriasis_vs_uninvolved · voom (fbe) | 20956 | 5.4e-14 | 5.1e-13 | 3.1e-12 | 3.1e-12 | 4.1e-14 | 2.899666 | 1.2e-14 | 5.9e-15 | app 873/1041, R 873/1041 | 0 | 27 |
| GSE83645 psoriasis_vs_uninvolved · welch (fbe) | 20956 | 5.4e-15 | 9.9e-14 | 1.8e-12 | 1.8e-12 | 7.1e-15 | – | – | – | app 818/855, R 818/855 | 0 | 23 |
| GSE83645 psoriasis_vs_uninvolved · mod (cpm) | 18591 | 8.1e-15 | 1.2e-13 | 3.0e-12 | 3.0e-12 | 7.1e-15 | 2.794949 | 1.2e-14 | 1.2e-14 | app 859/748, R 859/748 | 0 | 28 |
| GSE121212 PSOL_vs_CTRL · mod (fbe) | 23168 | 6.6e-15 | 1.3e-13 | 8.5e-12 | 8.5e-12 | 5.3e-15 | 5.381442 | 1.4e-14 | 1.0e-14 | app 1751/2640, R 1751/2640 | 0 | 28 |
| GSE121212 PSOL_vs_CTRL · voom (fbe) | 23168 | 6.0e-14 | 8.1e-13 | 8.5e-12 | 8.5e-12 | 5.9e-14 | 5.120704 | 2.3e-14 | 7.6e-15 | app 1967/3424, R 1967/3424 | 0 | 27 |
| GSE121212 PSOL_vs_CTRL · welch (fbe) | 23168 | 8.9e-15 | 1.3e-13 | 8.6e-12 | 8.6e-12 | 5.3e-15 | – | – | – | app 1750/2640, R 1750/2640 | 0 | 23 |
| GSE121212 PSOL_vs_PSONL · mod (fbe) | 23168 | 6.7e-15 | 1.4e-13 | 8.5e-12 | 8.5e-12 | 5.3e-15 | 4.786261 | 2.3e-14 | 1.8e-14 | app 1401/1846, R 1401/1846 | 0 | 28 |
| GSE121212 PSOL_vs_PSONL · mod · Patient (fbe) | 23168 | 1.1e-14 | 2.7e-13 | 4.4e-12 | 4.4e-12 | 5.3e-15 | 4.487269 | 2.7e-14 | 1.5e-14 | app 1420/1825, R 1420/1825 | 0 | 28 |
| GSE121212 PSOL_vs_PSONL · voom (fbe) | 23168 | 5.9e-14 | 3.2e-13 | 1.0e-11 | 8.5e-12 | 4.1e-14 | 4.667625 | 2.5e-14 | 8.3e-15 | app 1664/2499, R 1664/2499 | 0 | 27 |
| GSE121212 PSOL_vs_PSONL · welch (fbe) | 23168 | 7.2e-15 | 1.5e-13 | 8.4e-12 | 8.4e-12 | 5.3e-15 | – | – | – | app 1401/1846, R 1401/1846 | 0 | 23 |
| GSE121212 PSOL_vs_PSONL · mod (cpm) | 21755 | 7.3e-15 | 1.3e-13 | 8.6e-12 | 8.6e-12 | 7.1e-15 | 4.301113 | 2.2e-14 | 1.7e-14 | app 1293/1891, R 1293/1891 | 0 | 28 |

The Welch runs GSE63310 ML_vs_Basal and ML_vs_LP contain genes that R's `t.test()` rejects (section 4.4). All other runs pass every check.

Calls at all six thresholds: 252 comparisons, 3 differences, all three in GSE63310 ML_vs_Basal | Welch at FDR ≤ 0.01 (section 4.4). The per-threshold counts are in `VAL_DE_checks.csv`.

### 4.2 Covariates, pairing and missing values

- GSE186063 + Pair: the app leaves out the 2 non-lesion samples with Pair = NA and fits 52 samples with 26 pair levels (residual df 25), exactly as `lmFit(logCPM[, complete], ~ Type + factor(Pair))`; d0 4.413064 in both; calls 997/1429 in both. The prior line does not report that 2 samples were dropped.
- GSE186063 + Sex + Age: Age is fitted as a slope (prior line: "Sex (categorical) + Age (continuous)", residual df 50), matching `~ Type + Sex + Age` with numeric Age.
- GSE41745 + Patient: Patient1 has only a non-lesional sample; the app, like limma, keeps it with its own Patient coefficient (residual df 1, d0 16.376311 in both, calls 867/690 in both). With one residual df the result rests almost entirely on the prior.
- GSE83645 + Patient, + Site, + Patient + Site and GSE121212 + Patient: all pass (table 4.1).
- voom and Welch ignore selected covariates. The DE panel says so ("Covariate adjustment ... applies only to the moderated t method; current ... results are UNADJUSTED"), and the GSE186063 voom run with Pair selected is identical to the unadjusted voom run; the methods text, however, states the adjustment.

### 4.3 Moderated F

| Test | Genes | F (rel) | p (rel) | FDR (rel) | d0 (rel) | median s0² (rel) | FDR≤0.05 genes | FDR≤0.01 genes |
|---|---|---|---|---|---|---|---|---|
| GSE63310 moderated F over 3 groups | 16624 | 1.9e-12 | 8.8e-13 | 9.7e-13 | 5.7e-14 | 1.1e-14 | app 11335, R 11335 | app 8789, R 8789 |
| GSE186063 moderated F over 3 groups | 21112 | 8.2e-13 | 9.2e-12 | 9.2e-12 | 2.6e-14 | 1.5e-14 | app 11824, R 11824 | app 9337, R 9337 |
| GSE121212 moderated F over 6 groups | 23168 | 8.8e-14 | 7.7e-12 | 7.7e-12 | 2.4e-14 | 1.0e-14 | app 20107, R 20107 | app 18408, R 18408 |

Groups used: GSE63310 Basal, LP, ML (9 samples, residual df 6); GSE186063 lesion, non-lesion, normal_skin (65, df 62); GSE121212 all six Group levels including the atopic-dermatitis groups (144, df 138), as in `topTable(coef = 2:k)` on `~ Group`.

### 4.4 Welch t on genes with constant log-CPM

A gene with zero counts in every sample of the contrast has a constant log-CPM (the library-size-scaled prior count cancels), so both group variances are 0. The app returns t = 0, p = 1 and keeps the gene in the BH family; `t.test()` stops with "data are essentially constant". This affects 1 gene in GSE63310 ML vs Basal and 15 genes in ML vs LP. With p = 1 imputed for these genes in R, FDR agrees to 9.5×10⁻¹³ and 9.3×10⁻¹³ (relative), and calls and histogram bins agree exactly. Without imputation the BH family is larger by 1 or 15 genes; in ML vs Basal this moves a BH plateau from 0.00999959 to 0.01000019 and changes 7 up and 7 down calls at FDR ≤ 0.01 (app 812/728, R 819/735 at |log2FC| ≥ 1). This is a convention difference, not classified as a defect, but the methods text does not state it.

### 4.5 Displays and exports

- p-value histogram: the plotted 20 bins equal `floor(20p)` counts of the R p-values in all runs and equal `hist()` counts (no p fell on an inner bin edge), apart from the constant-gene Welch cases above.
- Volcano: plotted x equals log2FC and y equals −log10 p exactly (deviation 0) in all 42 runs.
- Top-50 DE table (gene, log2FC to 2 decimals, p and FDR as formatted by `fmtP`): identical to the R values sorted by adjusted p in all 42 runs.
- CSV export (`DE_<B>_vs_<A>.csv`): columns gene, gene_id, biotype, mean_log2_A, mean_log2_B, log2FC, p, FDR; one row per kept gene in all runs; values agree with R at the printed precision.
- Selection table and its CSV (GSE83645 + Patient, all 20,956 genes selected): 300 rendered rows and all CSV values agree with R at the printed precision; the `mean_log2` column is not limma's AveExpr.

## 5. Methods text

| Methods-text check | Configurations | Passed |
|---|---|---|
| methods text: covariate adjustment stated iff applied | 42 | 41 |
| methods text: no sentence for analyses not run (ORA, ANOVA/Kruskal-Wallis, k-means) | 42 | 0 |
| methods text: sample count stated = samples in DE fit | 42 | 16 |
| methods text: significance thresholds (FDR, abs(log2FC)) stated | 42 | 0 |
| methods text: transform stated matches DE input | 42 | 32 |

Configuration-level results are in `VAL_DE_checks.csv` (rows named "methods text: ...").

## 6. Orthogonal agreement (not parity)

edgeR `estimateDisp` + `glmQLFit` + `glmQLFTest`, and DESeq2 `DESeq` + `results` (default independent filtering and Cook's cutoff, unshrunk log2FC), on the same filterByExpr genes, sample subsets and designs as the app's moderated-t runs (Age scaled for DESeq2). Statistics: app t vs QL sign(logFC)·√F and vs DESeq2 Wald statistic (Spearman ρ). Significant sets: FDR ≤ 0.05 and |log2FC| ≥ 1, each method with its own log2FC and FDR (J = Jaccard); FDR-only Jaccard also shown. Rule fixed before inspection: the app default (moderated t) is flagged if J < 0.5 against both QL and DESeq2, or ρ < 0.8 against both.

| Dataset | Contrast | Covariates | App method | App sig | QL sig | DESeq2 sig | ρ vs QL | ρ vs DESeq2 | J vs QL | J vs DESeq2 | J QL–DESeq2 | J vs QL (FDR only) | J vs DESeq2 (FDR only) | Flag |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE63310 | LP_vs_Basal | none | mod | 6137 | 6212 | 6021 | 0.999 | 0.997 | 0.961 | 0.946 | 0.955 | 0.961 | 0.94 | False |
| GSE63310 | LP_vs_Basal | none | voom | 6148 | 6212 | 6021 | 0.999 | 0.997 | 0.963 | 0.943 | 0.955 | 0.959 | 0.941 | False |
| GSE63310 | LP_vs_Basal | none | welch | 4386 | 6212 | 6021 | 0.988 | 0.985 | 0.705 | 0.715 | 0.955 | 0.606 | 0.602 | False |
| GSE63310 | ML_vs_Basal | none | mod | 6284 | 6366 | 6186 | 0.999 | 0.997 | 0.959 | 0.945 | 0.946 | 0.963 | 0.939 | False |
| GSE63310 | ML_vs_Basal | none | voom | 6267 | 6366 | 6186 | 0.998 | 0.997 | 0.955 | 0.935 | 0.946 | 0.957 | 0.934 | False |
| GSE63310 | ML_vs_Basal | none | welch | 4704 | 6366 | 6186 | 0.989 | 0.986 | 0.736 | 0.747 | 0.946 | 0.643 | 0.637 | False |
| GSE63310 | ML_vs_LP | none | mod | 2488 | 2608 | 2697 | 0.996 | 0.995 | 0.905 | 0.814 | 0.822 | 0.922 | 0.8 | False |
| GSE63310 | ML_vs_LP | none | voom | 2499 | 2608 | 2697 | 0.997 | 0.993 | 0.889 | 0.814 | 0.822 | 0.916 | 0.809 | False |
| GSE63310 | ML_vs_LP | none | welch | 0 | 2608 | 2697 | 0.985 | 0.982 | 0 | 0.0 | 0.822 | 0 | 0.0 | False |
| GSE54456 | Psoriasis_vs_normal | none | mod | 2806 | 2941 | 2994 | 0.998 | 0.997 | 0.885 | 0.894 | 0.932 | 0.955 | 0.943 | False |
| GSE54456 | Psoriasis_vs_normal | none | voom | 2993 | 2941 | 2994 | 0.997 | 0.997 | 0.899 | 0.923 | 0.932 | 0.956 | 0.94 | False |
| GSE54456 | Psoriasis_vs_normal | none | welch | 2806 | 2941 | 2994 | 0.998 | 0.997 | 0.885 | 0.894 | 0.932 | 0.955 | 0.943 | False |
| GSE186063 | lesion_vs_nonlesion | none | mod | 2414 | 2772 | 2774 | 0.992 | 0.993 | 0.763 | 0.769 | 0.944 | 0.894 | 0.895 | False |
| GSE186063 | lesion_vs_nonlesion | none | voom | 2822 | 2772 | 2774 | 0.99 | 0.991 | 0.795 | 0.802 | 0.944 | 0.894 | 0.892 | False |
| GSE186063 | lesion_vs_nonlesion | none | welch | 2413 | 2772 | 2774 | 0.992 | 0.993 | 0.763 | 0.769 | 0.944 | 0.891 | 0.895 | False |
| GSE186063 | lesion_vs_nonlesion | Diag | mod | 2422 | 2732 | 2775 | 0.993 | 0.994 | 0.774 | 0.771 | 0.944 | 0.898 | 0.899 | False |
| GSE186063 | lesion_vs_nonlesion | Sex+Age | mod | 2419 | 2745 | 2708 | 0.993 | 0.994 | 0.775 | 0.771 | 0.932 | 0.903 | 0.904 | False |
| GSE186063 | lesion_vs_nonlesion | Pair | mod | 2426 | 2689 | 2818 | 0.998 | 0.998 | 0.864 | 0.845 | 0.924 | 0.962 | 0.963 | False |
| GSE41745 | lesional_vs_nonlesional | none | mod | 935 | 769 | 1280 | 0.996 | 0.995 | 0.796 | 0.695 | 0.581 | 0.769 | 0.601 | False |
| GSE41745 | lesional_vs_nonlesional | none | voom | 809 | 769 | 1280 | 0.994 | 0.994 | 0.822 | 0.598 | 0.581 | 0.807 | 0.516 | False |
| GSE41745 | lesional_vs_nonlesional | none | welch | 0 | 769 | 1280 | 0.972 | 0.968 | 0 | 0.0 | 0.581 | 0 | 0.0 | False |
| GSE41745 | lesional_vs_nonlesional | Patient | mod | 1557 | 1098 | 779 | 0.997 | 0.998 | 0.688 | 0.488 | 0.671 | 0.541 | 0.35 | False |
| GSE83645 | psoriasis_vs_uninvolved | none | mod | 1851 | 1869 | 1919 | 0.99 | 0.99 | 0.782 | 0.756 | 0.892 | 0.876 | 0.858 | False |
| GSE83645 | psoriasis_vs_uninvolved | none | voom | 1914 | 1869 | 1919 | 0.988 | 0.989 | 0.754 | 0.717 | 0.892 | 0.853 | 0.818 | False |
| GSE83645 | psoriasis_vs_uninvolved | none | welch | 1673 | 1869 | 1919 | 0.976 | 0.977 | 0.661 | 0.641 | 0.892 | 0.558 | 0.543 | False |
| GSE83645 | psoriasis_vs_uninvolved | Patient | mod | 2146 | 2271 | 2304 | 0.996 | 0.997 | 0.85 | 0.856 | 0.932 | 0.931 | 0.924 | False |
| GSE83645 | psoriasis_vs_uninvolved | Site | mod | 880 | 947 | 1123 | 0.99 | 0.99 | 0.827 | 0.742 | 0.805 | 0.87 | 0.7 | False |
| GSE83645 | psoriasis_vs_uninvolved | Patient+Site | mod | 1369 | 1556 | 1610 | 0.995 | 0.995 | 0.801 | 0.81 | 0.888 | 0.896 | 0.818 | False |
| GSE121212 | PSOL_vs_PSONL | none | mod | 3247 | 4226 | 4216 | 0.996 | 0.995 | 0.742 | 0.73 | 0.907 | 0.936 | 0.9 | False |
| GSE121212 | PSOL_vs_PSONL | none | voom | 4163 | 4226 | 4216 | 0.995 | 0.994 | 0.87 | 0.849 | 0.907 | 0.943 | 0.894 | False |
| GSE121212 | PSOL_vs_PSONL | none | welch | 3247 | 4226 | 4216 | 0.995 | 0.995 | 0.742 | 0.73 | 0.907 | 0.935 | 0.899 | False |
| GSE121212 | PSOL_vs_PSONL | Patient | mod | 3245 | 4153 | 4301 | 0.997 | 0.995 | 0.772 | 0.727 | 0.876 | 0.958 | 0.911 | False |
| GSE121212 | PSOL_vs_CTRL | none | mod | 4391 | 5466 | 5517 | 0.997 | 0.996 | 0.781 | 0.76 | 0.914 | 0.946 | 0.908 | False |
| GSE121212 | PSOL_vs_CTRL | none | voom | 5391 | 5466 | 5517 | 0.996 | 0.995 | 0.903 | 0.879 | 0.914 | 0.947 | 0.9 | False |
| GSE121212 | PSOL_vs_CTRL | none | welch | 4390 | 5466 | 5517 | 0.997 | 0.996 | 0.781 | 0.76 | 0.914 | 0.945 | 0.906 | False |

No contrast is flagged. Rank agreement of the default method is ρ ≥ 0.9897 against QL and ≥ 0.9903 against DESeq2 in every contrast. The set agreement is lower than the rank agreement mainly through the fold-change cut-off: in GSE121212 PSO_lesional vs PSO_non_lesional, 886 of the 4,014 genes called by both QL and DESeq2 are not called by the app default; 831 of these have app FDR ≤ 0.05 but |log2FC| < 1 (median |log2FC| 0.809 in the app vs 1.236 in QL; median of (mean_A + mean_B)/2 −1.95 log-CPM against 2.49 for all genes). The corresponding numbers are 998 of 5,246 (969 with FDR ≤ 0.05) for PSO_lesional vs CTRL_healthy and 461 of 2,693 (414) for GSE186063. The prior count of 2 in log-CPM moderates the fold changes of low-count genes, so the same |log2FC| ≥ 1 cut-off is more stringent than on GLM fold changes. This is how limma-trend on edgeR log-CPM behaves; it is reported as a property of the default, not as a defect. voom, which models counts, is closer to QL and DESeq2 in these contrasts (J 0.870 and 0.849 for PSO_lesional vs PSO_non_lesional). Welch calls nothing at FDR ≤ 0.05 for GSE63310 ML vs LP and GSE41745 (n = 3 v 3 and 2 v 3).

## 8. Methodological context for multi-group fits

| Dataset | Run | Samples in two-group fit | Residual df | d0 | Up/down FDR≤0.05 (two-group vs full) | Up/down FDR≤0.01 (two-group vs full) |
|---|---|---|---|---|---|---|
| GSE63310 | LP_vs_Basal · mod | 6 | 4 vs 6 | 5.128 vs 4.196 | 2780/3357 vs 2763/3305 | 2442/3002 vs 2446/2945 |
| GSE63310 | LP_vs_Basal · voom | 6 | 4 vs 6 | 5.481 vs 4.419 | 2810/3338 vs 2781/3301 | 2464/2950 vs 2419/2921 |
| GSE63310 | ML_vs_Basal · mod | 6 | 4 vs 6 | 4.411 vs 4.196 | 2969/3315 vs 2962/3307 | 2671/2954 vs 2704/2965 |
| GSE63310 | ML_vs_Basal · voom | 6 | 4 vs 6 | 4.727 vs 4.419 | 2966/3301 vs 2973/3329 | 2651/2893 vs 2699/2970 |
| GSE63310 | ML_vs_LP · mod | 6 | 4 vs 6 | 3.047 vs 4.196 | 1299/1189 vs 1418/1387 | 663/445 vs 993/854 |
| GSE63310 | ML_vs_LP · voom | 6 | 4 vs 6 | 3.578 vs 4.419 | 1287/1212 vs 1331/1290 | 760/528 vs 931/771 |
| GSE186063 | lesion_vs_nonlesion · mod | 54 | 52 vs 62 | 4.361 vs 4.290 | 1015/1399 vs 1016/1399 | 1008/1378 vs 1006/1367 |
| GSE186063 | lesion_vs_nonlesion · voom | 54 | 52 vs 62 | 4.194 vs 4.136 | 1114/1708 vs 1115/1703 | 1088/1631 vs 1089/1618 |
| GSE121212 | PSOL_vs_PSONL · mod | 54 | 52 vs 138 | 4.786 vs 5.009 | 1401/1846 vs 1401/1845 | 1394/1841 vs 1399/1839 |
| GSE121212 | PSOL_vs_PSONL · voom | 54 | 52 vs 138 | 4.668 vs 4.808 | 1664/2499 vs 1657/2481 | 1646/2444 vs 1639/2394 |
| GSE121212 | PSOL_vs_CTRL · mod | 64 | 62 vs 138 | 5.381 vs 5.009 | 1751/2640 vs 1750/2640 | 1749/2637 vs 1749/2635 |
| GSE121212 | PSOL_vs_CTRL · voom | 64 | 62 vs 138 | 5.121 vs 4.808 | 1967/3424 vs 1966/3419 | 1961/3397 vs 1962/3391 |

## 9. Feature status (DE tab)

| Feature | Status | Basis |
|---|---|---|
| filterByExpr filter | validated | Kept sets identical in count and identity on 6 datasets (16,624–23,168 genes). |
| CPM filter (fmode cpm) | validated | Kept sets identical on 6 datasets; downstream moderated t passes. |
| TMM normalization | validated | Factors within 4.4e-16 of calcNormFactors on 12 dataset/filter combinations. |
| log-CPM (prior count 2) | validated | 17,360,763 values compared; max deviation 6.2e-15. |
| Moderated t, limma-trend, two groups | validated_with_caveat | Parity on 15 runs without covariates (t <= 1.3e-12, p <= 8.7e-12 rel., calls exact). Caveats: two-group-only fit in multi-group data; normal approximation branch when d0 is infinite. |
| Covariate adjustment and fixed-effect blocking | validated_with_caveat | Parity on 8 runs incl. Pair with NA drop, continuous Age, Patient with an unpaired patient. |
| Missing covariate values | validated | Pair NA: 2 samples dropped, identical to limma on the complete cases. |
| voom | defect_found | 9 of 10 runs pass (log2FC <= 1.1e-13, t <= 1.8e-12, p <= 6.8e-11 rel.); covariates not supported in this version. |
| Welch t | validated_with_caveat | Parity on 9 runs (t <= 2.9e-12, p <= 8.6e-12 rel.); constant genes get p = 1 and stay in the BH family (convention, section 4.4). |
| BH FDR | validated | Relative deviation <= 6.8e-11 in all runs except one voom run (1.2e-4) and the Welch convention cases (section 4.4). |
| Moderated F | validated | 3 designs (3, 3, 6 groups); F, p, FDR <= 9.2e-12 rel.; call counts exact. |
| Significance calls (count pills) | validated | 252 threshold comparisons exact (3 Welch differences explained by the convention). |
| eBayes prior line | validated_with_caveat | d0, s0², df values correct; does not report samples dropped for missing covariate values. |
| MDS (plotMDS) | validated | 6 datasets, coordinates <= 8.3e-14, variance explained <= 2.4e-15. |
| PCA | validated | 6 datasets, top-2000 sets identical, scores <= 3.1e-12, axis percentages correct. |
| p-value histogram | validated | Bins exact in all runs (Welch constant-gene convention aside). |
| Volcano plot values | validated | Plotted x, y exactly equal to log2FC and -log10 p in 42 runs. |
| MA plot | defect_found | x axis is mean of group means, labeled mean log2 expression. |
| Top-50 DE table | validated | Rendered rows identical to R ordering and formatting in 42 runs. |
| DE CSV export | defect_found | Values correct at printed precision. |
| DE selection table and CSV | validated_with_caveat | Values correct; mean_log2 column is not AveExpr. |
| Methods text for DE | revised | Corrected in later versions. |

## 10. Not assessed

- Behaviour in real browsers (all runs used JavaScriptCore with the kit stub; in the stub, textContent does not unescape &amp;, so entity rendering in the methods text was not assessed).
- Per-gene s2.prior for voom and moderated F (d0 and median s0² compared; per-gene s2.post compared for moderated-t runs only).
- Voom or Welch with covariates (not implemented in the app; only the reporting of this was checked).
- Random-effect blocking (limma duplicateCorrelation); the app offers fixed-effect blocking only.
- p-value QQ plot, hit-count stability plot and voom mean-variance trend plot values.
- The d0 = infinity path on real data (did not occur; the affected p-value functions were tested on unit inputs only).
- All-pairs contrast panel, contrast-vs-contrast panel, sample exclusion via the included toggle, normalization other than TMM, log-scale and already-normalized input modes.
- DESeq2 on the unfiltered matrix (orthogonal runs used the filterByExpr genes).

## 11. Reproduction

Scripts (in `ci/`): `de_configs.json` (datasets and runs, shared by app driver and R), `app_de_driver.py` (runs the app through `harness.py` and dumps every DE output), `reference_DE.R` (R reference, writes binary log-CPM and %.17g tables), `compare_DE.py` (gene-by-gene comparison, writes the checks CSV, exits non-zero on any failure except the documented Welch convention rows), `ortho_DE.R` and `ortho_compare.py` (edgeR QL and DESeq2 agreement), `diag_voom_span.R`, `diag_numeric_block.R`, `run_VAL_DE.sh` (runs the pipeline). Inputs: the validation kit directory and a data directory holding `GSE63310_counts.tsv` and the five `<GSE>_raw_counts_GRCh38.p13_NCBI.tsv` files.

```bash
bash ci/run_VAL_DE.sh benchside/index.html kit data out
```

Outputs: `VAL_DE_checks.csv` (every check: name, dataset, reference, n_compared, max_deviation, tolerance, pass, note) and `VAL_DE_orthogonal.csv`.
