# Validation

benchsiDE reimplements its statistical methods in JavaScript. Each method is checked against its reference implementation (edgeR, limma, RRHO, RNASeqPower, stats, scipy or statsmodels) on inputs that anyone can obtain. Two independent checks are provided, and both run in continuous integration in Chromium, Firefox and WebKit.

1. **In-app self-test** (16 components). A deterministic synthetic dataset is regenerated in the browser and compared with embedded reference values. Any user can run it with **Validate statistics**.
2. **External validation on public data** (33 checks). The public demo dataset in this repository (`demo/`, GEO GSE63310) is analyzed by the application and by R, and the results are compared.

Deviations are the maximum over all compared values. They are absolute unless marked relative. Counts of significant genes and kept genes must agree exactly.

## 1. In-app self-test

Input: 2,000 genes by 6 samples (two groups of three) generated from a seeded pseudo-random generator with integer arithmetic, so the dataset is identical in every JavaScript engine. It contains 200 genes with a planted four-fold change and unequal library sizes. References were computed on the same matrix and are embedded in `index.html`.

| Component | Reference | Maximum deviation | Tolerance |
|---|---|---|---|
| TMM factors vs edgeR | `edgeR::calcNormFactors(m, method="TMM")` | 4.4e-15 | 1.0e-9 |
| moderated t vs limma | `limma::eBayes(lmFit(lg, model.matrix(~grp)))$t` | 3.6e-8 | 1.0e-6 |
| eBayes prior d0 | `limma fit$df.prior` | 1.0e-8 | 1.0e-6 |
| DE calls (mod) | `sum(p.adjust(p,"BH")<=0.05 & |logFC|>=1) in R` | 0 | 0 (exact) |
| BH FDR | `stats::p.adjust(p, method="BH")` | 3.7e-11 | 1.0e-9 |
| voom t vs limma | `limma::voom(counts, design) + lmFit + eBayes` | 4.8e-12 | 1.0e-6 |
| DE calls (voom) | `R voom pipeline hit counts` | 0 | 0 (exact) |
| CAMERA p vs limma | `limma::camera(lg, index=sets, design)` | 3.7e-10 | 1.0e-8 |
| inter-gene correlation | `limma::interGeneCorrelation via camera()$Correlation` | 2.4e-14 | 1.0e-8 |
| hypergeometric ORA p | `scipy.stats.hypergeom.sf(11, 2000, 30, 40)` | 3.7e-11 | 1.0e-6 |
| OLS fit vs statsmodels | `statsmodels OLS(y, add_constant(x)).fit()` | 4.9e-13 | 1.0e-8 |
| filterByExpr kept set vs edgeR | `edgeR::filterByExpr(m2, group=grp)` | 0 | 5.0e-1 |
| MDS coordinates vs plotMDS | `limma::plotMDS(lg, plot=FALSE)` | 3.3e-12 | 1.0e-8 |
| FRY directional+mixed p vs limma | `limma::fry(lg, index=sets, design, contrast=2)` | 7.9e-9 | 1.0e-5 |
| moderated F vs limma (2-group F=t²) | `limma (fit$t)^2` | 3.6e-9 | 1.0e-7 |
| power model vs RNASeqPower | `RNASeqPower::rnapower()` | 8.6e-10 | 1.0e-8 |

The self-test panel lists every compared value beside its reference, and **Download evidence (JSON)** saves the record together with the browser's engine string.

## 2. External validation on GSE63310

Input: `demo/GSE63310_counts.tsv` and `demo/GSE63310_design.tsv` (mouse mammary epithelium, 9 samples in three cell types; the dataset of the limma-voom workflow of Law et al., 2016). References: `validation/reference.R`, output `validation/gse63310_reference.json`. Comparison: `validation/compare.js`, run by `tests/validation.spec.js`.

Settings mirrored in R: filterByExpr on all samples (16,624 genes kept); TMM normalization; log2(CPM + 1) with effective library sizes; the contrast LP vs Basal fitted on the six samples of those groups. Gene-level statistics are compared at 41 genes (the first 25 and 16 spread across the matrix). The demo design has no covariates, so two test covariates are defined in the reference script to exercise the covariate code: `lcov`, the log2 raw library size (continuous), and `blk`, a two-level block alternating across samples. Gene sets for FRY and over-representation are defined in the script (a seeded random set of 60 genes, the 40 genes with the largest fold change, and a contiguous block of 80 genes).

| Check | Reference | Values compared | Maximum deviation | Tolerance |
|---|---|---|---|---|
| FilterByExpr: kept-gene count and index checksum (exact) | edgeR::filterByExpr | 2 | 0 | 0 (exact) |
| TMM normalization factors | edgeR::calcNormFactors(method = "TMM") | 9 | 4.9e-15 | 1.0e-9 |
| Log2(CPM + 1), gene 1, all samples | edgeR TMM effective library sizes | 9 | 4.0e-15 | 1.0e-9 |
| Moderated t | limma lmFit, eBayes | 41 | 4.0e-8 | 1.0e-6 |
| Moderated t p-values (relative) | limma eBayes, topTable | 41 | 2.2e-8 | 1.0e-6 |
| BH-adjusted p-values | stats::p.adjust(method = "BH") | 41 | 7.5e-10 | 1.0e-8 |
| Empirical-Bayes prior (d0, s0^2) | limma fit$df.prior, fit$s2.prior | 2 | 1.7e-8 | 1.0e-6 |
| Moderated t: significant genes up/down (exact) | limma eBayes, topTable | 2 | 0 | 0 (exact) |
| Hypergeometric ORA: overlap and hit counts (exact) | stats::phyper | 2 | 0 | 0 (exact) |
| Hypergeometric ORA p (relative) | stats::phyper | 1 | 1.6e-12 | 1.0e-6 |
| Voom: moderated t | limma::voom, lmFit, eBayes | 41 | 3.3e-8 | 1.0e-6 |
| Voom: prior d0 | limma::voom, lmFit, eBayes | 1 | 4.4e-8 | 1.0e-6 |
| Voom: significant genes up/down (exact) | limma::voom, lmFit, eBayes | 2 | 0 | 0 (exact) |
| Continuous covariate: moderated t | limma lmFit(~group + lcov), eBayes | 41 | 1.2e-8 | 1.0e-6 |
| Continuous covariate: prior d0 | limma lmFit(~group + lcov), eBayes | 1 | 9.6e-9 | 1.0e-6 |
| Continuous covariate: significant genes (exact) | limma lmFit(~group + lcov), eBayes | 2 | 0 | 0 (exact) |
| Blocking factor: moderated t | limma lmFit(~group + blk), eBayes | 41 | 2.1e-8 | 1.0e-6 |
| Blocking factor: prior d0 | limma lmFit(~group + blk), eBayes | 1 | 1.0e-8 | 1.0e-6 |
| Blocking factor: significant genes (exact) | limma lmFit(~group + blk), eBayes | 2 | 0 | 0 (exact) |
| Moderated F, three groups (relative) | limma topTable(coef = 2:3) | 41 | 4.7e-9 | 1.0e-7 |
| Moderated F: prior d0 | limma topTable(coef = 2:3) | 1 | 1.4e-8 | 1.0e-6 |
| RemoveBatchEffect (display values) | limma::removeBatchEffect | 90 | 5.3e-15 | 1.0e-9 |
| MDS coordinates and variance explained | limma::plotMDS | 20 | 1.6e-14 | 1.0e-8 |
| FRY without covariate: p and mixed p (relative) | limma::fry(design = ~lcov + group) | 6 | 5.4e-9 | 1.0e-5 |
| FRY without covariate: direction (exact) | limma::fry(design = ~lcov + group) | 3 | 0 | 0 (exact) |
| FRY with continuous covariate: p and mixed p (relative) | limma::fry(design = ~lcov + group) | 6 | 4.8e-9 | 1.0e-5 |
| FRY with continuous covariate: direction (exact) | limma::fry(design = ~lcov + group) | 3 | 0 | 0 (exact) |
| OLS regression: intercept, slope, SE, p | stats::lm | 4 | 6.2e-16 | 1.0e-8 |
| Density: nrd0 bandwidth | stats::bw.nrd0; exact Gaussian kernel sum | 1 | 8.9e-16 | 1.0e-12 |
| Density: curve values at 9 grid points | stats::bw.nrd0; exact Gaussian kernel sum | 9 | 2.3e-10 | 1.0e-8 |
| RRHO grid dimensions (exact) | RRHO::RRHO(alternative = "enrichment") | 2 | 0 | 0 (exact) |
| RRHO -log10 p, diagonal and first row | RRHO::RRHO(alternative = "enrichment") | 79 | 1.0e-10 | 1.0e-8 |
| Power model: n per group and power (relative) | RNASeqPower::rnapower | 2 | 1.4e-9 | 1.0e-8 |

Reference versions: R 4.5.3, edgeR 4.8.2, limma 3.66.0, RRHO 1.50.0, RNASeqPower 1.50.0.

## Reproducing the results

```bash
Rscript validation/reference.R      # regenerate the R reference values (optional)
npm install
npx playwright install --with-deps
npm test                            # self-test, demo pins and external validation in three browsers
```

## Conventions reproduced deliberately

- **limma voom span.** limma 3.66 and later select the lowess span adaptively (`adaptive.span = TRUE`), overriding `span = 0.5`. The application uses the same rule.
- **Zero-count genes in voom.** Genes with zero counts in every sample of the contrast are excluded from the mean-variance trend fit, as in limma.
- **TMM column sums.** TMM trims genes by rank, so its result can change when floating-point ties change. Library sizes are summed sequentially in plain double precision, matching R's `colSums`. A more accurate summation algorithm would move results away from edgeR's.
- **RRHO population size.** The RRHO package computes its hypergeometric probabilities with a population of N + 1. The application reproduces this for numerical agreement. Under-enrichment, which the package does not report, is shown as negative values.

## Checks without a numerical reference

Display and interface behavior is checked by test scripts, not against an external package: every plot is drawn on two-group and three-group designs at all three font scales and with 50-character sample and group names, and its layout is checked for margin and overlap violations. Session files are saved and restored, and the restored state is compared with the original.

## Scope

Agreement is shown for the inputs and settings listed above. The implementations are general, but agreement on other inputs is not established by these checks. The self-test can be run on any installation, and a discrepancy can be reported with its evidence file (see `CONTRIBUTING.md`).
