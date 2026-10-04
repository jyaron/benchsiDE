# Validation

benchsiDE reimplements its statistical methods in JavaScript. Each method is checked against its reference implementation (edgeR, limma, RRHO, RNASeqPower, metafor, stats, scipy or statsmodels) on inputs that anyone can obtain. Two independent checks are provided, and both run in continuous integration in Chromium, Firefox and WebKit.

1. **In-app self-test** (21 components). A deterministic synthetic dataset is regenerated in the browser and compared with embedded reference values. Any user can run it with **Validate statistics**.
2. **External validation on public data** (70 checks). The public demo dataset in this repository (`demo/`, GEO GSE63310) is analyzed by the application and by R, and the results are compared.

Deviations are the maximum over all compared values. They are absolute unless marked relative. Counts of significant genes and kept genes must agree exactly.

## 1. In-app self-test

Input: 2,000 genes by 6 samples (two groups of three) generated from a seeded pseudo-random generator with integer arithmetic, so the dataset is identical in every JavaScript engine. It contains 200 genes with a planted four-fold change and unequal library sizes. References were computed on the same matrix and are embedded in `index.html`.

| Component | Reference | Maximum deviation | Tolerance |
|---|---|---|---|
| TMM factors vs edgeR | `edgeR::calcNormFactors(m, method="TMM")` | 4.4e-15 | 1.0e-09 |
| moderated t vs limma | `limma::eBayes(lmFit(lg, model.matrix(~grp)))$t` | 6.4e-13 | 1.0e-06 |
| eBayes prior (trend): d0 and median s0² | `limma eBayes(fit, trend=TRUE): df.prior, median(s2.prior)` | 2.1e-14 | 1.0e-06 |
| DE calls (mod) | `sum(p.adjust(p,"BH")<=0.05 & |logFC|>=1) in R` | 0 | 0 (exact) |
| BH FDR | `stats::p.adjust(p, method="BH")` | 1.0e-18 | 1.0e-09 |
| voom t vs limma | `limma::voom(counts, design) + lmFit + eBayes` | 1.0e-12 | 1.0e-06 |
| DE calls (voom) | `R voom pipeline hit counts` | 0 | 0 (exact) |
| CAMERA p vs limma | `limma::camera(lg, index=sets, design)` | 3.2e-14 | 1.0e-08 |
| inter-gene correlation | `limma::interGeneCorrelation via camera()$Correlation` | 7.2e-16 | 1.0e-08 |
| hypergeometric ORA p | `scipy.stats.hypergeom.sf(11, 2000, 30, 40)` | 3.7e-11 | 1.0e-06 |
| OLS fit vs statsmodels | `statsmodels OLS(y, add_constant(x)).fit()` | 6.9e-15 | 1.0e-08 |
| filterByExpr kept set vs edgeR | `edgeR::filterByExpr(m2, group=grp)` | 0 | 5.0e-01 |
| MDS coordinates vs plotMDS | `limma::plotMDS(lg, plot=FALSE)` | 4.5e-10 | 1.0e-08 |
| FRY directional+mixed p vs limma | `limma::fry(lg, index=sets, design, contrast=2)` | 4.2e-11 | 1.0e-05 |
| moderated F vs limma (2-group F=t²) | `limma (fit$t)^2` | 8.4e-14 | 1.0e-07 |
| power model vs RNASeqPower | `RNASeqPower::rnapower()` | 8.6e-10 | 1.0e-08 |
| meta-analysis, random effects (DL) vs metafor | `metafor::rma(yi, vi, method = "DL")` | 3.8e-13 | 1.0e-08 |
| meta-analysis, 95% prediction interval vs metafor | `metafor::predict(rma(yi, vi, method = "DL"), predtype = "Riley")` | 3.0e-13 | 1.0e-08 |
| meta-analysis, REML with Hartung-Knapp vs metafor | `metafor::rma(yi, vi, method = "REML", test = "knha")` | 5.8e-13 | 1.0e-06 |
| meta-analysis, DerSimonian-Laird with Hartung-Knapp vs metafor | `metafor::rma(yi, vi, method = "DL", test = "knha")` | 1.1e-13 | 1.0e-08 |
| meta-analysis, fixed effect vs metafor | `metafor::rma(yi, vi, method = "FE")` | 5.8e-15 | 1.0e-08 |

The self-test panel lists every compared value beside its reference, and **Download evidence (JSON)** saves the record together with the browser's engine string.

## 2. External validation on GSE63310

Input: `demo/GSE63310_counts.tsv` and `demo/GSE63310_design.tsv` (mouse mammary epithelium, 9 samples in three cell types; the dataset of the limma-voom workflow of Law et al., 2016). References: `validation/reference.R`, output `validation/gse63310_reference.json`. Comparison: `validation/compare.js`, run by `tests/validation.spec.js`.

Settings mirrored in R: filterByExpr on all samples (16,624 genes kept); TMM normalization; edgeR log-CPM (`cpm(log = TRUE, prior.count = 2)`); moderated t with `eBayes(trend = TRUE)`; the contrast LP vs Basal fitted on the six samples of those groups. Gene-level statistics are compared at 41 genes (the first 25 and 16 spread across the matrix). The demo design has no covariates, so two test covariates are defined in the reference script to exercise the covariate code: `lcov`, the log2 raw library size (continuous), and `blk`, a two-level block alternating across samples. Gene sets for FRY and over-representation are defined in the script (a seeded random set of 60 genes, the 40 genes with the largest fold change, and a contiguous block of 80 genes).

| Check | Reference | Values compared | Maximum deviation | Tolerance |
|---|---|---|---|---|
| FilterByExpr: kept-gene count and index checksum (exact) | edgeR::filterByExpr | 2 | 0 | 0 (exact) |
| TMM normalization factors | edgeR::calcNormFactors(method = "TMM") | 9 | 4.9e-15 | 1.0e-09 |
| Log-CPM (edgeR, prior.count 2), gene 1, all samples | edgeR::cpm(log = TRUE, prior.count = 2) | 9 | 8.9e-15 | 1.0e-09 |
| Moderated t | limma lmFit, eBayes(trend = TRUE) | 41 | 6.6e-13 | 1.0e-06 |
| Moderated t p-values (relative) | limma eBayes(trend = TRUE), topTable | 41 | 9.9e-13 | 1.0e-06 |
| BH-adjusted p-values | stats::p.adjust(method = "BH") | 41 | 1.4e-13 | 1.0e-08 |
| Empirical-Bayes prior, trend (d0, median s0^2) | limma eBayes(trend = TRUE): df.prior, median(s2.prior) | 2 | 3.0e-13 | 1.0e-06 |
| Moderated t: significant genes up/down (exact) | limma eBayes(trend = TRUE), topTable | 2 | 0 | 0 (exact) |
| Hypergeometric ORA: overlap and hit counts (exact) | stats::phyper | 2 | 0 | 0 (exact) |
| Hypergeometric ORA p (relative) | stats::phyper | 1 | 1.2e-11 | 1.0e-06 |
| Voom: moderated t | limma::voom, lmFit, eBayes | 41 | 2.5e-13 | 1.0e-06 |
| Voom: prior d0 | limma::voom, lmFit, eBayes | 1 | 3.9e-13 | 1.0e-06 |
| Voom: significant genes up/down (exact) | limma::voom, lmFit, eBayes | 2 | 0 | 0 (exact) |
| Continuous covariate: moderated t | limma lmFit(~group + lcov), eBayes(trend = TRUE) | 41 | 2.9e-12 | 1.0e-06 |
| Continuous covariate: prior d0 | limma lmFit(~group + lcov), eBayes(trend = TRUE) | 1 | 2.6e-13 | 1.0e-06 |
| Continuous covariate: significant genes (exact) | limma lmFit(~group + lcov), eBayes(trend = TRUE) | 2 | 0 | 0 (exact) |
| Blocking factor: moderated t | limma lmFit(~group + blk), eBayes(trend = TRUE) | 41 | 1.2e-12 | 1.0e-06 |
| Blocking factor: prior d0 | limma lmFit(~group + blk), eBayes(trend = TRUE) | 1 | 2.1e-13 | 1.0e-06 |
| Blocking factor: significant genes (exact) | limma lmFit(~group + blk), eBayes(trend = TRUE) | 2 | 0 | 0 (exact) |
| Moderated F, three groups (relative) | limma topTable(coef = 2:3) | 41 | 1.8e-13 | 1.0e-07 |
| Moderated F: prior d0 | limma topTable(coef = 2:3) | 1 | 2.4e-13 | 1.0e-06 |
| RemoveBatchEffect (display values) | limma::removeBatchEffect | 90 | 8.0e-15 | 1.0e-09 |
| MDS coordinates and variance explained | limma::plotMDS | 20 | 1.1e-14 | 1.0e-08 |
| FRY without covariate: p and mixed p (relative) | limma::fry(design = ~lcov + group) | 6 | 2.4e-11 | 1.0e-05 |
| FRY without covariate: direction (exact) | limma::fry(design = ~lcov + group) | 3 | 0 | 0 (exact) |
| FRY with continuous covariate: p and mixed p (relative) | limma::fry(design = ~lcov + group) | 6 | 2.5e-11 | 1.0e-05 |
| FRY with continuous covariate: direction (exact) | limma::fry(design = ~lcov + group) | 3 | 0 | 0 (exact) |
| OLS regression: intercept, slope, SE, p | stats::lm | 4 | 2.3e-14 | 1.0e-08 |
| Co-expression: Pearson r | stats::cor | 40 | 1.3e-14 | 1.0e-10 |
| Co-expression: slope | stats::lm | 40 | 4.2e-15 | 1.0e-10 |
| Co-expression: p for r = 0 (relative) | stats::cor.test | 40 | 5.5e-13 | 1.0e-06 |
| Hub neighbourhood: partial correlation given group | stats::lm residuals, cor | 40 | 1.1e-14 | 1.0e-10 |
| Hub neighbourhood: p for partial r = 0 (relative) | stats::lm (gene ~ gene 1 + group) | 40 | 8.7e-13 | 1.0e-06 |
| Connectivity: sum of |r|^6 | stats::cor | 200 | 2.9e-14 | 1.0e-09 |
| Module preservation: 7 statistics, observed and 5 random sets (WGCNA 1.74) | WGCNA 1.74 .coreCalcForExpr | 42 | 1.0e-14 | 1.0e-10 |
| Robust rank aggregation: RRA score (relative) | RobustRankAggreg::aggregateRanks | 300 | 5.0e-15 | 1.0e-09 |
| Trait association: moderated t (limma-trend, ~ trait) | limma lmFit + eBayes(trend = TRUE) | 50 | 1.6e-14 | 1.0e-08 |
| Trait association: p (relative) | limma | 50 | 6.7e-13 | 1.0e-06 |
| Trait association with blocks: duplicateCorrelation consensus | limma::duplicateCorrelation | 1 | 9.4e-16 | 1.0e-08 |
| Trait association with blocks: moderated t (lmFit block, correlation) | limma lmFit(block, correlation) | 50 | 1.1e-13 | 1.0e-06 |
| Trait gene sets: CAMERA p (relative) | limma::camera | 3 | 6.3e-14 | 1.0e-06 |
| Module eigengene and its trait t (WGCNA moduleEigengenes, lm) | WGCNA 1.74 moduleEigengenes, stats::lm | 10 | 3.3e-14 | 1.0e-08 |
| Density: nrd0 bandwidth | stats::bw.nrd0; exact Gaussian kernel sum | 1 | 9.4e-16 | 1.0e-12 |
| Density: curve values at 9 grid points | stats::bw.nrd0; exact Gaussian kernel sum | 9 | 1.7e-10 | 1.0e-08 |
| RRHO grid dimensions (exact) | RRHO::RRHO(alternative = "enrichment") | 2 | 0 | 0 (exact) |
| RRHO -log10 p, diagonal and first row | RRHO::RRHO(alternative = "enrichment") | 80 | 9.6e-11 | 1.0e-08 |
| Full design: moderated t | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 41 | 1.5e-12 | 1.0e-06 |
| Full design: moderated t p-values (relative) | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 41 | 1.2e-12 | 1.0e-06 |
| Full design: prior d0 and residual df | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 2 | 2.4e-13 | 1.0e-06 |
| Full design: significant genes up/down (exact) | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 2 | 0 | 0 (exact) |
| Full design, continuous covariate: moderated t | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 41 | 1.5e-12 | 1.0e-06 |
| Full design, continuous covariate: significant genes (exact) | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 2 | 0 | 0 (exact) |
| Full design, blocking factor: moderated t | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 41 | 2.5e-12 | 1.0e-06 |
| Full design, blocking factor: significant genes (exact) | `limma::lmFit` + `eBayes(trend = TRUE)` (design `~ group`, all samples, B − A coefficient) | 2 | 0 | 0 (exact) |
| Full design, voom: moderated t | `limma::voom`, `lmFit`, `eBayes` (design `~ group`, all samples) | 41 | 1.1e-13 | 1.0e-06 |
| Full design, voom: significant genes up/down (exact) | `limma::voom`, `lmFit`, `eBayes` (design `~ group`, all samples) | 2 | 0 | 0 (exact) |
| FRY full design: p and mixed p (relative) | `limma::fry` (design and sample order as in the application) | 6 | 4.4e-11 | 1.0e-05 |
| FRY full design: direction (exact) | `limma::fry` (design and sample order as in the application) | 3 | 0 | 0 (exact) |
| Power model: n per group and power (relative) | RNASeqPower::rnapower | 2 | 1.4e-09 | 1.0e-08 |

The checks without the prefix "full design" use the option *variance from: these two groups only* (the model fitted to the two contrasted groups); the "full design" checks use the default, in which every group is fitted and the two groups are compared by a contrast. FRY's robust standardization depends on the basis of the residual space and therefore on sample order, so the R reference for the full-design FRY check uses the application's sample order (group A, group B, then the remaining groups).

Reference versions: R 4.5.3, edgeR 4.8.2, limma 3.66.0, RRHO 1.50.0, RNASeqPower 1.50.0; metafor 5.0.1 for the meta-analysis checks.

## Meta-analysis on public data

The meta-analysis was checked on two public psoriasis datasets, using the count matrices generated by NCBI for GEO: GSE186063 (lesional vs non-lesional skin) and GSE54456 (psoriatic vs normal skin). Each dataset was filtered (filterByExpr), TMM-normalized and tested (limma-trend moderated t on log-CPM) within itself. For each of the 19,308 genes matched in both, the application's pooled estimate was compared with `metafor::rma(yi, vi)` given the same per-dataset log2 fold changes and squared standard errors.

| Quantity | Random effects (DL): maximum deviation | Fixed effect: maximum deviation |
|---|---|---|
| Pooled log2FC (absolute) | 5.3e-15 | 8.9e-15 |
| Standard error (relative) | 1.2e-14 | 8.2e-16 |
| p-value (relative) | 1.2e-12 | 9.1e-13 |
| Cochran's Q (relative) | 6.1e-15 | 6.1e-15 |
| Q p-value (absolute) | 4.0e-14 | 4.0e-14 |
| I² (percentage points) | 3.6e-13 | 3.6e-13 |
| tau² (absolute) | 1.2e-14 | not applicable |

metafor 5.0.1. The six-gene synthetic case in the self-test (table 1) covers 3 to 5 datasets per gene.


## Full-matrix check of limma-trend on public data

GSE54456 (NCBI-generated counts; 171 samples, psoriatic vs normal skin), all 19,518 genes kept by filterByExpr, compared with `eBayes(lmFit(cpm(y, log = TRUE, prior.count = 2), ~group), trend = TRUE)` in R.

| Quantity | Maximum deviation |
|---|---|
| Prior df (app 4.848480, R 4.848480) | 0.0e+00 |
| log2FC (absolute) | 1.4e-14 |
| Moderated t (absolute) | 2.6e-13 |
| p-value (relative) | 6.8e-12 |
| Significant genes, FDR ≤ 0.05 and \|log2FC\| ≥ 1 | 1293 up, 1513 down in both |

## GSEA with sample permutation

Enrichment scores for the observed ranking match `fgsea::calcGseaStat(gseaParam = 1)` (fgsea 1.36.2) to 2.8e-16 for 10 Hallmark sets. NES, nominal p and FDR recomputed independently in Python from the permutation null of one run: maximum deviation 0, 0 and 0. Calibration on null data is reported in `validation/freeze/CALIBRATION_FREEZE.md`.

## Calibration and known-biology controls

The parity checks above show that each statistic matches its reference implementation. They cannot show that a default or threshold is appropriate. `validation/freeze/CALIBRATION_FREEZE.md` and `validation/study/VAL_CAL.md` report the controls used for that question: random splits of normal skin (no true difference) to measure false positives, and public psoriasis data to confirm that well-established changes are recovered.

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

## Performance changes in 0.23.0-beta

Three computations were reimplemented for speed. Each was checked against the previous implementation before release:

| Change | Check | Result |
|---|---|---|
| MDS leading-fold-change distance: top-500 selection by quickselect instead of a full sort | MDS coordinates, previous vs new, on GSE54456 (171 samples), a mouse time course (18 samples) and an 8-sample two-group set | bit-identical |
| FRY: eigenvalues by Householder tridiagonalization and implicit QL instead of Jacobi rotation | eigenvalues vs numpy.linalg.eigvalsh on 24 symmetric, Gram and rank-deficient matrices (1 to 170 rows) | max relative deviation 6e-15 |
| FRY: Gram matrices on typed arrays | FRY p and mixed p, previous vs new, Hallmark on GSE54456 and on the 8-sample set | max relative difference 1.2e-10; directional p identical |
| GSEA permutation loop moved into a separate function | GSEA ES, NES and p, previous vs new | identical |

The external validation against R (table above) gives the same deviations as before.

## Validation study on public psoriasis data

A validation study of 3,403 checks on six public datasets, including null-data calibration, is reported in `validation/study/`.

## Hub across datasets: calibration on public psoriasis data

Discovery dataset GSE54456 (psoriasis vs normal skin, 171 samples); comparison datasets GSE121212 (PSO lesional vs non-lesional, 54), GSE186063 (lesion vs non-lesion, 54) and GSE83645 (psoriasis vs uninvolved, 25). Neighbourhoods of 15 partners among the 2,000 most variable genes. Parity of the statistics with WGCNA 1.74, RobustRankAggreg and scipy is in the external validation table above and in `validation/study/hub_across_datasets.md`.

| Check | Expected under the null | Result | Consequence |
|---|---|---|---|
| Preservation of real neighbourhoods with gene correspondence shuffled (25 hubs) | Zsummary ≈ 0 | median -0.53 (all samples), -0.63 (within groups); none > 2 | the Z statistic is centred |
| Preservation of random 16-gene sets (40) | — | median Zsummary 3.3; 100% above the published Z = 2 threshold (57% within groups) | published thresholds not used; empirical baselines instead |
| Preservation of random sets of up-regulated genes against random sets (25 × 3 datasets) | — | 87% "preserved" (all samples; 76%–92% by dataset), 51% within groups (8%–84%) | DE-matched baseline added; the verdict requires preservation beyond it |
| Random up-regulated sets against the DE-matched baseline (20 × 3 datasets per mode) | 5% at p ≤ 0.05 | 3.3% (all samples), 3.3% (within groups) | calibrated |
| Random gene sets against the random-set baseline (30) | 5% at p ≤ 0.05 | 3% | calibrated |
| Fisher z on random splits of healthy skin, random gene pairs | 5% at p < 0.05 | 6 v 6: 4.9%; 10 v 10: 5.0%; 40 v 40: 4.9% | the test is calibrated on non-DE genes |
| Fisher z on random pairs of up-regulated genes | 5% if nothing rewires | GSE186063 40% (non-DE pairs 3%); GSE121212 9% (non-DE 6%) | raw Fisher z not used for calls; DE-matched baseline added |
| DE-matched rewiring p on random up-regulated pairs | 5% at p < 0.05 | A (this session): 4.7% (raw 39%), GSE121212: 5.0% (raw 7%), GSE186063: 4.1% (raw 38%) | calibrated |
| Consensus hubs with every dataset's co-expression destroyed (per-gene sample permutation, 2 runs) | 0 at FDR ≤ 0.05 | 0, 0 (real data: 145) | calibrated |
| Consensus hubs with only the discovery dataset intact (2 runs) | 0 | 0, 0 | one dataset cannot produce a consensus call |

Ten strongest hubs, three comparison datasets each (30 tests): preserved beyond DE-matched sets in 15 of 30 with correlation across all samples and in 30 of 30 within groups.

## Sample traits: parity and calibration

Parity with limma 3.66 on GSE186063 (Age; lesion and non-lesion skin; patient pairs as blocks), all genes: duplicateCorrelation consensus correlation to 6.7e-16; moderated t to 4.6e-14 and p to 5.9e-10 (relative) for three models (block means with group removed within blocks; ~ Age + group + Sex; ~ rank(Age)); block means computed independently in R agree to 7.1e-15; CAMERA on the trait coefficient (50 Hallmark sets) p to 3.7e-8 (relative) with and without the mixed-model transformation; module eigengenes against WGCNA moduleEigengenes to 4.6e-13 and their trait t against lm to 1e-10.

| Check | Expected | Result | Consequence |
|---|---|---|---|
| Random subject-level trait, GSE83645 lesional skin (20 biopsies from 5 patients), 40 replicates, samples treated as independent | 5% of genes at p < 0.05 | 25.0%; any gene at FDR ≤ 0.05 in 90% of replicates | repeated samples must be declared |
| Same, mixed model (duplicateCorrelation) | 5% | 8.4%; any FDR hit in 70% | not used for subject-level traits: a consensus correlation is too low for many genes with only 5 subjects |
| Same, block means (default for subject-level traits) | 5% | 6.3%; any FDR hit in 5% | calibrated |
| Random sample-level trait, GSE83645 lesional, 40 replicates: independent / mixed model | 5% | 4.7% / 5.0%; any FDR hit 5% / 2% | calibrated |
| GSE186063 Age shuffled across patients (lesion + non-lesion, group-adjusted), 40 replicates: independent / mixed model / block means | any FDR hit in ≤ 5% (FDR control under the complete null) | 18% / 8% / 2% | block means is the default |
| Random trait, GSE54456 (171 samples, group-adjusted), 25 replicates: genes / CAMERA (Hallmark) | 5% at p < 0.05 | 4.2% / 4.0% | calibrated |
| Random trait, module eigengene test (300 replicates; random 20-gene modules / hub neighbourhoods) | 5% | 6.0% / 6.3% | calibrated |
| Planted trait effects in healthy skin (GSE54456), 200 genes, 5 replicates each | FDP ≤ 5%; unbiased slopes | n = 10, 0.6 log₂/SD: 41% recovered, FDP 1.2%; n = 20, 0.4 log₂/SD: 57% recovered, FDP 2.4%; n = 40, 0.3 log₂/SD: 65% recovered, FDP 0.4%; n = 81, 0.2 log₂/SD: 66% recovered, FDP 0.6% | mean slope bias < 0.005 log₂/SD |

## Release candidate 1.0.0-rc.1: full check on the frozen build

The build was frozen (application script sha256 `cadbda80bd67db0e…` before the fixes listed below) and checked by three independent tracks in the JavaScriptCore test engine (`tests/engine/`). The reports and every check are in `validation/freeze/`.

| Track | Checks | Report |
|---|---|---|
| Parity with R (limma 3.66.0, edgeR 4.8.2, metafor 5.0.1, fgsea 1.36.2, WGCNA 1.74 functions) on GSE63310, GSE54456, GSE121212, GSE186063, GSE41745, GSE83645 and GSE143688 | 1,096 | PARITY_FREEZE.md |
| Null calibration and binomial-thinning simulation (seqgendiff 1.2.4) on four skin pools and one non-skin pool (GSE57945 non-IBD ileum) | 499 | CALIBRATION_FREEZE.md |
| Regression against the 0.23.0 baseline; every tab, export, report and session round trip | 227 | REGRESSION_FREEZE.md |

**Parity.** Filtering, TMM and log-CPM agreed in all 16 sessions (log-CPM 6.2e-15). Gene-level DE (moderated t, voom, Welch; all-groups and two-groups fits; categorical, continuous, patient and pair covariates) agreed to 3.4e-12 in t and 6.7e-11 (relative) in p, with exact significant-gene counts in 77 of 78 runs; the moderated F to 1.4e-11; all 12 interaction fits to 4.6e-13 in t with exact counts. ORA, CAMERA and GSEA enrichment scores agreed to 6.9e-11, 1.5e-11 and 4.6e-14. The meta-analysis matched metafor for REML, DerSimonian-Laird and fixed effect, with z, Hartung-Knapp and truncated Hartung-Knapp tests (pooled log₂FC 4.9e-10, prediction intervals 6.1e-9, counts exact). Hub preservation statistics matched WGCNA's own functions in 15 hub × dataset combinations (5.8e-13). With 100 design columns, FRY p-values depend on the parametrisation of the design, in limma as in the application.

**Null calibration.** The moderated t exceeded the binomial bound only at 3 samples per group in the four skin pools (14–21 of 100 splits with a discovery; the application warns at this size) and stayed within it in the non-skin pool and at 5 and 10 per group everywhere. FRY, CAMERA, sample-permutation GSEA and the Discovery screen were within the bound (one FRY cell 12/100 against a bound of 11). ORA of a top-200 list exceeded it in every cell (25–73 of 100), as stated in the application. The default meta-analysis (REML with the Hartung-Knapp test) stayed within the bound in all 12 configurations of 2–5 datasets at 3, 5 and 10 per group; the z test exceeded it with 4–5 datasets of 3 per group, and the fixed-effect model in every configuration with 3 per group.

**Simulation (true effects inserted into real null counts).** Gene-level tests controlled the FDR (largest mean false-discovery proportion 0.079, moderated t at 3 v 3); power at |log₂FC| ≥ 1 was 0.28, 0.69 and 0.91 at 3, 5 and 10 per group in skin. With effects shared by all datasets, every meta-analysis variant controlled the FDR. When effects differed between datasets (τ = 0.4), the observed FDR was 0.080–0.222 for REML or DerSimonian-Laird with the z test and 0.176–0.298 for the fixed-effect model; with the Hartung-Knapp test it was not significantly above 0.05 (largest mean false-discovery proportion 0.132, from very few calls), and power was at most 0.030. Its power, however, was 0 with 2 datasets, at most 0.001 with 3, 0.05–0.37 with 4 and 0.26–0.57 with 5, against 0.29–0.75 for REML with the z test. REML with Hartung-Knapp is kept as the default because it is the only variant that controlled the FDR in every condition; the application states its power with 2–4 datasets at the point of use and offers the z test as an exploratory alternative.

**Regression.** run_checks.py passed, and the snapshot reproduced 46 of 47 baseline fields bit for bit (the 47th is the version string). With the meta-analysis set to DerSimonian-Laird and the z test, every baseline field was reproduced exactly. All tabs, 185 plot exports, 84 CSV/GMT exports, the report, the methods text and 56 session comparisons ran without error. Both builds carry the Content-Security-Policy, and the subresource-integrity hash matches the Plotly file.

Corrections made after the check are covered by the engine check `rc_fixes`.

**Not assessed here.** Real-browser rendering and Content-Security-Policy enforcement (covered by `tests/privacy.spec.js` and `tests/freeze.spec.js` in CI on Chromium, Firefox and WebKit); exported image pixels; paired nulls, co-expression FDR, RRHO and signature transfer on the frozen build (code paths unchanged since VAL_CAL).

## Release candidate 1.0.0-rc.2: covariates with voom, and TREAT

voom with covariates and TREAT were checked against limma 3.66.0 in 80 configurations on GSE63310 and GSE186063 (two methods, up to two covariates including a blocking factor and a continuous covariate, both fit options, ordinary and TREAT tests). The maximum deviations were 1.4 × 10⁻¹² in log₂FC, 1.1 × 10⁻¹¹ in t and 1.4 × 10⁻¹⁰ (relative) in p, and the number of significant genes was identical in every configuration. The analysis of the limma/Glimma/edgeR workflow article on GSE63310 (voom adjusted for sequencing lane, TREAT with |log₂FC| > 1) gives the limma values exactly: 3,647, 3,831 and 2,782 genes for LP vs Basal, ML vs Basal and both. On 300 null splits of healthy skin with a random block covariate, the tests behaved as documented above for the moderated t and voom without covariates (above the bound at three samples per group, within it from five), and TREAT stayed within the bound at every size; every count was reproduced by limma in R. Details: `validation/rc2/RC2_DE.md`.
