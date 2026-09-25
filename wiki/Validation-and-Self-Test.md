# Validation and self-test

## How benchsiDE is validated

Each statistical routine is compared with its reference implementation on inputs that anyone can obtain. The results are recorded in [VALIDATION.md](../VALIDATION.md), with the maximum deviation observed for each check. Two checks are provided:

1. **In-app self-test.** Sixteen components on a deterministic synthetic dataset (below).
2. **External validation on public data.** The demo dataset in the repository (GEO GSE63310) is analyzed by the application and by R. `validation/reference.R` produces the reference values, and `tests/validation.spec.js` compares them with the application's results (33 checks).

| Component | Reference |
|---|---|
| TMM normalization | edgeR `calcNormFactors` |
| filterByExpr | edgeR `filterByExpr` |
| Moderated t, moderated F, covariate models | limma `lmFit`, `eBayes` |
| voom | limma `voom` (3.66) |
| CAMERA, FRY | limma `camera`, `fry` |
| MDS | limma `plotMDS` |
| Batch removal (display) | limma `removeBatchEffect` |
| RRHO | Bioconductor RRHO |
| Hypergeometric test, regression, density | R `stats`, scipy, statsmodels |
| Power model | RNASeqPower |

## In-browser self-test

**Validate statistics** runs the suite in your browser:

1. A deterministic synthetic dataset is regenerated from a seeded pseudo-random generator.
2. The dataset is processed by the same functions that the tabs use.
3. Sixteen components are compared with reference values computed in R, scipy and statsmodels and embedded in the application.
4. The panel lists every compared value with its reference and tolerance, the reference package versions and the browser engine.

**Download evidence (JSON)** saves the full record. The analysis state is saved before the test and restored afterwards.

## Continuous integration

The repository's test suite runs the self-test, the demo-dataset analysis with pinned results, and the external validation in Chromium, Firefox and WebKit on every commit.

## Platform note

TMM's trimming depends on exact ranks of log-ratios. Differences in the last bits of floating-point sums between platforms can change a tie and shift a factor slightly. The implementation reproduces R's summation order.
