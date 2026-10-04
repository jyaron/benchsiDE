# Differential expression in 1.0.0-rc.2: covariates with voom, and TREAT

## Changes

- **voom with covariates.** The Differential Expression tab fits voom with the design `~ 0 + group + covariates`. As in limma, voom estimates the precision weights from the fit of that design, each gene is fitted by weighted least squares, and the B − A contrast is taken with `contrasts.fit`. When precision weights are combined with a non-orthogonal design, `contrasts.fit` computes the contrast standard error from each gene's coefficient standard errors and the coefficient correlation of the unweighted design, which is an approximation (limma documentation). benchsiDE applies the same rule, here and in the Interaction card, so that its results equal the standard limma workflow.
- **TREAT** (McCarthy & Smyth, Bioinformatics 25:765–771, 2009), as limma `treat()`: p = P(T > (|b| − τ)/SE) + P(T > (|b| + τ)/SE) with T on the posterior degrees of freedom; t = sign(b)·max(0, (|b| − τ)/SE) when |b| > τ, otherwise 0. Available for the moderated t and voom.
- **Covariate terms that cannot be estimated** (aliased with the groups or with other covariates) are dropped, as `lmFit` does. The summary line and methods text now state this, and say that the comparison is not adjusted for a covariate when all of its terms are dropped.

## Parity with limma 3.66.0 / edgeR 4.8.2

Reference: `de_reference.R` (filterByExpr, TMM, `~ 0 + group + covariates`, `lmFit` on log-CPM with `eBayes(trend = TRUE)`/`treat(trend = TRUE)`, or `voom` + `lmFit` + `contrasts.fit` + `eBayes`/`treat`). 80 configurations: two datasets × method (moderated t, voom) × covariates × fit (all groups, two groups) × test (ordinary, TREAT τ = 1).

- GSE63310 (9 samples, 16,624 genes): contrasts LP vs Basal and ML vs LP; covariates none, sequencing lane (the article's), lane + a continuous covariate. The continuous covariate (`rin`) holds synthetic values and exists only to exercise a numeric covariate.
- GSE186063 (65 samples, 3 groups): lesion vs non-lesion; covariates Pair (blocking, normal skin has no pair and is left out), Age (continuous), Sex + Age, Pair + Age (Age is aliased with Pair and dropped, as in limma).

Maximum deviations over all 80 configurations and all genes: log₂FC 1.4e-12 (absolute), t 1.1e-11 (absolute), p 1.4e-10 (relative). Prior degrees of freedom agreed to 7 × 10⁻¹⁴ (relative), and residual degrees of freedom and sample counts exactly. The number of significant genes (FDR ≤ 0.05 with |log₂FC| ≥ 1, or TREAT FDR ≤ 0.05) was identical in all 80 configurations. Per-configuration values: `de_parity_80.csv`.

Under TREAT, no gene with |log₂FC| ≤ τ reached FDR ≤ 0.01, 0.05 or 0.10 in any configuration, as the definition requires.

## Interaction contrasts with voom and a covariate

Of the 12 interaction fits checked on GSE143688 for 1.0.0-rc.1 (`validation/freeze/PARITY_FREEZE.md`), two use voom with a covariate (Day) and change with the contrasts.fit rule. Both were recomputed against `voom(~ 0 + cell + day)`, `lmFit`, `contrasts.fit`, `eBayes`:

| Fit | max abs. diff. log₂FC | max abs. diff. t | max rel. diff. p | up/down app (R) |
|---|---|---|---|---|
| (Aldara_KO − Control_KO) − (Aldara_WT − Control_WT), voom, all groups + Day | 8.0e-14 | 3.1e-13 | 8.4e-12 | 314/455 (314/455) |
| (Aldara_KO − Control_KO) − (Aldara_WT − Control_WT), voom, four groups + Day | 5.0e-12 | 5.1e-12 | 1.1e-11 | 314/451 (314/451) |

The other ten fits are unchanged (no covariate, or the moderated t, for which contrasts.fit is exact).

## The workflow article (Law et al., F1000Research 5:1408, 2016)

The article analyses GSE63310 with `~ 0 + group + lane`, voom, `contrasts.fit` and `treat(lfc = 1)`. With limma 3.66.0, the LP vs Basal, ML vs Basal and shared significant genes number 3,647, 3,831 and 2,782; benchsiDE gives the same numbers. The article, computed with Bioconductor 3.8, reports 3,648, 3,834 and 2,784. This analysis is part of the external validation (`validation/compare.js`, 11 new checks, 70 in total).

## Null calibration

Healthy skin from GSE54456 (81 samples): 100 random splits into two groups of n = 3, 5 and 10, each with a random two-level block covariate (balanced overall, assigned independently of the groups; in 13 of the n = 3 splits it coincided with the groups and was dropped). Entries: splits with at least one gene at FDR ≤ 0.05 (TREAT: TREAT FDR ≤ 0.05). Under the global null at most about 5 splits in 100 are expected; the failure bound is the 99th percentile of Binomial(100, 0.05) = 11. Every count was reproduced exactly by limma in R on the same splits (2,100 analyses).

| Test | n = 3 | n = 5 | n = 10 |
|---|---|---|---|
| moderated t | 23 / 100 | 7 / 100 | 2 / 100 |
| voom | 17 / 100 | 5 / 100 | 1 / 100 |
| moderated t + block | 17 / 100 | 5 / 100 | 3 / 100 |
| voom + block | 15 / 100 | 5 / 100 | 4 / 100 |
| TREAT, moderated t, τ = 1 | 10 / 100 | 1 / 100 | 0 / 100 |
| TREAT, voom + block, τ = 1 | 9 / 100 | 1 / 100 | 0 / 100 |
| TREAT, voom + block, τ = 0.58 | 10 / 100 | 2 / 100 | 0 / 100 |

Mean fraction of genes with p < 0.05 (nominal 0.05 for the ordinary tests; TREAT p-values are conservative when the true |log₂FC| is below τ):

| Test | n = 3 | n = 5 | n = 10 |
|---|---|---|---|
| moderated t | 0.057 | 0.045 | 0.042 |
| voom | 0.056 | 0.045 | 0.042 |
| moderated t + block | 0.059 | 0.046 | 0.040 |
| voom + block | 0.057 | 0.046 | 0.041 |
| TREAT, moderated t, τ = 1 | 0.003 | 0.001 | 0.000 |
| TREAT, voom + block, τ = 1 | 0.004 | 0.001 | 0.000 |
| TREAT, voom + block, τ = 0.58 | 0.010 | 0.004 | 0.001 |

With five or more samples per group every test stayed within the bound. With three per group, the ordinary moderated t and voom, with or without a covariate, exceeded it (15–23 of 100), as documented for limma at this size in skin (VALIDATION.md). TREAT stayed within the bound at every size.

## Files

- `de_reference.R`: the R reference function.
- `parity_configs.json`, `de_parity_80.csv`: configurations and results.
- `null_splits.json`, `null_calibration.csv`: the null splits and results.
- `make_engine_reference.R`, `engine_configs.json`: the stored reference for the engine check `covariates_treat` (`tests/engine/de_rc2_reference.json.gz`).
