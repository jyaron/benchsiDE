# Differential expression

## Methods

| Method | Model | Requires | Notes |
|---|---|---|---|
| Moderated t (limma-trend) | Linear model on log₂ CPM per gene, empirical-Bayes variance moderation toward a mean–variance trend (limma `lmFit` + `eBayes(trend = TRUE)`) | Any input | Default. Supports covariates. |
| voom | Precision weights from the mean–variance trend, then the moderated t (limma `voom`) | Raw counts | Recommended for counts, particularly when library sizes vary. |
| Welch t | Unequal-variance t-test per gene | Any input | No variance moderation; low power at small n. Provided for comparison. |

All three report log₂ fold change, t, p and Benjamini–Hochberg FDR. A gene is called significant when FDR is at or below the chosen level and |log₂FC| is at or above the chosen threshold.

## Which samples enter the model

By default (*variance from: all groups*) every included sample of every group is fitted, with design `~ 0 + group` plus any covariates, and the two chosen groups are compared by the contrast B − A; the residual variance and the empirical-Bayes prior therefore use all samples, as in limma's standard analysis. The option *these two groups only* fits the two groups alone. Welch's test always uses the two groups only. Samples with a missing covariate value are left out of a covariate-adjusted fit. Covariates named like identifiers (patient, subject, donor, pair, individual, id) are categorical by default; each covariate can be set to continuous or categorical.

## voom details

The implementation follows limma 3.66, including two details that affect results:

1. Genes with zero counts in all samples of the contrast are excluded from the trend fit.
2. The lowess span is chosen adaptively (`chooseLowessSpan`, the limma 3.66 default), not fixed at 0.5.

The fitted trend is shown in the **voom mean–variance trend** panel.

## Diagnostics

- **p-value histogram**: a uniform distribution with a spike near zero is expected. A U-shape or a spike near one suggests a model problem.
- **Q–Q plot**: observed against expected −log₁₀ p.
- **Stability curves**: the number of hits across FDR and fold-change thresholds, with the current thresholds marked. A result that depends on one exact threshold is fragile.

## Volcano and MA plots

Top genes can be labelled, ranked by p-value (conventional), by |log₂FC|, or by π score (|log₂FC| × −log₁₀ p). The ranking used is stated in figure captions.

## DESeq2

benchsiDE does not reimplement DESeq2. **DESeq2 script (R)** exports an R script that runs DESeq2 on the original files for the current contrast.

## References

Robinson MD and Oshlack A (2010) Genome Biol 11:R25. Smyth GK (2004) Stat Appl Genet Mol Biol 3:Article3. Law CW et al. (2014) Genome Biol 15:R29. Ritchie ME et al. (2015) Nucleic Acids Res 43:e47. Benjamini Y and Hochberg Y (1995) J R Stat Soc B 57:289.

## Interaction contrasts

The **Interaction** card under the volcano plot tests whether a difference between two groups changes between two conditions: the difference of differences (B2 − A2) − (B1 − A1). Typical questions are whether a treatment response differs between genotypes (Aldara vs control in knockout mice compared with wild type) or between time points.

1. Combine the two factors into one grouping under **Review → group by** (for example *Treatment · Genotype*), so that each group is one combination of levels.
2. Choose the reference condition (A1 → B1, for example Control_WT → Aldara_WT) and the compared condition (A2 → B2, for example Control_KO → Aldara_KO).
3. Choose the test (moderated t or voom), whether the variance is estimated from all groups or from the four groups only, and any covariates (for example Day), then **Run**.

The model is ~ 0 + group (+ covariates); the interaction and the two simple effects are linear contrasts of the group coefficients. A positive interaction log₂FC means the B2 − A2 difference is larger than the B1 − A1 difference. The table gives both simple effects, so a gene that responds in one condition only can be recognized. Significant genes can be passed to the Enrichment tab, the CSV export contains every gene, and the methods paragraph describes the model.

An interaction is estimated with the variance of four group means rather than two, so it needs more samples than a simple comparison; the card warns when the smallest of the four groups has fewer than three samples.

**Agreement with limma.** On GSE143688 (74 samples, eight groups), the moderated-t interaction adjusted for day agrees with `lmFit`, `contrasts.fit` and `eBayes(trend = TRUE)` to 5 × 10⁻¹² in log₂FC and 4.5 × 10⁻¹¹ in t. For voom, benchsiDE computes each gene's contrast standard error exactly from that gene's weighted fit and agrees to 5 × 10⁻¹² with limma when the interaction is fitted as a model coefficient (`~ treatment * genotype + day`). limma's `contrasts.fit` does not refit each gene and is therefore approximate when precision weights are combined with a non-orthogonal design (as its documentation states); in this example its t-statistics differed from the exact values by up to 0.31.
