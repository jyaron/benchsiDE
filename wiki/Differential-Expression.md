# Differential expression

## Methods

| Method | Model | Requires | Notes |
|---|---|---|---|
| Moderated t (limma-trend) | Linear model on log₂ CPM per gene, empirical-Bayes variance moderation toward a mean–variance trend (limma `lmFit` + `eBayes(trend = TRUE)`) | Any input | Default. Supports covariates. |
| voom | Precision weights from the mean–variance trend, then the moderated t (limma `voom`) | Raw counts | Recommended for counts, particularly when library sizes vary. Supports covariates. |
| Welch t | Unequal-variance t-test per gene | Any input | No variance moderation; low power at small n; no covariates. Provided for comparison. |

All three report log₂ fold change, t, p and Benjamini–Hochberg FDR. A gene is called significant when FDR is at or below the chosen level and |log₂FC| is at or above the chosen threshold.

## Testing a fold-change threshold (TREAT)

With **test the threshold (TREAT)** ticked, the moderated t and voom test whether |log₂FC| is greater than the value in the |log₂FC| box, rather than whether it differs from zero (McCarthy & Smyth, Bioinformatics 25:765–771, 2009; limma `treat()`). The p-value is P(T > (|b| − τ)/SE) + P(T > (|b| + τ)/SE) on the posterior degrees of freedom, and genes are called when the TREAT FDR is at or below the chosen level. The t, p and FDR shown in the table, the volcano plot, the CSV export (columns `t_treat`, `p_treat`, `FDR_treat`), the Venn diagram and the all-contrasts table are then TREAT values, and every caption says so. The limma authors recommend TREAT over the default rule, which keeps genes whose *estimated* |log₂FC| reaches the threshold and so admits genes whose true change may be smaller. TREAT is stricter: on the demo dataset (LP vs Basal, voom adjusted for lane, FDR 0.05, threshold 1) it calls 3,647 genes against 6,133 with the default rule. Gene-set tests that rank genes (GSEA, the barcode plot) and the cross-dataset analyses always use the ordinary moderated t; over-representation analysis of significant genes uses the genes called under the current rule.

## Covariates with voom

Covariates can be used with the moderated t and with voom. With voom, the design is `~ 0 + group + covariates`; voom estimates the precision weights from the fit of that design, and the B − A contrast is taken with limma's `contrasts.fit`. When precision weights are combined with a non-orthogonal design (any covariate not balanced across groups), `contrasts.fit` computes the standard error approximately, from each gene's weighted coefficient standard errors and the coefficient correlation of the unweighted design; its documentation states this. benchsiDE uses the same rule so that its results equal the standard limma workflow. The Welch t does not support covariates; with covariates selected it shows an unadjusted warning.

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

**Agreement with limma.** On GSE143688 (74 samples, eight groups), the moderated-t interaction adjusted for day agrees with `lmFit`, `contrasts.fit` and `eBayes(trend = TRUE)` to 5 × 10⁻¹² in log₂FC and 4.5 × 10⁻¹¹ in t. The voom interaction adjusted for day agrees with `voom`, `lmFit`, `contrasts.fit` and `eBayes` to 5 × 10⁻¹² in log₂FC and t. With precision weights and a covariate, `contrasts.fit` uses limma's documented approximation for the standard error (see *Covariates with voom*), and benchsiDE follows it.
