# Differential expression

## Methods

| Method | Model | Requires | Notes |
|---|---|---|---|
| Moderated t (eBayes) | Linear model on log₂(CPM + 1) per gene, empirical-Bayes variance moderation (limma `lmFit` + `eBayes`) | Any input | Default. Supports covariates. |
| voom | Precision weights from the mean–variance trend, then the moderated t (limma `voom`) | Raw counts | Recommended for counts, particularly when library sizes vary. |
| Welch t | Unequal-variance t-test per gene | Any input | No variance moderation; low power at small n. Provided for comparison. |

All three report log₂ fold change, t, p and Benjamini–Hochberg FDR. A gene is called significant when FDR is at or below the chosen level and |log₂FC| is at or above the chosen threshold.

## Which samples enter the model

A contrast between two groups is fitted on the samples of those two groups. The genome-wide moderated F test (Group Patterns tab) uses all groups. See [Limitations](Limitations.md) for the consequences.

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
