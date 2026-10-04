# Interaction contrasts (design-freeze change)

Dataset: GSE143688, author-supplied counts, all 74 samples: genotype (wild type, Il36r fl/fl, Il36r−/−, keratinocyte deletion) × treatment (Aldara, control) × day (3, 7). Grouping factor: Cell = treatment_genotype (eight groups). Contrast: (Aldara_KO − Control_KO) − (Aldara_WT − Control_WT). Reference: `interaction_reference.R` (limma 3.66.0, edgeR 4.8.2); the engine check `interaction` in `tests/engine/run_checks.py` repeats the comparison.

| Fit | Reference in R | max abs. diff. log₂FC | max abs. diff. t | max rel. diff. p | prior df diff. | genes at FDR ≤ 0.05, \|log₂FC\| ≥ 1 (app / R) |
|---|---|---|---|---|---|---|
| moderated t, all eight groups, + Day | lmFit(~0 + cell + day), contrasts.fit, eBayes(trend = TRUE) | 5.0e-12 | 4.5e-11 | 1.3e-11 | 4.1e-14 | 871 / 871 |
| voom, four groups, + Day | voom(~0 + cell + day), lmFit, contrasts.fit, eBayes | 5.0e-12 | 5.1e-12 | 1.1e-11 | 2.9e-14 | 765 / 765 |
| moderated t, four groups | lmFit(~0 + cell), contrasts.fit, eBayes(trend = TRUE) | 8.9e-15 | 7.7e-14 | 6.8e-12 | 3.9e-14 | 701 / 701 |
| voom, all eight groups | voom(~0 + cell), lmFit, contrasts.fit, eBayes | 9.2e-14 | 1.9e-13 | 8.6e-12 | 3.2e-14 | 644 / 644 |

The filtered gene set (15,492 genes) and TMM factors were identical to edgeR (max. difference 2.2e-16). With voom and a covariate, limma's `contrasts.fit` computes the contrast standard error approximately when precision weights are combined with a non-orthogonal design (its documentation states this); from 1.0.0-rc.2 benchsiDE follows the same rule, so that results equal the standard limma workflow (the exact per-gene value, obtained by fitting the interaction as a coefficient, gives 760 significant genes).

Example result: Il22 is induced by Aldara in wild-type skin (log₂FC 5.89) but not in Il36r−/− skin (0.27); interaction log₂FC −5.62, FDR 9.1e-13.
