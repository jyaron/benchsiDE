# Limitations

## Statistical scope

- **Small groups.** With 3 or fewer samples in a group, FDR control is approximate: in benchsiDE's calibration on random 3-versus-3 splits of healthy skin, 10–25% of splits gave at least one gene at FDR ≤ 0.05. The app shows a warning at this size.
- **Discovery screen.** CAMERA with estimated inter-gene correlation is calibrated but has low power at small group sizes, and the screen does not use DE covariates.
- **Covariate adjustment** is available for the moderated t method only.
- **No negative-binomial GLM methods.** DESeq2 and edgeR's quasi-likelihood methods are not implemented; a DESeq2 script can be exported.
- **No mixed models.** Repeated measures are handled by fixed-effect blocking only.
- **GSEA** uses gene-set permutation, so its p-values have a floor of 1/(permutations + 1).

## Discovery

- Gene families are defined by symbol stems, which group paralogs of unrelated function in some cases.
- HGNC groups are applied to non-human data by symbol, an approximation.
- With small groups, few or no modules survive FDR; the nominal tier is exploratory.

## Cross-dataset

- Only mouse–human orthology is built in. Other species pairs require a custom map.
- Comparisons are between within-dataset statistics; batch-corrected integration of raw data is outside the tool's scope.

## Scale

Normalization and differential expression complete in seconds for typical bulk RNA-seq matrices. Above about three million matrix cells (genes × samples), the Discovery screen, GSEA and k-means become noticeably slower, and the application displays a warning. Very large matrices can exceed browser memory.

## Annotation releases

Embedded annotation and gene-set libraries are fixed at their download dates (recorded in the labels). Newer releases can be supplied as files.
- **Meta-analysis on small studies.** With the degrees-of-freedom adjustment, false shared genes on null data are reduced but not eliminated: in benchsiDE's calibration (three datasets of 3 versus 3 healthy-skin samples, six replicates), genes at meta FDR ≤ 0.05 remained in 3 of 6 replicates under random effects (2–28 genes) and in 6 of 6 under fixed effect (1–28 genes). Use the random-effects model and treat shared signatures from small studies as candidates.
- **Hub across datasets.** Preservation and rewiring are calibrated empirically within each dataset (100 random sets, 1,000 random pairs), so the smallest attainable p is about 0.01 and 0.001; the DE-matched baselines match on log₂FC decile only, not on expression level or variance. Rewiring needs 6 or more samples per group and assumes approximately bivariate-normal log expression. Consensus hubs depend on the gene universe and on β; robust rank aggregation assumes the datasets are independent.
- **Sample traits.** Associations are correlational, and no dataset with a recorded per-sample clinical severity score was available for a real positive control (validation used recorded age, shuffles and planted effects). For a trait that varies within subjects, the mixed model uses one consensus intra-block correlation (limma duplicateCorrelation); with few subjects this can understate the correlation for some genes.
