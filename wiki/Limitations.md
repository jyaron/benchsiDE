# Limitations

## Statistical scope

- **Contrasts are fitted on the samples of the two groups compared.** limma's usual practice fits one model to all samples and estimates the variance from all of them. The two approaches give the same fold changes but different variance estimates and degrees of freedom.
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
