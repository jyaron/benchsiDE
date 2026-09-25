# Glossary

| Term | Definition |
|---|---|
| BH / FDR | Benjamini–Hochberg procedure; controls the expected proportion of false discoveries among significant results |
| CAMERA | Competitive gene-set test that adjusts for inter-gene correlation (limma) |
| Competitive test | Tests whether a set is more affected than other genes |
| CPM | Counts per million: counts divided by library size, times 10⁶ |
| eBayes | Empirical-Bayes moderation of gene-wise variances toward a common prior (limma) |
| filterByExpr | edgeR's design-aware low-expression filter |
| FRY | Fast rotation gene-set test; self-contained (limma) |
| GSEA | Gene set enrichment analysis on a ranked gene list |
| Hidden module | A module significant as a whole while fewer than 25% of its members are individually significant |
| Homology family | A group of mouse and human genes connected through shared homology classes |
| Leading edge | The GSEA set members that drive the enrichment score |
| log₂FC | Log₂ fold change between two groups |
| MDS | Multidimensional scaling; limma's plotMDS uses leading log₂FC distances |
| Moderated F | Empirical-Bayes F test for any difference among groups |
| Moderated t | Empirical-Bayes t test (limma) |
| Module score | Per-sample mean z-score of a gene set's members |
| NES | Normalized enrichment score (GSEA) |
| ORA | Over-representation analysis; hypergeometric test of list overlap |
| RLE | Relative log expression: each gene's value minus its median across samples |
| RRHO | Rank–rank hypergeometric overlap |
| Self-contained test | Tests whether a set is differentially expressed at all |
| Silhouette | Measure of how well each item fits its cluster relative to the nearest other cluster (−1 to 1) |
| TMM | Trimmed mean of M-values normalization (edgeR) |
| voom | Precision weights from the mean–variance relationship of log-counts (limma) |
