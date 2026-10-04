# Quality control

All QC displays are on the Overview & QC tab and use only the included samples.

| Display | What it shows | What to look for |
|---|---|---|
| Library size | Column totals | Samples far below the others |
| Detected genes | Genes with CPM ≥ 1 per sample | Low complexity or degraded samples |
| PCA | Top 2,000 variable genes, log₂ scale | Separation by group; outliers; batch structure |
| PC loadings | Genes driving each component | Whether a component reflects biology or artifact |
| MDS | limma plotMDS convention (leading log₂FC, top 500 genes per pair) | Same structure as PCA from a different distance |
| Sample dendrogram | Correlation distance, average linkage | Samples joining the wrong branch |
| Sample–sample correlation | Pearson correlation on log₂ values | A sample correlating better with another group |
| RLE | Relative log expression per sample | Boxes not centred on zero indicate a normalization problem |
| Expression density | Gaussian kernel density per sample (nrd0 bandwidth) | Curves that do not overlay |
| Library complexity and saturation | Detected genes against depth | Samples that would gain genes with more sequencing |
| Sex check (mouse, human) | Expression of sex-specific genes | Mismatch with recorded sex; sample swaps |
| Biotype composition | Share of expressed genes by biotype | Unexpected rRNA or mitochondrial content |

## Automatic warnings

The top of the tab lists flags such as low-depth samples and samples whose correlation profile resembles another group.

## Excluding samples

Untick a sample in the sample bar to exclude it from every statistic. Exclusions are stored in session files and reported in the methods text. Record the reason for any exclusion.

## Sample dendrogram

Average-linkage clustering on 1 − Pearson correlation of the 2,000 most variable genes (identical to R `hclust(as.dist(1 - cor(x)), "average")`). **colour leaves by** marks each sample by group or by another design column, so you can see which samples the tree places together. **cut into** divides the tree into clusters (as R `cutree`) and compares them with the chosen column: a cross-table, the adjusted Rand index (1 = the clusters reproduce the column exactly, 0 = no better than chance), and a list of samples that sit in a cluster dominated by another group. A low index, or samples that cluster with the wrong group, can mean outliers, swapped labels, or a batch, sex or quality effect larger than the condition effect: colour by those columns to check.
