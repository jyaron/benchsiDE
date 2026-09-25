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
