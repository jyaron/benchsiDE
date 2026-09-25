# Heatmaps

## Heatmap tab

| Setting | Options |
|---|---|
| Rows | Top 25, 50, 100 or 200 most variable genes, or top DE genes of the current contrast |
| Row order | Hierarchical clustering, or k-means with k chosen by silhouette |
| Columns | Sample labels or group labels |

Values are z-scores per gene across the included samples. The colour scale runs from blue (low) to red (high). The figure height grows with the number of rows so that every gene label remains legible.

With k-means ordering, **Cluster enrichment** tests each cluster against the loaded gene-set library (see [Gene-set testing](Gene-Set-Testing.md)), and **Export clusters** saves the assignment.

## Module heatmaps

From Discovery cards and enrichment gene panels, **Heatmap (all N)** exports every member gene of a module as a z-scored heatmap, with a module-score strip above the matrix and group dividers through both. There is no cap on the number of genes.
