# Patterns and clustering

## Moderated F test

**Genome-wide moderated F (any group differs)** tests, for each gene, whether any group mean differs from the others, using all samples and empirical-Bayes variance moderation. It is the appropriate first test in designs with more than two groups.

## Pattern clustering

1. The most variable genes (500, 1,000 or 2,000) are selected.
2. Each gene's group-mean profile is z-scored.
3. Profiles are clustered by k-means with seeded k-means++ initialization, so repeated runs give identical results.

## Choosing k

With **auto (silhouette)**, k from 2 to 10 is evaluated and the k with the highest mean silhouette width is chosen. The per-k silhouette scores are reported.

| Best mean silhouette | Interpretation |
|---|---|
| Above 0.5 | Well-separated clusters |
| 0.3 to 0.5 | Moderate structure |
| Below 0.3 | Profiles form a continuum; treat clusters as a descriptive partition |

The silhouette criterion replaced an earlier elbow criterion, which selected four or five clusters regardless of the true structure in benchmarks with known cluster numbers.

## Trajectory heatmap

For numeric designs, the Trajectory heatmap shows every selected gene's profile over the ordered groups.

## Gene-set module score

Paste a list of genes into **Gene-set module score** to plot their mean z-score per sample.

## Reference

Rousseeuw PJ (1987) J Comput Appl Math 20:53.
