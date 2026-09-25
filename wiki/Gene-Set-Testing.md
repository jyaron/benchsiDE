# Gene-set testing

Three tests are offered. They answer different questions, so they can disagree without contradiction.

| Test | Hypothesis type | Question | Uses |
|---|---|---|---|
| ORA (over-representation) | Competitive | Are the significant genes over-represented in the set, relative to all kept genes? | The current DE hit list (up, down or both) |
| GSEA (preranked) | Competitive | Are the set's genes concentrated at one end of the genome-wide ranking? | All genes ranked by moderated t |
| FRY | Self-contained | Is the set, as a whole, differentially expressed between the groups? | Per-sample expression of the set's genes |

A competitive test compares the set with the other genes; a self-contained test compares the set with no change. A self-contained test can be significant for a set that is no more affected than the average gene, which is often the case when a large fraction of the transcriptome responds.

## ORA

One-sided hypergeometric test of the overlap between the hit list and each set. The universe is all genes kept after filtering. FDR by Benjamini–Hochberg across the tested sets.

## GSEA

Weighted Kolmogorov–Smirnov running enrichment score with weight p = 1 (hits weighted by |t|), computed on the moderated t ranking. Significance by gene-set permutation with a seeded generator (500 permutations for up to 200 sets, 250 for up to 1,000 sets, 100 above that). The normalized enrichment score (NES) divides by the mean of same-sign null scores. The smallest attainable p-value is 1/(permutations + 1), which limits FDR for large libraries. The leading edge is the set of member genes up to the running-score extremum.

## FRY

The fast rotation test of limma (`fry`, default settings), a closed-form approximation to ROAST. It reports a directional p-value (coordinated shift up or down) and a mixed p-value (any perturbation, including mixed directions). Results match limma to machine precision.

## Set size limits

Only sets with between 5 and 2,000 members present in the filtered data are tested by default. Change the limits with **set size**.

## Supporting displays

| Display | Content |
|---|---|
| Ridgeline | log₂FC distribution within each top set |
| Barcode plot | Positions of set members in the ranking, with the enrichment worm (limma `tricubeMovingAverage`) |
| Set-vs-rest ECDF | Cumulative distributions of the statistic for members and non-members, with the Kolmogorov–Smirnov distance |
| Jaccard heatmap | Pairwise overlap of top sets, to identify redundant results |
| Leading-edge matrix | Which genes drive which GSEA sets |

## Enrichment across contrasts

With three or more groups, **Enrichment across all contrasts** runs FRY independently for every pairwise contrast, with BH correction within each contrast. The heatmap shows signed −log₁₀ FDR for the 25 sets with the strongest evidence; the CSV export contains all sets.

## Cluster enrichment

On the Heatmap tab, when rows are ordered by k-means, each cluster is tested by ORA against the loaded library. Because the clusters were selected from the data, these results are descriptive.

## References

Subramanian A et al. (2005) PNAS 102:15545. Wu D et al. (2010) Bioinformatics 26:2176. Wu D and Smyth GK (2012) Nucleic Acids Res 40:e133.
