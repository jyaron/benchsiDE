# Cross-dataset comparison

## Principle

Each dataset is filtered, normalized and tested within itself. Raw or normalized values are never merged, pooled or renormalized across datasets, and no p-value is computed on expression pooled across datasets. Comparisons operate on within-dataset statistics of matched genes. Merging studies would require assumptions about batch structure that the tool does not make.

## Adding datasets

Up to eight datasets can be added to the session. Each has its own optional design file and species setting.

## Gene matching

| Situation | Matching |
|---|---|
| Same species | Case-normalized gene symbol |
| Mouse and human | MGI strict one-to-one orthologs |
| Custom ID map loaded | The map, applied before the rules above |

Identifiers are decoded with each dataset's own species annotation. If matching yields few genes, the slot reports the identifier types detected on each side and the remedy.

## Pairwise comparison

| Output | Statistic |
|---|---|
| Fold-change scatter | log₂FC in X against log₂FC in Y, coloured by quadrant; Pearson r |
| t concordance | Pearson correlation of moderated t-statistics |
| Hit-list overlap | One-sided hypergeometric (Fisher) test of shared up and shared down hits |
| RRHO map | Rank–rank hypergeometric overlap, computed as in the Bioconductor RRHO package; under-enrichment shown as a signed extension |
| Bland–Altman plot | Mean against difference of log₂FC |
| Shared and discordant lists | Genes significant in both datasets, ranked |

## All datasets

**Compute across all datasets** produces the correlation matrix and a consensus table. The **forest plot** shows one gene's log₂FC with confidence intervals in every dataset.

## Signature transfer

A module from dataset A (up- or down-regulated DE genes, or a Discovery module) is scored in each other dataset as a mean z-score. Discrimination is summarized as the Mann–Whitney AUC between that dataset's groups and calibrated against 200 seeded size-matched random signatures scored identically.

## References

Plaisier SB et al. (2010) Nucleic Acids Res 38:e169. Wu D and Smyth GK (2012) Nucleic Acids Res 40:e133.
