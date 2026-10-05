# Randomized validation against R (benchsiDE 1.0.0-rc.2)

Fixed reference comparisons test the configurations someone chose. This folder tests configurations drawn at random: analysis method, fit, covariates (including covariates that cannot be estimated, have missing values, or identify subjects), TREAT and its threshold, FDR, |log2FC| threshold, excluded samples, input type and filtering, on public datasets. Every configuration is run in the application (JavaScriptCore, `tests/engine/harness.py`) and in an independent R implementation written from the documented rules, not from the application's code (`reference.R`, limma 3.66.0, edgeR 4.8.2). The comparison also checks that the texts describe the model that was fitted.

## How to run

From the repository root, after `python3 tests/engine/fetch_data.py`:

```
sh validation/fuzz/run_all.sh 20261004      # one seed; RSCRIPT=/path/to/Rscript if Rscript is not on PATH
python3 validation/fuzz/stat_tests.py       # Gene Explorer tests (scipy) and meta-analysis (metafor)
```

Results go to `validation/fuzz/work/` (not kept in the repository).

| File | Role |
|---|---|
| `gen_configs.py` | draws the configurations from a seed; writes the matrices, designs and `configs.json` |
| `run_app.py`, `run_ds.py`, `run_sets.py` | run them in the application: DE tab; comparison-dataset DE; FRY and ORA |
| `reference.R` | the R implementation (DE; with `sets`, FRY and ORA) |
| `compare.py`, `compare_ds.py`, `compare_sets.py` | statistics, calls and texts against R |
| `stat_tests.py`, `meta_reference.R` | Gene Explorer tests against scipy/statsmodels; meta-analysis against metafor 5.0.1 |

## Results (three seeds: 20261004, 77031, 31337)

Thirteen inputs: GSE41745, GSE54456 (40-sample subset, filterByExpr and CPM filtering), GSE63310 (TMM, and library-size CPM with CPM filtering), GSE83645, GSE121212, GSE143688 (counts and log2-CPM values), GSE171012 (bulk samples), GSE186063 (counts, CPM values and log2(CPM + 1) values); seven synthetic covariates added to each design (continuous, three-level, aliased with the groups, with missing values, a subject identifier, single-level, partly aliased).

| Analysis | Compared | Result |
|---|---|---|
| Differential expression (moderated t, voom, Welch t) | 1,146 configurations; 950 fitted by both, 196 refused by both | calls identical in all 950; max abs. difference in t 2.2 × 10⁻¹¹, max rel. difference in p 9.6 × 10⁻¹⁰; no text inconsistency |
| Comparison datasets (moderated t with covariates) | 125 configurations | calls identical in all 125; max abs. difference in t 2.2 × 10⁻¹² |
| Over-representation (hypergeometric) | 31,650 set tests | max rel. difference in p 8.0 × 10⁻¹¹ |
| FRY | 791 configurations, Hallmark sets | 724 equal to limma::fry (rel. difference ≤ 10⁻⁶); see below for the other 67 |
| Gene Explorer: ANOVA, Kruskal–Wallis, Spearman trend, pairwise Welch with Holm | 400 random designs, including ties and constant groups | no disagreement; max rel. difference 3.5 × 10⁻¹² |
| Meta-analysis (fixed effect; DerSimonian–Laird or REML; z, Hartung–Knapp or truncated Hartung–Knapp) | 3,000 random cases × 6 variants, k = 2–12, including τ² near 0, large heterogeneity and one very imprecise dataset | no disagreement in estimate, SE, p, τ², I², confidence and prediction intervals or Q (2,668 cases per variant; 332 cases with identical estimates in every dataset excluded, where Q and the Hartung–Knapp SE are rounding-level values instead of exactly 0) |

**FRY.** FRY's standardization uses the largest squared residual effect of each gene, which depends on the basis chosen for the residual space and therefore on the order of the samples; limma::fry itself gives different p-values when the samples of the same design are reordered. In each of the 67 configurations that differ, limma's own result changed when the samples were reordered (by up to 85%). In the seed examined in detail (31337), all 23 were fitted on all groups of datasets with 5 to 8 groups; none of the 133 two-group fits differed. The application's differences from limma (up to 77%) are of the same size, and at FDR ≤ 0.05 the application called 6 of 1,191 sets differently from limma, against 12 for limma with the samples reordered.

## Defects found and fixed in this round

1. Library-size CPM normalization used library sizes over all genes instead of the kept genes (edgeR keep.lib.sizes = FALSE).
2. Covariates that could not be estimated were listed as adjusted for in the CSV headers, the Interaction card, the contrast comparison, the comparison-dataset methods text and the meta-analysis sample table.
3. The methods text did not state that excluded samples took part in filtering and normalization.
4. The log-scale input note gave the wrong back-transformation for log2(x + 1) input.
5. Pairwise Welch tests between two constant groups were reported as p = 1 and counted in the Holm correction.
6. FRY refused designs with an aliased covariate instead of dropping it, and depended on an internal sample order instead of the file order.

## Coverage map

| Result-producing path | Fixed comparison with R or Python (earlier records) | Randomized (this folder) |
|---|---|---|
| Input types, filtering (filterByExpr, CPM), TMM, library-size CPM, log-CPM | yes, 7 datasets | yes |
| DE: moderated t, voom, Welch t; covariates; all-group and two-group fits; TREAT; excluded samples | yes | yes, with texts |
| Moderated F | yes | no |
| Interaction contrasts | yes (12 fits) | no |
| Comparison-dataset DE | yes | yes, with texts |
| Meta-analysis (6 variants), prediction interval, leave-one-out | yes (metafor), null calibration | yes (per-gene engine; leave-one-out not randomized) |
| Gene Explorer tests | yes | yes |
| ORA | yes (scipy) | yes |
| FRY | yes | yes |
| CAMERA, GSEA (sample and gene permutation) | yes (limma, fgsea), null calibration | no |
| Co-expression test, hub network, WGCNA preservation, robust rank aggregation, rewiring | yes | no |
| Sample traits (slope, duplicateCorrelation, per-subject means, CAMERA, eigengenes) | yes | no |
| Power calculation | yes (RNASeqPower) | no |

Not compared with a reference implementation (no standard implementation exists, or the output is descriptive): homology-family scores, the Discovery screen (null-calibrated only), clustering and expression patterns, QC flags, dendrogram and heatmaps, the HTML report text. Texts were checked against the fitted model for the DE tab (summary, methods, CSV header, captions, counts), comparison datasets (methods text, sample table) and FRY; the texts of the other analyses were not checked automatically.
