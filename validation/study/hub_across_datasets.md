# Hub across datasets: parity and calibration

Discovery dataset GSE54456; comparison datasets GSE121212, GSE186063 and GSE83645 (public GEO count matrices, NCBI-generated). Hub IL36A, 15 partners.

## Parity

| Component | Reference | Values | Maximum deviation |
|---|---|---|---|
| Seven preservation statistics, observed and 50 random sets, then Z-scores | WGCNA 1.74 `.coreCalcForExpr`, `moduleEigengenes`, `signedKME`, sourced from the CRAN source and run in R (stats::cor) | 7 statistics × 51 sets | statistics 1.7e-14; Z 4.2e-14 |
| RRA score, 2,000 genes × 4 connectivity lists | RobustRankAggreg `aggregateRanks(method = "RRA", exact = FALSE)` | 2,000 | 2.6e-15 relative; identical order |
| Fisher z: r per group and p | scipy.stats.pearsonr and the normal tail | 15 partners | 5.4e-14 relative |
| Member forest: pooled log₂FC with and without the discovery dataset | DerSimonian–Laird in numpy on the same inputs | 16 members | 1.8e-15 |

WGCNA's compiled dependencies could not be built in the validation environment, so its R code was sourced directly; `validation/reference.R` downloads the WGCNA source at run time for the permanent check and does not redistribute it.

## Calibration

Discovery dataset GSE54456 (psoriasis vs normal skin, 171 samples); comparison datasets GSE121212 (PSO lesional vs non-lesional, 54), GSE186063 (lesion vs non-lesion, 54) and GSE83645 (psoriasis vs uninvolved, 25). Neighbourhoods of 15 partners among the 2,000 most variable genes. Parity of the statistics with WGCNA 1.74, RobustRankAggreg and scipy is in the external validation table above and in `validation/study/hub_across_datasets.md`.

| Check | Expected under the null | Result | Consequence |
|---|---|---|---|
| Preservation of real neighbourhoods with gene correspondence shuffled (25 hubs) | Zsummary ≈ 0 | median -0.53 (all samples), -0.63 (within groups); none > 2 | the Z statistic is centred |
| Preservation of random 16-gene sets (40) | — | median Zsummary 3.3; 100% above the published Z = 2 threshold (57% within groups) | published thresholds not used; empirical baselines instead |
| Preservation of random sets of up-regulated genes against random sets (25 × 3 datasets) | — | 87% "preserved" (all samples; 76%–92% by dataset), 51% within groups (8%–84%) | DE-matched baseline added; the verdict requires preservation beyond it |
| Random up-regulated sets against the DE-matched baseline (20 × 3 datasets per mode) | 5% at p ≤ 0.05 | 3.3% (all samples), 3.3% (within groups) | calibrated |
| Random gene sets against the random-set baseline (30) | 5% at p ≤ 0.05 | 3% | calibrated |
| Fisher z on random splits of healthy skin, random gene pairs | 5% at p < 0.05 | 6 v 6: 4.9%; 10 v 10: 5.0%; 40 v 40: 4.9% | the test is calibrated on non-DE genes |
| Fisher z on random pairs of up-regulated genes | 5% if nothing rewires | GSE186063 40% (non-DE pairs 3%); GSE121212 9% (non-DE 6%) | raw Fisher z not used for calls; DE-matched baseline added |
| DE-matched rewiring p on random up-regulated pairs | 5% at p < 0.05 | A (this session): 4.7% (raw 39%), GSE121212: 5.0% (raw 7%), GSE186063: 4.1% (raw 38%) | calibrated |
| Consensus hubs with every dataset's co-expression destroyed (per-gene sample permutation, 2 runs) | 0 at FDR ≤ 0.05 | 0, 0 (real data: 145) | calibrated |
| Consensus hubs with only the discovery dataset intact (2 runs) | 0 | 0, 0 | one dataset cannot produce a consensus call |

Ten strongest hubs, three comparison datasets each (30 tests): preserved beyond DE-matched sets in 15 of 30 with correlation across all samples and in 30 of 30 within groups.

