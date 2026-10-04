# benchsiDE 0.22.0-beta: validation of gene-set analyses, Discovery and heatmap clustering

Build under test: `benchside/index.html`, APP_VERSION 0.22.0-beta (2026-09-30), from the frozen repository archive (sha256 in MANIFEST.json, 17511676...). The build was not modified. All app numbers were produced by executing the app's own script in macOS JavaScriptCore with the kit's `harness.py` DOM/Plotly stub. References: R 4.5.3, limma 3.66.0, edgeR 4.8.2, fgsea 1.36.2, cluster (conda r-cluster), and Python numpy/scipy/statsmodels. Gene-set libraries: the app's built-in Hallmark, GO Biological Process and Reactome (human and mouse), decoded with `decodeBuiltin` and written to GMT files that were used unchanged by every reference.

Every number below was produced in this task. The full list of 652 checks (636 pass, 16 fail) is in `VAL_ENR_checks.csv`; the rerunnable pipeline is `VAL_ENR_runner.py` with the R scripts listed at the end.

## 1. Scope and inputs

| Dataset | Contrast (test vs reference) | Samples in contrast | Filtered genes | Notes |
|---|---|---|---|---|
| GSE54456 | Psoriasis_skin vs normal_skin (Cond) | 90 vs 81 | 19518 | no patient IDs |
| GSE186063 | lesion vs non-lesion (Type) | 26 vs 28 | 21112 | Pair inferred; normal_skin (ankylosing spondylitis donors) excluded from the contrast |
| GSE121212 | PSO_lesional vs PSO_non_lesional (Group) | 27 vs 27 | 23168 | AD and control groups present in the session; Patient pairing available |
| GSE41745 | lesional vs non_lesional (Cond) | 2 vs 3 | 18036 | small n |
| GSE63310 (mouse) | LP vs Basal (celltype) | 3 vs 3 | 16624 | 3 cell types |

Session settings: built-in annotation, TMM, filterByExpr, moderated t (limma-trend), FDR 0.05, |log2FC| >= 1, set size 5-2000 (ORA/FRY/GSEA) and 3-100 (Discovery), 1000 Discovery permutations, GSEA permutation "auto". The app's log2 matrix was exported and agreed with `edgeR::cpm(log=TRUE, prior.count=2)` after `filterByExpr` and TMM to at most 7.1e-15 on every dataset; the moderated t agreed with `eBayes(trend=TRUE)` to at most 1.1e-12. Gene-set tests in R therefore ran on the same matrix and the same ranking.

An earlier internal audit had already reported the gene-permutation GSEA calibration problem, the CAMERA correlation choice and the family-stem fix; those points are not repeated except where this validation adds evidence. an earlier internal audit listed heatmap clustering, auto-k, the ORA universe, and the Discovery size and hidden-module rules as not assessed; they are assessed here.

## 2. Summary of all checks

| Check | Rows | Values compared | Max deviation | Tolerance | Failed rows |
|---|---|---|---|---|---|
| Calibration: Discovery CAMERA (estimated inter-gene correlation, app setting) on label-shuffled data | 3 | 8,799 | 0 | 0 | 0 |
| Calibration: Discovery CAMERA path on label-shuffled data (Hallmark + families + HGNC) | 3 | 8,799 | 0.8793 | 2 | 0 |
| Calibration: Discovery permutation path on label-shuffled data (Hallmark + families + HGNC) | 9 | 26,028 | 1.789 | 2 | 0 |
| Calibration: ORA on annotation-biased random queries, app universe | 2 | 400 | 17.39 | 0.05 | 1 |
| Discovery >2-group module F statistic | 3 | 8,676 | 2.09e-12 | 1.00e-08 | 0 |
| Discovery >2-group random-set permutation p within Monte-Carlo error | 3 | 8,676 | 3.258 | 4 | 0 |
| Discovery BH FDR across candidates | 15 | 81,928 | 2.46e-11 | 1.00e-06 | 0 |
| Discovery CAMERA ignores DE covariates (Patient) without notice | 1 | 3,009 | 0 | 0 | 0 |
| Discovery CAMERA inter-gene correlation | 15 | 81,928 | 1.63e-14 | 1.00e-08 | 0 |
| Discovery CAMERA p on label-shuffled data | 3 | 8,799 | 1.15e-11 | 1.00e-06 | 0 |
| Discovery CAMERA p, inter.gene.cor=NA | 15 | 81,928 | 2.46e-11 | 1.00e-06 | 0 |
| Discovery CAMERA with an included one-sample third group | 1 | 50 | 0.01002 | 1.00e-06 | 1 |
| Discovery candidate list (families + HGNC + library, size 3-100) | 18 | 90,604 | 0 | 0 | 0 |
| Discovery curated HGNC groups mapped to genes | 5 | 7,469 | 0 | 0 | 0 |
| Discovery gene-family stems (familyStems rule) | 5 | 25,569 | 0 | 0 | 0 |
| Discovery heatmap-export caption names the test used (CAMERA mode) | 1 | 1 | 0 | 0 | 0 |
| Discovery leave-one-out direction text, top 20 modules | 3 | 60 | 0 | 0 | 0 |
| Discovery module score effect (delta mean z) | 15 | 81,928 | 1.11e-15 | 1.00e-08 | 0 |
| Discovery nominal tier and hidden-module counts in summary text | 15 | 81,928 | 0 | 0 | 0 |
| Discovery panel-figure caption names the test used (CAMERA mode) | 1 | 1 | 1 | 0 | 1 |
| Discovery per-module count of individually significant genes | 15 | 81,928 | 0 | 0 | 0 |
| Enrichment matrix: FRY BH FDR per pairwise contrast | 6 | 159,804 | 4.82e-10 | 1.00e-06 | 0 |
| Enrichment matrix: direction calls | 6 | 159,804 | 0 | 0 | 0 |
| Enrichment matrix: plotted signed -log10 FDR and top-25 selection | 6 | 1,350 | 0 | 1.00e-12 | 0 |
| FRY BH FDR directional/mixed | 15 | 87,404 | 3.89e-10 | 1.00e-06 | 0 |
| FRY direction | 15 | 43,702 | 0 | 0 | 0 |
| FRY directional p | 15 | 43,702 | 3.89e-10 | 1.00e-06 | 0 |
| FRY mixed p | 15 | 43,702 | 4.03e-10 | 1.00e-06 | 0 |
| FRY with pairing covariate Pair: directional and mixed p | 2 | 3,470 | 14.52 | 1.00e-06 | 2 |
| FRY with pairing covariate Patient: directional and mixed p | 2 | 3,474 | 8.70e-11 | 1.00e-06 | 0 |
| GSEA enrichment score, genes permutation run | 21 | 60,624 | 9.88e-15 | 1.00e-08 | 0 |
| GSEA enrichment score, samples permutation run | 9 | 26,780 | 3.50e-13 | 1.00e-08 | 0 |
| GSEA leading-edge size | 15 | 43,702 | 0 | 0 | 0 |
| GSEA plot title in sample-permutation mode | 1 | 1 | 1 | 0 | 1 |
| GSEA preranked (gene permutation) p within Monte-Carlo error | 15 | 43,702 | 3.297 | 4 | 0 |
| GSEA sample permutation ES identical to independent R ES | 9 | 26,780 | 3.50e-13 | 1.00e-08 | 0 |
| GSEA sample permutation NES/p/FDR distributional agreement | 9 | 26,780 | 5 | 4 | 0 |
| GSEA sample permutation seeded replay: NES, p, FDR | 9 | 26,780 | 0.008264 | 1.00e-08 | 0 |
| Heatmap cluster enrichment: BH FDR within cluster | 5 | 3,550 | 7.52e-11 | 1.00e-06 | 0 |
| Heatmap cluster enrichment: hypergeometric p | 5 | 3,550 | 8.44e-11 | 1.00e-06 | 0 |
| Heatmap cluster enrichment: set sizes and overlaps | 5 | 3,550 | 0 | 0 | 0 |
| Heatmap gene selection (top DE by p, N=200) | 5 | 1,000 | 0 | 0 | 0 |
| Heatmap gene selection (top DE by p, N=50) | 5 | 250 | 0 | 0 | 0 |
| Heatmap gene selection (top variable, N=200) | 5 | 1,000 | 0 | 0 | 0 |
| Heatmap gene selection (top variable, N=50) | 5 | 250 | 0 | 0 | 0 |
| Heatmap hierarchical clustering: displayed order is a leaf order of hclust(average, 1 - Pearson r) | 5 | 2,480 | 0 | 0 | 0 |
| Heatmap k-means total within-cluster SS relative to R kmeans(nstart=100), k=2..10 | 5 | 180 | 0.3433 | 0.01 | 5 |
| Heatmap row z-scores | 20 | 197,000 | 1.91e-14 | 1.00e-12 | 0 |
| Heatmap/auto-k mean silhouette, partitions with singleton clusters | 5 | 56 | 0.08 | 1.00e-08 | 5 |
| Heatmap/auto-k mean silhouette, partitions without singleton clusters | 5 | 124 | 3.83e-09 | 1.00e-08 | 0 |
| ORA: BH FDR | 45 | 131,106 | 2.00e-10 | 1.00e-06 | 0 |
| ORA: hypergeometric p | 45 | 131,106 | 2.00e-10 | 1.00e-06 | 0 |
| ORA: overlap gene lists | 45 | 131,106 | 0 | 0 | 0 |
| ORA: query size | 45 | 45 | 0 | 0 | 0 |
| ORA: set sizes m and overlaps k | 45 | 131,106 | 0 | 0 | 0 |
| built-in library: set names, sizes and gene contents | 6 | 18,446 | 0 | 0 | 0 |
| ranking statistic: moderated t (limma-trend) used for ORA query and GSEA | 5 | 98,458 | 1.08e-12 | 1.00e-08 | 0 |
| set mapping to filtered genes, size 5-2000 | 15 | 43,702 | 0 | 0 | 0 |
| upstream log-CPM matrix vs edgeR cpm(log,prior.count=2) after filterByExpr+TMM | 5 | 8,285,846 | 7.11e-15 | 1.00e-08 | 0 |

## 4. Results by feature

### 4.1 Built-in libraries against MSigDB v2024.1

The six libraries were compared with the MSigDB release 2024.1 GMT files downloaded from data.broadinstitute.org (h.all, c5.go.bp and c2.cp.reactome v2024.1.Hs; mh.all, m5.go.bp and m2.cp.reactome v2024.1.Mm). Set names and gene contents are identical for all sets: 50, 7,608 and 1,736 human sets and 50, 7,713 and 1,289 mouse sets. Gene order within sets differs, which has no effect. The set counts in the app labels match.

### 4.2 Over-representation analysis

For all five datasets, three libraries and three queries (DE up, down, up+down), the query size equals the number of genes with limma-trend FDR <= 0.05 and |log2FC| >= 1; set sizes after mapping, overlaps and overlap gene lists are identical to an independent R mapping; hypergeometric p agrees with `phyper(k-1, m, N-m, q, lower.tail=FALSE)` to at most 2.0e-10 (relative) and BH FDR to 2.0e-10. The smallest ORA p was 2.2e-70, so no value was below 1e-300. The concern is the universe.

### 4.3 FRY

Directional p, mixed p, both BH FDR columns and the direction agree with `limma::fry(sort="none")` for every set of every library and dataset (maximum relative deviation 4.0e-10; smallest p 6.1e-224). Sets at FDR <= 0.05 (directional): GSE54456 42/50, 6,294/7,149, 1,464/1,667; GSE63310 30/50, 4,774/7,172, 681/1,216. These large fractions are expected for a self-contained test, as the app states.

### 4.4 CAMERA

CAMERA is used only in the Discovery tab. On all 81,928 candidate modules (15 dataset-library runs, two groups each), p with `inter.gene.cor = NA`, the inter-gene correlation, and BH FDR agree with `limma::camera` to 2.5e-11, 1.6e-14 and 2.5e-11. With limma's default fixed correlation of 0.01, the same candidates give the following numbers of modules at FDR <= 0.05 (last column), against none with the estimated correlation:

| Dataset (2 groups) | Library | Candidates | FDR<=0.05 | p<=0.01 (nominal tier) | expected by chance | hidden (nominal) | min p | FDR<=0.05 with inter.gene.cor=0.01 |
|---|---|---|---|---|---|---|---|---|
| GSE54456 | hallmark_hs | 2933 | 0 | 47 | 29 | 0 | 8.45e-04 | 114 |
| GSE54456 | react_hs | 4456 | 0 | 100 | 45 | 0 | 1.03e-04 | 412 |
| GSE54456 | gobp_hs | 9341 | 0 | 266 | 93 | 0 | 1.01e-04 | 731 |
| GSE41745 | hallmark_hs | 2603 | 0 | 4 | 26 | 0 | 5.72e-03 | 97 |
| GSE41745 | react_hs | 4120 | 0 | 4 | 41 | 0 | 6.96e-03 | 432 |
| GSE41745 | gobp_hs | 9023 | 0 | 15 | 90 | 8 | 3.26e-03 | 498 |
| GSE186063 | hallmark_hs | 2996 | 0 | 18 | 30 | 0 | 1.18e-03 | 137 |
| GSE186063 | react_hs | 4516 | 0 | 51 | 45 | 0 | 1.18e-03 | 523 |
| GSE186063 | gobp_hs | 9395 | 0 | 93 | 94 | 0 | 1.18e-03 | 846 |
| GSE121212 | hallmark_hs | 3009 | 0 | 7 | 30 | 0 | 1.58e-03 | 115 |
| GSE121212 | react_hs | 4515 | 0 | 12 | 45 | 0 | 1.58e-03 | 532 |
| GSE121212 | gobp_hs | 9410 | 0 | 18 | 94 | 0 | 1.54e-03 | 676 |
| GSE63310 | hallmark_mm | 2671 | 0 | 0 | 27 | 0 | 1.70e-02 | 3 |
| GSE63310 | react_mm | 3807 | 0 | 0 | 38 | 0 | 1.70e-02 | 26 |
| GSE63310 | gobp_mm | 9133 | 0 | 0 | 91 | 0 | 1.27e-02 | 15 |

On label-shuffled GSE54456 (no true effect), the estimated correlation gave 0 modules at FDR <= 0.05, and the fixed 0.01 gave 78-97 (section 5). This confirms the audit's choice; the price is that no module reached FDR <= 0.05 in any of the 15 dataset-library runs on real contrasts, including IFN-alpha response, which is ranked first by the fixed-correlation test in GSE54456, GSE186063 and GSE41745.

### 4.5 GSEA

Enrichment score. For every set of every library and dataset, in both sample-permutation and gene-permutation runs, ES agrees with `fgsea::calcGseaStat(gseaParam=1)` on the limma-trend ranking to at most 3.5e-13. Leading-edge sizes follow the fgsea peak rule for every set.

Sample-label permutation (GSE54456, GSE186063, GSE121212; GSE41745 and GSE63310 have fewer than 1000 labelings and fall back to gene permutation). Two comparisons were made.

(a) Exact replay of the app's own arithmetic. `VAL_ENR_gsea_replay.py` re-implements the permutation loop from the algorithm description (mulberry32 seed 0x6E5A11 + nT*7919 + nG, Fisher-Yates shuffle, moderated t with the trend prior held fixed, float32 null storage, NES, p = (1 + #)/(1 + #same sign), pooled-null FDR made monotone):

| Dataset | Library | Sets | Permutations | max abs dES (observed) | sets with identical NES and p (within 1e-9) | sets differing by one tied null draw |
|---|---|---|---|---|---|---|
| GSE54456 | hallmark_hs | 50 | 1000 | 1.1e-15 | 50 | 0 |
| GSE54456 | gobp_hs | 7149 | 250 | 6.0e-15 | 7146 | 3 |
| GSE54456 | react_hs | 1667 | 500 | 1.4e-15 | 1667 | 0 |
| GSE186063 | hallmark_hs | 50 | 1000 | 8.9e-16 | 50 | 0 |
| GSE186063 | gobp_hs | 7213 | 250 | 3.1e-15 | 7211 | 2 |
| GSE186063 | react_hs | 1685 | 500 | 1.2e-15 | 1684 | 1 |
| GSE121212 | hallmark_hs | 50 | 1000 | 8.9e-16 | 50 | 0 |
| GSE121212 | gobp_hs | 7229 | 250 | 2.7e-15 | 7227 | 2 |
| GSE121212 | react_hs | 1687 | 500 | 1.7e-15 | 1687 | 0 |

All remaining differences were traced (instrumented copy of `gseaPhenoGen` run in the test harness only) to a single null permutation per set in which the maximum and the negative minimum of the running sum are equal to rounding, so the sign is resolved differently by the closed-form and the incremental sums. NES then differs by at most 0.0047 and FDR shifts slightly through the pooled null. This is not a defect.

(b) Independent R implementation (own permutations, `set.seed(1)`, limma's `s2.prior` and `df.prior` held fixed, ES by the fgsea formula, Subramanian 2005 NES and FDR). Agreement is distributional:

| Dataset | Library | Sets | Permutations | Pearson r (NES) | max abs dNES | sets with p outside 4 MC SE | FDR<=0.05 app / R | FDR<=0.25 app / R (Jaccard) |
|---|---|---|---|---|---|---|---|---|
| GSE54456 | hallmark_hs | 50 | 1000 | 0.9997 | 0.052 | 0 | 0 / 0 | 23 / 23 (1.00) |
| GSE54456 | gobp_hs | 7149 | 250 | 0.9995 | 0.198 | 0 | 0 / 0 | 1321 / 1431 (0.86) |
| GSE54456 | react_hs | 1667 | 500 | 0.9996 | 0.114 | 0 | 0 / 0 | 404 / 387 (0.91) |
| GSE186063 | hallmark_hs | 50 | 1000 | 0.9998 | 0.065 | 0 | 0 / 0 | 19 / 19 (1.00) |
| GSE186063 | gobp_hs | 7213 | 250 | 0.9994 | 0.162 | 1 | 0 / 0 | 0 / 0 (0.00) |
| GSE186063 | react_hs | 1685 | 500 | 0.9996 | 0.113 | 0 | 0 / 0 | 248 / 271 (0.87) |
| GSE121212 | hallmark_hs | 50 | 1000 | 0.9998 | 0.080 | 0 | 0 / 0 | 19 / 19 (1.00) |
| GSE121212 | gobp_hs | 7229 | 250 | 0.9993 | 0.193 | 1 | 0 / 0 | 0 / 0 (0.00) |
| GSE121212 | react_hs | 1687 | 500 | 0.9994 | 0.136 | 0 | 0 / 0 | 0 / 0 (0.00) |

The method is conservative on these data: no set reached FDR <= 0.05 in any of the nine runs, and the app and the reference agree on this. On GSE54456 Hallmark the leading sets have FDR 0.0557 (section 7). The app states that sample permutation is conservative and quotes FDR 0.25.

Gene-permutation fallback:

### 4.6 Enrichment across contrasts

| Dataset | Library | Contrasts | Sets | cells FDR<=0.05 (app = R pairwise fry) | same, full-design fry |
|---|---|---|---|---|---|
| GSE63310 | hallmark_mm | 3 | 50 | 62 = 62 | 61 |
| GSE63310 | react_mm | 3 | 1216 | 1627 = 1627 | 1784 |
| GSE63310 | gobp_mm | 3 | 7172 | 10481 = 10481 | 11341 |
| GSE121212 | hallmark_hs | 15 | 50 | 473 = 473 | 471 |
| GSE121212 | react_hs | 15 | 1687 | 14691 = 14691 | 14856 |
| GSE121212 | gobp_hs | 15 | 7229 | 58259 = 58259 | 59902 |

FDR (BH within contrast), direction and the plotted signed -log10 FDR of the top 25 sets agree exactly with `limma::fry` on each pair's samples (maximum relative FDR deviation 4.8e-10; 0 direction differences; plotted values identical).

### 4.7 Heatmap cluster enrichment

For every k-means cluster of every heatmap case (top 50 and 200 variable or DE genes, five datasets, Hallmark), set sizes, overlaps, hypergeometric p and within-cluster BH FDR agree with `scipy.stats.hypergeom.sf` to 8.4e-11. The display rule (sets with at least two overlapping genes; significant sets, else the top five, labelled "(none at FDR <= 0.05: top sets shown)") was confirmed in the code. The universe concern applies.

### 4.8 Discovery

- Family stems (`familyStems`) and HGNC curated groups (`curatedGroupRows`) are identical to independent Python implementations on all five datasets; candidate lists (families + HGNC + library, size 3-100) are identical in all 18 runs (15 two-group runs and 3 all-group runs; for example 9,341 candidates in GSE54456 with GO BP).
- Module effect (difference of mean per-gene z between groups, population SD) agrees to 1.1e-15.
- Per-module counts of individually significant genes (limma-trend q <= 0.05), the nominal-tier count (p <= 0.01 when no module passes FDR), the hidden-module count (< 25% individually significant) and the chance expectation stated in the summary text were recomputed and agree in all 15 runs.
- Leave-one-out direction text for the top 20 modules agrees with a recomputation on GSE54456, GSE41745 and GSE63310.
- Permutation path (more than two groups): F statistics agree with a one-way F on module scores to 2.1e-12; permutation p agrees with an independent R random-set null within 4 Monte-Carlo SE for every module (maximum z 3.26).

## 5. Calibration (label shuffles)

Group labels were permuted before analysis (three seeds per dataset, mulberry32 seeds 11-13; families + HGNC + Hallmark). No module reached FDR <= 0.05 in any of the 12 runs; the numbers at p <= 0.05 were 0.42-1.79 times the chance expectation.

| Dataset | Seed | Groups | Path | Candidates | p<=0.01 (expected) | p<=0.05 (expected) | FDR<=0.05 |
|---|---|---|---|---|---|---|---|
| GSE121212 | 11 | 6 | permutation | 3009 | 22 (30.1) | 159 (150.5) | 0 |
| GSE121212 | 12 | 6 | permutation | 3009 | 31 (30.1) | 107 (150.5) | 0 |
| GSE121212 | 13 | 6 | permutation | 3009 | 31 (30.1) | 144 (150.5) | 0 |
| GSE63310 | 11 | 3 | permutation | 2671 | 18 (26.7) | 106 (133.6) | 0 |
| GSE63310 | 12 | 3 | permutation | 2671 | 29 (26.7) | 111 (133.6) | 0 |
| GSE63310 | 13 | 3 | permutation | 2671 | 15 (26.7) | 103 (133.6) | 0 |
| GSE54456 | 11 | 2 | CAMERA | 2933 | 0 (29.3) | 61 (146.7) | 0 |
| GSE54456 | 12 | 2 | CAMERA | 2933 | 16 (29.3) | 129 (146.7) | 0 |
| GSE54456 | 13 | 2 | CAMERA | 2933 | 13 (29.3) | 105 (146.7) | 0 |
| GSE186063 | 11 | 3 | permutation | 2996 | 40 (30.0) | 268 (149.8) | 0 |
| GSE186063 | 12 | 3 | permutation | 2996 | 12 (30.0) | 111 (149.8) | 0 |
| GSE186063 | 13 | 3 | permutation | 2996 | 37 (30.0) | 122 (149.8) | 0 |

CAMERA on the shuffled GSE54456 labels, same candidates:

| Seed | Candidates | estimated correlation (app): p<=0.01 / FDR<=0.05 | fixed 0.01 (limma default): p<=0.01 / FDR<=0.05 |
|---|---|---|---|
| 11 | 2933 | 0 / 0 | 204 / 97 |
| 12 | 2933 | 16 / 0 | 201 / 78 |
| 13 | 2933 | 13 / 0 | 161 / 84 |

## 6. Heatmap clustering

### 6.1 Gene selection, z-scores and hierarchical clustering

For every dataset and N = 50 and 200: the top-variable selection equals the N largest population variances over the included samples, and the top-DE selection equals the N smallest limma-trend p-values; row z-scores (SD with n - 1) agree to 1e-12 with the plotted matrix. For all 20 hierarchical cases, the displayed row order is a valid leaf order of `hclust(as.dist(1 - cor(t(Z))), method = "average")`: every one of the 2,480 merges occupies a contiguous block of rows. The app does not expose merge heights, so heights were not compared.

### 6.2 k-means and auto-k

| Dataset | N genes | Source | auto-k (app) | app mean silhouette | singletons at k | cluster::silhouette, same partition | k by cluster::silhouette | app/R within-SS at k | ARI vs R at k | auto-k, rows reversed |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE121212 | 200 | de | 2 | 0.6743 | 0 | 0.6743 | 2 | 1.0000 | 1.000 | 2 |
| GSE121212 | 200 | var | 8 | 0.4042 | 1 | 0.3992 | 8 | 1.0834 | 0.915 | 8 |
| GSE121212 | 50 | de | 2 | 0.7091 | 0 | 0.7091 | 2 | 1.0000 | 1.000 | 2 |
| GSE121212 | 50 | var | 5 | 0.6671 | 0 | 0.6671 | 5 | 1.0000 | 1.000 | 6 |
| GSE186063 | 200 | de | 2 | 0.6762 | 0 | 0.6762 | 2 | 1.0000 | 1.000 | 2 |
| GSE186063 | 200 | var | 7 | 0.5074 | 0 | 0.5074 | 7 | 1.0248 | 0.981 | 7 |
| GSE186063 | 50 | de | 2 | 0.7339 | 0 | 0.7339 | 2 | 1.0000 | 1.000 | 2 |
| GSE186063 | 50 | var | 6 | 0.7705 | 1 | 0.7505 | 6 | 1.0000 | 1.000 | 5 |
| GSE41745 | 200 | de | 2 | 0.8435 | 0 | 0.8435 | 2 | 1.0000 | 1.000 | 2 |
| GSE41745 | 200 | var | 4 | 0.5484 | 0 | 0.5484 | 4 | 1.0000 | 1.000 | 4 |
| GSE41745 | 50 | de | 2 | 0.3623 | 0 | 0.3623 | 2 | 1.0000 | 1.000 | 2 |
| GSE41745 | 50 | var | 5 | 0.7099 | 1 | 0.6899 | 5 | 1.0000 | 1.000 | 4 |
| GSE54456 | 200 | de | 2 | 0.7702 | 0 | 0.7702 | 2 | 1.0000 | 1.000 | 2 |
| GSE54456 | 200 | var | 7 | 0.5455 | 0 | 0.5455 | 7 | 1.0394 | 0.957 | 6 |
| GSE54456 | 50 | de | 2 | 0.8460 | 1 | 0.8260 | 2 | 1.0000 | 1.000 | 2 |
| GSE54456 | 50 | var | 4 | 0.7314 | 0 | 0.7314 | 4 | 1.0000 | 1.000 | 4 |
| GSE63310 | 200 | de | 2 | 0.8513 | 0 | 0.8513 | 2 | 1.0000 | 1.000 | 2 |
| GSE63310 | 200 | var | 2 | 0.8400 | 0 | 0.8400 | 2 | 1.0000 | 1.000 | 2 |
| GSE63310 | 50 | de | 2 | 0.8770 | 0 | 0.8770 | 2 | 1.0000 | 1.000 | 2 |
| GSE63310 | 50 | var | 3 | 0.8777 | 1 | 0.8577 | 2 | 1.0000 | 1.000 | 2 |

Mean silhouettes agree with `cluster::silhouette` to 3.8e-9 when no cluster is a singleton.

## 7. Positive controls

Ranks are positions in the app's sorted tables (ORA by p; FRY by directional p; GSEA by FDR then |NES|; Discovery by p over all candidates). Discovery candidates are limited to 3-100 genes, so most Hallmark sets are not candidates. The seven Hallmark interferon, inflammatory, IL-6/JAK/STAT3, TNF-alpha/NF-kB, E2F and G2M sets were within the top 8 of 50 by ORA in GSE54456, GSE186063 and GSE121212 (top 13 in GSE41745); by FRY and GSEA they were within the top 10 in GSE54456 and within the top 14 (GSE186063), 20 (GSE121212) and 22 (GSE41745) elsewhere. Keratinization and cornified-envelope sets rank high by ORA of up-regulated genes (Reactome keratinization 3rd-23rd of about 1,650) but low by directional FRY (747th-1,140th), because these sets contain strongly up- and strongly down-regulated keratinocyte genes. The mixed FRY p ranks Reactome keratinization 44th, 65th, 259th and 151st (GSE54456, GSE186063, GSE121212, GSE41745), with mixed FDR <= 1.4e-3 in all four. GO IL-17 signaling ranked low by all methods in all datasets. Discovery placed no positive control below FDR 0.25.

| Dataset | Gene set | ORA up rank | FRY rank | GSEA rank (FDR, scheme) | Discovery rank / candidates (FDR) | CAMERA 0.01 rank |
|---|---|---|---|---|---|---|
| GSE54456 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 3/50 | 3/50 | 6 (0.0557, samples) | 25/2933 (0.473) | 1 |
| GSE54456 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 1/50 | 1/50 | 4 (0.0557, samples) | not a candidate (>100 genes) |  |
| GSE54456 | HALLMARK_INFLAMMATORY_RESPONSE | 5/50 | 4/50 | 1 (0.0557, samples) | not a candidate (>100 genes) |  |
| GSE54456 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | 8/50 | 10/50 | 8 (0.0588, samples) | not a candidate (>100 genes) |  |
| GSE54456 | HALLMARK_IL6_JAK_STAT3_SIGNALING | 7/50 | 2/50 | 3 (0.0557, samples) | 17/2933 (0.473) | 6 |
| GSE54456 | HALLMARK_E2F_TARGETS | 4/50 | 6/50 | 5 (0.0557, samples) | not a candidate (>100 genes) |  |
| GSE54456 | HALLMARK_G2M_CHECKPOINT | 2/50 | 5/50 | 2 (0.0557, samples) | not a candidate (>100 genes) |  |
| GSE54456 | GOBP_KERATINIZATION | 46/7149 | 814/7149 | 237 (0.11, samples) | 2086/9341 (0.528) | 61 |
| GSE54456 | GOBP_CORNIFIED_ENVELOPE_ASSEMBLY | 2545/7149 | 3224/7149 | 2198 (0.37, samples) | 4164/9341 (0.674) | 4528 |
| GSE54456 | GOBP_INTERLEUKIN_17_MEDIATED_SIGNALING_PATHWAY | 1756/7149 | 941/7149 | 934 (0.193, samples) | 2043/9341 (0.528) | 3330 |
| GSE54456 | GOBP_RESPONSE_TO_TYPE_I_INTERFERON | 270/7149 | 123/7149 | 51 (0.102, samples) | 59/9341 (0.256) | 20 |
| GSE54456 | GOBP_RESPONSE_TO_TYPE_II_INTERFERON | 243/7149 | 105/7149 | 85 (0.106, samples) | not a candidate (>100 genes) |  |
| GSE54456 | GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE | 28/7149 | 110/7149 | 12 (0.102, samples) | 257/9341 (0.346) | 40 |
| GSE54456 | GOBP_MITOTIC_CELL_CYCLE | 57/7149 | 904/7149 | 451 (0.117, samples) | not a candidate (>100 genes) |  |
| GSE54456 | REACTOME_KERATINIZATION | 23/1667 | 994/1667 | 595 (0.383, samples) | not a candidate (>100 genes) |  |
| GSE54456 | REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE | 8/1667 | 290/1667 | 85 (0.162, samples) | not a candidate (>100 genes) |  |
| GSE54456 | REACTOME_INTERLEUKIN_17_SIGNALING | 750/1667 | 187/1667 | 55 (0.162, samples) | 348/4456 (0.5) | 328 |
| GSE54456 | REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | 13/1667 | 74/1667 | 183 (0.188, samples) | 89/4456 (0.42) | 13 |
| GSE54456 | REACTOME_INTERFERON_GAMMA_SIGNALING | 18/1667 | 53/1667 | 108 (0.164, samples) | 79/4456 (0.42) | 63 |
| GSE54456 | REACTOME_TNF_SIGNALING | 798/1667 | 317/1667 | 418 (0.256, samples) | 405/4456 (0.509) | 472 |
| GSE54456 | REACTOME_CELL_CYCLE | 6/1667 | 115/1667 | 35 (0.162, samples) | not a candidate (>100 genes) |  |
| GSE54456 | REACTOME_ANTIMICROBIAL_PEPTIDES | 3/1667 | 1/1667 | 3 (0.162, samples) | 1/4456 (0.357) | 9 |
| GSE186063 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 2/50 | 5/50 | 6 (0.108, samples) | 50/2996 (0.732) | 1 |
| GSE186063 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 1/50 | 4/50 | 2 (0.108, samples) | not a candidate (>100 genes) |  |
| GSE186063 | HALLMARK_INFLAMMATORY_RESPONSE | 3/50 | 14/50 | 12 (0.13, samples) | not a candidate (>100 genes) |  |
| GSE186063 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | 6/50 | 12/50 | 13 (0.13, samples) | not a candidate (>100 genes) |  |
| GSE186063 | HALLMARK_IL6_JAK_STAT3_SIGNALING | 8/50 | 9/50 | 14 (0.133, samples) | 125/2996 (0.732) | 10 |
| GSE186063 | HALLMARK_E2F_TARGETS | 7/50 | 1/50 | 1 (0.108, samples) | not a candidate (>100 genes) |  |
| GSE186063 | HALLMARK_G2M_CHECKPOINT | 5/50 | 2/50 | 3 (0.108, samples) | not a candidate (>100 genes) |  |
| GSE186063 | GOBP_KERATINIZATION | 35/7213 | 1241/7213 | 748 (0.45, samples) | 1652/9395 (0.663) | 44 |
| GSE186063 | GOBP_CORNIFIED_ENVELOPE_ASSEMBLY | 2021/7213 | 811/7213 | 946 (0.467, samples) | 843/9395 (0.658) | 1769 |
| GSE186063 | GOBP_INTERLEUKIN_17_MEDIATED_SIGNALING_PATHWAY | 1354/7213 | 4688/7213 | 1344 (0.496, samples) | 6025/9395 (0.832) | 5059 |
| GSE186063 | GOBP_RESPONSE_TO_TYPE_I_INTERFERON | 251/7213 | 880/7213 | 548 (0.437, samples) | 844/9395 (0.658) | 94 |
| GSE186063 | GOBP_RESPONSE_TO_TYPE_II_INTERFERON | 104/7213 | 310/7213 | 390 (0.437, samples) | not a candidate (>100 genes) |  |
| GSE186063 | GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE | 16/7213 | 363/7213 | 177 (0.437, samples) | not a candidate (>100 genes) |  |
| GSE186063 | GOBP_MITOTIC_CELL_CYCLE | 115/7213 | 369/7213 | 269 (0.437, samples) | not a candidate (>100 genes) |  |
| GSE186063 | REACTOME_KERATINIZATION | 17/1685 | 1140/1685 | 715 (0.648, samples) | not a candidate (>100 genes) |  |
| GSE186063 | REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE | 10/1685 | 662/1685 | 513 (0.455, samples) | not a candidate (>100 genes) |  |
| GSE186063 | REACTOME_INTERLEUKIN_17_SIGNALING | 200/1685 | 642/1685 | 503 (0.454, samples) | 1262/4516 (0.694) | 432 |
| GSE186063 | REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | 9/1685 | 193/1685 | 191 (0.217, samples) | 252/4516 (0.601) | 28 |
| GSE186063 | REACTOME_INTERFERON_GAMMA_SIGNALING | 20/1685 | 292/1685 | 245 (0.247, samples) | 388/4516 (0.601) | 64 |
| GSE186063 | REACTOME_TNF_SIGNALING | 944/1685 | 390/1685 | 442 (0.432, samples) | 471/4516 (0.601) | 353 |
| GSE186063 | REACTOME_CELL_CYCLE | 7/1685 | 15/1685 | 11 (0.215, samples) | not a candidate (>100 genes) |  |
| GSE186063 | REACTOME_ANTIMICROBIAL_PEPTIDES | 3/1685 | 28/1685 | 34 (0.215, samples) | 48/4516 (0.601) | 37 |
| GSE121212 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 4/50 | 18/50 | 17 (0.217, samples) | 472/3009 (0.713) | 3 |
| GSE121212 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 1/50 | 14/50 | 12 (0.217, samples) | not a candidate (>100 genes) |  |
| GSE121212 | HALLMARK_INFLAMMATORY_RESPONSE | 3/50 | 16/50 | 15 (0.217, samples) | not a candidate (>100 genes) |  |
| GSE121212 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | 5/50 | 20/50 | 18 (0.217, samples) | not a candidate (>100 genes) |  |
| GSE121212 | HALLMARK_IL6_JAK_STAT3_SIGNALING | 8/50 | 6/50 | 7 (0.211, samples) | 210/3009 (0.713) | 8 |
| GSE121212 | HALLMARK_E2F_TARGETS | 7/50 | 3/50 | 2 (0.211, samples) | not a candidate (>100 genes) |  |
| GSE121212 | HALLMARK_G2M_CHECKPOINT | 6/50 | 4/50 | 3 (0.211, samples) | not a candidate (>100 genes) |  |
| GSE121212 | GOBP_KERATINIZATION | 73/7229 | 1186/7229 | 205 (0.335, samples) | 2092/9410 (0.654) | 583 |
| GSE121212 | GOBP_CORNIFIED_ENVELOPE_ASSEMBLY | 5466/7229 | 2609/7229 | 1388 (0.366, samples) | 3860/9410 (0.682) | 4581 |
| GSE121212 | GOBP_INTERLEUKIN_17_MEDIATED_SIGNALING_PATHWAY | 2032/7229 | 3826/7229 | 2580 (0.508, samples) | 4908/9410 (0.744) | 5445 |
| GSE121212 | GOBP_RESPONSE_TO_TYPE_I_INTERFERON | 183/7229 | 340/7229 | 184 (0.335, samples) | 423/9410 (0.654) | 16 |
| GSE121212 | GOBP_RESPONSE_TO_TYPE_II_INTERFERON | 115/7229 | 689/7229 | 159 (0.335, samples) | not a candidate (>100 genes) |  |
| GSE121212 | GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE | 29/7229 | 756/7229 | 233 (0.335, samples) | not a candidate (>100 genes) |  |
| GSE121212 | GOBP_MITOTIC_CELL_CYCLE | 205/7229 | 1131/7229 | 417 (0.335, samples) | not a candidate (>100 genes) |  |
| GSE121212 | REACTOME_KERATINIZATION | 20/1687 | 1103/1687 | 473 (0.359, samples) | not a candidate (>100 genes) |  |
| GSE121212 | REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE | 11/1687 | 580/1687 | 196 (0.293, samples) | not a candidate (>100 genes) |  |
| GSE121212 | REACTOME_INTERLEUKIN_17_SIGNALING | 668/1687 | 311/1687 | 363 (0.323, samples) | 704/4515 (0.642) | 272 |
| GSE121212 | REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | 9/1687 | 501/1687 | 329 (0.318, samples) | 885/4515 (0.642) | 70 |
| GSE121212 | REACTOME_INTERFERON_GAMMA_SIGNALING | 31/1687 | 528/1687 | 448 (0.356, samples) | 1037/4515 (0.642) | 170 |
| GSE121212 | REACTOME_TNF_SIGNALING | 1058/1687 | 626/1687 | 632 (0.421, samples) | 1228/4515 (0.646) | 515 |
| GSE121212 | REACTOME_CELL_CYCLE | 7/1687 | 148/1687 | 31 (0.293, samples) | not a candidate (>100 genes) |  |
| GSE121212 | REACTOME_ANTIMICROBIAL_PEPTIDES | 13/1687 | 182/1687 | 18 (0.293, samples) | 386/4515 (0.642) | 142 |
| GSE41745 | HALLMARK_INTERFERON_ALPHA_RESPONSE | 2/50 | 3/50 | 21 (0.00342, genes) | 2/2603 (0.99) | 1 |
| GSE41745 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 1/50 | 4/50 | 1 (0.00342, genes) | not a candidate (>100 genes) |  |
| GSE41745 | HALLMARK_INFLAMMATORY_RESPONSE | 4/50 | 2/50 | 12 (0.00342, genes) | not a candidate (>100 genes) |  |
| GSE41745 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | 7/50 | 20/50 | 5 (0.00342, genes) | not a candidate (>100 genes) |  |
| GSE41745 | HALLMARK_IL6_JAK_STAT3_SIGNALING | 5/50 | 1/50 | 22 (0.00342, genes) | 1/2603 (0.99) | 10 |
| GSE41745 | HALLMARK_E2F_TARGETS | 13/50 | 13/50 | 9 (0.00342, genes) | not a candidate (>100 genes) |  |
| GSE41745 | HALLMARK_G2M_CHECKPOINT | 9/50 | 7/50 | 7 (0.00342, genes) | not a candidate (>100 genes) |  |
| GSE41745 | GOBP_KERATINIZATION | 5/6816 | 716/6816 | 556 (0.0591, genes) | 238/9023 (0.87) | 1 |
| GSE41745 | GOBP_CORNIFIED_ENVELOPE_ASSEMBLY | 1790/6816 | 2832/6816 | 2090 (0.191, genes) | 3009/9023 (0.87) | 1800 |
| GSE41745 | GOBP_INTERLEUKIN_17_MEDIATED_SIGNALING_PATHWAY | 2083/6816 | 5842/6816 | 6636 (0.959, genes) | 8682/9023 (0.995) | 8727 |
| GSE41745 | GOBP_RESPONSE_TO_TYPE_I_INTERFERON | 70/6816 | 387/6816 | 631 (0.0591, genes) | 34/9023 (0.87) | 30 |
| GSE41745 | GOBP_RESPONSE_TO_TYPE_II_INTERFERON | 52/6816 | 4/6816 | 315 (0.0591, genes) | 3/9023 (0.87) | 14 |
| GSE41745 | GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE | 10/6816 | 62/6816 | 475 (0.0591, genes) | 8/9023 (0.87) | 2 |
| GSE41745 | GOBP_MITOTIC_CELL_CYCLE | 258/6816 | 798/6816 | 76 (0.0591, genes) | not a candidate (>100 genes) |  |
| GSE41745 | REACTOME_KERATINIZATION | 3/1618 | 747/1618 | 105 (0.0301, genes) | not a candidate (>100 genes) |  |
| GSE41745 | REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE | 2/1618 | 311/1618 | 201 (0.0301, genes) | not a candidate (>100 genes) |  |
| GSE41745 | REACTOME_INTERLEUKIN_17_SIGNALING | 362/1618 | 175/1618 | 462 (0.0459, genes) | 221/4120 (0.864) | 400 |
| GSE41745 | REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | 9/1618 | 71/1618 | 262 (0.0301, genes) | 14/4120 (0.864) | 2 |
| GSE41745 | REACTOME_INTERFERON_GAMMA_SIGNALING | 10/1618 | 53/1618 | 157 (0.0301, genes) | 9/4120 (0.864) | 36 |
| GSE41745 | REACTOME_TNF_SIGNALING | 585/1618 | 378/1618 | 227 (0.0301, genes) | 274/4120 (0.864) | 307 |
| GSE41745 | REACTOME_CELL_CYCLE | 18/1618 | 136/1618 | 8 (0.0301, genes) | not a candidate (>100 genes) |  |
| GSE41745 | REACTOME_ANTIMICROBIAL_PEPTIDES | 8/1618 | 13/1618 | 360 (0.0301, genes) | 1/4120 (0.864) | 3 |

## 8. Not assessed

- ORA sources other than DE lists: pasted gene lists and Pattern-tab clusters.
- Gene-set analyses with voom or Welch as the DE method, and ORA/GSEA with covariates.
- User-uploaded GMT parsing; GO Cellular Component, GO Molecular Function and rat libraries.
- Plotted values of the GSEA running-score curve, barcode, ECDF, ridgeline, Jaccard and leading-edge plots, and the CSV exports of the enrichment, Discovery and heatmap tables.
- Heatmap with batch-adjusted display (`dispM` with BATCHDISP), heatmap from volcano/MA selections, and hierarchical merge heights.
- Discovery for mouse with HGNC groups (matched by upper-cased symbol): mapping was checked, biological correctness of the orthology approximation was not.
- Discovery more-than-two-group path with GO BP or Reactome as library (Hallmark only).
- Behaviour in real browsers.

## 9. Files

- `VAL_ENR_checks.csv`: one row per check (name, dataset, reference, values compared, maximum deviation, tolerance, pass, note).
- `VAL_ENR_positive_controls.csv`: ranks and FDR of the positive-control sets per method and dataset.
- `VAL_ENR_runner.py`: runs the app headless, exports matrices and results, calls the R scripts, and writes a checks CSV.
- `VAL_ENR_reference.R`: upstream log-CPM, limma-trend, set mapping, phyper ORA (two universes), fry, camera, calcGseaStat, fgseaSimple, independent sample-permutation GSEA.
- `VAL_ENR_fry_cov.R`, `VAL_ENR_enrmat.R`, `VAL_ENR_discovery.R`, `VAL_ENR_heatmap.R`: FRY with covariates, enrichment matrix, Discovery CAMERA and permutation null, heatmap hclust/kmeans/silhouette.
- `VAL_ENR_gsea_replay.py`: exact replay of the seeded sample-permutation GSEA.

Runner status: `VAL_ENR_runner.py` was executed end to end from the kit files in a clean work directory. It completed library decoding, all core app runs, all R references, the seeded replay, the covariate FRY checks and the GSE63310 enrichment-matrix reference; all 148 reference CSV files it produced are byte-identical to those used in this report. It was stopped during the GSE121212 GO BP enrichment-matrix reference (15 contrasts, run with and without the full design) because the shared machine was heavily loaded (load average about 107 on 14 cores), so its Discovery and heatmap stages were not executed by the runner itself; those stages were executed in this task with the same JavaScript bodies and R scripts, and their results are those reported here.
