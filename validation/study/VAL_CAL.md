# benchsiDE calibration: false-positive control and power (VAL_CAL)

Build under test: `benchside/index.html`, APP_VERSION 0.22.0-beta (2026-09-30), sha256 `17511676c3096bda2f4ba964dea8ce5b78faa300a29ad0d1007c210da7c7ae65`. All application results were produced by the application's own functions, executed in macOS JavaScriptCore through the validation-kit harness (`harness.py`, DOM/Plotly stub); no statistic was re-implemented on the application side. Count matrices were checked against the MANIFEST sha256 values. References: R 4 with edgeR 4.8.2, limma 3.66.0, fgsea 1.36.2, metafor 5.0.1; Python with numpy/scipy. The index.html file was not modified.

## 1. Scope

an earlier internal audit calibrated the moderated t, voom, preranked and sample-permutation GSEA, FRY, CAMERA and Discovery on GSE54456 normal skin only (5v5, 10v10, 40v40; three seeds each), the meta-analysis on four 10v10 null datasets, and signature transfer on ten null splits. This report extends that work in five directions: (i) four null pools from three studies (GSE54456 normal skin, GSE121212 healthy controls, GSE121212 psoriasis non-lesional skin, GSE186063 non-lesional skin); (ii) the smallest designs in use (3 v 3), with 40 seeds per cell for gene-level tests; (iii) every inferential output, including ORA, co-expression FDR, the three-group moderated F and Discovery path, paired designs, meta-analysis with small studies, consensus, RRHO and signature transfer; (iv) R references run on the same splits, so that each miscalibration is classified as a property of the method or of the implementation; (v) positive controls on all five human datasets, including GSE41745 (2 lesional v 3 non-lesional), and on subsamples of 2 to 10 per group.

## 2. Methods

**Null pools.** GSE54456 normal_skin (81 samples), GSE121212 CTRL_healthy (37), GSE121212 PSO_non_lesional (27), GSE186063 non-lesion (28; psoriasis and psoriatic-arthritis patients). For each split, 2n samples were drawn without replacement from one pool, assigned to two groups of n, and loaded as a session containing only those samples (filterByExpr, TMM and log-CPM therefore computed on the split, as a user with those samples would do). Sizes: 3v3, 5v5, 10v10 and the largest balanced split (40v40, 18v18, 13v13, 14v14). Gene-level tests (moderated t, voom, Welch, moderated t with a random continuous covariate): 40 seeds per pool and size (640 sessions). Gene-set, Discovery and co-expression outputs: 10 seeds at 3v3, 8 at 5v5, 5 at the larger sizes (112 sessions). Seeds are fixed (numpy `default_rng`); sample lists are stored with every result.

**Outputs and settings.** Application defaults unless stated: filterByExpr, TMM, FDR 0.05, |log2FC| 1, gene-set size 5-2000, GSEA permutation mode `auto`. Each session ran: DE (mod, voom, welch; mod adjusted for a random N(0,1) covariate and for a random balanced binary covariate); for the built-in Hallmark (50 sets) and Reactome (1,736 sets) libraries: ORA of DE up, down and all hits, ORA of a pasted list of the 200 genes with the smallest moderated-t p, FRY, GSEA in `auto` mode, GSEA forced to gene permutation, and CAMERA. CAMERA is exposed only through the Discovery screen; it was run there with only the gene-set library as source and the enrichment size limits (5-2000). Discovery was also run with its defaults (gene families + HGNC groups, sizes 3-100). Co-expression: for 4 randomly chosen query genes per session, the query gene's values were permuted across samples 10 times and `drawCoexp` was re-run, so every correlation with the query is null. Three-group null: 3 groups of n from one pool (8 seeds), `computeModF` and Discovery (random-set permutation path, used when there are more than two groups). Paired null: GSE121212 psoriasis patients with both biopsies (26 pairs); within each patient the two samples were assigned at random to groups X and Y, and the patient was entered as covariate (3, 5, 10 and 26 pairs; 10, 10, 8, 8 seeds). Meta-analysis and consensus: 3, 4 or 5 null datasets (session plus comparison datasets added with `analyzeDatasetB`, disjoint samples when a pool is used twice), each 3v3, 5v5 or 10v10, 6 seeds; `runMeta` (REM default, REM with the prediction-interval criterion off, FEM) and `runMulti`. The same sessions supplied RRHO maps (`runCompare`, session v first comparison dataset) and null-signature transfers (the top three Discovery modules of the null session scored in every null comparison dataset). A real psoriasis signature (GSE54456 psoriasis v normal, 1,293 up and 1,513 down genes) was transferred to 72 null splits and to four lesional contrasts.

**Positive controls.** GSE54456 (Psoriasis_skin 90 v normal_skin 81), GSE121212 PSO_lesional v PSO_non_lesional (27 v 27; unpaired and paired by Patient) and v CTRL_healthy (27 v 37), GSE186063 lesion v non-lesion (26 v 28; unpaired and paired by the inferred Pair), GSE83645 psoriasis v uninvolved (20 v 5; unpaired and paired by Patient), GSE41745 lesional v non_lesional (2 v 3; unpaired and paired by Patient). Recovery was scored for the 20 genes used in an earlier internal audit (19 induced; KRT77, which is reduced in lesional skin, scored as down) and for 7 Hallmark and 8 Reactome sets expected to be up-regulated in lesional skin (interferon-gamma and -alpha, inflammatory response, IL6-JAK-STAT3, TNF-alpha via NF-kB, E2F, G2M; cornified envelope, keratinization, antimicrobial peptides, interferon signalling (three sets), metal sequestration by antimicrobial proteins, cell-cycle checkpoints). Low-n behaviour: GSE54456 subsamples 2v2, 2v3, 3v3, 5v5, 10v10 and paired subsamples of 3 and 5 patients from GSE121212 and GSE186063 (5 seeds each).

**Criteria.** Under a global null, BH at 0.05 should give at least one discovery in at most about 5% of splits. A cell fails when the number of splits with at least one discovery exceeds the 99th percentile of Binomial(R, 0.05) (R = splits in the cell); 95% Clopper-Pearson intervals are reported. Mean per-gene fraction p < 0.05 fails when mean - 3 SE (across seeds) exceeds 0.05. Parity checks use tolerance 0 for counts and 1e-8 (absolute) or 1e-6 (relative, p-values) for continuous values. Power checks require at least 75% of the classic genes that pass the filter, or at least 5 expected sets, at FDR 0.05.

**Reference runs.** `cal_reference.R` repeated every one of the 112 gene-set sessions and the 160 3v3 gene-level sessions in R (same samples, filterByExpr, TMM, `cpm(log=TRUE, prior.count=2)`, `eBayes(trend=TRUE)`, voom, `t.test`, `fry`, `camera(inter.gene.cor=NA)`, `fgsea` on the same moderated-t ranking), and recomputed the meta-analysis p-values of four null combinations with `metafor::rma(method="DL")` from the per-dataset estimates exported by the app. `robust_check.R` tested `eBayes(robust=TRUE)` on the 160 3v3 splits.

## 3. Parity of the calibration runs with R

| Check | Values compared | Max deviation | Tolerance | Pass |
|---|---|---|---|---|
| Null-split DE: genes at FDR<=0.05, moderated t vs limma-trend | 272 | 0 | 0 | yes |
| Null-split DE: genes at FDR<=0.05, voom vs limma voom | 272 | 0 | 0 | yes |
| Null-split DE: genes at FDR<=0.05, Welch vs stats::t.test | 272 | 0 | 0 | yes |
| Null-split DE: fraction p<0.05 (all genes), moderated t vs limma-trend | 272 | 5e-16 | 1e-08 | yes |
| Null-split DE: fraction p<0.05, Welch vs t.test | 272 | 4.86e-16 | 1e-08 | yes |
| Null-split filterByExpr kept-gene count | 112 | 0 | 0 | yes |
| Null-split FRY Hallmark: sets at FDR<=0.05 and at p<=0.05 | 224 | 0 | 0 | yes |
| Null-split FRY Reactome: sets at FDR<=0.05 and at p<=0.05 | 224 | 0 | 0 | yes |
| Null-split CAMERA (Discovery, library source) Hallmark: sets at FDR<=0.05 and p<=0.05 | 224 | 0 | 0 | yes |
| Null-split CAMERA (Discovery, library source) Reactome: sets at FDR<=0.05 and p<=0.05 | 224 | 0 | 0 | yes |
| Null-split testable set counts (size 5-2000) | 224 | 0 | 0 | yes |
| Meta-analysis REM p-values vs metafor, all genes of 4 null combinations (relative) | 68732 | 2.2e-13 | 1e-06 | yes |
| Meta-analysis REM genes at meta FDR<=0.05 vs metafor+BH (4 null combinations) | 4 | 0 | 0 | yes |
| Preranked GSEA (gene permutation) false-positive behaviour vs fgsea, Hallmark | 112 | 0.00589 | 0.02 | yes |

Every FDR count, kept-gene count and set count produced by the app on the null splits equals the R result exactly, and the meta-analysis p-values agree with metafor to 2.2e-13 (relative) on all 68,732 gene-level tests. Preranked GSEA behaves as fgsea does on the same rankings. All miscalibration reported below is therefore a property of the methods and defaults, not of the JavaScript implementation.

## 4. Null calibration

### 4.1 Gene-level tests (40 random splits per cell)

Entries: splits with at least one gene at FDR <= 0.05 out of 40 (expected by chance: at most 2 of 40; failure threshold 6) and the mean fraction of genes with p < 0.05 (nominal 0.05). Values in bold exceed the failure threshold.

| Pool | n per group | mod | voom | Welch | mod + random covariate | mean frac p<0.05 (mod / voom / Welch) | frac p<0.001 (mod; nominal 0.001) |
|---|---|---|---|---|---|---|---|
| GSE54456_normal | 3 | **10/40** (2837) | **10/40** (2875) | 0/40 (0) | **8/40** (2673) | 0.0751 / 0.0740 / 0.0461 | 0.00271 |
| GSE54456_normal | 5 | 0/40 (0) | 0/40 (0) | 1/40 (1) | 3/40 (42) | 0.0522 / 0.0518 / 0.0448 | 0.00105 |
| GSE54456_normal | 10 | 0/40 (0) | 0/40 (0) | 1/40 (1) | 0/40 (0) | 0.0490 / 0.0486 / 0.0466 | 0.00059 |
| GSE54456_normal | 40 | 1/40 (1) | 1/40 (1) | 1/40 (1) | 1/40 (1) | 0.0459 / 0.0453 / 0.0457 | 0.00070 |
| GSE121212_CTRL | 3 | **7/40** (82) | 3/40 (21) | 0/40 (0) | 6/40 (36) | 0.0423 / 0.0420 / 0.0311 | 0.00094 |
| GSE121212_CTRL | 5 | 2/40 (7) | 1/40 (3) | 1/40 (1) | 2/40 (36) | 0.0448 / 0.0448 / 0.0386 | 0.00069 |
| GSE121212_CTRL | 10 | 1/40 (2) | 0/40 (0) | 1/40 (1) | 0/40 (0) | 0.0538 / 0.0545 / 0.0516 | 0.00097 |
| GSE121212_CTRL | 18 | 2/40 (306) | 1/40 (84) | 1/40 (78) | 2/40 (239) | 0.0497 / 0.0500 / 0.0489 | 0.00104 |
| GSE121212_PSOnl | 3 | 4/40 (53) | 2/40 (51) | 0/40 (0) | 5/40 (40) | 0.0473 / 0.0489 / 0.0341 | 0.00090 |
| GSE121212_PSOnl | 5 | 1/40 (1) | 1/40 (1) | 0/40 (0) | 3/40 (3) | 0.0493 / 0.0498 / 0.0428 | 0.00110 |
| GSE121212_PSOnl | 10 | 3/40 (20) | 3/40 (23) | 3/40 (5) | 2/40 (32) | 0.0606 / 0.0640 / 0.0585 | 0.00160 |
| GSE121212_PSOnl | 13 | 0/40 (0) | 0/40 (0) | 0/40 (0) | 0/40 (0) | 0.0528 / 0.0547 / 0.0513 | 0.00083 |
| GSE186063_nonles | 3 | 5/40 (38) | 4/40 (27) | 2/40 (2) | 6/40 (43) | 0.0497 / 0.0509 / 0.0333 | 0.00115 |
| GSE186063_nonles | 5 | 3/40 (4) | 3/40 (4) | 2/40 (2) | 5/40 (16) | 0.0544 / 0.0541 / 0.0467 | 0.00099 |
| GSE186063_nonles | 10 | 0/40 (0) | 0/40 (0) | 0/40 (0) | 0/40 (0) | 0.0439 / 0.0445 / 0.0420 | 0.00053 |
| GSE186063_nonles | 14 | 0/40 (0) | 0/40 (0) | 0/40 (0) | 0/40 (0) | 0.0416 / 0.0417 / 0.0408 | 0.00034 |

Numbers in parentheses: total genes at FDR <= 0.05 over the 40 splits. Per-tertile fractions of p < 0.05 (mean expression tertiles, 40 seeds) lie within 3 SE of 0.05 in every cell for all three methods (largest: GSE54456 3v3 high tertile, moderated t 0.0840, SE 0.0141; Welch at 3v3 is conservative in every tertile, 0.030-0.050). The per-split KS distance from uniformity is not informative on its own, because the p-values of genes in one split are strongly dependent (mean KS D 0.039-0.147 across cells in the 112-split set).

At 3 per group, the moderated t and voom gave FDR discoveries in 10 of 40 GSE54456 splits (95% CI 0.127-0.412) and the moderated t in 7 of 40 GSE121212 control splits (0.073-0.328); limma in R gave the same count in every one of the 160 3v3 splits (Section 3), and `eBayes(robust = TRUE)` did not reduce it (12/40, 10/40, 5/40, 5/40 for GSE54456, GSE121212 CTRL, GSE121212 PSO_nl, GSE186063). The mean fraction p < 0.001 was 2.7 times nominal (0.0027) for GSE54456 at 3v3. From 5 per group upward every cell is within the binomial bound (moderated t, voom and Welch 0-3 of 40 splits; random-covariate model 0-5 of 40). Welch never exceeded the bound but has no power at 2-3 per group (Section 5).

### 4.2 Gene-set tests (112 null splits; Hallmark 50 sets, Reactome 1,670-1,671 testable sets)

| Test | Library | Splits with >= 1 set at FDR <= 0.05 | Total false sets (mean per split) | Mean fraction of sets p < 0.05 |
|---|---|---|---|---|
| GSEA, `auto` mode at 3v3 and 5v5 (automatic gene-permutation fallback) | Hallmark | **71/72** | 1617 (22.46) | 0.520 |
| GSEA, `auto` mode at >= 10 per group (sample permutation) | Hallmark | 3/40 | 6 (0.15) | 0.044 |
| GSEA, gene permutation selected (all sizes) | Hallmark | **111/112** | 2511 (22.42) | 0.519 |
| FRY (directional) | Hallmark | 7/112 | 47 (0.42) | 0.053 |
| CAMERA (Discovery, library source) | Hallmark | 1/112 | 2 (0.02) | 0.029 |
| ORA of a pasted top-200 list | Hallmark | **44/112** | 179 (1.60) | 0.086 |
| GSEA, `auto` mode at 3v3 and 5v5 (automatic gene-permutation fallback) | Reactome | 8/72 | 2921 (40.57) | 0.252 |
| GSEA, `auto` mode at >= 10 per group (sample permutation) | Reactome | 5/40 | 87 (2.17) | 0.051 |
| GSEA, gene permutation selected (all sizes) | Reactome | **16/112** | 6128 (54.71) | 0.255 |
| FRY (directional) | Reactome | 6/112 | 1137 (10.15) | 0.058 |
| CAMERA (Discovery, library source) | Reactome | 0/112 | 0 (0.00) | 0.024 |
| ORA of a pasted top-200 list | Reactome | **65/112** | 1206 (10.77) | 0.028 |

ORA of DE-hit lists could be run in only 12 of 112 null splits (the app requires at least 3 hits); in those, Hallmark gave false sets in 3 of 27 runs and Reactome in 13 of 27 (up to 149 sets in one GSE186063 3v3 split). FRY gave discoveries in 7 (Hallmark) and 6 (Reactome) of 112 splits, within the binomial limit of 11; the two largest counts (737 and 296 Reactome sets) occurred in the two splits in which the moderated t also called thousands of genes (GSE121212 CTRL 18v18 seed 0: 2,366 genes; GSE186063 3v3 seed 2: 4,308 genes), that is, splits whose groups differ systematically. FRY is a self-contained test and is expected to reject there; the remaining splits had 1-22 sets. GSEA preranked results (app) and fgsea agree (Section 3).

### 4.3 Discovery, co-expression and moderated F

| Pool | n | Discovery: splits with modules at FDR <= 0.05 | Nominal tier p <= 0.01: mean (range) v expected | Co-expression: permutations with >= 1 BH hit | Co-expression: mean frac p < 0.05 |
|---|---|---|---|---|---|
| GSE54456_normal | 3 | 0/10 | 1.4 (0-6) v 29.0 | 21/400 | 0.0498 |
| GSE54456_normal | 5 | 0/8 | 3.5 (0-12) v 29.0 | 17/320 | 0.0489 |
| GSE54456_normal | 10 | 0/5 | 11.8 (6-19) v 28.9 | 8/200 | 0.0568 |
| GSE54456_normal | 40 | 0/5 | 8.4 (1-19) v 29.1 | 7/200 | 0.0514 |
| GSE121212_CTRL | 3 | 0/10 | 0.9 (0-3) v 27.5 | 23/400 | 0.0493 |
| GSE121212_CTRL | 5 | 0/8 | 2.6 (0-5) v 27.5 | 14/320 | 0.0490 |
| GSE121212_CTRL | 10 | 0/5 | 7.4 (2-12) v 27.3 | 8/200 | 0.0450 |
| GSE121212_CTRL | 18 | 0/5 | 7.0 (0-13) v 27.4 | 7/200 | 0.0466 |
| GSE121212_PSOnl | 3 | 0/10 | 0.5 (0-2) v 27.3 | 14/400 | 0.0501 |
| GSE121212_PSOnl | 5 | 0/8 | 4.1 (0-11) v 27.4 | 18/320 | 0.0525 |
| GSE121212_PSOnl | 10 | 0/5 | 9.8 (5-15) v 27.2 | 9/200 | 0.0512 |
| GSE121212_PSOnl | 13 | 0/5 | 16.8 (10-27) v 27.2 | 13/200 | 0.0493 |
| GSE186063_nonles | 3 | 0/10 | 0.5 (0-4) v 28.8 | 27/400 | 0.0492 |
| GSE186063_nonles | 5 | 0/8 | 3.1 (0-6) v 28.7 | 14/320 | 0.0518 |
| GSE186063_nonles | 10 | 0/5 | 5.6 (3-8) v 28.5 | 8/200 | 0.0530 |
| GSE186063_nonles | 14 | 0/5 | 8.8 (4-17) v 28.7 | 5/200 | 0.0482 |

Pooled co-expression: 213 of 4480 permutations (0.0475; 95% CI 0.0415-0.0542) had at least one gene at BH FDR <= 0.05 against a permuted query gene; the mean fraction p < 0.05 was within 0.045-0.057 in every cell. Individual permutations of heavy-tailed query genes reached up to 0.369 of genes at p < 0.05 (the 90th percentile was 0.074-0.131), because the t-test on r assumes normality and one permuted outlier sample can align with sample-level outliers in many genes; the BH column absorbed this in about 95% of permutations. The Discovery nominal list was below the chance expectation in every cell and the app stated that it was not distinguishable from chance in all 112 sessions.

Three-group null (moderated F and Discovery random-set permutation path), 8 splits per cell:

| Pool | n per group | modF: splits with FDR hits (total genes) | modF mean frac p<0.05 | Discovery (families + HGNC): splits with FDR hits | Discovery mean frac p<0.05 | Discovery, Hallmark library: splits with FDR hits |
|---|---|---|---|---|---|---|
| GSE121212_CTRL | 3 | 1/8 (89) | 0.0734 | 0/8 | 0.0475 | 1/8 |
| GSE121212_CTRL | 5 | 2/8 (3) | 0.0500 | 0/8 | 0.0472 | 0/8 |
| GSE121212_CTRL | 10 | 0/8 (0) | 0.0319 | 0/8 | 0.0434 | 0/8 |
| GSE121212_CTRL | 12 | 0/8 (0) | 0.0535 | 0/8 | 0.0469 | 1/8 |
| GSE121212_PSOnl | 3 | 0/8 (0) | 0.0294 | 0/8 | 0.0438 | 0/8 |
| GSE121212_PSOnl | 5 | 0/8 (0) | 0.0533 | 0/8 | 0.0500 | 0/8 |
| GSE121212_PSOnl | 9 | 0/8 (0) | 0.0523 | 0/8 | 0.0447 | 0/8 |
| GSE186063_nonles | 3 | 0/8 (0) | 0.0578 | 0/8 | 0.0493 | 0/8 |
| GSE186063_nonles | 5 | 0/8 (0) | 0.0596 | 0/8 | 0.0571 | 1/8 |
| GSE186063_nonles | 9 | 0/8 (0) | 0.0429 | 0/8 | 0.0368 | 0/8 |
| GSE54456_normal | 3 | 1/8 (1) | 0.0669 | 0/8 | 0.0498 | 0/8 |
| GSE54456_normal | 5 | 1/8 (107) | 0.0506 | 0/8 | 0.0514 | 1/8 |
| GSE54456_normal | 10 | 0/8 (0) | 0.0408 | 0/8 | 0.0486 | 0/8 |
| GSE54456_normal | 20 | 0/8 (0) | 0.0446 | 0/8 | 0.0505 | 0/8 |

In total 5 of 112 three-group splits gave moderated-F discoveries (two of them with 89 and 107 genes) and 0 gave Discovery modules at FDR <= 0.05; both are within nominal limits.

### 4.4 Paired null (GSE121212 psoriasis patients, labels swapped at random within patient)

| Pairs | Splits | mod paired: splits with FDR hits | FRY paired Hallmark | FRY paired Reactome | GSEA paired Hallmark (mode) | GSEA paired Reactome | GSEA unpaired Hallmark (mode) | CAMERA unpaired Hallmark |
|---|---|---|---|---|---|---|---|---|
| 10 | 8 | 0/8 (0) | 0/8 (0) | 0/8 (0) | **8/8 (177)** (genes) | 4/8 (1599) | 0/8 (0) (samples) | 0/8 (0) |
| 26 | 8 | 0/8 (0) | 0/8 (0) | 0/8 (0) | **8/8 (176)** (genes) | 2/8 (961) | 1/8 (5) (samples) | 0/8 (0) |

With 3 or 5 pairs, a random within-patient swap puts every lesional sample in the same group with probability 1/4 or 1/16, and that labeling is the true lesional contrast. Of the 20 swaps at 3-5 pairs, 5 had a single orientation and all of them gave thousands of FDR hits; the 15 mixed-orientation swaps gave 0 splits with hits (paired moderated t). This design is therefore a valid null only for 10 or more pairs. Selecting the patient covariate makes GSEA fall back to gene permutation (`sample permutation is not available with covariates`), which gave false Hallmark sets in all 16 null splits at 10 and 26 pairs. The DE tab correctly warns that voom ignores the covariate ("Covariate adjustment (\"Patient\") applies only to the moderated t method; current voom ... results are UNADJUSTED").

### 4.5 Meta-analysis, consensus, RRHO and signature transfer

Meta-analysis on k independent null datasets (6 seeds per cell). Entries: runs with at least one gene at meta FDR <= 0.05 (total genes) / runs with at least one shared-signature gene (total). Failure threshold 2 of 6.

| k | n per group | REM (default; PI criterion on only when k >= 5) | FEM | Mean frac p<0.05 REM / FEM | Consensus genes (>= 2 datasets) |
|---|---|---|---|---|---|
| 3 | 3 | **6**/6 (81) / **4**/6 (17) | **6**/6 (709) / **6**/6 (102) | 0.0474 / 0.0875 | 0 |
| 3 | 5 | **5**/6 (86) / **3**/6 (12) | **5**/6 (461) / **4**/6 (58) | 0.0581 / 0.0929 | 0 |
| 3 | 10 | 1/6 (2) / 0/6 (0) | **5**/6 (20) / 0/6 (0) | 0.0413 / 0.0618 | 0 |
| 4 | 3 | **6**/6 (66) / **5**/6 (11) | **6**/6 (639) / **6**/6 (35) | 0.0495 / 0.0829 | 0 |
| 4 | 5 | **4**/6 (103) / 1/6 (5) | **6**/6 (858) / **3**/6 (8) | 0.0616 / 0.0938 | 0 |
| 4 | 10 | 0/6 (0) / 0/6 (0) | 1/6 (4) / 0/6 (0) | 0.0393 / 0.0609 | 0 |
| 5 | 3 | **3**/6 (6) / 2/6 (2) | **6**/6 (356) / **4**/6 (13) | 0.0411 / 0.0743 | 0 |
| 5 | 5 | 1/6 (1) / 0/6 (0) | **5**/6 (32) / 2/6 (3) | 0.0291 / 0.0539 | 0 |
| 5 | 10 | 1/6 (633) / 0/6 (0) | 2/6 (992) / 1/6 (1) | 0.0815 / 0.1086 | 0 |

Tail of the REM p-value distribution (all genes, prediction-interval criterion off, which does not affect p):

| Combination | Genes | frac p<0.05 | frac p<1e-3 | frac p<1e-4 | frac p<1e-5 | metafor genes at FDR<=0.05 (app) |
|---|---|---|---|---|---|---|
| 3 datasets of 3v3 (seed 0) | 17441 | 0.0673 | 0.00436 | 0.001548 | 0.000631 | 22 (22) |
| 3 datasets of 10v10 (seed 0) | 17299 | 0.0509 | 0.00306 | 0.000289 | 0.000116 | 2 (2) |
| 5 datasets of 3v3 (seed 0) | 17232 | 0.0168 | 0.00029 | 0.000058 | 0.000000 | 0 (0) |
| 5 datasets of 5v5 (seed 0) | 16760 | 0.0391 | 0.00179 | 0.000537 | 0.000060 | 0 (0) |

The REM z-test treats each dataset's standard error as known. With 3-5 samples per group the moderated-t standard errors rest on about 8 total degrees of freedom, so the pooled z has heavy tails: p < 1e-4 occurred 15 times more often than nominal for three 3v3 datasets, and every 3v3 combination with 3 or 4 datasets produced meta-FDR discoveries and, in 9 of 12 runs, shared-signature genes. With 10 per group the REM was within limits. FEM exceeded the limit in 7 of 9 cells, including three datasets of 10v10 (5 of 6 runs). metafor gives identical p-values. The consensus vote count was empty in all 54 runs.

RRHO on 54 pairs of independent null datasets: the maximum |signed -log10 p| on the 41 x 41 grid had median 54.9 and maximum 382.3, and exceeded a Bonferroni bound for the grid (4.53) in all 54 pairs (positive, concordant direction alone: 32 of 54). The correlation of t-statistics between two independent null datasets ranged from -0.43 to 0.33 (SD 0.174), against 1/sqrt(N) = 0.0075 expected if genes were independent. The hypergeometric overlap of FDR hit lists gave p < 0.05 in 0 of 54 pairs (most lists were empty).

Signature transfer, null signatures: the top three Discovery modules of each null session, scored in every null comparison dataset (463 transfers), gave empirical p <= 0.05 in 23 (0.0497; 95% CI 0.032-0.074); verdicts: does not replicate 440, partial replication 14, signature replicates 9. The real GSE54456 psoriasis signatures scored in 72 null splits (up and down, 144 transfers) gave empirical p <= 0.05 in 13 (0.0903; 95% CI 0.049-0.149; binomial 99th percentile 14), within the limit but above 0.05 in point estimate, consistent with the random-gene-set null understating the score variance of a co-regulated signature. All 8 positive transfers replicated (Section 5).

## 5. Positive controls

### 5.1 Classic psoriasis genes (20 genes; called = FDR <= 0.05, |log2FC| >= 1, expected direction)

| Contrast (test v reference) | n | Method | Genes called / passing filter | Not passing filter | Median rank by p | DE calls up / down | eBayes d0 |
|---|---|---|---|---|---|---|---|
| GSE54456 | 90 v 81 | mod | 20/20 | - | 29 | 1293 / 1513 | 4.85 |
| GSE54456 | 90 v 81 | voom | 20/20 | - | 56 | 1340 / 1653 | 4.65 |
| GSE54456 | 90 v 81 | welch | 20/20 | - | 76 | 1293 / 1513 | - |
| GSE121212_PSO_les_vs_nl | 27 v 27 | mod | 19/19 | IL17A | 1312 | 1106 / 1543 | 4.15 |
| GSE121212_PSO_les_vs_nl | 27 v 27 | voom | 19/19 | IL17A | 2325 | 1187 / 1711 | 4.07 |
| GSE121212_PSO_les_vs_nl | 27 v 27 | welch | 19/19 | IL17A | 1726 | 1106 / 1542 | - |
| GSE121212_PSO_les_vs_nl | 27 v 27 | mod+Patient | 19/19 | IL17A | 2579 | 1107 / 1529 | 3.88 |
| GSE121212_PSO_les_vs_CTRL | 27 v 37 | mod | 18/18 | IL17A,IL23A | 52.5 | 1411 / 2178 | 4.78 |
| GSE121212_PSO_les_vs_CTRL | 27 v 37 | voom | 18/18 | IL17A,IL23A | 347 | 1493 / 2459 | 4.57 |
| GSE121212_PSO_les_vs_CTRL | 27 v 37 | welch | 18/18 | IL17A,IL23A | 254.5 | 1410 / 2178 | - |
| GSE186063_les_vs_nl | 26 v 28 | mod | 19/19 | IL17A | 51 | 882 / 1159 | 4.18 |
| GSE186063_les_vs_nl | 26 v 28 | voom | 19/19 | IL17A | 122 | 949 / 1283 | 4.04 |
| GSE186063_les_vs_nl | 26 v 28 | welch | 19/19 | IL17A | 56 | 882 / 1159 | - |
| GSE186063_les_vs_nl | 26 v 28 | mod+Pair | 19/19 | IL17A | 98 | 874 / 1171 | 4.31 |
| GSE83645_pso_vs_uninv | 20 v 5 | mod | 20/20 | - | 79 | 932 / 919 | 2.93 |
| GSE83645_pso_vs_uninv | 20 v 5 | voom | 18/20 | - | 153 | 873 / 1041 | 2.90 |
| GSE83645_pso_vs_uninv | 20 v 5 | welch | 20/20 | - | 531.5 | 818 / 855 | - |
| GSE83645_pso_vs_uninv | 20 v 5 | mod+Patient | 20/20 | - | 146.5 | 1026 / 1120 | 3.02 |
| GSE41745_les_vs_nl | 2 v 3 | mod | 15/16 | CXCL8,IL17A,IL23A,CCL20 | 20.5 | 634 / 301 | 7.34 |
| GSE41745_les_vs_nl | 2 v 3 | voom | 15/16 | CXCL8,IL17A,IL23A,CCL20 | 31 | 627 / 182 | 7.52 |
| GSE41745_les_vs_nl | 2 v 3 | welch | 0/16 | CXCL8,IL17A,IL23A,CCL20 | 373.5 | 0 / 0 | - |
| GSE41745_les_vs_nl | 2 v 3 | mod+Patient | 15/16 | CXCL8,IL17A,IL23A,CCL20 | 13.5 | 867 / 690 | 16.38 |

All classic genes that pass filterByExpr were called by the moderated t and voom in every contrast with at least 5 samples per group (voom missed CXCL8 and IL23A in GSE83645, 18/20). IL17A (and IL23A in GSE121212 v controls; CXCL8, IL17A, IL23A, CCL20 in GSE41745) are removed by filterByExpr because of low counts. In GSE41745 (2 v 3) the moderated t called 15 of the 16 genes that pass the filter (KRT77 not called), and Welch called none: with 2 and 3 samples its Welch-Satterthwaite df is too small for any gene to reach FDR 0.05.

### 5.2 Expected pathways (sets at FDR <= 0.05 among 7 Hallmark / 8 Reactome expected up-regulated sets)

| Contrast | Library | ORA (DE up) | FRY | GSEA (mode) | CAMERA | All sets at FDR<=0.05: ORA / FRY / GSEA / CAMERA |
|---|---|---|---|---|---|---|
| GSE54456 | Hallmark | 7/7 | 7/7 | 0/7 (samples) | 1/7 | 20 / 42 / 0 / 1 |
| GSE54456 | Reactome | 8/8 | 8/8 | 0/8 (samples) | 0/8 | 91 / 1464 / 0 / 0 |
| GSE121212_PSO_les_vs_nl | Hallmark | 7/7 | 7/7 | 0/7 (samples) | 0/7 | 20 / 41 / 0 / 0 |
| GSE121212_PSO_les_vs_nl | Reactome | 8/8 | 8/8 | 0/8 (samples) | 0/8 | 133 / 1279 / 0 / 0 |
| GSE121212_PSO_les_vs_nl | paired, Hallmark | 7/7 | 7/7 | 7/7 (genes) | - | 20 / 43 / 33 / - |
| GSE121212_PSO_les_vs_nl | paired, Reactome | 8/8 | 8/8 | 8/8 (genes) | - | 119 / 1354 / 451 / - |
| GSE121212_PSO_les_vs_CTRL | Hallmark | 7/7 | 7/7 | 0/7 (samples) | 0/7 | 25 / 46 / 0 / 0 |
| GSE121212_PSO_les_vs_CTRL | Reactome | 8/8 | 8/8 | 0/8 (samples) | 0/8 | 209 / 1392 / 0 / 0 |
| GSE186063_les_vs_nl | Hallmark | 7/7 | 7/7 | 0/7 (samples) | 0/7 | 21 / 35 / 0 / 0 |
| GSE186063_les_vs_nl | Reactome | 8/8 | 8/8 | 0/8 (samples) | 0/8 | 94 / 1182 / 0 / 0 |
| GSE186063_les_vs_nl | paired, Hallmark | 7/7 | 7/7 | 7/7 (genes) | - | 20 / 38 / 31 / - |
| GSE186063_les_vs_nl | paired, Reactome | 8/8 | 8/8 | 8/8 (genes) | - | 96 / 1235 / 453 / - |
| GSE83645_pso_vs_uninv | Hallmark | 7/7 | 7/7 | 0/7 (samples) | 0/7 | 15 / 20 / 0 / 0 |
| GSE83645_pso_vs_uninv | Reactome | 8/8 | 6/8 | 0/8 (samples) | 0/8 | 58 / 691 / 0 / 0 |
| GSE83645_pso_vs_uninv | paired, Hallmark | 7/7 | 7/7 | 7/7 (genes) | - | 14 / 21 / 23 / - |
| GSE83645_pso_vs_uninv | paired, Reactome | 8/8 | 6/8 | 8/8 (genes) | - | 82 / 864 / 379 / - |
| GSE41745_les_vs_nl | Hallmark | 7/7 | 7/7 | 7/7 (genes) | 0/7 | 20 / 21 / 33 / 0 |
| GSE41745_les_vs_nl | Reactome | 8/8 | 7/8 | 8/8 (genes) | 0/8 | 143 / 354 / 514 / 0 |
| GSE41745_les_vs_nl | paired, Hallmark | 7/7 | 0/7 | 7/7 (genes) | - | 24 / 0 / 35 / - |
| GSE41745_les_vs_nl | paired, Reactome | 8/8 | 0/8 | 8/8 (genes) | - | 140 / 0 / 426 / - |

ORA of the up-regulated DE hits recovered all 7 Hallmark sets in every contrast, including GSE41745, and 8 of 8 Reactome sets in all contrasts. FRY recovered all expected sets except in the paired GSE41745 analysis (residual df 1: no set reaches FDR 0.05) and 6-7 of 8 Reactome sets in GSE83645/GSE41745; with large n FRY marks most of the library (for example 1,464 of 1,670 Reactome sets in GSE54456), as an earlier internal audit describes. Sample-permutation GSEA, the default when at least 7 samples are in each group, found no set at FDR <= 0.05 in any of the five contrasts in which it was used (Hallmark FDR <= 0.25: GSE54456 23, GSE121212_PSO_les_vs_nl 0, GSE121212_PSO_les_vs_CTRL 16, GSE186063_les_vs_nl 19, GSE83645_pso_vs_uninv 12); the 35 expected Hallmark set results in these five contrasts had nominal p 0.0019-0.23 and FDR 0.056-0.33. CAMERA found 1 Hallmark set at FDR <= 0.05 in GSE54456 and none elsewhere. Gene-permutation GSEA (automatic for GSE41745 and for every paired analysis) recovered the expected sets but also marked 23-35 of 50 Hallmark and 379-514 Reactome sets, comparable to the mean of 22 false Hallmark sets per null split in Section 4.2.

Discovery (families + HGNC groups, CAMERA): GSE54456: 0 of 2916 modules at FDR <= 0.05, 45 at p <= 0.01 (29 expected by chance); GSE121212_PSO_les_vs_nl: 0 of 2741 modules at FDR <= 0.05, 3 at p <= 0.01 (27 expected by chance); GSE121212_PSO_les_vs_CTRL: 0 of 2778 modules at FDR <= 0.05, 8 at p <= 0.01 (28 expected by chance); GSE186063_les_vs_nl: 0 of 2864 modules at FDR <= 0.05, 13 at p <= 0.01 (29 expected by chance); GSE83645_pso_vs_uninv: 0 of 3005 modules at FDR <= 0.05, 7 at p <= 0.01 (30 expected by chance); GSE41745_les_vs_nl: 0 of 2584 modules at FDR <= 0.05, 2 at p <= 0.01 (26 expected by chance). The SPRR, S100A, SERPINB and DEFB families were not among the top six modules except S100A* in GSE41745 (p 0.0083, FDR 0.99).

### 5.3 Low-n behaviour (subsamples, 5 seeds each)

| Subsample | Method | Classic genes called, mean (min) | Genes passing filter (mean) | DE calls up / down (mean) |
|---|---|---|---|---|
| GSE121212_pairsub_3 | mod | 6.2 (1) | 19.6 | 449 / 597 |
| GSE121212_pairsub_3 | mod+Patient | 7.0 (0) | 19.6 | 452 / 628 |
| GSE121212_pairsub_3 | voom | 5.4 (1) | 19.6 | 453 / 563 |
| GSE121212_pairsub_3 | welch | 1.4 (0) | 19.6 | 30 / 15 |
| GSE121212_pairsub_5 | mod | 17.8 (15) | 17.8 | 981 / 1138 |
| GSE121212_pairsub_5 | mod+Patient | 17.6 (15) | 17.8 | 1025 / 1265 |
| GSE121212_pairsub_5 | voom | 16.4 (13) | 17.8 | 1009 / 1190 |
| GSE121212_pairsub_5 | welch | 11.4 (0) | 17.8 | 683 / 660 |
| GSE186063_pairsub_3 | mod | 8.8 (0) | 18.6 | 221 / 202 |
| GSE186063_pairsub_3 | mod+Pair | 11.8 (0) | 18.6 | 403 / 374 |
| GSE186063_pairsub_3 | voom | 7.4 (0) | 18.6 | 170 / 169 |
| GSE186063_pairsub_3 | welch | 0.0 (0) | 18.6 | 0 / 0 |
| GSE186063_pairsub_5 | mod | 13.4 (5) | 18.0 | 666 / 657 |
| GSE186063_pairsub_5 | mod+Pair | 15.8 (9) | 18.0 | 899 / 959 |
| GSE186063_pairsub_5 | voom | 12.6 (3) | 18.0 | 647 / 577 |
| GSE186063_pairsub_5 | welch | 5.0 (0) | 18.0 | 181 / 185 |
| GSE54456_sub_10v10 | mod | 18.6 (18) | 18.8 | 1143 / 1397 |
| GSE54456_sub_10v10 | voom | 18.8 (18) | 18.8 | 1167 / 1442 |
| GSE54456_sub_10v10 | welch | 18.6 (18) | 18.8 | 1133 / 1352 |
| GSE54456_sub_2v2 | mod | 17.2 (15) | 19.4 | 558 / 356 |
| GSE54456_sub_2v2 | voom | 16.6 (11) | 19.4 | 536 / 317 |
| GSE54456_sub_2v2 | welch | 0.0 (0) | 19.4 | 0 / 0 |
| GSE54456_sub_2v3 | mod | 17.8 (15) | 20.0 | 873 / 639 |
| GSE54456_sub_2v3 | voom | 16.8 (15) | 20.0 | 795 / 738 |
| GSE54456_sub_2v3 | welch | 0.0 (0) | 20.0 | 0 / 0 |
| GSE54456_sub_3v3 | mod | 19.0 (18) | 19.8 | 921 / 957 |
| GSE54456_sub_3v3 | voom | 19.0 (18) | 19.8 | 904 / 941 |
| GSE54456_sub_3v3 | welch | 1.4 (0) | 19.8 | 35 / 20 |
| GSE54456_sub_5v5 | mod | 19.0 (18) | 19.0 | 1206 / 1041 |
| GSE54456_sub_5v5 | voom | 19.0 (18) | 19.0 | 1206 / 1029 |
| GSE54456_sub_5v5 | welch | 18.6 (17) | 19.0 | 827 / 469 |

| Subsample | Hallmark ORA | Hallmark FRY | Hallmark GSEA (mode) | Hallmark CAMERA | Reactome ORA | Reactome FRY |
|---|---|---|---|---|---|---|
| GSE121212_pairsub_3 | 3.8 | 1.2 | 7.0 (genes) | 0.0 | 3.2 | 1.6 |
| GSE121212_pairsub_5 | 7.0 | 7.0 | 7.0 (genes) | 0.0 | 7.6 | 7.2 |
| GSE186063_pairsub_3 | 2.8 | 1.4 | 7.0 (genes) | 0.0 | 4.0 | 1.6 |
| GSE186063_pairsub_5 | 5.0 | 4.4 | 6.8 (genes) | 0.8 | 6.4 | 5.0 |
| GSE54456_sub_10v10 | 7.0 | 7.0 | 0.0 (samples) | 0.4 | 8.0 | 7.0 |
| GSE54456_sub_2v2 | 4.2 | 1.4 | 7.0 (genes) | 0.0 | 5.8 | 0.4 |
| GSE54456_sub_2v3 | 6.8 | 2.6 | 7.0 (genes) | 0.0 | 8.0 | 3.0 |
| GSE54456_sub_3v3 | 7.0 | 5.0 | 7.0 (genes) | 0.0 | 8.0 | 4.4 |
| GSE54456_sub_5v5 | 7.0 | 7.0 | 7.0 (genes) | 0.0 | 8.0 | 7.2 |

Entries in the second table are the mean number of expected sets (of 7 Hallmark, 8 Reactome) at FDR <= 0.05, unpaired analysis. At 2-3 per group the moderated t and voom still recover most of the strongly induced genes in GSE54456 (17.2 of 19.4 at 2v2), but at 3 v 3 the same method gave FDR discoveries in 10-25% of null splits, depending on the pool (Section 4.1), so a list from a 3 v 3 experiment is not protected at the stated FDR. **Warnings at low n.** For GSE41745 (2 v 3) the app reported the eBayes prior (d0 7.34, total df 10.34, residual df 3) and, for GSEA, "only 10 distinct sample labelings (sample permutation needs 1000; at least 7 per group)" with the exploratory-ranking chip; Discovery stated that its nominal list was not distinguishable from chance. No message in the DE tab, the session notes or the methods text says that 2-3 samples per group give an inflated false-discovery risk, and Welch returned 0 calls without explanation. The only sample-size message in `analyzeGen` is for groups with fewer than 2 samples.

## 7. Feature status (calibration and power)

| Feature | Status | Basis |
|---|---|---|
| DE: moderated t (limma-trend) | validated_with_caveat | Null: within binomial limits at >= 5 per group in 4 pools (0-3/40 splits); 3v3: 10/40 and 7/40 splits with FDR hits; counts identical to limma on 272 splits; positive controls: all classic genes passing the filter called with >= 5 per group. |
| DE: voom | validated_with_caveat | As moderated t; identical to limma voom on 272 splits; 18-20/20 classic genes in positive controls. |
| DE: Welch t | validated_with_caveat | Never above nominal (0-3/40 splits); conservative at 3v3 (mean frac p<0.05 0.031-0.046); no power at 2-3 per group (0 calls in GSE41745). |
| DE: covariate-adjusted moderated t (random covariate) | validated_with_caveat | Continuous random covariate: within limits at >= 5 per group; 8/40 at GSE54456 3v3. Confounded binary covariate correctly refused. |
| DE: paired analysis (patient covariate) | validated | Paired null with 10 and 26 pairs: 0/16 splits with FDR hits; mixed-orientation 3-5 pair swaps 0/15. |
| Moderated F (3 groups) | validated | 5/112 three-group null splits with FDR hits (limit per cell respected); mean frac p<0.05 0.029-0.073. |
| ORA (over-representation) | defect_found | Recovers all expected Hallmark/Reactome sets in every positive control including GSE41745; ORA of correlated lists is anti-conservative; display exception with small queries. |
| FRY | validated | 7 (Hallmark) and 6 (Reactome) of 112 null splits with FDR sets, within limits; counts identical to limma::fry; recovers expected sets in all contrasts except paired GSE41745 (residual df 1). |
| CAMERA (via Discovery, library source) | validated_with_caveat | 1/112 and 0/112 null splits with FDR sets; identical to limma::camera; very low power. |
| GSEA, sample permutation | defect_found | Null: 3/40 (Hallmark) and 5/40 (Reactome) splits at >= 10 per group; no set at FDR <= 0.05 in five positive controls. |
| GSEA, gene permutation (automatic fallback for n <= 6 or covariates) | defect_found | 71/72 small-n null splits and 16/16 paired null splits with false Hallmark sets; same as fgsea. |
| Discovery screen, 2 groups (CAMERA) | defect_found | 0/112 null splits with FDR modules; nominal tier below chance and correctly labelled; 0 modules at FDR in six positive controls. |
| Discovery screen, > 2 groups (random-set permutation) | validated | 0/112 three-group null splits with FDR modules; mean frac p<0.05 0.037-0.057. |
| Co-expression p and BH FDR | validated_with_caveat | Permuted-query null: 213/4480 permutations (0.048) with a BH hit; mean frac p<0.05 0.045-0.057; individual heavy-tailed query genes up to 0.37 at p<0.05. |
| Meta-analysis (REM/FEM, shared signature) | defect_found | metafor parity 2.2e-13; null with 3-5 studies of 3-5 per group gives meta-FDR and signature genes; within limits with 10 per group (REM). |
| Consensus (vote count) | validated | 0 consensus genes in 54 null meta runs. |
| Signature transfer | defect_found | Null signatures: 23/463 (0.050) empirical p<=0.05; real signature into null splits 13/144; 8/8 positive transfers replicate; verdict ignores direction. |
| RRHO map | defect_found | Null pairs reach signed -log10 p up to 382; package parity not re-run here. |
| Compare: log2FC/t correlation and hit-overlap test | validated_with_caveat | Overlap test p<0.05 in 0/54 null pairs, but most null hit lists were empty, so the test was rarely exercised; t-statistic correlation between independent null datasets ranges -0.44 to 0.33 and log2FC correlation -0.40 to 0.27 (no inferential claim is made by the app). |

## 8. Not assessed

- Gene Explorer single-gene tests (pairwise Welch with Holm, ANOVA, Kruskal-Wallis, Spearman trend): not run under the null.
- Enrichment matrix (per-contrast FRY across contrasts), all-pairs DE and the Multiple-contrasts tab: not run (same moderated t/FRY code paths, but their multiplicity handling across contrasts was not tested).
- Gene-family scatter (runFamilies), hubs/network, heatmap cluster ORA, patterns/auto-k: not inferential in the sense tested here or not run.
- Power calculator: not assessed.
- Mouse data (GSE63310) and cross-species ortholog matching: not used in this calibration.
- GSE83645 was used only as a positive control; it has no healthy or uninvolved pool large enough for null splits (5 uninvolved samples).
- RRHO package parity on the null pairs (only the app map was examined; parity is reported in VALIDATION.md).
- Discovery and CAMERA with covariates or pairing: discoveryRun does not use covariates; not separately assessed or documented here.
- Behaviour in real browsers (all runs used JavaScriptCore with the kit DOM/Plotly stub).
- DESeq2 comparison: not run (not a method the app implements).

## 9. Deviations from the planned design

- Gene-set, Discovery and co-expression nulls used 10/8/5 seeds per cell (112 sessions), fewer than the 40 used for gene-level tests; their per-cell binomial limits are correspondingly wide. Co-expression used 4 random query genes per session with 10 permutations each (4,480 permutations).
- The GSE121212 paired null has 26 complete pairs: 26 PSO patients have both a lesional and a non-lesional sample in the matrix and design.
- Three-group nulls used the largest balanced size available (GSE54456 20, GSE121212 CTRL 12, PSO_nl 9, GSE186063 9 per group) in addition to 3, 5 and 10 where available.
- fgsea was run on Hallmark only (Reactome was too slow for a CI reference); preranked Reactome results are app-only.
- CAMERA was run through the Discovery screen, the only place the app exposes it, with the enrichment size limits (5-2000) instead of the Discovery defaults (3-100).

## 10. Files and CI

- `VAL_CAL_checks.csv`: all checks (check_id, name, dataset, reference, n_compared, max_deviation, tolerance, pass, kind, note).
- `calibration_suite.py` + `calibration_expected.json`: CI suite (78 metrics, fixed seeds; the gene-level null cells reuse exactly the splits of Section 4.1). Two runs on 0.22.0-beta: 64 s on 8 threads on an otherwise idle 14-core machine, and 3,330 s when the same machine had a 1-minute load average of about 96 from other jobs; both runs gave 67 PASS, 11 XFAIL, 0 FAIL and no value differing from the frozen values. Use a dedicated CI runner to stay within 15 minutes. Usage: `python3 calibration_suite.py --build benchside/index.html --kit <kit dir> --counts <dir with the five NCBI count matrices>`; exit status 1 if a metric outside its acceptance range is not a documented xfail. Each xfail carries its defect id.
- `cal_reference.R`: R reference for the null splits and the meta-analysis (edgeR/limma/fgsea/metafor); `robust_check.R`: eBayes(robust=TRUE) on the 3v3 splits; `cal_driver.js`: the JavaScript driver used for every app run (also embedded in the suite); `experiments.py`, `decode_gmt.py`: experiment definitions and GMT export of the built-in libraries.
