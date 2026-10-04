# benchsiDE validation report (1.0.0-rc.1)

This report summarizes the end-to-end validation of benchsiDE on public data. Detailed results for each part are in `validation/study/`: VAL_DE.md (differential expression against limma, edgeR and DESeq2), VAL_ENR.md (gene-set tests, Discovery, clustering), VAL_XDS.md (cross-dataset comparison and meta-analysis), VAL_CAL.md (false-positive calibration on null data and positive controls) and VAL_UNT.md (Gene Explorer, Patterns, QC, power, input handling, sessions, exports). The study tested 0.22.0-beta.

## Release candidate 1.0.0-rc.1

The frozen 1.0.0-rc.1 build was checked again in full: R parity (1,096 checks, `validation/freeze/PARITY_FREEZE.md`), null calibration and binomial-thinning simulation (499 checks, including a non-skin null pool, GSE57945; `CALIBRATION_FREEZE.md`) and regression against the 0.23.0 baseline (227 checks; `REGRESSION_FREEZE.md`). The results are summarized in VALIDATION.md ("Release candidate 1.0.0-rc.1"). Interaction contrasts were checked against limma on GSE143688 (`INTERACTION.md`), REML and Hartung-Knapp meta-analysis against metafor (`META_REML_HK.md`), and every literal-valued check in the earlier study was traced (`CHECK_PROVENANCE.md`). The sections below describe the earlier study.

## Data

GSE54456 (90 psoriasis, 81 normal skin), GSE121212 (psoriasis and atopic dermatitis, lesional and non-lesional, paired by patient), GSE186063 (lesional and non-lesional; pairs inferred from GEO metadata and labelled as inferred), GSE83645 (5 patients, uninvolved skin and lesional sites), GSE41745 (2 lesional, 3 non-lesional after one sample missing from the count matrix), and GSE63310 (mouse mammary cell types; the external-validation dataset).

## Checks

| Part | Checks |
|---|---|
| DE vs R | 1481 |
| Cross-dataset | 343 |
| False-positive calibration | 787 |
| Gene-set tests | 652 |
| Previously untested analyses | 140 |

The external validation (48 checks against R on GSE63310; 59 since the hub-neighbourhood, hub-across-datasets and sample-trait checks were added) and the in-app self-test (19 components; 21 from 1.0.0-rc.1) pass in full.

## Feature status

| Analysis | Status | Basis |
|---|---|---|
| Filtering (filterByExpr, CPM) | Validated | edgeR parity; external validation |
| TMM normalization, edgeR log-CPM | Validated | edgeR parity |
| Moderated t (limma-trend), full design or two groups | Validated | limma parity on GSE63310 and four psoriasis datasets; FDR approximate at ≤ 3 per group (warning shown) |
| voom | Validated | limma parity |
| Covariates and blocking factors | Validated | limma parity including paired designs |
| Welch t (per gene) | Validated | scipy and R t.test parity; constant genes untestable as in R (from 1.0.0-rc.1) |
| Moderated F (Patterns) | Validated | limma parity |
| ORA | Validated | scipy parity; anti-conservative for co-regulated lists (stated) |
| FRY | Validated | limma parity, including covariates and full design; with very many covariate columns p depends on the design parametrisation, as in limma |
| GSEA, sample permutation | Validated | calibrated on null splits; conservative; reported at FDR ≤ 0.25 |
| GSEA, gene permutation | Exploratory | not calibrated; ranking only, no FDR reported |
| CAMERA / Discovery screen | Validated, low power | limma parity; calibrated; no module at FDR ≤ 0.05 on six psoriasis contrasts |
| Discovery, more than two groups | Exploratory | permutation test with a p-value floor |
| Meta-analysis, REML + Hartung-Knapp (default from 1.0.0-rc.1) | Validated, low power with 2–4 datasets | metafor parity; within the null bound in all 12 configurations (2–5 datasets, 3–10 per group) and FDR controlled under heterogeneous effects; power 0–0.1% with 2–3 datasets (stated in the app) |
| Meta-analysis, DerSimonian-Laird or REML with z test; fixed effect | Validated, exploratory | metafor parity; above the null bound with 4–5 datasets of 3 per group (fixed effect at 3 per group throughout) and FDR 8–30% when true effects differ between datasets |
| Consensus vote count | Descriptive | no pooled test |
| t and log₂FC correlation between datasets | Validated | numpy parity |
| RRHO map | Descriptive | RRHO package parity; p-values assume independent genes |
| Signature transfer | Validated | positive control separates groups; direction required |
| Mouse–human orthologs | Validated | MGI one-to-one pairs |
| Interaction contrasts (from 1.0.0-rc.1) | Validated | limma parity on GSE143688, 12 fits (t 4.6e-13; counts exact) |
| Hub module preservation | Validated | WGCNA parity in 15 hub × dataset combinations |
| Heatmap clustering, k-means, silhouette auto-k | Descriptive | row-order invariant; silhouette as cluster::silhouette |
| Gene Explorer tests (ANOVA, Kruskal-Wallis, Welch/Holm) | Validated | scipy/statsmodels parity |
| Co-expression r, slope, p, FDR | Validated | R parity |
| Power calculator | Validated | RNASeqPower parity; pilot CV now biological CV |
| QC (MDS, RLE, densities, saturation) | Validated / descriptive | limma plotMDS parity; densities vs R KDE |
| Excel/Prism export | Checked | 185 plot exports in four configurations (1.0.0-rc.1 regression) |
| Browser behaviour | Checked in part | self-test 19/19 in Chrome 154 (macOS); other browsers not assessed |

*Validated*: agrees numerically with the reference implementation and behaves as expected on null and positive-control data. *Descriptive*: the output is a display or summary without a calibrated test. *Exploratory*: a test whose error rate is not controlled; results are labelled in the app.

## Real-browser check

The in-app self-test was run by the author in Chrome 154 on macOS with the GSE54456 matrix loaded: all 19 components reproduced their references (111 ms). The automated browser runs could not be performed in the validation environment.

## Not assessed

Behaviour in other browsers (Firefox, Safari, Edge) and the browser run of the external validation and of the full analysis workflow; DESeq2 and edgeR quasi-likelihood (used for comparison only, not implemented); Discovery on designs with more than two groups beyond the permutation floor; every plot's Excel export. Each part's own limits are listed in its report.

## Reproducing

The R reference script for the external validation is `validation/reference.R`; `validation/compare.js` runs the comparison in the page and is used by `tests/validation.spec.js` in CI. The study's scripts are referenced from each report in `validation/study/`.
