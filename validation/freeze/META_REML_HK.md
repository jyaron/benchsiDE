# Meta-analysis: REML and Hartung-Knapp (design-freeze change)

## Parity with metafor 5.0.1

Inputs: the per-gene log2FC and variance of the four psoriasis cohorts (GSE54456 session; GSE121212, GSE186063 and GSE83645 as comparison datasets, paired cohorts adjusted for patient), 19,518 genes measured in at least two cohorts. Each gene was refitted in R with `rma(yi, vi, method, test, control = list(threshold = 1e-12))`.

| Model / test | max abs. diff. pooled log2FC | max abs. diff. SE | max rel. diff. p | max abs. diff. tau² | genes at FDR ≤ 0.05, app / R |
|---|---|---|---|---|---|
| REML, Hartung-Knapp | 4.9e-10 | 4.0e-10 | 2.6e-9 | 1.1e-9 | 7,010 / 7,010 |
| REML, z | 4.9e-10 | 5.6e-10 | 9.5e-8 | 1.1e-9 | 12,909 / 12,909 |
| DerSimonian-Laird, Hartung-Knapp | 5.3e-15 | 4.4e-16 | 7.4e-13 | 3.6e-15 | 7,007 / 7,007 |

Confidence intervals agree to 1.8e-9 and prediction intervals to 6.1e-9 (REML) and 1.9e-12 (DL). The REML differences reflect the convergence tolerance of the two iterative solutions; for the gene with the largest tau² difference the restricted log-likelihood is identical to 12 digits. During development, Fisher scoring converged slowly on flat likelihoods and, in five genes, stopped at an interior local maximum below the likelihood at tau² = 0; the final implementation uses Newton steps with the observed information, step halving and a check of the boundary.

## Null calibration

Normal-skin samples of GSE54456 (81) were split at random into K datasets of n vs n (disjoint samples), each analysed with the default moderated t, and combined by every model. A replicate has discoveries if at least one gene reaches meta FDR ≤ 0.05 (no fold-change threshold; genes measured in all K datasets). The bound is the 99th percentile of Binomial(R, 0.05). Source: `meta_null_variants.csv`.

Replicates with discoveries:

|   K |   n |   R |   bound |   FEM_false |   REML_adhoc |   REML_false |   REML_hk |   REM_adhoc |   REM_false |   REM_hk |
|----:|----:|----:|--------:|------------:|-------------:|-------------:|----------:|------------:|------------:|---------:|
|   2 |   3 |  60 |       7 |          10 |            0 |            7 |         1 |           0 |           7 |        1 |
|   3 |   3 | 100 |      11 |          23 |            0 |           15 |         3 |           0 |          14 |        3 |
|   3 |   5 |  60 |       7 |           6 |            0 |            4 |         2 |           0 |           4 |        2 |
|   4 |   3 |  60 |       7 |          15 |            0 |            7 |         3 |           0 |           8 |        3 |
|   4 |   5 |  60 |       7 |           5 |            0 |            4 |         1 |           0 |           5 |        1 |
|   5 |   3 |  60 |       7 |          17 |            0 |            8 |         2 |           0 |           8 |        2 |

Fraction of genes with p < 0.05 (mean over replicates; nominal 0.05):

|   K |   n |   FEM_false |   REML_adhoc |   REML_false |   REML_hk |   REM_adhoc |   REM_false |   REM_hk |
|----:|----:|------------:|-------------:|-------------:|----------:|------------:|------------:|---------:|
|   2 |   3 |       0.053 |        0     |        0.041 |     0.052 |       0     |       0.041 |    0.052 |
|   3 |   3 |       0.054 |        0     |        0.04  |     0.054 |       0     |       0.04  |    0.054 |
|   3 |   5 |       0.064 |        0     |        0.049 |     0.053 |       0     |       0.049 |    0.053 |
|   4 |   3 |       0.054 |        0.002 |        0.04  |     0.048 |       0.002 |       0.039 |    0.048 |
|   4 |   5 |       0.046 |        0.001 |        0.033 |     0.045 |       0.001 |       0.033 |    0.045 |
|   5 |   3 |       0.06  |        0.007 |        0.044 |     0.058 |       0.007 |       0.044 |    0.058 |

Variants: REM = DerSimonian-Laird, REML = restricted maximum likelihood, FEM = fixed effect; `false` = Wald z, `hk` = Hartung-Knapp-Sidik-Jonkman, `adhoc` = truncated Hartung-Knapp.

With the z test, DerSimonian-Laird exceeds the bound with three, four and five datasets of 3 vs 3, and REML with three and five; the fixed-effect model exceeds it in every configuration with 3 samples per group. With the Hartung-Knapp test both estimators stay within the bound everywhere, and the fraction of p < 0.05 is close to nominal. The truncated form is very conservative. REML with Hartung-Knapp is the default from this version.
