# Signature transfer with homology families and pooled replication (1.0.0-rc.3)

## What changed
Signature transfer (Compare datasets) now (1) carries signature genes without a one-to-one human–mouse ortholog through their homology family, each family scored once as the mean z of its members present in the scored dataset; (2) reports Hedges' g of the signature score in each comparison dataset; (3) pools g over the independent comparison datasets with the meta-analysis model selected in the meta-analysis card (default REML with the Hartung-Knapp test); (4) offers a signed up-and-down signature; (5) calls replication only when both a self-contained test (Welch t-test of the score between the groups) and a competitive test (against 1,000 random signatures) give p <= 0.05, in the expected direction. Random signatures are gene sets of the same size and the same numbers of up- and down-regulated genes drawn from the session's genes, transferred and scored exactly as the signature (the same random gene set in every dataset, so the pooled comparison keeps the between-dataset correlation of a fixed gene list).

## Numerical agreement
- Per-sample signature scores recomputed independently in Python (numpy) from the dataset's log-expression: maximum difference 3.2e-15 (engine check `signature_transfer`).
- Hedges' g and its variance against metafor 5.0.1 `escalc(measure = "SMD")`: relative difference 2.6e-13 and 9.5e-14 (GSE83645, up-and-down signature).
- Pooled estimate over four cohorts against `rma(method = "REML", test = "knha")`: estimate, SE, p, CI, tau^2, I^2 and Q agree to 5e-14 (relative).

## Calibration
Bound: 99th percentile of Binomial(R, 0.05) (18 of 200; 11 of 100). Counted calls: "replicates", "partial replication" or "opposite direction".

| Null | Replicates | Calls per comparison dataset | Pooled replication calls |
|---|---|---|---|
| Random mouse gene sets (30 up) from GSE143688, imiquimod vs control (covariates Line, Day), scored in the four human psoriasis cohorts | 200 | 9, 7, 9, 7 | 6 |
| Random mouse gene sets (75 up, 75 down), same setting | 200 | 10, 8, 12, 13 | 11 |
| Real GSE143688 up-regulated signature scored in three comparison datasets of 10 vs 10 random GSE54456 normal-skin samples | 100 | 5, 2, 2 | 2 |

Each test alone is not calibrated under one of the nulls: in the psoriasis cohorts, random gene sets differ between the groups by the self-contained test in 59–161 of 200 replicates (the contrast changes most genes); in healthy splits, the co-regulated real signature beats random signatures by the competitive test in 18–26 of 100 (its genes vary together between samples). Requiring both is within the bound in every setting.

## Positive control
GSE143688 imiquimod signature (FDR <= 0.05, |log2FC| >= 1, adjusted for line and day) scored in GSE54456, GSE121212, GSE186063 and GSE83645 with families:
- de_up: GSE54456 g = 4.03, replicates; GSE121212 g = 2.81, replicates; GSE186063 g = 2.3, replicates; GSE83645 g = 2.11, replicates; pooled g = 2.88, p = 0.00734, empirical p = 0.000999, replicates: True
- de_dn: GSE54456 g = -2.67, replicates; GSE121212 g = -2.63, replicates; GSE186063 g = -1.72, partial replication; GSE83645 g = -1.27, partial replication; pooled g = -2.16, p = 0.00734, empirical p = 0.000999, replicates: True
- de_both: GSE54456 g = 3.46, replicates; GSE121212 g = 2.96, replicates; GSE186063 g = 2.03, replicates; GSE83645 g = 1.69, differs, but no more than random genes; pooled g = 2.60, p = 0.00775, empirical p = 0.00599, replicates: True

Of the up-regulated mouse genes without a one-to-one ortholog, 265 belong to a homology family (among them Serpinb3a/b/c, Ccl7/Ccl12, Gsdma/Gsdma3, Ifi20x). Without families this list cannot be transferred at all (no gene matched); with families it replicates in GSE54456, GSE121212 and GSE186063 and the pooled effect over the four cohorts replicates (g = 2.99, 95% CI 1.25 to 4.72).

## Power with two datasets
With two comparison datasets the Hartung-Knapp pooled test has 1 degree of freedom and rarely reaches p <= 0.05; the application states this. Per-dataset results remain valid.
