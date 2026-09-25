# Normalization

## TMM (default for counts)

The trimmed mean of M-values method of edgeR (`calcNormFactors`, method `TMM`, default settings: 30% trimming of log-ratios, 5% trimming of mean expression, precision weighting). One sample is chosen as the reference (the sample whose upper-quartile-scaled library size is closest to the mean). A scaling factor is computed for every other sample, and factors are rescaled to a geometric mean of one. The effective library size of each sample is its column total multiplied by its factor.

Normalized expression for display and for the moderated t-test is log₂(CPM + 1), with CPM computed from effective library sizes.

TMM factors match edgeR to machine precision on the validation datasets. On datasets with fractional counts, the column totals are accumulated in plain sequential double precision, which reproduces R's result; compensated summation would be more accurate in isolation but breaks parity, because TMM's trimming depends on exact ranks.

## Library-size CPM only

Scales by column totals without TMM factors. Available for comparison, or when the TMM assumption (most genes not differentially expressed) is known to fail.

## Normalized input

CPM/TPM/FPKM input is log₂-transformed as log₂(value + 1). log₂ input is used as supplied.

## Batch-adjusted display

The **display batch-adjusted** option (PCA card) removes a chosen factor's effect for visualization, as in `limma::removeBatchEffect`. It changes plots only; no test uses adjusted values. For statistical adjustment use covariates (see [Covariates and paired designs](Covariates-and-Paired-Designs.md)).
