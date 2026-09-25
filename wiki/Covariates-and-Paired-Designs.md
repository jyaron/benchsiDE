# Covariates and paired designs

## Covariate adjustment

With a design file containing more than one factor, the moderated t method can adjust for additional factors. The model per gene is

log₂ expression ~ group + covariate₁ + covariate₂ + …

and the test is on the group coefficient. Variance moderation uses the residual degrees of freedom of the adjusted model.

| Covariate type | Encoding |
|---|---|
| Categorical | Indicator columns |
| Numeric with more than two distinct values | Continuous slope |
| Numeric with two values | Indicator |

Covariate-adjusted results match limma `lmFit` with the corresponding design matrix and `eBayes`.

## Paired and blocked designs

Supply the subject (or block) identifier as a categorical covariate. This fits the subject as a fixed-effect block, which is the standard analysis for paired samples.

## Refused designs

benchsiDE refuses to fit, and states the reason, when:

- a covariate is completely confounded with the groups being compared,
- a covariate is constant within the contrast, or
- the adjusted model has no residual degrees of freedom.

## Other methods

voom and Welch t results are unadjusted. When a covariate is selected with these methods, the results carry an explicit warning.
