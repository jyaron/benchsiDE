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

Supply the subject (or block) identifier as a categorical covariate. This fits the subject as a fixed-effect block, which is the standard analysis for paired samples. If the design file has no subject column, add one in the review step: open **Edit the design**, type a column name (for example `Patient`), click **Add column**, and enter each sample's value in the table. Added columns are saved with the session.

## Design check and confounding

The **Design check** (review step, and under the sample list after analysis) tabulates each design factor against the groups. For each pair of groups it reports whether the factor can be adjusted for (the factor varies within the comparison) or is fully confounded with it (the groups do not share any level of the factor). A confounded factor cannot be adjusted for, because its effect cannot be separated from the group effect.

Example, GSE186063: all normal-skin samples come from ankylosing spondylitis patients, while lesional and non-lesional samples come from psoriasis and psoriatic arthritis patients. Diagnosis can therefore be adjusted for in lesion vs non-lesion, but not in any comparison with normal skin.

Options for a confounded comparison:

| Option | Where |
|---|---|
| Restrict the samples to one level of the confounding factor (for example psoriasis patients only) | "Include samples with …" checkboxes under the sample list |
| Group by a combined factor (for example Type × Diag) and compare within one level | "Group by factor" in the review step |
| Compare only groups in which the factor varies, with the factor as a covariate | Differential Expression tab, "adjust for" |
| Report the comparison as confounded | Methods and figure legends |

## Refused designs

benchsiDE refuses to fit, and states the reason, when:

- a covariate is completely confounded with the groups being compared,
- a covariate is constant within the contrast, or
- the adjusted model has no residual degrees of freedom.

## Other methods

voom and Welch t results are unadjusted. When a covariate is selected with these methods, the results carry an explicit warning.
