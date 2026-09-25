# Design types

The **Design type** setting on the review step determines how group labels are interpreted.

| Setting | Use for | Effect |
|---|---|---|
| Auto-detect | Most cases | Numeric if every group label is a number; categorical otherwise |
| Time course / dose–response (ordered numeric) | Time points, doses, concentrations | Numeric x-axis spaced by value; group means connected by lines; Spearman trend statistics |
| Categorical (genotypes / treatments / conditions) | Genotypes, treatments, cell types | Categorical axis; no connecting lines; no trend statistics |

## Why the distinction matters

A trend statistic assumes that groups have a meaningful order and spacing. Knockout lines or treatments have neither, so fitting a trend across them produces a number without meaning. Conversely, treating time points as unordered categories discards information.

If numeric labels are codes rather than quantities (for example genotype IDs 1, 2, 3), choose **Categorical** explicitly. The review step previews the consequences of the current choice before analysis.

## Group order

For categorical designs, groups can be reordered with the arrows in the sample bar. The first group is the default baseline for fold changes and pairwise comparisons.
