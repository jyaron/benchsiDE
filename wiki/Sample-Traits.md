# Sample traits

The **Sample traits** tab tests which genes, gene sets and modules are associated with a numeric value measured on each sample: a clinical score (PASI), immunohistochemistry counts of the same section, epidermal thickness, a serum level. Unlike a covariate, which is adjusted for, a trait is the variable of interest.

## Setting up

- Add the trait as a column of the design file, or with **Edit the design → Add column** in the review step. Any numeric column with three or more distinct values is offered. Blank or NA leaves only that sample out.
- **Transform:** log₂(x + 1) for counts, logit for proportions or percentages, rank for skewed values or outliers.
- **Samples / adjust for group:** when the trait differs between groups (PASI is 0 in controls), either keep one group or select *adjust for group*. Otherwise the result largely repeats the group comparison, and the tab says so.
- **Repeated samples:** select the column identifying the subject or animal when it contributes several samples. If the trait is constant within each subject (PASI, age), each subject is summarized by the mean of its samples (block means). If it varies between a subject's samples (IHC of each section), a mixed model with limma's duplicateCorrelation is used. *auto* chooses between the two.

## Results

| Output | Method |
|---|---|
| Gene table, volcano, gene plot | limma linear model on the trait, eBayes with intensity trend; effect per unit and per SD of the trait; Cook's distance flags genes where one sample (or subject) determines the slope |
| Gene sets | CAMERA on the trait coefficient (Wu & Smyth 2012) over the loaded library |
| Module × trait matrix | module eigengene (WGCNA) of each hub neighbourhood, gene set or Discovery module tested against every trait in the same model; partial r, BH over the matrix |
| Trait × trait | Spearman correlation between traits (per subject when blocks are set) |

## Why block means for subject-level traits

In calibration on GSE83645 (4 lesional biopsies from each of 5 patients), a random patient-level trait gave 25% of genes at p < 0.05 when the biopsies were treated as independent, and 8.4% with the mixed model, which uses one consensus correlation for all genes; with block means it was 6.3%. For a trait that varies within subjects, the mixed model was calibrated (5.0%). Full results: VALIDATION.md, section "Sample traits". The limitation of a consensus correlation is discussed by Hoffman & Roussos (Bioinformatics 37:192, 2021).

## Limits

- Associations are correlational. With traits that differ between groups, adjust for group or restrict to one group.
- Block means require the trait to be constant within blocks; covariates that vary within blocks cannot enter that model and are dropped (the tab lists them).
- The module × trait matrix is limited to 300 modules; use the CAMERA table for large libraries.
- Discovery modules were selected for differing between groups; their association with a trait that also differs between groups is partly circular.
