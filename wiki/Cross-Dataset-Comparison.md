# Cross-dataset comparison

## Principle

Each dataset is filtered, normalized and tested within itself. Raw or normalized values are never merged, pooled or renormalized across datasets, and no p-value is computed on expression pooled across datasets. Comparisons operate on within-dataset statistics of matched genes. Merging studies would require assumptions about batch structure that the tool does not make.

## Adding datasets

Any number of datasets can be added to the session, within memory. Each dataset uses about 8 bytes per gene per sample (roughly 30 MB for 20,000 genes × 170 samples). The slot list shows the total, warns above 1.5 GB, and refuses a dataset that would take the total above 3 GB. Each has its own optional design file and species setting.

## Gene matching

| Situation | Matching |
|---|---|
| Same species | Case-normalized gene symbol |
| Mouse and human | MGI strict one-to-one orthologs |
| Custom ID map loaded | The map, applied before the rules above |

Identifiers are decoded with each dataset's own species annotation. If matching yields few genes, the slot reports the identifier types detected on each side and the remedy.

## Pairwise comparison

| Output | Statistic |
|---|---|
| Fold-change scatter | log₂FC in X against log₂FC in Y, coloured by quadrant; Pearson r |
| t concordance | Pearson correlation of moderated t-statistics |
| Hit-list overlap | One-sided hypergeometric (Fisher) test of shared up and shared down hits |
| RRHO map | Rank–rank hypergeometric overlap, computed as in the Bioconductor RRHO package; under-enrichment shown as a signed extension |
| Shared and discordant lists | Genes significant in both datasets, ranked |

## All datasets

**Compute across all datasets** produces the correlation matrix and a consensus table. The **forest plot** shows one gene's log₂FC with confidence intervals in every dataset; its axis always includes 0 and is symmetric about it.

## Gene Explorer in comparison datasets

When comparison datasets are loaded, Gene Explorer shows the selected gene in each of them, one panel per dataset, using that dataset's groups, normalization and statistics (ANOVA, Welch tests with Holm adjustment, and the dataset's selected contrast). Genes are matched as in the table above. A symbol can be typed to show a gene that is not in the session.

## Dataset names

Each dataset has an editable name (GEO files default to their accession, for example `GSE54456`). The session's dataset can be renamed too. Names are used in every plot, table, export and the methods text, and outputs already on screen are redrawn with the new name.

## Meta-analysis

The meta-analysis card pools each dataset's own gene-level estimate: the log₂ fold change and its standard error from that dataset's test. Expression values are never combined.

| Quantity | Definition |
|---|---|
| Fixed effect | inverse-variance weighted mean of the log₂FCs |
| Random effects (default) | weights 1/(SE² + τ²); τ² by restricted maximum likelihood (REML, default) or by the DerSimonian–Laird method of moments |
| Heterogeneity | Cochran's Q with its χ² p-value; I² = (Q − df)/Q, truncated at 0 |
| Test | Hartung–Knapp–Sidik–Jonkman (default): the SE of the pooled log₂FC is rescaled by the weighted spread of the datasets around it and referred to a t distribution on k − 1 df; the truncated form keeps that SE at or above the model-based SE; or the Wald z test. Benjamini–Hochberg across genes |
| Prediction interval | 95% interval for the effect in a new study, μ ± t(0.975, k − 2)·√(τ² + SE²) (Riley et al. 2011; metafor `predict(..., predtype = "Riley")`); needs 3 or more datasets |
| Shared signature | meta FDR, minimum \|pooled log₂FC\|, the same direction in every dataset; optionally a prediction interval excluding 0 and a maximum I² (both off by default) |
| Robustness | the signature recomputed with each dataset left out (three or more datasets) |

The gene universe is the session's genes matched in at least the chosen number of datasets (default: a majority). Each dataset's variance is taken as (log₂FC / z)², where z is the normal quantile with the same two-sided p as the dataset's moderated t on its own degrees of freedom; treating small-study standard errors as known made the pooled test anti-conservative on null data. Pairwise comparisons and the consensus table are likewise restricted to genes present in the session dataset; the consensus counts significance among the datasets in which a gene is measured.  With the df adjustment switched off, results agree with `metafor::rma` (VALIDATION.md). With few datasets, τ² and I² are imprecise; the random-effects model is the conservative choice when datasets differ in platform, tissue or species.

## Signature transfer

A module from dataset A (up- or down-regulated DE genes, or a Discovery module) is scored in each other dataset as a mean z-score. Discrimination is summarized as the Mann–Whitney AUC between that dataset's groups and calibrated against 200 seeded size-matched random signatures scored identically.

## Choice of model and test

The default is REML random effects with the Hartung–Knapp–Sidik–Jonkman test. The choice rests on two kinds of evidence from the release-candidate calibration (`validation/freeze/CALIBRATION_FREEZE.md`):

- **Null data** (2–5 datasets built from different healthy or non-lesional skin pools, 3, 5 or 10 samples per group, 120 replicates each). REML or DerSimonian–Laird with the z test exceeded the 99% binomial bound with 4–5 datasets of 3 per group (13–22 of 120 replicates with a discovery, bound 12); the fixed-effect model exceeded it in every configuration with 3 per group. REML with Hartung–Knapp stayed within the bound in all 12 configurations.
- **Simulated true effects** (binomial thinning of real null counts). When the true effect is the same in every dataset, all models control the FDR. When it differs between datasets (τ = 0.4), the observed FDR was 8–22% with REML or DerSimonian–Laird and the z test and 18–30% with the fixed-effect model, whereas with Hartung–Knapp it was not significantly above 5%; in that setting Hartung–Knapp detected at most 3% of true effects.

The cost is power. Hartung–Knapp refers the pooled estimate to a t distribution on k − 1 degrees of freedom, and the FDR is controlled over thousands of genes, so in the simulations it detected essentially no true effect with 2–3 datasets (at most 0.1%), 5–37% with 4 and 26–57% with 5. REML with the z test detected 29–66% with 2–3 datasets and 48–75% with 4–5. The application states this in the meta-analysis summary when 2–4 datasets are combined.

| Datasets | Recommended reading |
|---|---|
| 2–3 | The Hartung–Knapp signature will usually be empty. Read the pooled estimates, the per-dataset results and the leave-one-out direction column; for an exploratory signature choose the z test and say so. |
| 4 | Hartung–Knapp signatures are conservative; the z test finds more but does not control the FDR when effects differ between datasets. |
| 5 or more | Hartung–Knapp is the recommended confirmatory analysis. |

The truncated Hartung–Knapp form is very conservative (no discovery in 1,440 null replicates and no power except with 5 datasets of 10 per group) and is intended as a sensitivity check.

The leave-one-dataset-out table reports both the signature recomputed without each dataset and the number of signature genes whose pooled estimate keeps its direction at nominal p ≤ 0.05. With Hartung–Knapp the first count falls sharply when 4 datasets become 3, because of the loss of a degree of freedom; the second does not depend on the correction across genes (on the four psoriasis cohorts it retained 84–94% of the 1,773 signature genes).

Sessions saved before version 1.0 restore the model they were run with (DerSimonian–Laird, z test).

## References

- DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials 7:177–188, 1986.
- Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med 20:3875–3889, 2001.
- Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Stat Med 21:1539–1558, 2002.
- IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol 14:25, 2014.
- Knapp G, Hartung J. Improved tests for a random effects meta-regression with a single covariate. Stat Med 22:2693–2710, 2003.
- Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. BMJ 342:d549, 2011.
- Sidik K, Jonkman JN. A simple confidence interval for meta-analysis. Stat Med 21:3153–3159, 2002.
- Viechtbauer W. Bias and efficiency of meta-analytic variance estimators in the random-effects model. J Educ Behav Stat 30:261–293, 2005.
- Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw 36(3), 2010.

Plaisier SB et al. (2010) Nucleic Acids Res 38:e169. Wu D and Smyth GK (2012) Nucleic Acids Res 40:e133.

## Hub across datasets

Open it with **Across datasets** in the hub neighbourhood panel (Co-expression tab). The neighbourhood is defined in the session (the discovery dataset) and tested in every comparison dataset, on the samples of that dataset's chosen contrast. Datasets flagged as duplicates of the session are excluded.

| Question | Method | How to read it |
|---|---|---|
| Is the neighbourhood preserved? | WGCNA module preservation (Langfelder et al. 2011): Zsummary = mean of Zdensity and Zconnectivity, each against 100 random gene sets of the same size in the test dataset | Use the verdict, not the published Z thresholds. **Beyond DE-matched** means the members stay co-expressed more than random sets of genes with the same differential expression (p ≤ 0.05); this is the evidence that the module, not only the disease difference, replicates. The percentile compares with neighbourhoods of random genes built the same way. |
| Is it a hub in every dataset? | Robust rank aggregation (Kolde et al. 2012) of connectivity ranks, over one common gene universe | Consensus hubs at FDR ≤ 0.05; the hub's own ranks are shown. With correlation across all samples, consistently DE genes tend to be hubs everywhere; use within-group mode for hubs of co-expression. |
| Does it rewire? | Fisher z test of the hub–partner correlation between the contrast groups, per dataset | Use the **DE-matched** FDR: correlations among DE genes can change between conditions generically. **Replicated in** counts independent comparison datasets; the discovery dataset and duplicates are not counted, and datasets are not pooled. Needs 6 or more samples per group. |
| How large are the members' effects? | The meta-analysis (run it first), pooled with and without the discovery dataset | The discovery dataset's effects can be inflated because the neighbourhood was selected there; the open diamonds exclude it. |

Calibration on public psoriasis data: in VALIDATION.md (section "Hub across datasets"). In that calibration, random gene sets exceeded the published Z = 2 threshold in 100% of cases, random sets of up-regulated genes were "preserved" against random genes in 87% of cases, and random pairs of up-regulated genes showed raw Fisher-z rewiring in 40% of cases in one dataset; each of these is why the corresponding verdict uses a matched baseline.

## Samples behind the meta-analysis

The group sizes of each dataset's contrast appear next to its checkbox. After the run, a table under the summary line lists, for each dataset, the samples compared (B vs A), further samples fitted only for the variance estimate (when the variance comes from all groups), samples left out for a missing covariate value, the number of distinct subjects when a patient or pair column is adjusted for, the dataset's median share of the pooled weight, and the number of genes it contributed. The **samples** column of the gene table gives the samples behind each gene's pooled estimate, which is smaller for genes measured in fewer datasets. A dataset's weight reflects its precision and, with random effects, the heterogeneity between datasets, so it is not proportional to its sample size.

## Pathway enrichment of meta-analysis genes

Below the meta-analysis summary, **pathway enrichment of** tests a list of meta-analysis genes for over-representation in the gene-set library chosen in the Enrichment tab. The list is either the shared signature or the top N genes ranked by meta p or by |pooled log₂FC|, in the up, down or both directions. **Run in Enrichment tab** opens the Enrichment tab with the list selected as the query (it also appears there as *Meta-analysis genes*). The universe is the genes tested in the meta-analysis that are annotated in the library. The hypergeometric test treats genes as independent; genes at the top of a ranking are co-regulated, so its p-values are optimistic, and the result describes what the list is made of rather than providing an independent test.

## Multi-gene forest plot

After the meta-analysis has run, type gene symbols into **multi-gene forest** (or leave the box blank for the genes of the shared signature with the largest pooled effects) and click **Draw**. Each gene is a row: coloured points are each dataset's own log₂FC, and the black diamond is the pooled estimate with its 95% confidence interval. The axis always includes 0. Genes not in the meta-analysis (not found, or measured in too few datasets) are listed under the plot.

## Homology-family dot plot

Type member symbols of either species into **family dot plot** (for example `SERPINB3, S100A7, LCE3`; a prefix such as `LCE3` matches every LCE3 member), or leave the box blank for the families significant in the same direction in the most datasets, then click **Draw**. Every loaded human or mouse dataset with a two-group contrast is shown, one coloured point per dataset; filled points are significant at the DE FDR, open points are not. Each family is scored within each dataset on that species' own members (mean z of the members found), so raw values are never merged; the FDR in this plot is corrected across all families scored in that dataset. The plot does not depend on the X/Y pair chosen for **Score families**. This is the view for gene families with no one-to-one ortholog, such as the S100A7, SERPINB3/4 and LCE3 families.
