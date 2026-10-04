# Changelog

## 1.0.0-rc.1 (2026-10-03)
Release candidate for version 1.0 (design freeze). The default meta-analysis changes in this release: re-run meta-analyses made with earlier versions before reporting them.

Meta-analysis
- CHANGED: the default meta-analysis is now REML random effects with the Hartung-Knapp-Sidik-Jonkman test (t distribution on k-1 degrees of freedom). On null data the previous default (DerSimonian-Laird, z test) gave genes at FDR ≤ 0.05 more often than chance allows with three to five small datasets; the new default stayed within the binomial bound in every tested configuration. On the four psoriasis cohorts the shared signature changes from 2,151 to 1,773 genes and keeps all 20 canonical psoriasis genes. DerSimonian-Laird, the z test and a truncated Hartung-Knapp test remain available. Sessions saved with earlier versions restore the model they were run with.
- NEW: REML estimate of the between-dataset variance (Newton iterations on the restricted likelihood, with the boundary at zero checked), agreeing with metafor::rma(method = "REML") to 1.1e-9 in tau² over 19,518 genes.
- CHANGED: confidence intervals in the meta-analysis table, CSV export, single-gene and multi-gene forest plots and the hub member forest plot use the standard error and critical value of the selected test (t on k-1 df with Hartung-Knapp); the prediction interval uses the same standard error, as metafor does.
- NEW: meta-analysis settings (model, test, minimum datasets, thresholds, prediction-interval rule) are saved in sessions.
- Validation: two self-test components added (REML with Hartung-Knapp, and DerSimonian-Laird with Hartung-Knapp, against metafor), 21 in total.

Differential expression
- NEW: interaction contrasts. A new card in the DE tab tests the difference of differences (B2 − A2) − (B1 − A1) between four groups of a combined factor (for example treatment × genotype or treatment × time), with the moderated t or voom, covariates, and the variance from all groups or the four groups only. It reports both simple effects, a volcano plot, a sortable table linked to the Gene Explorer, a CSV export, a methods paragraph and over-representation analysis of the significant genes. Results agree with limma to 5e-11 on GSE143688 (VALIDATION_REPORT.md).

Privacy
- NEW: both builds carry a Content-Security-Policy that forbids the page from opening any network connection or loading content from another server; the online build may load only Plotly.js from its CDN, now pinned by a subresource-integrity hash. The offline build allows no server. A new CI test confirms in Chromium, Firefox and WebKit that a full analysis makes no outside request and that attempts to reach another server are refused. The Privacy page describes how to check this in the browser.

Fixes from the release-candidate check
- FIXED: Welch t-test: a gene whose values are constant up to rounding in both compared groups (for example zero counts in every compared sample) received a small p-value computed from floating-point noise. Such genes are now untestable (t = 0, p = 1), with the rule used by R's t.test.
- FIXED: meta-analysis I² under REML is now computed from the fitted τ² (100 τ²/(τ² + s²), as metafor), so it is consistent with the displayed τ²; for DerSimonian-Laird and the fixed-effect model it is unchanged (from Cochran's Q). Agreement with metafor: 8e-7 percentage points over 19,401 genes.
- FIXED: the meta-analysis methods text stated 'a Wald z test' when the Hartung-Knapp test was selected, and the p-column tooltip always read 'meta z test'; both now name the selected test.
- FIXED: the set-overlap matrix and ridge plot titles contained click instructions, which appeared in exported figures; the instructions are now shown under the plots only.
- FIXED: with angled labels and the side margin at its limit, long labels could still run past the figure edge and be hidden; they are now shortened to the space available.
- FIXED: figures that only show a message (no biotype column, too few groups or contrasts) no longer offer image or data export.
- FIXED: the module × trait CSV export wrote a header-only file before any analysis was run; it now asks for the analysis first.
- FIXED: the multi-dataset summary CSV wrote 'NaN' for datasets in which a gene is not matched; these cells are now empty, as in the meta-analysis CSV.
- FIXED: the multi-gene forest plot axis is now symmetric about 0, like the other forest plots.
- NEW: with the Hartung-Knapp test and 2-4 datasets, the meta-analysis summary states the test's power in benchsiDE's simulations and the trade-off with the z test, which has more power but does not control the FDR when true effects differ between datasets. The help text of the test selector says the same, and describes the truncated form as a sensitivity check.
- CHANGED: the leave-one-dataset-out table adds, for each omitted dataset, the number of signature genes whose pooled estimate keeps its direction at nominal p ≤ 0.05; this does not depend on the correction across genes, which with the Hartung-Knapp test and 2-3 remaining datasets leaves very little power. On the four psoriasis cohorts 84-94% of the signature genes are retained by this measure.

Testing
- NEW: tests/engine/, a command-line test engine that runs the application on public GEO datasets (downloaded and checksum-verified by fetch_data.py) and checks interaction contrasts against limma, the cross-dataset functions, figure layout, quality-control warnings, Discovery captions and session round trips. CI adds a privacy test and interface tests for the new features in Chromium, Firefox and WebKit.

## 0.23.0-beta (2026-10-01)
Results change in this release for multi-group designs, meta-analysis, ORA and GSEA: re-run analyses made with earlier versions before reporting them.

Differential expression
- CHANGED: with more than two groups, the moderated t and voom now estimate the variance from all groups (design ~ group, the comparison as a contrast), as in limma's standard analysis. Earlier versions fitted only the two compared groups. The previous behaviour is available under "variance from: these two groups only". Covariates and blocking factors are fitted in the same full model.
- FIXED: a subject or patient column with numeric codes (1, 2, 3, ...) was fitted as a continuous slope instead of a blocking factor. Columns named like identifiers (patient, subject, donor, pair, individual, id) are now categorical, and each covariate has a continuous/categorical selector.
- FIXED: a covariate column that duplicated a group indicator (for example a patient with a single sample) made the design singular; aliased columns are now dropped, as lmFit does.
- FIXED: p-values for very large degrees of freedom used a normal approximation from 10⁶ df upward; the exact t and F tails are now used for all finite df.
- FIXED: the methods text described the analysis incorrectly in several places (contrast direction, samples and groups fitted, voom transform, covariate coding, ANOVA and pattern sentences when those analyses had not been run). It is now generated from the fitted model.
- NEW: the DE table export has a settings header, columns named by group, and AveExpr, SE, t and the significance call. The MA plot and the selection table use AveExpr.
- NEW: a warning is shown when the smaller group has 3 or fewer samples; at that size FDR control is approximate.

Gene-set analysis
- FIXED: ORA used every filtered gene as the universe, including genes not annotated in the library, which overstated enrichment. The universe is now the filtered genes annotated in at least one set of the library; pattern and heatmap clusters use the clustered genes as the universe.
- CHANGED: GSEA with sample permutation reports results at FDR ≤ 0.25 (shaded), the GSEA convention. With gene permutation (preranked, used when sample permutation is not possible) no FDR is reported, because that null gave false sets on data with no true difference; the table is a ranking only.
- FIXED: FRY with covariates treated a missing covariate value as a level; such samples are now left out, as in the DE test. FRY and the enrichment matrix use the same full model as the DE test.
- FIXED: the gene-set overlap panel failed (NaN) when two sets had no genes after mapping.
- FIXED: the GSEA plot title said "Preranked" for sample-permutation runs.

Discovery, clustering and patterns
- FIXED: the CAMERA residual degrees of freedom counted samples from groups that were not tested.
- FIXED: nomenclature stems that are not gene families (open reading frames, LINC, MIR, SNORD, KIAA, FAM, TMEM, ZNF and similar) were screened as families.
- NEW: the Discovery summary states the permutation floor, the low power of CAMERA at small group sizes, and that DE covariates are not used in the screen; the figure caption names the test that was run.
- FIXED: k-means results depended on row order; rows are now ordered canonically and 20 restarts are used. Singleton clusters now have silhouette 0, as in cluster::silhouette.
- NEW: a warning when auto-k is used with three groups, where silhouette favours k = 2 whatever the number of patterns.

Cross-dataset comparison and meta-analysis
- CHANGED: the meta-analysis now carries each dataset's evidence as its own moderated t supports, by converting the t (on its degrees of freedom) to the normal quantile with the same tail probability. Treating small-study standard errors as known gave false shared genes on data with no true difference. On six sets of three 3-versus-3 null splits of healthy skin (GSE54456), genes at meta FDR ≤ 0.05 fell from 9–113 to 0–28 (random effects) and from 65–435 to 1–28 (fixed effect).
- CHANGED: the prediction-interval criterion is off by default; leave-one-out robustness is computed without it; "gene must be matched" defaults to a majority of datasets; the summary reports genes not tested because they were measured in too few datasets.
- FIXED: the mouse-human ortholog table included genes from many-to-many homology classes as if they were one-to-one; mitochondrial genes (mt-Co1 → MT-CO1 and others) were missing.
- FIXED: regrouping a comparison dataset by a different design column did not re-run its filter and normalization.
- FIXED: the consensus table counted only genes measured in every dataset; significance is now counted among the datasets in which each gene is measured.
- FIXED: duplicate-dataset detection missed renamed or partial copies; it now matches sample IDs or identical library sizes.
- FIXED: signature transfer called a signature "replicated" when it separated the groups in the opposite direction.
- NEW: the RRHO map is labelled as descriptive (its p-values assume independent genes).
- FIXED: the homology-families card used "session" instead of the session's dataset name.

Gene Explorer, QC, power and input
- FIXED: the Kruskal-Wallis p-value lost accuracy below about 10⁻¹⁰ and returned 0 below about 10⁻¹⁶.
- FIXED: with a batch-adjusted display, Gene Explorer tests were computed on adjusted values; tests now always use unadjusted values.
- FIXED: the power calculator used the first two groups rather than the DE comparison, and included Poisson noise in the biological CV.
- FIXED: input units were detected from the first rows only; log₂(x + 1) input is now converted back as 2^x − 1, so the CPM filter is correct for such files.
- NEW: warnings for non-integer count input, missing or non-numeric cells, and duplicated gene IDs; a clearer message for decimal-comma files.
- FIXED: design files with quoted fields (R write.csv) were refused; GEO prefixes are removed when most values carry them; samples without a design label are excluded by default.
- FIXED: restoring a session made with different filter settings did not re-run the analysis.
- FIXED: the Excel/Prism export of the RLE plot was empty.
- FIXED: the UpSet plot did not state that only the largest intersections are shown, and its export lacked set labels. QC detection plots now say they count filtered genes.

Validation: the external validation has 12 new checks for the full-design analysis (moderated t, p-values, prior, covariate and blocking fits, voom and FRY against limma), bringing it to 48. The self-test references were regenerated for limma-trend and edgeR log-CPM. CHANGED: the header controls were crowded and wrapped unevenly. Export size, fonts, colours and the three axis-label settings are now in one **Figure settings** menu; the action buttons share one size and are grouped (figures and report; sessions and data; validation). The dataset description is shortened with an ellipsis when space is short (full text on hover).
- CHANGED: the validation of statistics opens as a centred dialog: an overall result, the browser and run time, the explanation folded under "What this test does", and a table with aligned numeric columns and a pass/FAIL badge per component; click a row for its value pairs. Close with the button, Esc or a click outside.
- NEW: pathway enrichment of meta-analysis genes. Under the meta-analysis, choose the shared signature or the top N genes (ranked by meta p or by |pooled log2FC|), up, down or both, and click Run in Enrichment tab; the same list is also available as a query source in the Enrichment tab. Over-representation is tested against the genes tested in the meta-analysis that are annotated in the library, so only genes that could have been selected are counted. Hypergeometric p-values agree with scipy.stats.hypergeom to 5.5e-11 (relative) on the four-cohort psoriasis meta-analysis.
- CHANGED: the homology-family dot plot now shows every loaded human or mouse dataset with a two-group contrast (the session and all comparison datasets, duplicates excluded), each family scored within each dataset on that species' own members, with BH correction across the families scored in that dataset. It no longer depends on the X/Y pair chosen for Score families, so it also works when the pair is changed or both datasets are of the same species. A blank selection draws the families significant in the same direction in the most datasets. Clicking a family shows its per-sample module scores in every dataset; the data export has the score change, FDR and members found for each dataset.
- FIXED: module heatmap export: long group names overlapped above their columns; each name is now wrapped to fit over its own block of samples (at underscores, spaces and case changes) and alternates between two rows if it still does not fit. The colour bar ticks overlapped on short figures; the bar now has a minimum height and fixed ticks. The module-score strip had been coloured on its own range while sharing the matrix's colour bar, so its colours did not match the bar; both now use one symmetric scale, stated in the caption. Every sample is labelled when there is room (vertical labels from about 8 px per column); otherwise the labels are left out and the caption says so.
- FIXED: with angled axis labels, a long label near the end of the axis could run past the edge of the figure, and Plotly then hides it entirely (for example the last group of a five-group bar plot had no label). The right margin (left margin for negative angles) now grows until every label fits, on screen and in exports.
- NEW: samples with almost no reads (under 1% of the median library size, or under 100,000 reads) are named separately in the notes above the dashboard, with their read count and a one-click link to exclude them; they appear as outliers in every plot and test.
- NEW: the meta-analysis reports the samples behind each dataset: the group sizes next to each dataset's checkbox, and after the run a table with the samples compared in each contrast, further samples that only inform the variance estimate, samples left out for a missing covariate value, distinct subjects when an identifier column is adjusted for, each dataset's median share of the pooled weight and the number of genes it contributed, with the total. Each gene's row in the table and CSV gives the number of samples behind its pooled estimate; the heatmap columns, forest plots, CSV header and methods text give the group sizes.
- NEW: multi-gene forest plot in the meta-analysis card: for a list of genes (or the top of the shared signature), each dataset's own log₂FC and the pooled estimate with its 95% confidence interval, genes ordered by pooled effect. Click a gene to open it in Gene Explorer; the data export gives one row per gene with every dataset's value, the pooled estimate, its CI and the meta FDR.
- NEW: homology-family dot plot in the homology-families card: chosen families (any member symbol of either species) or the strongest concordant families, with each dataset's family score change; filled points are significant, open points are not. Click a family for its per-sample scores; the data export gives scores, FDRs and members for both species.
- FIXED: changing the grouping column of a comparison dataset re-read its genes with whatever species was currently selected for adding datasets. After a dataset of another species had been added (for example a mouse dataset next to human ones), regrouping a human dataset left none of its genes matched, and it silently dropped out of the meta-analysis, forest plots, concordance and every other cross-dataset view. The dataset's own species is now used.
- NEW: load timing. After Analyze, a line above the dashboard gives the time this browser took to read and parse each dropped file, to put the review step on screen, and from Analyze to the first display of the results, with the total (time spent by the user on the review step is not counted); **copy** puts the timings and browser details on the clipboard as JSON.
- NEW: timing benchmark (Validate statistics → Timing benchmark…). Times a fixed analysis in the current browser on deterministic synthetic datasets (20,000 genes; 6, 24, 96 and 192 samples; one warm-up and 3–10 repetitions; median, minimum and maximum per step) and shows every run of every step (warm-up included) next to the medians. Results download with the browser and hardware details as JSON, and every run as CSV (one row per dataset, run and step, in milliseconds as measured). The same datasets can be exported as TSV files for other tools. `npm run bench` runs it in Chromium, Firefox and WebKit through Playwright.
- FIXED: the browser validation test checked the demo fixture against a wrong gene count and would have failed in every browser. CI now also runs on every branch and from the Actions tab (Run workflow), and keeps the test report.
- NEW: figure preview. The camera button on every plot, and the module-heatmap, gene-panel, GSEA and summary-figure exports, now open a preview window first: the figure as it will be saved, its pixel size, a choice of SVG or PNG (where both apply) and an editable file name. Save writes the file (and its caption file, where there is one); Close, Esc or a click outside the window cancels. Untick Figure settings → Preview to save at once as before.
- NEW: Sample traits links to Gene Explorer: an arrow next to every gene in the table, an Open in Gene Explorer button under the gene plot (the back button returns to the Traits tab), + Add to comparison, and clickable gene lists for a module × trait cell or a CAMERA gene set, ordered by each gene's association with the trait.
- NEW: the sample dendrogram colours each sample by group or by any design column (2–12 levels), with a legend, and can cut the tree into clusters (as many as the levels by default, or 2–10). Branches are coloured by cluster, a table cross-tabulates clusters against the chosen column with the adjusted Rand index, and samples placed in a cluster dominated by another group are listed. Heights, leaf order and cluster assignments are identical to R hclust (average linkage) and cutree; the adjusted Rand index matches scikit-learn.
- NEW: Sample traits tab. Numeric per-sample measurements (clinical scores, IHC counts, thickness, serum levels) can be analysed as variables of interest: genes associated with the trait (limma, effect per unit and per SD, Cook's distance flag), gene sets (CAMERA on the trait coefficient), a module × trait matrix (module eigengenes of hub neighbourhoods, gene sets or Discovery modules) and trait × trait correlations. Options: transformation (log₂(x + 1), logit, rank), group selection or adjustment, covariates, and repeated samples per subject (block means for subject-level traits, mixed model with duplicateCorrelation for traits that vary within subjects). A warning appears when the trait mostly differs between groups, when one sample has high leverage, and when a subject column is available but not selected. Settings are saved with the session; the methods text describes the model.
- Validation: six new external checks against limma (trait moderated t and p, duplicateCorrelation, blocked lmFit, CAMERA) and WGCNA (module eigengene), bringing the external validation to 59. Calibration results are in VALIDATION.md.
- NEW: hub across datasets (Compare tab, opened with **Across datasets** in the neighbourhood panel). For the selected hub neighbourhood, in every comparison dataset: (1) preservation (WGCNA Zsummary, with a verdict calibrated against random gene sets, gene sets with the same differential expression, and neighbourhoods of random genes); (2) consensus hubs across all datasets by robust rank aggregation of connectivity ranks over one common gene universe; (3) rewiring between the contrast groups (Fisher z, calibrated against gene pairs with the same differential expression, with replication counted across independent datasets); (4) the members in the meta-analysis, pooled with and without the discovery dataset. Results export to CSV.
- Validation: three new external checks (soft connectivity against R, the seven preservation statistics against WGCNA 1.74's own code, and RRA scores against RobustRankAggreg), bringing the external validation to 53. Calibration results on public psoriasis data are in VALIDATION.md.
- NEW: hub neighbourhood panel (Co-expression tab). Clicking **neighborhood** for a hub opens one view of its most correlated genes: an ego network (the hub in the centre; edge width = |r|, edge colour = sign; nodes coloured by log₂FC and outlined when differentially expressed), a heatmap of the members' z-scores with samples grouped by condition, a per-sample neighbourhood score by group, and the gene-set enrichment of the neighbourhood. Correlation can be computed across all samples or within groups (partial correlation given group, which removes correlation produced only by shared group differences). The neighbourhood can be scored in every comparison dataset through signature transfer and exported for Cytoscape (node and edge tables) or as a GMT gene set.
- NEW: Gene Explorer can draw each gene as a bar (mean ± SEM or SD), box (median, quartiles, 1.5 × IQR whiskers) or violin (kernel density with median and quartiles) plot, with the sample points on every type. The caption under the plot states what is drawn and warns when the smallest group is too small for a box or violin summary.
- CHANGED: the multi-gene comparison is drawn as grouped bars by default (one bar per gene within each group, mean ± SEM or SD, with sample points), with grouped box, grouped violin and the previous profile lines as options; profile lines remain the default for ordered (time-course) designs. A **group by** selector arranges the chart either by condition (genes side by side within each condition) or by gene (conditions or time points side by side within each gene, coloured as in the rest of the app). The z-score option now scales each gene across samples, not across group means. The Excel export gives a Prism Grouped table in the arrangement shown (rows: the x-axis categories; data sets: the coloured series), with the samples as replicate sub-columns.
- NEW: once a gene pair is plotted in the Co-expression tab, buttons below the statistics open either gene, or both side by side, in Gene Explorer; the Back button returns to the pair.
- FIXED: after **Load different data**, differential-expression results of the previous dataset were kept until a contrast was run on the new one. In that window the hub table, the connectivity-vs-DE scatter, the network colouring and the signature-transfer gene lists could read the previous dataset's log₂FC and FDR at the new dataset's gene positions. Loading or re-analyzing data now clears the DE results and the hub screen.
- Validation: two new external checks (partial correlation given group and its p-value, against residual correlation and lm in R).
- FIXED: numbers in tables (co-expression r, slope, p and FDR, and other result tables) could break across lines when the table was narrow (for example "1.00" / "0"). Numbers now never break; only the first column may break inside a long identifier, and a narrow table scrolls instead. The co-expression tables use a compact layout with right-aligned numbers.
- FIXED: in the self-test panel, short cells (values, tolerance, pass) wrapped one letter per line in narrow windows. The browser test now expects the 19 self-test components and pins the demo results of the full-design voom fit (verified in R).

## 0.22.0-beta (2026-09-30)
Results change in this release: re-run analyses made with earlier versions before reporting them.
- CHANGED: count data are now log-transformed as edgeR log-CPM (prior count 2, scaled by library size) instead of log₂(CPM + 1). log₂(CPM + 1) shrank the fold changes of lowly expressed genes toward 0 (GSE54456: IL17A log₂FC 0.95 instead of 3.6), so they failed the fold-change cut-off. Already-normalized input (CPM, TPM) still uses log₂(x + 1).
- FIXED: the default moderated t was labeled limma-trend but used a constant prior variance. It now fits the prior as a function of average expression, matching `eBayes(trend = TRUE)`. The moderated F-test and the comparison-dataset test use the same model.
- FIXED: when the empirical-Bayes prior degrees of freedom were infinite, p-values used 10⁶ degrees of freedom instead of limma's cap at the pooled residual degrees of freedom.
- FIXED: digamma and trigamma functions were accurate to about 2×10⁻⁹; now to about 10⁻¹². Agreement with limma improved from about 10⁻⁸ to 10⁻¹².
- CHANGED: GSEA now estimates significance by permuting sample labels (Subramanian et al. 2005), with FDR from the pooled normalized enrichment scores. Gene permutation (preranked) gave false positives on data with no true difference (12–37 of 50 Hallmark sets at FDR ≤ 0.05 on random splits of normal skin). Gene permutation remains available for ranking and is used automatically, with a warning, when there are fewer than about 7 samples per group or covariates are selected.
- FIXED: the GSEA table had no FDR column.
- NEW: comparison datasets keep every design column. Each slot lets you choose the column that defines the groups and columns to adjust for or pair by, fitted as in the session's covariate-adjusted analysis.
- CHANGED: comparison datasets are filtered with filterByExpr (previously CPM ≥ 1 in 2 samples), as the session is.
- FIXED: samples with a missing value (NA or blank) in a selected covariate were treated as a separate level. They are now left out of that fit, as in limma.
- NEW: the co-expression table and export include the Benjamini–Hochberg FDR across all genes tested against the query gene.
- FIXED: gene families in Discovery split names at the wrong digit: S100A8 was grouped under "S1". Families are now formed by removing the number after the last letter (S100A8 → S100A).
- FIXED: with two groups, Discovery reported effects as the first group relative to the second, regardless of the DE contrast. Effects now follow the DE tab's contrast.
- FIXED: hub, network, forest-plot and homology-family highlights used FDR ≤ 0.05 whatever FDR was set in the DE tab. They now use the DE tab's FDR.
- CHANGED: when no Discovery module passes FDR, the nominal list now states how many modules would reach p ≤ 0.01 by chance.
- Validation: references regenerated for the new defaults; self-test 19/19 and external validation 36/36. New checks: full-matrix limma-trend on GSE54456, GSEA enrichment scores against fgsea, and calibration and known-biology controls (VALIDATION.md).

## 0.21.0-beta (2026-09-30)
- NEW: meta-analysis across datasets (Compare datasets tab). Each dataset's own log₂ fold change and standard error are pooled per gene by a random-effects (DerSimonian–Laird) or fixed-effect model, with Cochran's Q, I², τ², a z test and FDR. Expression values are never combined. Outputs:
  - a shared signature: genes significant after pooling, with a minimum pooled effect, changed in the same direction in every dataset, and (by default with 5 or more datasets) with a 95% prediction interval that excludes 0; a maximum I² is optional and off by default;
  - a meta volcano coloured by heterogeneity;
  - a heatmap of each dataset's and the pooled log₂FC;
  - a sortable table with 95% confidence and prediction intervals;
  - a leave-one-dataset-out robustness table;
  - CSV export of all genes, and GMT export of the signature.

  Results agree with metafor::rma on all 19,464 genes of a two-dataset public example (VALIDATION.md). Two self-test components were added (18 in total).
- CHANGED: forest plots show the pooled estimate as a diamond once the meta-analysis has been run.
- NEW: network view for hub genes: the top hubs and their strongest partners, drawn in a circular layout (ordered by clustering, curved edges) or a force-directed layout. Edges are Pearson correlations above a chosen threshold, limited to each gene's strongest edges. Nodes are sized by connectivity and coloured by log₂FC or connectivity. Click a node to open the gene. The edge and node tables export as CSV for Cytoscape.
- NEW: dataset names. Every comparison dataset, and the session's own dataset, has an editable name that is used in all plots, tables, exports and the methods text. GEO files default to their accession (for example GSE54456). The session name is saved with the session.
- Methods text: sentences for the meta-analysis and network view, with citations.
- NEW: genes selected with a box or lasso on the volcano or MA plot are listed in a sortable table (log₂FC, group means, p, FDR, significance call). The selection is kept when switching between the two plots and is saved with the session. It can be exported to CSV, copied as a gene list, drawn as a heatmap (new heatmap gene source), or sent to the Gene Explorer comparison (up to 10 genes).
- NEW: colour options (header, **colours**): group palette (default, Okabe–Ito, Tol bright, ColorBrewer Dark2, greyscale), up/down colours, heatmap scale (RdBu, PuOr, PiYG, BrBG, blue–white–red) and continuous scale (viridis, cividis, magma, inferno, greyscale), plus a colour picker for each group. The choice applies to every plot and export and is saved with the session; the default output is unchanged.
- FIXED: table headers broke inside words (for example "GSE205748_raw_co / unts", "sig in" shown as "g in"). Headers now break only at spaces and after _ . - /, and wide comparison tables scroll horizontally instead.
- FIXED: the forest plot could draw over the table below it; every plot with a computed height now gets a container of that height. The bottom forest-plot row no longer sits on a horizontal zero line.
- FIXED: gene names in the consensus and meta-analysis tables were squeezed onto several lines.
- NEW: a comparison dataset with the same samples as the session or as another comparison dataset (the same data loaded twice) is flagged in its slot and left out of the meta-analysis by default.
- CHANGED: no fixed limit on the number of comparison datasets (previously 8). The slot list shows the memory in use, warns above 1.5 GB and refuses a dataset that would take the total above 3 GB. Comparison datasets no longer keep an unused copy of their raw counts, roughly halving their memory; all comparison results are unchanged.
- FIXED: large datasets made the page unresponsive. With 171 samples and 19,518 genes (GEO GSE54456), clicking Analyze froze the browser for about 85 s, and each change to the included samples froze it for about 40 s. The main cause was the MDS plot: every sample pair sorted all gene distances to find the 500 largest. It now selects them in linear time, and the result is identical to before. The same analysis now takes about 5 s, and changing the included samples takes about 4 s.
- FIXED: FRY gene-set testing was slow with many samples. The eigenvalues it needs are now computed by tridiagonal QL instead of Jacobi rotation, and the per-set matrix products use typed arrays. On GSE54456: Hallmark 7.0 s → 0.6 s, GO Biological Process (7,608 sets) 250 s → 7 s. p-values are unchanged to within 1e-10 (relative); agreement with limma::fry is unchanged (VALIDATION.md).
- NEW: progress bar for long computations: analysis, TMM normalization, MDS, the quality-control plots, FRY, GSEA and the multi-contrast enrichment matrix. The bar shows the fraction of work actually completed (for example, sample pairs done in MDS, samples normalized in TMM, gene sets tested in FRY, weighted by set size), the current step, elapsed time and an estimate of time remaining. The page repaints about 10 times per second while it runs. On large datasets, sample-inclusion changes and session restores also show it.

## 0.20.0-beta (2026-09-29)
- FIXED: long group names were cut off at the edge of plots and exports. Axis tick labels are now never cut off. By default, long labels wrap onto several lines, or are angled when wrapping is not enough, and the figure margin grows to fit them. New header controls set the label orientation (auto, horizontal, 45°, vertical), the label font size, and the treatment of long labels (wrap, full length, or shortened with an ellipsis). The settings apply to every plot and export and are saved with the session.
- NEW: every plot has a **Download data** button in its toolbar. It saves an Excel workbook: a sheet in GraphPad Prism table layout (Column, Grouped or XY, chosen from the plot type), a sheet listing every plotted value with its series, label and error values, and a notes sheet with axis titles, the Prism table type and import instructions. Where points carry sample names, a sheet of sample IDs is included.
- CHANGED: forest plots always include 0 on the x axis, and the axis is symmetric about 0, so effect sizes can be compared by magnitude.
- CHANGED: homology-family comparison works between any two uploaded datasets of different species (for example two comparison datasets), not only between the session and one dataset. The table, scatter plot, per-sample panels and CSV export are labeled with the chosen datasets.
- NEW: Gene Explorer shows the selected gene in every comparison dataset. Each dataset gets its own panel with its own groups, normalization and statistics (ANOVA, Welch tests with Holm adjustment, and the dataset's selected contrast). Matching is by symbol, by ortholog, or by custom ID map. A symbol can be typed to look up a gene that is not in the session. Values are never merged across datasets.
- CHANGED: the co-expression tables show the slope (OLS, gene on query gene) and the p-value for r = 0 beside r. Click a column header to sort by gene, r, slope or p. A selector sets how many genes are shown (30, 100 or 300 per direction), and **Export all (CSV)** writes every gene.
- NEW: **Design check**, shown in the review step and under the sample list. For each design factor it shows sample counts by group and reports, for each pair of groups, whether the factor can be adjusted for or is fully confounded with the comparison. For example, in GSE186063 all normal-skin samples are from ankylosing spondylitis patients, so diagnosis cannot be adjusted for in comparisons with normal skin, but it can in lesion vs non-lesion. The check lists the options for a confounded comparison.
- NEW: samples can be included or excluded by factor level under the sample list (for example, analyze psoriasis patients only).
- NEW: every pair of design factors is offered as a combined grouping factor (for example Type × Diag, giving groups such as "lesion · Psoriasis").
- NEW: design columns can be added and edited in the review step, for example a patient column for paired samples. Added columns are available as covariates and are saved with the session.
- FIXED: design values exported from GEO with a characteristic prefix (for example "diagnosis: Psoriasis") are read without the prefix.
- The error for a covariate that is confounded with the comparison now refers to the Design check.
- Validation: the external validation has three new checks, for the co-expression table (Pearson r, slope and p against stats::cor, stats::lm and stats::cor.test), bringing it to 36.

## 0.19.5-beta (2026-09-25)
- FIXED: running **Validate statistics** left the header summary, the methods filter sentence and the QC notes describing the self-test's synthetic dataset (6 samples, 2,000 genes) in place of the loaded data. The statistics themselves were restored correctly. The self-test also reset the normalization selector, the batch-adjusted display setting and a custom group order. It now saves and restores every display the analysis step rewrites; content panels are moved aside as live elements, so their click handlers keep working.
- New CI test: the self-test must leave the analysis displays unchanged and a control inside a content panel must remain functional.

## 0.19.4-beta (2026-09-24)
Validation release.
- VALIDATION.md rewritten. It now reports only checks that can be reproduced from this repository: the in-app self-test (16 components on a deterministic synthetic dataset) and a new external validation on the public demo dataset GSE63310 (33 checks against edgeR, limma, RRHO, RNASeqPower and stats). The previous record cited datasets that are not distributed with benchsiDE; it has been removed, together with the corresponding figures in this changelog.
- New: `validation/reference.R` generates the R reference values; `validation/compare.js` runs the comparison inside the application; `tests/validation.spec.js` runs it in CI in three browsers.
- FIXED: FRY with a continuous covariate failed with "No residual degrees of freedom". FRY coded every covariate as categorical indicators, while the differential-expression model fits a numeric covariate with more than two distinct values as a slope. FRY now uses the same rule; results match limma::fry(design = ~covariate + group).
- FIXED: RRHO maps were floored at -log10 p = 300, so strongly concordant comparisons were flattened to a constant. The hypergeometric tails are now kept on the log scale; values match the RRHO package beyond this range.

## 0.19.3-beta (2026-09-24)
Licensing release. No statistical output changed.
- Licenses of all embedded third-party data verified against the providers' current terms: MGI and RGD data are CC BY 4.0, HGNC data are CC0 1.0, MSigDB collections are CC BY 4.0, Plotly.js is MIT.
- The license notice now travels with the application file: a comment at the top of `index.html`, and a **licenses** panel in the header listing each resource, its license and how it was modified (column subsets, re-encoded).
- Loading built-in annotation now shows the data license under the citation.
- THIRD_PARTY_NOTICES.md rewritten; the MGI and RGD redistribution items are resolved.

## 0.19.2-beta (2026-09-24)
Documentation and language release. No statistical output changed.
- Interface text: all em-dashes removed from interface text, tooltips, headings, messages, methods text, report text and gene-set library labels (the characters, their escape sequences and their HTML entities), and replaced with the punctuation each sentence requires. Conversational and promotional wording removed. Empty table cells in the report now show an en-dash.
- The auto-generated methods paragraph was edited by hand; its content and citations are unchanged.
- The default SVG export file name is now `benchside_plot`.
- New user guide (`docs/user-guide.md`, compiled to `dist/benchside-guide.html`): eight step-by-step tutorials, with the demo-dataset results at each step reproduced independently in R.
- New knowledgebase (`wiki/`, 29 pages, published to the repository wiki with `scripts/publish_wiki.sh`). It replaces GUIDE.md.
- Test suite corrected: it expected 11 self-test components (now 16) and old panel text, and it pinned the demo results under the CPM filter while the default is filterByExpr. Both filter modes are now pinned, with values reproduced in R.
- New workflows: `pages.yml` deploys www.benchside.org from main; `release.yml` builds, tests and publishes a release when a version tag is pushed. `dist/` is no longer committed; its files are built in CI.
- README, CONTRIBUTING, CITATION.cff, THIRD_PARTY_NOTICES.md and package.json brought up to date.

## 0.19.1-beta (2026-09-23)
Tool-wide margin and overlap audit.
- Method: all 51 plot targets exercised in the engine across 8 scenarios (2-group and 5-group designs; screen, paper 1.5x and poster 2x fonts; a hostile dataset with 50-character sample and group names), producing 320 captured layouts, and checked against Plotly.js layout semantics: axis tick labels, axis titles and annotations do not auto-expand margins; legends and colorbars do; an auto-positioned title's baseline sits at the centre of the top margin.
- Found at screen size: 7 plots where the title collided with a legend drawn above the plot (stability curves, family concordance, co-expression, set ECDF, library complexity, sex check, category plot); voom-trend annotation clipped by an 8 px top margin; category labels wider than their margins on the ridgeline, Jaccard, signature-transfer and sample-correlation plots; small left/bottom deficits on volcano, patterns and p-value histogram. Poster fonts and long names multiplied these.
- Fix, central: a layout guard in the single render wrapper every plot and export already passes through. It only raises margins (never lowers), measures text with the browser's real font metrics, uses the actual container width for tick-angle decisions, runs after font scaling (so poster fonts get poster margins), pins titles to the container top and stacks legends/annotations beneath them, truncates labels that exceed the margin cap with a middle ellipsis (head and distinguishing tail kept; uniqueness enforced; full names stay in hover), and grows the figure and its container when margins would crush the plot area below 160 px.
- Tables: long unbreakable identifiers now wrap inside cells instead of pushing tables past the card edge.
- Real bug found by the audit harness: FRY result rows carried no member indices, so every FRY run threw inside the set-overlap panel after the table rendered. FRY rows now carry their members (which also enables barcode/ECDF/ridgeline for FRY); the Jaccard panel guards the degenerate case.
- Result: 0 layout flags and 0 crushed plot areas in all 8 scenarios; statistical results unchanged (self-test 16/16).

## 0.19.0-beta (2026-09-22)
Cross-species results are now displayed graphically in addition to the family table:
- Family concordance scatter (Compare tab, Homology families): one point per many-to-many homology family: x = module-score group difference in the session, y = the same in the chosen comparison dataset, each computed entirely within its own dataset. Point size scales with matched members; red/blue = concordant significant (both FDR<=0.05), amber diamonds = significant but discordant; diagonal = identical effect; the strongest concordant families are labeled (staggered, quadrant-inward). Clicking any point opens the per-family evidence plot.
- Per-family evidence plot: side-by-side per-sample module-score boxes (session groups | dataset groups, z within each dataset, never merged), with member counts, deltas, and FDR per side. Also reachable from a "plot" link on every family-table row.
- Validated: scatter coordinates equal the table's numpy/scipy-validated deltas exactly; per-sample plotted scores reproduce both deltas to 1e-9; clear empties all plots; layout verified by spec-reconstruction render (label staggering and quadrant-inward offsets fixed from the render check before ship).

## 0.18.2-beta (2026-09-22)
- User guide added: GUIDE.md at the repository root (canonical, versioned with the code) and a standalone dependency-free HTML rendering at dist/benchside-guide.html, following the same single-file/offline contract as the tool, intended to be served at benchside.org/guide.html. Covers quick start, input formats, a tab-by-tab reference, which-test-when guidance, the validation/evidence rationale, reproducibility machinery, troubleshooting (each entry from a real user-reported failure mode this cycle), and FAQ. A "guide" link joins the app header; the in-app tooltip layer remains the contextual documentation.

## 0.18.1-beta (2026-09-22)
Four reported UI defects on the homology-families panel and self-test, all reproduced in the engine before fixing:
- Family scoring always used the FIRST cross-species dataset; with several loaded, the others were unreachable. A dataset selector now chooses which one is scored (the table header names it).
- No way to reset the family panel, and a failed or changed re-run left the previous table on screen. Every run now clears prior output first (including error paths), and a Clear button empties the panel.
- The family table overflowed the frame: rows with 29 gene chips and long family labels had no width or height limits. The table now lives in a scrolling container (max 440 px), member chips cap at 6 per cell with a "+N more (CSV has all)" note, and family labels wrap.
- The statistics self-test report could only be dismissed via a "close" link buried at the bottom; a top-right X now closes it from anywhere.
Family-score regression values unchanged; full regression passes.

## 0.18.0-beta (2026-09-22)
- Per-sample expression density overlays (Overview QC, collapsible): exact-sum Gaussian KDE with the nrd0 bandwidth and grid conventions of stats::density defaults (n=512, cut=3), one curve per included sample colored by group. Validated at machine precision against an exact R computation (bandwidth 4e-15, curve values 3e-10; R's own binned density() agrees only to ~1e-4, its documented binning approximation). Computed lazily and cached (~8 ms/sample).
- Enrichment-across-contrasts matrix (Enrichment tab): the FRY directional test run independently for every pairwise contrast, displayed as a sets-by-contrasts heatmap of signed -log10 FDR (red = up, blue = down; BH within each contrast; top 25 sets by best FDR shown, full matrix in the CSV). Clicking a cell opens that set's barcode and ECDF in that contrast. Each column is the single-contrast FRY computation (validated against limma::fry; see VALIDATION.md). Contrast selectors are restored after the sweep. Not offered for two-group designs, where the single-contrast table gives the same result.

## 0.17.0-beta (2026-09-22)
- DEG Venn diagrams (all-pairs card, multi-contrast panes): pick 2-3 contrasts and a direction mode (significant / up-only / down-only); circles drawn with exact disjoint region counts from the same contrast engine and thresholds as the all-pairs table. Clicking a region count lists that exact region's genes (clickable, paged); a CSV exports every gene with its membership pattern. Areas are deliberately not proportional, and the plot states this. Region counts were checked against an independent set partition; hidden when fewer than 2 contrasts exist.

## 0.16.0-beta (2026-09-22)
- Many-to-many homology families are no longer lost. 987 families (MGI homology classes merged by shared membership, e.g., SERPINB3+SERPINB4 <-> Serpinb3a/b/c/d) are embedded alongside the strict 1:1 table. A new "Homology families" panel in the Compare tab scores each family as a module WITHIN each dataset (mean z of member genes; Welch t on that dataset's own contrast; BH across tested families) and reports cross-species direction concordance; no arbitrary 1:1 gene pick is ever made, and raw values are never merged. Gene-level matching policy is unchanged (strict 1:1 or custom map). Family scores were checked against an independent numpy/scipy recomputation. Methods text updated.
- Logo: "benchsi" and "DE" are now two distinct looks (regular-weight vs blue extra-bold) per design direction.

## 0.15.2-beta (2026-09-22)
- FIXED: cross-species ortholog matching returned 0 genes whenever identifier spaces differed (reported with a human Entrez-keyed dataset, GEO GSE117405, and a mouse Ensembl-keyed dataset). Two causes: comparison datasets were decoded with the session's annotation map (which cannot decode another species' identifiers), and an undecoded session (annotation never enabled) offered no symbols for the symbol-keyed ortholog table. Comparison datasets are now decoded with their own species' built-in annotation (cached per species, independent of the session's). TP53 to Trp53 spot-checked.
- Low-match diagnostic: when under 1% of session genes match, the slot names both identifier spaces with examples (e.g. "bare numeric NCBI Entrez IDs" vs "Ensembl mouse gene IDs") and states the specific remedy (enable species annotation and re-run; or load a custom ID map) instead of failing silently.

## 0.15.1-beta (2026-09-22)
- The self-test now shows its evidence instead of asserting a verdict. The report opens with the methodology (deterministic seeded dataset and its exact recipe; the same functions that power the UI; no test-only statistics path; externally computed references with full version provenance: R 4.5.3 / edgeR 4.8.2 / limma 3.66.0 / statmod 1.5.2 / RNASeqPower 1.50.0 / scipy 1.17.1 / statsmodels 0.14.6). Every component row states what is executed, the exact reference command, how many values are compared, the max deviation, and the tolerance. It expands to a side-by-side table of browser value vs reference value vs delta at full precision. A "Download evidence (JSON)" button exports all 110 compared value pairs with labels, per-component descriptions, engine user-agent, and provenance: a machine-readable validation certificate per run.

## 0.15.0-beta (2026-09-22)
Reviewer-concern release (statistical-scope ceilings from the publishability review).
- Continuous covariates: a numeric covariate with more than two distinct values is now fitted as a continuous slope in the moderated-t model (limma lmFit(~group + x) semantics); validated against limma (see VALIDATION.md). Two-level numeric covariates keep the indicator path (identical t either way). The adjustment note names each covariate's treatment (continuous / 2-level / categorical).
- Paired and blocked designs documented and validated: supplying the subject or block identifier as a covariate performs fixed-effect blocking (limma ~group + subject); validated against limma (see VALIDATION.md). Stated in the methods text.
- Custom identifier mapping for cross-dataset comparison: a two-column TSV (session ID -> target ID) can be loaded and is applied before symbol/ortholog fallback: the workaround for multi-mapped ortholog families the strict 1:1 MGI table excludes (Serpinb3a-class), non-model species, and probe-to-symbol maps. Ambiguous source ids are dropped, not guessed. Slot labels and methods text report custom-map usage per dataset. Validated: renamed genes unreachable by symbol matching are recovered exactly and concordance controls are unchanged.
- Signature-transfer calibration anchored to the competitive gene-set testing literature in the methods text (Wu & Smyth 2012).

## 0.14.3-beta (2026-09-22)
- Ridgeline hover fixed: violin outlines hovered on every KDE vertex, producing an unreadable tooltip storm. Violins no longer hover; each set instead carries one median marker with a single clean tooltip (set name, gene count, median log2FC, FDR, click hint). Clicks on violins and markers both resolve to the correct set.
- Enrichment pane reordered so click results land where you are: ridgeline first, then the barcode/ECDF panels its clicks draw into (with a gentle scroll-into-view), then the Jaccard heatmap with its shared-genes list directly beneath, then the leading-edge matrix. Previously ridge clicks painted results two plots up, above the Jaccard.

## 0.14.2-beta (2026-09-22)
- Project URLs are live: benchside.org (hosted instance) and github.com/jyaron/benchsiDE (source). Added to the header (small links), the methods text, the report meta line, CITATION.cff (repository-code + url), README, and the manuscript availability statement.

## 0.14.1-beta (2026-09-22)
- FIXED: plot click handlers were re-registered on every redraw (Plotly .on accumulates listeners), so one gene click pushed multiple history entries and the back button needed one extra click per redraw. All plot click bindings now go through a bind-once/swap-handler helper; engine-verified (4 redraws + 1 click = exactly 1 history entry).
- Contrast comparison lists are complete: scrollable tables with incremental paging (100 rows at a time) for shared, X-only, and Y-only genes. No more fixed 25/12 caps.
- UpSet intersections are inspectable: clicking a bar lists the genes with exactly that membership pattern (clickable, paged) with a CSV download of the full list.

## 0.14.0-beta (2026-09-22)
Interactivity release: every click-through uses Plotly events, so rendered figures and exports carry no link styling.
- Navigation: back button in the Gene Explorer returns to the view you clicked through from (20-step history); all gene click-throughs route through it.
- Contrast comparison (DE tab): OLS regression line with 95% CI band (validated olsCI; slope/intercept/r/slope-p stated with a plain-language reading), ranked shared-hits table with BOTH contrasts' p and FDR per gene, X-only/Y-only strongest contrast-specific genes, top-10 shared genes labeled on the scatter, points and rows click through to the gene.
- Volcano grid: per-panel up/down counts in the panel titles; significant points carry gene names and click through.
- UpSet: intersection-count labels no longer clip at the plot top.
- Hub screen rebuilt around the driver question: a connectivity-vs-|t| scatter (upper-right = coordinated AND responsive; red = FDR<=0.05), DE columns and row highlighting in the table, and an expandable per-hub correlation neighborhood (top-10 partners with r, DE-colored, clickable) with a one-line reading of what a red vs grey neighborhood means. Explainer rewritten to state the question the screen answers.
- Jaccard heatmap: rows similarity-ordered so redundant blocks are visible; clicking a cell lists the shared genes (clickable).
- Ridgeline: sets sorted by median log2FC, labels carry median and FDR, clicking a ridge opens its barcode + ECDF.
- Heatmap: clicking a row opens that gene in the Gene Explorer (tick labels visually unchanged).

## 0.13.1-beta (2026-09-21)
- Signature transfer rebuilt for interpretability. It now states its question ("does this gene list separate each dataset's groups better than a random list of the same size?") and answers per dataset with a verdict line: Mann-Whitney AUC (scipy-parity) for magnitude, calibrated by an empirical p against 200 seeded size-matched random signatures. Calibration uses the score's Welch |t|, not AUC: at small n, AUC saturates: 8% of random size-matched signatures (16/200, measured) tie a perfect real signature at AUC 1.0, flooring the AUC-null empirical p at 0.08 so even a perfect positive control cannot reach 0.05; the t null keeps resolution (real t 25.2 vs random max 8.1). Verdicts: replicates / partial / does not replicate. Controls: concordant copy replicates at p = 0.005 (floor); scrambled-design control does not (p = 0.91). Cross-dataset box magnitudes explicitly labeled as non-comparable (within-dataset z-units).

## 0.13.0-beta (2026-09-21)
- Cross-species comparison: mouse and human datasets can now be compared through MGI strict 1:1 orthologs (18,782 pairs from HOM_MouseHumanSequence.rpt, embedded; homology classes with multiple members in either species excluded rather than guessed). A per-dataset species selector controls the matching mode; slot labels and methods text state which mode was used. Case-normalized symbol matching loses 12.8% of true 1:1 orthologs (Trp53/TP53-class renames); validated on a pseudo-human rebuild of the mouse test data: +835 genes recovered, Trp53 matched by ortholog and not by symbol.
- Pairwise comparison readability: FC-FC scatter now colors concordant hits (red/blue by direction) and discordant hits (amber), labels the top 12 concordant genes by combined |t| rank, and is accompanied by a ranked shared-gene table (per-dataset log2FC and FDR, discordant genes listed separately, click-through to gene kinetics).
- Note: a renamed-subset control recovers r = 0.969 rather than 1.000 because dataset B is renormalized within its own (smaller) gene universe: the expected consequence of the never-merge rule, stated here so it is not mistaken for drift.

## 0.12.0-beta (2026-09-21)
Visualization-completeness release: 15 new displays, each validated against a reference before shipping.
- Enrichment: barcode plot (limma tricubeMovingAverage port, worm parity 3.5e-13; full-set ticks), set-vs-rest ECDF with KS distance (scipy parity), set-Jaccard redundancy heatmap (numpy-exact), leading-edge membership matrix (uses gseaES's own leadStart/leadEnd; full sets kept separately as allIdx), per-set log2FC ridgeline.
- DE: p-value Q-Q plot; hit-count stability curves (counts cross-checked against the validated contrast).
- QC/overview: PCA scree + PC1-3 pair panels; library-complexity bars (numpy-exact fractions); depth-saturation curves (seeded binomial thinning, exact endpoint at 100%); per-sample MA vs gene-median; sex-marker sample-identity check (mouse/human marker panels, display-only).
- Patterns: trajectory heatmap for numeric designs (peak-group ordering verified monotone); module volcano and module-score overview heatmap for Discovery.
- Multi-contrast: volcano grid (small multiples, shared thresholds; significant counts identical to the table).
- Compare: rank-rank scatter (Spearman parity with scipy), Bland-Altman of matched log2FCs (bias digit-exact), CAT concordance-at-the-top curves.
- Hubs: scale-free topology diagnostic; a positive slope is labelled as not scale-free.
- Fix found during validation: an initial leading-edge computation double-filtered the already-sliced overlap (mispairing pos with ovIdx); removed in favor of gseaES's own boundary. A saturation shortcut that over-counted detection at low fractions was removed before ship.

## 0.11.1-beta (2026-09-21)
- Comparison-dataset limit raised from 4 to 8 (the cap was UI legibility, not memory or statistics; ~2 MB/dataset).
- Per-gene forest plots across datasets: each dataset's own log2FC with its moderated-t 95% CI (limma topTable confint convention, validated to 1.9e-9 against R); red = significant within that dataset; deliberately no pooled estimate. Reachable from the multi-dataset summary (gene search + per-row links in the consensus table).
- UpSet-style intersection plot for the all-pairs runner: exact membership-pattern counts (disjoint columns), verified against independent recomputation.

## 0.11.0-beta (2026-09-21)
- Branding: benchsiDE (lowercase b) throughout the interface, methods text, reports, and exports.
- Compare tab scales to multiple datasets: up to 4 comparison datasets alongside the session (5 total). Each slot has its own contrast selectors and the same interactive group editor as onboarding (pattern-select and batch assign) for datasets without a design file. Pairwise view (FC scatter, concordance, RRHO) works for any dataset pair, with non-session pairs gene-matched through the session's symbols. New multi-dataset summary: pairwise log2FC-correlation matrix and a consensus gene table (significant within >=2 datasets at the current FDR, vote count only; pooled p-values deliberately not computed). Signature transfer now scores the module in every loaded dataset. Comparison datasets are not stored in session files (raw matrices; re-drop to restore).
- Validated: concordant-copy control r = 1.000; permuted-label null r = 0.107; slot group-editing rebuilds contrasts (null design corrected to the real split recovers r = 1.000); regression and 16-component self-test unchanged.

## 0.10.1-beta (2026-09-21)
- Onboarding: interactive group assignment in the review step for datasets without a design file: per-row checkboxes (click a sample name to toggle), select-all, select-by-name-pattern, and batch "assign selected to group". Reset re-infers from sample names. Per-row inputs remain editable; the design-type preview updates live. Manual grouping produces results identical to the design-file path (verified).

## 0.10.0-beta (2026-09-21)
Feature-completeness release. Every addition validated against its reference implementation before shipping (parity numbers in VALIDATION.md).
- Statistics: genome-wide moderated F (limma parity 1.2e-8); multi-covariate DE models y ~ group + cov1 + cov2 + ... (parity 4e-9, confounded designs refused); FRY self-contained gene-set test (directional p to 8 s.f., mixed p to 6+); limma removeBatchEffect as a display-only toggle (6e-15, statistics provably untouched); edgeR filterByExpr as the recommended default filter for counts (exact kept-set identity; CPM threshold retained); PCA loadings view (8e-16 vs SVD); TPM/FPKM/log input declaration with gene-length column handling; power/design guidance card (RNASeqPower parity <=8e-10).
- QC: MDS at plotMDS parity (1e-13); RLE boxplots; p-value histogram with shape interpretation; voom mean-variance trend plot; sample dendrogram (scipy average-linkage identity); consolidated QC flag card with correlation-outlier detection (verified to catch a deliberately mislabeled sample).
- Multi-contrast: all-pairs DE runner with hit-set overlap matrix; contrast-comparison scatter with concordance statistics.
- Cross-dataset comparison (new tab): second dataset slot with its own filter/TMM/log2 pipeline; symbol-based gene matching; FC/t concordance + directional Fisher overlaps; RRHO maps at exact Bioconductor RRHO parity (9.5e-11, package's N+1 convention reproduced and documented); module signature transfer scored entirely within dataset B. Raw values are never merged across datasets.
- Self-test extended from 11 to 16 components (filterByExpr, MDS, FRY, moderated F, power model).

## 0.9.6-beta (2026-09-21)
- Network hub genes (Co-expression tab): ranks genes by soft connectivity k_i = sum_j a_ij with a_ij = |r|^beta (or signed ((r+1)/2)^beta) over the top-M variable genes, per WGCNA's connectivity definition (Zhang & Horvath 2005; Langfelder & Horvath 2008). Table shows k, normalized k, strongest partner, and the current DE contrast's log2FC/FDR; full ranking exports as CSV; gene names click through to the Gene Explorer. Connectivity values verified against an independent implementation to <4x10^-13 end-to-end; a small-n instability caution is shown with every result and mirrored in the methods text.

## 0.9.5-beta (2026-09-18)
- GSEA running-enrichment-score ("mountain") plots: overlaid curves for the top-N sets (or any sets chosen via per-row "curve" links, up to 8) with per-set hit-position rugs, drawn from the exact running score the ES statistic is computed from (piecewise-linear, 2m+2 vertices per set; curve values verified against an independent implementation to <6x10^-16). SVG export at current size. ORA mode hides the plot.

## 0.9.4-beta (2026-09-18)
- Sample-inclusion bar: group blocks now wrap instead of overflowing the viewport on large sample counts, and each group has an "all" checkbox to include/exclude every member at once (per-sample checkboxes stay in sync).

## 0.9.3-beta (2026-09-18)
- Fixed line-ending handling in all four text loaders (matrix, design, annotation, GMT): bare carriage returns (classic Mac / some Excel exports) are now normalized to newlines instead of deleted. Previously a CR-only design file collapsed to one line and was silently ignored, so grouping fell back to per-sample name inference.
- A design file that matches sample names but yields no factors now raises an explanatory error instead of silently falling back.

## 0.9.2-beta (2026-09-18)
- Built-in annotation now also maps NCBI (Entrez) gene IDs for mouse (MGI_EntrezGene.rpt, 57,006 mappings), human (HGNC, 42,391), and rat (RGD GENES_RAT.txt, 63,039), alongside the existing Ensembl tables. MGI/RGD accessions are carried for deep links; Entrez-keyed genes link directly to their NCBI Gene page.
- Fixed matrix-column classification for numeric gene-ID columns: a column of all-distinct integers with an identifier-style header (or in first position) is now treated as the ID column rather than expression data. Previously an Entrez-keyed matrix gained a phantom sample and lost its IDs. Count columns are unaffected (verified on tied small-integer matrices).

## 0.9.1-beta (2026-09-18)
- Renamed to BenchsiDE (was: RNA-seq Explorer). Legacy session files still load.
- First public release candidate.

## 0.9.0-beta (2026-09-18)
- CAMERA competitive gene-set test for 2-group discovery screens (exact limma parity;
  inter-gene correlation shown per module); permutation test retained for >2 groups.
- In-app numerical self-test: 11 components vs embedded edgeR/limma/scipy/statsmodels
  references, one click, any browser.
- Curated HGNC gene groups as a discovery source (orthology-approximate for mouse/rat).
- Version pinning in header, sessions (with mismatch warning), reports, methods text.
- Discovery results included in HTML report; report DE header aware of voom.
- Fixed: fitFDist infinite-prior branch now matches limma (arithmetic mean of floored
  variances); discovered by the self-test.

## Earlier (development history)
- TMM normalization (edgeR parity), limma moderated t / voom (exact parity incl.
  adaptive lowess span), covariate adjustment, ORA + preranked GSEA, module discovery
  with leave-one-out robustness, k-means heatmaps with silhouette auto-k, publication
  figure exports, sessions/reports, species-aware annotation with deep links,
  deterministic seeding throughout.
