# VAL_UNT: validation of analyses not covered by earlier checks

Build: benchsiDE 0.22.0-beta (2026-09-30); `benchside/index.html` from the supplied repository archive (unmodified). Inputs: the validation kit and the five NCBI count matrices (sha256 identical to `MANIFEST.json`). The application script was executed in macOS JavaScriptCore with the kit's DOM/Plotly stub (`harness.py`); all HTML control defaults were set from the page before each run. References: scipy 1.17.1, statsmodels 0.14.6, numpy, openpyxl 3.1.5; R with edgeR 4.8.2, limma 3.66.0, cluster, RNASeqPower, babelgene. Every number below was produced in this task.

## Scope

This report covers the items in the earlier internal audit's 'Not assessed' list and the 16 inventory rows marked 'reviewed, not separately tested', plus the analyses listed in the task: Patterns and auto-k, single-gene tests, co-expression, hubs and network edges, QC plots, power calculator, input-type detection, design parsing, removeBatchEffect-style display, Venn/UpSet, Excel/Prism export, session round trip, ORA and Discovery size limits, hidden-module threshold, consensus rule and mouse-human ortholog mapping. 140 checks were run: 114 pass and 26 fail (`VAL_UNT_checks.csv`). Tolerances: 1e-8 absolute for log-scale statistics, 1e-6 relative for p-values, 0 for counts, sets and calls.

Coverage of the 16 'reviewed, not separately tested' rows: pairwise brackets; co-expression minimum n; auto k; top-1,000 genes; heatmap top 50; ORA universe; GSEA ranking statistic; hub network (2,000 genes, beta 6, unsigned); network edges (|r| >= 0.7, 3 per gene, 25 hubs + 5 neighbours); Compare concordance and overlaps; Compare minimum N; meta inclusion; signature-transfer minimum genes; low-depth flag; p-value histogram notes; combined factors. All were tested (checks below).

One reference setting was corrected during the work: the first limma reference used `keep.lib.sizes = TRUE`, which gave TMM factors differing by 1.0e-3 and log2FC by 4.7e-5. The app recomputes library sizes on kept genes (the edgeR user's guide setting, `keep.lib.sizes = FALSE`); with that setting log2FC agree to 3.9e-14 and FDR to 3.0e-11. This was a reference error, not an app defect.

## Results by feature

### Gene Explorer single-gene tests
All genes, GSE63310 (3 cell types, 16,624 genes) and GSE121212 Group (6 groups, 144 samples, 23,168 genes): ANOVA F and p agree with `scipy.stats.f_oneway` (F <= 1.0e-12, p <= 8.2e-12 relative). Pairwise Welch p and Holm-adjusted p agree with scipy and statsmodels in both modes (<= 5.4e-10). Brackets, stars, the pairwise table, bar means/SEM and points matched exactly for 60 genes per dataset. Spearman trend (GSE186063, Age as numeric grouping, 21,112 genes): rho exact, p <= 8.6e-12. Kruskal-Wallis H is exact, but small p-values are wrong: 4,893 genes of GSE121212 exceed the tolerance and 757 are reported as p = 0. Genes constant within both groups get Welch p = 1 where scipy returns NaN (1 and 16 genes on GSE63310; excluded from the comparison).

### Patterns
On GSE63310 and GSE121212, k-means started from the app's own k-means++ centres (captured by instrumenting a copy of `kmeans`) gives the same clusters as `stats::kmeans(algorithm = "Lloyd")` for k = 2..10; mean silhouettes agree with `cluster::silhouette` to 7.1e-11 and 2.1e-10 (the app stores distances as Float32). The app's 4-start optimum is up to 1.6% (GSE121212, k = 9) above the best of 50 random starts. Planted-k benchmark (`val/planted_k_benchmark.csv`):

| groups | planted_k | noise_sd | seed | auto_k | silhouette | ARI |
|---|---|---|---|---|---|---|
| 6.0 | 3.0 | 0.25 | 1.0 | 3.0 | 0.941 | 1.0 |
| 6.0 | 3.0 | 0.25 | 2.0 | 3.0 | 0.931 | 1.0 |
| 6.0 | 3.0 | 0.5 | 1.0 | 3.0 | 0.881 | 1.0 |
| 6.0 | 3.0 | 0.5 | 2.0 | 3.0 | 0.863 | 1.0 |
| 6.0 | 3.0 | 1.0 | 1.0 | 3.0 | 0.763 | 1.0 |
| 6.0 | 3.0 | 1.0 | 2.0 | 3.0 | 0.729 | 1.0 |
| 6.0 | 4.0 | 0.25 | 1.0 | 4.0 | 0.935 | 1.0 |
| 6.0 | 4.0 | 0.25 | 2.0 | 4.0 | 0.918 | 1.0 |
| 6.0 | 4.0 | 0.5 | 1.0 | 4.0 | 0.87 | 1.0 |
| 6.0 | 4.0 | 0.5 | 2.0 | 4.0 | 0.837 | 1.0 |
| 6.0 | 4.0 | 1.0 | 1.0 | 4.0 | 0.736 | 1.0 |
| 6.0 | 4.0 | 1.0 | 2.0 | 4.0 | 0.678 | 1.0 |
| 6.0 | 5.0 | 0.25 | 1.0 | 5.0 | 0.914 | 1.0 |
| 6.0 | 5.0 | 0.25 | 2.0 | 5.0 | 0.922 | 1.0 |
| 6.0 | 5.0 | 0.5 | 1.0 | 5.0 | 0.829 | 1.0 |
| 6.0 | 5.0 | 0.5 | 2.0 | 5.0 | 0.845 | 1.0 |
| 6.0 | 5.0 | 1.0 | 1.0 | 4.0 | 0.666 | 0.783 |
| 6.0 | 5.0 | 1.0 | 2.0 | 5.0 | 0.694 | 1.0 |
| 6.0 | 6.0 | 0.25 | 1.0 | 6.0 | 0.918 | 1.0 |
| 6.0 | 6.0 | 0.25 | 2.0 | 6.0 | 0.921 | 1.0 |
| 6.0 | 6.0 | 0.5 | 1.0 | 6.0 | 0.835 | 1.0 |
| 6.0 | 6.0 | 0.5 | 2.0 | 6.0 | 0.842 | 1.0 |
| 6.0 | 6.0 | 1.0 | 1.0 | 5.0 | 0.682 | 0.829 |
| 6.0 | 6.0 | 1.0 | 2.0 | 6.0 | 0.684 | 1.0 |
| 6.0 | 8.0 | 0.25 | 1.0 | 8.0 | 0.91 | 1.0 |
| 6.0 | 8.0 | 0.25 | 2.0 | 8.0 | 0.923 | 1.0 |
| 6.0 | 8.0 | 0.5 | 1.0 | 8.0 | 0.819 | 1.0 |
| 6.0 | 8.0 | 0.5 | 2.0 | 8.0 | 0.845 | 1.0 |
| 6.0 | 8.0 | 1.0 | 1.0 | 7.0 | 0.653 | 0.867 |
| 6.0 | 8.0 | 1.0 | 2.0 | 7.0 | 0.642 | 0.864 |
| 3.0 | 3.0 | 0.25 | 7.0 | 2.0 | 0.694 | 0.524 |
| 3.0 | 3.0 | 0.5 | 7.0 | 2.0 | 0.602 | 0.403 |
| 3.0 | 4.0 | 0.25 | 7.0 | 2.0 | 0.891 | 0.365 |
| 3.0 | 4.0 | 0.5 | 7.0 | 2.0 | 0.865 | 0.365 |
| 6.0 | 1.0 | 0.15 | 11.0 | 9.0 | 0.152 | nan |
| 6.0 | 1.0 | 0.15 | 12.0 | 10.0 | 0.155 | nan |

Six groups: exact at noise sd <= 0.5; at sd 1.0, 4 of 10 under-called by one. Three groups: auto k = 2 in every case. Null data: silhouette 0.152 and 0.155 with the weak-structure label shown. The silhouette of singleton clusters is 1 instead of 0. The page default is a fixed k = 6; auto is an option.

### Heatmap
Top-50 variable genes and z-scores exact; every cluster of scipy's average-linkage tree (1 - Pearson r) is contiguous in the app's row order (49 of 49).

### Co-expression, hubs, network
Three query genes against all genes (49,869 pairs): r <= 6.7e-16, p and BH FDR <= 8.9e-12 relative, slope <= 4.4e-15, top-30 tables identical, scatter regression, CI band and Spearman rho <= 7.1e-15. Fewer than 4 samples are refused. Soft connectivity (2,000 genes, unsigned beta 6 and signed beta 12) agrees with numpy to 6.6e-15 relative, with identical ranking and strongest partners. The network (25 hubs + 5 neighbours, |r| >= 0.7, 3 strongest per gene) has the same 65 nodes and 136 edges as the numpy reconstruction.

### QC plots
Library sizes, RLE quantiles and complexity shares are exact; density curves agree with an exact Gaussian KDE (nrd0, 512 points, cut 3) to 2.9e-10; saturation points lie within 2.43 SD of the binomial-thinning expectation. QC flags equal a reimplementation of the three rules on GSE63310 (none) and GSE121212 (low depth GSM3427959; correlation outliers GSM3427893, GSM3427900, GSM3427959) and on a crafted low-depth matrix. The median for the low-depth rule is the upper middle value for even n (19,847,506 vs 19,803,126.5 on GSE121212). Sex check on GSE186063 (XIST vs 5 chrY genes): plotted values exact, 65/65 samples on the side matching the Sex column. Detected genes and saturation count filtered genes only.

### Power calculator
`rnapowerJS` agrees with `RNASeqPower::rnapower` over 1,664 inputs (sample size 1.2e-9, detectable fold change 1.3e-9, power 7.4e-9 relative). Pilot inputs on GSE54456 reproduce the documented rule exactly (CV 0.4199 / 0.3574, depth 364.84), but the CV definition overstates sample sizes:

| cv_definition | cv_Psoriasis | cv_normal | n80_1.25x | n80_1.5x | n80_2x |
|---|---|---|---|---|---|
| app (raw-count CV) | 0.4199 | 0.3574 | 49 | 15 | 6 |
| library-normalized CV | 0.3389 | 0.279 | 32 | 10 | 4 |
| normalized, Poisson removed | 0.3233 | 0.2641 | 29 | 9 | 3 |
| edgeR common BCV | 0.3819 | 0.3398 | 43 | 13 | 5 |

### Input-type and design handling
Tab, comma, semicolon, CRLF and BOM versions of the count matrix give identical log-CPM; Ensembl IDs with version suffixes map like Entrez IDs (16,219 of 16,624 kept genes). log2(CPM + 1), log2 CPM with negative values, CPM and TPM are classified correctly. The six kit matrices are all classified as counts.

### Batch-adjusted display
`dispM` equals `limma::removeBatchEffect(x, batch, design = model.matrix(~group))`: GSE186063 with Sex 5.3e-14; GSE83645 with Patient 5.7e-14 and Site 6.0e-14.

### Venn, UpSet, all-pairs
GSE121212 Group, 15 contrasts: kept genes identical to `filterByExpr`; up/down counts and 34,517 hit-set memberships identical to limma-trend run on each pair of groups (`contrast_sets_reference.R`); all 21 Venn regions (3 contrasts; any/up/down) and the 14 UpSet bars and dot matrix identical to set algebra on the reference sets.

### Excel/Prism export
For p_gene, p_lib, p_pat, p_heat, p_upset, p_venn, p_coscatter, p_pca and p_saturation the generated .xlsx was opened with openpyxl: every 'All values' entry and the Prism-layout sheet equal the plotted values. RLE export fails.

### Session save/restore
With default filtering, save then restore into a fresh load reproduces DE (voom, one sample excluded, custom group order, FDR 0.1, |log2FC| 0.5), patterns, heatmap order, ORA and gene tests exactly. Non-default filtering is not reproduced.

### Enrichment, Discovery and consensus sweeps
ORA: 47,534 set tests (Hallmark and GO BP, up and down, four size settings) agree with `scipy.stats.hypergeom` and BH (<= 1.6e-10 relative). Significant sets at FDR <= 0.05 on GSE54456 by size limits (`val/ora_size_sweep.csv`): Hallmark up 20 and down 3 for any max >= 200 (3 and 0 at max 100), insensitive to min; GO BP up 265 to 739 and down 23 to 189 over min 1-25 and max 100 to unlimited (default 5-2,000: 739 up, 142 down).

Discovery (GSE54456, families + HGNC groups + Hallmark, CAMERA path; `val/discovery_size_sweep.csv`):

| min_size | max_size | modules_screened | families | hgnc | hallmark | sig_fdr05 | nominal_p01 | camera | hidden_lt_0.1 | hidden_lt_0.25 | hidden_lt_0.5 | jaccard_vs_default |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 100 | 2933 | 1388 | 1528 | 17 | 0 | 47 | True | 0 | 0 | 0 | 1.0 |
| 5 | 100 | 1503 | 596 | 890 | 17 | 0 | 28 | True | 0 | 0 | 0 | 0.596 |
| 3 | 50 | 2887 | 1380 | 1500 | 7 | 0 | 45 | True | 0 | 0 | 0 | 0.957 |
| 3 | 200 | 2979 | 1390 | 1539 | 50 | 0 | 51 | True | 0 | 0 | 0 | 0.922 |
| 10 | 500 | 676 | 207 | 419 | 50 | 0 | 14 | True | 0 | 0 | 0 | 0.196 |
| 3 | 1000 | 2981 | 1392 | 1539 | 50 | 0 | 51 | True | 0 | 0 | 0 | 0.922 |

No module passes FDR at any size setting; the nominal list changes with min size (Jaccard 0.596 for min 5, 0.196 for 10-500). On null splits of the 81 normal-skin samples into 3 or 4 random groups (permutation path), 0 of about 2,960 modules passed FDR in 6 of 6 splits, and 4.2-6.9% had p <= 0.05. Hidden-module counts at thresholds 0.1/0.25/0.5 (`val/discovery_hidden_sweep.csv`): 0/0/0 on GSE54456, GSE186063 and GSE83645; 0/0/2-3 on GSE41745. The flag did not fire at the default threshold on any kit dataset.

Consensus (5 datasets, 16,150 genes matched in all; `val/consensus_sweep.csv`): the app's 12,425 consensus genes and concordant/mixed flags equal a numpy vote count. Selected settings at FDR 0.05: >= 2 datasets 12,425 (77%); >= 3 8,522; >= 4 3,530; all 5 693; >= 2 with |log2FC| >= 1 2,090. Of the 20 classic psoriasis genes, 16 are matched in all datasets and all 16 qualify at >= 2; CCL20, CXCL8, IL17A and IL23A are excluded by GSE41745.

### Ortholog mapping Details: `val/ortholog_comparison.csv`.

## Feature status

| feature | status | basis |
|---|---|---|
| Gene Explorer: one-way ANOVA | validated | All genes of GSE63310 (16,624) and GSE121212 Group (23,168) vs scipy f_oneway; F <= 1.0e-12, p <= 8.2e-12 relative; display strings checked for 60 genes each |
| Gene Explorer: Kruskal-Wallis | defect_found | H exact; p loses precision below ~1e-10 and returns 0 below ~1e-16 |
| Gene Explorer: Spearman trend | validated | GSE186063 with Age as numeric grouping, 21,112 genes vs scipy spearmanr; rho exact, p <= 8.6e-12 relative |
| Gene Explorer: pairwise brackets (Welch, Holm, stars) | validated | Both modes, all genes, two datasets vs scipy/statsmodels (<= 5.4e-10 relative); brackets, stars and tables for 120 genes. Statistics under batch display mislabelled |
| Gene Explorer with batch-adjusted display | defect_found | Tests use adjusted data although labelled unadjusted |
| Patterns: k-means and silhouette | validated_with_caveat | Same k-means++ starts give identical clusters to stats::kmeans (Lloyd); silhouette matches cluster::silhouette to 2.1e-10 without singletons; singleton rule differs |
| Patterns: auto-k | validated_with_caveat | Planted k recovered in 20/20 six-group scenarios (noise <= 0.5) and within 1 at noise 1.0; fails for three groups. Default selection is k = 6, not auto |
| Patterns: gene selection and profiles | validated | Top-1,000 variable genes and z-profiles exact (GSE63310, GSE121212) |
| Heatmap: gene selection and hierarchical clustering | validated | Top-50 genes and z-scores exact; every scipy average-linkage cluster contiguous in the app row order (GSE63310) |
| Co-expression (r, slope, p, BH FDR, scatter fit) | validated | 3 query genes x 16,623 genes vs scipy/statsmodels; r <= 6.7e-16, p and FDR <= 8.9e-12 relative; n >= 4 rule checked |
| Hubs: soft connectivity and partners | validated | Unsigned beta 6 and signed beta 12, 2,000 genes vs numpy (<= 6.6e-15 relative), identical ranking and partners |
| Network view edges | validated | 65 nodes, 136 edges identical to numpy reconstruction of the stated rule |
| QC: library size, RLE, density, complexity, flags, sex check | validated | Plotted values equal independent computations (GSE63310, GSE121212, GSE186063); sex markers concordant with metadata in 65/65 samples |
| QC: detected genes and saturation | validated_with_caveat | Values exact for the implemented definition, which counts filtered genes only |
| Power calculator: rnapower | validated | 1,664 grid inputs vs RNASeqPower::rnapower, <= 7.4e-9 relative |
| Power calculator: pilot inputs | defect_found | Computed as documented (exact), but raw-count CV overstates n; groups not the DE contrast |
| Input-type detection and parsing | defect_found | Separators, CRLF, BOM, Ensembl version suffixes handled; log2(x+1) filter, FPKM, sorted matrices, NA/non-numeric cells, duplicated IDs fail |
| Design-file parsing | defect_found | Orientation, extra samples, GEO prefixes, numeric factors, combined-factor rule correct; missing samples become group 'NA'; quoted files refused |
| Batch-adjusted display (removeBatchEffect) | validated | GSE186063 (Sex) and GSE83645 (Patient, Site) vs limma::removeBatchEffect, <= 6.0e-14 |
| Venn / UpSet / all-pairs (DE tab) | validated_with_caveat | 15 contrasts on GSE121212: hit sets, Venn regions and UpSet counts identical to limma-trend reference; top-14 truncation undisclosed |
| Excel/Prism export | defect_found | 9 plots: exported values equal plotted values (openpyxl); RLE export throws; UpSet labels |
| Session save/restore | defect_found | Exact round trip with default filtering; non-default filter/normalization not reproduced; several settings not stored |
| ORA (hypergeometric, size limits, universe) | validated_with_caveat | 47,534 set tests vs scipy hypergeom (<= 1.6e-10 relative); size sweep reported; cluster-ORA universe concern |
| GSEA ranking statistic | validated | Observed-label ranking vector equals the DE tab's moderated t to 1.2e-11 (3 contrasts) |
| Discovery: size limits, hidden-module threshold, >2-group null | validated_with_caveat | 0 modules at FDR <= 0.05 on 6 null 3-/4-group splits; size and hidden-threshold sweeps reported; the hidden-module flag never fired on any kit dataset |
| Compare: concordance, overlaps, minimum N | validated | 4 dataset pairs: N, r of log2FC and t, overlaps and hypergeometric p identical to scipy; N >= 50 rule checked |
| Consensus genes | validated_with_caveat | Vote count reproduced exactly (12,425 genes); rule restricts to genes matched in all datasets |
| Meta inclusion (minK, finite SE) | validated | Rule reproduced on synthetic inputs |
| Signature transfer minimum genes | validated | < 3 genes refused |
| Mouse-human ortholog mapping | defect_found | 94.1% of pairs confirmed by HCOP; not 1:1 as documented; mitochondrial genes unmatched |

## Not assessed

- Behaviour in real browsers (all tests ran in JavaScriptCore with a DOM/Plotly stub; rendering, file dialogs and downloads were not exercised)
- IMQ mouse data: not among the supplied inputs; ortholog mapping was assessed on the full app table, GSE63310 (mouse) and GSE54456 (human)
- Heatmap k-means row ordering (uses the same autoK/kmeans functions validated in Patterns, but its own output was not compared)
- Heatmap cluster ORA (not run)
- Hub, discovery, consensus, Venn and pattern CSV exports (only the Excel/Prism export path was opened and checked)
- QC p-value histogram rule: tested on four synthetic p-value vectors for category only; bin counts not compared
- Sex check on a single-sex dataset (GSE63310: chrY markers are removed by the filter and no plot is drawn; behaviour noted, not validated)
- Signature-transfer AUC/verdict, RRHO, meta-analysis pooling, GSEA permutation p/FDR, FRY/CAMERA, voom: covered by other validation tracks, not re-tested here
- Session restore of comparison datasets (not stored in the session file; not tested further)
- Hub/network sensitivity to beta, M and |r| thresholds (only the defaults and one signed setting were checked)

## Files

- `VAL_UNT_checks.csv`: every check (id, name, dataset, reference, n compared, maximum deviation, tolerance, pass, note).
- `val_unt_runner.py`: CI runner; reproduces the core checks from the kit files (one run: 42 checks, 38 pass).
- `validation/ref/patterns_reference.R`, `contrast_sets_reference.R`, `removebatch_reference.R`, `power_reference.R`, `ortholog_reference.R`: R references called by the runner.
- Sweep and benchmark tables: `planted_k_benchmark.csv`, `power_cv_definitions.csv`, `ora_size_sweep.csv`, `discovery_size_sweep.csv`, `discovery_hidden_sweep.csv`, `discovery_null_multigroup.csv`, `consensus_sweep.csv`, `ortholog_comparison.csv`.
