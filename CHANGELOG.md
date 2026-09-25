# Changelog

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
