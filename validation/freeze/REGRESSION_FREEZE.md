# benchsiDE 1.0.0-rc.1: engine regression against the 0.23.0 baseline

Build under test: `benchside/index.html`, APP_VERSION `1.0.0-rc.1 (2026-10-03)`, sha256 of the application script `cadbda80bd67db0e3cb15fcbcf138c23fb75cb24ff17cf6c28a79c68886c7d5b` (matches the frozen value; check BLD-01). The build was not modified. Baseline: `baseline_0.23.0.json.gz` written by `tests/engine/snapshot.py` on `0.23.0-beta (2026-10-01)`.

All application numbers were produced by running the application script in macOS JavaScriptCore (`jsc`) with the repository's stub page and stub Plotly (`tests/engine/harness.py`). This is not a browser: rendering, real font metrics, image pixels and Content-Security-Policy enforcement are not tested here. References: R 4.5.3 with metafor 5.0.1 (per-gene `rma`), limma 3.66.0 results stored in `validation/freeze/interaction_reference_GSE143688.json.gz`, scipy, and the earlier validation documents. All count matrices matched their checksums in `datasets.json` (DAT-01 to DAT-08).

**Result: 211 of 227 checks pass.** The 16 failing rows correspond to the 7 defects listed below. None of them is in a statistical computation. The checks are in `REGRESSION_FREEZE_checks.csv`.

## 1. Engine checks and snapshot comparison

`tests/engine/run_checks.py` passed all five checks on the frozen build: core PASS; interaction PASS; cross_dataset PASS; figures_and_qc PASS; discovery_text PASS. External validation 59/59, self-test 21/21, interaction against limma max |ΔlogFC| 5.00e-12, default meta-analysis signature 1773 genes with 20/20 canonical genes (ENG-01 to ENG-06).

`tests/engine/snapshot.py` was run unchanged and every field compared with the baseline (47 fields: version string; per dataset n, up, down, log2FC vector and q vector for 5 datasets; meta method, counts, signature, gene list and the mu, se, p, q, k vectors; homology-family keys and d, p, q, k; meta-signature enrichment sets and m, k, p, q). 46 fields are identical to the last bit. The only difference is the version string (0.23.0-beta (2026-10-01) → 1.0.0-rc.1 (2026-10-03)), explained by the release (SNP-01).

The snapshot script sets the model to DerSimonian-Laird but does not set the test. In the stub the test selector returns an empty value (z), because the stub ignores `selected`. That is why the script reproduces the 0.23.0 analysis. In a browser the same script would run DerSimonian-Laird with Hartung-Knapp, giving 1778 signature genes instead of 2151.

To make the comparison independent of stub defaults, every page control was then set to its browser default (parsed from the HTML). The meta-analysis was run four times on the same four human cohorts as the snapshot:

| Model / test | genes tested | signature genes | fields differing from baseline |
|---|---|---|---|
| DL + z | 19401 | 2151 | none |
| REML + hk | 19401 | 1773 | meta.method, meta.sig, meta.up, meta.sigGenes, meta.rows.mu, meta.rows.se, meta.rows.p, meta.rows.q, enr.sets, enr.k, enr.p, enr.q |
| DL + hk | 19401 | 1778 | meta.sig, meta.up, meta.sigGenes, meta.rows.p, meta.rows.q, enr.sets, enr.k, enr.p, enr.q |
| REML + z | 19401 | 2157 | meta.method, meta.sig, meta.up, meta.sigGenes, meta.rows.mu, meta.rows.se, meta.rows.p, meta.rows.q |

- With the method set back to DerSimonian-Laird and the z test, the 0.23.0 results are reproduced exactly: every DE, meta, family and enrichment field is identical, max deviation 0 (VDZ rows).
- With the new default (REML + Hartung-Knapp), the DE and homology-family fields are unchanged. Only the meta fields and the enrichment of the meta signature differ. The signature changes from 2151 to 1773 genes, the figures stated in the 1.0.0-rc.1 CHANGELOG (VRH rows). The enrichment difference follows from the signature: enrichment is identical under DL + z.
- The new default was refitted gene by gene in metafor (19400 genes; `rma(method='REML', test='knha')`):
  - pooled log2FC agrees to 9.60e-11, the Hartung-Knapp standard error to 5.57e-11, p to 1.27e-09 (relative), tau² to 1.06e-09 and the CI to 1.91e-10.
  - The prediction interval agrees with metafor's `predtype='Riley'` interval (t on k−2 df) to 4.83e-09. It does not agree with metafor's default interval (t on k−1 df), so the claim in META_REML_HK.md refers to the Riley form.
  - One gene (TMEM126B) was not fitted by metafor's default Fisher scoring. For that gene the app's tau² is at the maximum of the restricted likelihood (MFR-TMEM).
- DerSimonian-Laird + z agrees with metafor to 7.99e-15 in the pooled estimate (MFR rows).
- With a minimum of two datasets, the counts of genes at meta FDR ≤ 0.05 equal those in META_REML_HK.md for REML + HK, REML + z and DL + HK (DOC rows).

A session saved before 1.0 (no meta settings, version 0.23.0) restores DerSimonian-Laird + z, shows the version warning, and reproduces the baseline meta results exactly (SES-PRE10). **Every difference from the baseline is explained by a change listed in the 1.0.0-rc.1 CHANGELOG. No difference is unexplained.**

## 2. Tabs, exports, report, methods and sessions

Four configurations were driven through the interface functions, with every control set to its browser default before the scenario values were applied:
- the four psoriasis cohorts, with the mouse dataset GSE143688 as a fourth comparison dataset (PSO_SETUP);
- GSE143688_all with the Cell factor, Day as a DE covariate, and the interaction card ((Aldara_KO − Control_KO) − (Aldara_WT − Control_WT), Day covariate, all eight groups);
- GSE171012 with Status on the 66 whole-skin samples;
- GSE171012 with Status on the 59 CD8 samples, for the Sample traits tab (see Deviations).

In each configuration every tab was opened twice. The following were then run:
- Gene Explorer: bar, box and violin plots; the multi-gene comparison in four chart types; gene set; gene in comparison datasets.
- DE: volcano selection, moderated F, power, PCA loadings, all pairs, volcano grid, contrast comparison and Venn diagram.
- Co-expression: hubs, neighbourhood and network.
- Heatmap: hierarchical and k-means.
- Patterns.
- Enrichment: ORA, GSEA with its running-score plot, FRY and the enrichment matrix.
- Discovery.
- Sample traits: association, CAMERA and the module × trait matrix.
- Compare tab: pair comparison, homology families, family dot plot, multi-dataset summary, meta-analysis, single-gene and multi-gene forest plots, signature transfer, enrichment of the meta signature, and the hub across datasets.

Every export button was then pressed. Every plot div was exported through the toolbar's image button (preview window, then Save) and its Excel/Prism button. The report, the summary figure, the module heatmap, the gene panels, the methods text (generated and copied) and the session file were also produced.

| Configuration | actions | plot divs | divs with data | Excel workbooks | image files | CSV/GMT files | runtime errors |
|---|---|---|---|---|---|---|---|
| PSO 4 cohorts + GSE143688 (comparison) | 106 | 56 | 53 | 56 | 60 | 23 | 0 |
| GSE143688_all (Cell; interaction) | 96 | 42 | 41 | 42 | 46 | 21 | 0 |
| GSE171012 (Status; whole skin) | 92 | 42 | 41 | 42 | 46 | 20 | 0 |
| GSE171012 (Status; CD8 samples, SortedCells trait) | 97 | 45 | 44 | 45 | 49 | 20 | 0 |

- No runtime error occurred in any tab, action or export (TAB, ACT rows).
- Every exported workbook was opened with Python's zipfile. Every plot div that has data produced a workbook with numeric cells (XLS rows).
- The plot divs without data are placeholders: no biotype column; patterns with fewer than three groups; UpSet with one contrast. They still offer both export buttons, and their workbook contains no values.
- Among CSV exports that can be reached in the interface, `trait_gene_sets_modules.csv` is written with a header only when no module or CAMERA result exists. The hub-across-datasets CSV was header-only in single-dataset runs, but its card cannot be opened without comparison datasets (the Across datasets button is hidden), so it is not counted.
- `consensus_genes.csv` writes the literal `NaN` for datasets in which a gene is not matched.
- The report contains the methods text, has no `undefined`/`NaN` tokens and has no external hyperlinks. The methods text describes the analyses run: REML + Hartung-Knapp for the meta-analysis, the interaction contrast, and the trait model (REP, MET rows).

Session round trip. The session file from each configuration was loaded in a fresh engine run with the same data, after every restorable setting had been set to a different value. The analyses were then re-run from the restored settings (DE, interaction, meta-analysis, trait association, ORA). In all 56 comparisons the results were identical, max deviation 0.0. These cover the log2FC/p/q vectors, the interaction, the meta rows including CIs and tau², the trait t/p/q, the enrichment, the volcano selection, the gene, the comparison list, the included samples and the groups. A session saved again after the restore equals the original (SES rows). The enrichment query source is not stored in sessions; the first source (DE up) was selected in both runs.

## 3. Figure content

- Hyperlinks: no figure text or layout contains a hyperlink or URL, and no caption file contains user-interface instructions or links (FLK, FCP).
- UI instructions: two figures carry instructions in their titles, which are part of the exported image and workbook: the set-overlap matrix ('… (click a cell for the shared genes)') and the ridge plot ('… (click a ridge for its barcode)'). This occurs in all four configurations.
- Gene names: individual-gene exports (Gene Explorer, gene in comparison datasets, co-expression pair, single-gene forest, trait gene) show only the decoded symbol in title and file name. A gene stored as `Pakap (ENSMUSG00000038729)` is exported as `Pakap` / `Pakap.svg` (GEN rows).
- Axes: all bar plots and the three forest plots include 0 on the value axis (ZER rows). The single-gene and hub forest plots are symmetric about 0, as the CHANGELOG states for forest plots. The multi-gene forest plot is not symmetric, though it includes 0.
- Angled labels: with the figures_and_qc method, labels fit in the Gene Explorer and comparison plots at 45°, −45° and auto, at 620, 1000 and 1460 px (ASW rows). They also fit in every final plot except the enrichment-across-contrasts matrix of GSE171012, which has 10 contrasts with labels of up to 75 characters. There the right margin reaches its cap of 45 % of the width (279 px) and the last label overruns the edge:
  - 77 px with the stub metric;
  - +19.5 px with the app's own fallback metric of 0.58 em per character;
  - -22.4 px at 0.5 em.

  Whether this happens in a browser depends on the font, so the finding is provisional.

## 4. Content-Security-Policy and subresource integrity

- Both builds carry exactly one CSP meta tag, placed before the first script, with `default-src 'none'` and `connect-src 'none'`.
- The online build's only external source is `https://cdn.plot.ly` in script-src, and its only external tag is the pinned Plotly 2.32.0 script.
- The offline build (present in `dist/`, not rebuilt) has no http(s) source in tags, CSS or CSP.
- The integrity attribute `sha384-7TVmlZWH60iKX5Uk7lSvQhjtcgw2tkFjuwLcXoRSR4zXTyWFJRm9aPAguMh7CIra` equals the sha384 of the file downloaded from cdn.plot.ly. The Plotly code inlined in the offline build has the same sha384 and the sha256 pinned in `build_offline.py`.
- The application script is byte-identical in the two builds (CSP, SRI, BLD-02).
- Browser enforcement of the policy was not tested (see Not assessed).

## Not assessed

- Real-browser rendering, layout and image pixels: Plotly.toImage is a stub, so exported SVG/PNG files and the composite figures (gene panels, module heatmap, summary figure) were checked only for being written and for their layout specification.
- Content-Security-Policy enforcement and refusal of outside requests in Chromium, Firefox and WebKit: this needs a browser; CI covers it.
- Angled-label fit under real fonts: text widths come from the stub (7 px per character); label-fit results are provisional.
- Timing benchmark (Validate statistics → Timing benchmark; brun/bjson/bcsv/bdata): not run.
- File input by drag-and-drop or file picker (FileReader path): data were passed to loadFromText directly; load timing and its clipboard copy were not tested.
- Unhandled promise rejections inside asynchronous handlers that do not return their promise cannot be intercepted in jsc. Errors were caught when they were thrown synchronously or by awaited promises.
- Sample traits on whole-skin GSE171012: the design has no numeric column, so the tab only shows its notice. Association, CAMERA and the module × trait matrix were exercised on the CD8 samples with a GEO annotation (see Deviations).
- Not applicable in some configurations, and not judged there: Venn/UpSet with one contrast (PSO session); homology families between two datasets of the same species; gene in comparison datasets without comparison datasets.

## Deviations

- The stub was extended for the scenario runs (run_checks.py and snapshot.py were run with the unmodified harness). The extensions are: browser-default values for every static control; select values updated when options are rewritten; querySelectorAll over rendered markup; canvas and Image stubs; download capture; and document.removeEventListener, without which the figure-preview Save failed in the stub.
- Sample traits: a SortedCells row (GEO GSE171012 series matrix, 'sorted cells', real values for all 59 CD8 samples, NA for whole skin) was added to the design. Traits were run on the 59 CD8 samples (Status groups, log2(x+1), Subject block) as a fourth configuration.
- Comparisons with the baseline used the four human cohorts, as snapshot.py does. In the full-interface PSO run the meta-analysis used the default selection of all five datasets.
- Session round trip: the enrichment query source is not stored in sessions and was set to the first option in both runs. For single-dataset configurations meta.mink was excluded from the re-save comparison, because it was empty in the saved session.
- Full-interface runs use 200 Discovery permutations (default 1000).

## Reproduction

Engine runs used at most 3 parallel jsc processes. Inputs: repository archive v371d5f9a, engine data vc0cbec56, validation kit v213ee441, and the GEO GSE171012 series matrix downloaded from ftp.ncbi.nlm.nih.gov. Plotly 2.32.0 was downloaded from cdn.plot.ly for the integrity check. R references: R 4.5.3, metafor 5.0.1.
