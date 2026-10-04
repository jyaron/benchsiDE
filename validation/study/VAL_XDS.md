# Validation of the Compare tab (cross-dataset analyses), benchsiDE 0.22.0-beta

Build under test: `benchside/index.html`, APP_VERSION "0.22.0-beta (2026-09-30)", SHA-256 17511676c3096bda2f4ba964dea8ce5b78faa300a29ad0d1007c210da7c7ae65 (equal to kit/MANIFEST.json). The build was not modified. All five count matrices matched their MANIFEST SHA-256.

Engine: macOS JavaScriptCore through `kit/harness.py` (DOM/Plotly stub). Matrix text was read inside JavaScript with the JSC `readFile()` builtin instead of being embedded with `json.dumps`; the text passed to the app is the same. References: R 4.5.3, edgeR 4.8.2, limma 3.66.0, metafor 5.0.1, RRHO 1.50.0; Python 3.11 with scipy and statsmodels for independent code. Ortholog references: MGI `HOM_MouseHumanSequence.rpt` and NCBI `gene_info` (human, mouse), both downloaded on 2026-09-30.

Machine-readable results: `VAL_XDS_checks.csv` (343 checks, 334 pass, 9 fail). Each row gives the comparison, the reference, the number of values compared, the maximum deviation and the tolerance. Supporting tables: `VAL_XDS_concordance.csv`, `VAL_XDS_classic_genes.csv`, `VAL_XDS_loo.csv`, `VAL_XDS_signature_transfer.csv`, `VAL_XDS_hallmark_ora.csv`, `VAL_XDS_mad3_top50.csv`.

Tolerances follow the task: 1e-8 absolute for continuous statistics on the log scale, 1e-6 relative for p-values, 0 for counts, sets and calls. Values that the app shows only as rounded text (r to 3 decimals, AUC to 2, p as `toExponential(2)`) were compared at half a unit of the last displayed digit. Where the full-precision value is available internally, it was also compared at the strict tolerance.

## 1. Design of the test

Session: GSE54456, psoriatic vs normal skin (Cond; 90 vs 81). Comparison datasets were added through `onAddDataset()` with `CMPB_TXT/CMPB_DES/CMPB_NAME`, re-grouped with `setDsFactor()` where needed, and paired through `d.covSel`:

| Study | Contrast (B vs A) | Grouping | Pairing factor | Genes kept | Matched to session |
|---|---|---|---|---|---|
| GSE54456 (session) | Psoriasis_skin vs normal_skin | Cond | none available | 19,518 | n/a |
| GSE186063 | lesion vs non-lesion | Type | Pair (inferred; 2 unpaired non-lesion samples dropped) | 21,112 | 19,308 |
| GSE121212 | PSO_lesional vs PSO_non_lesional | Group | Patient | 23,168 | 19,161 |
| GSE83645 | psoriasis vs uninvolved | Cond | Patient | 20,956 | 19,440 |
| GSE41745 | lesional vs non_lesional (2 vs 3) | Cond | Patient | 18,036 | 16,294 |

Two complete runs were made: an unpaired run, and a run with every available pairing factor. In each run the meta-analysis was computed four times: random effects (DL) and fixed effect (FE), each with the gene required in all 5 studies (the default, minK = 5) and in at least 2 (minK = 2). Additional runs:
- GSE121212 as the session, paired by Patient, with the other four studies as comparison datasets (meta-analysis parity with the session on the covariate path).
- GSE121212 as the session for the all-pairs, Venn and UpSet outputs. These outputs operate on session contrasts, not on comparison datasets.
- GSE54456 (human) with GSE63310 (mouse) for orthologs and homology families, and GSE63310 as the session with GSE54456 as a comparison.
- Duplicate-detection probes.
- A re-grouping probe.

The R reference matches genes across studies by NCBI GeneID, independently of the app's symbol matching; all five matrices use the same NCBI GRCh38.p13 annotation. As in the app, the session's genes form the universe.

## 2. Results by feature

### 2.1 Per-dataset statistics of comparison datasets (check 1)

Reference: `xds_reference.R`. Steps: filterByExpr on all samples, grouped by the first design factor; TMM; `cpm(log = TRUE, prior.count = 2)`; `lmFit` on the contrast samples with `~group` or `~group + pairing factor`; `eBayes(trend = TRUE)`; BH.

All genes were compared for each of the 9 study × mode combinations (19,518 to 23,168 genes each). The kept-gene sets were identical. Maximum deviations over all nine combinations: log2FC 3.0e-14, SE 6.0e-13, moderated t 3.4e-12, p 2.0e-11 (relative), BH q 2.0e-11 (relative), prior df 6.5e-12. DE calls (FDR ≤ 0.05, |log2FC| ≥ 1) were identical in every case; the counts are in the table below. GSE63310 as a mouse comparison slot (LP vs Basal, 16,624 genes) also agreed: t 1.1e-12, p 1.6e-12 (relative).

| Study, mode | d0 | residual df | Up / down |
|---|---|---|---|
| GSE54456 | 4.848 | 169 | 1,293 / 1,513 |
| GSE186063 unpaired / paired | 4.361 / 4.413 | 52 / 25 | 1,015 / 1,399; 997 / 1,429 |
| GSE121212 unpaired / paired | 4.786 / 4.487 | 52 / 25 | 1,401 / 1,846; 1,420 / 1,825 |
| GSE83645 unpaired / paired | 2.931 / 3.019 | 23 / 19 | 932 / 919; 1,026 / 1,120 |
| GSE41745 unpaired / paired | 7.340 / 16.376 | 3 / 1 | 634 / 301; 867 / 690 |

With pairing, GSE41745 has one residual degree of freedom, so its moderated variances are dominated by the prior. This is a property of the design (2 vs 3 samples), not of the app.

### 2.2 Gene matching by symbol (check 9)

For each of the four comparison datasets, the app's symbol match (`matchAtoB`) was compared with GeneID identity for all 19,518 session genes. They agreed exactly: 19,308, 19,161, 19,440 and 16,294 matches. No session gene matched a different GeneID, no shared GeneID was missed, and neither the session nor any dataset had duplicate decoded symbols, including case-insensitive duplicates.

Of the session genes, 1,421 kept a bare NCBI GeneID as their name because the built-in table has no symbol for them. Of these, 1,071 are in current NCBI gene_info (946 ncRNA, 930 with LOC symbols), and 89 have an official nomenclature symbol (for example RSC1A1 and PRKCZ-AS1). This does not affect same-species matching, because both sides keep the same numeric ID. Of the 18,097 decoded symbols, 38 differ from the current NCBI symbol; most are mitochondrial genes (app MT-TF, NCBI TRNF).

### 2.3 Meta-analysis (check 2)

Every tested gene was compared in all 8 configurations: 16,150 genes with minK = 5 and 19,518 with minK = 2. Reference: `metafor::rma(yi, vi, method = "DL"|"FE")` and `predict(predtype = "Riley")` for DL with k ≥ 3, followed by BH and the app's signature rule. All 116 numerical and set checks passed.

| Quantity | Maximum deviation (8 configurations) |
|---|---|
| tested gene set, k per gene, n_up / n_down | 0 (exact) |
| pooled log2FC | 9.7e-14 |
| pooled SE | 2.1e-14 |
| z | 3.0e-12 |
| p (relative) | 8.2e-11; 83 to 100 genes per FE configuration have p < 1e-300 (0 in both app and R) |
| Cochran Q | 1.7e-11 |
| Q p (relative) | 8.3e-12 |
| I² (percentage points) | 1.2e-11 |
| τ² (DL) | 2.5e-13 |
| 95% PI (Riley, k − 2 df) | 2.1e-12; PI availability (k ≥ 3) identical |
| BH q (relative) | 8.2e-11 |
| shared signature membership | 0 (exact) |
| leave-one-out: signature size and number retained, for each omitted study | 0 (exact; recomputed with metafor on the remaining four studies) |

With GSE121212 as the paired session, the 16,150 genes agreed with metafor: pooled log2FC 4.8e-14, PI 2.2e-13, p 9.3e-11 (relative). The signature was the same 1,087 genes as with GSE54456 as the session.

Signature sizes (FDR ≤ 0.05, |pooled log2FC| ≥ 1, same direction in every study; PI criterion automatic, on for REM with 5 studies):

| Configuration | Genes tested | Meta FDR ≤ 0.05 | Signature (up) |
|---|---|---|---|
| paired, REM, minK 5 | 16,150 | 11,335 | 1,087 (607) |
| paired, REM, minK 2 | 19,518 | 13,508 | 1,316 (728) |
| paired, FE, minK 5 | 16,150 | 14,179 | 1,677 (816) |
| unpaired, REM, minK 5 | 16,150 | 11,369 | 1,143 (621) |
| unpaired, FE, minK 5 | 16,150 | 13,890 | 1,694 (831) |

Exports and displays (all eight configurations):
- **meta_analysis.csv:** every row agreed with the internal values within its rounding (6 significant digits; p and FDR 5; I² 2 decimals). The per-dataset log2FC, SE and FDR columns agreed with R, and the rows were sorted by p.
- **shared_signature.gmt:** the UP and DOWN sets equal the signature split by sign.
- **Heatmap:** the top 50 genes, their order, the per-dataset cells and the asterisks (per-dataset q ≤ 0.05) agreed with R.
- **Signature table:** the order of the top 100 by p agreed. The CI and PI cells agreed to the displayed 2 decimals.
- **Summary sentence:** the counts of tested genes, genes at meta FDR ≤ 0.05 and signature genes up and down agreed exactly.

### 2.4 Concordance and overlap (check 3), RRHO (check 4), CAT

All 10 study pairs were checked in both modes.
- **Inputs:** the matched log2FC and t vectors (1,442,328 values over both modes) agreed with R to 3.4e-12.
- **Exact counts:** matched-gene counts, hit-set sizes, overlaps and the numbers of concordant and discordant genes significant in both studies all agreed.
- **Displayed statistics:** Pearson r of log2FC and of t, and Spearman ρ of t, agreed to the displayed 3 decimals.
- **Hypergeometric p:** the displayed value agreed to its 3 significant digits. At full precision (`fisherExact2`) the 10 values with p ≥ 1e-300 agreed with `phyper` to 1.7e-10 (relative). The other 30 values are below 1e-300 and are shown as "<1e-300".

Paired run, Pearson r of log2FC: 0.733 to 0.913 (`VAL_XDS_concordance.csv`).

**RRHO:** the over-enrichment cells (z > 0) agreed with `RRHO::RRHO(alternative = "enrichment", log10.ind = TRUE)$hypermat` to 1.3e-9 (33,458 cells over the two modes). The signed under-enrichment extension (z ≤ 0) agreed with −log10 `phyper(lower.tail = TRUE)` under the package's N + 1 population to 4.0e-11.

**CAT:** the fractions agreed with R to 5.0e-16; the underlying counts are identical.

**cross_dataset_matched.csv:** agreed with the internal values to the 4 exported decimals.

### 2.5 Venn and UpSet (check 5)

These operate on session contrasts. The test used GSE121212 as the session (6 Group levels, 15 contrasts). Up and down hit counts for all 15 contrasts agreed exactly with R.
- **Venn:** region counts for 2 and 3 contrasts and for each direction setting (any, up, down; 30 regions) agreed exactly with set algebra on R hit sets, as did the displayed numbers, the set sizes and `venn_regions.csv`. For example, for the three contrasts PSO_lesional vs PSO_non_lesional, PSO_lesional vs CTRL_healthy and AD_lesional vs CTRL_healthy, the set sizes were 3,247, 4,391 and 1,749 and the triple intersection 1,293.
- **UpSet:** the counts of the 14 largest of 499 membership patterns agreed exactly.

### 2.6 Consensus table (check 6)

Compared with independent vote counting from the R results: 16,150 genes matched in all five studies; 12,425 consensus genes (unpaired) and 13,458 (paired). Gene sets, the number of studies significant per gene, the concordance flag, per-dataset log2FC (3.0e-14), row order, and `consensus_genes.csv` (4 decimals) all agreed. Displayed to 2 decimals, the pairwise log2FC correlation matrix agreed with `cor`.

### 2.7 Signature transfer (check 7)

Modules: the GSE54456 up-regulated (1,293) and down-regulated (1,513) DE genes, scored in each of the four other studies. Independent code computed the module score as the mean z of the matched rows of edgeR log-CPM, with z using the population SD as the app does; AUC from `scipy.stats.mannwhitneyu`; and the Welch t of the score. The app's seeded null was reproduced exactly by reimplementing mulberry32 and Floyd sampling in Python.
- **Genes and scores:** matched gene counts were identical, and per-sample scores agreed to 7.9e-15.
- **Displayed statistics:** AUC, t and empirical p agreed to their displayed precision. The empirical p values were identical as counts: (1 + #null)/201 equals 1/201 in 14 cases and 9/201 in 2.
- **Verdicts:** identical, "signature replicates" in all 16 cases.

An independent null of 1,000 random sets from a different generator gave p ≤ 0.039 in every case, so the verdict does not depend on the seed. The up module gave AUC 0.97 to 1.00; the down module gave AUC 0.00 to 0.02, the expected inverted direction.

### 2.8 Forest plot (check 8)

`drawForest(g)` was called for all 19,518 session genes in the paired run (93,721 per-study rows):
- **Points:** agreed with the R log2FC to 3.0e-14.
- **95% CI half-widths:** agreed with `qt(0.975, df.total) × SE` to 1.1e-12.
- **Pooled diamond (μ ± 1.96 SE):** agreed with metafor DL to 6.1e-14. It is drawn exactly for the 16,150 genes in the meta-analysis.
- **Red markers:** exactly the studies with q ≤ 0.05.

### 2.9 Duplicate-dataset detection (check 10)

- **Exact copies:** a second copy of the session file and a second copy of GSE186063 were flagged (`dupOf` = session and GSE186063) and left unticked in the meta-analysis selector.
- **Partial and renamed copies:** a copy of GSE54456 restricted to 100 of its 171 samples, and a copy of GSE186063 with renamed sample IDs, were not flagged and were ticked by default.
- **Effect on the meta-analysis:** with GSE54456 + GSE186063, the signature had 2,242 genes. Adding the renamed GSE186063 copy changed it to 2,188 genes (PI off by default with 3 datasets) or 678 (PI on), and the copy was counted as independent evidence.

### 2.10 Orthologs and homology families (check 11)

The mapping was verified on every pair, not on a sample.
- **Embedded table against MGI:** the 18,782 pairs of the embedded `ORTHO_MH` table equal all MGI homology classes with exactly one mouse and one human member.
- **Pairs used in matching:** all 13,123 pairs used to match GSE54456 to GSE63310 agree with those classes, and both symbols are current in NCBI gene_info.
- **Homology families:** the 987 families in `ORTHO_FAM` equal an independent union-find over MGI classes that share a member, keeping families with more than one member in a species.
- **Family statistics:** for 492 tested families (GSE54456 Psoriasis_skin vs normal_skin; GSE63310 LP vs Basal), the module-score differences agreed to 2.9e-12, Welch p and within-dataset BH q to 1.5e-11 (relative), and the concordance flags exactly.

The class-level 1:1 rule is not one-to-one at the gene level.

### 2.11 Known psoriasis biology

Classic lesional genes (18 genes, the task list) in the default shared signature (paired, REM, minK = 5, PI on): 13 of 18. DEFB4A was excluded by the prediction interval. CXCL8, IL17A, IL23A and CCL20 were not tested, because filterByExpr removed them in GSE41745 (raw counts 0 to 16 per sample at 10.7 to 13.3 million reads). With minK = 2: 15 of 18; DEFB4A, CXCL8 and IL23A were excluded by the PI. With FE: 14 of 18, with the same four genes untested. Every tested classic gene had a positive pooled effect. None of the ten sex-chromosome marker genes was in any signature tested (four configurations), and no mitochondrial or ribosomal-protein gene was in the default signature.

The leading Hallmark sets for the default up-signature were interferon-γ response (55 of 186 genes, q 1.5e-32), G2M checkpoint, interferon-α response, E2F targets, IL6–JAK–STAT3 and inflammatory response (hypergeometric test; universe = tested genes; app's built-in Hallmark library). The down-signature reached FDR ≤ 0.05 only for KRAS_SIGNALING_DN (q 0.027) and myogenesis (q 0.046).

**Leave-one-study-out stability** (paired, REM, minK = 5): the share of the 1,087-gene signature retained was 51% without GSE54456, 68% without GSE186063, 64% without GSE121212, 65% without GSE83645 and 67% without GSE41745. With the PI criterion off, the same omissions retained 88%, 96%, 91%, 97% and 96% of a 1,594-gene signature. Several classic genes are lost in particular omissions: without GSE41745, S100A7, S100A8, S100A9, S100A12, SERPINB4 and PI3 drop out.

**Published meta-analysis signature:** the MAD-5 gene list (Tian et al. 2012, PLoS One 7:e44274, Table S2) could not be retrieved. The supplement redirects to storage.googleapis.com, which is on the sandbox denylist. MAD-5 is therefore not assessed. The article's Table 1 (top 25 up and top 25 down genes of MAD-3, with GeneIDs) was retrieved from the article PDF and compared by GeneID; this is a partial list from a related signature, not MAD-5.
- **Default signature:** 29 of 50 present (18 up, 11 down), none in the opposite direction. Fifteen were excluded only by the PI, 4 were untested (absent from GSE41745), 1 failed the direction rule (KRT77: +0.10, q 0.81 in GSE41745) and 1 is not in the session (COL6A4P1).
- **FE signature:** 44 of 50.

## 4. Caveats

- Consensus inclusion uses per-study FDR only, without a fold-change threshold, so 13,458 of 16,150 genes qualify in the paired run. The previous audit kept this rule and documented it as vote counting.
- The same-direction rule counts non-significant near-zero estimates. KRT77 fails it on a log2FC of +0.10 (q 0.81) in GSE41745. The previous audit measured this rule as not material.
- The forest-plot diamond shows the 95% CI, not the PI that the signature rule uses.
- The rank–rank Spearman ρ uses positional ranks (no tie averaging); it agreed with R to the displayed 3 decimals.
- With pairing, GSE41745 leaves 1 residual df. GSE186063 pairs are inferred from GEO titles (kit notes), and its normal_skin samples (ankylosing spondylitis patients) were not used.

## 5. Not assessed

- The MAD-5 gene list (Table S2 of Tian et al. 2012): its download host, storage.googleapis.com, is on the sandbox denylist. Only the MAD-3 top-50 list from the article's Table 1 was compared.
- The HGNC HCOP ortholog file: the download returned an error page. Orthologs were verified against MGI only, as a single source.
- Loading a design file after the comparison matrix (the `hookCmpFile` design branch, which uses only the first factor and does not refresh covariates). It needs a FileReader, which the harness does not provide.
- The custom ID map path (`parseCustMap`).
- voom and Welch as comparison-dataset methods; comparison datasets always use the moderated t.
- Comparison datasets supplied as log-scale or already-normalized input.
- Gene Explorer panels for comparison datasets (`drawGeneB`).
- Session save and restore of Compare state.
- Behavior in real browsers. All runs used JavaScriptCore with a DOM stub.
- Cross-species pairwise concordance statistics, apart from the matched-gene count. The pairwise statistics use the same code path as same-species pairs, which was validated.

## 6. Scripts

| File | Purpose |
|---|---|
| `scripts/xds_reference.R` | per-dataset R references |
| `scripts/xds_meta_reference.R` | metafor, all configurations and leave-one-out |
| `scripts/xds_concordance_reference.R` | Pearson, Spearman, hypergeometric overlap, RRHO, CAT |
| `scripts/xds_venn_reference.R` | all-pairs, Venn and UpSet references |
| `scripts/xds_mouse_reference.R` | GSE63310 slot and the re-grouping probe |
| `scripts/js/*.js` | app driver and run bodies |
| `scripts/xds_compare.py` | CI runner, writes `VAL_XDS_checks_ci.csv` |
| `scripts/README_XDS.md` | how to run |

A test run of `xds_compare.py` gave 348 checks with 0 parity failures. One known-defect probe, gene-level 1:1 orthologs, fails as expected. That check is flagged `probe` so that it does not fail CI.
