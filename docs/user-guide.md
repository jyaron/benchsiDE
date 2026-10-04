# benchsiDE User Guide

This guide takes you from an expression matrix to reportable results, one step at a time. Each tutorial names the exact controls you will use. Concepts, statistical definitions and troubleshooting are covered in more depth in the [Knowledgebase](../wiki/Home.md).

Version covered: 1.0.0-rc.2. The application is at <https://www.benchside.org>; source code and releases are at <https://github.com/jyaron/benchsiDE>.

## Contents

1. [Before you start](#1-before-you-start)
2. [Preparing your input files](#2-preparing-your-input-files)
3. [Tutorial A: your first analysis with the demo dataset](#3-tutorial-a-your-first-analysis-with-the-demo-dataset)
4. [Tutorial B: a two-group experiment with a batch covariate](#4-tutorial-b-a-two-group-experiment-with-a-batch-covariate)
5. [Tutorial C: a time course](#5-tutorial-c-a-time-course)
6. [Tutorial D: gene-set testing](#6-tutorial-d-gene-set-testing)
7. [Tutorial E: the Discovery screen](#7-tutorial-e-the-discovery-screen)
8. [Tutorial E2: sample traits](#8-tutorial-e2-sample-traits-clinical-scores-ihc-counts-and-other-per-sample-measurements)
9. [Tutorial F: comparing datasets and species](#9-tutorial-f-comparing-datasets-and-species)
10. [Tutorial G: figures, reports, sessions and methods text](#10-tutorial-g-figures-reports-sessions-and-methods-text)
11. [Tutorial H: verifying the statistics in your own browser](#11-tutorial-h-verifying-the-statistics-in-your-own-browser)
12. [Checklist before you report results](#12-checklist-before-you-report-results)

---

## 1. Before you start

### 1.1 What you need

- A current version of Chrome, Firefox, Safari or Edge. Nothing needs to be installed.
- An expression matrix: genes in rows, samples in columns. Raw counts are preferred.
- Optionally, a design file that assigns each sample to a group.

### 1.2 Choose how to run benchsiDE

| Option | When to use it | How |
|---|---|---|
| Hosted page | Everyday use | Open <https://www.benchside.org>. The page loads once; all computation then runs locally. |
| Downloaded file | Reproducible, version-pinned work | Download `index.html` from a GitHub release and open it by double-clicking. |
| Offline build | Air-gapped machines, controlled-access data | Download `benchside-offline.html` from a GitHub release. It contains the plotting library and makes no network requests. |

### 1.3 Confirm that your data stays local (optional)

1. Open the browser developer tools (F12, or Cmd+Option+I on macOS) and select the **Network** tab.
2. Load your data and run an analysis.
3. Observe that no request carries your data. With the offline build you can disconnect from the network entirely.

Both builds also carry a Content-Security-Policy that the browser enforces: the page may not open any network connection, and may load nothing from another server except, in the hosted build, the plotting library from its CDN, pinned by an integrity hash. The offline build allows no server at all. Any blocked attempt appears in the browser console as a policy violation. See the [Privacy and architecture](../wiki/Privacy-and-Architecture.md) page.

---

## 2. Preparing your input files

### 2.1 Expression matrix

- Tab-separated (`.tsv`, `.txt`) or comma-separated (`.csv`).
- First column: gene identifiers. Accepted: gene symbols, Ensembl gene IDs (`ENSG…`, `ENSMUSG…`, `ENSRNOG…`; version suffixes such as `.12` are tolerated), or NCBI/Entrez gene IDs (plain integers). The standard GEO files named `*_raw_counts_*_NCBI.tsv` load without modification.
- Remaining columns: one per sample, with the sample name in the header row. Extra annotation columns (gene name, biotype) are detected by their content and set aside.
- Values: raw counts are best, because they enable TMM normalization, filterByExpr and voom. CPM/TPM/FPKM and log-scale values are also accepted.

Example (first rows of the included demo file):

```
EntrezID   GSM1545535  GSM1545536  GSM1545538  ...
497097     1           2           342         ...
27395      431         771         1368        ...
```

### 2.2 Design file (optional but recommended)

Either orientation is accepted.

Samples in rows:

```
sample       celltype   batch
GSM1545535   LP         b1
GSM1545536   ML         b1
GSM1545538   Basal      b2
```

Factors in rows:

```
            GSM1545535  GSM1545536  GSM1545538
celltype    LP          ML          Basal
batch       b1          b1          b2
```

Rules:

- Sample names must match the matrix column headers exactly.
- The first factor is used for grouping by default; you can select a different one on the review step. Other factors become available as covariates.
- Any line-ending convention (Windows, Unix, classic Mac) is accepted.

If you have no design file, benchsiDE infers groups from sample names; you can correct them on the review step.

### 2.3 Gene-set libraries (optional)

Built-in libraries (from MSigDB v2024.1, CC BY 4.0): Hallmark, GO Biological Process, GO Cellular Component, GO Molecular Function and Reactome, each for human and mouse. You can also load any `.gmt` file.

### 2.4 Annotation and identifier maps (optional)

For mouse, human and rat, built-in annotation converts Ensembl and Entrez IDs to gene symbols. For other species, or for probe IDs, supply a two-column TSV (identifier, symbol).

---

### 2.5 Large datasets

benchsiDE runs entirely in the browser, so computation time grows with the size of the data. As a guide, GEO GSE54456 (171 samples, 19,518 genes after filtering) takes about 5 s to analyze, and FRY on GO Biological Process (7,608 sets) takes about 7 s on a recent laptop. Long computations show a progress bar with the current step, elapsed time and an estimate of time remaining; the page remains responsive while they run. Chrome, Edge, Firefox and Safari are supported; a 64-bit browser with at least 8 GB of system memory is recommended for datasets of this size.

After **Analyze →**, the line above the dashboard gives the load time in your browser: reading and parsing each file, and the analysis up to the first display (time spent on the review step is not counted). **copy** puts these timings and your browser details on the clipboard. For a standard comparison, **Validate statistics → Timing benchmark…** runs a fixed analysis on synthetic datasets of increasing size; the same datasets can be exported to time other tools on the same computer.

## 3. Tutorial A: your first analysis with the demo dataset

The repository includes a public dataset (GEO GSE63310: mouse mammary basal, luminal progenitor (LP) and mature luminal (ML) cells, three replicates each). The numbers quoted below were reproduced independently in R (edgeR and limma), so you can use them to confirm that each step went as intended.

### Step 1. Open the application

Open benchsiDE. The start screen shows numbered cards: **1 · Expression matrix**, **2 · Experimental design** and **2b · Species & gene symbols**.

### Step 2. Load the matrix

Drag `demo/GSE63310_counts.tsv` onto card 1, or click the card and select the file. A summary appears showing 9 samples and the detected input type (raw counts).

### Step 3. Load the design

Drag `demo/GSE63310_design.tsv` onto card 2. The sample-to-group table fills with Basal, LP and ML.

### Step 4. Set the species

On card 2b, choose **Mouse** in the species selector and click **Use built-in annotation**. The Entrez IDs are converted to gene symbols. If some identifiers have no entry, a link lets you download the unmapped IDs for inspection.

### Step 5. Review the setup

Card **3 · Review & analyze** shows the settings that will be used. Check each one:

1. **Design type**: leave **Auto-detect**. Cell types are categories, so the design is treated as categorical. For time courses or dose series choose **Time course / dose–response (ordered numeric)**.
2. **Normalization**: **TMM (edgeR method)**.
3. **Values are**: **auto-detect** should read raw counts.
4. **Filter**: **edgeR filterByExpr (recommended for counts)**.
5. **Group table**: confirm each sample's group. To reassign several samples at once, type a pattern into **match**, click **select matching**, enter a group name under **assign selected to**, and click **Apply**. **Reset to inferred** restores the original assignment.

### Step 6. Analyze

Click **Analyze →**. With filterByExpr, 16,624 genes are retained. All tabs are populated, and the sample bar below the header lists every sample with a checkbox.

### Step 7. Inspect quality control (Overview & QC tab)

Work down the cards:

1. **Library size** and **Detected genes per sample**: look for samples with unusually low values.
2. **PCA**: samples should cluster by cell type. Click **PC loadings (driver genes)** to see which genes drive each component.
3. **Sample dendrogram**, **Sample–sample correlation** and **MDS**: three independent views of sample similarity. A true outlier appears in all three.
4. **RLE** and **expression density per sample**: after normalization, boxes should be centred on zero and density curves should overlay.

Warnings about low depth or outlier samples appear at the top of the tab.

### Step 8. Exclude a sample (if needed)

Untick a sample in the sample bar. Every statistic in every tab recomputes immediately. Tick it again to restore it. Each group block also has an **all** checkbox. For this demo, keep all nine samples.

### Step 9. Run differential expression (Differential Expression tab)

1. Set the baseline group (left selector) to **Basal** and the comparison group (right selector) to **LP**.
2. Set the FDR selector to **0.05** and the **|log₂FC| ≥** threshold to **1**.
3. Choose the method **voom**.
4. Read the result: 2,781 genes up and 3,301 down in LP relative to Basal (variance estimated from all three groups, the default; with **two groups only** it is 2,810 and 3,338).
5. Switch between **Volcano** and **MA plot**. Set **label top** to the number of genes to annotate, and choose whether to rank them by **p-value**, **|log₂FC|** or **π score**.
6. To list a group of genes, set **drag to** to **select (box)** or **select (lasso)** and drag over them on the volcano or MA plot (Shift adds to the selection; double-click clears it). The table below the plot lists the selected genes with log₂FC, group means, p, FDR and the significance call; click a column header to sort. The selection is kept when you switch between volcano and MA plot. From the selection you can export a CSV, copy the gene symbols, draw a heatmap of the selected genes, or compare up to 10 of them in Gene Explorer.
7. Open the **p-value histogram (diagnostic)** and **voom mean–variance trend** panels to check the model assumptions.
8. Click **Export full table (CSV)** to save all genes with log₂FC, t, p and FDR.


**Variance from all groups.** With more than two groups, the default (*variance from: all groups*) fits every group and compares the two chosen groups by a contrast, as limma does. This uses all samples to estimate the variance and is usually more powerful. Choose *these two groups only* if other groups are known to be much more variable.

### Step 10. Look at individual genes (Gene Explorer tab)

Click any gene in the top-genes table. The Gene Explorer opens with per-group values. Use the toggles to switch between **log₂** and **linear** scale and between **± SEM** and **± SD**. Choose **pairwise vs first group** or **all pairs (sig. only)** to add Holm-corrected pairwise tests to the plot. **← back** returns you to the previous view.

The **bar / box / violin** toggle sets how each group is drawn; the sample points are shown on every type. Bars show the mean ± SEM or SD; boxes show the median and quartiles with whiskers to the most extreme sample within 1.5 × IQR; violins show a kernel density estimate with the median and quartiles. With fewer than about 5 (box) or 8 (violin) samples per group these summaries rest on very few values, and the caption says so.

**+ Add to comparison** collects genes for the multi-gene comparison below the plot. It is drawn as grouped bars by default, or as grouped box or violin plots, or as profile lines through the group means. **Group by** sets the arrangement: by condition (genes side by side within each condition) or by gene (conditions or time points side by side within each gene). **z-score each gene** puts genes with different expression levels on one axis. The Excel export of this plot is a Prism Grouped table in the arrangement shown, with the samples as replicate sub-columns.

### Step 11. Compare all contrasts

In the Differential Expression tab, click **Run all pairs (current method)**. The table and UpSet plot show the hits for every pairwise contrast. Select two or three contrasts and click **Draw Venn** for exact region counts; click a region to list its genes.

### Step 12. Save your work

Click **Save session** in the header. Keep the session file with your data; loading it later restores every setting.

---

## 4. Tutorial B: a two-group experiment with a batch covariate

Use this workflow when samples were processed in batches, or when a known variable (sex, RNA quality) should be adjusted for.

### Step 1. Include the covariate in the design file

```
sample   condition   batch   RIN
S1       Ctrl        b1      7.8
S2       Ctrl        b2      8.4
S3       Treated     b1      7.1
S4       Treated     b2      8.0
```

### Step 2. Load and analyze

Load the matrix and design, confirm on the review step that the grouping factor is `condition`, and click **Analyze →**.

### Step 3. Check whether the batch effect is visible

On the Overview & QC tab, look at the PCA. To see the data with the batch effect removed, set **display batch-adjusted** to the batch factor. This setting changes the display only; no statistic uses adjusted values.

### Step 4. Fit the adjusted model

1. In the Differential Expression tab, choose **Moderated t (eBayes)** or, for raw counts, **voom**. Covariate adjustment is available for both; the Welch t shows an explicit unadjusted warning.
2. Select the covariate(s) under **adjust for**. A numeric covariate with more than two distinct values (such as RIN) is fitted as a continuous slope; a categorical covariate is fitted as indicator terms.
3. The note under the plot states which covariates were fitted and the residual degrees of freedom.

If a covariate is confounded with the groups (for example, every control in batch 1 and every treated sample in batch 2), benchsiDE refuses to fit the model and explains why. No statistical method can separate the two effects in that case.

**Example: reproducing the limma/Glimma/edgeR workflow.** Law et al. (F1000Research 5:1408, 2016) analyse the demo dataset with sequencing lane as a covariate and TREAT. Add a `lane` column to the design file (GSM1545535, GSM1545536 and GSM1545538: L004; GSM1545539 to GSM1545542: L006; GSM1545544 and GSM1545545: L008), load the demo, and in the Differential Expression tab choose **voom**, adjust for **lane**, set FDR to **0.05** and |log₂FC| to **1**, and tick **test the threshold (TREAT)**. In Compare all contrasts, a Venn diagram of LP vs Basal and ML vs Basal then shows 3,647, 3,831 and 2,782 genes, the values limma 3.66 gives for the article's code (the article, computed with an earlier Bioconductor release, reports 3,648, 3,834 and 2,784). This analysis is part of the repository's external validation.

### Step 4b. Test the fold-change threshold (TREAT)

The default rule calls a gene when FDR is at or below the chosen level and the estimated |log₂FC| reaches the threshold. Ticking **test the threshold (TREAT)** instead tests whether |log₂FC| is greater than the threshold (McCarthy & Smyth 2009), as limma's `treat()`. This is stricter and is the approach recommended by the limma authors when a minimum fold change matters. The table, volcano plot, CSV, Venn diagram and methods text then report TREAT p-values and FDR, and say so.

The **Design check** (review step, and under the sample list after analysis) shows this before you fit anything. It tabulates each factor against the groups and lists, for each pair of groups, whether the factor can be adjusted for or is confounded. For a confounded comparison you can:

1. restrict the samples to one level of the factor with the **Include samples with …** checkboxes under the sample list;
2. group by a combined factor (for example `Type × Diag`) in the review step and compare within one level;
3. compare only groups in which the factor varies, with the factor under **adjust for**; or
4. report the comparison as confounded.

Example: in GEO GSE186063 every normal-skin sample is from an ankylosing spondylitis patient, and lesional and non-lesional samples are from psoriasis and psoriatic arthritis patients. Diagnosis can be adjusted for in lesion vs non-lesion, but a comparison of lesion with normal skin is also a comparison of diagnoses.


Each covariate has a **continuous / categorical** selector. Patient, subject, donor and pair columns are categorical by default, even when they are numbered, so they are fitted as blocking factors. Use continuous only for measurements such as age or RIN.

### Step 5. Paired designs

For paired samples (the same subject before and after treatment, or lesional and non-lesional biopsies from each patient), add a column identifying the subject and select it under **adjust for**. This fits subject as a blocking factor. If the design file has no subject column, add one on the review step: open **Edit the design**, type the column name, click **Add column**, and fill in each sample's subject.

---

## 5. Tutorial C: a time course

### Step 1. Encode time as a number

In the design file, use numeric labels such as `0`, `2`, `4`, `7`, `14`. On the review step, set **Design type** to **Time course / dose–response (ordered numeric)**. The review step previews how the choice will be used.

### Step 2. Analyze and read trends

After **Analyze →**, gene plots use a true numeric time axis, group means are connected, and Spearman trend statistics are reported. With a categorical design none of these appear, because ordering categories would be meaningless.

### Step 3. Find genes that change at any time point (Group Patterns tab)

Click **Compute** under **Genome-wide moderated F (any group differs)**. This single test uses all samples and is the appropriate first test in a multi-group design.

### Step 4. Group genes by temporal profile

1. Choose the number of genes (**500**, **1000** or **2000** most variable).
2. Leave the number of patterns at **auto (silhouette)**. The note reports the silhouette score for each k; a best score below about 0.3 means the profiles form a continuum rather than distinct clusters.
3. Click a pattern line to list its genes; **Export cluster genes (CSV)** saves the assignment.
4. The **Trajectory heatmap** shows every selected gene over time.

### Step 5. Contrast individual time points

In the Differential Expression tab, compare each time point with baseline, then use **Run all pairs (current method)** and **Compare two contrasts** to see which responses persist and which resolve.

With few replicates at baseline, genome-wide FDR can be unattainable even for large effects. See [Statistical power](../wiki/Statistical-Power.md) in the Knowledgebase.

---


### Step 6. Interaction: does a response differ between conditions?

For a factorial design (for example treatment × genotype, or treatment × time point), combine the two factors into one grouping under **Review → group by**, so that each group is one combination such as *Aldara_WT*. In the Differential Expression tab, the **Interaction** card tests the difference of differences (B2 − A2) − (B1 − A1):

1. Set the reference condition, for example Control_WT → Aldara_WT, and the compared condition, Control_KO → Aldara_KO.
2. Choose the test, whether the variance is estimated from all groups, and any covariates (for example Day), and click **Run**.
3. Read the table: the interaction log₂FC is the change in the response, and the two simple effects show whether a gene responds in one condition only. A positive value means the response is larger in the compared condition.

An interaction is estimated less precisely than a single comparison, so it needs larger groups; the card warns when a group has fewer than three samples.

## 6. Tutorial D: gene-set testing

### Step 1. Choose a library (Enrichment tab)

Under **1 · Gene-set library (GMT)**, choose a built-in library for your species or drop your own `.gmt` file.

### Step 2. Choose the question you are asking

| Method | Question | Input |
|---|---|---|
| **ORA (gene list)** | Are my significant genes over-represented in this set? | The current DE hit list |
| **GSEA** | Are the genes of this set concentrated at the top or bottom of the ranking? | All genes ranked by the moderated t-statistic; significance from permuting sample labels |
| **FRY (self-contained)** | Is the set as a whole differentially expressed? | Per-sample expression of the set's genes |

**GSEA significance.** By default GSEA permutes the sample labels (Subramanian et al. 2005). This keeps the correlation between genes, so the FDR is valid for sets whose members are co-regulated. It needs at least 1,000 distinct labelings (about 7 samples per group) and no covariates; otherwise, or if you choose **permute: genes**, the app uses gene permutation (preranked); it then reports no FDR, and the table is a ranking by NES and nominal p only. On random splits of normal skin, gene permutation reported 12–37 of the 50 Hallmark sets at FDR ≤ 0.05; sample permutation reported none. Sample permutation is conservative when a set is strongly co-regulated, so results are shaded at FDR ≤ 0.25, the threshold the GSEA documentation uses; FRY and CAMERA (Discovery tab) remain the recommended tests for inference.

FRY is usually the most sensitive of the three when replicates are few. FRY finding signal where ORA or GSEA do not is expected, not a contradiction: the tests answer different questions (see [Gene-set testing](../wiki/Gene-Set-Testing.md)).

### Step 3. Run and read the results

Click **Run enrichment**. The table lists each set with its statistic, p-value and FDR. Supporting plots follow:

1. **Ridgeline**: the distribution of log₂FC within each top set.
2. **Barcode plot** and **set-vs-rest ECDF**: click **bc** in a table row, or a ridgeline row, to draw them for that set.
3. **Jaccard heatmap**: overlap between the top sets. Click a cell to list shared genes.
4. **Leading-edge matrix** (GSEA): which genes drive which sets.

Each method has its own table and plots:

| Method | Table columns | Plots |
|---|---|---|
| ORA | size, overlap, ratio, p, FDR, overlapping genes | lollipop chart, ridgeline, Jaccard heatmap |
| GSEA | size, NES, p, FDR (sample permutation only), leading-edge genes | lollipop chart, running-score plot, ridgeline, Jaccard heatmap, leading-edge matrix |
| FRY | size, direction, p and FDR (directional), p and FDR (mixed) | ridgeline, Jaccard heatmap |

FRY's directional p asks whether the set's genes change together in one direction; the mixed p asks whether they change in either direction. A set can have a large directional p and a very small mixed p when its genes change strongly in opposite directions. With many samples or a strong treatment, FRY finds most sets significant; that is a property of a self-contained test, and CAMERA or GSEA show which sets change more than the other genes.

### Step 4. Test every contrast at once

With three or more groups, click **Enrichment across all contrasts**. The heatmap shows signed −log₁₀ FDR for each set in each pairwise contrast (red: up, blue: down). Click a cell to open the barcode plot for that set and contrast. **Export matrix CSV** saves all sets, not only the 25 displayed.

---

## 7. Tutorial E: the Discovery screen

The Discovery screen looks for groups of related genes that move together even when few of them pass a per-gene test.

### Step 1. Read the explanation

Open **How discovery works and how to read the results** at the top of the Discovery tab.

### Step 2. Choose candidate modules

Candidate modules can be gene families (genes sharing a symbol prefix), sets from the loaded GMT library, and HGNC gene groups. Set the size limits and the number of permutations.

### Step 3. Run

Click **Run discovery**. For two-group designs, each module is tested with CAMERA, which adjusts for correlation between genes. For designs with more groups, seeded permutation tests are used and the inter-gene correlation is displayed for each module.

### Step 4. Read the result correctly

- A module tagged **hidden module** is significant as a module while fewer than 25% of its members are individually significant.
- When no module survives FDR, which is common with small groups, the results are shown as a nominal tier and labelled accordingly. Treat these as ranked hypotheses to test in new data, not as findings.
- Each card offers **Panel figure (top 9)** and **Heatmap (all N)** exports. **CSV** exports the full result table.

---

## 8. Tutorial E2: sample traits (clinical scores, IHC counts and other per-sample measurements)

1. Add the measurement as a column of the design file, one value per sample (blank or NA where it was not measured), or add it in the review step with **Edit the design → Add column**.
2. Open the **Sample traits** tab and choose the trait. Choose a transformation if needed: log₂(x + 1) for counts, logit for percentages, rank when a few samples have extreme values.
3. If the trait differs between groups (PASI is 0 in healthy controls), tick **adjust for group** or keep one group only (for example lesional samples). Otherwise the genes found are mostly those that differ between the groups.
4. If a subject contributes several samples, choose the subject column under **repeated samples per**. Leave the analysis on **auto**: subject-level traits are analysed on subject means, and traits measured per sample with a mixed model.
5. Click **Test association**. Click a gene in the table or volcano plot to see expression against the trait with the fitted slope. Genes marked in the Cook's D column depend on one sample.
6. The arrow next to a gene, or **Open in Gene Explorer** under the plot, opens the gene in Gene Explorer; **← back** returns to the Traits tab.
7. With a gene-set library loaded (Enrichment tab), **Gene sets (CAMERA)** tests whole sets. **Module × trait matrix** tests hub neighbourhoods, gene sets or Discovery modules against every trait at once, and shows how the traits correlate with each other. Click a cell, or a gene-set name, to list its genes with links to Gene Explorer.

## 9. Tutorial F: comparing datasets and species

The Compare datasets tab relates your session (dataset A) to any number of other datasets. There is no fixed limit on the number of datasets; the limit is memory. Each dataset uses about 8 bytes per gene per sample (roughly 30 MB for 20,000 genes × 170 samples). The slot list shows the total, warns above 1.5 GB, and refuses a dataset that would take the total above 3 GB. Each dataset is filtered, normalized and tested on its own; raw values are never merged across datasets.

### Step 1. Add a dataset

1. Set **species of next dataset** (same as session, mouse or human).
2. Click **Add dataset (matrix)** and select the file.
3. Optionally click **Design for last added** and select its design file.

The slot reports how many genes matched. If the count is low, the message names the identifier types on both sides and the remedy.

Each dataset has a name field. A GEO file such as `GSE54456_raw_counts_GRCh38.p13_NCBI.tsv` is named `GSE54456` by default. Type any other name (for example `Li et al. 2014`) and press Enter; it is used in every plot, table, export and the methods text. The session's own dataset is named in the field above the slots (default "A (this session)"), and that name is saved with the session.

### Step 2. How genes are matched

- Same species: by gene symbol, case-insensitive.
- Mouse and human: by MGI strict one-to-one orthologs, which handles renamed genes (for example Trp53 and TP53).
- **Custom ID map**: a two-column file that overrides both, for other species or special cases.

### Step 3. Compare a pair

Choose X and Y and click **Compare pair**. You will see:

1. The fold-change scatter, coloured by quadrant, with ranked lists of shared and discordant genes.
2. t-statistic and rank–rank concordance, and a Fisher test of hit-list overlap.
3. An RRHO map of rank–rank overlap.
4. A Bland–Altman plot.

**Export matched table** saves every matched gene with its statistics in both datasets.

### Step 4. Summarize across all datasets

Click **Compute across all datasets** for the correlation matrix and consensus table. Type a gene symbol in the **forest plot** field and click **Draw** for its effect in every dataset. The axis always includes 0, so effect sizes are compared on a common scale. Gene Explorer also shows the selected gene in every comparison dataset, each with its own groups and statistics.

### Step 5. Meta-analysis: shared signatures

The meta-analysis card combines the datasets' own results into one pooled estimate per gene. It does not merge expression values.

1. Tick the datasets to include and set each dataset's contrast in its slot. If a comparison dataset's design file has several columns, the slot lets you choose which column defines the groups and which to **adjust for / pair by**. For paired samples (lesional and non-lesional skin from the same patients), tick the patient column: the fold change and its standard error then come from the within-patient comparison, as in the session's own analysis. Without it, a paired dataset is analyzed as unpaired and its standard errors are larger (1.6-fold in a simulated 8-patient example), which gives it less weight in the pooled estimate.
2. Choose the model and the test. The default is **random effects with REML** and the **Hartung–Knapp** test, the appropriate choice when datasets differ in platform, tissue, cohort or species and when there are few datasets: on null data it kept false positives at the nominal rate with two to five small datasets, where the DerSimonian–Laird model with the z test did not. DerSimonian–Laird, the z test and a fixed-effect model (one common effect) remain available.
3. Choose whether a gene must be matched in all datasets or in at least *k*.
4. Set the shared-signature criteria: meta FDR and minimum |pooled log₂FC|; optionally, that the 95% prediction interval excludes 0 and a maximum I² (both off by default).
5. Click **Run meta-analysis**.

Results:

| Output | Content |
|---|---|
| Meta volcano | pooled log₂FC against −log₁₀ meta p; colour = I²; the shared signature is outlined |
| Heatmap | the signature genes (up to 100, ranked by |pooled log₂FC| or by meta p): log₂FC in each dataset and pooled; * = FDR ≤ 0.05 in that dataset |
| Table | pooled log₂FC with 95% confidence and prediction intervals, p, FDR, I², τ², direction agreement and each dataset's log₂FC; sortable by any column |
| Leave-one-dataset-out | the signature recomputed with each dataset omitted, and how much of it is retained |
| Exports | all tested genes (CSV); the signature as up and down gene sets (GMT), which can be loaded as a gene-set library |

A gene enters the shared signature only if it is significant after pooling, has a pooled effect at least the chosen size, and changes in the same direction in every dataset. With the prediction-interval option, the 95% prediction interval, the range expected for the effect in a new comparable study, must also exclude 0. This asks whether the direction replicates while allowing the size of the change to differ between studies. I² is not a good default filter: large studies measure fold changes so precisely that a log₂ fold change of 9 in one study and 7 in another (512-fold and 128-fold) counts as heterogeneous, so an I² limit removes the largest, best-established effects (for psoriasis, S100A7–9, SERPINB3/4, PI3). Click any gene for its forest plot, which then shows the pooled estimate as a diamond. With two or three datasets, I² and τ² are imprecise; read the leave-one-out table before relying on the signature.

**Power with few datasets.** The Hartung–Knapp test refers each pooled estimate to a t distribution on k − 1 degrees of freedom, where k is the number of datasets, and the FDR is controlled across thousands of genes. In benchsiDE's simulations it detected almost no true effects with two or three datasets, 5–37% with four and 26–57% with five, while keeping false discoveries near the nominal rate even when effects differed between datasets. The z test detected far more but did not control the FDR when effects differed between datasets. The summary states this whenever two to four datasets are combined. With two or three datasets, report the pooled estimates and the per-dataset results; if you need a signature for exploration, choose **test: z (normal)** and say so in the methods.

**Samples behind the result.** Each dataset's checkbox shows its group sizes. After the run, a table lists for each dataset the samples compared, other samples used only for the variance estimate, samples left out for a missing covariate value, the number of subjects when a patient or pair column is adjusted for, the dataset's median share of the pooled weight, and the number of genes it contributed. The gene table and CSV give the number of samples behind each gene's pooled estimate.

**Leave-one-dataset-out.** For each omitted dataset the table gives the signature recomputed without it and, in the last column, the number of signature genes whose pooled estimate keeps its direction at nominal p ≤ 0.05. The last column does not depend on the correction across genes, so it is the better guide when the Hartung–Knapp test loses a degree of freedom.

**More views of the result.**

- **heatmap**: 25, 50 or 100 signature genes, ranked by |pooled log₂FC| or meta p.
- **multi-gene forest**: type gene symbols, or leave the box blank for the signature genes with the largest pooled effects, choose how many, and click **Draw**. Each gene shows its estimate in every dataset and the pooled estimate with its interval; the axis is symmetric about 0.
- **table**: the shared signature or all tested genes, 50 to 500 rows.
- **pathway enrichment of**: choose the shared signature or the top N genes (ranked by meta p or by |pooled log₂FC|), and up, down or both, then click **Run in Enrichment tab**. The background is the genes tested in the meta-analysis. As with any top-ranked list, these p-values are optimistic for co-regulated genes.

### Step 6. Score homology families (cross-species)

Many genes have no one-to-one ortholog (for example the mouse Serpinb3a to Serpinb3d genes and human SERPINB3 and SERPINB4). These families are excluded from gene-level matching but are not lost:

1. Choose the two datasets in the family panel's **X** and **Y** selectors: the session and a comparison dataset, or two comparison datasets. They must be of different species.
2. Click **Score families**. Each family is scored as a module in each dataset.
3. The **family concordance scatter** shows one point per family (x: effect in dataset X, y: effect in dataset Y). Red and blue points are concordant and significant in both datasets; amber diamonds are significant in both but opposite in direction.
4. Click a point, or **plot** in a table row, for per-sample scores in both datasets.
5. **Clear** resets the panel.

**Family dot plot across all datasets.** Under the scatter, type member symbols of either species into **family dot plot** (for example `SERPINB3, S100A7, LCE3`; a prefix such as `LCE3` matches every LCE3 member), choose how many families, and click **Draw**. Leave the box blank for the families significant in the same direction in the most datasets. Every loaded human or mouse dataset with a two-group contrast is drawn, one coloured point per dataset; filled points are significant at the DE FDR within that dataset, open points are not. Each family is scored in each dataset on that species' own members, so raw values are never merged. Click a family for its per-sample scores in every dataset. The plot does not depend on the X and Y chosen above.

### Step 7. Hub neighbourhoods across datasets

1. In the Co-expression tab, under **Network hub genes**, choose the number of variable genes and the soft power β, and click **Find hubs**.
2. Click **neighborhood** in a hub's table row to open its **Hub neighbourhood** card: the hub and its strongest partners, with correlation computed across all samples or within groups.
3. Click **Across datasets**. The Compare datasets tab opens the **Hub across datasets** card. Click **Analyze this neighbourhood** to answer four questions in every comparison dataset:
   1. Is the neighbourhood preserved? WGCNA module-preservation Z-scores, and a permutation p-value against random gene sets matched for differential expression.
   2. Is it a hub in every dataset? The hub's connectivity rank in each dataset.
   3. Does it rewire between conditions? The change in each hub–partner correlation between the two contrast groups (Fisher z test), per dataset and pooled.
   4. How do the members behave in the meta-analysis? Their pooled effects, with and without the dataset in which the hubs were found.
4. **Find consensus hubs** ranks genes by their connectivity across all datasets (robust rank aggregation). **Export tables (CSV)** saves every table.

A neighbourhood found in one dataset should be judged by its preservation in the others, not by its strength in the dataset where it was found.

### Step 8. Transfer a signature

Under **Signature transfer**, choose a module from dataset A (up- or down-regulated DE genes, or a Discovery module) and click **Score**. The module is scored in each other dataset and compared with random signatures of the same size.

---

## 10. Tutorial G: figures, reports, sessions and methods text

### 10.1 Figure exports

- Saving a figure opens a preview first: check it, choose SVG or PNG and the file name, then click **Save** (or **Close** to cancel). Untick **Figure settings → Preview** to save at once.
- Open **Figure settings** in the header. **Export size** sets the shape of exported figures (as shown; single column 4:3 or 1:1; double column 7:4 or wide).
- **Fonts** sets the type size: screen, paper (1.5×) or poster (2×). Margins adjust automatically.
- **Axis labels**, **Label size** and **Long labels** set the tick labels: orientation (auto, horizontal, 45°, vertical), label size, and long-label handling (wrap, full length, or shortened). With wrap or full length, labels are never cut off; the figure margin grows instead.
- The camera icon on any plot saves it as SVG, the recommended format for publication.
- **Colours** (header): the group palette (default, Okabe–Ito and Tol bright, both colour-blind safe, ColorBrewer Dark2, greyscale), the up/down colours, the heatmap scale (RdBu, PuOr, PiYG, BrBG, blue–white–red) and the continuous scale (viridis, cividis, magma, inferno, greyscale). Any single group can be given its own colour. The choice applies to every plot and export and is saved with the session.
- The **Download data** icon on any plot saves an Excel workbook of the plotted values, with a sheet laid out as a GraphPad Prism table (Column, Grouped or XY) so the figure can be redrawn in the same style as your other figures. In Prism, create a table of the type named on the Notes sheet and paste or import the Prism sheet.
- **Summary figure** (header) composes a four-panel overview with a caption.
- In the Gene Explorer, add genes with **+ Add to comparison**, then choose **Export panel figure** for a multi-panel figure with statistics.

### 10.2 Methods text

The **Auto-generated methods text** card on the Overview & QC tab describes exactly what you ran, with citations. It updates as you change settings. Click **Copy to clipboard** and paste it into your manuscript.

### 10.3 Report

**Generate report** (header) produces a single HTML file with results, figures, methods and citations.

### 10.4 Sessions

**Save session** writes a JSON file containing all settings: group assignments, design type, method, thresholds, covariates and Discovery settings. **Load session** restores them. The file records the application version, and a warning appears if you load it into a different version. **Load different data** returns to the start screen.

### 10.5 DESeq2 cross-check

**DESeq2 script (R)** (Differential Expression tab) downloads an R script, pre-filled with your contrast and sample selection, that runs DESeq2 on your original files.

---

## 11. Tutorial H: verifying the statistics in your own browser

### 11.1 Self-test (in the application)

1. Click **Validate statistics** in the header.
2. benchsiDE regenerates fixed synthetic datasets and runs them through the same functions the tabs use; there is no separate test-only code path.
3. The 21 components (filterByExpr, TMM, MDS, the moderated t and its empirical Bayes prior, voom, significant-gene calls, Benjamini–Hochberg FDR, the moderated F, CAMERA and its inter-gene correlation, FRY, over-representation p-values, a regression fit, the power model, and five meta-analysis models and tests) are compared with reference values computed in R (edgeR, limma, metafor, RNASeqPower), scipy and statsmodels. The badge reads **✓ all 21 components pass** when every value is within its tolerance. Click a row to see the compared values side by side.
4. Click **Download evidence (JSON)** to save the record, which includes your browser's engine string.
5. If any component fails, open an issue on GitHub and attach the evidence file.

### 11.2 Timing benchmark

**Validate statistics → Timing benchmark…** runs a fixed analysis on synthetic datasets of increasing size in your browser and reports the time for each step. Close other tabs first. The datasets can be exported to time other tools on the same computer.

### 11.3 Validation on public data (repository)

The repository's validation on public data is not run from the application. It compares the application's results on GEO GSE63310 with R results (70 checks), and is run automatically in Chromium, Firefox and WebKit on every change to the repository (the Actions tab on GitHub shows the result). To run it yourself, see [VALIDATION.md](../VALIDATION.md). The full validation record, with parity on seven public datasets and the false-positive calibration, is in the repository's `validation/` folder.

---

## 12. Checklist before you report results

- [ ] QC reviewed; any excluded sample is justified and stated in the methods.
- [ ] Design type matches the experiment (numeric only for ordered variables).
- [ ] The DE method matches the data (voom or moderated t for counts).
- [ ] Covariates are included where the design requires them, and are not confounded with the groups.
- [ ] Gene-set results state which test was used and what question it answers.
- [ ] Discovery results are reported as hypotheses unless they survive FDR and are replicated.
- [ ] A meta-analysis states the model and test, the number of datasets and the samples in each; with two or three datasets, any signature from the z test is reported as exploratory.
- [ ] The methods text has been copied from the application after the final analysis.
- [ ] The session file and application version are archived with the data.
