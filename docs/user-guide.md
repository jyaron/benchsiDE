# benchsiDE User Guide

This guide takes you from an expression matrix to reportable results, one step at a time. Each tutorial names the exact controls you will use. Concepts, statistical definitions and troubleshooting are covered in more depth in the [Knowledgebase](../wiki/Home.md).

Version covered: 0.19.x. The application is at <https://www.benchside.org>; source code and releases are at <https://github.com/jyaron/benchsiDE>.

## Contents

1. [Before you start](#1-before-you-start)
2. [Preparing your input files](#2-preparing-your-input-files)
3. [Tutorial A: your first analysis with the demo dataset](#3-tutorial-a-your-first-analysis-with-the-demo-dataset)
4. [Tutorial B: a two-group experiment with a batch covariate](#4-tutorial-b-a-two-group-experiment-with-a-batch-covariate)
5. [Tutorial C: a time course](#5-tutorial-c-a-time-course)
6. [Tutorial D: gene-set testing](#6-tutorial-d-gene-set-testing)
7. [Tutorial E: the Discovery screen](#7-tutorial-e-the-discovery-screen)
8. [Tutorial F: comparing datasets and species](#8-tutorial-f-comparing-datasets-and-species)
9. [Tutorial G: figures, reports, sessions and methods text](#9-tutorial-g-figures-reports-sessions-and-methods-text)
10. [Tutorial H: verifying the statistics in your own browser](#10-tutorial-h-verifying-the-statistics-in-your-own-browser)
11. [Checklist before you report results](#11-checklist-before-you-report-results)

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
4. Read the result: 2,810 genes up and 3,338 down in LP relative to Basal.
5. Switch between **Volcano** and **MA plot**. Set **label top** to the number of genes to annotate, and choose whether to rank them by **p-value**, **|log₂FC|** or **π score**.
6. Open the **p-value histogram (diagnostic)** and **voom mean–variance trend** panels to check the model assumptions.
7. Click **Export full table (CSV)** to save all genes with log₂FC, t, p and FDR.

### Step 10. Look at individual genes (Gene Explorer tab)

Click any gene in the top-genes table. The Gene Explorer opens with per-group values. Use the toggles to switch between **log₂** and **linear** scale and between **± SEM** and **± SD**. Choose **pairwise vs first group** or **all pairs (sig. only)** to add Holm-corrected pairwise tests to the plot. **← back** returns you to the previous view.

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

1. In the Differential Expression tab, choose **Moderated t (eBayes)**. Covariate adjustment is available for this method only; the other methods show an explicit unadjusted warning.
2. Select the covariate(s) under **adjust for**. A numeric covariate with more than two distinct values (such as RIN) is fitted as a continuous slope; a categorical covariate is fitted as indicator terms.
3. The note under the plot states which covariates were fitted and the residual degrees of freedom.

If a covariate is confounded with the groups (for example, every control in batch 1 and every treated sample in batch 2), benchsiDE refuses to fit the model and explains why. No statistical method can separate the two effects in that case.

### Step 5. Paired designs

For paired samples (the same subject before and after treatment), add a column identifying the subject and select it under **adjust for**. This fits subject as a blocking factor.

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

## 6. Tutorial D: gene-set testing

### Step 1. Choose a library (Enrichment tab)

Under **1 · Gene-set library (GMT)**, choose a built-in library for your species or drop your own `.gmt` file.

### Step 2. Choose the question you are asking

| Method | Question | Input |
|---|---|---|
| **ORA (gene list)** | Are my significant genes over-represented in this set? | The current DE hit list |
| **GSEA (preranked)** | Are the genes of this set concentrated at the top or bottom of the ranking? | All genes ranked by the moderated t-statistic |
| **FRY (self-contained)** | Is the set as a whole differentially expressed? | Per-sample expression of the set's genes |

FRY is usually the most sensitive of the three when replicates are few. FRY finding signal where ORA or GSEA do not is expected, not a contradiction: the tests answer different questions (see [Gene-set testing](../wiki/Gene-Set-Testing.md)).

### Step 3. Run and read the results

Click **Run enrichment**. The table lists each set with its statistic, p-value and FDR. Supporting plots follow:

1. **Ridgeline**: the distribution of log₂FC within each top set.
2. **Barcode plot** and **set-vs-rest ECDF**: click **bc** in a table row, or a ridgeline row, to draw them for that set.
3. **Jaccard heatmap**: overlap between the top sets. Click a cell to list shared genes.
4. **Leading-edge matrix** (GSEA): which genes drive which sets.

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

## 8. Tutorial F: comparing datasets and species

The Compare datasets tab relates your session (dataset A) to up to eight other datasets. Each dataset is filtered, normalized and tested on its own; raw values are never merged across datasets.

### Step 1. Add a dataset

1. Set **species of next dataset** (same as session, mouse or human).
2. Click **Add dataset (matrix)** and select the file.
3. Optionally click **Design for last added** and select its design file.

The slot reports how many genes matched. If the count is low, the message names the identifier types on both sides and the remedy.

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

Click **Compute across all datasets** for the correlation matrix and consensus table. Type a gene symbol in the **forest plot** field and click **Draw** for its effect in every dataset.

### Step 5. Score homology families (cross-species)

Many genes have no one-to-one ortholog (for example the mouse Serpinb3a to Serpinb3d genes and human SERPINB3 and SERPINB4). These families are excluded from gene-level matching but are not lost:

1. Choose the comparison dataset in the family panel's **Dataset** selector.
2. Click **Score families**. Each family is scored as a module in each dataset.
3. The **family concordance scatter** shows one point per family (x: effect in your session, y: effect in the other dataset). Red and blue points are concordant and significant in both datasets; amber diamonds are significant in both but opposite in direction.
4. Click a point, or **plot** in a table row, for per-sample scores in both datasets.
5. **Clear** resets the panel.

### Step 6. Transfer a signature

Under **Signature transfer**, choose a module from dataset A (up- or down-regulated DE genes, or a Discovery module) and click **Score**. The module is scored in each other dataset and compared with random signatures of the same size.

---

## 9. Tutorial G: figures, reports, sessions and methods text

### 9.1 Figure exports

- Set the export shape with the **export** selector (as shown; single column 4:3 or 1:1; double column 7:4 or wide).
- Set the type size with the font selector: **screen fonts**, **paper fonts (1.5×)** or **poster fonts (2×)**. Margins adjust automatically, and labels too long for the margin are shortened in the middle (the full name remains in the hover text).
- The camera icon on any plot saves it as SVG, the recommended format for publication.
- **Summary figure** (header) composes a four-panel overview with a caption.
- In the Gene Explorer, add genes with **+ Add to comparison**, then choose **Export panel figure** for a multi-panel figure with statistics.

### 9.2 Methods text

The **Auto-generated methods text** card on the Overview & QC tab describes exactly what you ran, with citations. It updates as you change settings. Click **Copy to clipboard** and paste it into your manuscript.

### 9.3 Report

**Generate report** (header) produces a single HTML file with results, figures, methods and citations.

### 9.4 Sessions

**Save session** writes a JSON file containing all settings: group assignments, design type, method, thresholds, covariates and Discovery settings. **Load session** restores them. The file records the application version, and a warning appears if you load it into a different version. **Load different data** returns to the start screen.

### 9.5 DESeq2 cross-check

**DESeq2 script (R)** (Differential Expression tab) downloads an R script, pre-filled with your contrast and sample selection, that runs DESeq2 on your original files.

---

## 10. Tutorial H: verifying the statistics in your own browser

1. Click **Validate statistics** in the header.
2. benchsiDE regenerates a fixed synthetic dataset and runs it through the same functions the tabs use.
3. Sixteen components are compared with reference values computed in R (edgeR, limma), scipy and statsmodels. The panel shows every compared value side by side with its tolerance.
4. Click **Download evidence (JSON)** to save the record, which includes your browser's engine string.
5. If any component fails, open an issue on GitHub and attach the evidence file.

---

## 11. Checklist before you report results

- [ ] QC reviewed; any excluded sample is justified and stated in the methods.
- [ ] Design type matches the experiment (numeric only for ordered variables).
- [ ] The DE method matches the data (voom or moderated t for counts).
- [ ] Covariates are included where the design requires them, and are not confounded with the groups.
- [ ] Gene-set results state which test was used and what question it answers.
- [ ] Discovery results are reported as hypotheses unless they survive FDR and are replicated.
- [ ] The methods text has been copied from the application after the final analysis.
- [ ] The session file and application version are archived with the data.
