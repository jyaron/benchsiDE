# Co-expression and hub genes

## Co-expression

Enter a query gene. benchsiDE computes the Pearson correlation between the query and every other gene across the included samples, on log₂ values, and lists the most positively and most negatively correlated genes (30, 100 or 300 per direction). Each row gives r, the slope of the ordinary least-squares regression of the gene on the query gene (log₂ units per log₂ unit), the p-value for r = 0 (t test with n − 2 degrees of freedom), and the Benjamini–Hochberg FDR across all genes tested against the query gene. With about 20,000 genes, use the FDR rather than p to decide which correlations to report. Click a column header to sort; **Export all (CSV)** writes every gene. Click a partner to plot the pair with the least-squares fit and its 95% confidence band. r, slope and p are validated against R `cor`, `lm` and `cor.test` (VALIDATION.md).

With few samples, correlations are imprecise: at n = 8, a correlation of 0.7 has a 95% confidence interval (Fisher z) of approximately −0.01 to 0.94. Treat co-expression lists as hypotheses.


To look at a pair in more detail, use **Open in Gene Explorer** below the pair statistics: it shows either gene by group, with its tests, or both genes side by side. **Back** in Gene Explorer returns to the co-expression view.

## Hub genes

**Find hubs** computes soft connectivity for the most variable genes (default 2,000): the sum over all other genes of the adjacency |r|^β (unsigned) or ((r + 1)/2)^β (signed), with soft power β (default 6). This follows the WGCNA connectivity definition. Genes with the highest connectivity are listed, with their DE statistics where available.

The **scale-free diagnostic** plots log frequency against log connectivity. A high R² with a negative slope indicates an approximately scale-free network at the chosen β. A positive slope means the network is not scale-free at this β, which is common for small variable-gene subsets. The diagnostic is shown for information and is not used to choose β.


## Hub neighbourhood

Click **neighborhood** next to a hub in the hub table. The panel shows the hub's most strongly correlated genes (10, 15 or 25) among the genes screened for connectivity, in four views:

| View | Question it answers |
|---|---|
| Ego network | What is the hub connected to, and do its partners form one module? Edge width is \|r\| and edge colour the sign; partner–partner edges above the chosen \|r\| are drawn faintly. Nodes are coloured by log₂FC and outlined when differentially expressed. |
| Heatmap | Do the members move together, and in which condition? Rows are z-scores (hub on top), samples are grouped by condition. |
| Neighbourhood score | How large is the effect? Each sample's mean z of the members, negatively correlated members sign-reversed, by group. Descriptive: the members were chosen from the same samples, so no p-value is given. |
| Gene-set enrichment | What is it? Hypergeometric test against the screened genes annotated in the loaded library, BH-FDR. |

**Correlation mode.** *Across all samples* is ordinary Pearson r. In a two-group design it is dominated by the group difference: any two genes that are both up in one group correlate whether or not they are co-regulated. *Within groups* uses the partial correlation given group (the correlation of the residuals after removing each group's mean), which keeps only co-variation within conditions (Baba, Shibata & Sibuya 2004; cf. Parsana et al. 2019). Its p-values use n − G − 1 degrees of freedom. A neighbourhood that holds in this mode reflects coordination beyond the shared group effect.

**Replication.** *Test in comparison datasets* scores the positively correlated members in every comparison dataset with the signature-transfer test (calibrated against size-matched random gene sets).

**Exports.** Cytoscape node and edge tables (CSV), the neighbourhood as a GMT gene set, and the Prism/Excel data of each plot.

## Network view

After **Find hubs**, **Draw network** shows the most connected genes as a network.

| Setting | Effect |
|---|---|
| top hubs | number of hubs shown (10, 25, 50) |
| partners per hub | each hub's most strongly correlated genes added as nodes (0 to 10) |
| edges \|r\| ≥ | minimum absolute correlation for an edge |
| strongest per gene | each gene keeps only its strongest edges (2 to 10, or all); an edge is drawn if it is among the strongest of either gene. Hubs are correlated with one another by construction, so without this limit most pairs are connected. |
| layout | **circular**: genes ordered by average-linkage clustering on 1 − \|r\|, edges drawn as curves; **force-directed**: strongly correlated genes pulled together (Fruchterman–Reingold, seeded, reproducible) |
| colour | log₂FC of the current contrast (outline = FDR ≤ 0.05), or connectivity |

Node size is connectivity. Edges are red for positive and blue for negative correlation, darker and thicker when stronger. Click a node to open the gene. **Export network (CSV)** writes an edge table (source, target, r) and a node table, which Cytoscape imports directly. Correlation edges describe co-variation across samples, not regulation.

## Reference

- Fruchterman TMJ, Reingold EM. Graph drawing by force-directed placement. Softw Pract Exp 21:1129–1164, 1991.

Zhang B and Horvath S (2005) Stat Appl Genet Mol Biol 4:Article17.
