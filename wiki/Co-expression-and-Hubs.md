# Co-expression and hub genes

## Co-expression

Enter a query gene. benchsiDE computes the Pearson correlation between the query and every other gene across the included samples, on log₂ values, and lists the 30 most positively and 30 most negatively correlated genes. Click a partner to plot the pair with an ordinary least-squares fit and its 95% confidence band.

With few samples, correlations are imprecise: at n = 8, a correlation of 0.7 has a 95% confidence interval (Fisher z) of approximately −0.01 to 0.94. Treat co-expression lists as hypotheses.

## Hub genes

**Find hubs** computes soft connectivity for the most variable genes (default 2,000): the sum over all other genes of the adjacency |r|^β (unsigned) or ((r + 1)/2)^β (signed), with soft power β (default 6). This follows the WGCNA connectivity definition. Genes with the highest connectivity are listed, with their DE statistics where available.

The **scale-free diagnostic** plots log frequency against log connectivity. A high R² with a negative slope indicates an approximately scale-free network at the chosen β. A positive slope means the network is not scale-free at this β, which is common for small variable-gene subsets. The diagnostic is shown for information and is not used to choose β.

## Reference

Zhang B and Horvath S (2005) Stat Appl Genet Mol Biol 4:Article17.
