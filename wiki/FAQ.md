# FAQ

**Is my data uploaded?**
No. All computation happens in the browser. See [Privacy and architecture](Privacy-and-Architecture.md).

**Are the results the same as R?**
For the implemented methods, yes, to the precision stated in [VALIDATION.md](../VALIDATION.md). The self-test lets you confirm this in your own browser.

**Why is DESeq2 not included?**
DESeq2 is an iterative estimation procedure (negative-binomial GLMs, dispersion shrinkage, outlier handling). A browser version that did not reproduce Bioconductor's results exactly would be misleading. The application exports a DESeq2 script for the current contrast instead.

**Should I use voom or the moderated t?**
For raw counts, voom. The moderated t on log₂-CPM is appropriate for normalized input and gives similar results when library sizes are similar.

**Why do ORA, GSEA and FRY give different results?**
They test different hypotheses. See [Gene-set testing](Gene-Set-Testing.md).

**Can I analyze a paired design?**
Yes. Add the subject identifier to the design file and select it as a covariate.

**Can I compare mouse and human data?**
Yes. See [Cross-dataset comparison](Cross-Dataset-Comparison.md) and [Orthologs and homology families](Orthologs-and-Homology-Families.md).

**Does it work offline?**
The offline build (`benchside-offline.html`) needs no network connection.

**How do I cite it?**
See [Citing](Citing.md).
