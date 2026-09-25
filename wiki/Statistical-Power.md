# Statistical power

## Small samples

With two or three replicates per group, genome-wide FDR control can be unattainable even for large effects: the smallest p-value the data can produce may be larger than the BH threshold for a single discovery. A contrast with few significant genes is therefore not evidence of no effect.

Useful instruments in this situation:

- voom or the moderated t, which borrow variance information across genes.
- The stability curves on the Differential Expression tab.
- Gene-set tests, particularly FRY, which aggregate evidence across genes.
- The moderated F test, which uses all samples.

## Power and design guidance

**Compute from this pilot** (Overview & QC tab) estimates the replicates needed per group for a chosen fold change, using the sequencing depth and biological coefficient of variation observed in the current data. The calculation follows the model of the RNASeqPower package (Hart SN et al. 2013, J Comput Biol 20:970).

Treat the estimate as approximate: it assumes independent genes and a negative-binomial count model, and it uses the pilot's dispersion as if it were known.
