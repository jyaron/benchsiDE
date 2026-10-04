# Sample traits: parity and calibration

Public data: GSE186063, GSE83645, GSE54456 (NCBI-generated GEO count matrices). Traits in the calibration are random values or shuffles of the recorded age; planted effects were added to healthy-skin samples. No dataset with a recorded PASI per RNA-seq sample was available, so no real-severity positive control was run.

Parity with limma 3.66 on GSE186063 (Age; lesion and non-lesion skin; patient pairs as blocks), all genes: duplicateCorrelation consensus correlation to 6.7e-16; moderated t to 4.6e-14 and p to 5.9e-10 (relative) for three models (block means with group removed within blocks; ~ Age + group + Sex; ~ rank(Age)); block means computed independently in R agree to 7.1e-15; CAMERA on the trait coefficient (50 Hallmark sets) p to 3.7e-8 (relative) with and without the mixed-model transformation; module eigengenes against WGCNA moduleEigengenes to 4.6e-13 and their trait t against lm to 1e-10.

| Check | Expected | Result | Consequence |
|---|---|---|---|
| Random subject-level trait, GSE83645 lesional skin (20 biopsies from 5 patients), 40 replicates, samples treated as independent | 5% of genes at p < 0.05 | 25.0%; any gene at FDR ≤ 0.05 in 90% of replicates | repeated samples must be declared |
| Same, mixed model (duplicateCorrelation) | 5% | 8.4%; any FDR hit in 70% | not used for subject-level traits: a consensus correlation is too low for many genes with only 5 subjects |
| Same, block means (default for subject-level traits) | 5% | 6.3%; any FDR hit in 5% | calibrated |
| Random sample-level trait, GSE83645 lesional, 40 replicates: independent / mixed model | 5% | 4.7% / 5.0%; any FDR hit 5% / 2% | calibrated |
| GSE186063 Age shuffled across patients (lesion + non-lesion, group-adjusted), 40 replicates: independent / mixed model / block means | any FDR hit in ≤ 5% (FDR control under the complete null) | 18% / 8% / 2% | block means is the default |
| Random trait, GSE54456 (171 samples, group-adjusted), 25 replicates: genes / CAMERA (Hallmark) | 5% at p < 0.05 | 4.2% / 4.0% | calibrated |
| Random trait, module eigengene test (300 replicates; random 20-gene modules / hub neighbourhoods) | 5% | 6.0% / 6.3% | calibrated |
| Planted trait effects in healthy skin (GSE54456), 200 genes, 5 replicates each | FDP ≤ 5%; unbiased slopes | n = 10, 0.6 log₂/SD: 41% recovered, FDP 1.2%; n = 20, 0.4 log₂/SD: 57% recovered, FDP 2.4%; n = 40, 0.3 log₂/SD: 65% recovered, FDP 0.4%; n = 81, 0.2 log₂/SD: 66% recovered, FDP 0.6% | mean slope bias < 0.005 log₂/SD |
