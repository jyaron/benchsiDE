# Sessions, reports and methods text

## Sessions

**Save session** writes a JSON file containing the analysis state: sample inclusion, group assignments and order, design type, normalization, filter, DE method, contrast and thresholds, covariates, font scale, annotation choice and Discovery settings. The file does not contain the expression data. To restore, load the same data files, then **Load session**.

The file records the application version. Loading a session saved by a different version displays a warning, because defaults or methods may have changed.

## Methods text

The Overview & QC tab contains **Auto-generated methods text**, which describes the analysis as performed: input dimensions, filter, normalization, methods, thresholds, covariates, gene-set tests and cross-dataset procedures, each with its citation. It updates when settings change. Copy it after the final analysis.

## Report

**Generate report** writes a single self-contained HTML file with the settings, sample table, DE summary and top genes, figures, Discovery results, methods text and references.

## DESeq2 companion script

**DESeq2 script (R)** writes an R script for the current contrast and sample selection, for users who want a DESeq2 analysis of the same data.
