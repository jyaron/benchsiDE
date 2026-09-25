# Input formats

## Expression matrix

| Property | Requirement |
|---|---|
| Delimiter | Tab or comma; detected automatically |
| Orientation | Genes in rows and samples in columns |
| Identifier column | First column; symbols, Ensembl gene IDs (version suffix tolerated) or NCBI/Entrez integers |
| Annotation columns | Detected by content (for example gene name, biotype) and excluded from the numeric matrix |
| Values | Raw counts, CPM/TPM/FPKM, or log-scale values |

### Value types

The review step classifies the values as raw counts, normalized linear values or log₂ values, and the classification can be overridden with **Values are**. The classification determines which methods are available:

| Units | TMM | filterByExpr | voom | Moderated t | Welch t |
|---|---|---|---|---|---|
| Raw counts | yes | yes | yes | yes | yes |
| CPM / TPM / FPKM | no | no | no | yes | yes |
| log₂ | no | no | no | yes | yes |

Non-integer counts (for example estimated counts from Salmon) are accepted on the counts path. Classification is automatic and can be overridden with **Values are**.

## Design file

Tab-separated, in either orientation (samples in rows with factors in columns, or factors in rows with samples in columns). Sample names must match the matrix headers exactly. The first factor is used for grouping unless another is selected; additional factors can be used as covariates. All line-ending conventions are accepted.

## Gene-set files (GMT)

One set per line: set name, description, then member genes, separated by tabs. Members can be symbols or identifiers; matching is case-insensitive. Sets are tested only if their number of members present in the filtered data lies within the size limits (default 5 to 2,000).

## Annotation files

Two columns, tab-separated: identifier and symbol. Common formats (HGNC, MGI, BioMart exports) are recognized automatically. Only tab is used as the delimiter, so descriptive fields that contain commas are read correctly.

## Session files

JSON files written by **Save session**. See [Sessions, reports and methods text](Sessions-Reports-and-Methods-Text.md).
