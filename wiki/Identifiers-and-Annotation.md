# Identifiers and annotation

## Built-in annotation

For mouse, human and rat, benchsiDE embeds identifier-to-symbol tables:

| Species | Source | Identifier types |
|---|---|---|
| Mouse | MGI (MRK_ENSEMBL.rpt and MGI_EntrezGene.rpt) | Ensembl, Entrez |
| Human | HGNC complete set | Ensembl, Entrez |
| Rat | RGD GENES_RAT | Ensembl, Entrez |

Select the species on the start screen and click **Use built-in annotation**. The source and its citation are added to the methods text. The tables are subsets of the source files, redistributed under the providers' licenses (MGI and RGD: CC BY 4.0; HGNC: CC0 1.0).

## Unmapped identifiers

If some identifiers are not found, the start screen reports the number and offers **Download unmapped IDs**. Typical causes are retired identifiers, non-coding features absent from the annotation release, or identifiers from another species.

## Custom annotation

Supply a two-column file (identifier, symbol) for species or identifier types not covered by the built-in tables. The file is cited as a user-supplied mapping.

## Gene links

The Gene Explorer links each gene to external resources appropriate to its species (Ensembl, NCBI Gene, MGI, RGD, GeneCards). Links open in a new tab; no expression data is sent.

## Identifiers in comparisons

Gene matching between datasets is described in [Cross-dataset comparison](Cross-Dataset-Comparison.md) and [Orthologs and homology families](Orthologs-and-Homology-Families.md).
