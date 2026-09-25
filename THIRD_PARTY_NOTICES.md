# Third-party material embedded in benchsiDE

benchsiDE application code is released under the MIT License (`LICENSE`). The application file embeds the third-party material listed below. The same notice is included in `index.html` itself (a comment at the top of the file, and the **licenses** panel in the application header), so it accompanies every copy of the file.

| Resource | Use in benchsiDE | License | Modification | Terms verified |
|---|---|---|---|---|
| Plotly.js 2.32.0 | Plotting. Loaded from the Plotly CDN by `index.html`; embedded in the offline build. | [MIT](https://github.com/plotly/plotly.js/blob/master/LICENSE) | Unmodified | Bundle header |
| MSigDB v2024.1: Hallmark, GO BP, GO CC, GO MF and Reactome collections (human and mouse) | Built-in gene-set libraries | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Re-encoded; set membership unchanged | Gene-set pages on gsea-msigdb.org. Collections with additional terms (KEGG, BioCarta, Database of Cell Signaling) are not embedded. |
| MGI (The Jackson Laboratory): MRK_ENSEMBL.rpt, MGI_EntrezGene.rpt, HOM_MouseHumanSequence.rpt | Mouse identifier-to-symbol mapping, accession links, mouse–human orthologs, homology families | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Subset: selected columns; re-encoded | <https://www.informatics.jax.org/mgihome/other/copyright.shtml> (2026-09-24) |
| RGD (Medical College of Wisconsin): GENES_RAT | Rat identifier-to-symbol mapping, accession links | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Subset: selected columns; re-encoded | <https://rgd.mcw.edu/wg/disclaimer/> (2026-09-24) |
| HGNC (HUGO Gene Nomenclature Committee at the University of Cambridge): complete set and gene groups | Human identifier-to-symbol mapping; curated discovery modules | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) | Subset: selected columns; re-encoded | <https://www.genenames.org/about/license/> (2026-09-24) |
| GEO GSE63310 (demo dataset, `demo/`) | Demonstration and tests; not embedded in the application | Public GEO data | Unmodified counts assembled per the published workflow | Cite the original study |

## Attribution

Each resource is credited in the application: in the annotation and library labels, in the **licenses** panel, and in the auto-generated methods text, which includes the citation of every resource used in an analysis. HGNC does not require attribution; it is given as requested by HGNC.

## Warranty

The MGI, RGD and HGNC data are provided by their sources "as is", without warranty of any kind.

## Statistical methods

The statistical algorithms are independent JavaScript implementations validated against their reference implementations. No source code from the reference packages is included.
