# benchsiDE

Validated, client-side RNA-seq differential expression and discovery analysis in a single HTML file.

- Application: <https://www.benchside.org>
- User guide: [docs/user-guide.md](docs/user-guide.md) (also at <https://www.benchside.org/guide.html>)
- Knowledgebase: [wiki/Home.md](wiki/Home.md) (also on the repository wiki)
- Validation record: [VALIDATION.md](VALIDATION.md)

benchsiDE performs normalization, differential expression, gene-set testing, module discovery and cross-dataset comparison in the browser. Expression data is read locally and is not transmitted. Each statistical routine is validated against its reference implementation, and the validation suite can be rerun in any browser from the application header.

## Properties

| Property | Detail |
|---|---|
| Local computation | No server component and no upload. The offline build makes no network requests. |
| Validated statistics | TMM normalization and filterByExpr (edgeR); moderated t, moderated F, voom, covariate models, CAMERA and FRY (limma); RRHO; hypergeometric, KS and OLS tests (scipy, statsmodels). Checked on a synthetic dataset and on the public demo dataset GSE63310; see [VALIDATION.md](VALIDATION.md). |
| In-browser verification | **Validate statistics** regenerates a deterministic dataset and compares 16 components with embedded reference values, listing every compared value and its tolerance. |
| Reproducibility | Seeded randomness throughout; version-stamped session files; self-contained HTML reports; methods text generated from the settings used, with citations; DESeq2 companion script export. |

## Quick start

1. Open <https://www.benchside.org>, or download `index.html` or `benchside-offline.html` from the latest [release](https://github.com/jyaron/benchsiDE/releases).
2. Drop an expression matrix (genes in rows, samples in columns) onto card 1, and optionally a design file onto card 2.
3. Review the detected settings and click **Analyze**.

The included demo dataset (`demo/GSE63310_counts.tsv`, `demo/GSE63310_design.tsv`; mouse mammary epithelium, GEO GSE63310) is the dataset of the published limma-voom workflow. Tutorial A of the user guide walks through it with the expected results at each step.

## Features

| Area | Content |
|---|---|
| Quality control | Library size, detected genes, PCA with loadings, MDS, dendrogram, sample correlation, RLE, expression density, library complexity, sex check, biotype composition, power estimation |
| Normalization and filtering | TMM or library-size CPM; filterByExpr or CPM threshold; display-only batch removal |
| Differential expression | Moderated t with covariates (categorical, continuous, blocking), voom, Welch t; moderated F; all-pairs contrasts, UpSet and Venn; diagnostic plots |
| Gene-set testing | ORA, preranked GSEA, FRY; built-in MSigDB Hallmark, GO (BP, CC, MF) and Reactome for human and mouse, or any GMT file; enrichment across all contrasts |
| Discovery | Module screen over gene families, HGNC gene groups and gene-set libraries; CAMERA (two groups) or seeded permutation (more groups) |
| Cross-dataset | Up to eight datasets; symbol, mouse–human ortholog or custom matching; fold-change and t concordance, RRHO, forest plots, signature transfer, homology-family scoring |
| Output | SVG export at screen, paper and poster type sizes; composite figures; HTML report; session files; methods text |

## Repository layout

| Path | Content |
|---|---|
| `index.html` | The application (loads Plotly.js from its CDN) |
| `scripts/build_offline.py` | Builds `dist/benchside-offline.html` with Plotly.js embedded (checksum-pinned) |
| `scripts/build_guide.py` | Builds `dist/benchside-guide.html` from `docs/user-guide.md` |
| `scripts/publish_wiki.sh` | Copies `wiki/` to the repository wiki |
| `docs/` | User guide |
| `wiki/` | Knowledgebase source |
| `demo/` | Public demo dataset |
| `validation/` | R reference script and comparison runner for the external validation |
| `tests/` | Playwright test suite (Chromium, Firefox, WebKit) |

## Running the tests

```bash
npm install
npx playwright install --with-deps
python3 scripts/build_offline.py
npm test
```

## Citation

See `CITATION.cff`. A manuscript is in preparation.

## License

The benchsiDE code is released under the MIT license. Embedded third-party resources retain their own terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
