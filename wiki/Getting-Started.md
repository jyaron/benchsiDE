# Getting started

## Three ways to run the application

| Build | Network use | Typical use |
|---|---|---|
| Hosted (<https://www.benchside.org>) | Loads the page and the plotting library once | Everyday analysis |
| `index.html` from a release | Loads the plotting library from its CDN | Version-pinned, archived analyses |
| `benchside-offline.html` from a release | None | Air-gapped systems, controlled-access data |

All three builds contain identical application code. The offline build additionally embeds Plotly.js 2.32.0.

## Minimum inputs

An expression matrix with genes in rows and samples in columns. A design file is optional; without one, groups are inferred from sample names and can be edited before analysis. See [Input formats](Input-Formats.md).

## The analysis flow

1. Load the matrix, and optionally the design and species annotation.
2. Review the detected settings: design type, normalization, input units and filter.
3. Click **Analyze**. Every tab is computed from the included samples.
4. Toggle samples in the sample bar to include or exclude them; all results recompute.

For a guided first run with the included demo dataset, follow Tutorial A in the [User Guide](../docs/user-guide.md).
