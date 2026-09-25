# Tests

`selftest.spec.js` runs in Chromium, Firefox and WebKit through Playwright. It contains:

1. The in-app numerical self-test (16 components compared with embedded edgeR, limma, scipy and statsmodels references), on both the CDN build and the offline build.
2. The demo dataset (GSE63310) analyzed with TMM normalization and voom, under both filter modes, against results pinned from an independent R analysis.
3. An interface test of the **Validate statistics** button.

`validation.spec.js` runs the external validation: the demo dataset is analyzed by the application and compared with the R reference values in `validation/gse63310_reference.json` (33 checks; see `VALIDATION.md`).

Run locally:

```bash
npm install
npx playwright install --with-deps
python3 scripts/build_offline.py
npm test
```

The self-test references were generated with edgeR, limma, scipy and statsmodels on the deterministic synthetic dataset that the self-test regenerates in the browser (seeded pseudo-random generator, integer arithmetic). The package versions are recorded in the self-test evidence file.
