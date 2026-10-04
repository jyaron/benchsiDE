# Tests

`selftest.spec.js` runs in Chromium, Firefox and WebKit through Playwright. It contains:

1. The in-app numerical self-test (21 components compared with embedded edgeR, limma, metafor, scipy and statsmodels references), on both the CDN build and the offline build.
2. The demo dataset (GSE63310) analyzed with TMM normalization and voom, under both filter modes, against results pinned from an independent R analysis.
3. An interface test of the **Validate statistics** button.

`validation.spec.js` runs the external validation: the demo dataset is analyzed by the application and compared with the R reference values in `validation/gse63310_reference.json` (59 checks; see `VALIDATION.md`).

Run locally:

```bash
npm install
npx playwright install --with-deps
python3 scripts/build_offline.py
npm test
```

The self-test references were generated with edgeR, limma, scipy and statsmodels on the deterministic synthetic dataset that the self-test regenerates in the browser (seeded pseudo-random generator, integer arithmetic). The package versions are recorded in the self-test evidence file.

## Timing benchmark

For comparisons with other tools on a user's own files, the line above the dashboard after **Analyze** gives the time to read and parse the dropped files and from **Analyze** to the first display, excluding the time spent on the review step.


`npm run bench` runs `tests/bench/benchmark.spec.js` in Chromium, Firefox and WebKit, one browser at a time, on the offline build. It times a fixed analysis (load and normalize, moderated t, voom, FRY, GSEA with 1,000 permutations, hub screen, plot rendering) on deterministic synthetic datasets of 20,000 genes and 6, 24, 96 and 192 samples: one warm-up and 5 timed repetitions each. Results are written to `bench-results/`: a JSON file per browser and a `_runs.csv` with every run of every step (warm-up flagged), in milliseconds as measured. `BENCH_REPS` and `BENCH_SIZES` (for example `BENCH_SIZES=6,24`) change the defaults. The same benchmark runs in any browser from **Validate statistics → Timing benchmark…**, which also exports the datasets for use in other tools. Timings on shared CI runners are not reproducible and are not reported.
