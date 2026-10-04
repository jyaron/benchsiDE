# Engine checks

These checks run the application's script in Apple's JavaScriptCore shell (`jsc`, present on every Mac) with a minimal stub page and a stub Plotly that records every plot specification. They test computation, the generated plot specifications and the export files. They do not test rendering, layout or browser behaviour; those are covered by the Playwright tests in `tests/`.

1. `python3 tests/engine/fetch_data.py` downloads the public GEO count matrices listed in `datasets.json` into `tests/engine/data/` (not committed) and verifies their checksums.
2. `python3 tests/engine/run_checks.py` runs every check; `--only NAME` runs one. Each line reports PASS or FAIL with the values compared.

Design files in `designs/` were built from the GEO sample annotations; their provenance is described in `datasets.json`.
