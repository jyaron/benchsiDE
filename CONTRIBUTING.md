# Contributing

benchsiDE is a single-file application: `index.html` contains all application code and embedded resources. This design keeps the application auditable, portable and free of a build step for users. Contributions must preserve it.

## Requirements for changes

1. **Statistical changes require reference validation.** A change to a statistical routine must include parity evidence against the reference implementation (edgeR or limma in R, or scipy or statsmodels), and the embedded self-test references must be regenerated if they are affected. Changes to numerical output without updated validation are not merged.
2. **Tests must pass in all three browser engines** (`npm test` runs Chromium, Firefox and WebKit).
3. **No new network dependencies.** The offline build must remain functional without network access.
4. **Documentation moves with the code.** Update `docs/user-guide.md`, the relevant `wiki/` pages and the methods text (including citations) in the same pull request as the feature.
5. **Exports carry no interface text.** Exported figures must not contain on-screen instructions or hyperlinks.

## Development

Open `index.html` directly; no build step is required. `scripts/build_offline.py` and `scripts/build_guide.py` produce the files in `dist/`. The Playwright suite in `tests/` is the required check for pull requests.

## Reporting a numerical discrepancy

Run **Validate statistics**, download the evidence file, and attach it to an issue together with the browser and operating system.
