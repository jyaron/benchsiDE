# Provenance scan of the validation check tables

Every call that recorded a check (`addcheck` / `addc`) in the five validation sessions was parsed. Thirty-two calls gave both the count (`n_compared`) and the deviation as literal numbers; they are listed in `check_provenance_literal_calls.csv` with their status in the saved tables.

| Group | Calls | Finding | Action |
|---|---|---|---|
| Pass flag computed from the application's output in the same call (input detection, design parsing, QC flags, minimum-size rules, p-value histogram notes, comparison and signature-transfer limits) | 14 | Count = number of cases; the pass flag is an expression on the observed output (for example the detected input scale or the groups formed). | None needed. |
| Defect records (deviation ≥ 1, pass false) | 14 | Observed failures, with the observed output in the note. | Covered by the freeze regression. |
| Sex-chromosome genes absent from the shared signature (VAL_XDS) | 1 | The first entry had literal values. The same session replaced it with four computed entries (one per meta-analysis configuration); only the computed entries are in the saved table. | Recomputed on the freeze build (engine check `cross_dataset`): none of the 10 genes is in the default REML/Hartung-Knapp signature, all 10 were tested, and no mitochondrial or ribosomal-protein gene is in it. |
| Power table cells (VAL_UNT) | 1 | Literal pass that refers to the separate computed RNASeqPower grid checks in the same table, which pass. | None needed; cross-referenced. |
| Discovery heatmap-export caption names CAMERA (VAL_ENR) | 1 | Literal pass based on inspection. | Recomputed on the freeze build (engine check `discovery_text`): both Discovery figure captions name the test. |
| Discovery CAMERA with DE covariates (VAL_ENR) | 1 | Literal pass recording that p-values were identical with and without the covariate, reported as a UI issue (no notice). | The notice was added in the audit; engine check `discovery_text` confirms it appears when covariates are selected. |

No check in the saved tables now reports a pass without a computation behind it or a recomputation on the freeze build.
