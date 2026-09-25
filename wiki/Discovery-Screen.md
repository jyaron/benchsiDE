# Discovery screen

The Discovery tab tests many candidate modules for coordinated change. Its purpose is to find groups of related genes whose members each change modestly, so that few pass a per-gene test, but which together show a consistent shift.

## Candidate modules

| Source | Definition |
|---|---|
| Gene families | Genes sharing a symbol stem (for example Skint1 to Skint11) |
| Loaded GMT library | Each set in the library |
| HGNC gene groups | Curated human gene groups; for non-human data, matched by case-normalized symbol (an orthology approximation) |

Modules are tested if their size in the filtered data lies within the size limits (default 3 to 100).

## Module score

For each sample, the module score is the mean z-score of the module's genes, with each gene standardized across the included samples.

## Tests

| Design | Test |
|---|---|
| Two groups | CAMERA (limma), a competitive test that inflates the variance by the estimated inter-gene correlation. Parametric; no permutation floor. |
| Three or more groups | Seeded permutation test against size-matched random modules. The inter-gene correlation is shown for each module as a caution, because random-set nulls understate it. |

FDR is controlled by Benjamini–Hochberg across all candidates. Overlapping modules are positively dependent, a condition under which BH remains valid.

## Result tiers

- **Significant**: FDR at or below 0.05.
- **Hidden module**: significant as a module while fewer than 25% of members are individually significant.
- **Nominal tier**: shown only when no module survives FDR; modules with nominal p at or below 0.01, labelled as nominal. These are ranked hypotheses, not findings.

## Robustness

A leave-one-out note appears on a result card when the direction of the module's effect depends on a single sample.

## Interpretation

Structural gene families (for example DEAD-box helicases) can move together because of global changes in cell composition or proliferation rather than a shared regulatory programme. Functionally coherent families are stronger candidates. Confirm results in an independent dataset before reporting them as findings.

## Reference

Wu D and Smyth GK (2012) Nucleic Acids Res 40:e133.
