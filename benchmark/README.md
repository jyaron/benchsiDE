# Timing benchmark, benchsiDE 1.0.0-rc.1

Records exported by **Validate statistics → Timing benchmark** on 4 October 2026, one computer, in Chrome 154, Firefox 157 and Safari 18.6 (versions as reported by each browser). Synthetic data: 20,000 genes, two groups, 6, 24, 96 and 192 samples; one warm-up and five measured repetitions per size.

- `benchside_benchmark_*.json`: every repetition, warm-up included, in milliseconds, with the browser's environment report.
- `timing_benchmark_medians.csv`: median, minimum and maximum per browser, size and step.
- `timing_benchmark.pdf`: run time against the number of samples (median; shaded: range of the five repetitions).

GSEA permutes genes with 6 samples (too few distinct sample labelings) and sample labels from 24 samples on. The plot-drawing step includes the browser's rendering of the next frame.
