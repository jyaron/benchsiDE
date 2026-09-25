# Filtering

Genes with too few counts to be tested reliably are removed before normalization factors and statistics are computed.

## edgeR filterByExpr (default for counts)

A reimplementation of `edgeR::filterByExpr` with default arguments (`min.count = 10`, `min.total.count = 15`, `large.n = 10`, `min.prop = 0.7`). A gene is kept if its CPM exceeds the CPM equivalent of 10 counts at the median library size in at least as many samples as the smallest group (adjusted for large groups), and its total count is at least 15. The filter is design-aware: it uses the group assignment.

The kept gene set is identical to edgeR's on the validation datasets.

## CPM threshold

Keeps genes with CPM at or above a threshold (default 1) in at least a given number of samples (default 2). This is the only filter available for non-count input.

## Which to use

Use filterByExpr for raw counts. The CPM threshold is retained for normalized input and for comparison with older analyses. The filter mode is recorded in the methods text and in session files.
