# Multiple contrasts

## All pairwise contrasts

**Run all pairs (current method)** computes every pairwise contrast with the current method and thresholds. FDR is controlled within each contrast; no joint model across contrasts is fitted.

Outputs:

- A summary table of up and down hits per contrast.
- An UpSet plot of hit-set intersections. Click a bar for its genes.
- A volcano grid with one panel per contrast.

## Venn diagrams

Select two or three contrasts and a direction (either, up only, down only). Region counts are exact set partitions. Circle areas are not proportional to counts, because exact area-proportional three-set diagrams are generally impossible; the plot states this. Click a region to list its genes; **Export regions CSV** writes every gene with its membership pattern.

## Comparing two contrasts

**Compare two contrasts** plots the log₂FC of one contrast against another, with a regression line, and lists genes shared by both contrasts or exclusive to one.

## Enrichment across contrasts

See [Gene-set testing](Gene-Set-Testing.md).

## Any-group test

For the question "does this gene differ among any groups", use the genome-wide moderated F test (Group Patterns tab), which uses all samples in a single test rather than many pairwise tests.
