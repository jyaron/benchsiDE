# Reference values for the engine check "covariates_treat" (benchsiDE 1.0.0-rc.2).
# Run from the repository root after tests/engine/fetch_data.py:  Rscript validation/rc2/make_engine_reference.R
# GSE63310: the lane column is the sequencing lane used by Law et al. (F1000Research 5:1408, 2016); the rin column
# holds synthetic values that exist only to exercise a continuous covariate.
source("validation/rc2/de_reference.R")
cfg <- fromJSON("validation/rc2/engine_configs.json", simplifyVector = FALSE)
td <- tempfile(); dir.create(td)
cnt1 <- read.delim("demo/GSE63310_counts.tsv", row.names = 1, check.names = FALSE)
d1 <- read.delim("tests/engine/designs/GSE63310_lane_rin_design.tsv", stringsAsFactors = FALSE); d1 <- d1[match(colnames(cnt1), d1$sample), ]
cnt2 <- read.delim("tests/engine/data/GSE186063.tsv", row.names = 1, check.names = FALSE)
d2 <- read.delim("tests/engine/designs/GSE186063_design_long.tsv", stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN")); d2 <- d2[match(colnames(cnt2), d2$sample), ]
de_reference(cnt1, d1, "celltype", cfg$GSE63310, file.path(td, "a"))
de_reference(cnt2, d2, "Type", cfg$GSE186063, file.path(td, "b"))
out <- list()
for (ds in c("GSE63310", "GSE186063")) {
  pre <- if (ds == "GSE63310") "a" else "b"
  out[[ds]] <- lapply(seq_along(cfg[[ds]]), function(i) { r <- read.csv(sprintf("%s/%s_%d.csv", td, pre, i - 1), colClasses = c(id = "character"))
    list(cfg = cfg[[ds]][[i]], ids = r$id, t = signif(r$t, 15), p = signif(r$p, 15)) })
}
con <- gzfile("tests/engine/de_rc2_reference.json.gz", "w"); writeLines(toJSON(out, digits = NA, auto_unbox = TRUE), con); close(con)
cat("written\n")
