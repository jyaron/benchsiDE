# Reference values for benchsiDE interaction contrasts on GSE143688 (all 74 samples).
# Run from the repository root after tests/engine/fetch_data.py:  Rscript validation/freeze/interaction_reference.R
suppressMessages({ library(edgeR); library(limma); library(jsonlite) })
cnt <- as.matrix(read.delim("tests/engine/data/GSE143688_all.tsv", row.names = 1, check.names = FALSE))
des <- read.csv("tests/engine/designs/GSE143688_all_design.csv", header = FALSE, row.names = 1, check.names = FALSE)
des <- as.data.frame(t(des), stringsAsFactors = FALSE); rownames(des) <- colnames(cnt)
cell <- factor(des$Cell); day <- factor(des$Day)
y <- DGEList(cnt, group = cell); keep <- filterByExpr(y, group = cell)
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE]); lc <- cpm(y, log = TRUE, prior.count = 2)
out <- list(limma = as.character(packageVersion("limma")), edgeR = as.character(packageVersion("edgeR")), genes = rownames(y))
# 1. moderated t (limma-trend), all eight groups, adjusted for Day: ~ 0 + cell + day with a contrast
d1 <- model.matrix(~ 0 + cell + day); colnames(d1) <- sub("^cell", "", colnames(d1))
con <- makeContrasts(contrasts = "(Aldara_KO-Control_KO)-(Aldara_WT-Control_WT)", levels = d1)
f1 <- eBayes(contrasts.fit(lmFit(lc, d1), con), trend = TRUE)
t1 <- topTable(f1, number = Inf, sort.by = "none")
out$mod_all_day <- list(logFC = signif(t1$logFC, 12), t = signif(t1$t, 12), p = signif(t1$P.Value, 12), df0 = f1$df.prior[1])
# 2. voom, the four groups only, adjusted for Day. limma's contrasts.fit is approximate when precision weights are
#    combined with a non-orthogonal design (?contrasts.fit), so the interaction is fitted as a coefficient instead.
sel <- cell %in% c("Control_WT", "Aldara_WT", "Control_KO", "Aldara_KO"); cf <- droplevels(cell[sel])
trt <- factor(sub("_.*", "", cf), levels = c("Control", "Aldara")); gen <- factor(sub(".*_", "", cf), levels = c("WT", "KO"))
d2 <- model.matrix(~ trt * gen + day[sel]); v2 <- voom(y[, sel], d2); f2 <- eBayes(lmFit(v2, d2))
t2 <- topTable(f2, coef = "trtAldara:genKO", number = Inf, sort.by = "none")
out$voom_four_day <- list(logFC = signif(t2$logFC, 12), t = signif(t2$t, 12), p = signif(t2$P.Value, 12), df0 = f2$df.prior[1])
con2 <- gzfile("validation/freeze/interaction_reference_GSE143688.json.gz", "w"); writeLines(toJSON(out, digits = NA), con2); close(con2)
cat("written:", length(out$genes), "genes\n")
