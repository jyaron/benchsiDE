# benchsiDE external validation: reference values on the public demo dataset.
#
# Input:  demo/GSE63310_counts.tsv, demo/GSE63310_design.tsv (GEO GSE63310; Law et al. 2016).
# Output: validation/gse63310_reference.json, consumed by tests/validation.spec.js.
# Run from the repository root:   Rscript validation/reference.R
#
# Every quantity mirrors the corresponding benchsiDE computation:
#   filter  edgeR::filterByExpr on all nine samples, grouped by cell type
#   norm    TMM factors on the kept genes; log2(CPM + 1) with effective library sizes
#   DE      contrast LP vs Basal fitted on the six samples of those two groups
#
# Two test covariates are defined here, because the demo design has none:
#   lcov  log2 raw library size (continuous)
#   blk   two-level block, alternating b1/b2 in matrix column order
# They exist only to exercise the covariate, blocking and batch-removal code paths.

suppressMessages({
  library(edgeR); library(limma); library(jsonlite); library(RRHO); library(RNASeqPower)
})

cnt <- as.matrix(read.delim("demo/GSE63310_counts.tsv", row.names = 1, check.names = FALSE))
des <- read.delim("demo/GSE63310_design.tsv")
grp <- factor(des$celltype[match(colnames(cnt), des$sample)])
lcov <- log2(colSums(cnt))
blk <- rep(c("b1", "b2"), length.out = ncol(cnt))

keep <- filterByExpr(cnt, group = grp)
m <- cnt[keep, ]
nf <- calcNormFactors(m, method = "TMM")
lg <- log2(t(t(m) / (colSums(m) * nf)) * 1e6 + 1)
N <- nrow(m)

# LP vs Basal, samples in the app's order (group A = Basal, then group B = LP)
selA <- which(grp == "Basal"); selB <- which(grp == "LP"); sel <- c(selA, selB)
g2 <- factor(c(rep("Basal", length(selA)), rep("LP", length(selB))), levels = c("Basal", "LP"))

fit_mod <- function(design, coef) {
  fit <- eBayes(lmFit(lg[, sel], design))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  list(fit = fit, tt = tt)
}
hits <- function(tt) c(sum(tt$adj.P.Val <= 0.05 & tt$logFC >= 1), sum(tt$adj.P.Val <= 0.05 & tt$logFC <= -1))

mod <- fit_mod(model.matrix(~g2), 2)
covc <- fit_mod(model.matrix(~g2 + lcov[sel]), 2)
covb <- fit_mod(model.matrix(~g2 + factor(blk[sel])), 2)

dge <- DGEList(m); dge$samples$norm.factors <- nf
v <- voom(dge[, sel], model.matrix(~g2))
fv <- eBayes(lmFit(v, model.matrix(~g2)))
ttv <- topTable(fv, coef = 2, number = Inf, sort.by = "none")

# moderated F over all three groups
fF <- eBayes(lmFit(lg, model.matrix(~grp)))
ttF <- topTable(fF, coef = 2:3, number = Inf, sort.by = "none")

# display-only batch removal, all samples
rbe <- removeBatchEffect(lg, batch = blk, design = model.matrix(~grp))

mds <- plotMDS(lg, plot = FALSE)

# gene sets (0-based row indices into the kept matrix, for the app)
set.seed(20260924)
sets <- list(
  SET_RANDOM = sort(sample(N, 60)),
  SET_TOP_LFC = order(-mod$tt$logFC)[1:40],
  SET_BLOCK = 1001:1080
)
fr <- fry(lg[, sel], index = sets, design = model.matrix(~g2), contrast = 2, sort = "none")
frc <- fry(lg[, sel], index = sets, design = model.matrix(~lcov[sel] + g2), contrast = 3, sort = "none")

# over-representation of LP-up hits in SET_TOP_LFC
up <- which(mod$tt$adj.P.Val <= 0.05 & mod$tt$logFC >= 1)
k <- length(intersect(up, sets$SET_TOP_LFC)); q <- length(up); ms <- length(sets$SET_TOP_LFC)
ora_p <- phyper(k - 1, ms, N - ms, q, lower.tail = FALSE)

# regression of gene 2 on gene 1 (all samples)
ols <- summary(lm(lg[2, ] ~ lg[1, ]))$coefficients

# density of sample 1: nrd0 bandwidth and exact-sum Gaussian KDE on the 512-point grid
x1 <- lg[, 1]; bw <- bw.nrd0(x1)
grid <- seq(min(x1) - 3 * bw, max(x1) + 3 * bw, length.out = 512)
probes <- c(1, 64, 128, 200, 256, 320, 400, 470, 512)
kde <- sapply(grid[probes], function(x) mean(dnorm((x - x1) / bw)) / bw)

# RRHO: LP vs Basal against ML vs Basal moderated t, step floor(N/40)
selM <- which(grp == "ML"); selC <- c(selA, selM)
gM <- factor(c(rep("Basal", length(selA)), rep("ML", length(selM))), levels = c("Basal", "ML"))
tML <- topTable(eBayes(lmFit(lg[, selC], model.matrix(~gM))), coef = 2, number = Inf, sort.by = "none")$t
step <- max(1, floor(N / 40))
rr <- RRHO(data.frame(g = paste0("g", 1:N), v = mod$tt$t), data.frame(g = paste0("g", 1:N), v = tML),
           stepsize = step, alternative = "enrichment", BY = FALSE, log10.ind = TRUE)

pw <- c(n = rnapower(depth = 20, cv = 0.4, effect = 2, alpha = 0.05, power = 0.8),
        power = rnapower(depth = 20, n = 3, cv = 0.4, effect = 2, alpha = 0.05))

idx <- c(1:25, seq(1000, N, by = 997))   # compared genes: first 25 plus a spread across the matrix
f12 <- function(x) as.numeric(format(x, digits = 15))
ref <- list(
  source = "Rscript validation/reference.R",
  packages = list(R = paste(R.version$major, R.version$minor, sep = "."),
                  edgeR = as.character(packageVersion("edgeR")), limma = as.character(packageVersion("limma")),
                  RRHO = as.character(packageVersion("RRHO")), RNASeqPower = as.character(packageVersion("RNASeqPower"))),
  covariates = list(lcov = f12(lcov), blk = blk),
  n_kept = N, kept_checksum = sum(which(keep)),
  tmm = f12(nf),
  compared_genes = idx - 1,
  logcpm_1 = f12(lg[1, ]),
  mod = list(t = f12(mod$tt$t[idx]), p = f12(mod$tt$P.Value[idx]), q = f12(mod$tt$adj.P.Val[idx]),
             df0 = f12(mod$fit$df.prior), s20 = f12(mod$fit$s2.prior), hits = hits(mod$tt)),
  voom = list(t = f12(ttv$t[idx]), df0 = f12(fv$df.prior), hits = hits(ttv)),
  cov_continuous = list(t = f12(covc$tt$t[idx]), df0 = f12(covc$fit$df.prior), hits = hits(covc$tt)),
  cov_block = list(t = f12(covb$tt$t[idx]), df0 = f12(covb$fit$df.prior), hits = hits(covb$tt)),
  modF = list(F = f12(ttF$F[idx]), df0 = f12(fF$df.prior)),
  rbe_rows = f12(rbe[idx[1:10], ]),
  mds = list(x = f12(mds$x), y = f12(mds$y), ve = f12(mds$var.explained[1:2])),
  sets = lapply(sets, function(s) s - 1),
  fry = list(p = f12(fr$PValue), pm = f12(fr$PValue.Mixed), dir = as.character(fr$Direction)),
  fry_cov = list(p = f12(frc$PValue), pm = f12(frc$PValue.Mixed), dir = as.character(frc$Direction)),
  ora = list(k = k, q = q, m = ms, N = N, p = f12(ora_p)),
  ols = f12(c(ols[1, 1], ols[2, 1], ols[2, 2], ols[2, 4])),
  kde = list(bw = f12(bw), probes = probes - 1, y = f12(kde)),
  rrho = list(step = step, dim = dim(rr$hypermat), diag = f12(diag(rr$hypermat)),
              row1 = f12(rr$hypermat[1, ]), max = f12(max(rr$hypermat, na.rm = TRUE))),
  rnapower = f12(pw)
)
writeLines(toJSON(ref, digits = I(15), auto_unbox = TRUE, pretty = FALSE), "validation/gse63310_reference.json")
cat("wrote validation/gse63310_reference.json; kept genes:", N, "\n")
