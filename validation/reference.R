# benchsiDE external validation: reference values on the public demo dataset.
#
# Input:  demo/GSE63310_counts.tsv, demo/GSE63310_design.tsv (GEO GSE63310; Law et al. 2016).
# Output: validation/gse63310_reference.json, consumed by tests/validation.spec.js.
# Run from the repository root:   Rscript validation/reference.R
#
# Every quantity mirrors the corresponding benchsiDE computation:
#   filter  edgeR::filterByExpr on all nine samples, grouped by cell type
#   norm    TMM factors on the kept genes; edgeR log-CPM, cpm(log = TRUE, prior.count = 2)
#   DE      moderated t with eBayes(trend = TRUE) (limma-trend); voom with its own weights
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
lg <- cpm(DGEList(m, norm.factors = nf), log = TRUE, prior.count = 2)
N <- nrow(m)

# LP vs Basal, samples in the app's order (group A = Basal, then group B = LP)
selA <- which(grp == "Basal"); selB <- which(grp == "LP"); sel <- c(selA, selB)
g2 <- factor(c(rep("Basal", length(selA)), rep("LP", length(selB))), levels = c("Basal", "LP"))

fit_mod <- function(design, coef) {
  fit <- eBayes(lmFit(lg[, sel], design), trend = TRUE)
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

# full design (application default, "variance from all groups"): ~ group over all samples, LP vs Basal coefficient
gF <- relevel(factor(grp), ref = "Basal"); XF <- model.matrix(~gF)
fitFull <- function(X) { f <- eBayes(lmFit(lg, X), trend = TRUE); list(fit = f, tt = topTable(f, coef = "gFLP", number = Inf, sort.by = "none")) }
fmod <- fitFull(XF)
fcovc <- fitFull(model.matrix(~gF + lcov))
fcovb <- fitFull(model.matrix(~gF + factor(blk)))
vF <- voom(dge, XF); fvF <- eBayes(lmFit(vF, XF)); ttvF <- topTable(fvF, coef = "gFLP", number = Inf, sort.by = "none")

# moderated F over all three groups
fF <- eBayes(lmFit(lg, model.matrix(~grp)), trend = TRUE)
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
# FRY's robust standardization depends on the basis of the residual space, hence on sample order:
# use the application's order (group A, group B, then the other groups in order of appearance)
ordF <- c(selA, selB, setdiff(seq_along(grp), c(selA, selB)))
frF <- fry(lg[, ordF], index = sets, design = XF[ordF, ], contrast = which(colnames(XF) == "gFLP"), sort = "none")

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
tML <- topTable(eBayes(lmFit(lg[, selC], model.matrix(~gM)), trend = TRUE), coef = 2, number = Inf, sort.by = "none")$t
step <- max(1, floor(N / 40))
rr <- RRHO(data.frame(g = paste0("g", 1:N), v = mod$tt$t), data.frame(g = paste0("g", 1:N), v = tML),
           stepsize = step, alternative = "enrichment", BY = FALSE, log10.ind = TRUE)

pw <- c(n = rnapower(depth = 20, cv = 0.4, effect = 2, alpha = 0.05, power = 0.8),
        power = rnapower(depth = 20, n = 3, cv = 0.4, effect = 2, alpha = 0.05))

idx <- c(1:25, seq(1000, N, by = 997))   # compared genes: first 25 plus a spread across the matrix
f12 <- function(x) as.numeric(format(x, digits = 15))
# co-expression table: every compared gene against gene 1 (Pearson r, OLS slope on gene 1, cor.test p)
cidx <- setdiff(idx, 1)
coex <- t(sapply(cidx, function(g) c(cor(lg[g, ], lg[1, ]), coef(lm(lg[g, ] ~ lg[1, ]))[2],
                                     cor.test(lg[g, ], lg[1, ])$p.value)))
# hub neighbourhood, within-group mode: partial correlation of each compared gene with gene 1 given
# cell type = correlation of residuals from ~ grp; p = t test of gene 1's coefficient in lm(gene ~ gene1 + grp)
pcor <- t(sapply(cidx, function(g) c(cor(resid(lm(lg[g, ] ~ grp)), resid(lm(lg[1, ] ~ grp))),
                                      summary(lm(lg[g, ] ~ lg[1, ] + grp))$coefficients[2, 4])))

# ---- hub across datasets ----
# soft connectivity k_i = sum_j |r_ij|^6 over the first 200 kept genes (unsigned), all samples
cidxK <- 1:200; Ck <- abs(cor(t(lg[cidxK, ]))); diag(Ck) <- 0; kconn <- rowSums(Ck^6)
# module-preservation statistics: WGCNA 1.74's own .coreCalcForExpr, sourced at run time from the CRAN source
# (not redistributed here); reference module = genes 1..12 over Basal and LP samples, test = the same genes over
# LP and ML samples, and five fixed random 12-gene sets as permutation draws
pres <- NULL
wg <- tryCatch({
  td <- tempfile(); dir.create(td)
  tf <- file.path(td, "WGCNA.tar.gz")
  ok <- FALSE
  for (u in c("https://cloud.r-project.org/src/contrib/WGCNA_1.74.tar.gz", "https://cloud.r-project.org/src/contrib/Archive/WGCNA/WGCNA_1.74.tar.gz"))
    if (!ok) ok <- !inherits(try(download.file(u, tf, quiet = TRUE, mode = "wb"), silent = TRUE), "try-error")
  if (!ok) stop("WGCNA source not available")
  untar(tf, exdir = td); td
}, error = function(e) NULL)
if (!is.null(wg)) {
  grab <- function(path, name) { t <- paste(readLines(path, warn = FALSE), collapse = "\n"); i <- regexpr(paste0("\n", gsub(".", "\\.", name, fixed = TRUE), " = function"), t)
    s0 <- i + 1; ch <- strsplit(substring(t, s0), "")[[1]]; d <- 0; st <- FALSE
    for (j in seq_along(ch)) { if (ch[j] == "{") { d <- d + 1; st <- TRUE } else if (ch[j] == "}") { d <- d - 1; if (st && d == 0) break } }
    substring(t, s0, s0 + j - 1) }
  e <- new.env()
  with(e, { indentSpaces <- function(indent = 0) ""; printFlush <- function(...) invisible(NULL); moduleColor.getMEprefix <- function() "ME"
    .checkAndScaleWeights <- function(weights, expr, ...) weights; colVars <- function(x, na.rm = FALSE, ...) apply(x, 2, var, na.rm = na.rm)
    ..minNSamples <- 4; ..minNGenes <- 4
    cor <- function(x, y = NULL, use = "everything", method = "pearson", quick = 0, weights.x = NULL, weights.y = NULL, ...) stats::cor(x, y, use = use, method = method) })
  Rd <- file.path(wg, "WGCNA", "R")
  for (nm in list(c("Functions.R", "moduleEigengenes"), c("Functions.R", "signedKME"), c("Functions.R", "prepComma"),
                  c("modulePreservation.R", ".coreCalcForExpr"), c("modulePreservation.R", ".clusterCoeff"), c("modulePreservation.R", ".MAR")))
    eval(parse(text = grab(file.path(Rd, nm[1]), nm[2])), envir = e)
  opt <- list(densityOnly = FALSE, calculatePermutation = FALSE, corFnc = "cor", corOptions = "use = 'p'", quickCor = 0, nType = 1,
              MEgold = "MEgold", MEgrey = "MEgrey", calculateClusterCoeff = FALSE)
  refS <- which(grp %in% c("Basal", "LP")); tstS <- which(grp %in% c("LP", "ML")); modG <- 1:12
  presSets <- list(seq(13, 13 + 11 * 37, by = 37), seq(41, 41 + 11 * 53, by = 53), seq(97, 97 + 11 * 29, by = 29), seq(201, 201 + 11 * 61, by = 61), seq(333, 333 + 11 * 17, by = 17))
  st7 <- function(tg) { s <- suppressWarnings(e$.coreCalcForExpr(t(lg[modG, refS]), t(lg[modG, refS]), t(lg[tg, tstS]), rep("m1", 12), NULL, NULL, NULL, opt))
    c(s$proVar[1, 2], s$meanSignAwareKME[1, 2], s$MeanSignAwareCorDat[1, 2], s$MeanAdj[1, 2], s$corkIM, s$corkME, s$ICORdat) }
  pres <- list(refS = refS - 1, tstS = tstS - 1, mod = modG - 1, sets = lapply(presSets, function(x) x - 1),
               obs = f12(st7(modG)), perm = lapply(presSets, function(x) f12(st7(x))), wgcna = "1.74")
}
# robust rank aggregation: four fixed permutations of 300 items, RobustRankAggreg::aggregateRanks (RRA, exact = FALSE)
set.seed(20261001); items <- paste0("g", 1:300)
rl <- list(items, sample(items), c(items[1:20], sample(items[21:300])), c(sample(items[1:40]), sample(items[41:300])))
rra <- if (requireNamespace("RobustRankAggreg", quietly = TRUE)) { ag <- RobustRankAggreg::aggregateRanks(glist = rl, N = 300, method = "RRA", exact = FALSE)
  list(lists = lapply(rl, function(x) as.integer(sub("g", "", x)) - 1), item = as.integer(sub("g", "", rownames(ag))) - 1, score = f12(ag$Score)) } else NULL


# ---- sample traits ----
# synthetic sample-level trait and block on the demo samples (in matrix column order); limma lmFit + eBayes(trend = TRUE)
trait_x <- c(3.1, 7.4, 1.2, 5.5, 9.0, 2.8, 6.3, 4.1, 8.7); block_x <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
tz <- (trait_x - mean(trait_x)) / sd(trait_x)
fT <- eBayes(lmFit(lg, cbind(1, tz)), trend = TRUE)
XG <- cbind(1, tz, as.numeric(grp == "LP"), as.numeric(grp == "ML"))
dcT <- duplicateCorrelation(lg, XG, block = block_x)
fTB <- eBayes(lmFit(lg, XG, block = block_x, correlation = dcT$consensus.correlation), trend = TRUE)
camT <- camera(lg, cbind(1, tz), contrast = 2, index = list(a = 1:30, b = seq(500, 2000, 50), c = 3000:3100), inter.gene.cor = NA, sort = FALSE)
trait <- list(x = trait_x, block = block_x, t_all = f12(fT$t[1:50, 2]), p_all = f12(fT$p.value[1:50, 2]), df0_all = f12(fT$df.prior[1]),
              rho = f12(dcT$consensus.correlation), t_blk = f12(fTB$t[1:50, 2]), p_blk = f12(fTB$p.value[1:50, 2]),
              cam_p = f12(camT$PValue), cam_cor = f12(camT$Correlation))
if (!is.null(wg)) {
  meT <- e$moduleEigengenes(t(lg[1:40, ]), colors = rep("m", 40), excludeGrey = FALSE)
  trait$me <- f12(meT$eigengenes[, 1]); trait$me_t <- f12(summary(lm(meT$eigengenes[, 1] ~ tz))$coefficients[2, 3])
}

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
             df0 = f12(mod$fit$df.prior), s20 = f12(median(mod$fit$s2.prior)), hits = hits(mod$tt)),
  full_mod = list(t = f12(fmod$tt$t[idx]), p = f12(fmod$tt$P.Value[idx]), df0 = f12(fmod$fit$df.prior), dfr = fmod$fit$df.residual[1], hits = hits(fmod$tt)),
  full_covc = list(t = f12(fcovc$tt$t[idx]), df0 = f12(fcovc$fit$df.prior), hits = hits(fcovc$tt)),
  full_covb = list(t = f12(fcovb$tt$t[idx]), df0 = f12(fcovb$fit$df.prior), hits = hits(fcovb$tt)),
  full_voom = list(t = f12(ttvF$t[idx]), df0 = f12(fvF$df.prior), hits = hits(ttvF)),
  full_fry = list(p = f12(frF$PValue), pm = f12(frF$PValue.Mixed), dir = as.character(frF$Direction)),
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
  coexp = list(genes = cidx - 1, r = f12(coex[, 1]), slope = f12(coex[, 2]), p = f12(coex[, 3])),
  pcor = list(genes = cidx - 1, r = f12(pcor[, 1]), p = f12(pcor[, 2])),
  kconn = f12(kconn),
  trait = trait,
  pres = pres,
  rra = rra,
  kde = list(bw = f12(bw), probes = probes - 1, y = f12(kde)),
  rrho = list(step = step, dim = dim(rr$hypermat), diag = f12(diag(rr$hypermat)),
              row1 = f12(rr$hypermat[1, ]), max = f12(max(rr$hypermat, na.rm = TRUE))),
  rnapower = f12(pw)
)
writeLines(toJSON(ref, digits = I(15), auto_unbox = TRUE, pretty = FALSE), "validation/gse63310_reference.json")
cat("wrote validation/gse63310_reference.json; kept genes:", N, "\n")
