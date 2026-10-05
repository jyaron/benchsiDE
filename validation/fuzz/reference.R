# Randomized comparison of benchsiDE differential expression: independent reference in limma/edgeR.
# Implements the documented rules of the application from its documentation, not from its code:
#  filtering (filterByExpr with the grouping factor, or CPM >= t in >= m samples) and TMM on all loaded samples;
#  log-CPM with prior count 2; samples without a group label and user-excluded samples left out of the fit;
#  "all": every group with an included sample is fitted, "pair": the two compared groups only;
#  samples with a missing value (blank or NA) in a selected covariate are left out (moderated t and voom);
#  a numeric covariate with more than two distinct values among the fitted samples is a slope unless its name
#  marks an identifier; otherwise it is categorical; non-estimable terms are dropped as lmFit does;
#  design ~0 + group + covariates; contrasts.fit(B - A); eBayes(trend = TRUE) / treat(trend = TRUE) on log-CPM,
#  or voom + eBayes / treat; Welch: per-gene unequal-variance t on log-CPM of the two groups (no covariates).
# Usage: Rscript validation/fuzz/reference.R <load name>    -> work/<load>.ref.json
suppressMessages({ library(limma); library(edgeR); library(jsonlite) })
args <- commandArgs(TRUE); W <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "work")
if (!length(W) || !dir.exists(W)) W <- "validation/fuzz/work"
spec <- fromJSON(file.path(W, "configs.json"), simplifyVector = FALSE)
name <- args[1]; L <- Filter(function(l) l$name == name, spec$loads)[[1]]; CF <- Filter(function(c) c$load == name, spec$configs)
SETMODE <- length(args) > 1 && args[2] == "sets"
if (SETMODE) {
  sl <- readLines(file.path(W, paste0(name, ".sets.jsonl"))); S0 <- fromJSON(substring(sl[startsWith(sl, "S ")][1], 3), simplifyVector = FALSE)
  SETS <- lapply(S0$sets, function(z) unlist(z$rows) + 1L); names(SETS) <- sapply(S0$sets, function(z) z$name)
  UNIV <- unlist(S0$universe) + 1L
  CF <- Filter(function(c) c$m %in% c("mod", "voom"), CF)
}
cnt <- as.matrix(read.delim(file.path(W, paste0(name, ".matrix.tsv")), row.names = 1, check.names = FALSE))
des <- read.delim(file.path(W, paste0(name, ".design.tsv")), colClasses = "character", na.strings = character(0), check.names = FALSE)
stopifnot(identical(des$sample, colnames(cnt)))
idre <- "patient|subject|donor|pair|individual|animal|block|batch|litter|replicate|(^|[_ .])id$|^id([_ .]|$)|_id$|sample.?id"
is_miss <- function(v) is.na(v) | trimws(v) == "" | toupper(trimws(v)) == "NA"
grp_all <- des[[L$factor]]; grp_all[is_miss(grp_all)] <- "NA"
units <- if (is.null(L$units)) "counts" else L$units
lib0 <- colSums(cnt)
if (units == "counts") {
  keep <- if (L$fmode == "fbe") filterByExpr(DGEList(cnt), group = factor(grp_all)) else rowSums(t(t(cnt) / lib0 * 1e6) >= 1) >= 2
  x <- DGEList(cnt[keep, , drop = FALSE])
  if (L$norm == "tmm") x <- calcNormFactors(x)
  lc <- cpm(x, log = TRUE, prior.count = 2)
} else {
  # normalized linear values: used as is, log2(x + 1); log values: used as is, linear 2^x - 1 when no value is
  # negative (log2(x + 1) input), else 2^x. Filter: linear value >= 1 in >= 2 samples. TMM and voom do not apply.
  lin <- if (units == "norm") cnt else if (min(cnt) >= 0) 2^cnt - 1 else 2^cnt
  keep <- rowSums(lin >= 1) >= 2
  lc <- if (units == "norm") log2(cnt[keep, , drop = FALSE] + 1) else cnt[keep, , drop = FALSE]
  x <- NULL
}
ids <- rownames(lc)
step <- max(1L, floor(length(ids) / 400)); samp <- seq(1L, length(ids), by = step)
res <- list(); res[["_load"]] <- list(nG = length(ids), ids = ids)
welch_vec <- function(a, b) {
  ma <- rowMeans(a); mb <- rowMeans(b); na <- ncol(a); nb <- ncol(b)
  va <- apply(a, 1, var); vb <- apply(b, 1, var); se2 <- va / na + vb / nb; se <- sqrt(se2)
  t <- (mb - ma) / se; df <- se2^2 / ((va / na)^2 / (na - 1) + (vb / nb)^2 / (nb - 1))
  p <- 2 * pt(-abs(t), df)
  const <- se < 10 * .Machine$double.eps * pmax(abs(ma), abs(mb)) | se2 == 0
  t[const] <- 0; p[const] <- 1
  list(lfc = mb - ma, t = t, p = p)
}
for (c in CF) {
  out <- list(id = c$id)
  r <- tryCatch({
    inc <- !(des$sample %in% unlist(c$excl)) & grp_all != "NA"
    A <- c$A; B <- c$B; covs <- unlist(c$cov); if (is.null(covs)) covs <- character(0)
    if (c$m == "welch") covs <- character(0)
    ok <- inc
    for (cv in covs) ok <- ok & !is_miss(des[[cv]])
    nA <- sum(ok & grp_all == A); nB <- sum(ok & grp_all == B)
    if (nA < 2 || nB < 2) stop(sprintf("fewer than 2 samples: %s %d, %s %d", A, nA, B, nB))
    groups <- if (c$fit == "all" && c$m != "welch") c(A, B, setdiff(sort(unique(grp_all[ok]), method = "radix"), c(A, B))) else c(A, B)
    use <- ok & grp_all %in% groups
    if (c$m == "welch") {
      w <- welch_vec(lc[, use & grp_all == A, drop = FALSE], lc[, use & grp_all == B, drop = FALSE])
      lfc <- w$lfc; t <- w$t; p <- w$p; dfr <- NA; df0 <- NA; kept <- character(0); dropped <- character(0); types <- list()
    } else {
      dd <- data.frame(g = factor(grp_all[use], levels = groups)); types <- list()
      for (cv in covs) {
        v <- des[[cv]][use]; lv <- unique(v); num <- suppressWarnings(as.numeric(v))
        cont <- all(!is.na(num)) && length(lv) > 2 && !grepl(idre, cv, ignore.case = TRUE, perl = TRUE)
        types[[cv]] <- if (cont) "continuous" else "categorical"
        dd[[cv]] <- if (cont) num else factor(v, levels = sort(lv, method = "radix"))
      }
      fml <- as.formula(paste("~ 0 + g", if (length(covs)) paste("+", paste(sprintf("`%s`", covs), collapse = " + ")) else ""))
      X <- model.matrix(fml, dd); colnames(X) <- gsub("`", "", colnames(X))
      q <- qr(X); keepc <- sort(q$pivot[seq_len(q$rank)])
      if (!all(c(1, 2) %in% keepc)) stop("contrast not estimable")
      dropped_cols <- colnames(X)[-keepc]; X <- X[, keepc, drop = FALSE]
      kept <- character(0); dropped <- character(0)
      for (cv in covs) {
        cols <- colnames(X)[startsWith(colnames(X), cv)]
        if (!length(cols)) dropped <- c(dropped, cv) else kept <- c(kept, cv)
      }
      if (nrow(X) - ncol(X) < 1) stop("no residual degrees of freedom")
      cm <- matrix(0, ncol(X), 1); cm[1, 1] <- -1; cm[2, 1] <- 1
      if (c$m == "voom") { if (is.null(x)) stop("voom needs raw counts"); fit <- lmFit(voom(x[, use], X), X); trend <- FALSE } else { fit <- lmFit(lc[, use], X); trend <- TRUE }
      fit <- contrasts.fit(fit, cm)
      f2 <- if (isTRUE(c$treat)) treat(fit, lfc = c$thr, trend = trend) else eBayes(fit, trend = trend)
      lfc <- f2$coefficients[, 1]; t <- f2$t[, 1]; p <- f2$p.value[, 1]; dfr <- f2$df.residual[1]; df0 <- f2$df.prior
    }
    qv <- p.adjust(p, "BH")
    call <- if (isTRUE(c$treat)) qv <= c$fdr else qv <= c$fdr & abs(lfc) >= c$thr
    if (SETMODE) {
      # FRY on the unweighted log-CPM with the same design and contrast (limma::fry defaults)
      # limma::fry with design ~ group + covariates, reference level A, contrast = the B coefficient (the
      # parametrization the application uses; FRY's robust variance depends on the parametrization, also in limma)
      gord <- unlist(S0$groupOrder); oth <- setdiff(gord[gord %in% groups], c(A, B))
      gcols <- colnames(X)[seq_along(groups)]; names(gcols) <- groups
      covcols <- colnames(X)[-seq_along(groups)]
      XF <- cbind(1, X[, gcols[oth], drop = FALSE], X[, covcols, drop = FALSE], X[, gcols[B], drop = FALSE])
      fr <- fry(lc[, use], index = SETS, design = XF, contrast = ncol(XF), sort = "none")
      # the same test with the samples in a different order: limma's FRY standardization uses the largest squared
      # residual effect, whose value depends on the (arbitrary) basis of the residual space, so this measures
      # limma's own variation for this design
      set.seed(1); po <- sample(nrow(XF))
      fr2 <- fry(lc[, use][, po, drop = FALSE], index = SETS, design = XF[po, , drop = FALSE], contrast = ncol(XF), sort = "none")
      up <- intersect(which(call & lfc > 0), UNIV); N <- length(UNIV); q <- length(up)
      ora <- if (q >= 3) t(sapply(SETS, function(st) { st <- intersect(st, UNIV); k <- length(intersect(st, up)); c(length(st), k, phyper(k - 1, length(st), N - length(st), q, lower.tail = FALSE)) })) else NULL
      SETRES <- list(fry = lapply(seq_along(SETS), function(i) list(names(SETS)[i], fr$NGenes[i], fr$PValue[i], fr$PValue.Mixed[i], fr$Direction[i], fr2$PValue[i], fr2$PValue.Mixed[i])),
                  ora = if (is.null(ora)) NULL else lapply(seq_len(nrow(ora)), function(i) list(rownames(ora)[i], ora[i, 1], ora[i, 2], ora[i, 3])),
                  kept = kept, dropped = dropped, n = sum(use))
    }
    gi <- which(call) - 1L
    if (SETMODE) SETRES else list(calls = ifelse(lfc[call] > 0, gi, -gi - 1L), lfc = unname(lfc[samp]), t = unname(t[samp]), p = unname(p[samp]),
         n = sum(use), nA = nA, nB = nB, dfr = dfr, df0 = df0, groups = groups, kept = kept, dropped = dropped, types = types,
         dropNA = sum(inc & grp_all %in% groups & !ok))
  }, error = function(e) list(err = conditionMessage(e)))
  res[[c$id]] <- c(out, r)
}
writeLines(toJSON(res, auto_unbox = TRUE, digits = NA, null = "null", na = "null"), file.path(W, paste0(name, if (SETMODE) ".setsref.json" else ".ref.json")))
cat(name, length(CF), "configurations\n")
