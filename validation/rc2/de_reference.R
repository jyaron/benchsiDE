# benchsiDE 1.0.0-rc.2: reference results for differential expression with covariates and TREAT,
# computed with the standard limma workflow (limma 3.66.0, edgeR 4.8.2):
#   filterByExpr, calcNormFactors (TMM);
#   design ~0 + group + covariates (aliased columns dropped);
#   moderated t: lmFit on log-CPM (prior.count 2), contrasts.fit(B - A), eBayes(trend = TRUE) or treat(trend = TRUE);
#   voom: voom(counts, design), lmFit, contrasts.fit(B - A), eBayes() or treat().
# "all": every group is fitted; "pair": only the two compared groups. Samples with a missing covariate are dropped.
suppressMessages({ library(limma); library(edgeR); library(jsonlite) })

de_reference <- function(counts, design_df, factor, cfg, out_prefix) {
  grp_all <- factor(design_df[[factor]])
  x <- DGEList(counts, group = grp_all)
  keep <- filterByExpr(x, group = grp_all)
  x <- x[keep, , keep.lib.sizes = FALSE]
  x <- calcNormFactors(x)
  lcpm <- cpm(x, log = TRUE, prior.count = 2)
  res <- list()
  for (i in seq_along(cfg)) {
    c <- cfg[[i]]; A <- c$A; B <- c$B; covs <- unlist(c$cov)
    use <- if (c$fit == "all") rep(TRUE, ncol(x)) else grp_all %in% c(A, B)
    for (cv in covs) { v <- design_df[[cv]]; use <- use & !is.na(v) & trimws(as.character(v)) != "" & toupper(trimws(as.character(v))) != "NA" }
    g <- droplevels(grp_all[use])
    dd <- data.frame(g = g)
    for (cv in covs) {
      v <- design_df[[cv]][use]
      num <- suppressWarnings(as.numeric(as.character(v)))
      dd[[cv]] <- if (!any(is.na(num)) && length(unique(num)) > 2) num else factor(v)
    }
    X <- model.matrix(as.formula(paste("~ 0 + g", if (length(covs)) paste("+", paste(covs, collapse = " + ")) else "")), dd)
    X <- X[, sort(qr(X)$pivot[seq_len(qr(X)$rank)]), drop = FALSE]
    cm <- matrix(0, ncol(X), 1, dimnames = list(colnames(X), "BvsA")); cm[paste0("g", B), 1] <- 1; cm[paste0("g", A), 1] <- -1
    if (c$m == "voom") { fit <- lmFit(voom(x[, use], X), X); trend <- FALSE }
    else { fit <- lmFit(lcpm[, use], X); trend <- TRUE }
    fit <- contrasts.fit(fit, cm)
    f2 <- if (is.null(c$tr)) eBayes(fit, trend = trend) else treat(fit, lfc = c$tr, trend = trend)
    write.csv(data.frame(id = rownames(x), lfc = f2$coefficients[, 1], t = f2$t[, 1], p = f2$p.value[, 1]),
              sprintf("%s_%d.csv", out_prefix, i - 1), row.names = FALSE)
    res[[i]] <- list(df_prior = f2$df.prior, df_residual = f2$df.residual[1], n = sum(use))
  }
  res
}
