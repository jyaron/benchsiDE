# metafor reference for validation/fuzz/stat_tests.py
suppressMessages({ library(metafor); library(jsonlite) })
a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); W <- file.path(dirname(normalizePath(f)), "work")
J <- fromJSON(file.path(W, "meta_fuzz.json"), simplifyVector = FALSE)
meth <- list(c("FE", "z"), c("DL", "z"), c("REML", "z"), c("DL", "knha"), c("REML", "knha"), c("REML", "adhoc"))
RES <- lapply(J$cases, function(cc) { y <- unlist(cc$y); v <- unlist(cc$v)
  lapply(meth, function(m) { r <- tryCatch(suppressWarnings(rma(y, v, method = m[1], test = m[2], control = list(threshold = 1e-12, maxiter = 1000, stepadj = 0.5))), error = function(e) NULL)
    if (is.null(r)) return(list(err = TRUE))
    pr <- if (m[1] != "FE" && length(y) >= 3) tryCatch(predict(r, predtype = "Riley"), error = function(e) NULL) else NULL
    list(mu = as.numeric(r$b), se = r$se, p = r$pval, tau2 = r$tau2, I2 = r$I2 / 100, lo = r$ci.lb, hi = r$ci.ub,
         plo = if (is.null(pr)) NA else pr$pi.lb, phi = if (is.null(pr)) NA else pr$pi.ub, Q = r$QE) }) })
writeLines(toJSON(RES, auto_unbox = TRUE, digits = NA, na = "null"), file.path(W, "meta_fuzz_ref.json"))
cat("metafor", as.character(packageVersion("metafor")), length(RES), "cases\n")
