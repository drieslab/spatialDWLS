# Timing and accuracy against known truth.
#   Rscript dev/bench.R <n_spots> [n_pure_niches]      (env NC = n_cell)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("dev/sim.R")
args <- commandArgs(TRUE)
n <- as.integer(args[[1]])
n_pure <- if (length(args) > 1) as.integer(args[[2]]) else 0L
nc <- as.numeric(Sys.getenv("NC", "50"))

s <- simulate_dwls(n_spots = n, n_pure = n_pure)
acc <- function(p) {
    ok <- !is.na(colSums(p))
    d <- p[rownames(s$truth), ok] - s$truth[, ok]
    sprintf("RMSE %.4f  NA spots %d", sqrt(mean(d^2)), sum(!ok))
}
cat(sprintf("spots %d  pure niches %d  n_cell %g\n", n, n_pure, nc))
t1 <- system.time(p1 <- dwls(s$expr, s$sig, s$clusters, n_cell = nc))
cat(sprintf("serial     %6.1fs  %s\n", t1[["elapsed"]], acc(p1)))
if (requireNamespace("BiocParallel", quietly = TRUE) &&
    .Platform$OS.type == "unix") {
    bp <- BiocParallel::MulticoreParam(8)
    t2 <- system.time(p2 <- dwls(s$expr, s$sig, s$clusters, n_cell = nc,
        BPPARAM = bp))
    cat(sprintf("8 workers  %6.1fs  %s  identical: %s\n", t2[["elapsed"]],
        acc(p2), isTRUE(all.equal(p1, p2))))
}
