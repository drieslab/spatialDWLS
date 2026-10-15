# Old (pre-extraction, Giotto gsource) DWLS vs spatialDWLS::dwls().
#   Rscript dev/equivalence.R ../GiottoSuite/Giotto [n_spots]
#
# The old solver is read from Giotto's git history at the last commit before
# the rewrite, with its two known bugs patched (negatives clipped only on the
# last row; single-type clusters left at zero), so the comparison isolates the
# rewrite. On the simulation below and the GiottoData mini Visium object the
# two agree to ~1e-14.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("dev/sim.R")
args <- commandArgs(TRUE)
giotto <- normalizePath(args[[1]])
n <- if (length(args) > 1) as.integer(args[[2]]) else 2000L
ref <- Sys.getenv("DWLS_OLD_REF", "cc39eb8f")  # gsource before the rewrite

old_src <- system2("git", c("-C", giotto, "show",
    paste0(ref, ":R/spatial_enrichment.R")), stdout = TRUE)
from <- grep("^#' @title enrich_deconvolution", old_src)
to <- grep("^#' @title Estimate cell type proportions with spatialDWLS",
    old_src) - 1L
core <- old_src[from:to]
core <- sub("for (i in dim(dwls_results)[1])",
    "for (i in seq_len(dim(dwls_results)[1]))", core, fixed = TRUE)
core <- sub("which(cluster_info == cluster_sort[i])] == 1",
    "which(cluster_info == cluster_sort[i])] <- 1", core, fixed = TRUE)

old <- new.env()
old$mixedsort <- function(x) GiottoUtils::mixedsort(x)
old$package_check <- function(...) invisible(TRUE)
eval(parse(text = core), envir = old)
if (!requireNamespace("Rfast", quietly = TRUE)) stop("the old code needs Rfast")
if (!requireNamespace("GiottoUtils", quietly = TRUE)) {
    stop("the old code needs GiottoUtils (local_seed)")
}

old_dwls <- function(expr, sig, clusters, n_cell = 50, cutoff = 2) {
    genes <- intersect(rownames(sig), rownames(expr))
    log_e <- as.matrix(expr[genes, ])
    nolog <- 2^log_e - 1
    cl <- unname(clusters[colnames(expr)])
    r1 <- old$enrich_deconvolution(nolog, log_e, cl, sig[genes, ], cutoff)
    old$spot_deconvolution(nolog, cl, sig[genes, ],
        ifelse(r1 >= 1 / n_cell, 1, 0))
}

s <- simulate_dwls(n_spots = n, n_pure = 3L)
for (nc in c(50, 5)) {
    t_old <- system.time(p_old <- old_dwls(s$expr, s$sig, s$clusters, nc))
    t_new <- system.time(p_new <- dwls(s$expr, s$sig, s$clusters, n_cell = nc))
    p_old <- p_old[rownames(p_new), colnames(p_new)]
    cat(sprintf(
        "n_cell %2g: old %.1fs  new %.1fs  max|diff| %.2e  NA old/new %d/%d\n",
        nc, t_old[["elapsed"]], t_new[["elapsed"]],
        max(abs(p_old - p_new), na.rm = TRUE),
        sum(is.na(p_old)), sum(is.na(p_new))))
}
