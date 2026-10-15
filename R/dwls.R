# spatialDWLS runs in two rounds, both built from per-cluster summaries plus
# per-spot solves:
#
#   1. cluster round. PAGE z-scores pick the candidate cell types of each
#      cluster; every spot is then solved against that cluster's types.
#   2. spot round. Types below 1 / n_cell in round 1 are dropped per spot and
#      the spot is solved again against what is left.
#
# The per-cluster summaries (cluster means, feature reference levels, best z
# per type, types present) are reductions over spots; a per-spot solve reads
# only its own column. `dwls_stream()` is organised around that split: it reads
# spots window by window through a callback, folds the summaries across
# windows, and hands each window to the helpers in plan.R and solve.R, which
# accept any subset of spots. `dwls()` is the one-window case.

#' spatialDWLS: cell type deconvolution of spatial spots
#'
#' @description
#' Estimate what fraction of each spot is made up of each cell type, given a
#' mean-expression reference built from single-cell data.
#'
#' * [dwls()] deconvolves an in-memory expression matrix.
#' * [dwls_stream()] does the same for a backend that reads spots in windows,
#'   so the expression never has to fit in memory at once.
#' * [dwls_signature()] builds the reference from a single-cell matrix.
#'
#' @references
#' Dong R, Yuan GC (2021). SpatialDWLS: accurate deconvolution of spatial
#' transcriptomic data. *Genome Biology* 22, 145.
#' \doi{10.1186/s13059-021-02362-7}
#'
#' Tsoucas D, Dong R, Chen H, Zhu Q, Guo G, Yuan GC (2019). Accurate
#' estimation of cell-type composition from gene expression data. *Nature
#' Communications* 10, 2975. \doi{10.1038/s41467-019-10802-z}
#' @keywords internal
"_PACKAGE"


#' Deconvolve spots into cell type proportions
#'
#' @description
#' Runs both rounds of spatialDWLS on an expression matrix. Each spot's
#' result is a composition: the proportions of a solved spot sum to one.
#'
#' @param expr log-scaled, normalised expression: features (rows) by spots
#'   (columns), as a base `matrix` or a `Matrix` (sparse is kept sparse).
#'   Spot IDs are taken from the column names.
#' @param sign_matrix mean expression per cell type: signature features (rows)
#'   by cell type (columns), on the reverse-logged scale, e.g. from
#'   [dwls_signature()]. Not a binary marker matrix.
#' @param clusters cluster of every spot. Either named by spot ID (matched
#'   against `colnames(expr)`; extra entries, such as the clusters of a whole
#'   dataset when `expr` is a subset, are ignored) or unnamed in the column
#'   order of `expr`. Spots are grouped by cluster to pick candidate cell types
#'   in round 1. Every spot needs a non-`NA` cluster.
#' @param n_cell numeric. expected number of cells per spot. In round 2 a cell
#'   type is kept for a spot only when its round-1 proportion is at least
#'   `1 / n_cell`.
#' @param cutoff numeric. PAGE enrichment z-score a cell type must reach in at
#'   least one spot of a cluster to be a candidate for that cluster. The two
#'   highest-scoring types are always candidates.
#' @param logbase numeric. log base of `expr`.
#' @param na_as_zero logical. spots that cannot be solved (no signal on the
#'   signature features) are `NA` by default; `TRUE` returns them as 0.
#' @param BPPARAM a `BiocParallel::BiocParallelParam` for the per-cluster
#'   plans and per-spot solves, or `NULL` (default) to run serially.
#' @param chunk_size integer. spots per parallel task.
#' @returns a numeric matrix, cell types (rows) by spots (columns).
#' @seealso [dwls_stream()] for expression read in windows.
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:30)
#' sig <- matrix(1, 30, 3, dimnames = list(genes, c("A", "B", "C")))
#' sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
#'
#' # 20 spots mixing A and B, 20 spots of C alone
#' mix <- cbind(rbind(a <- runif(20, 0.3, 0.7), 1 - a, 0),
#'              matrix(c(0, 0, 1), 3, 20))
#' counts <- matrix(rpois(30 * 40, sig %*% mix * 5), 30,
#'                  dimnames = list(genes, paste0("s", 1:40)))
#' clusters <- rep(c("mixed", "pure"), each = 20)
#'
#' prop <- dwls(log2(counts + 1), sig, clusters, n_cell = 5)
#' round(prop[, c(1:3, 21:23)], 2)
#' @export
dwls <- function(expr, sign_matrix, clusters, n_cell = 50, cutoff = 2,
    logbase = 2, na_as_zero = FALSE, BPPARAM = NULL, chunk_size = 256L) {
    spots <- colnames(expr)
    if (is.null(spots)) {
        stop("[dwls] `expr` needs spot IDs as column names.", call. = FALSE)
    }
    if (anyDuplicated(spots)) {
        stop("[dwls] `expr` has duplicated spot IDs (column names), e.g. ",
            spots[anyDuplicated(spots)], ".", call. = FALSE)
    }
    if (is.null(names(clusters))) {
        if (length(clusters) != length(spots)) {
            stop("[dwls] unnamed `clusters` must have one value per ",
                "expression column.", call. = FALSE)
        }
        names(clusters) <- spots
    }
    # One window: the whole matrix, restricted to the signature features when
    # the driver asks for them.
    dwls_stream(
        read_window = function(i, features) expr[features, , drop = FALSE],
        n_windows = 1L,
        feat_ids = rownames(expr),
        clusters = clusters,
        sign_matrix = sign_matrix,
        n_cell = n_cell, cutoff = cutoff, logbase = logbase,
        na_as_zero = na_as_zero, BPPARAM = BPPARAM, chunk_size = chunk_size
    )
}


#' Deconvolve spots read in windows
#'
#' @description
#' The spatialDWLS driver for expression backends that read spots in windows,
#' such as an on-disk store. [dwls()] calls it with a single window.
#'
#' The reader is called once per window in each of four passes: cluster and
#' feature means, PAGE z-scores, the cluster round of solves and the spot
#' round. Only the signature features are ever requested, and no more than one
#' window is held at a time (apart from the cell types x spots result).
#'
#' @param read_window `function(i, features)` returning the log-scaled
#'   expression of window `i` (1..`n_windows`) for exactly `features`, in that
#'   order, as a features x spots `matrix` or `Matrix` with spot IDs as column
#'   names. It must return the same spots for the same `i` on every pass,
#'   in any order.
#' @param n_windows integer. number of windows.
#' @param feat_ids character. every feature the backend holds; the signature
#'   is restricted to these.
#' @param clusters cluster of every spot, named by spot ID. Extra entries are
#'   ignored.
#' @inheritParams dwls
#' @returns a numeric matrix, cell types (rows) by spots (columns), spots in
#'   window order.
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:30)
#' sig <- matrix(1, 30, 3, dimnames = list(genes, c("A", "B", "C")))
#' sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
#' mix <- cbind(rbind(a <- runif(20, 0.3, 0.7), 1 - a, 0),
#'              matrix(c(0, 0, 1), 3, 20))
#' expr <- log2(matrix(rpois(30 * 40, sig %*% mix * 5), 30,
#'                     dimnames = list(genes, paste0("s", 1:40))) + 1)
#' clusters <- setNames(rep(c("mixed", "pure"), each = 20), colnames(expr))
#'
#' # read the spots ten at a time
#' windows <- split(colnames(expr), rep(1:4, each = 10))
#' prop <- dwls_stream(
#'     read_window = function(i, features) expr[features, windows[[i]]],
#'     n_windows = length(windows), feat_ids = rownames(expr),
#'     clusters = clusters, sign_matrix = sig, n_cell = 5
#' )
#' all.equal(prop, dwls(expr, sig, clusters, n_cell = 5))
#' @export
dwls_stream <- function(read_window, n_windows, feat_ids, clusters,
    sign_matrix, n_cell = 50, cutoff = 2, logbase = 2, na_as_zero = FALSE,
    BPPARAM = NULL, chunk_size = 256L) {
    if (!is.function(read_window)) {
        stop("[dwls] `read_window` must be a function(i, features).",
            call. = FALSE)
    }
    if (is.null(names(clusters))) {
        stop("[dwls] `clusters` must be named by spot ID.", call. = FALSE)
    }
    if (length(n_cell) != 1L || !is.finite(n_cell) || n_cell <= 0) {
        stop("[dwls] `n_cell` must be a single positive number.",
            call. = FALSE)
    }
    if (!is.null(BPPARAM) && !methods::is(BPPARAM, "BiocParallelParam")) {
        stop("[dwls] `BPPARAM` must be a BiocParallelParam or NULL.",
            call. = FALSE)
    }
    opts <- list(chunk_size = as.integer(chunk_size), BPPARAM = BPPARAM)

    sig <- as.matrix(sign_matrix)
    .dwls_check_names(colnames(sig), "cell type (column) names of `sign_matrix`")
    .dwls_check_names(rownames(sig), "feature (row) names of `sign_matrix`")
    genes <- intersect(rownames(sig), feat_ids)
    if (length(genes) == 0L) {
        stop("[dwls] no features shared between `sign_matrix` and the ",
            "expression values.", call. = FALSE)
    }
    sig <- sig[genes, , drop = FALSE]
    # Clusters of spots that are never read are dropped after pass 1, so a
    # superset of labels (e.g. the full metadata of a subset) is fine.
    levels <- sort(unique(as.character(clusters[!is.na(clusters)])))

    # Features are restricted *before* the reverse-log transform, which then
    # touches only the non-zeros of a sparse window.
    load_window <- function(i) {
        L <- read_window(i, genes)
        if (!methods::is(L, "Matrix")) L <- as.matrix(L)
        if (!identical(rownames(L), genes)) {
            stop("[dwls] `read_window` must return the requested features, ",
                "in order.", call. = FALSE)
        }
        idx <- match(colnames(L), names(clusters))
        if (anyNA(idx)) {
            stop("[dwls] `clusters` has no entry for ", sum(is.na(idx)),
                " spot(s) of the expression values.", call. = FALSE)
        }
        cl <- clusters[idx]
        if (anyNA(cl)) {
            stop("[dwls] `clusters` is NA for ", sum(is.na(cl)),
                " spot(s) of the expression values; every spot needs a ",
                "cluster.", call. = FALSE)
        }
        list(log = L, nolog = .dwls_unlog(L, logbase),
            clusters = as.character(cl))
    }
    # A single window is read once and reused by every pass.
    cache <- if (n_windows == 1L) load_window(1L) else NULL
    get_window <- function(i) if (is.null(cache)) load_window(i) else cache

    ## pass 1: per-cluster and per-feature means ##
    k_sum <- matrix(0, length(genes), length(levels),
        dimnames = list(genes, levels))
    k_n <- stats::setNames(numeric(length(levels)), levels)
    g_sum <- numeric(length(genes))
    n_spots <- 0L
    for (i in seq_len(n_windows)) {
        w <- get_window(i)
        k_sum <- k_sum + .dwls_cluster_sums(w$nolog, w$clusters, levels)
        k_n <- k_n + tabulate(factor(w$clusters, levels), length(levels))
        g_sum <- g_sum + Matrix::rowSums(w$nolog)
        n_spots <- n_spots + ncol(w$log)
    }
    keep <- k_n > 0
    levels <- levels[keep]
    cluster_mean <- sweep(k_sum[, keep, drop = FALSE], 2L, k_n[keep], "/")
    # PAGE's per-feature reference level: the mean on the reverse-logged
    # scale, logged back
    gene_mean <- log(g_sum / n_spots + 1, base = logbase)

    ## pass 2: round-1 candidates, from the best PAGE z per cluster ##
    sig1 <- sig[rowSums(sig) > 0, , drop = FALSE]
    marker1 <- .dwls_marker_matrix(sig1, min_markers = TRUE)
    max_z <- matrix(-Inf, ncol(sig1), length(levels),
        dimnames = list(colnames(sig1), levels))
    for (i in seq_len(n_windows)) {
        w <- get_window(i)
        z <- .dwls_page_z(w$log, marker1, gene_mean)
        max_z <- pmax(max_z, .dwls_group_rowmax(z, w$clusters, levels))
    }
    plans1 <- .dwls_plan_round1(sig1, marker1, max_z, cluster_mean, cutoff,
        BPPARAM = opts$BPPARAM)

    ## pass 3: cluster round ##
    binary <- vector("list", n_windows)
    present <- matrix(FALSE, ncol(sig1), length(levels),
        dimnames = list(colnames(sig1), levels))
    for (i in seq_len(n_windows)) {
        w <- get_window(i)
        p <- .dwls_solve_round(w$nolog, w$clusters, plans1,
            cell_types = colnames(sig1), round = 1L, opts = opts)
        b <- p >= (1 / n_cell)
        b[is.na(b)] <- FALSE
        present <- present | .dwls_group_rowany(b, w$clusters, levels)
        binary[[i]] <- b
    }
    marker2 <- .dwls_marker_matrix(sig, min_markers = FALSE)
    plans2 <- .dwls_plan_round2(sig, marker2, present, cluster_mean,
        BPPARAM = opts$BPPARAM)

    ## pass 4: spot round, on the types above 1 / n_cell ##
    out <- vector("list", n_windows)
    for (i in seq_len(n_windows)) {
        w <- get_window(i)
        # pass 3's per-spot types are matched to this pass's spots by name: a
        # reader may return a window's spots in another order the second time
        if (!identical(colnames(w$log), colnames(binary[[i]]))) {
            if (!setequal(colnames(w$log), colnames(binary[[i]]))) {
                stop("[dwls] `read_window(", i, ", )` returned different ",
                    "spots on different passes.", call. = FALSE)
            }
            binary[[i]] <- binary[[i]][, colnames(w$log), drop = FALSE]
        }
        out[[i]] <- .dwls_solve_round(w$nolog, w$clusters, plans2,
            cell_types = colnames(sig), round = 2L, binary = binary[[i]],
            marker = marker2, opts = opts)
    }
    res <- do.call(cbind, out)
    if (isTRUE(na_as_zero)) res[is.na(res)] <- 0
    res
}


#' Build a DWLS signature matrix
#'
#' @description
#' Mean expression per cell type over the signature features: the reference
#' [dwls()] expects.
#'
#' @param expr single-cell expression, features (rows) by cells (columns).
#'   Pass the scale the spots will be solved on, i.e. normalised and *not*
#'   logged (`dwls()` reverse-logs the spot expression before solving).
#' @param sign_gene character. signature features, e.g. cell type markers.
#'   Features missing from `expr` are dropped.
#' @param cell_type cell type of each column of `expr`.
#' @returns a numeric matrix, signature features (rows) by cell type
#'   (columns, in order of first appearance in `cell_type`).
#' @examples
#' set.seed(1)
#' x <- matrix(rpois(200, 3), 20, 10,
#'             dimnames = list(paste0("g", 1:20), paste0("c", 1:10)))
#' dwls_signature(x, sign_gene = paste0("g", 1:5),
#'                cell_type = rep(c("A", "B"), each = 5))
#' @export
dwls_signature <- function(expr, sign_gene, cell_type) {
    if (ncol(expr) != length(cell_type)) {
        stop("[dwls_signature] `cell_type` needs one value per column of ",
            "`expr`.", call. = FALSE)
    }
    if (!is.character(sign_gene)) {
        stop("[dwls_signature] `sign_gene` must be a character vector.",
            call. = FALSE)
    }
    genes <- intersect(rownames(expr), sign_gene)
    types <- unique(as.character(cell_type))
    out <- vapply(types, function(ct) {
        Matrix::rowMeans(expr[genes, cell_type == ct, drop = FALSE])
    }, numeric(length(genes)))
    out <- matrix(out, nrow = length(genes), dimnames = list(genes, types))
    out
}


# Names that index results must exist and be unique: a duplicated cell type
# would silently overwrite another type's proportions.
.dwls_check_names <- function(x, what) {
    if (is.null(x) || anyNA(x) || any(!nzchar(x))) {
        stop("[dwls] the ", what, " are missing.", call. = FALSE)
    }
    if (anyDuplicated(x)) {
        stop("[dwls] the ", what, " must be unique; duplicated: ",
            paste(unique(x[duplicated(x)]), collapse = ", "), ".",
            call. = FALSE)
    }
    invisible(x)
}

# logbase^x - 1, keeping a sparse matrix sparse (0 maps to 0).
.dwls_unlog <- function(x, logbase) {
    if (methods::is(x, "sparseMatrix")) {
        x <- methods::as(x, "CsparseMatrix")
        x@x <- logbase^x@x - 1
        return(x)
    }
    logbase^x - 1
}
