#' The dampened WLS problem of each spot
#' @description
#' The pieces of the spatialDWLS fit as data: the spot expression on the
#' signature features as the solver sees it (reverse-logged), the signature on
#' those features, and the dampening constant `J` of each spot. Use it with
#' [dwls_expression_loss()] to score cell type proportions from any source,
#' e.g. a model that predicts the cell type of each cell in a spot, against
#' the spot's expression with spatialDWLS's weights.
#'
#' `J` is chosen per cluster as the first round of [dwls()] does it: by
#' cross-validation on the cluster's mean expression, against all cell types.
#' Every spot of a cluster gets that cluster's `J`.
#' @inheritParams dwls
#' @returns a list:
#' \describe{
#'   \item{`y`}{reverse-logged expression, signature features by spots. Sparse
#'     input stays sparse.}
#'   \item{`S`}{the signature on those features, features by cell types.}
#'   \item{`J`}{the dampening constant of each spot, named by spot ID. The
#'     weights are capped at `2^(J - 1)`.}
#'   \item{`clusters`}{the cluster of each spot, named by spot ID.}
#' }
#' @seealso [dwls_expression_loss()]
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:30)
#' sig <- matrix(1, 30, 3, dimnames = list(genes, c("A", "B", "C")))
#' sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
#' mix <- cbind(rbind(a <- runif(20, 0.3, 0.7), 1 - a, 0),
#'              matrix(c(0, 0, 1), 3, 20))
#' counts <- matrix(rpois(30 * 40, sig %*% mix * 5), 30,
#'                  dimnames = list(genes, paste0("s", 1:40)))
#' clusters <- rep(c("mixed", "pure"), each = 20)
#' d <- dwls_design(log2(counts + 1), sig, clusters)
#' table(d$J)
#' @export
dwls_design <- function(expr, sign_matrix, clusters, logbase = 2) {
    spots <- colnames(expr)
    if (is.null(spots)) {
        stop("[dwls_design] `expr` needs spot IDs as column names.",
            call. = FALSE)
    }
    .dwls_check_names(spots, "spot IDs (column names of `expr`)")
    .dwls_check_names(colnames(sign_matrix), "cell types (column names of `sign_matrix`)")
    if (is.null(names(clusters))) {
        if (length(clusters) != length(spots)) {
            stop("[dwls_design] unnamed `clusters` must have one value per ",
                "expression column.", call. = FALSE)
        }
        names(clusters) <- spots
    }
    if (!all(spots %in% names(clusters))) {
        stop("[dwls_design] `clusters` has no cluster for some spots, e.g. ",
            setdiff(spots, names(clusters))[1], ".", call. = FALSE)
    }
    clusters <- clusters[spots]
    if (anyNA(clusters)) {
        stop("[dwls_design] `clusters` is NA for some spots.", call. = FALSE)
    }
    feats <- intersect(rownames(sign_matrix), rownames(expr))
    if (length(feats) < ncol(sign_matrix)) {
        stop("[dwls_design] fewer signature features in `expr` (",
            length(feats), ") than cell types (", ncol(sign_matrix), ").",
            call. = FALSE)
    }
    S <- as.matrix(sign_matrix[feats, , drop = FALSE])
    y <- .dwls_unlog(expr[feats, , drop = FALSE], logbase)

    lev <- sort(unique(as.character(clusters)))
    cmean <- sweep(.dwls_cluster_sums(y, clusters, lev), 2,
        tabulate(factor(clusters, levels = lev), length(lev)), "/")
    J <- vapply(lev, function(cl) {
        as.numeric(.dwls_cluster_constant(S, cmean[, cl]))
    }, numeric(1))
    list(
        y = y,
        S = S,
        J = stats::setNames(unname(J[as.character(clusters)]), spots),
        clusters = clusters
    )
}


#' Dampened WLS loss of cell type proportions against spot expression
#' @description
#' How well cell type proportions explain each spot's expression, measured the
#' way spatialDWLS fits it: the spot is predicted as `c * S %*% p`, with the
#' scale `c` solved in closed form, and the residual is weighted with the
#' dampened weights `(1 / S %*% p)^2`, divided by their minimum and capped at
#' `2^(J - 1)`. The weighted residual is divided by the weighted total of the
#' spot, so spots of any depth count alike, and averaged over spots.
#'
#' The gradient treats the weights as fixed, as the iterations of
#' spatialDWLS do; the scale `c` needs no term of its own because it is the
#' optimum for the current proportions. This makes the loss usable for
#' training a model whose predictions are pooled into spot proportions.
#' @param design a [dwls_design()] result.
#' @param prop cell type proportions, cell types (rows, matched to
#'   `colnames(design$S)` by name) by spots (columns, matched to the spot IDs
#'   of `design` by name). Columns need not sum to one.
#' @param grad logical. also return the gradient.
#' @returns a list: `value`, the mean loss over the spots of `prop`, and with
#'   `grad = TRUE` also `grad`, its gradient with respect to `prop` (same
#'   shape and names as `prop`). A spot with no expression on the signature
#'   features has loss 0 and gradient 0.
#' @seealso [dwls_design()]
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:30)
#' sig <- matrix(1, 30, 3, dimnames = list(genes, c("A", "B", "C")))
#' sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
#' mix <- cbind(rbind(a <- runif(20, 0.3, 0.7), 1 - a, 0),
#'              matrix(c(0, 0, 1), 3, 20))
#' counts <- matrix(rpois(30 * 40, sig %*% mix * 5), 30,
#'                  dimnames = list(genes, paste0("s", 1:40)))
#' d <- dwls_design(log2(counts + 1), sig, rep(c("mixed", "pure"), each = 20))
#' dimnames(mix) <- list(colnames(sig), colnames(counts))
#' flat <- mix; flat[] <- 1 / 3
#' dwls_expression_loss(d, mix, grad = FALSE)$value    # the true mixtures
#' dwls_expression_loss(d, flat, grad = FALSE)$value   # much worse
#' @export
dwls_expression_loss <- function(design, prop, grad = TRUE) {
    S <- design$S
    if (is.null(rownames(prop)) || !all(colnames(S) %in% rownames(prop))) {
        stop("[dwls_expression_loss] `prop` needs the cell types of the ",
            "design as row names.", call. = FALSE)
    }
    spots <- colnames(prop)
    if (is.null(spots) || !all(spots %in% colnames(design$y))) {
        stop("[dwls_expression_loss] `prop` needs spot IDs of the design ",
            "as column names.", call. = FALSE)
    }
    p <- as.matrix(prop[colnames(S), , drop = FALSE])
    y <- as.matrix(design$y[, spots, drop = FALSE])
    m <- S %*% p                                          # features x spots
    w <- .dwls_dampened_weight_matrix(m, 2^(design$J[spots] - 1), max(S))
    cc <- colSums(w * y * m) / pmax(colSums(w * m * m), 1e-12)
    total <- pmax(colSums(w * y * y), 1e-12)
    r <- y - sweep(m, 2, cc, "*")
    value <- mean(colSums(w * r * r) / total)
    if (!isTRUE(grad)) return(list(value = value))
    g <- crossprod(S, sweep(w * r, 2, -2 * cc / total, "*")) / length(spots)
    if (nrow(prop) > ncol(S) || !identical(rownames(prop), colnames(S))) {
        full <- matrix(0, nrow(prop), ncol(prop), dimnames = dimnames(prop))
        full[colnames(S), ] <- g
        g <- full
    } else {
        dimnames(g) <- dimnames(prop)
    }
    list(value = value, grad = g)
}

# .dwls_dampened_weights() for every spot at once: column j holds the weights
# of spot j for its prediction m[, j], capped at cap[j]. Predictions are
# floored at 1e-6 of the largest signature value so a feature no proposed
# type expresses gets a large but finite weight.
.dwls_dampened_weight_matrix <- function(m, cap, s_max) {
    w <- (1 / pmax(m, 1e-6 * s_max))^2
    w <- sweep(w, 2, apply(w, 2, min), "/")
    sweep(w, 2, cap, pmin)
}
