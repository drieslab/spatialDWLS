# Per-cluster summaries and plans: everything a spot needs that is not in its
# own column. Internal; see dwls.R for how the passes use them.

# per-cluster summaries ####
#
# All three are reductions over spots that combine across windows: sums add,
# maxima take pmax, presence ORs. `levels` fixes the cluster columns so a
# window missing a cluster still lines up.

# (Reverse-logged) expression summed per feature per cluster.
.dwls_cluster_sums <- function(nolog_expr, clusters, levels) {
    f <- factor(clusters, levels = levels)
    ind <- Matrix::sparseMatrix(
        i = seq_along(f), j = as.integer(f), x = 1,
        dims = c(length(f), length(levels)),
        dimnames = list(NULL, levels)
    )
    as.matrix(nolog_expr %*% ind)
}

# Which features mark which cell type: each feature goes to the type where its
# mean expression is highest. Round 1 additionally gives any type left with two
# or fewer markers its top two features.
.dwls_marker_matrix <- function(sig, min_markers = FALSE) {
    m <- matrix(0, nrow(sig), ncol(sig), dimnames = dimnames(sig))
    m[cbind(seq_len(nrow(sig)), max.col(sig, ties.method = "first"))] <- 1
    if (isTRUE(min_markers)) {
        for (i in which(colSums(m) <= 2)) {
            m[rank(-sig[, i]) <= 2, i] <- 1
        }
    }
    m
}

# PAGE z-scores, cell types x spots, for any subset of spots.
#
# Same statistic as the old element-wise loop: within each spot, the mean
# fold change of a type's markers against the mean and sd of all features'
# fold changes, scaled by sqrt(markers). Built from column sums and one
# cross-product instead of materialising the dense fold-change matrix.
.dwls_page_z <- function(log_expr, marker, gene_mean) {
    n_g <- nrow(log_expr)
    col_sum <- Matrix::colSums(log_expr)
    col_mean <- col_sum / n_g - mean(gene_mean)
    cross_mean <- as.vector(Matrix::crossprod(gene_mean, log_expr))
    col_ss <- Matrix::colSums(log_expr^2) - 2 * cross_mean + sum(gene_mean^2)
    col_sd <- sqrt(pmax(col_ss - n_g * col_mean^2, 0) / (n_g - 1))

    # Every feature counts toward the spot's mean and sd, markers or not.
    full <- matrix(0, n_g, ncol(marker),
        dimnames = list(rownames(log_expr), colnames(marker)))
    hit <- intersect(rownames(marker), rownames(log_expr))
    full[hit, ] <- marker[hit, , drop = FALSE]
    marker <- full
    n_m <- colSums(marker)
    sig_mean <- as.matrix(Matrix::crossprod(marker, log_expr))
    sig_mean <- (sig_mean - as.vector(crossprod(marker, gene_mean))) / n_m

    z <- sweep(sig_mean, 2L, col_mean, "-") * sqrt(n_m)
    z <- sweep(z, 2L, col_sd, "/")
    dimnames(z) <- list(colnames(marker), colnames(log_expr))
    z
}

.dwls_group_rowmax <- function(x, clusters, levels) {
    out <- matrix(-Inf, nrow(x), length(levels),
        dimnames = list(rownames(x), levels))
    for (cl in intersect(levels, clusters)) {
        out[, cl] <- suppressWarnings(MatrixGenerics::rowMaxs(
            x[, clusters == cl, drop = FALSE], na.rm = TRUE
        ))
    }
    out
}

.dwls_group_rowany <- function(x, clusters, levels) {
    out <- matrix(FALSE, nrow(x), length(levels),
        dimnames = list(rownames(x), levels))
    for (cl in intersect(levels, clusters)) {
        out[, cl] <- rowSums(x[, clusters == cl, drop = FALSE]) > 0
    }
    out
}


# per-cluster plans ####

# A plan is everything a spot of that cluster needs to be solved: which cell
# types are in play, on which features, and the dampening constant chosen on
# the cluster mean. The constant's cross-validation is the costly part, so
# clusters are planned in parallel like the spots.

.dwls_plan_round1 <- function(sig1, marker1, max_z, cluster_mean, cutoff,
    BPPARAM = NULL) {
    plans <- .dwls_lapply(colnames(max_z), function(cl) {
        row_max <- max_z[, cl]
        names(row_max) <- colnames(sig1)
        ct <- names(row_max)[which(row_max > cutoff)]
        if (length(ct) < 2L) {
            second <- sort(row_max, decreasing = TRUE)[2L]
            ct <- names(row_max)[which(row_max >= second)]
        }
        genes <- .dwls_marker_genes(marker1, ct, rownames(sig1))
        S <- sig1[genes, ct, drop = FALSE]
        list(ct = ct, genes = genes, S = S,
            J = .dwls_cluster_constant(S, cluster_mean[genes, cl]))
    }, BPPARAM)
    names(plans) <- colnames(max_z)
    plans
}

.dwls_plan_round2 <- function(sig, marker2, present, cluster_mean,
    BPPARAM = NULL) {
    plans <- .dwls_lapply(colnames(present), function(cl) {
        ct <- rownames(present)[present[, cl]]
        if (length(ct) <= 1L) {
            return(list(ct = ct, genes = character(), S = NULL, J = NA))
        }
        genes <- .dwls_marker_genes(marker2, ct, rownames(sig))
        S <- sig[genes, ct, drop = FALSE]
        list(ct = ct, genes = genes, S = S,
            J = .dwls_cluster_constant(S, cluster_mean[genes, cl]))
    }, BPPARAM)
    names(plans) <- colnames(present)
    plans
}

# Features marking any of `ct`, in the order of `all_genes`.
.dwls_marker_genes <- function(marker, ct, all_genes) {
    hit <- rowSums(marker[, ct, drop = FALSE]) > 0
    intersect(all_genes, rownames(marker)[hit])
}

.dwls_cluster_constant <- function(S, B) {
    .dwls_dampening_constant(S, B, .dwls_ols(S, B))
}
