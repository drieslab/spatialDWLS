# Per-spot solves and the dampened WLS solver. Internal.

# per-spot solves ####

# Solves one round for the spots in `nolog_expr`, which may be any subset of
# spots as long as `clusters` (and `binary`, round 2) line up with its
# columns. Returns cell types x spots.
.dwls_solve_round <- function(nolog_expr, clusters, plans, cell_types, round,
    binary = NULL, marker = NULL, opts) {
    out <- matrix(0, length(cell_types), ncol(nolog_expr),
        dimnames = list(cell_types, colnames(nolog_expr)))

    # Each task carries its own slice of the expression, so a parallel worker
    # is sent the features and spots it solves and not the whole matrix. The
    # slice is densified once per cluster, then cut into tasks.
    tasks <- .dwls_tasks(clusters, opts$chunk_size)
    for (cl in unique(vapply(tasks, `[[`, "", "cluster"))) {
        in_cl <- which(vapply(tasks, function(t) t$cluster == cl, TRUE))
        idx <- unlist(lapply(tasks[in_cl], `[[`, "idx"))
        dense <- as.matrix(nolog_expr[plans[[cl]]$genes, idx, drop = FALSE])
        colnames(dense) <- as.character(idx)
        for (i in in_cl) {
            sel <- as.character(tasks[[i]]$idx)
            tasks[[i]]$B <- dense[, sel, drop = FALSE]
            if (round == 2L) {
                tasks[[i]]$binary <- binary[, tasks[[i]]$idx, drop = FALSE]
            }
        }
    }
    solve_task <- function(task) {
        plan <- plans[[task$cluster]]
        if (round == 1L) {
            .dwls_task_round1(plan, task$B)
        } else {
            .dwls_task_round2(plan, task$B, task$binary, marker, cell_types)
        }
    }
    res <- .dwls_lapply(tasks, solve_task, opts$BPPARAM)

    for (i in seq_along(tasks)) {
        r <- res[[i]]
        out[rownames(r), tasks[[i]]$idx] <- r
    }
    # Clip the constrained solver's round-off. Every row: the loop this
    # replaces ran `for (i in nrow)` and only ever touched the last type.
    out[] <- pmax(out, 0)
    out
}

.dwls_task_round1 <- function(plan, B) {
    r <- matrix(0, length(plan$ct), ncol(B),
        dimnames = list(plan$ct, NULL))
    for (k in seq_len(ncol(B))) {
        if (sum(B[, k]) > 0) {
            r[, k] <- .dwls_spot(plan$S, B[, k], plan$J)
        }
    }
    r
}

.dwls_task_round2 <- function(plan, B, binary, marker, cell_types) {
    r <- matrix(0, length(cell_types), ncol(binary),
        dimnames = list(cell_types, NULL))

    if (length(plan$ct) == 0L) {
        r[] <- NA_real_
        return(r)
    }
    if (length(plan$ct) == 1L) {
        # Only one type anywhere in the cluster: every spot is that type.
        # This branch used to read `== 1` and left the cluster all zero.
        r[plan$ct, ] <- 1
        return(r)
    }

    # Spots share a handful of type combinations; look each one's markers up
    # once rather than once per spot.
    genes_for <- new.env(parent = emptyenv())
    for (k in seq_len(ncol(binary))) {
        b <- B[, k]
        ct_k <- rownames(binary)[binary[, k]]
        if (sum(b) == 0 || length(ct_k) == 0L) {
            r[, k] <- NA_real_
            next
        }
        if (length(ct_k) == 1L) {
            r[ct_k, k] <- 1
            next
        }
        key <- paste(ct_k, collapse = "\r")
        genes_k <- genes_for[[key]]
        if (is.null(genes_k)) {
            genes_k <- .dwls_marker_genes(marker, ct_k, rownames(marker))
            genes_for[[key]] <- genes_k
        }
        if (sum(b[genes_k]) == 0) {
            r[, k] <- NA_real_
            next
        }
        # genes_k and ct_k are subsets of the cluster's, so this is a slice
        # of the plan's signature, not a fresh lookup
        sol <- .dwls_spot(plan$S[genes_k, ct_k, drop = FALSE], b[genes_k],
            plan$J)
        r[names(sol), k] <- sol
    }
    r
}

# Chunks of spots within a cluster: the unit of parallel work.
.dwls_tasks <- function(clusters, chunk_size = 256L) {
    by_cluster <- split(seq_along(clusters), clusters)
    tasks <- list()
    for (cl in names(by_cluster)) {
        idx <- by_cluster[[cl]]
        for (part in split(idx, ceiling(seq_along(idx) / chunk_size))) {
            tasks[[length(tasks) + 1L]] <- list(cluster = cl, idx = part)
        }
    }
    tasks
}

# Serial unless a parallel BiocParallel param is passed. BiocParallel is only
# suggested, so it is looked up when, and only when, it is asked for.
.dwls_lapply <- function(X, FUN, BPPARAM = NULL) {
    if (is.null(BPPARAM) || methods::is(BPPARAM, "SerialParam") ||
        length(X) < 2L) {
        return(lapply(X, FUN))
    }
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        stop("[dwls] `BPPARAM` needs the BiocParallel package.",
            call. = FALSE)
    }
    BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM)
}


# solver ####

# Dampened weighted least squares for one spot (Tsoucas et al. 2019): start
# from the non-negative OLS solution, then iterate weights 1 / (S x)^2,
# capped at 2^(J - 1) times the smallest, averaging each step 1:4 with the
# previous solution until it moves less than 0.01.
.dwls_spot <- function(S, B, J) {
    solution <- .dwls_ols(S, B)
    change <- 1
    iterations <- 0L
    while (change > 0.01 && iterations < 1000L) {
        new <- .dwls_wls(S, B, solution, J)
        averaged <- (new + 4 * solution) / 5
        change <- sum(abs(averaged - solution))
        solution <- averaged
        iterations <- iterations + 1L
    }
    solution / sum(solution)
}

# Non-negative least squares as a QP.
.dwls_ols <- function(S, B) {
    k <- ncol(S)
    D <- crossprod(S)
    d <- crossprod(S, B)
    out <- .dwls_qp(D, d, k)
    names(out) <- colnames(S)
    out
}

# One dampened WLS step. The weights enter as `crossprod(S * w, S)`; the
# previous version built diag(w), a dense features x features matrix, on
# every iteration of every spot.
.dwls_wls <- function(S, B, gold, J) {
    w <- .dwls_dampened_weights(S, gold, 2^(J - 1))
    Sw <- S * w
    D <- crossprod(Sw, S)
    d <- crossprod(Sw, B)
    D <- (D + t(D)) / 2
    # one eigendecomposition serves both the scaling (D is PSD, so its
    # largest eigenvalue is its 2-norm) and the positive-definiteness check
    ev <- eigen(D, symmetric = TRUE, only.values = TRUE)$values
    sc <- max(abs(ev))
    out <- .dwls_qp(.dwls_near_pd(D / sc, ev / sc), d / sc, ncol(S))
    names(out) <- colnames(S)
    out
}

.dwls_dampened_weights <- function(S, gold, multiplier) {
    ws <- as.vector((1 / (S %*% gold))^2)
    ws <- ws / min(ws)
    ws[which(ws > multiplier)] <- multiplier
    ws
}

# Matrix::nearPD(), but only when it would change something. nearPD() drops
# eigenvalues below 1e-6 of the largest and then floors them at 1e-8; a
# matrix already above the first threshold comes back as it went in. That is
# nearly every one here, and nearPD() was most of the per-iteration cost.
# `ev` are the eigenvalues of symmetric `D`, when the caller already has them.
.dwls_near_pd <- function(D, ev = NULL) {
    if (is.null(ev)) {
        D <- (D + t(D)) / 2
        ev <- eigen(D, symmetric = TRUE, only.values = TRUE)$values
    }
    if (min(ev) > 1e-6 * max(ev)) return(D)
    as.matrix(Matrix::nearPD(D)$mat)
}

.dwls_qp <- function(D, d, k) {
    A <- diag(k)
    b0 <- rep(0, k)
    out <- tryCatch(
        quadprog::solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b0)$solution,
        error = function(e) NULL
    )
    if (is.null(out)) {
        D <- as.matrix(Matrix::nearPD(D)$mat)
        out <- tryCatch(
            quadprog::solve.QP(
                Dmat = D, dvec = d, Amat = A, bvec = b0
            )$solution,
            error = function(e) {
                stop("[DWLS] quadratic program could not be solved, even ",
                    "after Matrix::nearPD(): ", conditionMessage(e),
                    call. = FALSE)
            }
        )
    }
    out
}

# Choose the dampening constant J by cross-validation on the cluster mean:
# for each candidate cap 2^(j - 1), refit on 100 random halves of the
# features and keep the j whose solutions vary least.
#
# Same halves and the same answer as before, but each refit is a k x k weighted
# least squares through .lm.fit() rather than stats::lm() with its formula and
# model frame, and the halves are drawn once instead of once per candidate.
# This was ~90% of the run time.
.dwls_dampening_constant <- function(S, B, gold) {
    n <- nrow(S)
    k <- ncol(S)
    ws <- as.vector((1 / (S %*% gold))^2)
    ws <- ws / min(ws)
    finite <- ws[is.finite(ws)]
    n_j <- if (length(finite)) ceiling(log2(max(finite))) else 0
    if (!is.finite(n_j) || n_j < 1) return(1L)

    halves <- .dwls_cv_halves(n)
    scale <- sum(gold)
    sd2 <- numeric(n_j)
    for (j in seq_len(n_j)) {
        w <- ws
        w[which(ws > 2^(j - 1))] <- 2^(j - 1)
        sw <- sqrt(w)
        Sw <- S * sw
        Bw <- B * sw
        sols <- vapply(halves, function(h) {
            fit <- stats::.lm.fit(Sw[h, , drop = FALSE], Bw[h])
            if (fit$rank < k) return(rep(NA_real_, k))
            coef <- fit$coefficients
            coef[fit$pivot] <- coef
            coef * scale / sum(coef)
        }, numeric(k))
        sd2[j] <- mean(apply(sols, 1L, stats::sd)^2)
    }
    # all NA only when every half was rank deficient; cap at the first level
    if (all(is.na(sd2))) return(1L)
    which.min(sd2)
}

# The 100 random halves, seeded 1..100 so the choice of J is reproducible.
#
# The caller's RNG is left exactly as it was: its kind and its seed (or its
# absence) are restored on exit. The kind is pinned to R's defaults while the
# halves are drawn because parallel workers run L'Ecuyer-CMRG, under which
# set.seed(i) draws different halves, and a cluster planned on a worker would
# otherwise get a different J than in the session.
.dwls_cv_halves <- function(n, n_rep = 100L) {
    .dwls_local_rng()
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    lapply(seq_len(n_rep), function(i) {
        set.seed(i)
        sample(n, size = n * 0.5)
    })
}

# Restore the RNG kind and `.Random.seed` when the *calling* function exits.
# The kind is restored first: changing it reseeds, which the seed restore then
# overwrites.
.dwls_local_rng <- function(envir = parent.frame()) {
    genv <- globalenv()
    kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
    seed <- if (had_seed) get(".Random.seed", envir = genv) else NULL
    restore <- function() {
        RNGkind(kind[[1L]], kind[[2L]], kind[[3L]])
        if (had_seed) {
            assign(".Random.seed", seed, envir = genv)
        } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
            rm(".Random.seed", envir = genv)
        }
    }
    # the call carries the closure itself, so it runs in the caller's frame
    # without `restore` having to be visible there
    do.call(on.exit, list(as.call(list(restore)), add = TRUE), envir = envir)
    invisible()
}
