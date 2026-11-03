test_that("dwls_design() holds the reverse-logged signature features and J per spot", {
    s <- .synthetic()
    d <- dwls_design(s$expr, s$sig, s$clusters)
    expect_identical(rownames(d$y), rownames(s$sig))
    expect_identical(colnames(d$y), colnames(s$expr))
    expect_equal(d$y, 2^s$expr - 1)
    expect_identical(d$S, s$sig)
    expect_identical(names(d$J), colnames(s$expr))
    expect_true(all(d$J >= 1 & d$J == round(d$J)))
    # one J per cluster, the cross-validated constant of the cluster mean
    for (cl in unique(s$clusters)) {
        sp <- names(s$clusters)[s$clusters == cl]
        expect_length(unique(d$J[sp]), 1L)
        expect_equal(unname(d$J[sp[1]]),
            as.numeric(.dwls_cluster_constant(s$sig, rowMeans(d$y[, sp]))))
    }
})

test_that("dwls_design() matches clusters by name and restricts to shared features", {
    s <- .synthetic()
    keep <- rownames(s$expr)[-(1:3)]
    shuffled <- rev(s$clusters)
    d <- dwls_design(s$expr[keep, ], s$sig, shuffled)
    expect_identical(rownames(d$y), keep)
    expect_identical(d$clusters, s$clusters)
    expect_error(dwls_design(s$expr, s$sig, s$clusters[1:5]), "no cluster")
    expect_error(dwls_design(unname(s$expr), s$sig, s$clusters), "spot IDs")
})

test_that("dwls_design() keeps sparse expression sparse", {
    s <- .synthetic()
    d <- dwls_design(Matrix::Matrix(s$expr, sparse = TRUE), s$sig, s$clusters)
    expect_s4_class(d$y, "sparseMatrix")
    expect_equal(as.matrix(d$y), 2^s$expr - 1)
})

test_that("the loss weights are spatialDWLS's dampened weights", {
    s <- .synthetic()
    d <- dwls_design(s$expr, s$sig, s$clusters)
    p <- c(A = 0.5, B = 0.3, C = 0.2)
    m <- d$S %*% p
    w <- .dwls_dampened_weight_matrix(m, 2^(d$J[1] - 1), max(d$S))
    expect_equal(as.vector(w), .dwls_dampened_weights(d$S, p, 2^(d$J[1] - 1)))
})

test_that("dwls_expression_loss() prefers the true mixtures and ignores scale", {
    s <- .synthetic()
    d <- dwls_design(s$expr, s$sig, s$clusters)
    n <- ncol(s$expr) / 2
    truth <- matrix(0, 3, 2 * n, dimnames = list(colnames(s$sig), colnames(s$expr)))
    truth["C", 1:n] <- 1
    a <- (2^s$expr[1:10, n + seq_len(n)] - 1)
    b <- (2^s$expr[11:20, n + seq_len(n)] - 1)
    truth["A", n + seq_len(n)] <- colSums(a) / (colSums(a) + colSums(b))
    truth["B", n + seq_len(n)] <- 1 - truth["A", n + seq_len(n)]
    flat <- truth; flat[] <- 1 / 3
    l_true <- dwls_expression_loss(d, truth, grad = FALSE)$value
    expect_lt(l_true, dwls_expression_loss(d, flat, grad = FALSE)$value / 5)
    expect_equal(dwls_expression_loss(d, 3 * truth, grad = FALSE)$value, l_true)
    # rows are matched by name, extra rows get a zero gradient
    shuffled <- rbind(truth[3:1, ], D = 0.1)
    l <- dwls_expression_loss(d, shuffled)
    expect_equal(l$value, dwls_expression_loss(d, truth + 0)$value)
    expect_identical(dimnames(l$grad), dimnames(shuffled))
    expect_true(all(l$grad["D", ] == 0))
})

test_that("dwls_expression_loss() gradient matches finite differences with the weights held", {
    s <- .synthetic()
    d <- dwls_design(s$expr, s$sig, s$clusters)
    set.seed(5)
    sp <- colnames(s$expr)[c(1:4, 25:30)]
    p <- matrix(stats::runif(3 * length(sp)), 3, dimnames = list(colnames(s$sig), sp))
    g <- dwls_expression_loss(d, p)$grad
    # hold the weights at p: same loss with a fixed weight matrix
    w0 <- .dwls_dampened_weight_matrix(d$S %*% p, 2^(d$J[sp] - 1), max(d$S))
    y <- as.matrix(d$y[, sp])
    fixed <- function(q) {
        m <- d$S %*% q
        cc <- colSums(w0 * y * m) / colSums(w0 * m * m)
        mean(colSums(w0 * (y - sweep(m, 2, cc, "*"))^2) / colSums(w0 * y * y))
    }
    expect_equal(fixed(p), dwls_expression_loss(d, p, grad = FALSE)$value)
    h <- 1e-6
    for (k in seq_len(8)) {
        i <- sample(3, 1); j <- sample(length(sp), 1)
        up <- dn <- p; up[i, j] <- up[i, j] + h; dn[i, j] <- dn[i, j] - h
        num <- (fixed(up) - fixed(dn)) / (2 * h)
        expect_equal(g[i, j], num, tolerance = 1e-5)
    }
})

test_that("a spot without signal has loss 0 and gradient 0", {
    s <- .synthetic()
    s$expr[, 1] <- 0
    d <- dwls_design(s$expr, s$sig, s$clusters)
    p <- matrix(1 / 3, 3, 1, dimnames = list(colnames(s$sig), colnames(s$expr)[1]))
    l <- dwls_expression_loss(d, p)
    expect_equal(l$value, 0)
    expect_true(all(l$grad == 0))
    expect_error(dwls_expression_loss(d, unname(p)), "cell types")
})
