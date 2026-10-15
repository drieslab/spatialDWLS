test_that("dwls() returns a composition per spot", {
    d <- .synthetic()
    prop <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)

    expect_true(is.matrix(prop))
    expect_identical(dimnames(prop), list(colnames(d$sig), colnames(d$expr)))
    expect_gte(min(prop), 0)
    expect_equal(unname(colSums(prop)), rep(1, ncol(prop)), tolerance = 1e-8)
})


test_that("clusters are matched by name, not position", {
    d <- .synthetic()
    in_order <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)
    shuffled <- dwls(d$expr, d$sig, rev(d$clusters), n_cell = 5)
    expect_equal(shuffled, in_order)

    # unnamed: positional
    expect_equal(dwls(d$expr, d$sig, unname(d$clusters), n_cell = 5),
        in_order)

    expect_error(dwls(d$expr, d$sig, unname(d$clusters)[-1]), "one value per")
    expect_error(dwls(d$expr, d$sig, d$clusters[-1]), "no entry")
})


test_that("a cluster with a single cell type is assigned to it", {
    d <- .synthetic()
    prop <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)
    pure <- prop[, d$clusters == "pure"]

    # Was all zero in Giotto <= 4.3.1: the round-2 branch for this case
    # compared (`==`) instead of assigning.
    expect_equal(unname(pure["C", ]), rep(1, ncol(pure)))
    expect_equal(unname(colSums(pure[c("A", "B"), ])), rep(0, ncol(pure)))
})


test_that("unsolvable spots are NA, or 0 with na_as_zero", {
    d <- .synthetic()
    d$expr[, "s25"] <- 0

    prop <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)
    expect_true(all(is.na(prop[, "s25"])))
    expect_false(anyNA(prop[, -25]))

    prop0 <- dwls(d$expr, d$sig, d$clusters, n_cell = 5, na_as_zero = TRUE)
    expect_identical(unname(prop0[, "s25"]), c(0, 0, 0))
})


test_that("dwls() leaves the session's RNG alone", {
    d <- .synthetic()

    set.seed(11)
    expected <- stats::runif(1)
    set.seed(11)
    dwls(d$expr, d$sig, d$clusters)
    expect_identical(stats::runif(1), expected)

    # a non-default kind survives, and so does having no seed at all
    old <- RNGkind()
    on.exit(RNGkind(old[[1]], old[[2]], old[[3]]), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    dwls(d$expr, d$sig, d$clusters)
    expect_identical(RNGkind()[[1]], "L'Ecuyer-CMRG")

    if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
    }
    dwls(d$expr, d$sig, d$clusters)
    expect_false(exists(".Random.seed", envir = globalenv(),
        inherits = FALSE))
})


test_that("the result does not depend on the session's RNG kind", {
    d <- .synthetic()
    default <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)

    old <- RNGkind()
    on.exit(RNGkind(old[[1]], old[[2]], old[[3]]), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    expect_equal(dwls(d$expr, d$sig, d$clusters, n_cell = 5), default)
})


test_that("parallel solves match serial ones", {
    skip_on_os("windows")
    skip_if_not_installed("BiocParallel")
    d <- .synthetic()

    serial <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)
    # chunk_size = 7 splits each cluster over several tasks; the workers run
    # L'Ecuyer-CMRG, which the dampening constant's halves must not depend on
    parallel <- dwls(d$expr, d$sig, d$clusters, n_cell = 5, chunk_size = 7L,
        BPPARAM = BiocParallel::MulticoreParam(2))
    expect_equal(parallel, serial)
})


test_that("sparse and dense input give the same result", {
    d <- .synthetic()
    dense <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)
    sparse <- dwls(Matrix::Matrix(d$expr, sparse = TRUE), d$sig, d$clusters,
        n_cell = 5)
    expect_equal(sparse, dense)
})


test_that("arguments are validated", {
    d <- .synthetic()
    expect_error(dwls(d$expr, d$sig, d$clusters, n_cell = 0), "n_cell")
    expect_error(dwls(d$expr, d$sig, d$clusters, n_cell = c(5, 10)), "n_cell")
    expect_error(dwls(d$expr, d$sig, d$clusters, BPPARAM = "two"), "BPPARAM")
    x <- d$expr
    colnames(x) <- NULL
    expect_error(dwls(x, d$sig, d$clusters), "column names")
    sig <- d$sig
    rownames(sig) <- paste0("other", seq_len(nrow(sig)))
    expect_error(dwls(d$expr, sig, d$clusters), "no features shared")
})
