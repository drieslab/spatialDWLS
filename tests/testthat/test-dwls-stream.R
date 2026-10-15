test_that("dwls_stream() gives the same answer over several windows", {
    d <- .synthetic(n_per = 30L)
    spots <- colnames(d$expr)
    # uneven windows that cut through both clusters
    cuts <- split(sample(spots), c(rep(1, 7), rep(2, 25), rep(3, 28)))
    reader <- function(i, features) d$expr[features, cuts[[i]], drop = FALSE]

    one <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)
    many <- dwls_stream(read_window = reader, n_windows = length(cuts),
        feat_ids = rownames(d$expr), clusters = d$clusters,
        sign_matrix = d$sig, n_cell = 5)

    expect_setequal(colnames(many), spots)
    expect_equal(many[, spots], one[, spots], tolerance = 1e-10)
})


test_that("dwls_stream() checks what the reader returns", {
    d <- .synthetic()
    expect_error(
        dwls_stream(
            read_window = function(i, features) d$expr[rev(features), ],
            n_windows = 1L, feat_ids = rownames(d$expr),
            clusters = d$clusters, sign_matrix = d$sig),
        "requested features"
    )
    expect_error(
        dwls_stream(read_window = "not a function", n_windows = 1L,
            feat_ids = rownames(d$expr), clusters = d$clusters,
            sign_matrix = d$sig),
        "must be a function"
    )
    expect_error(
        dwls_stream(
            read_window = function(i, features) d$expr[features, ],
            n_windows = 1L, feat_ids = rownames(d$expr),
            clusters = unname(d$clusters), sign_matrix = d$sig),
        "named by spot ID"
    )
})


test_that("a window may return its spots in a different order each pass", {
    d <- .synthetic(n_per = 30L)
    spots <- colnames(d$expr)
    cuts <- split(spots, rep(1:3, each = 20))
    calls <- 0L
    # every other read of a window comes back reversed
    reader <- function(i, features) {
        calls <<- calls + 1L
        s <- cuts[[i]]
        if (calls %% 2L == 0L) s <- rev(s)
        d$expr[features, s, drop = FALSE]
    }
    one <- dwls(d$expr, d$sig, d$clusters, n_cell = 5)
    many <- dwls_stream(read_window = reader, n_windows = length(cuts),
        feat_ids = rownames(d$expr), clusters = d$clusters,
        sign_matrix = d$sig, n_cell = 5)
    expect_equal(many[, spots], one[, spots], tolerance = 1e-10)

    # different spots on the last pass than on the one before is an error,
    # not a silent mix-up (3 windows x 4 passes: calls 10-12 are the last)
    calls <- 0L
    bad <- function(i, features) {
        calls <<- calls + 1L
        s <- if (calls > 9L) cuts[[(i %% 3L) + 1L]] else cuts[[i]]
        d$expr[features, s, drop = FALSE]
    }
    expect_error(dwls_stream(read_window = bad, n_windows = 3L,
        feat_ids = rownames(d$expr), clusters = d$clusters,
        sign_matrix = d$sig, n_cell = 5), "different spots")
})
