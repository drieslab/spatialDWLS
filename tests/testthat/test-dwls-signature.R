test_that("dwls_signature() averages each type over the signature features", {
    x <- matrix(seq_len(40), 4, 10,
        dimnames = list(paste0("g", 1:4), paste0("c", 1:10)))
    ct <- rep(c("B", "A"), each = 5)

    sig <- dwls_signature(x, sign_gene = c("g3", "g1", "missing"),
        cell_type = ct)
    expect_identical(dimnames(sig), list(c("g1", "g3"), c("B", "A")))
    expect_equal(sig[, "B"], rowMeans(x[c("g1", "g3"), 1:5]))
    expect_equal(sig[, "A"], rowMeans(x[c("g1", "g3"), 6:10]))

    # sparse input
    expect_equal(dwls_signature(Matrix::Matrix(x, sparse = TRUE),
        c("g3", "g1"), ct), sig)

    expect_error(dwls_signature(x, "g1", ct[-1]), "one value per column")
    expect_error(dwls_signature(x, 1, ct), "character")
})
