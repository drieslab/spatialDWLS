# Synthetic Visium-like data with known composition.
#  - K cell types over G genes; each type has `n_marker` up-regulated markers
#  - reference single cells -> signature (mean of normalised, not logged,
#    expression: the scale dwls() solves spots on)
#  - spots: niches with Dirichlet-like compositions, 2..max_cells cells per
#    spot, Poisson counts, library-normalised and log2(x + 1)
#  - the niche is the cluster; `n_pure` niches hold a single cell type
simulate_dwls <- function(n_spots = 4000L, G = 2000L, K = 12L, n_marker = 50L,
                          n_niche = 10L, max_cells = 10L, depth = 0.5,
                          seed = 42L, n_pure = 0L) {
    set.seed(seed)
    base <- stats::rlnorm(G, 0, 1)
    mu <- matrix(base, G, K)
    markers <- split(sample(G, K * n_marker), rep(seq_len(K), each = n_marker))
    for (k in seq_len(K)) {
        mu[markers[[k]], k] <- mu[markers[[k]], k] * stats::runif(n_marker, 4, 10)
    }
    mu <- sweep(mu, 2, colSums(mu), "/") * 5000
    rownames(mu) <- paste0("g", seq_len(G))
    colnames(mu) <- paste0("ct", seq_len(K))

    norm <- function(cnt) {
        cnt <- Matrix::Matrix(cnt, sparse = TRUE)
        lib <- Matrix::colSums(cnt)
        lib[lib == 0] <- 1
        Matrix::t(Matrix::t(cnt) / lib) * 1e4
    }
    ref_ct <- rep(colnames(mu), each = 80)
    ref <- matrix(stats::rpois(G * length(ref_ct), mu[, ref_ct] * 0.3), G,
        dimnames = list(rownames(mu), paste0("ref", seq_along(ref_ct))))
    sig <- spatialDWLS::dwls_signature(norm(ref),
        sign_gene = rownames(mu)[unlist(markers)], cell_type = ref_ct)

    niche_p <- t(vapply(seq_len(n_niche), function(i) {
        a <- rep(0.05, K)
        a[sample(K, 4)] <- c(4, 2, 1, 0.5)
        p <- stats::rgamma(K, a)
        p / sum(p)
    }, numeric(K)))
    for (i in seq_len(n_pure)) {
        niche_p[i, ] <- 0
        niche_p[i, i] <- 1
    }
    niche <- sample(n_niche, n_spots, replace = TRUE)
    n_c <- sample(2:max_cells, n_spots, replace = TRUE)
    comp <- t(vapply(seq_len(n_spots), function(s) {
        stats::rmultinom(1, n_c[s], niche_p[niche[s], ])
    }, numeric(K)))
    lambda <- mu %*% t(comp) * depth
    ids <- paste0("spot", seq_len(n_spots))
    cnt <- matrix(stats::rpois(length(lambda), lambda), G,
        dimnames = list(rownames(mu), ids))
    expr <- norm(cnt)
    expr@x <- log2(expr@x + 1)

    truth <- t(comp / n_c)
    dimnames(truth) <- list(colnames(mu), ids)
    list(
        expr = expr, counts = cnt, sig = sig, truth = truth,
        clusters = stats::setNames(paste0("n", niche), ids),
        n_cells = stats::setNames(n_c, ids)
    )
}
