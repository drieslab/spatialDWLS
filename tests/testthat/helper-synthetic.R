# Small synthetic problem: three cell types with disjoint markers, a cluster of
# spots that are pure type C and a cluster of A/B mixtures.
#
# Pure C rather than pure A on purpose. Round 1 picks a cluster's candidate
# types by PAGE z-score, which is relative to the whole dataset: had half the
# spots been pure A, A's markers would sit *below* their dataset mean in the
# mixtures and A would not be a candidate there at all.
.synthetic <- function(n_per = 20L, seed = 3L) {
    set.seed(seed)
    genes <- paste0("g", 1:30)
    sig <- matrix(1, 30, 3, dimnames = list(genes, c("A", "B", "C")))
    sig[1:10, "A"] <- 40
    sig[11:20, "B"] <- 40
    sig[21:30, "C"] <- 40
    mix <- cbind(
        matrix(rep(c(0, 0, 1), n_per), 3),
        rbind(stats::runif(n_per, 0.3, 0.7), 0, 0)
    )
    mix[2, (n_per + 1):(2 * n_per)] <- 1 - mix[1, (n_per + 1):(2 * n_per)]
    counts <- sig %*% mix * 5
    counts[] <- stats::rpois(length(counts), counts)
    ids <- paste0("s", seq_len(ncol(counts)))
    dimnames(counts) <- list(genes, ids)
    list(
        expr = log2(counts + 1),
        sig = sig,
        clusters = stats::setNames(rep(c("pure", "mixed"), each = n_per), ids)
    )
}
