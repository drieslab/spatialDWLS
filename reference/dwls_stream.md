# Deconvolve spots read in windows

The spatialDWLS driver for expression backends that read spots in
windows, such as an on-disk store.
[`dwls()`](https://drieslab.github.io/spatialDWLS/reference/dwls.md)
calls it with a single window.

The reader is called once per window in each of four passes: cluster and
feature means, PAGE z-scores, the cluster round of solves and the spot
round. Only the signature features are ever requested, and no more than
one window is held at a time (apart from the cell types x spots result).

## Usage

``` r
dwls_stream(
  read_window,
  n_windows,
  feat_ids,
  clusters,
  sign_matrix,
  n_cell = 50,
  cutoff = 2,
  logbase = 2,
  na_as_zero = FALSE,
  BPPARAM = NULL,
  chunk_size = 256L
)
```

## Arguments

- read_window:

  `function(i, features)` returning the log-scaled expression of window
  `i` (1..`n_windows`) for exactly `features`, in that order, as a
  features x spots `matrix` or `Matrix` with spot IDs as column names.
  It must return the same spots for the same `i` on every pass.

- n_windows:

  integer. number of windows.

- feat_ids:

  character. every feature the backend holds; the signature is
  restricted to these.

- clusters:

  cluster of every spot, named by spot ID.

- sign_matrix:

  mean expression per cell type: signature features (rows) by cell type
  (columns), on the reverse-logged scale, e.g. from
  [`dwls_signature()`](https://drieslab.github.io/spatialDWLS/reference/dwls_signature.md).
  Not a binary marker matrix.

- n_cell:

  numeric. expected number of cells per spot. In round 2 a cell type is
  kept for a spot only when its round-1 proportion is at least
  `1 / n_cell`.

- cutoff:

  numeric. PAGE enrichment z-score a cell type must reach in at least
  one spot of a cluster to be a candidate for that cluster. The two
  highest-scoring types are always candidates.

- logbase:

  numeric. log base of `expr`.

- na_as_zero:

  logical. spots that cannot be solved (no signal on the signature
  features) are `NA` by default; `TRUE` returns them as 0.

- BPPARAM:

  a
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  for the per-cluster plans and per-spot solves, or `NULL` (default) to
  run serially.

- chunk_size:

  integer. spots per parallel task.

## Value

a numeric matrix, cell types (rows) by spots (columns), spots in window
order.

## Examples

``` r
set.seed(1)
genes <- paste0("g", 1:30)
sig <- matrix(1, 30, 3, dimnames = list(genes, c("A", "B", "C")))
sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40
mix <- cbind(rbind(a <- runif(20, 0.3, 0.7), 1 - a, 0),
             matrix(c(0, 0, 1), 3, 20))
expr <- log2(matrix(rpois(30 * 40, sig %*% mix * 5), 30,
                    dimnames = list(genes, paste0("s", 1:40))) + 1)
clusters <- setNames(rep(c("mixed", "pure"), each = 20), colnames(expr))

# read the spots ten at a time
windows <- split(colnames(expr), rep(1:4, each = 10))
prop <- dwls_stream(
    read_window = function(i, features) expr[features, windows[[i]]],
    n_windows = length(windows), feat_ids = rownames(expr),
    clusters = clusters, sign_matrix = sig, n_cell = 5
)
all.equal(prop, dwls(expr, sig, clusters, n_cell = 5))
#> [1] TRUE
```
