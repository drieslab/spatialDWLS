# Deconvolve spots into cell type proportions

Runs both rounds of spatialDWLS on an expression matrix. Each spot's
result is a composition: the proportions of a solved spot sum to one.

## Usage

``` r
dwls(
  expr,
  sign_matrix,
  clusters,
  n_cell = 50,
  cutoff = 2,
  logbase = 2,
  na_as_zero = FALSE,
  BPPARAM = NULL,
  chunk_size = 256L
)
```

## Arguments

- expr:

  log-scaled, normalised expression: features (rows) by spots (columns),
  as a base `matrix` or a `Matrix` (sparse is kept sparse). Spot IDs are
  taken from the column names.

- sign_matrix:

  mean expression per cell type: signature features (rows) by cell type
  (columns), on the reverse-logged scale, e.g. from
  [`dwls_signature()`](https://drieslab.github.io/spatialDWLS/reference/dwls_signature.md).
  Not a binary marker matrix.

- clusters:

  cluster of every spot. Either named by spot ID (matched against
  `colnames(expr)`; extra entries, such as the clusters of a whole
  dataset when `expr` is a subset, are ignored) or unnamed in the column
  order of `expr`. Spots are grouped by cluster to pick candidate cell
  types in round 1. Every spot needs a non-`NA` cluster.

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

a numeric matrix, cell types (rows) by spots (columns).

## See also

[`dwls_stream()`](https://drieslab.github.io/spatialDWLS/reference/dwls_stream.md)
for expression read in windows.

## Examples

``` r
set.seed(1)
genes <- paste0("g", 1:30)
sig <- matrix(1, 30, 3, dimnames = list(genes, c("A", "B", "C")))
sig[1:10, "A"] <- sig[11:20, "B"] <- sig[21:30, "C"] <- 40

# 20 spots mixing A and B, 20 spots of C alone
mix <- cbind(rbind(a <- runif(20, 0.3, 0.7), 1 - a, 0),
             matrix(c(0, 0, 1), 3, 20))
counts <- matrix(rpois(30 * 40, sig %*% mix * 5), 30,
                 dimnames = list(genes, paste0("s", 1:40)))
clusters <- rep(c("mixed", "pure"), each = 20)

prop <- dwls(log2(counts + 1), sig, clusters, n_cell = 5)
round(prop[, c(1:3, 21:23)], 2)
#>     s1   s2   s3 s21 s22 s23
#> A 0.42 0.44 0.52   0   0   0
#> B 0.58 0.56 0.48   0   0   0
#> C 0.00 0.00 0.00   1   1   1
```
