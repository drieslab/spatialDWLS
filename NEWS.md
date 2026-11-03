# spatialDWLS 0.1.2

## new features
* `dwls_design()` returns the dampened WLS problem of each spot as data: the
  reverse-logged expression on the signature features, the signature on those
  features and the dampening constant `J` of each spot (cross-validated per
  cluster, as in the first round of `dwls()`).
* `dwls_expression_loss()` scores cell type proportions from any source
  against the spot expression with spatialDWLS's dampened weights, and returns
  the gradient with the weights held fixed. With it a model that types single
  cells can be trained from spot expression.

# spatialDWLS 0.1.1

## bug fixes
* A cluster label with no spots in the expression no longer stops the run
  with "quadratic program could not be solved". Named `clusters` may now be a
  superset, e.g. the clusters of a whole dataset when `expr` is a subset.
* In a cluster with a single candidate cell type, a spot with no signal on
  the signature features is now `NA`, as it is in every other cluster; it
  was returned as proportion 1 of that type.
* `dwls_stream()` matches the second round's spots to the first round's by
  name, so a `read_window()` that returns a window's spots in a different
  order on a later pass no longer gives silently wrong results.
* Duplicated or missing cell type names, feature names or spot IDs are now an
  error; a duplicated cell type silently overwrote another's proportions.
  `NA` clusters get an error that says so.
* Parallel runs no longer send every task's expression slice to each worker.

## changes
* Noted for 0.1.0: the PAGE step that picks candidate cell types reverse-logs
  with `logbase`. Giotto <= 4.3.1 used base 2 there whatever `logbase` was, so
  results differ from it slightly for `logbase` other than 2 (identical at the
  default).

# spatialDWLS 0.1.0

First release as a standalone package. The algorithm was previously
`runDWLSDeconv()` inside Giotto (<= 4.3.1), which now calls this package.

Results match Giotto 4.3.1 to floating point (~1e-15), except where these bugs
were fixed:

* Clusters with a single candidate cell type in the spot round now assign that
  type proportion 1; they were left at 0.
* Small negative values from the constrained solver are now clipped to 0 for
  every cell type, not only the last one.
* Running DWLS no longer leaves the session's random number generator in a
  fixed state.

Faster (12-22x serial on 600-20,000 spots):

* The dampening constant's cross-validation uses `.lm.fit()` on k x k
  weighted least squares instead of `lm()`, and draws its 100 subsets once.
* Weighted solves use `crossprod(S * w, S)` instead of a features x features
  `diag(w)`, and call `Matrix::nearPD()` only when the matrix needs it.
* The reverse-log transform keeps sparse input sparse and touches only the
  signature features.

New:

* `BPPARAM` for parallel per-cluster plans and per-spot solves.
* `dwls_stream()` for expression read in windows.
* `na_as_zero`.
* `dwls_signature()`, the matrix form of Giotto's
  `makeSignMatrixDWLSfromMatrix()`.
