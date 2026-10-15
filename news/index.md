# Changelog

## spatialDWLS 0.1.0

First release as a standalone package. The algorithm was previously
`runDWLSDeconv()` inside Giotto (\<= 4.3.1), which now calls this
package.

Results match Giotto 4.3.1 to floating point (~1e-15), except where
these bugs were fixed:

- Clusters with a single candidate cell type in the spot round now
  assign that type proportion 1; they were left at 0.
- Small negative values from the constrained solver are now clipped to 0
  for every cell type, not only the last one.
- Running DWLS no longer leaves the session’s random number generator in
  a fixed state.

Faster (12-22x serial on 600-20,000 spots):

- The dampening constant’s cross-validation uses
  [`.lm.fit()`](https://rdrr.io/r/stats/lmfit.html) on k x k weighted
  least squares instead of [`lm()`](https://rdrr.io/r/stats/lm.html),
  and draws its 100 subsets once.
- Weighted solves use `crossprod(S * w, S)` instead of a features x
  features `diag(w)`, and call
  [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html) only
  when the matrix needs it.
- The reverse-log transform keeps sparse input sparse and touches only
  the signature features.

New:

- `BPPARAM` for parallel per-cluster plans and per-spot solves.
- [`dwls_stream()`](https://drieslab.github.io/spatialDWLS/reference/dwls_stream.md)
  for expression read in windows.
- `na_as_zero`.
- [`dwls_signature()`](https://drieslab.github.io/spatialDWLS/reference/dwls_signature.md),
  the matrix form of Giotto’s `makeSignMatrixDWLSfromMatrix()`.
