# Build a DWLS signature matrix

Mean expression per cell type over the signature features: the reference
[`dwls()`](https://drieslab.github.io/spatialDWLS/reference/dwls.md)
expects.

## Usage

``` r
dwls_signature(expr, sign_gene, cell_type)
```

## Arguments

- expr:

  single-cell expression, features (rows) by cells (columns). Pass the
  scale the spots will be solved on, i.e. normalised and *not* logged
  ([`dwls()`](https://drieslab.github.io/spatialDWLS/reference/dwls.md)
  reverse-logs the spot expression before solving).

- sign_gene:

  character. signature features, e.g. cell type markers. Features
  missing from `expr` are dropped.

- cell_type:

  cell type of each column of `expr`.

## Value

a numeric matrix, signature features (rows) by cell type (columns, in
order of first appearance in `cell_type`).

## Examples

``` r
set.seed(1)
x <- matrix(rpois(200, 3), 20, 10,
            dimnames = list(paste0("g", 1:20), paste0("c", 1:10)))
dwls_signature(x, sign_gene = paste0("g", 1:5),
               cell_type = rep(c("A", "B"), each = 5))
#>      A   B
#> g1 4.2 4.0
#> g2 2.6 2.8
#> g3 3.2 3.0
#> g4 2.6 3.6
#> g5 3.0 4.2
```
