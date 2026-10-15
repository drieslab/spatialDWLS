# spatialDWLS

<!-- badges: start -->
[![R-CMD-check](https://github.com/drieslab/spatialDWLS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/drieslab/spatialDWLS/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/drieslab/spatialDWLS/actions/workflows/pkgdown.yaml/badge.svg)](https://drieslab.github.io/spatialDWLS/)
<!-- badges: end -->

Cell type deconvolution of spatial transcriptomics spots with dampened
weighted least squares.

`spatialDWLS` estimates what fraction of each spot is made up of each cell
type, given a mean-expression reference built from single-cell data. Each
spot's result is a composition: the proportions of a solved spot sum to one.

The method (Dong and Yuan 2021) runs in two rounds:

1. **Cluster round.** For each cluster of spots, cell types enriched in it
   (PAGE z-score above `cutoff`) become candidates, and every spot is solved
   against them.
2. **Spot round.** Per spot, cell types whose round-1 proportion is below
   `1 / n_cell` are dropped and the spot is solved again against the rest.

Each solve is DWLS (Tsoucas et al. 2019): non-negative weighted least squares
whose weights are iterated and capped by a cross-validated dampening constant.

This package is the standalone engine. It works on plain matrices, so it has
no dependency on Giotto; [Giotto](https://github.com/giotto-suite/Giotto)'s
`runDWLSDeconv()` calls it.

**Website:** <https://drieslab.github.io/spatialDWLS/>, with two worked
examples: [plain matrices](https://drieslab.github.io/spatialDWLS/articles/spatialDWLS.html)
and [through Giotto](https://drieslab.github.io/spatialDWLS/articles/giotto.html).

## Installation

```r
# install.packages("remotes")
remotes::install_github("drieslab/spatialDWLS")
```

## Usage

```r
library(spatialDWLS)

# reference: mean normalised (not logged) expression per cell type
sig <- dwls_signature(sc_norm, sign_gene = markers, cell_type = sc_labels)

# spots: log-normalised expression, features x spots; one cluster per spot
prop <- dwls(spot_lognorm, sig, clusters = spot_clusters, n_cell = 10)
```

`prop` is a cell types x spots matrix. Spots with no signal on the signature
features cannot be solved and are `NA` (`na_as_zero = TRUE` returns 0).

- **Sparse input** stays sparse; only the signature features are densified,
  one cluster at a time.
- **Parallel:** pass a `BiocParallel` param, e.g.
  `BPPARAM = BiocParallel::MulticoreParam(8)`. Results are identical to the
  serial run.
- **Larger than memory:** `dwls_stream()` takes a `read_window(i, features)`
  callback and reads spots one window at a time. This is how Giotto's
  disk-backed stores (GiottoDisk) run it.
- **Reproducible:** the dampening constant's cross-validation is seeded
  internally; your session's RNG state is left untouched.

## Citation

Dong R, Yuan GC (2021). SpatialDWLS: accurate deconvolution of spatial
transcriptomic data. *Genome Biology* 22, 145.
<https://doi.org/10.1186/s13059-021-02362-7>

Tsoucas D, Dong R, Chen H, Zhu Q, Guo G, Yuan GC (2019). Accurate estimation
of cell-type composition from gene expression data. *Nature Communications*
10, 2975. <https://doi.org/10.1038/s41467-019-10802-z>
