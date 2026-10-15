# spatialDWLS: cell type deconvolution of spatial spots

Estimate what fraction of each spot is made up of each cell type, given
a mean-expression reference built from single-cell data.

- [`dwls()`](https://drieslab.github.io/spatialDWLS/reference/dwls.md)
  deconvolves an in-memory expression matrix.

- [`dwls_stream()`](https://drieslab.github.io/spatialDWLS/reference/dwls_stream.md)
  does the same for a backend that reads spots in windows, so the
  expression never has to fit in memory at once.

- [`dwls_signature()`](https://drieslab.github.io/spatialDWLS/reference/dwls_signature.md)
  builds the reference from a single-cell matrix.

## References

Dong R, Yuan GC (2021). SpatialDWLS: accurate deconvolution of spatial
transcriptomic data. *Genome Biology* 22, 145.
[doi:10.1186/s13059-021-02362-7](https://doi.org/10.1186/s13059-021-02362-7)

Tsoucas D, Dong R, Chen H, Zhu Q, Guo G, Yuan GC (2019). Accurate
estimation of cell-type composition from gene expression data. *Nature
Communications* 10, 2975.
[doi:10.1038/s41467-019-10802-z](https://doi.org/10.1038/s41467-019-10802-z)

## See also

Useful links:

- <https://drieslab.github.io/spatialDWLS/>

- <https://github.com/drieslab/spatialDWLS>

- Report bugs at <https://github.com/drieslab/spatialDWLS/issues>

## Author

**Maintainer**: Ruben Dries <rubendries@gmail.com>
([ORCID](https://orcid.org/0000-0001-7650-7754))

Authors:

- Ruben Dries <rubendries@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-7650-7754))

- Rui Dong (original spatialDWLS implementation)

- Guo-Cheng Yuan <guo-cheng.yuan@mssm.edu>

Other contributors:

- Jiaji George Chen <jiajic@bu.edu> \[contributor\]
