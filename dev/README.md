# Development scripts

Not part of the package build (`dev/` is in `.Rbuildignore`).

| Script | What it does |
|---|---|
| `sim.R` | `simulate_dwls()`: Visium-like spots with known cell type proportions, a matching signature, and a cluster column. Matrix-only, no Giotto. |
| `equivalence.R` | Runs the pre-extraction Giotto DWLS (taken from Giotto's git history) and `spatialDWLS::dwls()` on the same inputs, and reports max abs difference and timing. |
| `bench.R` | Timing (serial / parallel) and accuracy against the known truth, at a given number of spots. |

```sh
# from the package root
Rscript dev/equivalence.R ../GiottoSuite/Giotto   # path to a Giotto checkout
Rscript dev/bench.R 20000                  # spots
NC=5 Rscript dev/bench.R 2000 3            # n_cell = 5, 3 single-type niches
```
