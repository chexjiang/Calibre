# Calibre

Calibre is an R package for inferring cell-type-specific enhancer–gene links
from paired single-cell RNA and ATAC measurements. It fits count-regression
models, generates residual-based synthetic null statistics, and calibrates a
shared empirical false-discovery-rate threshold.

## Repository contents

- `R/`: Calibre implementation and supporting methods.
- `man/`: package documentation.
- `config/`: YAML configuration files for simulation and validation analyses.
- `docs/`: methodological notes and analysis documentation.
- `DESCRIPTION`, `NAMESPACE`, and `LICENSE`: R package metadata.

## Installation

From a local clone, install the package with:

```bash
R CMD INSTALL .
```

The required R dependencies are listed in `DESCRIPTION`. Optional integrations
are listed under `Suggests`.

## Core API

The main functions are:

- `Calibre()`: infer enhancer–gene links for a Seurat multiome object.
- `CalibreFilter()`: apply the calibrated FDR threshold to a Calibre result.
- `CalibreDiff()`: compare enhancer–gene associations between conditions.
- `CalibreDynamic()`: model associations that vary along pseudotime.
- `Calibre_single()`: fit one enhancer–gene model.
- `calibre_simulation_methods()`: list the configured simulation methods.
- `simulate_calibre_data()`: generate paired RNA/ATAC simulation data.
- `run_calibre_simulation()`: evaluate selected methods on simulated data.

Example:

```r
library(Calibre)

result <- Calibre(
  object = multiome_object,
  gene_peak_df = candidate_pairs,
  celltype_col = "celltype",
  family = "poisson",
  returnObject = TRUE
)

significant_links <- CalibreFilter(
  result,
  stat = "zval",
  fdrcutoff = 0.05
)
```

See the files under `man/` and the methodological notes under `docs/` for
function arguments, return values, and analysis details.

## License

Calibre is released under the MIT License. See `LICENSE`.
