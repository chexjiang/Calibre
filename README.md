# EpiLink

EpiLink builds cell-type-specific enhancer-gene maps from paired single-cell
RNA and ATAC measurements. This directory is the clean, reproducible rewrite;
`incoming/` is an immutable historical snapshot and is never used as an output
location.

## Layout

- `R/`: tested method implementation.
- `config/datasets.yml`: dataset-specific inputs and analysis settings.
- `scripts/run_epilink.R`: command-line real-data runner.
- `simulation/` and `scripts/run_simulation.R`: configured Poisson/NB/ZINB study.
- `tests/testthat/`: unit and small end-to-end tests.
- `incoming/`: original code, data, and results (ignored by Git).
- `results/`: newly generated results (ignored by Git).

## Test

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Run

Run from the project root:

```bash
Rscript scripts/run_epilink.R pbmc Bcell poisson
Rscript scripts/run_epilink.R pbmc Bcell nb
Rscript scripts/run_epilink.R bmmc Myeloid poisson
Rscript scripts/run_epilink.R bmmc Myeloid nb
```

Run the model-based simulation profile and generate PDF plots:

```bash
Rscript scripts/run_simulation.R model_based
Rscript scripts/plot_model_based_simulation.R
```

The empirical-null benchmark is run with `scripts/run_empirical_null.R` and
its PDF is generated with `scripts/plot_empirical_null.R`.

Formal EpiLink runs use the standardized Pearson residual bootstrap together
with the multiplicative finite-cell correction and the plus-one empirical-FDP
numerator. Use `diagnostic` for high-replicate global-null calibration and
`power_pilot` for the small-effect sample-size pilot. See
`docs/simulation_calibration.md` for the current findings and limitations.
The `trend` profile adds the Signac benchmark and ten enhancer-effect sizes for
the combined FDR/power figure.

The runner validates all configured inputs before loading large matrices. New
results are written below `results/`; existing files are not overwritten unless
`--overwrite` is supplied.
