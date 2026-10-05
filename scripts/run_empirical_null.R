#!/usr/bin/env Rscript

stop(
  "The empirical-null pipeline is now method-separated. Run: ",
  "(1) Rscript scripts/prepare_empirical_null.R; ",
  "(2) Rscript scripts/run_empirical_null_method.R <method> --workers=10; ",
  "(3) Rscript scripts/summarize_empirical_null.R; ",
  "(4) Rscript scripts/plot_empirical_null.R",
  call. = FALSE
)
