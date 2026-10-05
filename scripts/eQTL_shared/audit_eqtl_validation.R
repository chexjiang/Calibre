#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(yaml)})

config <- yaml::read_yaml("config/eqtl_validation.yml")
rows <- list()
for (dataset in names(config$benchmarks)) {
  for (celltype in names(config$benchmarks[[dataset]])) {
    for (method in names(config$methods)) {
      raw_dir <- file.path("results/real", dataset, celltype, method)
      eqtl1 <- file.path("results/validation/eQTL/SCENT", dataset, celltype,
                         method, "summary.tsv")
      eqtl2 <- file.path("results/validation/eQTL/pgBoost", dataset, celltype,
                         method, "ranked_curve_primary.tsv")
      rows[[length(rows) + 1L]] <- data.table(
        dataset, celltype, method,
        all_pairs_available = file.exists(file.path(raw_dir, "all_pairs.tsv.gz")),
        significant_pairs_available = file.exists(file.path(raw_dir, "significant_pairs.tsv.gz")),
        eqtl1_complete = file.exists(eqtl1),
        eqtl2_complete = file.exists(eqtl2)
      )
    }
  }
}
audit <- rbindlist(rows)
dir.create("results/validation/eQTL", recursive = TRUE, showWarnings = FALSE)
fwrite(audit, "results/validation/eQTL/audit.tsv", sep = "\t")
cat("Configured combinations:", nrow(audit), "\n")
cat("eQTL1 complete:", sum(audit$eqtl1_complete), "\n")
cat("eQTL2 complete:", sum(audit$eqtl2_complete), "\n")
cat("Missing because upstream all_pairs is absent:",
    sum(!audit$all_pairs_available), "\n")
