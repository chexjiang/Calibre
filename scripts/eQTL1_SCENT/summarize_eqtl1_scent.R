#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})
source("R/eqtl_validation.R")

config <- read_eqtl_validation_config("config/eqtl_validation.yml")
out_root <- "results/validation/eQTL/SCENT"
count_files <- Sys.glob(file.path(out_root, "*", "*", "*",
                                  "per_gene_counts.tsv.gz"))

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, probability))
}

rows <- lapply(count_files, function(path) {
  parts <- strsplit(path, .Platform$file.sep, fixed = TRUE)[[1L]]
  method <- parts[[length(parts) - 1L]]
  celltype <- parts[[length(parts) - 2L]]
  dataset <- parts[[length(parts) - 3L]]
  counts <- fread(path)
  old_path <- file.path(dirname(path), "summary.tsv")
  old <- if (file.exists(old_path)) fread(old_path) else data.table()
  pooled <- odds_ratio_from_counts(counts)
  boot <- bootstrap_scent_odds(
    counts, config$reference$bootstrap_replicates,
    config$reference$bootstrap_seed
  )
  finite_boot <- boot[is.finite(boot)]
  result <- cbind(data.frame(
    dataset = dataset, celltype = celltype,
    reference_group = config$benchmarks[[dataset]][[celltype]], method = method,
    n_selected_pairs = if (nrow(old)) old$n_selected_pairs[[1L]] else NA_integer_,
    n_bootstrap = length(boot)
  ), pooled, data.frame(
    bootstrap_median_or = if (length(finite_boot)) median(finite_boot) else NA_real_,
    bootstrap_ci_low = safe_quantile(boot, 0.025),
    bootstrap_ci_high = safe_quantile(boot, 0.975),
    bootstrap_probability_or_gt_1 = if (length(finite_boot))
      mean(finite_boot > 1) else NA_real_,
    n_genes = nrow(counts)
  ))
  fwrite(result, old_path, sep = "\t")
  fwrite(data.table(bootstrap_replicate = seq_along(boot), odds_ratio = boot),
         file.path(dirname(path), "gene_bootstrap_odds_ratios.tsv.gz"), sep = "\t")
  result
})

summary <- rbindlist(rows, fill = TRUE)
summary[, padj := p.adjust(pval, method = "BH")]
fwrite(summary, file.path(out_root, "all_results.tsv"), sep = "\t")
message("Recomputed ", nrow(summary), " SCENT-style summaries from cached counts")
