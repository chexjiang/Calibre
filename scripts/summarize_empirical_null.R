#!/usr/bin/env Rscript

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml is required", call. = FALSE)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)
settings <- yaml::read_yaml("config/empirical_null.yml")
output_dir <- settings$output_dir
setup_dir <- file.path(output_dir, "setup")
pairs <- readRDS(file.path(setup_dir, "candidate_pairs.rds"))
registry <- empirical_null_methods()
n_pairs <- nrow(pairs)
n_methods <- nrow(registry)
tested_count <- raw_count <- operational_count <- matrix(
  0L, nrow = n_pairs, ncol = n_methods,
  dimnames = list(as.character(pairs$pair_id), registry$method_id)
)
replicate_rows <- list()
row_index <- 0L

for (method_index in seq_len(n_methods)) {
  method_id <- registry$method_id[[method_index]]
  for (replicate in seq_len(settings$n_replicates)) {
    path <- empirical_null_result_path(
      output_dir, replicate, method_id, settings$replicates_per_folder
    )
    row_index <- row_index + 1L
    if (!file.exists(path)) {
      replicate_rows[[row_index]] <- data.frame(
        replicate = replicate, method_id = method_id,
        method = registry$method[[method_index]], status = "missing",
        n_tested = 0L, raw_rejections = NA_integer_,
        operational_rejections = NA_integer_, raw_type1 = NA_real_,
        operational_rejection_rate = NA_real_, raw_any_rejection = NA,
        operational_any_rejection = NA, runtime_seconds = NA_real_,
        error = NA_character_, stringsAsFactors = FALSE
      )
      next
    }
    result <- utils::read.delim(path, stringsAsFactors = FALSE)
    tested <- is.finite(result$pval)
    indices <- match(result$pair_id[tested], pairs$pair_id)
    tested_count[indices, method_index] <- tested_count[indices, method_index] + 1L
    raw_count[indices, method_index] <- raw_count[indices, method_index] +
      as.integer(as.logical(result$raw_rejected[tested]))
    operational_count[indices, method_index] <-
      operational_count[indices, method_index] +
      as.integer(as.logical(result$operational_rejected[tested]))
    raw_rejections <- sum(as.logical(result$raw_rejected[tested]))
    operational_rejections <- sum(
      as.logical(result$operational_rejected[tested])
    )
    errors <- unique(stats::na.omit(result$error))
    replicate_rows[[row_index]] <- data.frame(
      replicate = replicate, method_id = method_id,
      method = registry$method[[method_index]],
      status = if (any(tested)) "complete" else "error",
      n_tested = sum(tested), raw_rejections = raw_rejections,
      operational_rejections = operational_rejections,
      raw_type1 = if (any(tested)) raw_rejections / sum(tested) else NA_real_,
      operational_rejection_rate = if (any(tested)) {
        operational_rejections / sum(tested)
      } else NA_real_,
      raw_any_rejection = if (any(tested)) raw_rejections > 0L else NA,
      operational_any_rejection = if (any(tested)) {
        operational_rejections > 0L
      } else NA,
      runtime_seconds = if ("runtime_seconds" %in% names(result)) {
        unique(result$runtime_seconds)[[1]]
      } else NA_real_,
      error = if (length(errors)) paste(errors, collapse = "; ") else NA_character_,
      stringsAsFactors = FALSE
    )
  }
}
replicate_metrics <- do.call(rbind, replicate_rows)
utils::write.table(
  replicate_metrics, file.path(output_dir, "replicate_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pair_rows <- lapply(seq_len(n_methods), function(method_index) {
  denominator <- tested_count[, method_index]
  data.frame(
    pair_id = pairs$pair_id, gene = pairs$gene, peak = pairs$peak,
    method_id = registry$method_id[[method_index]],
    method = registry$method[[method_index]],
    n_replicates = denominator,
    raw_type1 = ifelse(
      denominator > 0L, raw_count[, method_index] / denominator, NA_real_
    ),
    operational_type1 = ifelse(
      denominator > 0L,
      operational_count[, method_index] / denominator, NA_real_
    ),
    stringsAsFactors = FALSE
  )
})
pair_type1 <- do.call(rbind, pair_rows)
utils::write.table(
  pair_type1, file.path(output_dir, "pair_type1_error.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_rows <- lapply(seq_len(n_methods), function(method_index) {
  method_id <- registry$method_id[[method_index]]
  x <- replicate_metrics[
    replicate_metrics$method_id == method_id &
      replicate_metrics$status == "complete", , drop = FALSE
  ]
  pair_x <- pair_type1[
    pair_type1$method_id == method_id, , drop = FALSE
  ]
  data.frame(
    method_id = method_id, method = registry$method[[method_index]],
    n_complete_replicates = nrow(x),
    n_missing_replicates = sum(
      replicate_metrics$method_id == method_id &
        replicate_metrics$status == "missing"
    ),
    n_error_replicates = sum(
      replicate_metrics$method_id == method_id &
        replicate_metrics$status == "error"
    ),
    mean_pairwise_raw_type1 = mean(pair_x$raw_type1, na.rm = TRUE),
    mean_pairwise_operational_type1 = mean(
      pair_x$operational_type1, na.rm = TRUE
    ),
    mean_raw_false_discoveries = mean(x$raw_rejections, na.rm = TRUE),
    mean_operational_false_discoveries = mean(
      x$operational_rejections, na.rm = TRUE
    ),
    complete_null_fdr_raw = mean(x$raw_any_rejection, na.rm = TRUE),
    complete_null_fdr_operational = mean(
      x$operational_any_rejection, na.rm = TRUE
    ),
    mean_runtime_seconds = mean(x$runtime_seconds, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
method_summary <- do.call(rbind, summary_rows)
utils::write.table(
  method_summary, file.path(output_dir, "method_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("Empirical-null summaries written to ", output_dir)
