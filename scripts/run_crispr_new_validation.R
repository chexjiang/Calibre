#!/usr/bin/env Rscript

raw_args <- commandArgs(trailingOnly = TRUE)
bootstrap_arg <- grep("^--n-bootstrap=", raw_args, value = TRUE)
seed_arg <- grep("^--seed=", raw_args, value = TRUE)
only_method_args <- grep("^--only-method=", raw_args, value = TRUE,
                         ignore.case = TRUE)
append <- "--append" %in% raw_args
dataset_args <- raw_args[!grepl("^--", raw_args)]
only_method <- if (length(only_method_args)) {
  tolower(sub("^--only-method=", "", only_method_args[[1L]], ignore.case = TRUE))
} else NULL
n_bootstrap <- if (length(bootstrap_arg)) {
  as.integer(sub("^--n-bootstrap=", "", bootstrap_arg[[1L]]))
} else 200L
seed <- if (length(seed_arg)) {
  as.integer(sub("^--seed=", "", seed_arg[[1L]]))
} else 20260921L
if (!is.finite(n_bootstrap) || n_bootstrap < 0L) {
  stop("--n-bootstrap must be a non-negative integer", call. = FALSE)
}

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})
source("R/validation.R")
source("R/crispr_new.R")

benchmark_file <- file.path(
  project_root, "reference", "CRISPR_new",
  "crispr_new_snp_gene_benchmark.tsv.gz"
)
if (!file.exists(benchmark_file)) {
  stop("Missing prepared CRISPR_new benchmark. Run: Rscript scripts/prepare_crispr_new_reference.R",
       call. = FALSE)
}
benchmark <- read_crispr_new_benchmark(benchmark_file)
method_specs <- crispr_new_method_specs()
methods_to_compute <- if (is.null(only_method)) {
  method_specs
} else {
  method_specs[tolower(method_specs$method_id) %in% only_method, , drop = FALSE]
}
if (!nrow(methods_to_compute)) {
  stop("Requested method is not available: ", paste(only_method, collapse = ", "),
       call. = FALSE)
}
distance_thresholds <- c(1000L, 10000L, 100000L)

real_root <- file.path(project_root, "results", "real")
available <- list.dirs(real_root, recursive = FALSE, full.names = FALSE)
available <- setdiff(available, c("data_explore", "SatpathyK562"))
if (!length(dataset_args)) {
  datasets <- available
} else {
  datasets <- vapply(dataset_args, function(dataset) {
    hit <- available[tolower(available) == tolower(dataset)]
    if (length(hit) != 1L) {
      stop("No results/real dataset named ", dataset, call. = FALSE)
    }
    hit
  }, character(1L))
}

for (dataset in unique(datasets)) {
  dataset_root <- file.path(real_root, dataset)
  celltypes <- list.dirs(dataset_root, recursive = FALSE, full.names = FALSE)
  celltypes <- celltypes[file.exists(file.path(
    dataset_root, celltypes, "epilink_poisson", "all_pairs.tsv.gz"
  ))]
  if (!length(celltypes)) {
    message("Skipping ", dataset, ": no completed EpiLink-Poisson cell types")
    next
  }

  validation_root <- file.path(
    project_root, "results", "validation", "CRISPR", "CRISPR_pgBoost"
  )
  dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)
  metric_rows <- list()
  curve_rows <- list()
  mapping_rows <- list()
  scored_for_all <- list()
  metric_index <- curve_index <- mapping_index <- 1L

  for (celltype in celltypes) {
    message("Mapping CRISPR_pgBoost links for ", dataset, "/", celltype, " ...")
    baseline_file <- file.path(
      dataset_root, celltype, "epilink_poisson", "all_pairs.tsv.gz"
    )
    baseline <- data.table::fread(
      baseline_file, select = c("gene", "peak"), data.table = FALSE,
      showProgress = FALSE
    )
    mapping <- map_crispr_new_to_pairs(baseline, benchmark)
    if (!nrow(mapping)) {
      message("Skipping ", dataset, "/", celltype,
              ": no CRISPR_pgBoost SNP-gene links overlap tested peaks")
      next
    }
    mapping_rows[[mapping_index]] <- data.frame(
      dataset = dataset, celltype = celltype,
      n_peak_snp_gene = nrow(mapping),
      n_snp_gene = nrow(unique(mapping[, c("snp", "gene")])),
      n_positive = length(unique(paste(
        mapping$snp[mapping$positive], mapping$gene[mapping$positive], sep = "\r"
      ))),
      stringsAsFactors = FALSE
    )
    mapping_index <- mapping_index + 1L

    for (method_index in seq_len(nrow(methods_to_compute))) {
      method <- methods_to_compute$method[[method_index]]
      method_id <- methods_to_compute$method_id[[method_index]]
      result_file <- file.path(dataset_root, celltype,
                               validation_result_method_id(method_id),
                               "all_pairs.tsv.gz")
      if (!file.exists(result_file)) next
      header <- names(data.table::fread(result_file, nrows = 0L,
                                        showProgress = FALSE))
      pvalue_column <- if (method_id == "signac_uncorrected" &&
                           "pvalue_uncorrected" %in% header) {
        "pvalue_uncorrected"
      } else if (method_id == "signac" && "pvalue_corrected" %in% header) {
        "pvalue_corrected"
      } else intersect(c("pval", "pvalue"), header)
      if (!length(pvalue_column)) {
        message("Skipping ", dataset, "/", celltype, "/", method_id,
                ": no p-value column")
        next
      }
      results <- data.table::fread(
        result_file, select = c("gene", "peak", pvalue_column[[1L]]),
        data.table = FALSE, showProgress = FALSE
      )
      scored <- score_crispr_new_links(mapping, results, pvalue_column[[1L]])
      if (!nrow(scored)) next
      scored$dataset <- dataset
      scored$celltype <- celltype
      scored$method <- method
      scored$method_id <- method_id
      scored_for_all[[paste(celltype, method_id, sep = "\r")]] <- scored

      for (distance_threshold in distance_thresholds) {
        evaluated <- crispr_new_rank_metrics(scored, distance_threshold)
        if (is.null(evaluated)) next
        ci <- bootstrap_crispr_new_metrics(
          scored, distance_threshold, n_bootstrap = n_bootstrap,
          seed = seed + metric_index
        )
        row <- evaluated$metrics
        for (name in names(ci)) row[[name]] <- unname(ci[[name]])
        row$dataset <- dataset
        row$celltype <- celltype
        row$method <- method
        row$method_id <- method_id
        metric_rows[[metric_index]] <- row
        metric_index <- metric_index + 1L

        curve <- evaluated$curve
        curve$dataset <- dataset
        curve$celltype <- celltype
        curve$method <- method
        curve$method_id <- method_id
        curve$distance_threshold <- distance_threshold
        curve_rows[[curve_index]] <- curve
        curve_index <- curve_index + 1L
      }
    }
  }

  # Match pgBoost's across-cell-type comparison by taking the best score for
  # each SNP-gene pair across all available cell types, separately by method.
  for (method_index in seq_len(nrow(methods_to_compute))) {
    method <- methods_to_compute$method[[method_index]]
    method_id <- methods_to_compute$method_id[[method_index]]
    keys <- names(scored_for_all)[endsWith(names(scored_for_all),
                                           paste0("\r", method_id))]
    if (!length(keys)) next
    pooled <- data.table::rbindlist(scored_for_all[keys], fill = TRUE)
    pooled <- pooled[, .(
      positive = any(positive),
      distance_to_tss = min(distance_to_tss),
      score = if (all(is.na(score))) NA_real_ else max(score, na.rm = TRUE)
    ), by = .(snp, gene)]
    pooled <- as.data.frame(pooled)
    for (distance_threshold in distance_thresholds) {
      evaluated <- crispr_new_rank_metrics(pooled, distance_threshold)
      if (is.null(evaluated)) next
      ci <- bootstrap_crispr_new_metrics(
        pooled, distance_threshold, n_bootstrap = n_bootstrap,
        seed = seed + metric_index
      )
      row <- evaluated$metrics
      for (name in names(ci)) row[[name]] <- unname(ci[[name]])
      row$dataset <- dataset
      row$celltype <- "All"
      row$method <- method
      row$method_id <- method_id
      metric_rows[[metric_index]] <- row
      metric_index <- metric_index + 1L

      curve <- evaluated$curve
      curve$dataset <- dataset
      curve$celltype <- "All"
      curve$method <- method
      curve$method_id <- method_id
      curve$distance_threshold <- distance_threshold
      curve_rows[[curve_index]] <- curve
      curve_index <- curve_index + 1L
    }
  }

  if (!length(metric_rows)) {
    message("Skipping ", dataset, ": no evaluable CRISPR_pgBoost rows")
    next
  }
  metrics <- data.table::rbindlist(metric_rows, fill = TRUE)
  curves <- data.table::rbindlist(curve_rows, fill = TRUE)
  mappings <- data.table::rbindlist(mapping_rows, fill = TRUE)
  metrics_path <- file.path(
    validation_root, paste0(dataset, "_enrichment_metrics.tsv")
  )
  curves_path <- file.path(
    validation_root, paste0(dataset, "_enrichment_recall_curves.tsv.gz")
  )
  if (append && file.exists(metrics_path)) {
    previous <- data.table::fread(
      metrics_path, data.table = FALSE
    )
    previous$method[previous$method == "Signac"] <- "Signac (corrected)"
    previous <- previous[previous$method != "Signac (uncorrected)", , drop = FALSE]
    metrics <- data.table::rbindlist(list(previous, metrics), fill = TRUE)
    if (file.exists(curves_path)) {
      previous_curves <- data.table::fread(curves_path, data.table = FALSE)
      previous_curves$method[previous_curves$method == "Signac"] <- "Signac (corrected)"
      previous_curves <- previous_curves[
        previous_curves$method != "Signac (uncorrected)", , drop = FALSE
      ]
      curves <- data.table::rbindlist(list(previous_curves, curves), fill = TRUE)
    }
  }
  data.table::fwrite(metrics, metrics_path, sep = "\t")
  data.table::fwrite(curves, curves_path, sep = "\t", compress = "gzip")
  data.table::fwrite(
    mappings, file.path(validation_root, paste0(dataset, "_mapping_summary.tsv")),
    sep = "\t"
  )
  message("Wrote CRISPR_pgBoost validation for ", dataset, " to ", validation_root)
}
