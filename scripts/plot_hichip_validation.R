#!/usr/bin/env Rscript

# Write HiChIP odds-ratio, p-value, and true-positive-count PDFs for each
# correction setting and requested dataset.

raw_args <- commandArgs(trailingOnly = TRUE)
datasets <- raw_args[!grepl("^--", raw_args)]

file_arguments <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arguments)) {
  script_path <- normalizePath(sub("^--file=", "", file_arguments[[1L]]))
  project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
  if (!file.exists(file.path(project_root, "R", "validation_plots.R"))) {
    project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
  }
}
setwd(project_root)

suppressPackageStartupMessages(library(ggplot2))
source("R/validation_plots.R")

real_root <- file.path(project_root, "results", "real")
validation_root <- file.path(project_root, "results", "validation", "HiChIP")
figure_dir <- file.path(validation_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
available_datasets <- list.dirs(real_root, recursive = FALSE, full.names = FALSE)
available_datasets <- setdiff(available_datasets, "data_explore")
if (!length(datasets)) {
  resolved <- available_datasets
} else {
  resolved <- vapply(datasets, function(dataset) {
    hit <- available_datasets[tolower(available_datasets) == tolower(dataset)]
    if (length(hit) != 1L) stop("No results/real dataset named ", dataset, call. = FALSE)
    hit
  }, character(1L))
}
if (!length(resolved)) stop("No datasets found under ", real_root, call. = FALSE)

for (dataset in unique(resolved)) {
  for (correction in c(FALSE, TRUE)) {
    suffix <- if (correction) "T" else "F"
    metrics_path <- file.path(validation_root, paste0(dataset, "_enrichment_", suffix, ".tsv"))
    if (!file.exists(metrics_path)) {
      message("Skipping ", dataset, " correction ", if (correction) "T" else "F",
              ": missing HiChIP metrics")
      next
    }
    metrics <- read.delim(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
    observed <- unique(as.character(metrics$celltype))
    metrics <- append_validation_all(metrics,
                                     celltypes = setdiff(observed, "All"))
    metrics <- metrics[as.character(metrics$celltype) %in%
                         unique(c(observed, "All")), , drop = FALSE]
    p_or <- plot_validation_odds_ratio(metrics, dataset,
                                       paste0("H3K27ac HiChIP (correction=", correction, ")"))
    p_pvalue <- plot_validation_pvalue(
      metrics, dataset, paste0("H3K27ac HiChIP (correction=", correction, ")")
    )
    p_tp <- plot_validation_true_positive_count(
      metrics, dataset, paste0("H3K27ac HiChIP (correction=", correction, ")")
    )
    ggplot2::ggsave(file.path(figure_dir,
                              paste0(dataset, "_odds_ratio_hichip_", suffix, ".pdf")),
                    p_or, width = 12, height = 8)
    ggplot2::ggsave(file.path(figure_dir,
                              paste0(dataset, "_pvalue_hichip_", suffix, ".pdf")),
                    p_pvalue, width = 12, height = 8)
    ggplot2::ggsave(file.path(figure_dir,
                              paste0(dataset, "_true_positive_count_hichip_",
                                     suffix, ".pdf")),
                    p_tp, width = 12, height = 8)
  }
  message("Wrote HiChIP correction F/T PDFs to ", figure_dir)
}
