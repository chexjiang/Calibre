#!/usr/bin/env Rscript

# Write the two canonical CRISPR figures for each correction variant and
# dataset.  The files are kept beside that dataset's other validation output:
#   validation/figures/odds_ratio_<type>_correction_{F,T}.pdf
# The metrics table retains the raw exact Fisher results; no confidence-interval
# figures are produced.

raw_args <- commandArgs(trailingOnly = TRUE)
correction_args <- grep("^--correction=", raw_args, value = TRUE,
                        ignore.case = TRUE)
reference_args <- grep("^--reference=", raw_args, value = TRUE,
                       ignore.case = TRUE)
dataset_args <- raw_args[!grepl("^--", raw_args)]

reference_mode <- if (length(reference_args)) {
  tolower(sub("^--reference=", "", reference_args[[1L]], ignore.case = TRUE))
} else {
  "flowfish"
}
if (reference_mode == "epcrispr") reference_mode <- "engreitz"
if (!reference_mode %in% c("flowfish", "engreitz")) {
  stop("--reference must be flowfish or engreitz", call. = FALSE)
}
validation_name <- if (reference_mode == "engreitz") {
  "CRISPR_EPCrispr"
} else {
  "CRISPR_FlowFISH"
}
plot_prefix <- if (reference_mode == "engreitz") {
  "crispr_epcrispr"
} else {
  "crispr_flowfish"
}

parse_correction <- function(value) {
  value <- toupper(sub("^--correction=", "", value))
  if (!value %in% c("T", "TRUE", "F", "FALSE")) {
    stop("--correction must be T, TRUE, F, or FALSE", call. = FALSE)
  }
  value %in% c("T", "TRUE")
}
corrections <- if (length(correction_args)) parse_correction(correction_args[[1L]]) else c(FALSE, TRUE)
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
validation_root <- file.path(project_root, "results", "validation", "CRISPR",
                             validation_name)
figure_dir <- file.path(validation_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
available_datasets <- list.dirs(real_root, recursive = FALSE, full.names = FALSE)
available_datasets <- setdiff(available_datasets, c("data_explore", "SatpathyK562"))
if (!length(dataset_args)) {
  datasets <- available_datasets
} else {
  datasets <- vapply(dataset_args, function(dataset) {
    hit <- available_datasets[tolower(available_datasets) == tolower(dataset)]
    if (length(hit) != 1L) {
      stop("No results/real dataset named ", dataset, call. = FALSE)
    }
    hit
  }, character(1L))
}
if (!length(datasets)) stop("No datasets found under ", real_root, call. = FALSE)

for (dataset in unique(datasets)) {
  for (correction in corrections) {
    correction_label <- if (correction) "correction_T" else "correction_F"
    suffix <- if (correction) "T" else "F"
    metrics_path <- file.path(
      validation_root, paste0(dataset, "_enrichment_", suffix, ".tsv")
    )
    if (!file.exists(metrics_path)) {
      message("Skipping ", dataset, " correction ", if (correction) "T" else "F",
              ": missing CRISPR metrics")
      next
    }
    metrics <- read.delim(metrics_path, stringsAsFactors = FALSE,
                          check.names = FALSE)
    if (!"celltype" %in% names(metrics)) {
      stop("CRISPR metrics must contain a celltype column: ", metrics_path,
           call. = FALSE)
    }

    p_or <- plot_validation_odds_ratio(
      metrics, dataset, paste0("CRISPR (correction=", correction, ")")
    )
    benchmark_title <- paste0("CRISPR (correction=", correction, ")")
    ggplot2::ggsave(file.path(figure_dir,
                              paste0(dataset, "_odds_ratio_", plot_prefix, "_correction_",
                                     suffix, ".pdf")),
                    p_or, width = 12, height = 8)
    message("Wrote CRISPR correction=", correction, " figures to ", figure_dir)
  }
}
