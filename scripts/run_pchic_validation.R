#!/usr/bin/env Rscript

# Compute PCHiC enrichment metrics for datasets discovered under results/real/.
# With no --correction argument, both correction_F and correction_T are written.

args <- commandArgs(trailingOnly = TRUE)
correction_args <- grep("^--correction=", args, value = TRUE, ignore.case = TRUE)
only_method_args <- grep("^--only-method=", args, value = TRUE, ignore.case = TRUE)
append <- "--append" %in% args
datasets <- args[!grepl("^--", args)]
only_method <- if (length(only_method_args)) {
  tolower(sub("^--only-method=", "", only_method_args[[1L]], ignore.case = TRUE))
} else NULL
parse_correction <- function(value) {
  value <- toupper(sub("^--correction=", "", value, ignore.case = TRUE))
  if (!value %in% c("T", "TRUE", "F", "FALSE")) {
    stop("--correction must be T, TRUE, F, or FALSE", call. = FALSE)
  }
  value %in% c("T", "TRUE")
}
corrections <- if (length(correction_args)) parse_correction(correction_args[[1L]]) else c(FALSE, TRUE)

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})
source("R/validation.R")
source("scripts/PCHiC_enrichment.R")

default_references <- load_pchic_reference(project_root)
real_root <- file.path(project_root, "results", "real")
available <- list.dirs(real_root, recursive = FALSE, full.names = FALSE)
if (!length(datasets)) {
  datasets <- available
} else {
  datasets <- vapply(datasets, function(dataset) {
    hit <- available[tolower(available) == tolower(dataset)]
    if (length(hit) != 1L) stop("No results/real dataset named ", dataset, call. = FALSE)
    hit
  }, character(1L))
}
if (!length(datasets)) stop("No datasets found under ", real_root, call. = FALSE)

for (dataset in unique(datasets)) {
  references <- if (tolower(dataset) == "bmmc") {
    load_pchic_reference(project_root, celltype_map = bmmc_pchic_celltype_map)
  } else if (tolower(dataset) == "pbmc") {
    load_pchic_reference(project_root, celltype_map = pbmc_pchic_celltype_map)
  } else if (tolower(dataset) %in% c("dogma_stim", "dogma_ctrl")) {
    load_pchic_reference(project_root, celltype_map = dogma_pchic_celltype_map)
  } else if (tolower(dataset) == "share") {
    load_pchic_reference(project_root, celltype_map = share_pchic_celltype_map)
  } else if (tolower(dataset) == "brain") {
    load_pchic_reference(project_root, celltype_map = brain_pchic_celltype_map)
  } else {
    default_references
  }
  if (!length(intersect(names(references),
                        list.dirs(file.path(real_root, dataset), recursive = FALSE,
                                  full.names = FALSE)))) {
    message("Skipping ", dataset, ": no PCHiC-mapped celltype results")
    next
  }
  for (correction in corrections) {
    coordinate_chain <- if (tolower(dataset) == "share") {
      rtracklayer::import.chain(file.path(
        project_root, "reference", "PCHiC", "hg19ToHg38.over.chain"
      ))
    } else NULL
    metrics <- compute_dataset_pchic(
      project_root, dataset, references, allow_missing = TRUE,
      correction = correction, coordinate_chain = coordinate_chain,
      method_ids = if (is.null(only_method)) pchic_method_specs$method_id else only_method
    )
    output_dir <- file.path(project_root, "results", "validation", "PCHiC")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(output_dir,
                             paste0(dataset, "_enrichment_",
                                    if (correction) "T" else "F", ".tsv"))
    if (append && file.exists(output_file)) {
      previous <- data.table::fread(output_file, data.table = FALSE)
      if ("method" %in% names(previous)) {
        previous$method[previous$method == "Signac"] <- "Signac (corrected)"
        previous <- previous[previous$method != "Signac (uncorrected)", , drop = FALSE]
      }
      metrics <- rbind(previous, metrics)
    }
    data.table::fwrite(metrics, output_file, sep = "\t")
    message("Wrote ", nrow(metrics), " PCHiC rows (correction=", correction,
            ") to ", output_file)
  }
}
