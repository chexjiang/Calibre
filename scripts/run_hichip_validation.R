#!/usr/bin/env Rscript

# Compute H3K27ac HiChIP enrichment metrics for datasets discovered under
# results/real/. With no --correction argument, both correction_F and
# correction_T are written.

raw_args <- commandArgs(trailingOnly = TRUE)
correction_args <- grep("^--correction=", raw_args, value = TRUE,
                        ignore.case = TRUE)
only_method_args <- grep("^--only-method=", raw_args, value = TRUE,
                         ignore.case = TRUE)
append <- "--append" %in% raw_args
datasets <- raw_args[!grepl("^--", raw_args)]
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

real_root <- file.path(project_root, "results", "real")
available_datasets <- list.dirs(real_root, recursive = FALSE, full.names = FALSE)
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

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(readxl)
})
source("R/validation.R")
source("scripts/HiCHIP_enrichment.R")

reference_file <- file.path(project_root, "reference", "HiChIP",
                            "41588_2017_BFng3963_MOESM4_ESM.xlsx")
promoters <- hichip_make_promoters()
default_references <- hichip_make_reference(reference_file, promoters)
message("HiChIP reference links: ",
        paste(names(default_references), vapply(default_references, length, integer(1L)),
              sep = "=", collapse = "; "))

for (dataset in unique(resolved)) {
  references <- if (tolower(dataset) == "bmmc") {
    hichip_make_reference(reference_file, promoters,
                          sheet_map = bmmc_hichip_sheet_map)
  } else if (tolower(dataset) == "pbmc") {
    hichip_make_reference(reference_file, promoters,
                          sheet_map = pbmc_hichip_sheet_map)
  } else if (tolower(dataset) %in% c("dogma_stim", "dogma_ctrl")) {
    hichip_make_reference(reference_file, promoters,
                          sheet_map = dogma_hichip_sheet_map)
  } else if (tolower(dataset) == "share") {
    hichip_make_reference(reference_file, promoters,
                          sheet_map = share_hichip_sheet_map)
  } else if (tolower(dataset) == "brain") {
    hichip_make_reference(reference_file, promoters,
                          sheet_map = brain_hichip_sheet_map)
  } else {
    default_references
  }
  if (!length(intersect(names(references),
                        list.dirs(file.path(real_root, dataset), recursive = FALSE,
                                  full.names = FALSE)))) {
    message("Skipping ", dataset, ": no HiChIP-mapped celltype results")
    next
  }
  for (correction in corrections) {
    coordinate_chain <- if (tolower(dataset) == "share") {
      rtracklayer::import.chain(file.path(
        project_root, "reference", "HiChIP", "hg19ToHg38.over.chain"
      ))
    } else NULL
    metrics <- compute_dataset_hichip(
      project_root, dataset, references, correction = correction,
      allow_missing = TRUE, coordinate_chain = coordinate_chain,
      method_ids = if (is.null(only_method)) hichip_method_specs$method_id else only_method
    )
    output_dir <- file.path(project_root, "results", "validation", "HiChIP")
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
    message("Wrote ", nrow(metrics), " HiChIP rows for ", dataset,
            " (correction=", correction, ") to ", output_file)
  }
}
