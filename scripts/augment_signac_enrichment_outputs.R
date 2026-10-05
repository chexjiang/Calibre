#!/usr/bin/env Rscript

# Normalize legacy Signac labels in the centralized CRISPR validation tables.
# The validation layout is now results/validation/CRISPR/<dataset>_*.tsv;
# this helper no longer writes into results/real/<dataset>/validation/.

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(getwd(), mustWork = TRUE)
validation_root <- file.path(project_root, "results", "validation", "CRISPR")
if (!dir.exists(validation_root)) {
  stop("Missing centralized CRISPR validation directory: ", validation_root,
       call. = FALSE)
}

files <- list.files(validation_root, pattern = "_enrichment_[FT][.]tsv$",
                    full.names = TRUE)
report <- list()

for (path in sort(files)) {
  tab <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  changed <- FALSE
  if ("method" %in% names(tab)) {
    changed <- any(tab$method == "Signac")
    tab$method[tab$method == "Signac"] <- "Signac (corrected)"
  }
  if (changed) data.table::fwrite(tab, path, sep = "\t", na = "NA")
  report[[length(report) + 1L]] <- data.frame(
    file = path, normalized_signac_label = changed,
    stringsAsFactors = FALSE
  )
}

if (length(report)) {
  data.table::fwrite(
    data.table::rbindlist(report, fill = TRUE),
    file.path(validation_root, "signac_enrichment_augmentation_summary.tsv"),
    sep = "\t"
  )
}
