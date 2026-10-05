#!/usr/bin/env Rscript

# Migrate existing native Signac outputs without rerunning LinkPeaks.
# The legacy `pvalue` column in these files is the two-sided value
# (min(2 * Signac's one-tail value, 1)); recover the one-tail value and
# recompute both BH-adjusted columns over all candidate pairs.

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

project_root <- normalizePath(".", mustWork = TRUE)
config <- yaml::read_yaml(file.path(project_root, "config", "datasets.yml"))

default_fdr <- if (!is.null(config$defaults$fdr)) as.numeric(config$defaults$fdr) else 0.10

dataset_fdr <- function(dataset) {
  value <- config[[dataset]]$fdr
  if (is.null(value)) default_fdr else as.numeric(value)
}

write_tsv_gz_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(tmp)) unlink(tmp)
  ok <- FALSE
  tryCatch({
    con <- gzfile(tmp, open = "wt")
    utils::write.table(
      x, con, sep = "\t", quote = FALSE, row.names = FALSE,
      na = "NA", col.names = TRUE
    )
    close(con)
    ok <- TRUE
  }, error = function(e) {
    if (exists("con") && isOpen(con)) close(con)
    if (file.exists(tmp)) unlink(tmp)
    stop(e)
  })
  if (!ok || !file.rename(tmp, path)) {
    if (file.exists(tmp)) unlink(tmp)
    stop("Could not replace ", path)
  }
}

write_tsv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(tmp)) unlink(tmp)
  utils::write.table(
    x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA", col.names = TRUE
  )
  if (!file.rename(tmp, path)) {
    if (file.exists(tmp)) unlink(tmp)
    stop("Could not replace ", path)
  }
}

all_pair_files <- list.files(
  file.path(project_root, "results", "real"),
  pattern = "^all_pairs\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE
)
all_pair_files <- all_pair_files[basename(dirname(all_pair_files)) == "signac"]
all_pair_files <- sort(all_pair_files)

if (!length(all_pair_files)) {
  stop("No results/real/*/*/signac/all_pairs.tsv.gz files found")
}

summary_rows <- vector("list", length(all_pair_files))

for (i in seq_along(all_pair_files)) {
  all_pairs_path <- all_pair_files[[i]]
  signac_dir <- dirname(all_pairs_path)
  dataset <- basename(dirname(dirname(dirname(all_pairs_path))))
  fdr <- dataset_fdr(dataset)

  pairs <- data.table::fread(all_pairs_path, data.table = FALSE)
  if (!"pvalue_corrected" %in% names(pairs)) {
    if (!"pvalue" %in% names(pairs)) {
      stop("Missing pvalue in ", all_pairs_path)
    }
    pairs$pvalue_corrected <- as.numeric(pairs$pvalue)
  }
  if (!"pvalue_uncorrected" %in% names(pairs)) {
    pairs$pvalue_uncorrected <- pmin(
      pmax(as.numeric(pairs$pvalue_corrected) / 2, 0), 1
    )
  }

  # Recompute both adjustments from all candidate pairs. This is important:
  # padj_uncorrected must not be obtained by dividing padj_corrected by two.
  pairs$padj_uncorrected <- stats::p.adjust(
    as.numeric(pairs$pvalue_uncorrected), method = "BH"
  )
  pairs$padj_corrected <- stats::p.adjust(
    as.numeric(pairs$pvalue_corrected), method = "BH"
  )

  legacy_columns <- c("pvalue", "padj")
  base_columns <- setdiff(names(pairs), c(
    legacy_columns, "pvalue_uncorrected", "pvalue_corrected",
    "padj_uncorrected", "padj_corrected"
  ))
  pairs <- pairs[, c(
    base_columns, "pvalue_uncorrected", "pvalue_corrected",
    "padj_uncorrected", "padj_corrected"
  ), drop = FALSE]

  corrected <- pairs[is.finite(pairs$padj_corrected) & pairs$padj_corrected <= fdr,
                     , drop = FALSE]
  uncorrected <- pairs[
    is.finite(pairs$padj_uncorrected) & pairs$padj_uncorrected <= fdr,
    , drop = FALSE
  ]

  write_tsv_gz_atomic(pairs, all_pairs_path)
  write_tsv_gz_atomic(
    corrected, file.path(signac_dir, "significant_pairs.tsv.gz")
  )
  write_tsv_gz_atomic(
    corrected, file.path(signac_dir, "significant_pairs_corrected.tsv.gz")
  )
  write_tsv_gz_atomic(
    uncorrected, file.path(signac_dir, "significant_pairs_uncorrected.tsv.gz")
  )

  metadata_path <- file.path(signac_dir, "run_metadata.tsv")
  if (file.exists(metadata_path)) {
    metadata <- data.table::fread(metadata_path, data.table = FALSE)
    if (nrow(metadata) == 1L) {
      metadata$n_significant <- nrow(corrected)
      metadata$n_significant_uncorrected <- nrow(uncorrected)
      write_tsv_atomic(metadata, metadata_path)
    }
  }

  summary_rows[[i]] <- data.frame(
    dataset = dataset,
    celltype = basename(dirname(signac_dir)),
    n_pairs = nrow(pairs),
    n_significant_corrected = nrow(corrected),
    n_significant_uncorrected = nrow(uncorrected),
    fdr = fdr,
    stringsAsFactors = FALSE
  )
  message(
    "[", i, "/", length(all_pair_files), "] ", dataset, "/",
    basename(dirname(signac_dir)), ": ", nrow(corrected), " corrected; ",
    nrow(uncorrected), " uncorrected"
  )
}

summary <- do.call(rbind, summary_rows)
summary_path <- file.path(project_root, "results", "real", "signac_migration_summary.tsv")
write_tsv_atomic(summary, summary_path)
message("Wrote migration summary: ", summary_path)
