#!/usr/bin/env Rscript

# Repair aggregate CRISPR enrichment tables when the per-celltype Signac rows
# are present but one aggregate Signac view is missing.  This reads the saved
# all_pairs files only; it does not rerun any association method.

suppressPackageStartupMessages(library(data.table))

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep(
  "^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."),
                             mustWork = TRUE)
setwd(project_root)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}

configuration <- yaml::read_yaml("config/datasets.yml")
real_root <- file.path(project_root, "results", "real")
validation_root <- file.path(project_root, "results", "validation", "CRISPR")

target_methods <- c(
  epilink_poisson = "EpiLink-Poisson",
  epilink_nb = "EpiLink-NB",
  poisson_glm = "Poisson-GLM",
  nb_glm = "NB-GLM",
  scmultimap = "scMultiMap",
  scmultimap_nobatch = "scMultiMap-noBatch",
  signac = "Signac",
  signac_uncorrected = "Signac (uncorrected)"
)

method_from_label <- function(label) {
  hit <- names(target_methods)[unname(target_methods) == label]
  if (length(hit)) hit[[1L]] else NULL
}

read_pairs <- function(path, full = FALSE) {
  data.table::fread(
    path,
    select = if (full) c("gene", "peak", "pvalue", "pvalue_uncorrected",
                          "pvalue_corrected") else c("gene", "peak"),
    data.table = FALSE, showProgress = FALSE
  )
}

common_universe <- function(dataset, celltype) {
  cached <- file.path(
    validation_root,
    paste0(dataset, "_", gsub("[^A-Za-z0-9]+", "_", celltype),
           "_tested_universe_F.rds")
  )
  if (!file.exists(cached)) {
    cached <- sub("_F[.]rds$", "_T.rds", cached)
  }
  if (file.exists(cached)) {
    return(readRDS(cached)[, c("gene", "peak"), drop = FALSE])
  }
  root <- file.path(real_root, dataset, celltype)
  method_ids <- names(target_methods)
  paths <- setNames(file.path(root, method_ids, "all_pairs.tsv.gz"), method_ids)
  paths <- paths[file.exists(paths)]
  if (!"signac" %in% names(paths) || length(paths) < 2L) {
    stop("Insufficient method results for ", dataset, "/", celltype,
         call. = FALSE)
  }
  tables <- lapply(paths, read_pairs, full = FALSE)
  keys <- lapply(tables, function(x) paste(x$gene, x$peak, sep = "\r"))
  common <- Reduce(intersect, keys)
  template <- tables[[1L]]
  template[paste(template$gene, template$peak, sep = "\r") %in% common,
           c("gene", "peak"), drop = FALSE]
}

selected_signac <- function(dataset, celltype, universe, method_id, cutoff) {
  path <- file.path(real_root, dataset, celltype, "signac", "all_pairs.tsv.gz")
  result <- read_pairs(path, full = TRUE)
  keys <- paste(universe$gene, universe$peak, sep = "\r")
  result_keys <- paste(result$gene, result$peak, sep = "\r")
  index <- match(keys, result_keys)
  pvalue <- validation_signac_pvalues(result, method_id)[index]
  padj <- stats::p.adjust(pvalue, method = "BH")
  universe[is.finite(padj) & padj < cutoff, c("gene", "peak"), drop = FALSE]
}

reference_for <- function(dataset) {
  settings_name <- names(configuration)[tolower(names(configuration)) ==
                                         tolower(dataset)]
  settings_name <- setdiff(settings_name, "defaults")
  settings <- utils::modifyList(configuration$defaults,
                                 configuration[[settings_name[[1L]]]])
  read_crispr_reference_units(settings$crispr_positive,
                              settings$crispr_negative)
}

aggregate_signac <- function(dataset, celltypes, method_id, cutoffs = c(.05, .10, .15, .20)) {
  reference <- reference_for(dataset)
  rows <- list()
  universes <- list()
  selections <- list()
  for (celltype in celltypes) {
    universe <- common_universe(dataset, celltype)
    universes[[celltype]] <- universe
    selections[[celltype]] <- list()
    for (cutoff in cutoffs) {
      selected <- selected_signac(dataset, celltype, universe, method_id, cutoff)
      selections[[celltype]][[format(cutoff, nsmall = 2L)]] <- selected
      row <- evaluate_crispr_reference_units(
        universe, selected, reference,
        method = if (method_id == "signac") "Signac (corrected)" else
          "Signac (uncorrected)",
        cutoff = cutoff, benchmark = "CRISPR-FlowFISH"
      )
      row$celltype <- celltype
      rows[[length(rows) + 1L]] <- row
    }
  }

  union_universe <- unique(do.call(rbind, universes))
  for (cutoff in cutoffs) {
    key <- format(cutoff, nsmall = 2L)
    union_selected <- unique(do.call(rbind, lapply(selections, `[[`, key)))
    union_selected <- union_selected[
      paste(union_selected$gene, union_selected$peak, sep = "\r") %in%
        paste(union_universe$gene, union_universe$peak, sep = "\r"),
      , drop = FALSE
    ]
    row <- evaluate_crispr_reference_units(
      union_universe, union_selected, reference,
      method = if (method_id == "signac") "Signac (corrected)" else
        "Signac (uncorrected)",
      cutoff = cutoff, benchmark = "CRISPR-FlowFISH"
    )
    row$celltype <- "All"
    rows[[length(rows) + 1L]] <- row
  }
  do.call(rbind, rows)
}

merge_rows <- function(path, rows) {
  previous <- data.table::fread(path, data.table = FALSE)
  previous$method[previous$method == "Signac"] <- "Signac (corrected)"
  previous <- previous[previous$method != unique(rows$method), , drop = FALSE]
  for (column in setdiff(names(previous), names(rows))) rows[[column]] <- NA
  rows <- rows[, names(previous), drop = FALSE]
  data.table::fwrite(
    data.table::rbindlist(list(previous, rows), fill = TRUE),
    path, sep = "\t", na = "NA"
  )
}

jobs <- list(
  list(dataset = "bmmc", celltypes = c("Bcell", "Dendritic", "Erythroid", "Myeloid", "T"),
       method_id = "signac_uncorrected",
       path = file.path(validation_root, "bmmc_FlowFISH_enrichment_T.tsv")),
  list(dataset = "Brain", celltypes = c("Astrocytes", "Excitatory neurons", "Oligodendroglia"),
       method_id = "signac_uncorrected",
       path = file.path(validation_root, "Brain_enrichment_F.tsv")),
  list(dataset = "Brain", celltypes = c("Astrocytes", "Excitatory neurons", "Oligodendroglia"),
       method_id = "signac_uncorrected",
       path = file.path(validation_root, "Brain_enrichment_T.tsv")),
  list(dataset = "SatpathyK562", celltypes = "K562", method_id = "signac",
       path = file.path(validation_root, "SatpathyK562_enrichment_F.tsv")),
  list(dataset = "SatpathyK562", celltypes = "K562", method_id = "signac",
       path = file.path(validation_root, "SatpathyK562_enrichment_T.tsv"))
)

for (job in jobs) {
  rows <- aggregate_signac(job$dataset, job$celltypes, job$method_id)
  merge_rows(job$path, rows)
  message("Updated ", job$path, " with ", nrow(rows), " ", unique(rows$method), " rows")
}
