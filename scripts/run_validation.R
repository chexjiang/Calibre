#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
partial <- "--partial" %in% args
append <- "--append" %in% args
correction_args <- grep("^--correction=", args, value = TRUE)
only_method_args <- grep("^--only-method=", args, value = TRUE)
args <- args[!args %in% c("--partial", "--append") &
             !grepl("^--correction=", args) &
             !grepl("^--only-method=", args)]
if (length(args) != 2L) {
  stop("Usage: Rscript scripts/run_validation.R <dataset> <celltype> ",
       "[--correction=T|F] [--partial] [--only-method=METHOD_ID] [--append]",
       call. = FALSE)
}
dataset_name <- args[[1L]]
celltype <- args[[2L]]
only_method <- if (length(only_method_args)) {
  sub("^--only-method=", "", only_method_args[[1L]])
} else NULL
correction <- FALSE
if (length(correction_args)) {
  value <- toupper(sub("^--correction=", "", correction_args[[1L]]))
  if (!value %in% c("T", "TRUE", "F", "FALSE")) {
    stop("--correction must be T, TRUE, F, or FALSE", call. = FALSE)
  }
  correction <- value %in% c("T", "TRUE")
}
correction_label <- if (correction) "correction_T" else "correction_F"

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)
configuration <- yaml::read_yaml("config/datasets.yml")
settings <- utils::modifyList(configuration$defaults, configuration[[dataset_name]])
if (tolower(dataset_name) == "share") {
  # SHARE candidate peaks are hg19.  Validation is performed against the
  # pooled hg38 CRISPR reference after lifting candidate ranges on the fly.
  settings$crispr_reference_build <- "hg38"
  settings$crispr_positive <- "reference/CRISPR/crispr_hg38_positive.tsv"
  settings$crispr_negative <- "reference/CRISPR/crispr_hg38_negative.tsv"
}
coordinate_chain <- if (tolower(dataset_name) == "share") {
  rtracklayer::import.chain("reference/CRISPR/hg19ToHg38.over.chain")
} else NULL
result_root <- file.path(settings$benchmark_output_dir, celltype)
if (tolower(dataset_name) == "brain2") {
  result_root <- file.path(result_root, "ALL")
}
safe_celltype <- gsub("[^A-Za-z0-9]+", "_", celltype)
output_dir <- file.path(project_root, "results", "validation", "CRISPR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# lock_dir <- file.path(output_dir, ".validation.lock")
# if (!dir.create(lock_dir, showWarnings = FALSE)) {
#   message("Validation update already running; skipping this update")
#   quit(save = "no", status = 0)
# }
# on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)

method_specs <- data.frame(
  method = c("EpiLink-Poisson", "EpiLink-NB", "Poisson-GLM", "NB-GLM",
             "scMultiMap", "scMultiMap-noBatch", "Signac",
             "Signac (uncorrected)"),
  method_id = c("epilink_poisson", "epilink_nb", "poisson_glm", "nb_glm",
                "scmultimap", "scmultimap_nobatch", "signac",
                "signac_uncorrected"),
  stringsAsFactors = FALSE
)
has_result <- function(method_id) {
  result_id <- validation_result_method_id(method_id)
  all_file <- file.path(result_root, result_id, "all_pairs.tsv.gz")
  object_file <- file.path(result_root, result_id, "epilink_object.rds")
  file.exists(all_file) && (!grepl("^epilink_", method_id) || file.exists(object_file))
}
available <- vapply(method_specs$method_id, has_result, logical(1L))
if (partial) {
  method_specs <- method_specs[available, , drop = FALSE]
  if (!nrow(method_specs)) {
    message("No completed method result is available yet")
    quit(save = "no", status = 0)
  }
} else if (any(!available)) {
  stop("missing method result: ", paste(method_specs$method_id[!available], collapse = ", "), call. = FALSE)
}
methods_to_compute <- if (is.null(only_method)) {
  method_specs
} else {
  method_specs[method_specs$method_id %in% only_method, , drop = FALSE]
}
if (!nrow(methods_to_compute)) {
  stop("Requested method is not available: ", only_method, call. = FALSE)
}

cutoffs <- c(0.05, 0.10, 0.15, 0.20)
full_result_ids <- unique(vapply(
  methods_to_compute$method_id, validation_result_method_id, character(1L)
))
result_cache <- list()
method_results <- lapply(method_specs$method_id, function(method_id) {
  result_id <- validation_result_method_id(method_id)
  if (!is.null(result_cache[[result_id]])) return(result_cache[[result_id]])
  result_path <- file.path(result_root, result_id, "all_pairs.tsv.gz")
  result <- data.table::fread(
    result_path,
    select = if (result_id %in% full_result_ids) NULL else c("gene", "peak"),
    data.table = FALSE, showProgress = FALSE
  )
  result_cache[[result_id]] <<- result
  result
})
names(method_results) <- method_specs$method_id

# Every method is evaluated on the intersection of its tested candidate pairs.
# This is especially important in partial runs, where available result files
# may not contain exactly the same rows.
pair_tables <- lapply(method_results, function(x) {
  unique(x[, c("gene", "peak"), drop = FALSE])
})
key_sets <- lapply(pair_tables, function(x) paste(x$gene, x$peak, sep = "\r"))
common_keys <- Reduce(base::intersect, key_sets)
tested_universe <- pair_tables[[1L]]
tested_universe <- tested_universe[
  paste(tested_universe$gene, tested_universe$peak, sep = "\r") %in% common_keys,
  , drop = FALSE
]
saveRDS(
  tested_universe,
  file.path(output_dir, paste0(
    dataset_name, "_", safe_celltype, "_tested_universe_",
    if (correction) "T" else "F", ".rds"
  ))
)

reference_units <- read_crispr_reference_units(
  settings$crispr_positive, settings$crispr_negative
)
epilink_ids <- methods_to_compute$method_id[
  grepl("^epilink_", methods_to_compute$method_id)
]
epilink_objects <- lapply(epilink_ids, function(method_id) readRDS(file.path(result_root, method_id, "epilink_object.rds")))
names(epilink_objects) <- epilink_ids

restrict_to_universe <- function(pairs) {
  if (!nrow(pairs)) return(tested_universe[FALSE, , drop = FALSE])
  universe_key <- paste(tested_universe$gene, tested_universe$peak, sep = "\r")
  keep <- paste(pairs$gene, pairs$peak, sep = "\r") %in% universe_key
  unique(pairs[keep, c("gene", "peak"), drop = FALSE])
}

metric_rows <- list()
row_index <- 1L
for (index in seq_len(nrow(methods_to_compute))) {
  method <- methods_to_compute$method[index]
  method_id <- methods_to_compute$method_id[index]
  result <- method_results[[method_id]]
  for (cutoff in cutoffs) {
    if (grepl("^epilink_", method_id)) {
      # The correction flag is explicit so both EpiLink variants can be
      # reproduced without overwriting one another.
      selected <- EpiLinkFilter(epilink_objects[[method_id]], stat = "zval",
                                fdrcutoff = cutoff, correction = correction,
                                conservative = FALSE, verbose = FALSE)
      selected <- restrict_to_universe(selected)
    } else if (method_id %in% validation_signac_method_ids) {
      # Incoming Signac path restricts to the EpiLink-Poisson universe and
      # recalculates BH before applying strict '<' cutoff.  The uncorrected
      # view divides the stored corrected Signac p-values by two first.
      signac <- restrict_to_universe(result)
      if (nrow(signac)) {
        keys <- paste(signac$gene, signac$peak, sep = "\r")
        result_keys <- paste(result$gene, result$peak, sep = "\r")
        result_index <- match(keys, result_keys)
        signac$pvalue <- validation_signac_pvalues(result, method_id)[result_index]
        signac$padj_incoming <- stats::p.adjust(signac$pvalue, method = "BH")
        selected <- signac[is.finite(signac$padj_incoming) & signac$padj_incoming < cutoff,
                           c("gene", "peak"), drop = FALSE]
      } else selected <- tested_universe[FALSE, , drop = FALSE]
    } else {
      selected <- result[is.finite(result$padj) & result$padj < cutoff,
                         c("gene", "peak"), drop = FALSE]
      selected <- restrict_to_universe(selected)
    }
    metric_rows[[row_index]] <- evaluate_crispr_reference_units(
      tested_universe, selected, reference_units,
      method = method, cutoff = cutoff, benchmark = "CRISPR-FlowFISH",
      coordinate_chain = coordinate_chain
    )
    row_index <- row_index + 1L
  }
}
metrics <- do.call(rbind, metric_rows)
if (!is.null(only_method) && identical(tolower(only_method), "signac")) {
  metrics$method[metrics$method == "Signac"] <- "Signac (corrected)"
}
if (append) {
  output_file <- file.path(
    output_dir, paste0(dataset_name, "_", safe_celltype,
                       "_enrichment_", if (correction) "T" else "F", ".tsv")
  )
  if (file.exists(output_file)) {
    previous <- utils::read.delim(output_file, stringsAsFactors = FALSE,
                                  check.names = FALSE)
    if ("method" %in% names(previous)) {
      previous$method[previous$method == "Signac"] <- "Signac (corrected)"
      previous <- previous[previous$method != "Signac (uncorrected)", , drop = FALSE]
    }
    metrics <- rbind(previous, metrics)
  }
}
output_file <- file.path(
  output_dir, paste0(dataset_name, "_", safe_celltype,
                     "_enrichment_", if (correction) "T" else "F", ".tsv")
)
utils::write.table(metrics, output_file,
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("Enrichment metrics (correction=", correction, ") written to ",
        output_file)
