#!/usr/bin/env Rscript

# Compute CRISPR enrichment metrics for datasets discovered under
# results/real/.  Each correction setting is written to the global
# results/validation/CRISPR/<validation type>/ directory as one flat table
# containing every matched cell type plus the union (celltype = "All").
# With no --correction argument, both correction_F and correction_T are written.

raw_args <- commandArgs(trailingOnly = TRUE)
correction_args <- grep("^--correction=", raw_args, value = TRUE,
                        ignore.case = TRUE)
reference_args <- grep("^--reference=", raw_args, value = TRUE,
                       ignore.case = TRUE)
only_method_args <- grep("^--only-method=", raw_args, value = TRUE,
                         ignore.case = TRUE)
append <- "--append" %in% raw_args
dataset_args <- raw_args[!grepl("^--", raw_args)]
only_method <- if (length(only_method_args)) {
  tolower(sub("^--only-method=", "", only_method_args[[1L]], ignore.case = TRUE))
} else NULL

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

parse_correction <- function(value) {
  value <- toupper(sub("^--correction=", "", value, ignore.case = TRUE))
  if (!value %in% c("T", "TRUE", "F", "FALSE")) {
    stop("--correction must be T, TRUE, F, or FALSE", call. = FALSE)
  }
  value %in% c("T", "TRUE")
}
corrections <- if (length(correction_args)) {
  parse_correction(correction_args[[1L]])
} else {
  c(FALSE, TRUE)
}

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}

real_root <- file.path(project_root, "results", "real")
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

configuration <- yaml::read_yaml("config/datasets.yml")
flowfish_reference <- file.path(
  project_root, "reference", "CRISPR", "41586_2021_3446_MOESM7_ESM.txt"
)
engreitz_reference <- file.path(
  project_root, "reference", "CRISPR_new",
  "EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"
)

settings_for_dataset <- function(dataset) {
  config_name <- names(configuration)[tolower(names(configuration)) == tolower(dataset)]
  config_name <- config_name[config_name != "defaults"]
  if (length(config_name) == 1L) {
    settings <- utils::modifyList(configuration$defaults, configuration[[config_name]])
    return(settings)
  }
  # Future datasets can be processed without adding a config entry when they
  # use the standard CRISPR reference files and already have results/real/<dataset>.
  message("No config entry for ", dataset,
          "; using the standard dataset defaults")
  utils::modifyList(configuration$defaults, list(
    benchmark_output_dir = file.path("results", "real", dataset)
  ))
}

method_specs <- data.frame(
  method = c("EpiLink-Poisson", "EpiLink-NB", "Poisson-GLM", "NB-GLM",
             "scMultiMap", "scMultiMap-noBatch", "Signac",
             "Signac (uncorrected)"),
  method_id = c("epilink_poisson", "epilink_nb", "poisson_glm", "nb_glm",
                "scmultimap", "scmultimap_nobatch", "signac",
                "signac_uncorrected"),
  stringsAsFactors = FALSE
)
methods_to_compute <- if (is.null(only_method)) {
  method_specs
} else {
  method_specs[tolower(method_specs$method_id) %in% only_method, , drop = FALSE]
}
if (!nrow(methods_to_compute)) {
  stop("Requested method is not available: ", paste(only_method, collapse = ", "),
       call. = FALSE)
}
cutoffs <- c(0.05, 0.10, 0.15, 0.20)

compute_dataset_metrics <- function(dataset, settings, correction) {
  dataset_root <- file.path(real_root, dataset)
  result_subdir <- if (tolower(dataset) == "brain2") "ALL" else NULL
  celltype_result_root <- function(celltype) {
    if (is.null(result_subdir)) {
      file.path(dataset_root, celltype)
    } else {
      file.path(dataset_root, celltype, result_subdir)
    }
  }
  celltypes <- list.dirs(dataset_root, recursive = FALSE, full.names = FALSE)
  celltypes <- celltypes[vapply(
    celltypes,
    function(celltype) file.exists(file.path(
      celltype_result_root(celltype), "epilink_poisson", "all_pairs.tsv.gz"
    )),
    logical(1L)
  )]
  if (tolower(dataset) == "brain") {
    expected_brain_celltypes <- c("Excitatory neurons", "Astrocytes", "Oligodendroglia")
    missing_brain <- setdiff(expected_brain_celltypes, celltypes)
    if (length(missing_brain)) {
      message("Brain CRISPR validation: missing celltype result directories: ",
              paste(missing_brain, collapse = ", "))
  }
  if (tolower(dataset) %in% c("dogma_stim", "dogma_ctrl")) {
    expected_dogma_celltypes <- c("Bcell", "Myeloid", "T")
    missing_dogma <- setdiff(expected_dogma_celltypes, celltypes)
    if (length(missing_dogma)) {
      message("Dogma CRISPR validation: missing celltype result directories: ",
              paste(missing_dogma, collapse = ", "))
    }
  }
  }
  if (!length(celltypes)) {
    stop("No completed celltype results found under ", dataset_root, call. = FALSE)
  }

  reference_by_celltype <- if (reference_mode == "engreitz") {
    list(K562 = read_epcrispr_reference_units(engreitz_reference,
                                              celltype = "K562"))
  } else {
    dataset_build <- if (is.null(settings$genome)) "hg38" else
      tolower(as.character(settings$genome))
    reference_chain <- if (dataset_build == "hg19") NULL else
      file.path(project_root, "reference", "CRISPR", "hg19ToHg38.over.chain")
    read_flowfish_reference_units(
      flowfish_reference,
      celltype_map = flowfish_celltype_map,
      coordinate_chain = reference_chain
    )
  }
  # Do not load or intersect result tables for cell types that have no
  # matching CRISPR reference.  Besides being biologically invalid, retaining
  # them unnecessarily makes multi-cell-type datasets very memory intensive.
  celltypes <- intersect(celltypes, names(reference_by_celltype))
  matched_celltypes <- intersect(celltypes, names(reference_by_celltype))
  if (!length(matched_celltypes)) {
    message("Skipping ", dataset,
            ": no celltype-matched CRISPR reference (", reference_mode, ")")
    return(NULL)
  }
  reference_units_all <- combine_crispr_reference_units(
    reference_by_celltype[matched_celltypes]
  )
  reference_label <- if (reference_mode == "engreitz") {
    "CRISPR-EPCrispr"
  } else {
    "CRISPR-FlowFISH"
  }
  # The reference is lifted above, so candidate and reference coordinates now
  # share the same build.  SHARE remains hg19 and is therefore also unlifted.
  coordinate_chain <- NULL

  available_for <- function(ct, method_id) {
    result_id <- validation_result_method_id(method_id)
    result_root <- celltype_result_root(ct)
    all_file <- file.path(result_root, result_id, "all_pairs.tsv.gz")
    object_file <- file.path(result_root, result_id, "epilink_object.rds")
    file.exists(all_file) &&
      (!grepl("^epilink_", method_id) || file.exists(object_file))
  }

  result_by_celltype <- list()
  object_by_celltype <- list()
  full_result_ids <- unique(vapply(
    methods_to_compute$method_id, validation_result_method_id, character(1L)
  ))
  for (ct in celltypes) {
    result_by_celltype[[ct]] <- list()
    object_by_celltype[[ct]] <- list()
    result_cache <- list()
    for (method_id in method_specs$method_id) {
      if (!available_for(ct, method_id)) next
      result_id <- validation_result_method_id(method_id)
      if (!is.null(result_cache[[result_id]])) {
        result <- result_cache[[result_id]]
      } else {
        result_path <- file.path(celltype_result_root(ct), result_id, "all_pairs.tsv.gz")
        # Keep the in-memory validation object small.  EpiLink selection is
        # taken from epilink_object.rds, while non-EpiLink methods only need
        # their adjusted p value and pair coordinates.
        header <- names(data.table::fread(result_path, nrows = 0L,
                                           showProgress = FALSE))
        select_cols <- intersect(c("gene", "peak", "padj", "pvalue",
                                   "pvalue_corrected", "pvalue_uncorrected"),
                                 header)
        result <- data.table::fread(
          result_path,
          select = select_cols,
          data.table = FALSE, showProgress = FALSE
        )
        result_cache[[result_id]] <- result
      }
      result_by_celltype[[ct]][[method_id]] <- result
      if (method_id %in% methods_to_compute$method_id && grepl("^epilink_", method_id)) {
        object_by_celltype[[ct]][[method_id]] <- readRDS(
          file.path(celltype_result_root(ct), method_id, "epilink_object.rds")
        )
      }
    }
  }

  # Restrict every method to exactly the same candidate gene-peak pairs.  This
  # prevents unequal tested universes from changing the contingency table.
  universe_by_celltype <- setNames(lapply(celltypes, function(ct) {
    pair_tables <- lapply(result_by_celltype[[ct]], function(x) {
      unique(x[, c("gene", "peak"), drop = FALSE])
    })
    key_sets <- lapply(pair_tables, function(x) paste(x$gene, x$peak, sep = "\r"))
    common_keys <- Reduce(base::intersect, key_sets)
    template <- pair_tables[[1L]]
    template[paste(template$gene, template$peak, sep = "\r") %in% common_keys,
             , drop = FALSE]
  }), celltypes)
  assessable_by_celltype <- setNames(lapply(matched_celltypes, function(ct) {
    reference_units_hit(
      universe_by_celltype[[ct]], reference_by_celltype[[ct]],
      coordinate_chain = coordinate_chain
    )
  }), matched_celltypes)

  restrict_to <- function(pairs, universe) {
    if (!nrow(pairs)) return(universe[FALSE, , drop = FALSE])
    keys <- paste(universe$gene, universe$peak, sep = "\r")
    keep <- paste(pairs$gene, pairs$peak, sep = "\r") %in% keys
    unique(pairs[keep, c("gene", "peak"), drop = FALSE])
  }

  selected_for_uncached <- function(ct, method_id, cutoff) {
    result <- result_by_celltype[[ct]][[method_id]]
    universe <- universe_by_celltype[[ct]]
    if (is.null(result)) return(NULL)
    if (grepl("^epilink_", method_id)) {
      selected <- EpiLinkFilter(
        object_by_celltype[[ct]][[method_id]], stat = "zval",
        fdrcutoff = cutoff, correction = correction,
        conservative = FALSE, verbose = FALSE
      )
      return(restrict_to(selected, universe))
    }
    if (method_id %in% validation_signac_method_ids) {
      signac <- restrict_to(result, universe)
      if (!nrow(signac)) return(universe[FALSE, , drop = FALSE])
      keys <- paste(signac$gene, signac$peak, sep = "\r")
      result_keys <- paste(result$gene, result$peak, sep = "\r")
      result_index <- match(keys, result_keys)
      signac$pvalue <- validation_signac_pvalues(result, method_id)[result_index]
      signac$padj_incoming <- stats::p.adjust(signac$pvalue, method = "BH")
      return(signac[is.finite(signac$padj_incoming) &
                      signac$padj_incoming < cutoff,
                    c("gene", "peak"), drop = FALSE])
    }
    restrict_to(result[is.finite(result$padj) & result$padj < cutoff,
                       c("gene", "peak"), drop = FALSE], universe)
  }

  selection_cache <- new.env(parent = emptyenv())
  selected_for <- function(ct, method_id, cutoff) {
    cache_key <- paste(ct, method_id, format(cutoff, nsmall = 2L), sep = "\r")
    if (!exists(cache_key, envir = selection_cache, inherits = FALSE)) {
      assign(cache_key, selected_for_uncached(ct, method_id, cutoff),
             envir = selection_cache)
    }
    get(cache_key, envir = selection_cache, inherits = FALSE)
  }

  metric_for_celltype <- function(ct, method, method_id, cutoff) {
    if (!available_for(ct, method_id)) return(NULL)
    universe <- universe_by_celltype[[ct]]
    selected <- selected_for(ct, method_id, cutoff)
    row <- evaluate_crispr_reference_units(
      universe, selected, reference_by_celltype[[ct]],
      method = method, cutoff = cutoff, benchmark = reference_label,
      assessable_units = assessable_by_celltype[[ct]],
      coordinate_chain = coordinate_chain
    )
    row$celltype <- ct
    row
  }

  metric_rows <- list()
  row_index <- 1L
  for (ct in matched_celltypes) {
    rows <- list()
    cell_row_index <- 1L
    for (index in seq_len(nrow(methods_to_compute))) {
      method <- methods_to_compute$method[index]
      method_id <- methods_to_compute$method_id[index]
      for (cutoff in cutoffs) {
        row <- metric_for_celltype(ct, method, method_id, cutoff)
        if (!is.null(row)) {
          rows[[cell_row_index]] <- row
          cell_row_index <- cell_row_index + 1L
        }
      }
    }
    if (length(rows)) {
      metric_rows[[row_index]] <- do.call(rbind, rows)
      row_index <- row_index + 1L
    }
  }

  for (index in seq_len(nrow(methods_to_compute))) {
    method <- methods_to_compute$method[index]
    method_id <- methods_to_compute$method_id[index]
    for (cutoff in cutoffs) {
      usable <- matched_celltypes[vapply(matched_celltypes, function(ct) {
        available_for(ct, method_id)
      }, logical(1L))]
      if (!length(usable)) next
      tested_union <- unique(do.call(rbind, universe_by_celltype[usable]))
      selected_union <- unique(do.call(rbind, lapply(
        usable, selected_for, method_id = method_id, cutoff = cutoff
      )))
      selected_union <- restrict_to(selected_union, tested_union)
      row <- evaluate_crispr_reference_units(
        tested_union, selected_union, reference_units_all,
        method = method, cutoff = cutoff, benchmark = reference_label,
        coordinate_chain = coordinate_chain
      )
      row$celltype <- "All"
      metric_rows[[row_index]] <- row
      row_index <- row_index + 1L
    }
  }

  if (!length(metric_rows)) {
    stop("No CRISPR enrichment rows were generated for ", dataset, call. = FALSE)
  }
  do.call(rbind, metric_rows)
}

for (dataset in unique(datasets)) {
  settings <- settings_for_dataset(dataset)
  for (correction in corrections) {
    correction_label <- if (correction) "correction_T" else "correction_F"
    metrics <- compute_dataset_metrics(dataset, settings, correction)
    if (is.null(metrics)) next
    metrics$reference <- reference_mode
    validation_root <- file.path(
      project_root, "results", "validation", "CRISPR", validation_name
    )
    dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)
    suffix <- if (correction) "T" else "F"
    output_file <- file.path(
      validation_root, paste0(dataset, "_enrichment_", suffix, ".tsv")
    )
    if (append && file.exists(output_file)) {
      previous <- utils::read.delim(output_file, stringsAsFactors = FALSE,
                                    check.names = FALSE)
      if ("method" %in% names(previous)) {
        previous$method[previous$method == "Signac"] <- "Signac (corrected)"
        previous <- previous[previous$method != "Signac (uncorrected)", , drop = FALSE]
      }
      metrics <- data.table::rbindlist(list(previous, metrics), fill = TRUE)
    }
    # Always replace the table, matching the PCHiC and HiChIP runners.
    utils::write.table(metrics, output_file,
                       sep = "\t", quote = FALSE, row.names = FALSE)
    message("Wrote ", nrow(metrics), " CRISPR rows for ", dataset,
            " (correction=", correction, ") to ", output_file)
  }
}
