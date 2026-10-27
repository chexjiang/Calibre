# Utilities for the Brain2/Microglia differential-association label-null pilot.

da_null_methods <- function() {
  data.frame(
    method_id = c(
      "calibre_poisson", "calibre_nb", "glm_poisson", "glm_nb",
      "scmultimap"
    ),
    method = c(
      "CalibreDiff-Poisson", "CalibreDiff-NB", "Poisson-GLM", "NB-GLM",
      "scMultiMap"
    ),
    stringsAsFactors = FALSE
  )
}

da_null_method_group <- function(method_id) {
  method_id <- calibre_output_method_id(method_id)
  if (method_id %in% c("calibre_poisson", "glm_poisson")) {
    return(c("calibre_poisson", "glm_poisson"))
  }
  if (method_id %in% c("calibre_nb", "glm_nb")) {
    return(c("calibre_nb", "glm_nb"))
  }
  method_id
}

da_null_folder <- function(replicate, replicates_per_folder = 100L) {
  sprintf("%02d", ((as.integer(replicate) - 1L) %/%
                     as.integer(replicates_per_folder)) + 1L)
}

da_null_result_path <- function(output_dir, replicate, method_id,
                                replicates_per_folder = 100L) {
  file.path(
    output_dir,
    da_null_folder(replicate, replicates_per_folder),
    sprintf("replicate_%04d_%s.tsv.gz", as.integer(replicate), method_id)
  )
}

da_null_pair_key <- function(x) paste(x$gene, x$peak, sep = "\r")

empty_da_null_result <- function(pairs, method_id, method, error = NA_character_) {
  data.frame(
    pair_id = pairs$pair_id,
    gene = pairs$gene,
    peak = pairs$peak,
    method_id = method_id,
    method = method,
    statistic = NA_real_,
    pval = NA_real_,
    padj = NA_real_,
    raw_rejected = FALSE,
    operational_rejected = FALSE,
    error = error,
    stringsAsFactors = FALSE
  )
}

align_da_null_result <- function(results, pairs, method_id, method, alpha,
                                 operational_keys = NULL) {
  out <- empty_da_null_result(pairs, method_id, method)
  if (is.null(results) || !nrow(results)) return(out)
  matched <- match(da_null_pair_key(pairs), da_null_pair_key(results))
  statistic_column <- intersect(
    c("zval", "test_stat", "zscore", "score", "covar"), names(results)
  )
  if (length(statistic_column)) {
    out$statistic <- as.numeric(results[[statistic_column[[1L]]]][matched])
  }
  out$pval <- as.numeric(results$pval[matched])
  if ("padj" %in% names(results)) {
    out$padj <- as.numeric(results$padj[matched])
  } else {
    out$padj <- stats::p.adjust(out$pval, method = "BH")
  }
  out$raw_rejected <- is.finite(out$pval) & out$pval <= alpha
  if (is.null(operational_keys)) {
    out$operational_rejected <- is.finite(out$padj) & out$padj <= alpha
  } else {
    out$operational_rejected <- da_null_pair_key(pairs) %in% operational_keys
  }
  out
}

fit_da_null_calibre_family <- function(
    object, pairs, family, settings, seed, n_cores = 1L
) {
  calibre_id <- paste0("calibre_", family)
  glm_id <- paste0("glm_", family)
  registry <- da_null_methods()
  calibre_label <- registry$method[match(calibre_id, registry$method_id)]
  glm_label <- registry$method[match(glm_id, registry$method_id)]

  fitted <- CalibreDiff(
    object = object,
    gene_peak_df = pairs[, c("gene", "peak"), drop = FALSE],
    celltype_col = settings$celltype_column,
    family = family,
    celltype = NULL,
    condition_col = settings$condition_column,
    condition = NULL,
    condition_diff = NULL,
    # This is an unrestricted cell-level label null. Passing donor_col would
    # make CalibreDiff enforce one condition per donor, which is incompatible
    # with the permuted labels. Donor effect and donor-robust mode are both off.
    donor_col = NULL,
    donor_effect = FALSE,
    donor_robust = FALSE,
    covariates = unlist(settings$covariates, use.names = FALSE),
    null_generation = "s",
    stat = "zval",
    bin = isTRUE(settings$binarize_atac),
    selection = settings$selection,
    returnObject = TRUE,
    fdrcutoff = settings$alpha,
    ncore = n_cores,
    theory = FALSE,
    conservative = isTRUE(settings$calibre_conservative),
    correction = isTRUE(settings$calibre_correction),
    verbose = FALSE,
    seed = seed
  )
  if (is.null(fitted)) {
    return(setNames(list(
      empty_da_null_result(pairs, calibre_id, calibre_label,
                           "CalibreDiff returned no fitted pairs"),
      empty_da_null_result(pairs, glm_id, glm_label,
                           "CalibreDiff returned no fitted pairs")
    ), c(calibre_id, glm_id)))
  }

  results <- fitted@misc$Calibre$res_df
  significant <- CalibreFilter(
    fitted, stat = "zval", fdrcutoff = settings$alpha,
    correction = isTRUE(settings$calibre_correction),
    conservative = isTRUE(settings$calibre_conservative),
    returnObject = FALSE, verbose = FALSE
  )
  significant_keys <- da_null_pair_key(significant)
  out <- list(
    calibre = align_da_null_result(
      results, pairs, calibre_id, calibre_label, settings$alpha,
      operational_keys = significant_keys
    ),
    glm = align_da_null_result(
      results, pairs, glm_id, glm_label, settings$alpha
    )
  )
  names(out) <- c(calibre_id, glm_id)
  out
}

make_da_null_object <- function(rna, atac, metadata, condition_labels,
                                celltype_column = "newCT",
                                condition_column = "condition") {
  metadata <- metadata[colnames(rna), , drop = FALSE]
  metadata[[condition_column]] <- factor(
    as.character(condition_labels[colnames(rna)]), levels = c("Ctrl", "AD")
  )
  object <- Seurat::CreateSeuratObject(
    counts = rna, meta.data = metadata, assay = "RNA"
  )
  # The assays are reduced to the fixed benchmark feature set. Restore the
  # original full-assay sequencing depths for covariates and scMultiMap.
  if ("nCount_RNA" %in% names(metadata)) {
    object@meta.data$nCount_RNA <- as.numeric(metadata$nCount_RNA)
  }
  object[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = atac, sep = InferPeakSeparator(rownames(atac))
  )
  if ("nCount_ATAC" %in% names(metadata)) {
    object@meta.data$nCount_ATAC <- as.numeric(metadata$nCount_ATAC)
  }
  object
}

scmultimap_da_difference <- function(object, pairs, condition_labels,
                                     irls = TRUE) {
  if (!requireNamespace("scMultiMap", quietly = TRUE)) {
    stop("scMultiMap is required", call. = FALSE)
  }
  metadata <- object@meta.data
  cell_names <- colnames(object)
  labels <- as.character(condition_labels[cell_names])
  control <- which(labels == "Ctrl")
  case <- which(labels == "AD")
  if (!length(control) || !length(case)) {
    stop("both Ctrl and AD must contain cells", call. = FALSE)
  }

  pairs <- pairs[, c("gene", "peak"), drop = FALSE]
  assays <- c("RNA", "ATAC")
  feature_indices <- lapply(seq_along(assays), function(i) {
    assay_counts <- object[[assays[[i]]]]$counts
    match(unique(pairs[[i]]), rownames(assay_counts))
  })
  names(feature_indices) <- assays
  depths <- list(
    RNA = as.numeric(metadata$nCount_RNA),
    ATAC = as.numeric(metadata$nCount_ATAC)
  )
  if (any(!is.finite(depths$RNA)) || any(!is.finite(depths$ATAC))) {
    stop("full-assay RNA/ATAC depths are required for scMultiMap")
  }

  run_group <- function(indices) {
    count_list <- lapply(seq_along(assays), function(i) {
      assay_counts <- object[[assays[[i]]]]$counts
      as.matrix(Matrix::t(assay_counts[
        feature_indices[[i]], indices, drop = FALSE
      ]))
    })
    names(count_list) <- assays
    depth_list <- lapply(depths, function(x) x[indices])
    irls_results <- lapply(seq_along(assays), function(i) {
      scMultiMap::scMultiMap_IRLS(
        count_list[[i]], depth_list[[i]], bsample = NULL,
        irls = irls, verbose = FALSE
      )
    })
    names(irls_results) <- assays
    suppressMessages(scMultiMap::scMultiMap_WLS(
      count_list, depth_list, pairs, irls_results, verbose = FALSE
    ))
  }

  control_result <- run_group(control)
  case_result <- run_group(case)
  pair_key <- da_null_pair_key(pairs)
  control_key <- da_null_pair_key(control_result)
  case_key <- da_null_pair_key(case_result)
  control_covar <- as.numeric(control_result$covar[match(pair_key, control_key)])
  case_covar <- as.numeric(case_result$covar[match(pair_key, case_key)])
  case_covar - control_covar
}
