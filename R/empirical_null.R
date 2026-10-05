# Utilities for the PBMC empirical-null type-I error benchmark.

#' Empirical-null benchmark method registry
#' @return A data frame with stable method identifiers and labels.
empirical_null_methods <- function() {
  data.frame(
    method_id = c(
      "epilink_poisson", "epilink_nb", "glm_poisson", "glm_nb",
      "scmultimap", "signac", "scent"
    ),
    method = c(
      "EpiLink-Poisson", "EpiLink-NB", "Poisson-GLM", "NB-GLM",
      "scMultiMap", "Signac", "SCENT"
    ),
    stringsAsFactors = FALSE
  )
}

#' Permute a complete ATAC matrix across cells
#'
#' One permutation is shared by every peak, preserving the empirical ATAC
#' matrix while breaking its cell-level alignment with the fixed RNA matrix.
#' @param atac Peak-by-cell count matrix.
#' @param seed Permutation seed.
#' @return A list containing the permuted counts and permutation audit table.
permute_empirical_atac <- function(atac, seed) {
  if (is.null(colnames(atac)) || anyDuplicated(colnames(atac))) {
    stop("ATAC cells must have unique names", call. = FALSE)
  }
  target_cells <- colnames(atac)
  permutation <- with_epilink_seed(as.integer(seed), sample.int(ncol(atac)))
  permuted <- atac[, permutation, drop = FALSE]
  source_cells <- colnames(permuted)
  colnames(permuted) <- target_cells
  list(
    counts = permuted,
    permutation = data.frame(
      target_cell = target_cells,
      source_cell = source_cells,
      source_index = permutation,
      stringsAsFactors = FALSE
    )
  )
}

empirical_pair_key <- function(x) paste(x$gene, x$peak, sep = "\r")

#' Choose the fixed candidate-pair universe before null permutations
#' @param object Unpermuted PBMC/Tnk Seurat object.
#' @param pairs Candidate gene-peak pairs.
#' @param n_pairs Number of pairs to sample.
#' @param selection Candidate filtering rule passed to FilterGenePeakPairs.
#' @param seed Sampling seed.
#' @return A fixed pair table with pair_id.
sample_empirical_pairs <- function(object, pairs, n_pairs, selection, seed) {
  filtered <- FilterGenePeakPairs(
    object, pairs, selection = selection, verbose = TRUE
  )
  filtered <- unique(filtered[, c("gene", "peak"), drop = FALSE])
  if (nrow(filtered) < n_pairs) {
    stop("Only ", nrow(filtered), " eligible pairs; requested ", n_pairs,
         call. = FALSE)
  }
  chosen <- with_epilink_seed(
    as.integer(seed), sample.int(nrow(filtered), as.integer(n_pairs))
  )
  out <- filtered[chosen, , drop = FALSE]
  rownames(out) <- NULL
  out$pair_id <- seq_len(nrow(out))
  out
}

empty_empirical_result <- function(pairs, method_id, method, error = NA_character_) {
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

align_empirical_result <- function(results, pairs, method_id, method, alpha,
                                   operational_keys = NULL) {
  out <- empty_empirical_result(pairs, method_id, method)
  if (is.null(results) || !nrow(results)) return(out)
  matched <- match(empirical_pair_key(pairs), empirical_pair_key(results))
  statistic_column <- intersect(
    c("zval", "test_stat", "zscore", "score", "cor", "covar"),
    names(results)
  )
  if (length(statistic_column)) {
    out$statistic <- as.numeric(results[[statistic_column[[1]]]][matched])
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
    out$operational_rejected <- empirical_pair_key(pairs) %in% operational_keys
  }
  out
}

#' Fit a BH-adjusted Poisson or NB GLM on one null replicate
fit_empirical_glm <- function(object, pairs, family, alpha, covariates,
                              n_cores = 1L, selection = "rate",
                              binarize_atac = TRUE, seed = 1L) {
  method_id <- paste0("glm_", family)
  registry <- empirical_null_methods()
  method <- registry$method[match(method_id, registry$method_id)]
  results <- EpiLink(
    object = object, gene_peak_df = pairs[, c("gene", "peak")],
    family = family, celltype = NULL, covariates = covariates,
    null_generation = "s", stat = "zval", bin = binarize_atac,
    selection = selection, returnObject = FALSE, fdrcutoff = alpha,
    ncore = n_cores, theory = TRUE, correction = FALSE,
    conservative = FALSE, verbose = FALSE, seed = seed
  )
  align_empirical_result(results, pairs, method_id, method, alpha)
}

#' Fit EpiLink and its corresponding BH-GLM on one null replicate
#'
#' A single EpiLink fit supplies the observed GLM p-values used by both rows,
#' avoiding a redundant second model fit.
fit_empirical_epilink_family <- function(object, pairs, family, alpha,
                                         covariates, n_cores, seed,
                                         selection = "rate", binarize_atac = TRUE,
                                         correction = TRUE, conservative = TRUE) {
  epilink_id <- paste0("epilink_", family)
  glm_id <- paste0("glm_", family)
  registry <- empirical_null_methods()
  epilink_label <- registry$method[match(epilink_id, registry$method_id)]
  glm_label <- registry$method[match(glm_id, registry$method_id)]
  fitted <- EpiLink(
    object = object, gene_peak_df = pairs[, c("gene", "peak")],
    family = family, celltype = NULL, covariates = covariates,
    null_generation = "s", stat = "zval", bin = binarize_atac,
    selection = selection, returnObject = TRUE, fdrcutoff = alpha,
    ncore = n_cores, theory = FALSE, correction = correction,
    conservative = conservative, verbose = FALSE, seed = seed
  )
  if (is.null(fitted)) {
    return(list(
      epilink = empty_empirical_result(pairs, epilink_id, epilink_label,
                                       "EpiLink returned no fitted pairs"),
      glm = empty_empirical_result(pairs, glm_id, glm_label,
                                   "EpiLink returned no fitted pairs")
    ))
  }
  results <- fitted@misc$EpiLink$res_df
  significant <- EpiLinkFilter(
    fitted, stat = "zval", fdrcutoff = alpha, correction = correction,
    conservative = conservative, returnObject = FALSE, verbose = FALSE
  )
  significant_keys <- empirical_pair_key(significant)
  list(
    epilink = align_empirical_result(
      results, pairs, epilink_id, epilink_label, alpha,
      operational_keys = significant_keys
    ),
    glm = align_empirical_result(results, pairs, glm_id, glm_label, alpha)
  )
}

#' Fit scMultiMap on one empirical-null replicate
fit_empirical_scmultimap <- function(object, pairs, alpha, irls = TRUE) {
  results <- suppressMessages(run_scMultiMap_compat(
    object, pairs[, c("gene", "peak")], gene_assay = "RNA",
    peak_assay = "ATAC", bsample = NULL, irls = irls, verbose = FALSE
  ))
  align_empirical_result(
    results, pairs, "scmultimap", "scMultiMap", alpha
  )
}

#' Fit SCENT on one empirical-null replicate.
#'
#' SCENT receives the same fixed candidate-pair universe and the same permuted
#' RNA/ATAC object as every other method. Its bootstrap p-value is used as the
#' method p-value, followed by BH correction over the fixed candidate pairs.
fit_empirical_scent <- function(object, pairs, alpha, seed,
                                n_cores = 1L, covariates = character()) {
  if (!requireNamespace("SCENT", quietly = TRUE)) {
    stop("SCENT is not installed")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is not installed")
  }
  rna <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
  atac <- Seurat::GetAssayData(object, assay = "ATAC", layer = "counts")
  cell_names <- colnames(rna)
  metadata <- object@meta.data[cell_names, , drop = FALSE]
  metadata$cell <- cell_names
  metadata$celltypes <- rep("all", length(cell_names))
  rownames(metadata) <- cell_names
  if (length(covariates)) {
    for (covariate in covariates) {
      if (grepl("^log\\(.+\\)$", covariate)) {
        source_name <- sub("^log\\((.+)\\)$", "\\1", covariate)
        if (!source_name %in% names(metadata)) {
          stop("SCENT covariate source is missing: ", source_name)
        }
        metadata[[source_name]] <- log(pmax(as.numeric(metadata[[source_name]]), 1))
      } else if (!covariate %in% names(metadata)) {
        stop("SCENT covariate is missing: ", covariate)
      }
    }
  }
  peak_info <- pairs[, c("gene", "peak"), drop = FALSE]
  scent_object <- SCENT::CreateSCENTObj(
    rna = Matrix::Matrix(rna, sparse = TRUE),
    atac = Matrix::Matrix(atac, sparse = TRUE),
    meta.data = metadata, peak.info = peak_info,
    covariates = as.character(covariates), celltypes = "celltypes"
  )
  scent_object <- with_epilink_seed(as.integer(seed), SCENT::SCENT_algorithm(
    scent_object, celltype = "all", ncores = as.integer(n_cores),
    regr = "poisson", bin = TRUE
  ))
  results <- scent_object@SCENT.result
  if (is.null(results) || !nrow(results)) {
    return(empty_empirical_result(pairs, "scent", "SCENT"))
  }
  if (ncol(results) == 7L) {
    names(results) <- c("gene", "peak", "beta", "se", "z", "p", "boot_p")
  }
  if (!"boot_p" %in% names(results)) {
    stop("SCENT result has no boot_basic_p column")
  }
  scent_results <- data.frame(
    gene = as.character(results$gene),
    peak = as.character(results$peak),
    statistic = as.numeric(results$z),
    pval = as.numeric(results$boot_p),
    stringsAsFactors = FALSE
  )
  align_empirical_result(scent_results, pairs, "scent", "SCENT", alpha)
}

#' Fit native Signac LinkPeaks on one empirical-null replicate
#'
#' LinkPeaks uses normalized RNA data and its usual matched-background null.
#' Results are then restricted to the fixed 1000 benchmark candidate pairs.
fit_empirical_signac <- function(object, pairs, alpha, seed, n_cores = 1L) {
  object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  old_plan <- NULL
  if (requireNamespace("future", quietly = TRUE) && n_cores > 1L &&
      .Platform$OS.type != "windows") {
    old_plan <- future::plan()
    old_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 8 * 1024^3)
    on.exit({
      options(future.globals.maxSize = old_size)
      future::plan(old_plan)
    }, add = TRUE)
    future::plan(future::multicore, workers = as.integer(n_cores))
  }
  run_linkpeaks <- function() Signac::LinkPeaks(
    object, peak.assay = "ATAC", expression.assay = "RNA",
    genes.use = unique(pairs$gene), pvalue_cutoff = 1, score_cutoff = 0,
    verbose = FALSE
  )
  linked <- with_epilink_seed(as.integer(seed), try(run_linkpeaks(), silent = TRUE))
  if (inherits(linked, "try-error")) {
    if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
    linked <- with_epilink_seed(as.integer(seed), run_linkpeaks())
  }
  links <- as.data.frame(Signac::Links(linked[["ATAC"]]))
  if (!nrow(links)) {
    return(empty_empirical_result(
      pairs, "signac", "Signac", "LinkPeaks returned no links"
    ))
  }
  links$pval <- pmin(2 * links$pvalue, 1)
  # BH is computed after restricting to the same fixed candidate universe used
  # by every method, rather than over incidental cis links returned by LinkPeaks.
  align_empirical_result(links, pairs, "signac", "Signac", alpha)
}

#' Fit one configured method on one empirical-null dataset
fit_empirical_null_method <- function(object, pairs, method_id, alpha,
                                      covariates, seed, n_cores = 1L,
                                      selection = "rate", binarize_atac = TRUE,
                                      epilink_correction = TRUE,
                                      epilink_conservative = TRUE,
                                      scmultimap_irls = TRUE,
                                      scent_cores = 1L) {
  registry <- empirical_null_methods()
  if (!method_id %in% registry$method_id) {
    stop("Unknown empirical-null method: ", method_id, call. = FALSE)
  }
  if (method_id %in% c("epilink_poisson", "epilink_nb")) {
    family <- sub("^epilink_", "", method_id)
    fitted <- fit_empirical_epilink_family(
      object, pairs, family = family, alpha = alpha,
      covariates = covariates, n_cores = n_cores, seed = seed,
      selection = selection, binarize_atac = binarize_atac,
      correction = epilink_correction,
      conservative = epilink_conservative
    )
    return(fitted$epilink)
  }
  if (method_id %in% c("glm_poisson", "glm_nb")) {
    family <- sub("^glm_", "", method_id)
    return(fit_empirical_glm(
      object, pairs, family = family, alpha = alpha,
      covariates = covariates, n_cores = n_cores,
      selection = selection, binarize_atac = binarize_atac, seed = seed
    ))
  }
  if (method_id == "scmultimap") {
    return(fit_empirical_scmultimap(
      object, pairs, alpha = alpha, irls = scmultimap_irls
    ))
  }
  if (method_id == "scent") {
    return(fit_empirical_scent(
      object, pairs, alpha = alpha, seed = seed,
      n_cores = scent_cores, covariates = covariates
    ))
  }
  fit_empirical_signac(
    object, pairs, alpha = alpha, seed = seed, n_cores = n_cores
  )
}

#' Return the result files produced by one independent method job
empirical_null_method_group <- function(method_id) {
  if (method_id %in% c("epilink_poisson", "glm_poisson")) {
    return(c("epilink_poisson", "glm_poisson"))
  }
  if (method_id %in% c("epilink_nb", "glm_nb")) {
    return(c("epilink_nb", "glm_nb"))
  }
  method_id
}

#' Fit one independent job, returning one or two method result tables
fit_empirical_null_method_group <- function(object, pairs, method_id, alpha,
                                            covariates, seed, n_cores = 1L,
                                            selection = "rate",
                                            binarize_atac = TRUE,
                                            epilink_correction = TRUE,
                                            epilink_conservative = TRUE,
                                            scmultimap_irls = TRUE,
                                            scent_cores = 1L) {
  group <- empirical_null_method_group(method_id)
  if (length(group) == 2L && grepl("^epilink_", group[[1]])) {
    family <- sub("^epilink_", "", group[[1]])
    fitted <- fit_empirical_epilink_family(
      object, pairs, family = family, alpha = alpha,
      covariates = covariates, n_cores = n_cores, seed = seed,
      selection = selection, binarize_atac = binarize_atac,
      correction = epilink_correction,
      conservative = epilink_conservative
    )
    return(setNames(list(fitted$epilink, fitted$glm), group))
  }
  one_result <- fit_empirical_null_method(
    object, pairs, method_id = group[[1]], alpha = alpha,
    covariates = covariates, seed = seed, n_cores = n_cores,
    selection = selection, binarize_atac = binarize_atac,
    epilink_correction = epilink_correction,
    epilink_conservative = epilink_conservative,
    scmultimap_irls = scmultimap_irls,
    scent_cores = scent_cores
  )
  setNames(list(one_result), group)
}

empirical_null_folder <- function(replicate, replicates_per_folder = 100L) {
  sprintf("%02d", ((as.integer(replicate) - 1L) %/%
                      as.integer(replicates_per_folder)) + 1L)
}

empirical_null_result_path <- function(output_dir, replicate, method_id,
                                       replicates_per_folder = 100L) {
  file.path(
    output_dir,
    empirical_null_folder(replicate, replicates_per_folder),
    sprintf("replicate_%04d_%s.tsv.gz", as.integer(replicate), method_id)
  )
}

#' Summarize one pair-level method result
summarize_empirical_method <- function(result, replicate, runtime_seconds) {
  tested <- is.finite(result$pval)
  data.frame(
    replicate = as.integer(replicate),
    method_id = result$method_id[[1]],
    method = result$method[[1]],
    n_pairs = nrow(result),
    n_tested = sum(is.finite(result$pval)),
    raw_type1 = if (any(tested)) mean(result$raw_rejected[tested]) else NA_real_,
    operational_rejection_rate = if (any(tested)) {
      mean(result$operational_rejected[tested])
    } else NA_real_,
    runtime_seconds = runtime_seconds,
    error = {
      values <- unique(stats::na.omit(result$error))
      if (length(values)) paste(values, collapse = "; ") else NA_character_
    },
    stringsAsFactors = FALSE
  )
}
