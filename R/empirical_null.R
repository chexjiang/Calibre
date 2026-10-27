# Utilities for the PBMC empirical-null type-I error benchmark.

#' Empirical-null benchmark method registry
#' @return A data frame with stable method identifiers and labels.
empirical_null_methods <- function() {
  data.frame(
    method_id = c(
      "calibre_poisson", "calibre_nb", "glm_poisson", "glm_nb",
      "scmultimap", "signac", "scent", "archr"
    ),
    method = c(
      "Calibre-Poisson", "Calibre-NB", "Poisson-GLM", "NB-GLM",
      "scMultiMap", "Signac", "SCENT", "ArchR"
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
  permutation <- with_calibre_seed(as.integer(seed), sample.int(ncol(atac)))
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
  chosen <- with_calibre_seed(
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
    c("zval", "test_stat", "zscore", "score", "cor", "Correlation", "covar"),
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
  results <- Calibre(
    object = object, gene_peak_df = pairs[, c("gene", "peak")],
    family = family, celltype = NULL, covariates = covariates,
    null_generation = "s", stat = "zval", bin = binarize_atac,
    selection = selection, returnObject = FALSE, fdrcutoff = alpha,
    ncore = n_cores, theory = TRUE, correction = FALSE,
    conservative = FALSE, verbose = FALSE, seed = seed
  )
  align_empirical_result(results, pairs, method_id, method, alpha)
}

#' Fit Calibre and its corresponding BH-GLM on one null replicate
#'
#' A single Calibre fit supplies the observed GLM p-values used by both rows,
#' avoiding a redundant second model fit.
fit_empirical_calibre_family <- function(object, pairs, family, alpha,
                                         covariates, n_cores, seed,
                                         selection = "rate", binarize_atac = TRUE,
                                         correction = TRUE, conservative = TRUE) {
  calibre_id <- paste0("calibre_", family)
  glm_id <- paste0("glm_", family)
  registry <- empirical_null_methods()
  calibre_label <- registry$method[match(calibre_id, registry$method_id)]
  glm_label <- registry$method[match(glm_id, registry$method_id)]
  fitted <- Calibre(
    object = object, gene_peak_df = pairs[, c("gene", "peak")],
    family = family, celltype = NULL, covariates = covariates,
    null_generation = "s", stat = "zval", bin = binarize_atac,
    selection = selection, returnObject = TRUE, fdrcutoff = alpha,
    ncore = n_cores, theory = FALSE, correction = correction,
    conservative = conservative, verbose = FALSE, seed = seed
  )
  if (is.null(fitted)) {
    empty_calibre <- empty_empirical_result(
      pairs, calibre_id, calibre_label, "Calibre returned no fitted pairs"
    )
    return(list(
      calibre = empty_calibre,
      epilink = empty_calibre,
      glm = empty_empirical_result(pairs, glm_id, glm_label,
                                   "Calibre returned no fitted pairs")
    ))
  }
  results <- fitted@misc$Calibre$res_df
  significant <- CalibreFilter(
    fitted, stat = "zval", fdrcutoff = alpha, correction = correction,
    conservative = conservative, returnObject = FALSE, verbose = FALSE
  )
  significant_keys <- empirical_pair_key(significant)
  calibre_result <- align_empirical_result(
      results, pairs, calibre_id, calibre_label, alpha,
      operational_keys = significant_keys
    )
  list(
    calibre = calibre_result,
    epilink = calibre_result,
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
  scent_object <- with_calibre_seed(as.integer(seed), SCENT::SCENT_algorithm(
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
  linked <- with_calibre_seed(as.integer(seed), try(run_linkpeaks(), silent = TRUE))
  if (inherits(linked, "try-error")) {
    if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
    linked <- with_calibre_seed(as.integer(seed), run_linkpeaks())
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

#' Cache only the selected cells' original fragments.
#'
#' Unlike the count-based methods, ArchR cannot consume a permuted peak-count
#' matrix.  It builds an Arrow file from a bgzip/tabix fragment file, so the
#' empirical null is generated by remapping the fourth (barcode) column of
#' the original fragments.  This cache avoids scanning the full PBMC fragment
#' file once for every replicate.
prepare_empirical_fragment_cache <- function(fragment_file, selected_cells,
                                             cache_file, cells_file = NULL,
                                             overwrite = FALSE) {
  selected_cells <- as.character(selected_cells)
  if (anyDuplicated(selected_cells)) {
    stop("Selected fragment barcodes must be unique", call. = FALSE)
  }
  cache_index <- paste0(cache_file, ".tbi")
  cache_ok <- file.exists(cache_file) && file.exists(cache_index) &&
    !isTRUE(overwrite)
  if (cache_ok && !is.null(cells_file) && file.exists(cells_file)) {
    cache_ok <- identical(as.character(readRDS(cells_file)), selected_cells)
  }
  if (cache_ok) return(cache_file)

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  plain_file <- paste0(cache_file, ".tmp-", Sys.getpid(), ".tsv")
  if (file.exists(cache_file)) unlink(cache_file)
  if (file.exists(cache_index)) unlink(cache_index)
  on.exit(if (file.exists(plain_file)) unlink(plain_file), add = TRUE)

  input <- if (grepl("[.]gz$", fragment_file, ignore.case = TRUE)) {
    gzfile(fragment_file, open = "rt")
  } else file(fragment_file, open = "rt")
  output <- file(plain_file, open = "wt")
  on.exit(try(close(input), silent = TRUE), add = TRUE)
  on.exit(try(close(output), silent = TRUE), add = TRUE)
  selected_lookup <- selected_cells
  repeat {
    lines <- readLines(input, n = 100000L)
    if (!length(lines)) break
    comment <- startsWith(lines, "#")
    keep <- comment
    row_index <- which(!comment)
    if (length(row_index)) {
      fields <- strsplit(lines[row_index], "\t", fixed = TRUE)
      has_barcode <- vapply(fields, length, integer(1)) >= 4L
      keep[row_index[has_barcode]] <- vapply(
        fields[has_barcode], function(x) x[[4L]] %in% selected_lookup,
        logical(1)
      )
    }
    writeLines(lines[keep], output)
  }
  close(input)
  close(output)
  gz_file <- Rsamtools::bgzip(plain_file, dest = cache_file, overwrite = TRUE)
  Rsamtools::indexTabix(gz_file, format = "bed")
  if (!is.null(cells_file)) saveRDS(selected_cells, cells_file)
  cache_file
}

#' Write one ArchR-compatible empirical-null fragment file.
#' @param source_to_target Named character vector, names are source barcodes
#'   and values are target barcodes.  This is the inverse orientation of the
#'   permutation matrix: target[j] receives source[permutation[j]].
write_empirical_fragment_replicate <- function(cache_file, source_to_target,
                                               output_file) {
  if (!file.exists(cache_file) || !file.exists(paste0(cache_file, ".tbi"))) {
    stop("Cached empirical fragments are missing or not indexed", call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  output_index <- paste0(output_file, ".tbi")
  if (file.exists(output_file)) unlink(output_file)
  if (file.exists(output_index)) unlink(output_index)
  plain_file <- paste0(output_file, ".tmp-", Sys.getpid(), ".tsv")
  on.exit(if (file.exists(plain_file)) unlink(plain_file), add = TRUE)
  input <- gzfile(cache_file, open = "rt")
  output <- file(plain_file, open = "wt")
  on.exit(try(close(input), silent = TRUE), add = TRUE)
  on.exit(try(close(output), silent = TRUE), add = TRUE)
  repeat {
    lines <- readLines(input, n = 100000L)
    if (!length(lines)) break
    comment <- startsWith(lines, "#")
    row_index <- which(!comment)
    if (length(row_index)) {
      fields <- strsplit(lines[row_index], "\t", fixed = TRUE)
      for (i in seq_along(fields)) {
        barcode <- fields[[i]][[4L]]
        replacement <- unname(source_to_target[barcode])
        if (length(replacement) == 1L && !is.na(replacement)) {
          fields[[i]][[4L]] <- replacement
          lines[row_index[[i]]] <- paste(fields[[i]], collapse = "\t")
        }
      }
    }
    writeLines(lines, output)
  }
  close(input)
  close(output)
  gz_file <- Rsamtools::bgzip(plain_file, dest = output_file, overwrite = TRUE)
  Rsamtools::indexTabix(gz_file, format = "bed")
  output_file
}

#' Run ArchR on one fragment-level empirical-null replicate.
fit_empirical_archr <- function(object, pairs, alpha, seed, n_cores = 1L,
                                fragment_file, sample_name, genome_name,
                                work_directory) {
  if (!exists("run_archr_real", mode = "function")) {
    stop("run_archr_real() is not available; source R/archr_benchmark.R")
  }
  archr_requirements()
  rna <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
  atac <- Seurat::GetAssayData(object, assay = "ATAC", layer = "counts")
  peak_names <- rownames(atac)
  if (is.null(peak_names) && methods::is(atac, "Matrix")) {
    peak_names <- methods::slot(atac, "Dimnames")[[1L]]
  }
  if (is.null(peak_names)) peak_names <- dimnames(atac)[[1L]]
  if (is.null(peak_names) || !length(peak_names)) {
    stop("ArchR empirical-null ATAC matrix has no peak names", call. = FALSE)
  }
  rownames(atac) <- peak_names
  colnames(atac) <- colnames(rna)
  linked <- run_archr_real(
    rna = rna, atac = atac, pairs = pairs,
    keep_cells = colnames(rna), fragment_file = fragment_file,
    sample_name = sample_name, genome_name = genome_name,
    output_dir = work_directory, fdr = alpha,
    n_cores = n_cores, seed = seed
  )
  archr_results <- linked$results
  # ArchR versions differ in whether addPeak2GeneLinks() emits Pval when
  # empirical/permuted p-values are disabled.  FDR is still its native
  # adjusted significance measure; use it as a finite audit p-value only when
  # no raw Pval is available, while retaining native FDR in padj.
  if (!any(is.finite(archr_results$pval)) && "FDR" %in% names(archr_results)) {
    archr_results$pval <- archr_results$FDR
  }
  significant_keys <- empirical_pair_key(linked$significant)
  aligned <- align_empirical_result(
    archr_results, pairs, "archr", "ArchR", alpha,
    operational_keys = significant_keys
  )
  aligned
}

#' Fit one configured method on one empirical-null dataset
fit_empirical_null_method <- function(object, pairs, method_id, alpha,
                                      covariates, seed, n_cores = 1L,
                                      selection = "rate", binarize_atac = TRUE,
                                      calibre_correction = TRUE,
                                      calibre_conservative = TRUE,
                                      scmultimap_irls = TRUE,
                                      scent_cores = 1L,
                                      archr_fragment_file = NULL,
                                      archr_sample_name = "pbmc",
                                      archr_genome_name = "hg38",
                                      archr_work_directory = NULL) {
  method_id <- calibre_output_method_id(method_id)
  registry <- empirical_null_methods()
  if (!method_id %in% registry$method_id) {
    stop("Unknown empirical-null method: ", method_id, call. = FALSE)
  }
  if (method_id %in% c("calibre_poisson", "calibre_nb")) {
    family <- sub("^calibre_", "", method_id)
    fitted <- fit_empirical_calibre_family(
      object, pairs, family = family, alpha = alpha,
      covariates = covariates, n_cores = n_cores, seed = seed,
      selection = selection, binarize_atac = binarize_atac,
      correction = calibre_correction,
      conservative = calibre_conservative
    )
    return(fitted$calibre)
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
  if (method_id == "archr") {
    return(fit_empirical_archr(
      object, pairs, alpha = alpha, seed = seed, n_cores = n_cores,
      fragment_file = archr_fragment_file,
      sample_name = archr_sample_name, genome_name = archr_genome_name,
      work_directory = archr_work_directory
    ))
  }
  fit_empirical_signac(
    object, pairs, alpha = alpha, seed = seed, n_cores = n_cores
  )
}

#' Return the result files produced by one independent method job
empirical_null_method_group <- function(method_id) {
  method_id <- calibre_output_method_id(method_id)
  if (method_id %in% c("calibre_poisson", "glm_poisson")) {
    return(c("calibre_poisson", "glm_poisson"))
  }
  if (method_id %in% c("calibre_nb", "glm_nb")) {
    return(c("calibre_nb", "glm_nb"))
  }
  method_id
}

#' Fit one independent job, returning one or two method result tables
fit_empirical_null_method_group <- function(object, pairs, method_id, alpha,
                                            covariates, seed, n_cores = 1L,
                                            selection = "rate",
                                            binarize_atac = TRUE,
                                            calibre_correction = TRUE,
                                            calibre_conservative = TRUE,
                                            scmultimap_irls = TRUE,
                                            scent_cores = 1L,
                                            archr_fragment_file = NULL,
                                            archr_sample_name = "pbmc",
                                            archr_genome_name = "hg38",
                                            archr_work_directory = NULL) {
  method_id <- calibre_output_method_id(method_id)
  group <- empirical_null_method_group(method_id)
  if (length(group) == 2L && grepl("^calibre_", group[[1]])) {
    family <- sub("^calibre_", "", group[[1]])
    fitted <- fit_empirical_calibre_family(
      object, pairs, family = family, alpha = alpha,
      covariates = covariates, n_cores = n_cores, seed = seed,
      selection = selection, binarize_atac = binarize_atac,
      correction = calibre_correction,
      conservative = calibre_conservative
    )
    return(setNames(list(fitted$calibre, fitted$glm), group))
  }
  one_result <- fit_empirical_null_method(
    object, pairs, method_id = group[[1]], alpha = alpha,
    covariates = covariates, seed = seed, n_cores = n_cores,
    selection = selection, binarize_atac = binarize_atac,
    calibre_correction = calibre_correction,
    calibre_conservative = calibre_conservative,
    scmultimap_irls = scmultimap_irls,
    scent_cores = scent_cores,
    archr_fragment_file = archr_fragment_file,
    archr_sample_name = archr_sample_name,
    archr_genome_name = archr_genome_name,
    archr_work_directory = archr_work_directory
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
