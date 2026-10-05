#!/usr/bin/env Rscript

# Run EpiLinkDiff on one real-data cell type.
#
# Usage:
#   Rscript scripts/run_real_diff.R <dataset> <celltype>
#       [epilink_poisson|epilink_nb|scmultimap]
#       [--donor-effect] [--donor-robust|--no-donor-robust] [--overwrite]
#
# Results are written to results/real/<dataset>_diff/<celltype>/<method>/,
# using the same all_pairs/significant_pairs/candidate_pairs/run_metadata/run.log
# layout as scripts/run_real.R.

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
donor_effect <- "--donor-effect" %in% args
has_donor_robust <- "--donor-robust" %in% args
no_donor_robust <- "--no-donor-robust" %in% args
if (has_donor_robust && no_donor_robust) {
  stop("--donor-robust and --no-donor-robust are mutually exclusive",
       call. = FALSE)
}
donor_robust <- if (has_donor_robust) TRUE else if (no_donor_robust) FALSE else NULL
args <- setdiff(args, "--overwrite")
args <- setdiff(args, c("--donor-effect", "--donor-robust", "--no-donor-robust"))
if (length(args) < 2L || length(args) > 3L) {
  stop(
    "Usage: Rscript scripts/run_real_diff.R <dataset> <celltype> ",
    "[epilink_poisson|epilink_nb|scmultimap] ",
    "[--donor-effect] [--donor-robust|--no-donor-robust] [--overwrite]",
    call. = FALSE
  )
}
dataset_name <- args[[1L]]
celltype <- args[[2L]]
method_id <- if (length(args) == 3L) args[[3L]] else NULL
if (!is.null(method_id) &&
    !method_id %in% c("epilink_poisson", "epilink_nb", "scmultimap")) {
  stop("unknown method: ", method_id,
       "; use epilink_poisson, epilink_nb, or scmultimap", call. = FALSE)
}

file_argument <- sub(
  "^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
)
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml is required", call. = FALSE)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}
configuration <- yaml::read_yaml("config/datasets.yml")
if (!dataset_name %in% names(configuration)) {
  stop("unknown dataset: ", dataset_name, call. = FALSE)
}
settings <- utils::modifyList(configuration$defaults, configuration[[dataset_name]])

condition_column <- settings$diff_condition_column
donor_column <- settings$diff_donor_column
if (is.null(condition_column)) {
  stop(
    "dataset ", dataset_name,
    " needs diff_condition_column in config/datasets.yml",
    call. = FALSE
  )
}
diff_covariates <- if (is.null(settings$diff_covariates)) {
  settings$covariates
} else {
  settings$diff_covariates
}
diff_covariates <- unlist(diff_covariates, use.names = FALSE)
if (is.null(diff_covariates)) diff_covariates <- character()
# The donor effect is represented by `donor_col`; including the same metadata
# column again as a covariate would be redundant/collinear.
if (!is.null(donor_column)) {
  diff_covariates <- setdiff(diff_covariates, donor_column)
}
diff_family <- if (identical(method_id, "scmultimap")) {
  NA_character_
} else if (!is.null(method_id)) {
  if (method_id == "epilink_nb") "nb" else "poisson"
} else if (is.null(settings$diff_family)) {
  settings$family
} else {
  settings$diff_family
}
diff_null_generation <- if (is.null(settings$diff_null_generation)) {
  settings$null_generation
} else {
  settings$diff_null_generation
}
diff_bin <- if (is.null(settings$diff_binarize_atac)) {
  isTRUE(settings$binarize_atac)
} else {
  isTRUE(settings$diff_binarize_atac)
}
diff_fdr <- if (is.null(settings$diff_fdr)) settings$fdr else settings$diff_fdr
diff_correction <- if (is.null(settings$diff_correction)) {
  isTRUE(settings$correction)
} else {
  isTRUE(settings$diff_correction)
}

output_root <- file.path("results", "real", paste0(dataset_name, "_diff"))
effective_donor_robust <- if (is.null(donor_robust)) {
  !is.null(donor_column)
} else {
  donor_robust
}
output_method <- if (method_id == "scmultimap") {
  "scmultimap"
} else {
  paste0(
    "epilink_diff_", diff_family,
    if (donor_effect) "_donor_effect" else "_no_donor_effect",
    if (effective_donor_robust) "_donor_robust" else "_no_donor_robust"
  )
}
output_dir <- file.path(output_root, celltype, output_method)
done_file <- file.path(output_dir, "run_metadata.tsv")
if (file.exists(done_file) && !overwrite) {
  stop(
    "completed output exists: ", output_dir,
    "; use --overwrite to replace",
    call. = FALSE
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(output_dir, "run.log")
log_message <- function(...) {
  line <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    paste0(..., collapse = "")
  )
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}
started <- Sys.time()
log_message("Loading ", dataset_name, " inputs for ", celltype)

rna <- readRDS(settings$rna)
atac <- readRDS(settings$atac)
metadata <- readRDS(settings$metadata)
pairs_all <- readRDS(settings$candidate_pairs)
pairs_all$peak <- CanonicalPeakNames(pairs_all$peak)

if (is.null(rownames(rna)) && !is.null(settings$rna_features) &&
    file.exists(settings$rna_features)) {
  feature_table <- utils::read.delim(
    settings$rna_features, stringsAsFactors = FALSE,
    check.names = FALSE, row.names = 1
  )
  if (nrow(feature_table) != nrow(rna)) {
    stop("RNA feature metadata length does not match RNA matrix", call. = FALSE)
  }
  rownames(rna) <- rownames(feature_table)
}
if (!identical(rownames(metadata), colnames(rna)) &&
    !is.null(settings$metadata_barcode_column)) {
  barcode_column <- settings$metadata_barcode_column
  if (!barcode_column %in% names(metadata)) {
    stop("metadata barcode column not found: ", barcode_column, call. = FALSE)
  }
  rownames(metadata) <- as.character(metadata[[barcode_column]])
}
if (is.null(rownames(rna))) stop("RNA matrix has no feature names", call. = FALSE)
rna_cells <- as.character(colnames(rna))
atac_cells <- as.character(colnames(atac))
if (!setequal(rna_cells, atac_cells)) {
  stop("RNA and ATAC contain different cell barcodes", call. = FALSE)
}
if (!identical(rna_cells, atac_cells)) {
  atac <- atac[, rna_cells, drop = FALSE]
}
metadata <- metadata[rna_cells, , drop = FALSE]
if (!condition_column %in% names(metadata)) {
  stop("condition column not found: ", condition_column, call. = FALSE)
}
if (!is.null(donor_column) && !donor_column %in% names(metadata)) {
  stop("donor column not found: ", donor_column, call. = FALSE)
}
if ((donor_effect || effective_donor_robust) && is.null(donor_column)) {
  stop("donor_col is required for the selected donor options", call. = FALSE)
}
if (!settings$celltype_column %in% names(metadata)) {
  stop("celltype column not found: ", settings$celltype_column, call. = FALSE)
}

celltype_values <- as.character(metadata[[settings$celltype_column]])
keep_mask <- !is.na(celltype_values) & celltype_values == celltype
keep_cells <- rownames(metadata)[keep_mask]
if (!length(keep_cells)) stop("no cells found for ", celltype, call. = FALSE)
rna <- rna[, keep_cells, drop = FALSE]
atac <- atac[, keep_cells, drop = FALSE]
metadata <- metadata[keep_cells, , drop = FALSE]

# Brain2 is stored with raw Status labels. Make the reference/comparison order
# explicit whenever the usual Ctrl/AD labels are present.
condition_values <- unique(as.character(metadata[[condition_column]]))
if (all(c("Ctrl", "AD") %in% condition_values)) {
  metadata[[condition_column]] <- factor(
    metadata[[condition_column]], levels = c("Ctrl", "AD")
  )
}

if (!requireNamespace("Seurat", quietly = TRUE) ||
    !requireNamespace("Signac", quietly = TRUE)) {
  stop("Seurat and Signac are required", call. = FALSE)
}
object <- Seurat::CreateSeuratObject(
  counts = rna, meta.data = metadata, assay = "RNA"
)
peak_sep <- InferPeakSeparator(rownames(atac))
object[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac, sep = peak_sep
)

# Use the same default candidate filtering as EpiLink. This intentionally does
# not impose scMultiMap-specific top-gene/top-peak restrictions.
pairs <- FilterGenePeakPairs(
  object, pairs_all, selection = settings$selection, verbose = TRUE
)
pairs <- unique(pairs[, c("gene", "peak"), drop = FALSE])
if (!nrow(pairs)) stop("no candidate pairs remain after filtering", call. = FALSE)
saveRDS(pairs, file.path(output_dir, "candidate_pairs.rds"))
log_message("Prepared ", ncol(object), " cells and ", nrow(pairs), " candidate pairs")
log_message(
  "condition=", condition_column, "; donor=",
  if (is.null(donor_column)) "none" else donor_column,
  "; donor_effect=", donor_effect,
  "; donor_robust=", effective_donor_robust,
  "; bin=", if (identical(method_id, "scmultimap")) "not_used" else diff_bin,
  "; family=", ifelse(is.na(diff_family), "scMultiMap", diff_family)
)

condition_levels <- if (is.factor(metadata[[condition_column]])) {
  levels(droplevels(metadata[[condition_column]]))
} else {
  unique(as.character(metadata[[condition_column]]))
}
if (length(condition_levels) != 2L) {
  stop("differential analysis requires exactly two condition labels", call. = FALSE)
}

if (identical(method_id, "scmultimap")) {
  if (!requireNamespace("scMultiMap", quietly = TRUE)) {
    stop("scMultiMap is required for method=scmultimap", call. = FALSE)
  }
  if (is.null(donor_column)) {
    stop("scMultiMap disease-control analysis requires diff_donor_column",
         call. = FALSE)
  }

  # scMultiMap's disease-control workflow estimates covariances separately in
  # the two groups, then permutes biological-sample labels while preserving
  # the observed number of disease samples. Use the configured donor
  # column as the biological-sample label, as in the vignette.
  sample_labels <- as.character(metadata[[donor_column]])
  sample_condition <- split(
    as.character(metadata[[condition_column]]), sample_labels
  )
  if (any(vapply(sample_condition, function(x) length(unique(x)) != 1L,
                 logical(1L)))) {
    stop("each donor must belong to exactly one condition for scMultiMap",
         call. = FALSE)
  }
  sample_ids <- names(sample_condition)
  sample_status <- vapply(sample_condition, function(x) unique(x)[[1L]],
                          character(1L))
  control_ids <- sample_ids[sample_status == condition_levels[[1L]]]
  disease_ids <- sample_ids[sample_status == condition_levels[[2L]]]
  if (!length(control_ids) || !length(disease_ids)) {
    stop("both conditions must contain at least one donor", call. = FALSE)
  }

  irls_enabled <- if (is.null(settings$scmultimap_irls)) TRUE else
    isTRUE(settings$scmultimap_irls)
  if (!irls_enabled) {
    # Some scMultiMap releases leave theta_median undefined when IRLS is
    # disabled. Apply the same local namespace patch as scripts/run_real.R.
    scmultimap_ns <- asNamespace("scMultiMap")
    scmultimap_irls_original <- get("scMultiMap_IRLS", envir = scmultimap_ns)
    scmultimap_irls_fixed <- scmultimap_irls_original
    scmultimap_body <- as.list(body(scmultimap_irls_fixed))
    irls_branch <- which(vapply(scmultimap_body, function(expr) {
      is.call(expr) && identical(expr[[1L]], as.name("if")) &&
        identical(expr[[2L]], as.name("irls"))
    }, logical(1L)))
    if (length(irls_branch) != 1L) {
      stop("could not patch scMultiMap_IRLS for irls=FALSE", call. = FALSE)
    }
    if_expr <- scmultimap_body[[irls_branch]]
    else_expr <- as.list(if_expr[[4L]])
    if_expr[[4L]] <- as.call(c(
      list(as.name("{")),
      list(quote(theta_median <- rep(NA_real_, n_cell))),
      else_expr[-1L]
    ))
    scmultimap_body[[irls_branch]] <- if_expr
    body(scmultimap_irls_fixed) <- as.call(scmultimap_body)
    lock_after_patch <- bindingIsLocked("scMultiMap_IRLS", scmultimap_ns)
    if (lock_after_patch) unlockBinding("scMultiMap_IRLS", scmultimap_ns)
    assign("scMultiMap_IRLS", scmultimap_irls_fixed, envir = scmultimap_ns)
    if (lock_after_patch) lockBinding("scMultiMap_IRLS", scmultimap_ns)
  }

  run_scmultimap <- function(sample_ids_subset) {
    subset_cells <- rownames(metadata)[sample_labels %in% sample_ids_subset]
    scMultiMap::scMultiMap(
      object[, subset_cells], pairs,
      gene_assay = "RNA", peak_assay = "ATAC",
      bsample = donor_column, p.adjust.method = "BH",
      irls = irls_enabled, verbose = TRUE
    )
  }
  align_scmultimap <- function(result) {
    result_key <- paste(result$gene, result$peak, sep = "\r")
    pair_key <- paste(pairs$gene, pairs$peak, sep = "\r")
    aligned <- result[match(pair_key, result_key), , drop = FALSE]
    aligned$gene <- pairs$gene
    aligned$peak <- pairs$peak
    aligned
  }

  log_message(
    "Running scMultiMap disease-control analysis with ",
    length(disease_ids), " disease and ", length(control_ids),
    " control donors"
  )
  control_raw <- run_scmultimap(control_ids)
  disease_raw <- run_scmultimap(disease_ids)
  control_result <- align_scmultimap(control_raw)
  disease_result <- align_scmultimap(disease_raw)
  observed_diff <- disease_result$covar - control_result$covar

  all_combinations <- utils::combn(sample_ids, length(disease_ids))
  n_permutation <- min(100L, ncol(all_combinations))
  set.seed(settings$seed)
  permutation_indices <- sample.int(
    ncol(all_combinations), n_permutation, replace = FALSE
  )
  null_diff <- matrix(
    NA_real_, nrow = nrow(pairs), ncol = n_permutation
  )
  for (i_permutation in seq_len(n_permutation)) {
    disease_ids_null <- all_combinations[
      , permutation_indices[[i_permutation]], drop = TRUE
    ]
    control_ids_null <- setdiff(sample_ids, disease_ids_null)
    null_control <- align_scmultimap(run_scmultimap(control_ids_null))
    null_disease <- align_scmultimap(run_scmultimap(disease_ids_null))
    null_diff[, i_permutation] <-
      null_disease$covar - null_control$covar
    log_message("Permutation ", i_permutation, "/", n_permutation)
  }

  # Match the scMultiMap disease-control vignette: two-sided p-values are
  # twice the smaller one-sided proportion, without a +1 correction.
  p_upper <- rowMeans(observed_diff > null_diff, na.rm = TRUE)
  p_lower <- rowMeans(observed_diff < null_diff, na.rm = TRUE)
  pval <- pmin(p_upper, p_lower) * 2
  pval[rowSums(is.finite(null_diff)) == 0L | !is.finite(observed_diff)] <- NA_real_
  results <- data.frame(
    gene = pairs$gene,
    peak = pairs$peak,
    pval = pval,
    test_stat = observed_diff,
    covar = observed_diff,
    covar_control = control_result$covar,
    covar_disease = disease_result$covar,
    cor_control = control_result$cor,
    cor_disease = disease_result$cor,
    pval_control = control_result$pval,
    pval_disease = disease_result$pval,
    padj = stats::p.adjust(pval, method = "BH"),
    stringsAsFactors = FALSE
  )
  significant <- results[
    is.finite(results$padj) & results$padj <= diff_fdr,
    , drop = FALSE
  ]
  saveRDS(
    list(control = control_raw, disease = disease_raw,
         observed_difference = observed_diff, null_difference = null_diff,
         n_permutation = n_permutation),
    file.path(output_dir, "scmultimap_diff_object.rds")
  )
  metadata_null_generation <- NA_character_
  metadata_binarize_atac <- NA
  metadata_covariates <- NA_character_
  metadata_correction <- NA
} else {
  fitted <- EpiLinkDiff(
    object = object,
    gene_peak_df = pairs,
    celltype_col = settings$celltype_column,
    family = diff_family,
    celltype = NULL,
    condition_col = condition_column,
    donor_col = donor_column,
    donor_effect = donor_effect,
    donor_robust = donor_robust,
    covariates = diff_covariates,
    null_generation = diff_null_generation,
    stat = settings$statistic,
    bin = diff_bin,
    selection = settings$selection,
    returnObject = TRUE,
    fdrcutoff = diff_fdr,
    ncore = settings$n_cores,
    theory = FALSE,
    conservative = settings$conservative,
    correction = diff_correction,
    verbose = TRUE,
    seed = settings$seed
  )
  if (is.null(fitted)) stop("EpiLinkDiff returned no fitted pairs", call. = FALSE)
  saveRDS(fitted, file.path(output_dir, "epilink_diff_object.rds"))
  results <- fitted@misc$EpiLink$res_df
  significant <- results[
    !is.na(results$significant) & results$significant,
    , drop = FALSE
  ]
  metadata_null_generation <- diff_null_generation
  metadata_binarize_atac <- diff_bin
  metadata_covariates <- paste(diff_covariates, collapse = ";")
  metadata_correction <- diff_correction
}

utils::write.table(
  results, gzfile(file.path(output_dir, "all_pairs.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  significant, gzfile(file.path(output_dir, "significant_pairs.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

metadata_run <- data.frame(
  dataset = dataset_name,
  celltype = celltype,
  method = output_method,
  condition_column = condition_column,
  donor_column = if (is.null(donor_column)) NA_character_ else donor_column,
  donor_effect = donor_effect,
  donor_robust = effective_donor_robust,
  condition_levels = paste(condition_levels, collapse = ";"),
  started = as.character(started),
  finished = as.character(Sys.time()),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  R = R.version.string,
  n_cells = length(keep_cells),
  n_candidates = nrow(pairs),
  n_tested = sum(is.finite(results$pval)),
  n_significant = nrow(significant),
  genome = if (is.null(settings$genome)) NA_character_ else as.character(settings$genome),
  fdr = diff_fdr,
  n_cores = settings$n_cores,
  seed = settings$seed,
  correction = metadata_correction,
  conservative = settings$conservative,
  binarize_atac = metadata_binarize_atac,
  null_generation = metadata_null_generation,
  covariates = metadata_covariates,
  n_permutation = if (identical(method_id, "scmultimap")) n_permutation else NA_integer_,
  stringsAsFactors = FALSE
)
utils::write.table(
  metadata_run, done_file, sep = "\t", quote = FALSE, row.names = FALSE
)
log_message(
  "Finished ", output_method, ": ", nrow(significant), "/", nrow(results),
  " significant"
)
