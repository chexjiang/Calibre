#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
condition_effect <- "--condition-effect" %in% args
donor_effect <- "--donor-effect" %in% args
if (condition_effect && donor_effect) {
  stop("--condition-effect and --donor-effect are mutually exclusive", call. = FALSE)
}
condition_subset <- NULL
condition_equals <- grep("^--condition=", args)
if (length(condition_equals) > 1L) {
  stop("provide at most one --condition argument", call. = FALSE)
}
if (length(condition_equals) == 1L) {
  condition_subset <- sub("^--condition=", "", args[[condition_equals]])
  args <- args[-condition_equals]
}
condition_flag <- which(args == "--condition")
if (length(condition_flag) > 1L) {
  stop("provide at most one --condition argument", call. = FALSE)
}
if (length(condition_flag) == 1L) {
  i <- condition_flag[[1L]]
  if (i == length(args)) stop("--condition requires AD or Ctrl", call. = FALSE)
  if (!is.null(condition_subset)) stop("duplicate --condition argument", call. = FALSE)
  condition_subset <- args[[i + 1L]]
  args <- args[-c(i, i + 1L)]
}
args <- setdiff(args, c("--overwrite", "--condition-effect", "--donor-effect"))
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/run_real.R <dataset> <celltype> <method> [--condition AD|Ctrl] [--donor-effect|--condition-effect] [--overwrite]", call. = FALSE)
}
dataset_name <- args[[1]]
celltype <- args[[2]]
method_id <- args[[3]]
allowed <- c("epilink_poisson", "epilink_nb", "scmultimap",
             "scmultimap_nobatch", "signac")
if (!method_id %in% allowed) stop("unknown method: ", method_id, call. = FALSE)
if (condition_effect && !method_id %in% c("epilink_poisson", "epilink_nb")) {
  stop("--condition-effect is only supported for EpiLink methods", call. = FALSE)
}
if (donor_effect && !method_id %in% c("epilink_poisson", "epilink_nb", "scmultimap")) {
  stop("--donor-effect is supported for EpiLink and scMultiMap", call. = FALSE)
}
if (!is.null(condition_subset) && condition_effect) {
  stop("--condition-effect cannot be used with --condition AD/Ctrl", call. = FALSE)
}
if (dataset_name == "Brain2" && is.null(condition_subset) &&
    method_id %in% c("epilink_poisson", "epilink_nb") &&
    !condition_effect && !donor_effect) {
  stop("Brain2 EpiLink runs require --donor-effect or --condition-effect", call. = FALSE)
}

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml is required", call. = FALSE)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)
configuration <- yaml::read_yaml("config/datasets.yml")
if (!dataset_name %in% names(configuration)) stop("unknown dataset: ", dataset_name, call. = FALSE)
settings <- utils::modifyList(configuration$defaults, configuration[[dataset_name]])

condition_column <- if (!is.null(settings$diff_condition_column)) {
  as.character(settings$diff_condition_column)
} else {
  "condition"
}
donor_column <- if (!is.null(settings$diff_donor_column)) {
  as.character(settings$diff_donor_column)
} else {
  "donor"
}
condition_dir <- if (is.null(condition_subset)) {
  "ALL"
} else if (tolower(condition_subset) %in% c("ctrl", "control")) {
  "CTRL"
} else if (tolower(condition_subset) == "ad") {
  "AD"
} else {
  stop("--condition must be AD or Ctrl", call. = FALSE)
}

# In condition-subset runs, condition is constant and is not a model term.
# The default subset model uses only log(nCount_RNA); --donor-effect adds the
# donor fixed effect (or donor bsample for scMultiMap).
if (!is.null(condition_subset)) {
  subset_covariates <- if (is.null(settings$subset_covariates)) {
    c("log(nCount_RNA)")
  } else {
    unlist(settings$subset_covariates, use.names = FALSE)
  }
  model_covariates <- subset_covariates
  if (donor_effect && method_id %in% c("epilink_poisson", "epilink_nb")) {
    model_covariates <- c(setdiff(model_covariates, condition_column), donor_column)
  }
} else {
  model_covariates <- unlist(settings$covariates, use.names = FALSE)
}
if (is.null(condition_subset) && condition_effect) {
  model_covariates <- model_covariates[model_covariates != donor_column]
  if (!condition_column %in% model_covariates) {
    model_covariates <- c(model_covariates, condition_column)
  }
} else if (is.null(condition_subset) && donor_effect) {
  model_covariates <- model_covariates[model_covariates != condition_column]
  if (!donor_column %in% model_covariates) {
    model_covariates <- c(model_covariates, donor_column)
  }
}

output_root <- file.path(settings$benchmark_output_dir, celltype)
if (dataset_name == "Brain2") {
  output_root <- file.path(output_root, condition_dir)
}
effect_suffix <- if (condition_effect) {
  "_condition_effect"
} else if (donor_effect) {
  "_donor_effect"
} else {
  ""
}
output_method <- paste0(method_id, effect_suffix)
output_dir <- file.path(output_root, output_method)
done_file <- file.path(output_dir, "run_metadata.tsv")
if (file.exists(done_file) && !overwrite) {
  stop("completed output exists: ", output_dir, "; use --overwrite to replace", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(output_dir, "run.log")
log_message <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}
started <- Sys.time()
n_significant_uncorrected <- NA_integer_
batch_used <- NULL
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
  stop(
    "RNA and ATAC contain different cell barcodes; RNA-only=",
    length(setdiff(rna_cells, atac_cells)),
    ", ATAC-only=", length(setdiff(atac_cells, rna_cells)),
    call. = FALSE
  )
}
# Match the ATAC columns to RNA by barcode. This also avoids false failures
# from harmless names/attributes on sparse-matrix dimnames.
if (!identical(rna_cells, atac_cells)) {
  atac <- atac[, rna_cells, drop = FALSE]
  atac_cells <- as.character(colnames(atac))
}
metadata <- metadata[rna_cells, , drop = FALSE]
celltype_values <- as.character(metadata[[settings$celltype_column]])
keep_mask <- !is.na(celltype_values) & celltype_values == celltype
if (!is.null(condition_subset)) {
  if (!condition_column %in% names(metadata)) {
    stop("condition column not found in metadata: ", condition_column, call. = FALSE)
  }
  condition_values <- as.character(metadata[[condition_column]])
  keep_mask <- keep_mask & !is.na(condition_values) &
    tolower(condition_values) == tolower(condition_subset)
}
keep_cells <- rownames(metadata)[keep_mask]
if (!length(keep_cells)) {
  suffix <- if (is.null(condition_subset)) "" else paste0(" and condition ", condition_subset)
  stop("no cells found for ", celltype, suffix, call. = FALSE)
}
rna <- rna[, keep_cells, drop = FALSE]
atac <- atac[, keep_cells, drop = FALSE]
metadata <- metadata[keep_cells, , drop = FALSE]

if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Signac", quietly = TRUE)) {
  stop("Seurat and Signac are required", call. = FALSE)
}
object <- Seurat::CreateSeuratObject(counts = rna, meta.data = metadata, assay = "RNA")
peak_sep <- InferPeakSeparator(rownames(atac))
genome_name <- if (is.null(settings$genome)) "hg38" else tolower(as.character(settings$genome))
if (!genome_name %in% c("hg19", "hg38")) {
  stop("unsupported genome in config: ", settings$genome,
       "; expected hg19 or hg38", call. = FALSE)
}
if (method_id == "signac") {
  if (genome_name == "hg19") {
    if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE) ||
        !requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
      stop("EnsDb.Hsapiens.v75 and BSgenome.Hsapiens.UCSC.hg19 are required for hg19 Signac runs")
    }
    annotation <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
    signac_genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else {
    if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) ||
        !requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("EnsDb.Hsapiens.v86 and BSgenome.Hsapiens.UCSC.hg38 are required for hg38 Signac runs")
    }
    annotation <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
    signac_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
    # Rename locally instead of setting seqlevelsStyle/genome, which can make
    # GenomeInfoDb contact UCSC for chromInfo and breaks reproducible offline runs.
    annotation_levels <- GenomeInfoDb::seqlevels(annotation)
    GenomeInfoDb::seqlevels(annotation) <- ifelse(
      annotation_levels == "MT", "chrM", paste0("chr", annotation_levels)
    )
  object[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = atac, sep = peak_sep, annotation = annotation
  )
} else {
  object[["ATAC"]] <- Signac::CreateChromatinAssay(counts = atac, sep = peak_sep)
}
pairs <- FilterGenePeakPairs(object, pairs_all, selection = settings$selection, verbose = TRUE)
pairs <- unique(pairs[, c("gene", "peak"), drop = FALSE])
saveRDS(pairs, file.path(output_dir, "candidate_pairs.rds"))
log_message("Prepared ", ncol(object), " cells and ", nrow(pairs), " candidate pairs")

write_results <- function(results, significant, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(results, gzfile(file.path(directory, "all_pairs.tsv.gz")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(significant, gzfile(file.path(directory, "significant_pairs.tsv.gz")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}

if (method_id %in% c("epilink_poisson", "epilink_nb")) {
  family <- if (method_id == "epilink_poisson") "poisson" else "nb"
  log_message("Running EpiLink-", family, " with ", settings$n_cores, " cores")
  fitted <- EpiLink(
    object = object, gene_peak_df = pairs, family = family,
    celltype = NULL, covariates = model_covariates,
    null_generation = settings$null_generation, stat = settings$statistic,
    bin = settings$binarize_atac,
    selection = settings$selection, returnObject = TRUE, fdrcutoff = settings$fdr,
    ncore = settings$n_cores, theory = FALSE,
    correction = settings$correction, conservative = settings$conservative,
    verbose = TRUE, seed = settings$seed
  )
  saveRDS(fitted, file.path(output_dir, "epilink_object.rds"))
  results <- fitted@misc$EpiLink$res_df
  significant <- EpiLinkFilter(
    fitted, stat = settings$statistic, fdrcutoff = settings$fdr,
    correction = settings$correction, conservative = settings$conservative,
    verbose = FALSE
  )
  write_results(results, significant, output_dir)
  glm_id <- paste0(family, "_glm", effect_suffix)
  glm_dir <- file.path(output_root, glm_id)
  glm_significant <- results[is.finite(results$padj) & results$padj <= settings$fdr, , drop = FALSE]
  write_results(results, glm_significant, glm_dir)
  utils::write.table(
    data.frame(derived_from = method_id, family = family, fdr = settings$fdr,
               elapsed_seconds = 0, n_tested = nrow(results),
               n_significant = nrow(glm_significant)),
    file.path(glm_dir, "run_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  n_tested <- nrow(results)
  n_significant <- nrow(significant)
} else if (method_id %in% c("scmultimap", "scmultimap_nobatch")) {
  use_batch <- method_id == "scmultimap"
  irls_enabled <- if (is.null(settings$scmultimap_irls)) TRUE else isTRUE(settings$scmultimap_irls)
  # A condition-subset run has no between-condition comparison.  By default
  # it therefore follows scMultiMap's no-batch analysis; --donor-effect opts
  # into donor-aware bootstrap resampling within the selected condition.
  bsample <- if (use_batch && (is.null(condition_subset) || donor_effect)) {
    settings$batch_column
  } else {
    NULL
  }
  batch_used <- bsample
  log_message("Running original scMultiMap workflow; irls=", irls_enabled,
              "; batch=", if (use_batch) {
                if (is.null(bsample)) "none" else as.character(bsample)
              } else "none (scmultimap_nobatch)")
  if (!irls_enabled) {
    # scMultiMap releases before the fix for this branch leave theta_median
    # undefined when IRLS is disabled, although the downstream WLS step only
    # needs mu and sigma_sq. Patch the namespace locally for this run.
    scmultimap_ns <- asNamespace("scMultiMap")
    scmultimap_irls_original <- get("scMultiMap_IRLS", envir = scmultimap_ns)
    scmultimap_irls_fixed <- scmultimap_irls_original
    scmultimap_body <- as.list(body(scmultimap_irls_fixed))
    irls_branch <- which(vapply(scmultimap_body, function(expr) {
      is.call(expr) && identical(expr[[1L]], as.name("if")) &&
        identical(expr[[2L]], as.name("irls"))
    }, logical(1L)))
    if (length(irls_branch) == 1L) {
      if_expr <- scmultimap_body[[irls_branch]]
      else_expr <- as.list(if_expr[[4L]])
      if_expr[[4L]] <- as.call(c(
        list(as.name("{")),
        list(quote(theta_median <- rep(NA_real_, n_cell))),
        else_expr[-1L]
      ))
      scmultimap_body[[irls_branch]] <- if_expr
      body(scmultimap_irls_fixed) <- as.call(scmultimap_body)
      if (bindingIsLocked("scMultiMap_IRLS", scmultimap_ns)) {
        unlockBinding("scMultiMap_IRLS", scmultimap_ns)
        lock_after_patch <- TRUE
      } else {
        lock_after_patch <- FALSE
      }
      assign("scMultiMap_IRLS", scmultimap_irls_fixed, envir = scmultimap_ns)
      if (lock_after_patch) lockBinding("scMultiMap_IRLS", scmultimap_ns)
    } else {
      stop("could not patch scMultiMap_IRLS for irls=FALSE", call. = FALSE)
    }
  }
  results <- scMultiMap::scMultiMap(
    object, pairs, gene_assay = "RNA", peak_assay = "ATAC",
    bsample = bsample, p.adjust.method = "BH",
    irls = irls_enabled, verbose = TRUE
  )
  significant <- results[is.finite(results$padj) & results$padj <= settings$fdr, , drop = FALSE]
  write_results(results, significant, output_dir)
  n_tested <- sum(is.finite(results$pval))
  n_significant <- nrow(significant)
} else {
  log_message("Computing peak sequence statistics")
  object[["ATAC"]] <- Signac::RegionStats(
    object[["ATAC"]], genome = signac_genome,
    sep = c("-", "-")
  )
  object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  if (requireNamespace("future", quietly = TRUE)) {
    old_future_max_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 8 * 1024^3)
    on.exit(options(future.globals.maxSize = old_future_max_size), add = TRUE)
    future::plan(future::multicore, workers = settings$n_cores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }
  log_message("Running native Signac::LinkPeaks with ", settings$n_cores, " workers")
  linkpeaks_call <- function() Signac::LinkPeaks(
    object, peak.assay = "ATAC", expression.assay = "RNA",
    genes.use = unique(pairs$gene), pvalue_cutoff = 1, score_cutoff = 0,
    verbose = TRUE
  )
  linkpeaks_result <- try(linkpeaks_call(), silent = TRUE)
  if (inherits(linkpeaks_result, "try-error")) {
    # Signac's future-multicore path can fail to serialize its large closure
    # when several LinkPeaks jobs start simultaneously. Retry with the same
    # statistics in a sequential plan; this is slower but deterministic.
    log_message("LinkPeaks parallel attempt failed; retrying sequentially")
    if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
    object <- linkpeaks_call()
  } else {
    object <- linkpeaks_result
  }
  links <- as.data.frame(Signac::Links(object[["ATAC"]]))
  # Preserve both views in one all_pairs table. Signac currently returns the
  # one-tail Gaussian probability, so pvalue_uncorrected is the raw value and
  # pvalue is the corrected two-sided value used by the legacy output.
  links$pvalue_uncorrected <- pmin(pmax(as.numeric(links$pvalue), 0), 1)
  links$pvalue <- pmin(2 * links$pvalue_uncorrected, 1)
  pair_key <- paste(pairs$gene, pairs$peak, sep = "\r")
  link_key <- paste(links$gene, links$peak, sep = "\r")
  results <- links[match(pair_key, link_key), , drop = FALSE]
  results$gene <- pairs$gene
  results$peak <- pairs$peak
  results$padj_uncorrected <- stats::p.adjust(
    results$pvalue_uncorrected, method = "BH"
  )
  results$padj <- stats::p.adjust(results$pvalue, method = "BH")
  significant <- results[is.finite(results$padj) & results$padj <= settings$fdr, , drop = FALSE]
  significant_uncorrected <- results[
    is.finite(results$padj_uncorrected) &
      results$padj_uncorrected <= settings$fdr, , drop = FALSE
  ]
  write_results(results, significant, output_dir)
  utils::write.table(
    significant, gzfile(file.path(output_dir, "significant_pairs_corrected.tsv.gz")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    significant_uncorrected,
    gzfile(file.path(output_dir, "significant_pairs_uncorrected.tsv.gz")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  saveRDS(object, file.path(output_dir, "signac_object.rds"))
  n_tested <- sum(is.finite(results$pvalue))
  n_significant <- nrow(significant)
  n_significant_uncorrected <- nrow(significant_uncorrected)
}

metadata_run <- data.frame(
  dataset = dataset_name, celltype = celltype, method = method_id,
  condition_subset = if (is.null(condition_subset)) NA_character_ else condition_dir,
  condition_column = condition_column,
  condition_effect = condition_effect,
  donor_effect = donor_effect,
  started = as.character(started), finished = as.character(Sys.time()),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  R = R.version.string, n_cells = length(keep_cells), n_candidates = nrow(pairs),
  n_tested = n_tested, n_significant = n_significant,
  n_significant_uncorrected = n_significant_uncorrected,
  genome = genome_name, fdr = settings$fdr, n_cores = settings$n_cores, seed = settings$seed,
  correction = settings$correction, conservative = settings$conservative,
  covariates = paste(model_covariates, collapse = ";"),
  batch_column = if (is.null(batch_used)) NA_character_ else as.character(batch_used)
)
utils::write.table(metadata_run, done_file, sep = "\t", quote = FALSE, row.names = FALSE)
log_message("Finished ", method_id, ": ", n_significant, "/", n_tested, " significant")

# Refresh the cell-type-aware FlowFISH table after each completed method. The
# validator accepts a partial method set; it reads the raw ESM7 reference and
# no longer depends on pooled crispr_hg38_positive/negative TSVs.
if (file.exists("scripts/run_crispr_validation.R") &&
    file.exists("reference/CRISPR/41586_2021_3446_MOESM7_ESM.txt")) {
  try(system2(
    file.path(R.home("bin"), "Rscript"),
    c("scripts/run_crispr_validation.R", dataset_name,
      "--reference=flowfish", "--correction=F"),
    stdout = FALSE, stderr = FALSE
  ), silent = TRUE)
}
