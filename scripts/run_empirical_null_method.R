#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
args <- setdiff(args, "--overwrite")
workers_arg <- grep("^--workers=", args, value = TRUE)
batch_arg <- grep("^--batch=", args, value = TRUE)
positional <- args[!grepl("^--(workers|batch)=", args)]
if (length(positional) != 1L || length(workers_arg) > 1L || length(batch_arg) > 1L) {
  stop(
    "Usage: Rscript scripts/run_empirical_null_method.R <method> ",
    "[--workers=10] [--batch=1..10] [--overwrite]", call. = FALSE
  )
}
method_id <- positional[[1]]
file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml is required", call. = FALSE)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)
settings <- yaml::read_yaml("config/empirical_null.yml")
registry <- empirical_null_methods()
if (!method_id %in% registry$method_id) {
  stop("Unknown method. Choose: ", paste(registry$method_id, collapse = ", "),
       call. = FALSE)
}
method_group <- empirical_null_method_group(method_id)
workers <- if (length(workers_arg)) {
  as.integer(sub("^--workers=", "", workers_arg))
} else as.integer(settings$parallel_replicates)
if (!is.finite(workers) || workers < 1L) stop("workers must be positive", call. = FALSE)

replicates <- seq_len(settings$n_replicates)
if (length(batch_arg)) {
  batch <- as.integer(sub("^--batch=", "", batch_arg))
  n_folders <- ceiling(settings$n_replicates / settings$replicates_per_folder)
  if (!is.finite(batch) || batch < 1L || batch > n_folders) {
    stop("batch must be between 1 and ", n_folders, call. = FALSE)
  }
  lower <- (batch - 1L) * settings$replicates_per_folder + 1L
  upper <- min(batch * settings$replicates_per_folder, settings$n_replicates)
  replicates <- seq.int(lower, upper)
}

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
if (!requireNamespace("Seurat", quietly = TRUE) ||
    !requireNamespace("Signac", quietly = TRUE)) {
  stop("Seurat and Signac are required", call. = FALSE)
}
if (method_id == "signac" && !requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  stop("EnsDb.Hsapiens.v86 is required for Signac", call. = FALSE)
}
if (method_id == "scent" && !requireNamespace("SCENT", quietly = TRUE)) {
  stop("SCENT is required for the SCENT empirical-null method", call. = FALSE)
}
if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)

dataset_configuration <- yaml::read_yaml("config/datasets.yml")
dataset_settings <- utils::modifyList(
  dataset_configuration$defaults,
  dataset_configuration[[settings$dataset]]
)
output_dir <- settings$output_dir
setup_dir <- file.path(output_dir, "setup")
required_setup <- c(
  "candidate_pairs.rds", "cell_names.rds", "permutation_indices.rds",
  "permutation_manifest.tsv"
)
if ("signac" %in% method_group) {
  required_setup <- c(required_setup, "signac_peak_metadata.rds")
}
missing_setup <- required_setup[!file.exists(file.path(setup_dir, required_setup))]
if (length(missing_setup)) {
  stop("Missing setup files; run scripts/prepare_empirical_null.R first",
       call. = FALSE)
}

message("Loading fixed PBMC/Tnk inputs for method ", method_id)
rna <- readRDS(dataset_settings$rna)
atac <- readRDS(dataset_settings$atac)
metadata <- readRDS(dataset_settings$metadata)
cell_names <- readRDS(file.path(setup_dir, "cell_names.rds"))
if (!identical(colnames(rna), colnames(atac))) {
  stop("RNA and ATAC cell order differs", call. = FALSE)
}
rna <- rna[, cell_names, drop = FALSE]
atac <- atac[, cell_names, drop = FALSE]
metadata <- metadata[cell_names, , drop = FALSE]
pairs <- readRDS(file.path(setup_dir, "candidate_pairs.rds"))
permutation_indices <- readRDS(file.path(setup_dir, "permutation_indices.rds"))
manifest <- utils::read.delim(file.path(setup_dir, "permutation_manifest.tsv"))
inference_seeds <- epilink_pair_seeds(
  settings$n_replicates, settings$seed + 2000L
)

annotation <- NULL
peak_metadata <- NULL
if (method_id == "signac") {
  annotation <- Signac::GetGRangesFromEnsDb(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  )
  annotation_levels <- GenomeInfoDb::seqlevels(annotation)
  GenomeInfoDb::seqlevels(annotation) <- ifelse(
    annotation_levels == "MT", "chrM", paste0("chr", annotation_levels)
  )
  peak_metadata <- readRDS(file.path(setup_dir, "signac_peak_metadata.rds"))
}

make_object <- function(atac_counts) {
  object_metadata <- metadata
  object_metadata$nCount_ATAC <- as.numeric(Matrix::colSums(atac_counts))
  object <- Seurat::CreateSeuratObject(
    counts = rna, meta.data = object_metadata, assay = "RNA"
  )
  object[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = atac_counts, sep = InferPeakSeparator(rownames(atac_counts)),
    annotation = annotation
  )
  if (!is.null(peak_metadata)) {
    object[["ATAC"]] <- Seurat::AddMetaData(
      object[["ATAC"]], metadata = peak_metadata
    )
  }
  object
}

result_path <- function(replicate, one_method) empirical_null_result_path(
  output_dir, replicate, one_method, settings$replicates_per_folder
)
group_complete <- function(replicate) all(vapply(
  method_group, function(one_method) {
    file.exists(result_path(replicate, one_method))
  }, logical(1)
))
pending <- replicates[overwrite | !vapply(replicates, group_complete, logical(1))]
if (!length(pending)) {
  message("All requested results already exist for ", method_id)
  quit(save = "no", status = 0)
}

run_one <- function(replicate) {
  paths <- vapply(method_group, function(one_method) {
    path <- result_path(replicate, one_method)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    path
  }, character(1))
  permutation <- permutation_indices[, replicate]
  atac_null <- atac[, permutation, drop = FALSE]
  colnames(atac_null) <- cell_names
  object <- make_object(atac_null)
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    fit_empirical_null_method_group(
      object, pairs, method_id = method_id, alpha = settings$alpha,
      covariates = unlist(settings$covariates, use.names = FALSE),
      seed = inference_seeds[replicate],
      n_cores = settings$per_replicate_cores,
      selection = settings$selection,
      binarize_atac = settings$binarize_atac,
      epilink_correction = settings$epilink_correction,
      epilink_conservative = settings$epilink_conservative,
      scmultimap_irls = settings$scmultimap_irls,
      scent_cores = if (is.null(settings$scent_cores)) {
        settings$per_replicate_cores
      } else settings$scent_cores
    ),
    error = function(e) {
      setNames(lapply(method_group, function(one_method) {
        empty_empirical_result(
          pairs, one_method,
          registry$method[match(one_method, registry$method_id)],
          conditionMessage(e)
        )
      }), method_group)
    }
  )
  elapsed <- proc.time()[["elapsed"]] - started
  for (one_method in names(result)) {
    one_result <- result[[one_method]]
    one_result$replicate <- as.integer(replicate)
    one_result$folder <- empirical_null_folder(
      replicate, settings$replicates_per_folder
    )
    one_result$permutation_seed <- manifest$permutation_seed[
      match(replicate, manifest$replicate)
    ]
    one_result$inference_seed <- inference_seeds[replicate]
    one_result$runtime_seconds <- elapsed
    path <- result_path(replicate, one_method)
    temporary_path <- paste0(path, ".tmp-", Sys.getpid())
    connection <- gzfile(temporary_path, open = "wt")
    utils::write.table(one_result, connection, sep = "\t", quote = FALSE,
                       row.names = FALSE)
    close(connection)
    if (file.exists(path)) unlink(path)
    if (!file.rename(temporary_path, path)) {
      stop("Could not finalize result file: ", path, call. = FALSE)
    }
  }
  data.frame(
    replicate = replicate,
    n_tested = sum(vapply(result, function(x) sum(is.finite(x$pval)), numeric(1))),
    n_errors = as.integer(any(vapply(result, function(x) {
      any(!is.na(x$error))
    }, logical(1)))),
    runtime_seconds = elapsed,
    error = {
      values <- unique(unlist(lapply(result, function(x) stats::na.omit(x$error))))
      if (length(values)) paste(values, collapse = "; ") else NA_character_
    },
    stringsAsFactors = FALSE
  )
}

safe_run_one <- function(replicate) {
  tryCatch(
    run_one(replicate),
    error = function(e) data.frame(
      replicate = replicate, n_tested = 0L, n_errors = 1L,
      runtime_seconds = NA_real_, error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  )
}

message("Running ", length(pending), " pending replicates for ", method_id,
        " in groups of at most ", workers)
status_rows <- list()
chunk_starts <- seq.int(1L, length(pending), by = workers)
for (chunk_index in seq_along(chunk_starts)) {
  from <- chunk_starts[[chunk_index]]
  chunk <- pending[from:min(from + workers - 1L, length(pending))]
  message("Round ", chunk_index, "/", length(chunk_starts),
          ": replicates ", paste(chunk, collapse = ","))
  rows <- if (length(chunk) > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      chunk, safe_run_one, mc.cores = min(workers, length(chunk)),
      mc.preschedule = FALSE, mc.set.seed = FALSE
    )
  } else lapply(chunk, safe_run_one)
  status_rows <- c(status_rows, rows)
  gc()
}
status <- do.call(rbind, status_rows)
status$method_id <- method_id
utils::write.table(
  status,
  file.path(output_dir, paste0("last_run_status_", method_id, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("Finished ", method_id, ": ", nrow(status), " replicates; errors=",
        sum(status$n_errors))
