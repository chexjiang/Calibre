#!/usr/bin/env Rscript

# Run ArchR only on the deterministic model-based simulation grid.  Existing
# methods are never re-run: the ArchR rows are checkpointed separately and are
# merged into replicate_metrics.rds/tsv.gz after every completed replicate.

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
workers_arg <- grep("^--workers=", args, value = TRUE)
archr_cores_arg <- grep("^--archr-cores=", args, value = TRUE)
max_combinations_arg <- grep("^--max-combinations=", args, value = TRUE)
args <- setdiff(args, c("--overwrite", workers_arg, archr_cores_arg,
                        max_combinations_arg))
if (length(args) || length(workers_arg) > 1L || length(archr_cores_arg) > 1L ||
    length(max_combinations_arg) > 1L) {
  stop(
    "Usage: Rscript scripts/run_archr_model_based.R ",
    "[--workers=8] [--archr-cores=1] [--max-combinations=1] [--overwrite]",
    call. = FALSE
  )
}
workers_override <- if (length(workers_arg)) {
  as.integer(sub("^--workers=", "", workers_arg))
} else NULL
archr_cores_override <- if (length(archr_cores_arg)) {
  as.integer(sub("^--archr-cores=", "", archr_cores_arg))
} else NULL
max_combinations <- if (length(max_combinations_arg)) {
  as.integer(sub("^--max-combinations=", "", max_combinations_arg))
} else NULL
if (!is.null(workers_override) &&
    (!is.finite(workers_override) || workers_override < 1L)) {
  stop("--workers must be a positive integer", call. = FALSE)
}
if (!is.null(archr_cores_override) &&
    (!is.finite(archr_cores_override) || archr_cores_override < 1L)) {
  stop("--archr-cores must be a positive integer", call. = FALSE)
}
if (!is.null(max_combinations) &&
    (!is.finite(max_combinations) || max_combinations < 1L)) {
  stop("--max-combinations must be a positive integer", call. = FALSE)
}

file_argument <- sub(
  "^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
)
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("the yaml package is required", call. = FALSE)
}
if (!requireNamespace("ArchR", quietly = TRUE)) {
  stop(
    "ArchR is not installed. Install ArchR and its Bioconductor dependencies " ,
    "before running this script.", call. = FALSE
  )
}
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}

configuration <- yaml::read_yaml("config/simulation.yml")
settings <- utils::modifyList(
  configuration$defaults, configuration$profiles$model_based
)
if (!is.null(workers_override)) settings$n_cores <- workers_override
archr_cores <- if (is.null(archr_cores_override)) 1L else archr_cores_override
output_dir <- settings$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

main_metrics_file <- file.path(output_dir, "replicate_metrics.rds")
main_metrics_tsv <- file.path(output_dir, "replicate_metrics.tsv.gz")
if (!file.exists(main_metrics_file)) {
  stop("Existing simulation result is missing: ", main_metrics_file, call. = FALSE)
}
main_metrics <- readRDS(main_metrics_file)
required_main_columns <- c(
  "method_id", "method", "dgp", "effect_size", "n_cells", "replicate",
  "task_id", "fdr", "power"
)
missing_main_columns <- setdiff(required_main_columns, names(main_metrics))
if (length(missing_main_columns)) {
  stop(
    "Existing replicate_metrics.rds is missing columns: ",
    paste(missing_main_columns, collapse = ", "), call. = FALSE
  )
}
base_metrics <- main_metrics[main_metrics$method_id != "archr", , drop = FALSE]
embedded_archr_metrics <- main_metrics[
  main_metrics$method_id == "archr", , drop = FALSE
]
archr_metrics_file <- file.path(output_dir, "archr_replicate_metrics.rds")
archr_metrics_tsv <- file.path(output_dir, "archr_replicate_metrics.tsv.gz")

thread_env <- c(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
do.call(Sys.setenv, as.list(thread_env))

make_grid <- function(settings, existing_metrics) {
  grid <- expand.grid(
    n_cells = unlist(settings$n_cells, use.names = FALSE),
    dgp = unlist(settings$data_generating_models, use.names = FALSE),
    effect_size = unlist(settings$effect_sizes, use.names = FALSE),
    replicate = seq_len(settings$n_replicates),
    stringsAsFactors = FALSE
  )
  grid$task_id <- seq_len(nrow(grid))
  grid$simulation_seed <- epilink_pair_seeds(nrow(grid), settings$seed)
  grid$inference_seed <- epilink_pair_seeds(
    nrow(grid), settings$seed + 1000000L
  )
  # Prefer the seeds already recorded by the completed six/seven-method run.
  # This guarantees that ArchR sees exactly the same simulated dataset even if
  # the YAML seed is changed later.
  key <- c("n_cells", "dgp", "effect_size", "replicate")
  seed_cols <- intersect(c(key, "task_id", "simulation_seed", "inference_seed"),
                         names(existing_metrics))
  if (length(setdiff(c("simulation_seed", "inference_seed"), seed_cols)) == 0L) {
    recorded <- existing_metrics[
      !duplicated(existing_metrics[, key, drop = FALSE]), seed_cols, drop = FALSE
    ]
    grid <- merge(grid, recorded, by = key, all.x = TRUE, sort = FALSE,
                  suffixes = c("_config", ""))
    if (anyNA(grid$simulation_seed) || anyNA(grid$inference_seed)) {
      stop("Existing metrics do not cover the configured model-based grid", call. = FALSE)
    }
    if ("task_id" %in% names(grid)) {
      grid$task_id <- grid$task_id
    }
    grid$task_id_config <- NULL
  }
  grid
}

grid <- make_grid(settings, main_metrics)
empty_archr_score <- function(message, runtime = NA_real_) {
  data.frame(
    discoveries = NA_integer_, true_positives = NA_integer_,
    false_positives = NA_integer_, fdr = NA_real_, power = NA_real_,
    method = "ArchR", method_id = "archr",
    finite_sample_correction = FALSE, plus_one = FALSE,
    n_fitted = NA_integer_, runtime_seconds = runtime,
    error = message, stringsAsFactors = FALSE
  )
}

write_archr_summary <- function(metrics, output_file) {
  if (is.null(metrics) || !nrow(metrics)) return(invisible(NULL))
  long <- rbind(
    data.frame(
      n_cells = metrics$n_cells, dgp = metrics$dgp,
      effect_size = metrics$effect_size, method = metrics$method,
      method_id = metrics$method_id, metric = "FDR", value = metrics$fdr,
      stringsAsFactors = FALSE
    ),
    data.frame(
      n_cells = metrics$n_cells, dgp = metrics$dgp,
      effect_size = metrics$effect_size, method = metrics$method,
      method_id = metrics$method_id, metric = "Power", value = metrics$power,
      stringsAsFactors = FALSE
    )
  )
  grouping <- interaction(
    long$n_cells, long$dgp, long$effect_size,
    long$method_id, long$metric, drop = TRUE
  )
  rows <- lapply(split(long, grouping), function(x) {
    values <- x$value[is.finite(x$value)]
    data.frame(
      n_cells = x$n_cells[[1]], dgp = x$dgp[[1]],
      effect_size = x$effect_size[[1]], method = x$method[[1]],
      method_id = x$method_id[[1]], metric = x$metric[[1]],
      mean = if (length(values)) mean(values) else NA_real_,
      se = if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_,
      n = length(values), stringsAsFactors = FALSE
    )
  })
  summary_metrics <- do.call(rbind, rows)
  summary_metrics <- summary_metrics[
    order(summary_metrics$n_cells, summary_metrics$dgp,
          summary_metrics$effect_size, summary_metrics$method_id,
          summary_metrics$metric), , drop = FALSE
  ]
  temporary <- paste0(output_file, ".tmp-", Sys.getpid())
  utils::write.table(summary_metrics, temporary, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  if (file.exists(output_file)) unlink(output_file)
  if (!file.rename(temporary, output_file)) {
    stop("Could not finalize ArchR summary: ", output_file, call. = FALSE)
  }
  invisible(summary_metrics)
}

write_main_metrics_checkpoint <- function(archr_metrics) {
  combined <- rbind(base_metrics, archr_metrics)
  combined <- combined[
    order(combined$n_cells, combined$dgp, combined$effect_size,
          combined$replicate, combined$method_id), , drop = FALSE
  ]
  rds_temporary <- paste0(main_metrics_file, ".tmp-", Sys.getpid())
  saveRDS(combined, rds_temporary)
  if (file.exists(main_metrics_file)) unlink(main_metrics_file)
  if (!file.rename(rds_temporary, main_metrics_file)) {
    stop("Could not finalize main replicate_metrics.rds", call. = FALSE)
  }
  tsv_temporary <- paste0(main_metrics_tsv, ".tmp-", Sys.getpid())
  connection <- gzfile(tsv_temporary, open = "wt")
  utils::write.table(combined, connection, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  close(connection)
  if (file.exists(main_metrics_tsv)) unlink(main_metrics_tsv)
  if (!file.rename(tsv_temporary, main_metrics_tsv)) {
    stop("Could not finalize main replicate_metrics.tsv.gz", call. = FALSE)
  }
  invisible(combined)
}

write_combined_summary_checkpoint <- function() {
  combined <- readRDS(main_metrics_file)
  long <- rbind(
    data.frame(
      n_cells = combined$n_cells, dgp = combined$dgp,
      effect_size = combined$effect_size, method = combined$method,
      method_id = combined$method_id, metric = "FDR", value = combined$fdr,
      stringsAsFactors = FALSE
    ),
    data.frame(
      n_cells = combined$n_cells, dgp = combined$dgp,
      effect_size = combined$effect_size, method = combined$method,
      method_id = combined$method_id, metric = "Power", value = combined$power,
      stringsAsFactors = FALSE
    )
  )
  grouping <- interaction(
    long$n_cells, long$dgp, long$effect_size,
    long$method_id, long$metric, drop = TRUE
  )
  rows <- lapply(split(long, grouping), function(x) {
    values <- x$value[is.finite(x$value)]
    data.frame(
      n_cells = x$n_cells[[1]], dgp = x$dgp[[1]],
      effect_size = x$effect_size[[1]], method = x$method[[1]],
      method_id = x$method_id[[1]], metric = x$metric[[1]],
      mean = if (length(values)) mean(values) else NA_real_,
      se = if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_,
      n = length(values), stringsAsFactors = FALSE
    )
  })
  summary_metrics <- do.call(rbind, rows)
  summary_metrics <- summary_metrics[
    order(summary_metrics$n_cells, summary_metrics$dgp,
          summary_metrics$effect_size, summary_metrics$method_id,
          summary_metrics$metric), , drop = FALSE
  ]
  summary_file <- file.path(output_dir, "summary_metrics_cells_100_500.tsv")
  temporary <- paste0(summary_file, ".tmp-", Sys.getpid())
  utils::write.table(summary_metrics, temporary, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  if (file.exists(summary_file)) unlink(summary_file)
  if (!file.rename(temporary, summary_file)) {
    stop("Could not finalize combined summary: ", summary_file, call. = FALSE)
  }
  invisible(summary_metrics)
}

if (isTRUE(overwrite)) {
  unlink(c(archr_metrics_file, archr_metrics_tsv), force = TRUE)
  archr_metrics <- NULL
} else if (file.exists(archr_metrics_file)) {
  archr_metrics <- readRDS(archr_metrics_file)
} else {
  archr_metrics <- embedded_archr_metrics
}
if (!is.null(archr_metrics) && nrow(archr_metrics)) {
  required_archr_columns <- c("task_id", "n_cells", "dgp", "effect_size", "replicate")
  missing <- setdiff(required_archr_columns, names(archr_metrics))
  if (length(missing)) {
    stop("ArchR checkpoint is missing columns: ", paste(missing, collapse = ", "),
         "; use --overwrite", call. = FALSE)
  }
}
completed_task_ids <- if (is.null(archr_metrics)) integer() else {
  unique(as.integer(archr_metrics$task_id))
}
message(
  "Running ArchR for ", nrow(grid), " datasets; outer workers=", settings$n_cores,
  ", ArchR threads per worker=", archr_cores
)

run_one <- function(task_id) {
  scenario <- grid[task_id, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    simulation <- simulate_epilink_data(
      n_cells = scenario$n_cells, n_pairs = settings$n_pairs,
      effect_fraction = settings$effect_fraction,
      effect_size = scenario$effect_size, dgp = scenario$dgp,
      atac_probability = settings$atac_probability,
      intercept = settings$intercept, nb_size = settings$nb_size,
      zero_inflation = settings$zero_inflation,
      seed = scenario$simulation_seed
    )
    archr_seed <- epilink_pair_seeds(7L, scenario$inference_seed)[7L]
    run_archr_simulation(
      simulation, q = settings$fdr, seed = archr_seed,
      n_cores = archr_cores,
      work_directory = file.path(
        output_dir, "archr_work",
        sprintf("task_%05d", as.integer(scenario$task_id))
      )
    )
  }, error = function(e) empty_archr_score(conditionMessage(e)))
  result$finite_sample_correction <- FALSE
  result$plus_one <- FALSE
  result$dgp <- scenario$dgp
  result$effect_size <- scenario$effect_size
  result$effect_fraction <- settings$effect_fraction
  result$n_cells <- scenario$n_cells
  result$n_pairs <- settings$n_pairs
  result$replicate <- scenario$replicate
  result$task_id <- scenario$task_id
  result$simulation_seed <- scenario$simulation_seed
  result$inference_seed <- scenario$inference_seed
  result$dataset_runtime_seconds <- proc.time()[["elapsed"]] - started
  result[, c(
    "discoveries", "true_positives", "false_positives", "fdr", "power",
    "method", "method_id", "finite_sample_correction", "plus_one",
    "n_fitted", "runtime_seconds", "error", "dgp", "effect_size",
    "effect_fraction", "n_cells", "n_pairs", "replicate", "task_id",
    "simulation_seed", "inference_seed", "dataset_runtime_seconds"
  ), drop = FALSE]
}

safe_run_one <- function(task_id) {
  tryCatch(run_one(task_id), error = function(e) {
    scenario <- grid[task_id, , drop = FALSE]
    result <- empty_archr_score(conditionMessage(e))
    result$dgp <- scenario$dgp
    result$effect_size <- scenario$effect_size
    result$effect_fraction <- settings$effect_fraction
    result$n_cells <- scenario$n_cells
    result$n_pairs <- settings$n_pairs
    result$replicate <- scenario$replicate
    result$task_id <- scenario$task_id
    result$simulation_seed <- scenario$simulation_seed
    result$inference_seed <- scenario$inference_seed
    result$dataset_runtime_seconds <- NA_real_
    result[, c(
      "discoveries", "true_positives", "false_positives", "fdr", "power",
      "method", "method_id", "finite_sample_correction", "plus_one",
      "n_fitted", "runtime_seconds", "error", "dgp", "effect_size",
      "effect_fraction", "n_cells", "n_pairs", "replicate", "task_id",
      "simulation_seed", "inference_seed", "dataset_runtime_seconds"
    ), drop = FALSE]
  })
}

record_result <- function(result) {
  archr_metrics <<- if (is.null(archr_metrics)) result else rbind(archr_metrics, result)
  archr_metrics <<- archr_metrics[
    order(archr_metrics$task_id), , drop = FALSE
  ]
  completed_task_ids <<- unique(c(completed_task_ids, result$task_id))
  write_main_metrics_checkpoint(archr_metrics)
  temporary <- paste0(archr_metrics_file, ".tmp-", Sys.getpid())
  saveRDS(archr_metrics, temporary)
  if (file.exists(archr_metrics_file)) unlink(archr_metrics_file)
  if (!file.rename(temporary, archr_metrics_file)) {
    stop("Could not finalize ArchR checkpoint: ", archr_metrics_file, call. = FALSE)
  }
  connection <- gzfile(paste0(archr_metrics_tsv, ".tmp-", Sys.getpid()), open = "wt")
  utils::write.table(
    archr_metrics, connection, sep = "\t", quote = FALSE, row.names = FALSE
  )
  close(connection)
  tsv_temporary <- paste0(archr_metrics_tsv, ".tmp-", Sys.getpid())
  if (file.exists(archr_metrics_tsv)) unlink(archr_metrics_tsv)
  if (!file.rename(tsv_temporary, archr_metrics_tsv)) {
    stop("Could not finalize ArchR table checkpoint: ", archr_metrics_tsv, call. = FALSE)
  }
  message(
    "ArchR replicate completed: task_id=", result$task_id,
    " (", length(completed_task_ids), "/", nrow(grid), ")"
  )
}

run_pending <- function(task_ids) {
  if (settings$n_cores <= 1L || .Platform$OS.type == "windows") {
    for (task_id in task_ids) record_result(safe_run_one(task_id))
    return(invisible(NULL))
  }
  active <- list()
  next_index <- 1L
  launch_next <- function() {
    if (next_index > length(task_ids)) return(invisible(NULL))
    task_id <- task_ids[[next_index]]
    next_index <<- next_index + 1L
    handle <- parallel::mcparallel(safe_run_one(task_id), mc.set.seed = FALSE)
    active[[as.character(handle$pid)]] <<- list(handle = handle, task_id = task_id)
    invisible(NULL)
  }
  while (length(active) < min(settings$n_cores, length(task_ids))) launch_next()
  while (length(active)) {
    handles <- lapply(active, `[[`, "handle")
    names(handles) <- names(active)
    collected <- parallel::mccollect(handles, wait = FALSE)
    if (!length(collected)) {
      Sys.sleep(0.2)
      next
    }
    for (pid in names(collected)) {
      result <- collected[[pid]]
      if (inherits(result, "try-error")) result <- safe_run_one(active[[pid]]$task_id)
      record_result(result)
      active[[pid]] <- NULL
      launch_next()
    }
  }
  invisible(NULL)
}

combinations <- unique(grid[, c("n_cells", "dgp", "effect_size"), drop = FALSE])
if (!is.null(max_combinations)) {
  combinations <- head(combinations, max_combinations)
  message("Limiting this run to the first ", nrow(combinations),
          " combination(s); rerun without --max-combinations to continue.")
}
for (combination_index in seq_len(nrow(combinations))) {
  combination <- combinations[combination_index, , drop = FALSE]
  task_ids <- which(
    grid$n_cells == combination$n_cells & grid$dgp == combination$dgp &
      grid$effect_size == combination$effect_size
  )
  pending <- setdiff(task_ids, completed_task_ids)
  if (!length(pending)) {
    message(
      "[model_based] combination ", combination_index, "/", nrow(combinations),
      " already complete: n_cells=", combination$n_cells,
      ", dgp=", combination$dgp, ", effect_size=", combination$effect_size
    )
    next
  }
  message(
    "[model_based] combination ", combination_index, "/", nrow(combinations),
    ": n_cells=", combination$n_cells, ", dgp=", combination$dgp,
    ", effect_size=", combination$effect_size,
    " (", length(pending), " datasets pending)"
  )
  run_pending(pending)
  write_archr_summary(
    archr_metrics, file.path(output_dir, "archr_summary_metrics_cells_100_500.tsv")
  )
  write_combined_summary_checkpoint()
}

write_archr_summary(
  archr_metrics, file.path(output_dir, "archr_summary_metrics_cells_100_500.tsv")
)
write_combined_summary_checkpoint()
utils::write.table(
  data.frame(
    method = "ArchR", method_id = "archr", profile = "model_based",
    n_datasets = length(unique(archr_metrics$task_id)),
    n_method_runs = nrow(archr_metrics),
    n_errors = sum(!is.na(archr_metrics$error)),
    outer_workers = settings$n_cores, archr_cores = archr_cores,
    ArchR_version = as.character(utils::packageVersion("ArchR")),
    cor_cutoff = 0.45, max_dist = 250000L,
    finished = as.character(Sys.time()), stringsAsFactors = FALSE
  ),
  file.path(output_dir, "archr_run_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("ArchR results written to ", output_dir)
