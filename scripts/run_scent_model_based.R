#!/usr/bin/env Rscript

# Run SCENT only on the deterministic datasets used by the existing model-based
# benchmark. The existing combined replicate files are extended with SCENT
# rows: they are checkpointed after every replicate, while the combined summary
# is refreshed after every complete (n_cells, dgp, effect_size) group. An
# interrupted run can therefore resume without discarding completed work.

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
include_reference <- "--include-reference" %in% args
workers_arg <- grep("^--workers=", args, value = TRUE)
scent_cores_arg <- grep("^--scent-cores=", args, value = TRUE)
args <- setdiff(args, c("--overwrite", "--include-reference",
                        workers_arg, scent_cores_arg))
if (length(args) || length(workers_arg) > 1L || length(scent_cores_arg) > 1L) {
  stop(
    "Usage: Rscript scripts/run_scent_model_based.R ",
    "[--workers=8] [--scent-cores=1] [--overwrite] [--include-reference]",
    call. = FALSE
  )
}
workers_override <- if (length(workers_arg)) {
  as.integer(sub("^--workers=", "", workers_arg))
} else NULL
scent_cores_override <- if (length(scent_cores_arg)) {
  as.integer(sub("^--scent-cores=", "", scent_cores_arg))
} else NULL
if (!is.null(workers_override) &&
    (!is.finite(workers_override) || workers_override < 1L)) {
  stop("--workers must be a positive integer", call. = FALSE)
}
if (!is.null(scent_cores_override) &&
    (!is.finite(scent_cores_override) || scent_cores_override < 1L)) {
  stop("--scent-cores must be a positive integer", call. = FALSE)
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
if (!requireNamespace("SCENT", quietly = TRUE)) {
  stop(
    "SCENT is not installed. Install with devtools::install_github(",
    "'immunogenomics/SCENT')", call. = FALSE
  )
}
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}

configuration <- yaml::read_yaml("config/simulation.yml")
sample_settings <- utils::modifyList(
  configuration$defaults, configuration$profiles$model_based
)
output_dir <- sample_settings$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
main_metrics_file <- file.path(output_dir, "replicate_metrics.rds")
main_metrics_tsv <- file.path(output_dir, "replicate_metrics.tsv.gz")
if (!file.exists(main_metrics_file)) {
  stop("Existing simulation result is missing: ", main_metrics_file)
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
    paste(missing_main_columns, collapse = ", ")
  )
}
base_metrics <- main_metrics[main_metrics$method_id != "scent", , drop = FALSE]
embedded_scent_metrics <- main_metrics[
  main_metrics$method_id == "scent", , drop = FALSE
]

if (isTRUE(include_reference)) {
  profiles <- c("model_based", "trend")
  output_files <- c(
    model_based = file.path(output_dir, "scent_replicate_metrics.rds"),
    trend = file.path(output_dir, "scent_reference_500_metrics.rds")
  )
} else {
  profiles <- "model_based"
  output_files <- c(model_based = file.path(output_dir, "scent_replicate_metrics.rds"))
}

thread_env <- c(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
do.call(Sys.setenv, as.list(thread_env))

make_grid <- function(settings) {
  cell_counts <- unlist(settings$n_cells, use.names = FALSE)
  grid <- expand.grid(
    n_cells = cell_counts,
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
  grid
}

empty_scent_score <- function(message, runtime = NA_real_) {
  data.frame(
    discoveries = NA_integer_, true_positives = NA_integer_,
    false_positives = NA_integer_, fdr = NA_real_, power = NA_real_,
    method = "SCENT", method_id = "scent",
    finite_sample_correction = FALSE, plus_one = FALSE,
    n_fitted = NA_integer_, runtime_seconds = runtime,
    error = message, stringsAsFactors = FALSE
  )
}

write_scent_summary <- function(metrics, output_file) {
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
  group <- interaction(
    long$n_cells, long$dgp, long$effect_size,
    long$method_id, long$metric, drop = TRUE
  )
  summary_rows <- lapply(split(long, group), function(x) {
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
  summary_metrics <- do.call(rbind, summary_rows)
  summary_metrics <- summary_metrics[
    order(summary_metrics$n_cells, summary_metrics$dgp,
          summary_metrics$effect_size, summary_metrics$method_id,
          summary_metrics$metric), , drop = FALSE
  ]
  summary_temporary <- paste0(output_file, ".tmp-", Sys.getpid())
  utils::write.table(
    summary_metrics, summary_temporary, sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  if (file.exists(output_file)) unlink(output_file)
  if (!file.rename(summary_temporary, output_file)) {
    stop("Could not finalize SCENT summary: ", output_file, call. = FALSE)
  }
  invisible(summary_metrics)
}

write_main_metrics_checkpoint <- function(scent_metrics) {
  if (is.null(scent_metrics) || !nrow(scent_metrics)) return(invisible(NULL))
  combined <- rbind(base_metrics, scent_metrics)
  combined <- combined[
    order(combined$n_cells, combined$dgp, combined$effect_size,
          combined$replicate, combined$method_id), , drop = FALSE
  ]
  rds_temporary <- paste0(main_metrics_file, ".tmp-", Sys.getpid())
  saveRDS(combined, rds_temporary)
  if (file.exists(main_metrics_file)) unlink(main_metrics_file)
  if (!file.rename(rds_temporary, main_metrics_file)) {
    stop("Could not finalize main replicate_metrics.rds")
  }
  tsv_temporary <- paste0(main_metrics_tsv, ".tmp-", Sys.getpid())
  connection <- gzfile(tsv_temporary, open = "wt")
  utils::write.table(combined, connection, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  close(connection)
  if (file.exists(main_metrics_tsv)) unlink(main_metrics_tsv)
  if (!file.rename(tsv_temporary, main_metrics_tsv)) {
    stop("Could not finalize main replicate_metrics.tsv.gz")
  }
  invisible(combined)
}

write_combined_summary_checkpoint <- function() {
  if (!file.exists(main_metrics_file)) return(invisible(NULL))
  combined <- readRDS(main_metrics_file)
  reference_file <- sample_settings$reference_500_file
  if (!is.null(reference_file) && file.exists(reference_file)) {
    reference <- readRDS(reference_file)
    scent_reference_file <- file.path(
      output_dir, "scent_reference_500_metrics.rds"
    )
    if (file.exists(scent_reference_file)) {
      reference <- rbind(
        reference[reference$method_id != "scent", , drop = FALSE],
        readRDS(scent_reference_file)
      )
    }
    combined <- rbind(combined, reference)
  }
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
  summary_rows <- lapply(split(long, grouping), function(x) {
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
  summary_metrics <- do.call(rbind, summary_rows)
  summary_metrics <- summary_metrics[
    order(summary_metrics$n_cells, summary_metrics$dgp,
          summary_metrics$effect_size, summary_metrics$method_id,
          summary_metrics$metric), , drop = FALSE
  ]
  summary_file <- file.path(output_dir, "summary_metrics_cells_100_500.tsv")
  summary_temporary <- paste0(summary_file, ".tmp-", Sys.getpid())
  utils::write.table(summary_metrics, summary_temporary, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  if (file.exists(summary_file)) unlink(summary_file)
  if (!file.rename(summary_temporary, summary_file)) {
    stop("Could not finalize summary_metrics_cells_100_500.tsv")
  }
  invisible(summary_metrics)
}

write_combined_summary_checkpoint <- function() {
  if (!file.exists(main_metrics_file)) return(invisible(NULL))
  combined <- readRDS(main_metrics_file)
  reference_file <- sample_settings$reference_500_file
  if (!is.null(reference_file) && file.exists(reference_file)) {
    reference <- readRDS(reference_file)
    scent_reference_file <- file.path(
      output_dir, "scent_reference_500_metrics.rds"
    )
    if (file.exists(scent_reference_file)) {
      reference <- rbind(
        reference[reference$method_id != "scent", , drop = FALSE],
        readRDS(scent_reference_file)
      )
    }
    combined <- rbind(combined, reference)
  }
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
  summary_rows <- lapply(split(long, grouping), function(x) {
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
  summary_metrics <- do.call(rbind, summary_rows)
  summary_metrics <- summary_metrics[
    order(summary_metrics$n_cells, summary_metrics$dgp,
          summary_metrics$effect_size, summary_metrics$method_id,
          summary_metrics$metric), , drop = FALSE
  ]
  summary_file <- file.path(output_dir, "summary_metrics_cells_100_500.tsv")
  summary_temporary <- paste0(summary_file, ".tmp-", Sys.getpid())
  utils::write.table(summary_metrics, summary_temporary, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  if (file.exists(summary_file)) unlink(summary_file)
  if (!file.rename(summary_temporary, summary_file)) {
    stop("Could not finalize summary_metrics_cells_100_500.tsv")
  }
  invisible(summary_metrics)
}

run_profile <- function(profile_name, destination) {
  settings <- utils::modifyList(
    configuration$defaults, configuration$profiles[[profile_name]]
  )
  if (!is.null(workers_override)) settings$n_cores <- workers_override
  if (!is.null(scent_cores_override)) settings$scent_cores <- scent_cores_override
  grid <- make_grid(settings)
  tsv_destination <- sub("[.]rds$", ".tsv.gz", destination)
  if (isTRUE(overwrite)) {
    if (file.exists(destination)) unlink(destination)
    if (file.exists(tsv_destination)) unlink(tsv_destination)
    metrics <- NULL
  } else if (file.exists(destination)) {
    metrics <- tryCatch(
      readRDS(destination),
      error = function(e) stop(
        "Could not read existing SCENT checkpoint ", destination,
        "; use --overwrite to restart. Error: ", conditionMessage(e),
        call. = FALSE
      )
    )
    required_columns <- c("task_id", "n_cells", "dgp", "effect_size", "replicate")
    missing_columns <- setdiff(required_columns, names(metrics))
    if (length(missing_columns)) {
      stop(
        "Existing SCENT checkpoint is missing columns: ",
        paste(missing_columns, collapse = ", "), "; use --overwrite",
        call. = FALSE
      )
    }
  } else {
    metrics <- if (identical(profile_name, "model_based") &&
                   nrow(embedded_scent_metrics)) {
      embedded_scent_metrics
    } else NULL
  }
  completed_task_ids <- if (is.null(metrics)) integer() else {
    unique(as.integer(metrics$task_id))
  }
  message(
    "Running SCENT for profile=", profile_name, ": ", nrow(grid),
    " datasets; outer workers=", settings$n_cores,
    ", SCENT bootstrap cores=", settings$scent_cores
  )

  run_one <- function(task_id) {
    scenario <- grid[task_id, , drop = FALSE]
    started <- proc.time()[["elapsed"]]
    score <- tryCatch({
      simulation <- simulate_epilink_data(
        n_cells = scenario$n_cells, n_pairs = settings$n_pairs,
        effect_fraction = settings$effect_fraction,
        effect_size = scenario$effect_size, dgp = scenario$dgp,
        atac_probability = settings$atac_probability,
        intercept = settings$intercept, nb_size = settings$nb_size,
        zero_inflation = settings$zero_inflation,
        seed = scenario$simulation_seed
      )

      # In the full six-method run SCENT receives the seventh deterministic
      # method seed. Reusing it makes this SCENT-only run reproducible.
      scent_seed <- epilink_pair_seeds(
        7L, scenario$inference_seed
      )[7L]
      run_scent_simulation(
        simulation, q = settings$fdr, seed = scent_seed,
        n_cores = settings$scent_cores
      )
    }, error = function(e) empty_scent_score(conditionMessage(e)))

    score$finite_sample_correction <- FALSE
    score$plus_one <- FALSE
    score$dgp <- scenario$dgp
    score$effect_size <- scenario$effect_size
    score$effect_fraction <- settings$effect_fraction
    score$n_cells <- scenario$n_cells
    score$n_pairs <- settings$n_pairs
    score$replicate <- scenario$replicate
    score$task_id <- scenario$task_id
    score$simulation_seed <- scenario$simulation_seed
    score$inference_seed <- scenario$inference_seed
    score$dataset_runtime_seconds <- proc.time()[["elapsed"]] - started
    score <- score[, c(
      "discoveries", "true_positives", "false_positives", "fdr", "power",
      "method", "method_id", "finite_sample_correction", "plus_one",
      "n_fitted", "runtime_seconds", "error", "dgp", "effect_size",
      "effect_fraction", "n_cells", "n_pairs", "replicate", "task_id",
      "simulation_seed", "inference_seed", "dataset_runtime_seconds"
    ), drop = FALSE]
    score
  }

  # Submit one (n_cells, dgp, effect_size) group at a time. Each checkpoint
  # then represents a complete DGP/cell-size/effect-size summary (20 reps).
  combinations <- unique(grid[, c("n_cells", "dgp", "effect_size"), drop = FALSE])
  completed <- length(completed_task_ids)
  for (combination_index in seq_len(nrow(combinations))) {
    combination <- combinations[combination_index, , drop = FALSE]
    task_ids <- which(
      grid$n_cells == combination$n_cells &
        grid$dgp == combination$dgp &
        grid$effect_size == combination$effect_size
    )
    pending_task_ids <- setdiff(task_ids, completed_task_ids)
    if (!length(pending_task_ids)) {
      message(
        "[", profile_name, "] combination ", combination_index, "/",
        nrow(combinations), " already complete: n_cells=", combination$n_cells,
        ", dgp=", combination$dgp, ", effect_size=", combination$effect_size
      )
      next
    }
    message(
      "[", profile_name, "] combination ", combination_index, "/",
      nrow(combinations), ": n_cells=", combination$n_cells,
      ", dgp=", combination$dgp, ", effect_size=", combination$effect_size,
      " (", length(pending_task_ids), " datasets pending)"
    )
    safe_run_one <- function(task_id) {
      tryCatch(run_one(task_id), error = function(e) {
        scenario <- grid[task_id, , drop = FALSE]
        score <- empty_scent_score(conditionMessage(e))
        score$finite_sample_correction <- FALSE
        score$plus_one <- FALSE
        score$dgp <- scenario$dgp
        score$effect_size <- scenario$effect_size
        score$effect_fraction <- settings$effect_fraction
        score$n_cells <- scenario$n_cells
        score$n_pairs <- settings$n_pairs
        score$replicate <- scenario$replicate
        score$task_id <- scenario$task_id
        score$simulation_seed <- scenario$simulation_seed
        score$inference_seed <- scenario$inference_seed
        score$dataset_runtime_seconds <- NA_real_
        score[, c(
          "discoveries", "true_positives", "false_positives", "fdr", "power",
          "method", "method_id", "finite_sample_correction", "plus_one",
          "n_fitted", "runtime_seconds", "error", "dgp", "effect_size",
          "effect_fraction", "n_cells", "n_pairs", "replicate", "task_id",
          "simulation_seed", "inference_seed", "dataset_runtime_seconds"
        ), drop = FALSE]
      })
    }
    record_result <- function(result) {
      metrics <<- if (is.null(metrics)) result else rbind(metrics, result)
      metrics <<- metrics[order(metrics$task_id), , drop = FALSE]
      completed_task_ids <<- unique(c(completed_task_ids, result$task_id))
      completed <<- length(completed_task_ids)
      if (identical(profile_name, "model_based")) {
        write_main_metrics_checkpoint(metrics)
      }
      message(
        "[", profile_name, "] replicate completed: task_id=", result$task_id,
        " (", completed, "/", nrow(grid), " datasets checkpointed)"
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
        handle <- parallel::mcparallel(
          safe_run_one(task_id), mc.set.seed = FALSE
        )
        active[[as.character(handle$pid)]] <<- list(
          handle = handle, task_id = task_id
        )
        invisible(NULL)
      }
      while (length(active) < min(settings$n_cores, length(task_ids))) {
        launch_next()
      }
      while (length(active)) {
        handles <- lapply(active, function(x) x$handle)
        names(handles) <- names(active)
        collected <- parallel::mccollect(handles, wait = FALSE)
        if (!length(collected)) {
          Sys.sleep(0.2)
          next
        }
        for (pid in names(collected)) {
          result <- collected[[pid]]
          if (inherits(result, "try-error")) {
            result <- safe_run_one(active[[pid]]$task_id)
          }
          record_result(result)
          active[[pid]] <- NULL
          launch_next()
        }
      }
      invisible(NULL)
    }
    run_pending(pending_task_ids)
    completed <- length(completed_task_ids)
    # Keep the dedicated SCENT checkpoint in sync with the main combined
    # replicate file after each complete DGP/cell-size/effect-size group.
    destination_temporary <- paste0(destination, ".tmp-", Sys.getpid())
    saveRDS(metrics, destination_temporary)
    if (file.exists(destination)) unlink(destination)
    if (!file.rename(destination_temporary, destination)) {
      stop("Could not finalize SCENT checkpoint: ", destination, call. = FALSE)
    }
    destination_tsv <- sub("[.]rds$", ".tsv.gz", destination)
    destination_tsv_temporary <- paste0(destination_tsv, ".tmp-", Sys.getpid())
    connection <- gzfile(destination_tsv_temporary, open = "wt")
    utils::write.table(metrics, connection, sep = "\t", quote = FALSE,
                       row.names = FALSE)
    close(connection)
    if (file.exists(destination_tsv)) unlink(destination_tsv)
    if (!file.rename(destination_tsv_temporary, destination_tsv)) {
      stop("Could not finalize SCENT table checkpoint: ", destination_tsv,
           call. = FALSE)
    }
    other_destinations <- unname(output_files)[unname(output_files) != destination]
    other_metrics <- lapply(other_destinations, function(path) {
      if (file.exists(path)) readRDS(path) else NULL
    })
    combined_metrics <- do.call(
      rbind, c(list(metrics), other_metrics[!vapply(other_metrics, is.null, logical(1))])
    )
    write_scent_summary(
      combined_metrics,
      file.path(output_dir, "scent_summary_metrics_cells_100_500.tsv")
    )
    write_combined_summary_checkpoint()
    message(
      "[", profile_name, "] completed combination ", combination_index,
      "/", nrow(combinations), " (", completed, "/", nrow(grid),
      " datasets checkpointed)"
    )
  }
  if (is.null(metrics)) {
    metrics <- data.frame()
  }
  metrics
}

all_metadata <- list()
for (profile_name in profiles) {
  metrics <- run_profile(profile_name, output_files[[profile_name]])
  existing_destinations <- unname(output_files)[vapply(
    unname(output_files), file.exists, logical(1)
  )]
  existing_metrics <- lapply(existing_destinations, readRDS)
  if (length(existing_metrics)) {
    write_scent_summary(
      do.call(rbind, existing_metrics),
      file.path(output_dir, "scent_summary_metrics_cells_100_500.tsv")
    )
    write_combined_summary_checkpoint()
  }
  all_metadata[[profile_name]] <- data.frame(
    profile = profile_name,
    n_datasets = length(unique(metrics$task_id)),
    n_method_runs = nrow(metrics),
    n_errors = sum(!is.na(metrics$error)),
    n_cores = utils::modifyList(
      configuration$defaults, configuration$profiles[[profile_name]]
    )$n_cores,
    scent_cores = utils::modifyList(
      configuration$defaults, configuration$profiles[[profile_name]]
    )$scent_cores,
    stringsAsFactors = FALSE
  )
}
utils::write.table(
  do.call(rbind, all_metadata),
  file.path(output_dir, "scent_run_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("SCENT-only results written to ", output_dir)
