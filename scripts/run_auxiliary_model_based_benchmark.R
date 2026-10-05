#!/usr/bin/env Rscript

# Re-run one auxiliary method on the exact deterministic datasets from the
# existing model_based simulation.  This keeps the main simulation
# checkpoint untouched while allowing alternative Signac p-value handling and
# an independent ArchR checkpoint workflow.

raw_args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% raw_args
workers_arg <- grep("^--workers=", raw_args, value = TRUE)
archr_cores_arg <- grep("^--archr-cores=", raw_args, value = TRUE)
raw_args <- setdiff(raw_args, c("--overwrite", workers_arg, archr_cores_arg))
if (length(raw_args) != 1L || length(workers_arg) > 1L ||
    length(archr_cores_arg) > 1L) {
  stop(
    "Usage: Rscript scripts/run_auxiliary_model_based_benchmark.R ",
    "<signac_uncorrected|archr> [--workers=8] [--archr-cores=1] [--overwrite]",
    call. = FALSE
  )
}
method_id <- raw_args[[1L]]
if (!method_id %in% c("signac_uncorrected", "archr")) {
  stop("method must be signac_uncorrected or archr", call. = FALSE)
}

file_argument <- sub(
  "^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
)
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

workers <- if (length(workers_arg)) as.integer(sub("^--workers=", "", workers_arg)) else {
  if (method_id == "archr") 1L else 8L
}
archr_cores <- if (length(archr_cores_arg)) {
  as.integer(sub("^--archr-cores=", "", archr_cores_arg))
} else 1L
if (!is.finite(workers) || workers < 1L || !is.finite(archr_cores) || archr_cores < 1L) {
  stop("workers and archr-cores must be positive integers", call. = FALSE)
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("the yaml package is required", call. = FALSE)
}
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}
if (method_id == "archr") archr_requirements()

configuration <- yaml::read_yaml("config/simulation.yml")
sample_settings <- utils::modifyList(
  configuration$defaults, configuration$profiles$model_based
)
output_dir <- sample_settings$output_dir
main_metrics_file <- file.path(output_dir, "replicate_metrics.rds")
summary_file <- file.path(output_dir, "summary_metrics_cells_100_500.tsv")
if (!file.exists(main_metrics_file)) {
  stop("Missing existing simulation checkpoint: ", main_metrics_file, call. = FALSE)
}
if (!file.exists(summary_file)) {
  stop("Missing existing simulation summary: ", summary_file, call. = FALSE)
}

main_metrics <- readRDS(main_metrics_file)
scenarios <- main_metrics[
  main_metrics$method_id == "signac" &
    main_metrics$n_cells %in% unlist(sample_settings$n_cells, use.names = FALSE),
  c("dgp", "effect_size", "effect_fraction", "n_cells", "n_pairs",
    "replicate", "task_id", "simulation_seed", "inference_seed"),
  drop = FALSE
]
scenarios <- scenarios[!duplicated(scenarios$task_id), , drop = FALSE]
scenarios <- scenarios[order(scenarios$task_id), , drop = FALSE]
if (!nrow(scenarios)) stop("No Signac scenarios found in main checkpoint", call. = FALSE)

registry <- epilink_simulation_methods()
method_spec <- if (method_id == "archr") {
  registry[registry$method_id == "archr", , drop = FALSE]
} else {
  registry[registry$method_id == "signac", , drop = FALSE]
}
if (nrow(method_spec) != 1L) stop("Could not resolve method specification")
method_seed_index <- match(method_spec$method_id, registry$method_id)

checkpoint_file <- file.path(
  output_dir, paste0(method_id, "_replicate_metrics.rds")
)
checkpoint_tsv <- sub("[.]rds$", ".tsv.gz", checkpoint_file)
if (file.exists(checkpoint_file) && !isTRUE(overwrite)) {
  checkpoint <- readRDS(checkpoint_file)
} else {
  checkpoint <- NULL
}

required_checkpoint_columns <- c("task_id", "n_cells", "dgp", "effect_size",
                                 "replicate", "method_id")
if (!is.null(checkpoint) &&
    length(setdiff(required_checkpoint_columns, names(checkpoint)))) {
  stop("Existing checkpoint is missing required columns; use --overwrite",
       call. = FALSE)
}
completed_task_ids <- if (is.null(checkpoint)) integer() else {
  unique(as.integer(checkpoint$task_id))
}

save_checkpoint <- function(metrics) {
  metrics <- metrics[order(metrics$task_id), , drop = FALSE]
  temporary <- paste0(checkpoint_file, ".tmp-", Sys.getpid())
  saveRDS(metrics, temporary)
  if (file.exists(checkpoint_file)) unlink(checkpoint_file)
  if (!file.rename(temporary, checkpoint_file)) {
    stop("Could not finalize checkpoint: ", checkpoint_file, call. = FALSE)
  }
  connection <- gzfile(paste0(checkpoint_tsv), open = "wt")
  utils::write.table(metrics, connection, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  close(connection)
  invisible(metrics)
}

run_one <- function(index) {
  scenario <- scenarios[index, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    simulation <- simulate_epilink_data(
      n_cells = as.integer(scenario$n_cells),
      n_pairs = as.integer(scenario$n_pairs),
      effect_fraction = as.numeric(scenario$effect_fraction),
      effect_size = as.numeric(scenario$effect_size),
      dgp = as.character(scenario$dgp),
      atac_probability = sample_settings$atac_probability,
      intercept = sample_settings$intercept,
      nb_size = sample_settings$nb_size,
      zero_inflation = sample_settings$zero_inflation,
      seed = as.integer(scenario$simulation_seed)
    )
    method_seed <- epilink_pair_seeds(
      nrow(registry), as.integer(scenario$inference_seed)
    )[[method_seed_index]]
    if (method_id == "archr") {
      score <- run_archr_simulation(
        simulation, q = sample_settings$fdr, seed = method_seed,
        n_cores = archr_cores
      )
      score$method <- "ArchR"
      score$method_id <- "archr"
    } else {
      object <- make_simulation_seurat(simulation)
      fitted <- fit_simulated_pairs(
        simulation, object, method_spec, q = sample_settings$fdr,
        seed = method_seed, signac_pvalue_correction = FALSE
      )
      score <- score_simulation_selection(fitted$selected, simulation$truth)
      score$method <- "Signac (uncorrected)"
      score$method_id <- "signac_uncorrected"
      score$finite_sample_correction <- FALSE
      score$plus_one <- FALSE
      score$n_fitted <- fitted$n_fitted
    }
    score$finite_sample_correction <- FALSE
    score$plus_one <- FALSE
    score$runtime_seconds <- proc.time()[["elapsed"]] - started
    score$error <- NA_character_
    score$dgp <- scenario$dgp
    score$effect_size <- scenario$effect_size
    score$effect_fraction <- scenario$effect_fraction
    score$n_cells <- scenario$n_cells
    score$n_pairs <- scenario$n_pairs
    score$replicate <- scenario$replicate
    score$task_id <- scenario$task_id
    score$simulation_seed <- scenario$simulation_seed
    score$inference_seed <- scenario$inference_seed
    score
  }, error = function(e) {
    data.frame(
      discoveries = 0L, true_positives = 0L, false_positives = 0L,
      fdr = 0, power = NA_real_,
      method = if (method_id == "archr") "ArchR" else "Signac (uncorrected)",
      method_id = method_id, finite_sample_correction = FALSE,
      plus_one = FALSE, n_fitted = 0L, runtime_seconds = NA_real_,
      error = conditionMessage(e), dgp = scenario$dgp,
      effect_size = scenario$effect_size,
      effect_fraction = scenario$effect_fraction,
      n_cells = scenario$n_cells, n_pairs = scenario$n_pairs,
      replicate = scenario$replicate, task_id = scenario$task_id,
      simulation_seed = scenario$simulation_seed,
      inference_seed = scenario$inference_seed,
      stringsAsFactors = FALSE
    )
  })
  result
}

if (!length(setdiff(scenarios$task_id, completed_task_ids))) {
  message("All ", method_id, " scenarios are already checkpointed")
} else {
  combinations <- unique(scenarios[, c("n_cells", "dgp", "effect_size"), drop = FALSE])
  for (combination_index in seq_len(nrow(combinations))) {
    combination <- combinations[combination_index, , drop = FALSE]
    group <- scenarios[
      scenarios$n_cells == combination$n_cells &
        scenarios$dgp == combination$dgp &
        scenarios$effect_size == combination$effect_size,
      , drop = FALSE
    ]
    pending <- group$task_id[!group$task_id %in% completed_task_ids]
    if (!length(pending)) next
    message("[", method_id, "] combination ", combination_index, "/",
            nrow(combinations), ": n_cells=", combination$n_cells,
            ", dgp=", combination$dgp, ", effect_size=", combination$effect_size,
            " (", length(pending), " datasets)")
    results <- if (workers <= 1L || .Platform$OS.type == "windows") {
      lapply(match(pending, scenarios$task_id), run_one)
    } else {
      parallel::mclapply(match(pending, scenarios$task_id), run_one,
                         mc.cores = workers, mc.set.seed = FALSE)
    }
    new_metrics <- do.call(rbind, results)
    checkpoint <- if (is.null(checkpoint)) new_metrics else
      rbind(checkpoint, new_metrics)
    completed_task_ids <- unique(c(completed_task_ids, new_metrics$task_id))
    save_checkpoint(checkpoint)
    message("[", method_id, "] checkpointed ", length(completed_task_ids),
            "/", nrow(scenarios), " datasets")
  }
}

if (is.null(checkpoint) || length(setdiff(scenarios$task_id, checkpoint$task_id))) {
  stop("Checkpoint is incomplete; rerun to finish before updating summary", call. = FALSE)
}

summary_metrics <- read.delim(summary_file, stringsAsFactors = FALSE,
                              check.names = FALSE)
summary_metrics <- summary_metrics[
  summary_metrics$method_id != method_id, , drop = FALSE
]
summary_rows <- do.call(rbind, lapply(
  c("FDR", "Power"), function(metric_name) {
    value_column <- if (metric_name == "FDR") "fdr" else "power"
    split_rows <- split(
      checkpoint,
      interaction(checkpoint$n_cells, checkpoint$dgp,
                  checkpoint$effect_size, drop = TRUE)
    )
    do.call(rbind, lapply(split_rows, function(x) {
      values <- x[[value_column]][is.finite(x[[value_column]])]
      data.frame(
        n_cells = x$n_cells[[1]], dgp = x$dgp[[1]],
        effect_size = x$effect_size[[1]],
        method = x$method[[1]], method_id = method_id,
        metric = metric_name,
        mean = if (length(values)) mean(values) else NA_real_,
        se = if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_,
        n = length(values), stringsAsFactors = FALSE
      )
    }))
  }
))
summary_metrics <- rbind(summary_metrics, summary_rows)
summary_metrics <- summary_metrics[
  order(summary_metrics$n_cells, summary_metrics$dgp,
        summary_metrics$effect_size, summary_metrics$method_id,
        summary_metrics$metric), , drop = FALSE
]
summary_temporary <- paste0(summary_file, ".tmp-", Sys.getpid())
utils::write.table(summary_metrics, summary_temporary, sep = "\t",
                   quote = FALSE, row.names = FALSE)
if (file.exists(summary_file)) unlink(summary_file)
if (!file.rename(summary_temporary, summary_file)) {
  stop("Could not finalize summary: ", summary_file, call. = FALSE)
}
message("Updated summary with ", method_id, " rows: ", summary_file)
