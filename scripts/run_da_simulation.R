#!/usr/bin/env Rscript

# Differential-association model-based simulation.
#
# Methods are run independently so that each method has its own resumable
# checkpoint. After every worker chunk, the aggregate replicate table and its
# summary are rewritten. This makes partial results visible while a long
# scMultiMap run is still in progress.

raw_args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% raw_args
workers_arg <- grep("^--workers=", raw_args, value = TRUE)
rerun_dgp_arg <- grep("^--rerun-dgp=", raw_args, value = TRUE)
raw_args <- setdiff(raw_args, c("--overwrite", workers_arg, rerun_dgp_arg))
if (length(raw_args) || length(workers_arg) > 1L ||
    length(rerun_dgp_arg) > 1L) {
  stop(
    "Usage: Rscript scripts/run_da_simulation.R [--workers=16] ",
    "[--overwrite] [--rerun-dgp=zinb]",
    call. = FALSE
  )
}

file_argument <- sub(
  "^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
)
project_root <- normalizePath(file.path(dirname(file_argument), ".."),
                              mustWork = TRUE)
setwd(project_root)
if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml is required")
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}

settings <- yaml::read_yaml("config/simulation_da.yml")
rerun_dgps <- if (length(rerun_dgp_arg)) {
  unique(trimws(strsplit(sub("^--rerun-dgp=", "", rerun_dgp_arg),
                         ",", fixed = TRUE)[[1L]]))
} else character()
available_dgps <- as.character(unlist(settings$data_generating_models,
                                       use.names = FALSE))
if (length(rerun_dgps) &&
    any(!nzchar(rerun_dgps) | !rerun_dgps %in% available_dgps)) {
  stop("--rerun-dgp must name DGPs in config: ",
       paste(available_dgps, collapse = ", "), call. = FALSE)
}
workers <- if (length(workers_arg)) {
  as.integer(sub("^--workers=", "", workers_arg))
} else as.integer(settings$n_cores)
if (!is.finite(workers) || workers < 1L) {
  stop("workers must be a positive integer", call. = FALSE)
}

output_dir <- settings$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
registry <- da_simulation_methods()
method_ids <- registry$method_id

grid <- expand.grid(
  dgp = as.character(unlist(settings$data_generating_models,
                             use.names = FALSE)),
  n_cells = as.integer(unlist(settings$n_cells, use.names = FALSE)),
  effect_size = as.numeric(unlist(settings$effect_sizes, use.names = FALSE)),
  replicate = seq_len(as.integer(settings$n_replicates)),
  stringsAsFactors = FALSE
)
if (any(grid$n_cells %% settings$n_donors != 0L)) {
  stop("Each n_cells value must be divisible by n_donors", call. = FALSE)
}
grid$task_id <- seq_len(nrow(grid))
grid$simulation_seed <- epilink_pair_seeds(nrow(grid), settings$seed)
grid$inference_seed <- epilink_pair_seeds(
  nrow(grid), as.integer(settings$seed) + 1000000L
)

method_checkpoint <- function(method_id) {
  file.path(output_dir, paste0(method_id, "_replicate_metrics.rds"))
}
method_checkpoint_tsv <- function(method_id) {
  file.path(output_dir, paste0(method_id, "_replicate_metrics.tsv.gz"))
}
aggregate_rds <- file.path(output_dir, "replicate_metrics.rds")
aggregate_tsv <- file.path(output_dir, "replicate_metrics.tsv.gz")
summary_file <- file.path(output_dir, "summary_metrics.tsv")
metadata_file <- file.path(output_dir, "run_metadata.tsv")

if (overwrite) {
  unlink(c(
    aggregate_rds, aggregate_tsv, summary_file, metadata_file,
    unlist(lapply(method_ids, method_checkpoint)),
    unlist(lapply(method_ids, method_checkpoint_tsv))
  ), force = TRUE)
}

write_gz_table <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  connection <- gzfile(temporary, open = "wt")
  utils::write.table(x, connection, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  close(connection)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) {
    stop("Could not finalize ", path, call. = FALSE)
  }
}

write_atomic_rds <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, temporary)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) {
    stop("Could not finalize ", path, call. = FALSE)
  }
}

read_checkpoint <- function(method_id) {
  path <- method_checkpoint(method_id)
  if (!file.exists(path)) return(NULL)
  x <- readRDS(path)
  required <- c("task_id", "dgp", "n_cells", "effect_size", "replicate",
                "method_id", "n_pairs")
  if (length(setdiff(required, names(x)))) {
    stop("Checkpoint for ", method_id,
         " is incompatible; rerun with --overwrite", call. = FALSE)
  }
  expected <- grid[match(x$task_id, grid$task_id), , drop = FALSE]
  if (anyNA(expected$task_id) || any(x$n_pairs != settings$n_pairs) ||
      any(x$dgp != expected$dgp) || any(x$n_cells != expected$n_cells) ||
      any(x$effect_size != expected$effect_size) ||
      any(x$replicate != expected$replicate)) {
    stop("Checkpoint for ", method_id,
         " was produced with a different DA grid; use --overwrite",
         call. = FALSE)
  }
  x
}

all_checkpoint_metrics <- function() {
  pieces <- lapply(method_ids, read_checkpoint)
  pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
  if (!length(pieces)) return(NULL)
  out <- do.call(rbind, pieces)
  out[order(out$task_id, out$method_id), , drop = FALSE]
}

refresh_aggregate <- function() {
  all_metrics <- all_checkpoint_metrics()
  if (is.null(all_metrics)) return(invisible(NULL))
  write_atomic_rds(all_metrics, aggregate_rds)
  write_gz_table(all_metrics, aggregate_tsv)
  summary <- summarize_da_metrics(all_metrics)
  utils::write.table(summary, summary_file, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  invisible(all_metrics)
}

write_metadata <- function(status = "running") {
  utils::write.table(
    data.frame(
      profile = "da_poisson",
      status = status,
      dgp = paste(settings$data_generating_models, collapse = ";"),
      n_donors = settings$n_donors,
      donors_per_condition = paste(
        settings$n_donors_control,
        settings$n_donors - settings$n_donors_control, sep = ";"
      ),
      n_cells = paste(settings$n_cells, collapse = ";"),
      n_pairs = settings$n_pairs,
      n_replicates = settings$n_replicates,
      effect_sizes = paste(settings$effect_sizes, collapse = ";"),
      effect_fraction = settings$effect_fraction,
      fdr = settings$fdr,
      scmultimap_permutations = settings$scmultimap_permutations,
      workers = workers,
      finished = as.character(Sys.time()),
      R = R.version.string,
      stringsAsFactors = FALSE
    ), metadata_file, sep = "\t", quote = FALSE, row.names = FALSE
  )
}

method_label <- function(method_id) {
  registry$method[match(method_id, registry$method_id)]
}

fit_one_method <- function(method_id, simulation, object, alpha, seed) {
  if (method_id == "epilink_poisson") {
    return(fit_da_epilink(
      simulation, object, family = "poisson", alpha = alpha, seed = seed,
      n_cores = settings$pair_cores,
      null_generation = settings$null_generation,
      correction = settings$epilink_correction,
      conservative = settings$epilink_conservative,
      donor_robust = settings$donor_robust
    ))
  }
  if (method_id == "epilink_nb") {
    return(fit_da_epilink(
      simulation, object, family = "nb", alpha = alpha, seed = seed,
      n_cores = settings$pair_cores,
      null_generation = settings$null_generation,
      correction = settings$epilink_correction,
      conservative = settings$epilink_conservative,
      donor_robust = settings$donor_robust
    ))
  }
  if (method_id == "glm_poisson") {
    return(fit_da_glm_theory(
      simulation, object, family = "poisson", alpha = alpha,
      n_cores = settings$pair_cores
    ))
  }
  if (method_id == "glm_nb") {
    return(fit_da_glm_theory(
      simulation, object, family = "nb", alpha = alpha,
      n_cores = settings$pair_cores
    ))
  }
  fit_da_scmultimap(
    simulation, object, alpha = alpha, seed = seed,
    n_permutations = settings$scmultimap_permutations,
    irls = settings$scmultimap_irls,
    permutation_cores = settings$pair_cores
  )
}

run_one <- function(method_id, scenario) {
  started <- proc.time()[["elapsed"]]
  label <- method_label(method_id)
  tryCatch({
    simulation <- simulate_da_data(
      n_donors = settings$n_donors,
      cells_per_donor = as.integer(scenario$n_cells / settings$n_donors),
      n_donors_control = settings$n_donors_control,
      n_pairs = settings$n_pairs,
      effect_fraction = settings$effect_fraction,
      effect_size = scenario$effect_size,
      atac_probability = settings$atac_probability,
      intercept = settings$intercept,
      atac_effect = settings$atac_effect,
      condition_effect = settings$condition_effect,
      dgp = scenario$dgp,
      nb_size = settings$nb_size,
      zero_inflation = settings$zero_inflation,
      seed = scenario$simulation_seed
    )
    object <- if (method_id %in% c(
      "epilink_poisson", "epilink_nb", "glm_poisson", "glm_nb",
      "scmultimap"
    )) {
      make_da_seurat(simulation)
    } else NULL
    method_seed <- epilink_pair_seeds(
      nrow(registry), scenario$inference_seed
    )[[match(method_id, registry$method_id)]]
    fitted <- fit_one_method(
      method_id, simulation, object, alpha = settings$fdr,
      seed = method_seed
    )
    score <- score_da_selection(fitted$selected, simulation$truth)
    score$method <- label
    score$method_id <- method_id
    score$n_fitted <- fitted$n_fitted
    score$runtime_seconds <- proc.time()[["elapsed"]] - started
    score$error <- NA_character_
    score$dgp <- scenario$dgp
    score$effect_size <- scenario$effect_size
    score$effect_fraction <- settings$effect_fraction
    score$n_donors <- settings$n_donors
    score$cells_per_donor <- as.integer(scenario$n_cells / settings$n_donors)
    score$n_cells <- simulation$n_cells
    score$n_pairs <- settings$n_pairs
    score$replicate <- scenario$replicate
    score$task_id <- scenario$task_id
    score$simulation_seed <- scenario$simulation_seed
    score$inference_seed <- scenario$inference_seed
    score
  }, error = function(e) {
    data.frame(
      discoveries = 0L, true_positives = 0L, false_positives = 0L,
      fdr = 0, power = NA_real_, method = label, method_id = method_id,
      n_fitted = 0L, runtime_seconds = NA_real_, error = conditionMessage(e),
      dgp = scenario$dgp, effect_size = scenario$effect_size,
      effect_fraction = settings$effect_fraction,
      n_donors = settings$n_donors,
      cells_per_donor = as.integer(scenario$n_cells / settings$n_donors),
      n_cells = scenario$n_cells,
      n_pairs = settings$n_pairs, replicate = scenario$replicate,
      task_id = scenario$task_id, simulation_seed = scenario$simulation_seed,
      inference_seed = scenario$inference_seed,
      stringsAsFactors = FALSE
    )
  })
}

write_metadata("running")
message(
  "DA simulation: ", nrow(grid), " datasets x ", length(method_ids),
  " methods; workers=", workers, "; DGPs=",
  paste(settings$data_generating_models, collapse = ","),
  "; n_cells=", paste(settings$n_cells, collapse = ","),
  "; n_pairs=", settings$n_pairs,
  "; scMultiMap permutations=", settings$scmultimap_permutations
)

for (method_id in method_ids) {
  checkpoint <- read_checkpoint(method_id)
  if (length(rerun_dgps) && !is.null(checkpoint)) {
    remove_rows <- checkpoint$dgp %in% rerun_dgps
    if (any(remove_rows)) {
      checkpoint <- checkpoint[!remove_rows, , drop = FALSE]
      write_atomic_rds(checkpoint, method_checkpoint(method_id))
      write_gz_table(checkpoint, method_checkpoint_tsv(method_id))
      refresh_aggregate()
      message("[", method_id, "] removed ", sum(remove_rows),
              " checkpointed datasets for rerun: ",
              paste(rerun_dgps, collapse = ","))
    }
  }
  completed <- if (is.null(checkpoint)) integer() else
    unique(as.integer(checkpoint$task_id))
  pending <- setdiff(grid$task_id, completed)
  if (length(rerun_dgps)) {
    pending <- pending[grid$dgp[pending] %in% rerun_dgps]
  }
  if (!length(pending)) {
    message("[", method_id, "] already complete")
    refresh_aggregate()
    next
  }
  message("[", method_id, "] pending datasets: ", length(pending))
  chunk_starts <- seq.int(1L, length(pending), by = workers)
  for (chunk_index in seq_along(chunk_starts)) {
    from <- chunk_starts[[chunk_index]]
    task_ids <- pending[from:min(from + workers - 1L, length(pending))]
    message(
      "[", method_id, "] chunk ", chunk_index, "/", length(chunk_starts),
      ": task_id=", paste(task_ids, collapse = ",")
    )
    rows <- if (length(task_ids) > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(
        task_ids,
        function(task_id) run_one(
          method_id, grid[grid$task_id == task_id, , drop = FALSE]
        ),
        mc.cores = min(workers, length(task_ids)), mc.set.seed = FALSE
      )
    } else lapply(
      task_ids,
      function(task_id) run_one(
        method_id, grid[grid$task_id == task_id, , drop = FALSE]
      )
    )
    new_metrics <- do.call(rbind, rows)
    checkpoint <- if (is.null(checkpoint)) new_metrics else
      rbind(checkpoint, new_metrics)
    checkpoint <- checkpoint[order(checkpoint$task_id), , drop = FALSE]
    write_atomic_rds(checkpoint, method_checkpoint(method_id))
    write_gz_table(checkpoint, method_checkpoint_tsv(method_id))
    refresh_aggregate()
    completed <- unique(checkpoint$task_id)
    message(
      "[", method_id, "] checkpointed ", length(completed), "/",
      nrow(grid), " datasets; summary refreshed"
    )
  }
}

refresh_aggregate()
write_metadata("complete")
message("DA simulation finished: ", output_dir)
