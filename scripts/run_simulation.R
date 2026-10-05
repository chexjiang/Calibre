#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
include_scent_flag <- "--include-scent" %in% args
args <- setdiff(args, c("--overwrite", "--include-scent"))
if (length(args) != 1L) {
  stop("Usage: Rscript scripts/run_simulation.R <diagnostic|power_pilot|trend|model_based|quick|full> [--include-scent] [--overwrite]", call. = FALSE)
}
profile_name <- args[[1]]
file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("yaml", quietly = TRUE)) stop("the yaml package is required", call. = FALSE)
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)
configuration <- yaml::read_yaml("config/simulation.yml")
if (!profile_name %in% names(configuration$profiles)) {
  stop("unknown simulation profile: ", profile_name, call. = FALSE)
}
settings <- utils::modifyList(configuration$defaults, configuration$profiles[[profile_name]])
if (include_scent_flag) settings$include_scent <- TRUE
output_dir <- settings$output_dir
raw_file <- file.path(output_dir, "replicate_metrics.tsv.gz")
result_file <- file.path(output_dir, "replicate_metrics.rds")
if ((file.exists(raw_file) || file.exists(result_file)) && !overwrite) {
  stop("simulation output already exists in ", output_dir, "; use --overwrite to replace it", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
grid$inference_seed <- epilink_pair_seeds(nrow(grid), settings$seed + 1000000L)

run_started <- Sys.time()
message(
  "Running ", nrow(grid), " datasets (", settings$n_pairs, " pairs each; cells=",
  paste(cell_counts, collapse = ","), ") with ",
  settings$n_cores, " worker(s); SCENT=", settings$include_scent
)
run_task <- function(task_id) {
  scenario <- grid[task_id, ]
  started <- proc.time()[["elapsed"]]
  simulation <- simulate_epilink_data(
    n_cells = scenario$n_cells, n_pairs = settings$n_pairs,
    effect_fraction = settings$effect_fraction,
    effect_size = scenario$effect_size, dgp = scenario$dgp,
    atac_probability = settings$atac_probability, intercept = settings$intercept,
    nb_size = settings$nb_size, zero_inflation = settings$zero_inflation,
    seed = scenario$simulation_seed
  )
  result <- run_epilink_simulation(
    simulation, method_ids = unlist(settings$methods, use.names = FALSE),
    q = settings$fdr, seed = scenario$inference_seed,
    epilink_correction = settings$epilink_correction,
    epilink_conservative = settings$epilink_conservative,
    include_scent = settings$include_scent, scent_cores = settings$scent_cores
  )
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
  result
}

task_ids <- seq_len(nrow(grid))
results <- if (settings$n_cores > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(task_ids, run_task, mc.cores = settings$n_cores, mc.set.seed = FALSE)
} else {
  lapply(task_ids, run_task)
}
metrics <- do.call(rbind, results)
metrics <- metrics[order(metrics$dgp, metrics$effect_size, metrics$replicate, metrics$method_id), ]
saveRDS(metrics, result_file)
utils::write.table(metrics, gzfile(raw_file), sep = "\t", quote = FALSE, row.names = FALSE)

run_metadata <- data.frame(
  profile = profile_name, started = as.character(run_started),
  finished = as.character(Sys.time()), R = R.version.string,
  n_datasets = nrow(grid), n_method_runs = nrow(metrics),
  n_errors = sum(!is.na(metrics$error)), include_scent = settings$include_scent,
  epilink_correction = settings$epilink_correction,
  epilink_conservative = settings$epilink_conservative,
  config = "config/simulation.yml"
)
utils::write.table(run_metadata, file.path(output_dir, "run_metadata.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("Simulation finished: ", output_dir)
