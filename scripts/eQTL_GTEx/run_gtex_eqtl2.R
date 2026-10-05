#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})
source("scripts/eQTL_GTEx/gtex_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args)) args[[1L]] else "config/gtex_eqtl_validation.yml"
config <- read_yaml(config_path)
ref <- config$reference
out_root <- "results/validation/eQTL/GTEx/eQTL2"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
eqtl <- as.data.table(readRDS(file.path(
  ref$processed_dir, "GTEx_v8_CAVIAR_PIPgt0.2_49tissues_aggregated.rds")))
top_fractions <- as.numeric(unlist(config$ranking$top_fractions))
recall_grid <- as.numeric(unlist(config$ranking$recall_grid))
rank_methods <- c("epilink_poisson", "epilink_nb", "scmultimap", "signac")
specs <- gtex_method_specs()
all_curves <- list()
all_native <- list()

evaluate_native <- function(selected, units) {
  called <- paste(units$gene, units$peak, sep = "\r") %chin%
    paste(selected$gene, selected$peak, sep = "\r")
  tp <- sum(called & units$positive)
  fp <- sum(called & !units$positive)
  positives <- sum(units$positive)
  base <- mean(units$positive)
  data.table(n_selected_units = sum(called), TP = tp, FP = fp,
             recall = if (positives) tp / positives else NA_real_,
             precision = if (tp + fp) tp / (tp + fp) else NA_real_,
             baseline_precision = base,
             enrichment = if ((tp + fp) && base > 0) tp / (tp + fp) / base else NA_real_,
             n_units = nrow(units), n_positive = positives)
}

run_scope <- function(dataset, celltype = NULL) {
  scope <- if (is.null(celltype)) "dataset_pool" else "panel"
  label <- if (is.null(celltype)) dataset else paste(dataset, celltype, sep = " / ")
  message("eQTL2 ", scope, ": ", label)
  universe <- read_scope_universe(dataset, config, celltype)
  units <- make_ranked_validation_units(universe, eqtl,
                                        as.numeric(ref$fine_map_pip_rank) + .Machine$double.eps)
  if (!nrow(units) || !any(units$positive)) {
    message("  SKIP: no assessable PIP>0.5 units")
    return(invisible(NULL))
  }
  scope_dir <- if (is.null(celltype)) file.path(out_root, "dataset_pool", dataset) else
    file.path(out_root, "panel", dataset, celltype)
  dir.create(scope_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(units, file.path(scope_dir, "validation_units.tsv.gz"), sep = "\t")

  for (method in rank_methods) {
    scores <- read_scope_scores(dataset, method, config, celltype)
    if (!nrow(scores)) next
    for (score_type in c("primary", "evidence")) {
      score_col <- paste0(score_type, "_score")
      scored <- scores[, .(gene, peak, score = eqtl_rank_score(get(score_col), method))]
      curve <- ranked_curve(scored, units, top_fractions)
      if (!nrow(curve)) next
      curve[, `:=`(scope = scope, dataset = dataset,
                   celltype = if (is.null(celltype)) "POOL" else celltype,
                   method = method, score_type = score_type,
                   score_column = if (score_type == "primary")
                     config$methods[[method]]$primary_score else
                     config$methods[[method]]$evidence_score)]
      setcolorder(curve, c("scope", "dataset", "celltype", "method", "score_type"))
      method_dir <- file.path(scope_dir, method)
      dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
      fwrite(curve, file.path(method_dir, paste0("ranked_curve_", score_type, ".tsv")), sep = "\t")
      all_curves[[paste(scope, dataset, celltype, method, score_type, sep = "::")]] <<- curve
    }
  }

  for (i in seq_len(nrow(specs))) {
    spec <- specs[i]
    selected <- read_scope_selected(dataset, spec, config, celltype)
    if (!nrow(selected)) next
    native <- evaluate_native(selected, units)
    native[, `:=`(scope = scope, dataset = dataset,
                  celltype = if (is.null(celltype)) "POOL" else celltype,
                  method = spec$method, n_selected_pairs = nrow(selected))]
    setcolorder(native, c("scope", "dataset", "celltype", "method"))
    method_dir <- file.path(scope_dir, spec$method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(native, file.path(method_dir, "native_significance.tsv"), sep = "\t")
    all_native[[paste(scope, dataset, celltype, spec$method, sep = "::")]] <<- native
  }
}

panels <- gtex_panels(config)
for (i in seq_len(nrow(panels))) run_scope(panels$dataset[[i]], panels$celltype[[i]])
for (dataset in names(config$benchmarks)) run_scope(dataset, NULL)

curves <- rbindlist(all_curves, fill = TRUE)
fwrite(curves, file.path(out_root, "all_curves.tsv.gz"), sep = "\t")
summaries <- curves[, as.list(summarize_ranked_curve(.SD, recall_grid)[1L, ]),
                    by = .(scope, dataset, celltype, method, score_type, score_column)]
fwrite(summaries, file.path(out_root, "all_summaries.tsv"), sep = "\t")
fwrite(rbindlist(all_native, fill = TRUE),
       file.path(out_root, "all_native_significance.tsv"), sep = "\t")
message("GTEx eQTL2 complete: ", out_root)

