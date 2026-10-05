#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(GenomicRanges)
  library(IRanges)
})
source("R/eqtl_validation.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else "config/eqtl_validation.yml"
only_dataset <- if (length(args) >= 2L) args[[2L]] else NULL
only_celltype <- if (length(args) >= 3L) args[[3L]] else NULL
config <- read_eqtl_validation_config(config_path)
out_root <- "results/validation/eQTL/OneK1K/eQTL2"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
tss <- load_eqtl_tss(config$reference$processed_dir)
top_fractions <- as.numeric(unlist(config$ranking$top_fractions))
recall_grid <- as.numeric(unlist(config$ranking$recall_grid))
all_curves <- list()
all_summaries <- list()
all_native <- list()

evaluate_native <- function(selected, units) {
  selected_key <- paste(selected$gene, selected$peak, sep = "\r")
  unit_key <- paste(units$gene, units$peak, sep = "\r")
  called <- unit_key %in% selected_key
  tp <- sum(called & units$positive)
  fp <- sum(called & !units$positive)
  positives <- sum(units$positive)
  base_rate <- mean(units$positive)
  data.frame(
    n_selected_units = sum(called), TP = tp, FP = fp,
    recall = if (positives) tp / positives else NA_real_,
    precision = if (tp + fp) tp / (tp + fp) else NA_real_,
    baseline_precision = base_rate,
    enrichment = if ((tp + fp) && base_rate > 0) tp / (tp + fp) / base_rate else NA_real_,
    n_units = nrow(units), n_positive = positives
  )
}

for (dataset in names(config$benchmarks)) {
  if (!is.null(only_dataset) && dataset != only_dataset) next
  for (celltype in names(config$benchmarks[[dataset]])) {
    if (!is.null(only_celltype) && celltype != only_celltype) next
    group <- config$benchmarks[[dataset]][[celltype]]
    cell_dir <- file.path("results/real", dataset, celltype)
    if (!dir.exists(cell_dir)) {
      message("SKIP missing result directory: ", dataset, "/", celltype)
      next
    }
    message("Ranked eQTL: ", dataset, "/", celltype, " [", group, "]")
    universe <- read_candidate_universe(dataset, celltype)
    eqtl <- load_eqtl_reference(group, config$reference$processed_dir)
    eqtl <- filter_distal_eqtls(eqtl, tss,
                                config$reference$min_tss_distance_bp)
    units <- make_ranked_validation_units(
      universe, eqtl, config$reference$fine_map_pip_rank
    )
    if (!nrow(units) || !any(units$positive)) {
      message("  SKIP: no assessable PIP >= 0.5 variant-gene units")
      next
    }
    unit_path <- file.path(out_root, dataset, celltype)
    dir.create(unit_path, recursive = TRUE, showWarnings = FALSE)
    fwrite(units, file.path(unit_path, "validation_units.tsv.gz"), sep = "\t")

    score_cache <- list()
    score_source <- c(poisson_glm = "epilink_poisson", nb_glm = "epilink_nb")
    for (method in names(config$methods)) {
      method_out <- file.path(unit_path, method)
      source_method <- if (method %in% names(score_source)) score_source[[method]] else method
      if (is.null(score_cache[[source_method]])) {
        score_cache[[source_method]] <- read_method_score_table(
          dataset, celltype, source_method, config
        )
      }
      score_table <- score_cache[[source_method]]
      if (is.null(score_table) && source_method != method) {
        score_table <- read_method_score_table(dataset, celltype, method, config)
      }
      for (score_type in c("primary", "evidence")) {
        if (is.null(score_table)) next
        score_column_internal <- paste0(score_type, "_score")
        scored <- score_table[, .(
          gene, peak,
          score = eqtl_rank_score(get(score_column_internal), method)
        )]
        curve <- ranked_curve(scored, units, top_fractions)
        if (!nrow(curve)) next
        curve[, `:=`(dataset = dataset, celltype = celltype,
                     reference_group = group, method = method,
                     score_type = score_type,
                     score_column = if (score_type == "primary")
                       config$methods[[method]]$primary_score else
                       config$methods[[method]]$evidence_score)]
        setcolorder(curve, c("dataset", "celltype", "reference_group", "method",
                             "score_type", "score_column"))
        dir.create(method_out, recursive = TRUE, showWarnings = FALSE)
        fwrite(curve, file.path(method_out,
                               paste0("ranked_curve_", score_type, ".tsv")),
               sep = "\t")
        summary <- summarize_ranked_curve(curve, recall_grid)
        summary <- cbind(data.frame(
          dataset = dataset, celltype = celltype, reference_group = group,
          method = method, score_type = score_type,
          score_column = curve$score_column[[1L]]
        ), summary)
        all_curves[[paste(dataset, celltype, method, score_type)]] <- curve
        all_summaries[[paste(dataset, celltype, method, score_type)]] <- summary
      }

    }

    # Native significance is evaluated separately from continuous ranking.
    # Signac has two native views (BH-corrected and uncorrected); these share
    # the same all_pairs score and therefore are not duplicated in ranking.
    for (spec_i in seq_len(nrow(eqtl1_method_specs(config)))) {
      spec <- eqtl1_method_specs(config)[spec_i]
      selected <- read_method_pairs(
        dataset, celltype, spec$source_method, config,
        significant = TRUE, score_type = "primary",
        significance_file = spec$significance_file
      )
      if (is.null(selected)) next
      # The significant file already encodes the method's two-sided/native
      # threshold. Do not discard significant negative-effect links.
      selected <- selected[is.finite(score), .(gene, peak)]
      native <- evaluate_native(selected, units)
      native <- cbind(data.frame(
        dataset = dataset, celltype = celltype, reference_group = group,
        method = spec$method,
        selection = "method_native_positive_significant_pairs",
        n_selected_pairs = nrow(selected)
      ), native)
      method_out <- file.path(unit_path, spec$method)
      dir.create(method_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(native, file.path(method_out, "native_significance.tsv"), sep = "\t")
      all_native[[paste(dataset, celltype, spec$method)]] <- native
    }
  }
}

configured_panels <- rbindlist(lapply(names(config$benchmarks), function(dataset) {
  data.table(dataset = dataset,
             celltype = names(config$benchmarks[[dataset]]))
}))

curve_files <- Sys.glob(file.path(out_root, "*", "*", "*", "ranked_curve_*.tsv"))
if (length(curve_files)) {
  existing_curves <- rbindlist(lapply(curve_files, fread), fill = TRUE)
  existing_curves <- merge(existing_curves, configured_panels,
                            by = c("dataset", "celltype"), all = FALSE)
  fwrite(existing_curves, file.path(out_root, "all_curves.tsv.gz"), sep = "\t")
  existing_summaries <- existing_curves[, {
    as.list(summarize_ranked_curve(.SD, recall_grid)[1L, ])
  }, by = .(dataset, celltype, reference_group, method, score_type, score_column)]
  fwrite(existing_summaries, file.path(out_root, "all_summaries.tsv"), sep = "\t")
}
native_files <- Sys.glob(file.path(out_root, "*", "*", "*", "native_significance.tsv"))
if (length(native_files)) {
  native <- rbindlist(lapply(native_files, fread), fill = TRUE)
  native <- merge(native, configured_panels,
                  by = c("dataset", "celltype"), all = FALSE)
  native <- native[method %in% eqtl1_method_specs(config)$method]
  fwrite(native,
         file.path(out_root, "all_native_significance.tsv"), sep = "\t")
}

message("Ranked eQTL validation complete: ", out_root)
