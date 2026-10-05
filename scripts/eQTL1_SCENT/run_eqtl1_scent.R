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
out_root <- "results/validation/eQTL/SCENT"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

message("Loading 1000 Genomes background SNPs...")
background <- load_background_variants()
tss <- load_eqtl_tss(config$reference$processed_dir)
all_summary <- list()

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
    message("SCENT-style: ", dataset, "/", celltype, " [", group, "]")
    universe <- read_candidate_universe(dataset, celltype)
    cache_dir <- file.path(out_root, dataset, celltype, "_cache")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    background_cache <- file.path(cache_dir, "background_peak_overlaps.rds")
    if (file.exists(background_cache)) {
      message("  Loading cached 1000G-to-peak overlaps")
      background_hits <- readRDS(background_cache)
    } else {
      message("  Computing 1000G-to-peak overlaps once for all methods")
      background_hits <- background_peak_overlap(universe, background)
      saveRDS(background_hits, background_cache, compress = FALSE)
    }
    eqtl <- load_eqtl_reference(group, config$reference$processed_dir)
    eqtl <- filter_distal_eqtls(eqtl, tss,
                                config$reference$min_tss_distance_bp)
    cis_cache <- file.path(cache_dir, "background_cis_counts.rds")
    if (file.exists(cis_cache)) {
      background_cis_counts <- readRDS(cis_cache)
    } else {
      benchmark_genes <- intersect(unique(universe$gene), unique(eqtl$gene))
      benchmark_tss <- tss[tss$gene %in% benchmark_genes, , drop = FALSE]
      benchmark_tss <- benchmark_tss[!duplicated(benchmark_tss$gene), , drop = FALSE]
      background_cis_counts <- count_positions_in_windows(
        background, benchmark_tss, config$reference$cis_window_bp
      )
      names(background_cis_counts) <- benchmark_tss$gene
      saveRDS(background_cis_counts, cis_cache, compress = FALSE)
    }
    method_rows <- list()
    method_specs <- eqtl1_method_specs(config)
    for (spec_i in seq_len(nrow(method_specs))) {
      spec <- method_specs[spec_i]
      method <- spec$method
      selected <- read_method_pairs(
        dataset, celltype, spec$source_method, config,
        significant = TRUE,
        significance_file = spec$significance_file
      )
      if (is.null(selected)) {
        message("  SKIP ", method, ": no significant pairs")
        next
      }
      # Native significance is two-sided; retain significant links regardless
      # of the sign of beta/zval/correlation.
      selected <- selected[is.finite(score), .(gene, peak)]
      counts <- scent_gene_counts(
        universe = universe, selected = selected, eqtl = eqtl,
        background = background, tss = tss,
        pip_threshold = config$reference$fine_map_pip_scent,
        window_bp = config$reference$cis_window_bp,
        background_peak_hits = background_hits,
        background_cis_counts = background_cis_counts
      )
      pooled <- odds_ratio_from_counts(counts)
      boot <- bootstrap_scent_odds(
        counts,
        replicates = config$reference$bootstrap_replicates,
        seed = config$reference$bootstrap_seed
      )
      result <- cbind(
        data.frame(dataset = dataset, celltype = celltype,
                   reference_group = group, method = method,
                   n_selected_pairs = nrow(selected),
                   n_bootstrap = length(boot)),
        pooled,
        data.frame(
          bootstrap_median_or = stats::median(boot, na.rm = TRUE),
          bootstrap_ci_low = unname(stats::quantile(boot, 0.025, na.rm = TRUE)),
          bootstrap_ci_high = unname(stats::quantile(boot, 0.975, na.rm = TRUE)),
          bootstrap_probability_or_gt_1 = mean(boot > 1, na.rm = TRUE),
          n_genes = nrow(counts)
        )
      )
      method_out <- file.path(out_root, dataset, celltype, method)
      dir.create(method_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(counts, file.path(method_out, "per_gene_counts.tsv.gz"), sep = "\t")
      fwrite(data.table(bootstrap_replicate = seq_along(boot), odds_ratio = boot),
             file.path(method_out, "gene_bootstrap_odds_ratios.tsv.gz"), sep = "\t")
      fwrite(result, file.path(method_out, "summary.tsv"), sep = "\t")
      method_rows[[method]] <- result
    }
    if (length(method_rows)) {
      cell_summary <- rbindlist(method_rows, fill = TRUE)
      dir.create(file.path(out_root, dataset, celltype), recursive = TRUE,
                 showWarnings = FALSE)
      fwrite(cell_summary, file.path(out_root, dataset, celltype, "summary.tsv"),
             sep = "\t")
      all_summary[[paste(dataset, celltype)]] <- cell_summary
    }
  }
}

summary_files <- Sys.glob(file.path(out_root, "*", "*", "*", "summary.tsv"))
summary_files <- summary_files[basename(dirname(summary_files)) != "_cache"]
if (length(summary_files)) {
  summary <- rbindlist(lapply(summary_files, fread), fill = TRUE)
  configured_panels <- rbindlist(lapply(names(config$benchmarks), function(dataset) {
    data.table(dataset = dataset,
               celltype = names(config$benchmarks[[dataset]]))
  }))
  summary <- merge(summary, configured_panels,
                   by = c("dataset", "celltype"), all = FALSE)
  summary <- summary[method %in% eqtl1_method_specs(config)$method]
  summary[, padj := p.adjust(pval, method = "BH")]
  fwrite(summary, file.path(out_root, "all_results.tsv"), sep = "\t")
}
message("SCENT-style validation complete: ", out_root)
