#!/usr/bin/env Rscript

# eQTL1 validation with two gene-bootstrap summaries:
#   1) pooled 2x2-table odds ratio;
#   2) SCENT-style mean of per-gene causal-variant enrichments.
#
# Both statistics are calculated from the same 1,000 gene-bootstrap draws.

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
out_root <- if (length(args) >= 4L) args[[4L]] else
  "results/validation/eQTL/OneK1K/eQTL1"

config <- read_eqtl_validation_config(config_path)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
stopifnot(config$reference$bootstrap_replicates == 1000L ||
          config$reference$bootstrap_replicates == 1000)

quantile_or_na <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, probability, names = FALSE, type = 6))
}
finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}
finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}
finite_probability_gt_one <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x > 1)
}

message("Loading 1000 Genomes background SNPs...")
background <- load_background_variants()
tss <- load_eqtl_tss(config$reference$processed_dir)

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

    message("eQTL1 pooled + SCENT mean: ", dataset, "/", celltype,
            " [", group, "]")
    universe <- read_candidate_universe(dataset, celltype)
    cache_dir <- file.path(out_root, dataset, celltype, "_cache")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    background_cache <- file.path(cache_dir, "background_peak_overlaps.rds")
    if (file.exists(background_cache)) {
      background_hits <- readRDS(background_cache)
    } else {
      message("  Computing 1000G-to-peak overlaps once")
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

      gene_enrichment <- per_gene_scent_enrichment(counts)
      causal_total <- counts$causal_annotated + counts$causal_not_annotated
      common_total <- counts$common_annotated + counts$common_not_annotated
      pooled_counts <- counts[
        is.finite(causal_total) & is.finite(common_total) &
          causal_total > 0 & common_total > 0,
        , drop = FALSE
      ]
      pooled <- odds_ratio_from_counts(pooled_counts[, c(
        "gene", "causal_annotated", "causal_not_annotated",
        "common_annotated", "common_not_annotated"
      ), drop = FALSE])
      boot <- bootstrap_scent_metrics(
        counts,
        replicates = config$reference$bootstrap_replicates,
        seed = config$reference$bootstrap_seed
      )

      result <- data.frame(
        dataset = dataset,
        celltype = celltype,
        reference_group = group,
        method = method,
        n_selected_pairs = nrow(selected),
        n_genes_all = nrow(counts),
        n_genes_bootstrap = nrow(gene_enrichment),
        n_bootstrap = nrow(boot),
        pooled_or_point = pooled$odds_ratio_ha[[1L]],
        pooled_or_fisher = pooled$odds_ratio_fisher[[1L]],
        pooled_or_pval = pooled$pval[[1L]],
        pooled_or_bootstrap_median = finite_median(boot$pooled_or),
        pooled_or_ci_low = quantile_or_na(boot$pooled_or, 0.025),
        pooled_or_ci_high = quantile_or_na(boot$pooled_or, 0.975),
        pooled_or_probability_gt_1 = finite_probability_gt_one(boot$pooled_or),
        scent_mean_point = finite_mean(gene_enrichment$scent_enrichment),
        scent_mean_bootstrap_median = finite_median(boot$scent_mean_enrichment),
        scent_mean_ci_low = quantile_or_na(boot$scent_mean_enrichment, 0.025),
        scent_mean_ci_high = quantile_or_na(boot$scent_mean_enrichment, 0.975),
        scent_mean_probability_gt_1 = finite_probability_gt_one(
          boot$scent_mean_enrichment
        ),
        stringsAsFactors = FALSE
      )

      method_out <- file.path(out_root, dataset, celltype, method)
      dir.create(method_out, recursive = TRUE, showWarnings = FALSE)
      fwrite(counts, file.path(method_out, "per_gene_counts.tsv.gz"), sep = "\t")
      fwrite(gene_enrichment,
             file.path(method_out, "per_gene_scent_enrichment.tsv.gz"),
             sep = "\t")
      fwrite(boot,
             file.path(method_out,
                       "gene_bootstrap_pooled_or_and_scent_mean.tsv.gz"),
             sep = "\t")
      fwrite(result, file.path(method_out, "summary.tsv"), sep = "\t")
      method_rows[[method]] <- result
    }

    if (length(method_rows)) {
      cell_summary <- rbindlist(method_rows, fill = TRUE)
      dir.create(file.path(out_root, dataset, celltype), recursive = TRUE,
                 showWarnings = FALSE)
      fwrite(cell_summary, file.path(out_root, dataset, celltype, "summary.tsv"),
             sep = "\t")
    }
  }
}

summary_files <- Sys.glob(file.path(out_root, "*", "*", "*", "summary.tsv"))
if (length(summary_files)) {
  summary <- rbindlist(lapply(summary_files, fread), fill = TRUE)
  configured_panels <- rbindlist(lapply(names(config$benchmarks), function(dataset) {
    data.table(dataset = dataset,
               celltype = names(config$benchmarks[[dataset]]))
  }))
  summary <- merge(summary, configured_panels,
                   by = c("dataset", "celltype"), all = FALSE)
  summary <- summary[method %in% eqtl1_method_specs(config)$method]
  summary[, pooled_or_padj := p.adjust(pooled_or_pval, method = "BH")]
  fwrite(summary, file.path(out_root, "all_results.tsv"), sep = "\t")
}
message("Completed pooled-OR and SCENT-style mean bootstrap analysis: ", out_root)
