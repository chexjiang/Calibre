#!/usr/bin/env Rscript

# Dataset-level eQTL validation plots. eQTL1 pools gene x cell-type counts;
# eQTL2 pools TP/selected/positive counts at each top fraction.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})
source("R/eqtl_validation.R")

config <- read_eqtl_validation_config("config/eqtl_validation.yml")
oneK1K_root <- "results/validation/eQTL/OneK1K"
eqtl1_plot_dir <- file.path(oneK1K_root, "eQTL1")
eqtl2_plot_dir <- file.path(oneK1K_root, "eQTL2")
eqtl1_root <- "results/validation/eQTL/OneK1K/eQTL1"
eqtl2_root <- "results/validation/eQTL/OneK1K/eQTL2"
dir.create(eqtl1_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(eqtl2_plot_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
  epilink_poisson = "EpiLink-Poisson", epilink_nb = "EpiLink-NB",
  scmultimap = "scMultiMap", signac = "Signac (corrected)",
  signac_corrected = "Signac (corrected)",
  signac_uncorrected = "Signac (uncorrected)"
)
plot_methods <- c("epilink_poisson", "epilink_nb", "scmultimap",
                  "signac_corrected", "signac_uncorrected")
rank_methods <- c("epilink_poisson", "epilink_nb", "scmultimap", "signac")
method_levels <- c("EpiLink-Poisson", "EpiLink-NB", "scMultiMap",
                   "Signac (corrected)", "Signac (uncorrected)")
method_colors <- c(
  "EpiLink-Poisson" = "#E47159", "EpiLink-NB" = "#C87878",
  "scMultiMap" = "#edd089", "Signac (corrected)" = "#54b2a2",
  "Signac (uncorrected)" = "#78A878"
)

panels <- rbindlist(lapply(names(config$benchmarks), function(dataset) {
  data.table(dataset = dataset,
             celltype = names(config$benchmarks[[dataset]]))
}))

finite_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, prob, names = FALSE, type = 6))
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

# eQTL1: pool the gene-level tables across cell types and recalculate both
# statistics. Gene x cell type is kept as the bootstrap unit.
eqtl1_rows <- list()
for (dataset in names(config$benchmarks)) {
  for (method in c("epilink_poisson", "epilink_nb", "scmultimap",
                   "signac_corrected", "signac_uncorrected")) {
    count_list <- list()
    for (celltype in names(config$benchmarks[[dataset]])) {
      path <- file.path(eqtl1_root, dataset, celltype, method,
                        "per_gene_counts.tsv.gz")
      if (!file.exists(path)) next
      counts <- fread(path, showProgress = FALSE)
      if (!nrow(counts)) next
      counts[, gene := paste(celltype, gene, sep = "::")]
      count_list[[celltype]] <- counts[, .(
        gene, causal_annotated, causal_not_annotated,
        common_annotated, common_not_annotated
      )]
    }
    if (!length(count_list)) next
    counts <- rbindlist(count_list, fill = TRUE)
    causal_total <- counts$causal_annotated + counts$causal_not_annotated
    common_total <- counts$common_annotated + counts$common_not_annotated
    pooled_counts <- counts[
      is.finite(causal_total) & is.finite(common_total) &
        causal_total > 0 & common_total > 0
    ]
    if (!nrow(pooled_counts)) next
    pooled <- odds_ratio_from_counts(pooled_counts)
    boot <- bootstrap_scent_metrics(
      counts, replicates = as.integer(config$reference$bootstrap_replicates),
      seed = as.integer(config$reference$bootstrap_seed)
    )
    gene_enrichment <- per_gene_scent_enrichment(counts)
    eqtl1_rows[[paste(dataset, method, sep = "::")]] <- data.table(
      dataset = dataset, method = method, n_celltypes = length(count_list),
      n_gene_celltype_units = nrow(counts),
      pooled_or_point = pooled$odds_ratio_ha[[1L]],
      pooled_or_bootstrap_median = finite_median(boot$pooled_or),
      pooled_or_ci_low = finite_quantile(boot$pooled_or, 0.025),
      pooled_or_ci_high = finite_quantile(boot$pooled_or, 0.975),
      scent_mean_point = finite_mean(gene_enrichment$scent_enrichment),
      scent_mean_bootstrap_median = finite_median(boot$scent_mean_enrichment),
      scent_mean_ci_low = finite_quantile(boot$scent_mean_enrichment, 0.025),
      scent_mean_ci_high = finite_quantile(boot$scent_mean_enrichment, 0.975)
    )
  }
}

if (length(eqtl1_rows)) {
  eqtl1 <- rbindlist(eqtl1_rows, fill = TRUE)
  eqtl1[, method_label := factor(unname(method_labels[method]),
                                 levels = method_levels)]

  make_eqtl1_forest <- function(dat, estimate, low, high, title, xlab,
                                filename, x_transform = c("log10", "pseudo_log")) {
    x_transform <- match.arg(x_transform)
    dat <- dat[is.finite(get(estimate)) & is.finite(get(low)) &
                 is.finite(get(high))]
    if (x_transform == "log10") {
      dat <- dat[get(estimate) > 0 & get(low) > 0 & get(high) > 0]
    }
    dat[, `:=`(estimate_plot = get(estimate), low_plot = get(low),
               high_plot = get(high))]
    x_scale <- if (x_transform == "log10") {
      scale_x_log10()
    } else {
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10))
    }
    p <- ggplot(dat, aes(y = method_label, colour = method_label)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
      geom_segment(aes(x = low_plot, xend = high_plot,
                       y = method_label, yend = method_label),
                   linewidth = 0.65) +
      geom_point(aes(x = estimate_plot), size = 1.8) +
      facet_wrap(~dataset, ncol = 3, scales = "free_y") +
      x_scale +
      scale_colour_manual(values = method_colors, drop = FALSE) +
      labs(title = title, x = xlab, y = NULL, colour = NULL) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            strip.background = element_rect(fill = "grey92", colour = "grey70"),
            panel.grid.minor = element_blank())
    ggsave(file.path(eqtl1_plot_dir, paste0(filename, ".pdf")), p,
           width = 12, height = 7, useDingbats = FALSE)
    ggsave(file.path(eqtl1_plot_dir, paste0(filename, ".png")), p,
           width = 12, height = 7, dpi = 220)
  }
  make_eqtl1_forest(
    eqtl1, "pooled_or_bootstrap_median", "pooled_or_ci_low",
    "pooled_or_ci_high", "Dataset-level eQTL1 pooled odds ratio",
    "Pooled OR (gene x cell-type bootstrap 95% CI; log10 scale)",
    "eQTL1_dataset_pooled_OR_forest"
  )
  make_eqtl1_forest(
    eqtl1, "scent_mean_bootstrap_median", "scent_mean_ci_low",
    "scent_mean_ci_high", "Dataset-level eQTL1 SCENT-style mean",
    "Mean per-gene enrichment (gene x cell-type bootstrap 95% CI; pseudo-log scale)",
    "eQTL1_dataset_SCENT_mean_forest", x_transform = "pseudo_log"
  )
}

# eQTL2 ranking curves: pool counts at each top fraction.
curve_path <- file.path(eqtl2_root, "all_curves.tsv.gz")
if (file.exists(curve_path)) {
  curves <- fread(curve_path, showProgress = FALSE)
  curves <- merge(curves, panels, by = c("dataset", "celltype"), all = FALSE)
  curves <- curves[score_type == "primary" & method %in% rank_methods]
  if (nrow(curves)) {
    pooled_curves <- curves[, .(
      TP = sum(TP, na.rm = TRUE),
      n_selected = sum(n_selected, na.rm = TRUE),
      n_positive = sum(n_positive, na.rm = TRUE),
      n_units = sum(n_units, na.rm = TRUE),
      n_celltypes = uniqueN(celltype)
    ), by = .(dataset, method, top_fraction)]
    pooled_curves[, `:=`(
      recall = fifelse(n_positive > 0, TP / n_positive, NA_real_),
      precision = fifelse(n_selected > 0, TP / n_selected, NA_real_),
      baseline_precision = fifelse(n_units > 0, n_positive / n_units, NA_real_)
    )]
    pooled_curves[, enrichment := fifelse(
      baseline_precision > 0, precision / baseline_precision, NA_real_
    )]
    pooled_curves[, method_label := factor(unname(method_labels[method]),
                                           levels = method_levels)]
    y_max <- finite_quantile(pooled_curves$enrichment, 0.98)
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
    p <- ggplot(pooled_curves,
                aes(recall, enrichment, colour = method_label,
                    group = method_label)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
      geom_line(linewidth = 0.7) + geom_point(size = 1) +
      scale_colour_manual(values = method_colors, drop = FALSE) +
      facet_wrap(~dataset, ncol = 3, scales = "free") +
      coord_cartesian(ylim = c(0, y_max)) +
      labs(title = "Dataset-level eQTL2 ranked enrichment",
           x = "Recall of fine-mapped variant-gene links",
           y = "Precision enrichment over pooled assessable base rate",
           colour = NULL) +
      theme_bw(base_size = 10) + theme(legend.position = "bottom")
    ggsave(file.path(eqtl2_plot_dir, "eQTL2_dataset_ranked_enrichment_recall.pdf"),
           p, width = 12, height = 7, useDingbats = FALSE)
    ggsave(file.path(eqtl2_plot_dir, "eQTL2_dataset_ranked_enrichment_recall.png"),
           p, width = 12, height = 7, dpi = 220)
  }
}

# eQTL2 native significant sets: exact pooling of TP/FP across cell types.
native_path <- file.path(eqtl2_root, "all_native_significance.tsv")
if (file.exists(native_path)) {
  native <- fread(native_path, showProgress = FALSE)
  native <- merge(native, panels, by = c("dataset", "celltype"), all = FALSE)
  native <- native[method %in% plot_methods]
  if (nrow(native)) {
    pooled_native <- native[, .(
      n_celltypes = uniqueN(celltype),
      n_selected_pairs = sum(n_selected_pairs, na.rm = TRUE),
      TP = sum(TP, na.rm = TRUE), FP = sum(FP, na.rm = TRUE),
      n_units = sum(n_units, na.rm = TRUE),
      n_positive = sum(n_positive, na.rm = TRUE)
    ), by = .(dataset, method)]
    pooled_native[, `:=`(
      recall = fifelse(n_positive > 0, TP / n_positive, NA_real_),
      precision = fifelse(TP + FP > 0, TP / (TP + FP), NA_real_),
      baseline_precision = fifelse(n_units > 0, n_positive / n_units, NA_real_)
    )]
    pooled_native[, enrichment := fifelse(
      baseline_precision > 0, precision / baseline_precision, NA_real_
    )]
    pooled_native[, method_label := factor(unname(method_labels[method]),
                                           levels = method_levels)]
    p <- ggplot(pooled_native, aes(recall, enrichment, colour = method_label)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
      geom_point(size = 2.5) +
      scale_colour_manual(values = method_colors, drop = FALSE) +
      facet_wrap(~dataset, ncol = 3, scales = "free") +
      labs(title = "Dataset-level eQTL2 native significance",
           x = "Recall (native significant set)",
           y = "Enrichment (native significant set)", colour = NULL) +
      theme_bw(base_size = 10) + theme(legend.position = "bottom")
    ggsave(file.path(eqtl2_plot_dir, "eQTL2_dataset_native_significance.pdf"),
           p, width = 12, height = 7, useDingbats = FALSE)
    ggsave(file.path(eqtl2_plot_dir, "eQTL2_dataset_native_significance.png"),
           p, width = 12, height = 7, dpi = 220)
  }
}

message("Dataset-level plots written to: ", oneK1K_root)
