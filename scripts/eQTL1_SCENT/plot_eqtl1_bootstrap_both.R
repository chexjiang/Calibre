#!/usr/bin/env Rscript

# Plot eQTL1 results from run_eqtl1_bootstrap_both.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) args[[1L]] else
  "results/validation/eQTL/OneK1K/eQTL1/all_results.tsv"
plot_dir <- if (length(args) >= 2L) args[[2L]] else
  "results/validation/eQTL/OneK1K/eQTL1"
requested_datasets <- if (length(args) >= 3L) args[3:length(args)] else
  c("pbmc", "bmmc", "Dogma_CTRL", "Dogma_STIM", "Arthritis", "SHARE")

if (!file.exists(input_path)) {
  stop("Missing input: ", input_path,
       ". Run run_eqtl1_bootstrap_both.R first.", call. = FALSE)
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
comparison_dir <- file.path(dirname(plot_dir), "figures")
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
x <- fread(input_path)

keep_datasets <- requested_datasets
x <- x[dataset %in% keep_datasets]
if (!nrow(x)) stop("No requested datasets found in ", input_path, call. = FALSE)

method_labels <- c(
  epilink_poisson = "EpiLink-Poisson",
  epilink_nb = "EpiLink-NB",
  poisson_glm = "Poisson GLM",
  nb_glm = "NB GLM",
  scmultimap = "scMultiMap",
  signac = "Signac (corrected)",
  signac_corrected = "Signac (corrected)",
  signac_uncorrected = "Signac (uncorrected)"
)
plot_methods <- c("epilink_poisson", "epilink_nb", "scmultimap",
                  "signac_corrected", "signac_uncorrected")
method_levels <- c("EpiLink-Poisson", "EpiLink-NB", "scMultiMap",
                   "Signac (corrected)", "Signac (uncorrected)")
method_colors <- c(
  "EpiLink-Poisson" = "#E47159",
  "EpiLink-NB" = "#C87878",
  "scMultiMap" = "#edd089",
  "Signac (corrected)" = "#54b2a2",
  "Signac (uncorrected)" = "#78A878"
)

x <- x[method %in% plot_methods]
if (!nrow(x)) stop("No requested methods found in ", input_path, call. = FALSE)

x[, method_label := factor(unname(method_labels[method]),
                           levels = method_levels)]
x[, dataset := factor(dataset, levels = keep_datasets)]
# Arthritis/B and SHAREseq/LCL use the same broad B-cell reference group as
# Bcell.  Collapse these labels in the figure so that the comparison is made
# in one Bcell facet rather than in separate B and LCL facets.
x[, celltype := fifelse(celltype %in% c("B", "LCL"), "Bcell", celltype)]
x[, celltype := factor(celltype,
                       levels = c("Bcell", "Tnk", "T", "Myeloid"))]

make_forest <- function(dat, estimate_col, low_col, high_col,
                        title, xlab, filename,
                        x_transform = c("log10", "pseudo_log")) {
  x_transform <- match.arg(x_transform)
  dat <- dat[is.finite(get(estimate_col)) &
             is.finite(get(low_col)) & is.finite(get(high_col))]
  if (x_transform == "log10") {
    dat <- dat[get(estimate_col) > 0 & get(low_col) > 0 &
               get(high_col) > 0]
  }
  dat[, estimate := get(estimate_col)]
  dat[, ci_low := get(low_col)]
  dat[, ci_high := get(high_col)]
  x_scale <- if (x_transform == "log10") {
    scale_x_log10()
  } else {
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10))
  }
  p <- ggplot(dat, aes(y = method_label, colour = method_label)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
    geom_segment(aes(x = ci_low, xend = ci_high,
                     y = method_label, yend = method_label),
                 linewidth = 0.65, na.rm = TRUE) +
    geom_point(aes(x = estimate), size = 1.8, na.rm = TRUE) +
    facet_grid(dataset ~ celltype, scales = "free_y", drop = TRUE) +
    x_scale +
    scale_colour_manual(values = method_colors, drop = FALSE) +
    labs(title = title, x = xlab, y = NULL, colour = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), p,
         width = 13, height = 10, useDingbats = FALSE)
  ggsave(file.path(plot_dir, paste0(filename, ".png")), p,
         width = 13, height = 10, dpi = 220)
}

make_forest(
  x,
  "pooled_or_bootstrap_median", "pooled_or_ci_low", "pooled_or_ci_high",
  "eQTL1 pooled odds ratio",
  "Pooled OR (gene-bootstrap 95% CI; log scale)",
  "eQTL1_pooled_OR_forest",
  x_transform = "log10"
)

make_forest(
  x,
  "scent_mean_bootstrap_median", "scent_mean_ci_low",
  "scent_mean_ci_high",
  "eQTL1 SCENT-style mean enrichment",
  "Mean per-gene enrichment (gene-bootstrap 95% CI; pseudo-log scale)",
  "eQTL1_SCENT_mean_forest",
  x_transform = "pseudo_log"
)

comparison <- x[is.finite(pooled_or_bootstrap_median) &
                is.finite(scent_mean_bootstrap_median) &
                pooled_or_bootstrap_median > 0 &
                scent_mean_bootstrap_median > 0]
p_compare <- ggplot(
  comparison,
  aes(x = pooled_or_bootstrap_median,
      y = scent_mean_bootstrap_median,
      colour = method_label)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 2, alpha = 0.85) +
  facet_grid(dataset ~ celltype, drop = TRUE) +
  scale_x_log10() +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_colour_manual(values = method_colors, drop = FALSE) +
  labs(
    title = "Pooled OR versus SCENT-style mean enrichment",
    x = "Pooled OR (bootstrap median; log scale)",
    y = "SCENT-style mean enrichment (bootstrap median; pseudo-log scale)",
    colour = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(comparison_dir, "eQTL1_pooled_vs_SCENT_mean.pdf"), p_compare,
       width = 13, height = 10, useDingbats = FALSE)
ggsave(file.path(comparison_dir, "eQTL1_pooled_vs_SCENT_mean.png"), p_compare,
       width = 13, height = 10, dpi = 220)

method_summary <- x[, .(
  n_dataset_celltype = .N,
  pooled_or_median_across_panels = median(pooled_or_bootstrap_median,
                                          na.rm = TRUE),
  scent_mean_median_across_panels = median(scent_mean_bootstrap_median,
                                           na.rm = TRUE)
), by = .(method, method_label)]
setorder(method_summary, -scent_mean_median_across_panels)
fwrite(method_summary,
       file.path(plot_dir, "eQTL1_method_summary.tsv"), sep = "\t")

message("Plots written to: ", plot_dir)
