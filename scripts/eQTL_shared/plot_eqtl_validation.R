#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

plot_dir <- "results/validation/eQTL/OneK1K"
eqtl1_plot_dir <- file.path(plot_dir, "eQTL1")
eqtl2_plot_dir <- file.path(plot_dir, "eQTL2")
dir.create(eqtl1_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(eqtl2_plot_dir, recursive = TRUE, showWarnings = FALSE)
method_labels <- c(
  epilink_poisson = "EpiLink-Poisson", epilink_nb = "EpiLink-NB",
  poisson_glm = "Poisson GLM", nb_glm = "NB GLM",
  scmultimap = "scMultiMap", signac = "Signac (corrected)",
  signac_corrected = "Signac (corrected)",
  signac_uncorrected = "Signac (uncorrected)"
)
plot_methods <- c("epilink_poisson", "epilink_nb", "scmultimap",
                  "signac_corrected", "signac_uncorrected")
rank_methods <- c("epilink_poisson", "epilink_nb", "scmultimap", "signac")
plot_levels <- c("EpiLink-Poisson", "EpiLink-NB", "scMultiMap",
                 "Signac (corrected)", "Signac (uncorrected)")
method_colours <- c(
  "EpiLink-Poisson" = "#E47159",
  "EpiLink-NB" = "#C87878",
  "scMultiMap" = "#edd089",
  "Signac (corrected)" = "#54b2a2",
  "Signac (uncorrected)" = "#78A878"
)

config <- yaml::read_yaml("config/eqtl_validation.yml")
configured_panels <- data.table::rbindlist(lapply(names(config$benchmarks), function(dataset) {
  data.table::data.table(dataset = dataset,
                         celltype = names(config$benchmarks[[dataset]]))
}))
collapse_bcell_label <- function(x) {
  data.table::fifelse(x %in% c("B", "LCL"), "Bcell", x)
}

scent_path <- "results/validation/eQTL/OneK1K/eQTL1/all_results.tsv"
if (file.exists(scent_path)) {
  x <- fread(scent_path)[method %in% plot_methods]
  x <- merge(x, configured_panels, by = c("dataset", "celltype"), all = FALSE)
  x[, method_label := factor(unname(method_labels[method]), levels = plot_levels)]
  x[, celltype := collapse_bcell_label(celltype)]
  x[, panel := paste(dataset, celltype, sep = " / ")]
  p <- ggplot(x, aes(method_label, pooled_or_bootstrap_median,
                     colour = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_errorbar(aes(ymin = pooled_or_ci_low, ymax = pooled_or_ci_high),
                  width = 0.15) +
    geom_point(size = 2) +
    facet_wrap(~panel, scales = "free_x") +
  scale_colour_manual(values = method_colours, drop = FALSE) +
    scale_y_log10() +
    labs(x = NULL, y = "Pooled eQTL1 odds ratio (gene bootstrap)",
         colour = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  ggsave(file.path(eqtl1_plot_dir, "eQTL1_SCENT_odds_ratios.pdf"), p,
         width = 12, height = 8)
  ggsave(file.path(eqtl1_plot_dir, "eQTL1_SCENT_odds_ratios.png"), p,
         width = 12, height = 8, dpi = 180)
}

curve_path <- "results/validation/eQTL/OneK1K/eQTL2/all_curves.tsv.gz"
if (file.exists(curve_path)) {
  x <- fread(curve_path)[score_type == "primary" & method %in% rank_methods]
  x <- merge(x, configured_panels, by = c("dataset", "celltype"), all = FALSE)
  rank_levels <- plot_levels[plot_levels != "Signac (uncorrected)"]
  rank_colours <- method_colours[names(method_colours) != "Signac (uncorrected)"]
  x[, method_label := factor(unname(method_labels[method]), levels = rank_levels)]
  x[, celltype := collapse_bcell_label(celltype)]
  x[, panel := paste(dataset, celltype, sep = " / ")]
  p <- ggplot(x, aes(recall, enrichment, colour = method_label,
                     group = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_line(linewidth = 0.7) +
  geom_point(size = 1) +
    scale_colour_manual(values = rank_colours, drop = FALSE) +
    facet_wrap(~panel, scales = "free") +
    coord_cartesian(ylim = c(0, stats::quantile(x$enrichment, 0.98, na.rm = TRUE))) +
    labs(x = "Recall of fine-mapped variant-gene links",
         y = "Precision enrichment over assessable base rate", colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(eqtl2_plot_dir, "eQTL2_ranked_enrichment_recall.pdf"), p,
         width = 12, height = 8)
  ggsave(file.path(eqtl2_plot_dir, "eQTL2_ranked_enrichment_recall.png"), p,
         width = 12, height = 8, dpi = 180)
}

native_path <- "results/validation/eQTL/OneK1K/eQTL2/all_native_significance.tsv"
if (file.exists(native_path)) {
  x <- fread(native_path)[method %in% plot_methods]
  x <- merge(x, configured_panels, by = c("dataset", "celltype"), all = FALSE)
  x[, method_label := factor(unname(method_labels[method]), levels = plot_levels)]
  x[, celltype := collapse_bcell_label(celltype)]
  x[, panel := paste(dataset, celltype, sep = " / ")]
  p <- ggplot(x, aes(recall, enrichment, colour = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
  geom_point(size = 2.5) +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    facet_wrap(~panel, scales = "free") +
    labs(x = "Recall (native significant set)",
         y = "Enrichment (native significant set)", colour = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(eqtl2_plot_dir, "eQTL2_native_significance.pdf"), p,
         width = 12, height = 8)
  ggsave(file.path(eqtl2_plot_dir, "eQTL2_native_significance.png"), p,
         width = 12, height = 8, dpi = 180)
}
