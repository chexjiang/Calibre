#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

root <- "results/validation/eQTL/GTEx"
fig <- file.path(root, "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
labels <- c(epilink_poisson = "EpiLink-Poisson", epilink_nb = "EpiLink-NB",
            scmultimap = "scMultiMap", signac_corrected = "Signac (corrected)",
            signac_uncorrected = "Signac (uncorrected)")
levels <- unname(labels)
cols <- c("EpiLink-Poisson" = "#E47159", "EpiLink-NB" = "#C87878",
          "scMultiMap" = "#edd089", "Signac (corrected)" = "#54b2a2",
          "Signac (uncorrected)" = "#78A878")

e1 <- fread(file.path(root, "eQTL1", "all_results.tsv"))
e1[, method_label := factor(unname(labels[method]), levels = levels)]
e1[, scope_label := fifelse(scope == "panel", paste(dataset, celltype, sep = " / "), dataset)]

save_forest <- function(x, estimate, low, high, title, ylab, file, facet = NULL,
                        y_transform = "pseudo") {
  x <- x[is.finite(get(estimate)) & is.finite(get(low)) & is.finite(get(high))]
  if (!nrow(x)) return(invisible(NULL))
  x[, `:=`(est = get(estimate), lo = get(low), hi = get(high))]
  dodge <- position_dodge(width = .65)
  p <- ggplot(x, aes(x = if (is.null(facet)) dataset else scope_label,
                     y = est, colour = method_label, group = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = .18) +
    geom_point(position = dodge, size = 1.7) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(title = title, x = NULL, y = ylab, colour = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "bottom",
          panel.grid.minor = element_blank())
  if (!is.null(facet)) p <- p + facet_wrap(as.formula(paste("~", facet)), scales = "free_x", ncol = 4)
  p <- p + if (y_transform == "log") scale_y_log10() else scale_y_continuous(trans = pseudo_log_trans(base = 10))
  ggsave(file.path(fig, paste0(file, ".pdf")), p, width = 15, height = 9, useDingbats = FALSE)
  ggsave(file.path(fig, paste0(file, ".png")), p, width = 15, height = 9, dpi = 220)
}

panel <- e1[scope == "panel"]
pool <- e1[scope == "dataset_pool"]
save_forest(panel, "scent_mean_bootstrap_median", "scent_mean_ci_low", "scent_mean_ci_high",
            "GTEx eQTL1: SCENT-style mean enrichment by panel", "Mean per-gene enrichment (95% gene-bootstrap CI)",
            "eQTL1_panel_SCENT_mean", facet = "dataset")
save_forest(pool, "scent_mean_bootstrap_median", "scent_mean_ci_low", "scent_mean_ci_high",
            "GTEx eQTL1: dataset-pooled SCENT-style mean", "Mean per-gene enrichment (95% gene-bootstrap CI)",
            "eQTL1_dataset_pool_SCENT_mean")
save_forest(panel, "pooled_or_bootstrap_median", "pooled_or_ci_low", "pooled_or_ci_high",
            "GTEx eQTL1: pooled odds ratio by panel", "Pooled odds ratio (95% gene-bootstrap CI)",
            "eQTL1_panel_pooled_OR", facet = "dataset", y_transform = "log")
save_forest(pool, "pooled_or_bootstrap_median", "pooled_or_ci_low", "pooled_or_ci_high",
            "GTEx eQTL1: dataset-pooled odds ratio", "Pooled odds ratio (95% gene-bootstrap CI)",
            "eQTL1_dataset_pool_OR", y_transform = "log")

tissue <- fread(file.path(root, "eQTL1", "all_tissue_results.tsv.gz"))
tissue[, method_label := factor(unname(labels[method]), levels = levels)]
tissue[, scope_label := fifelse(scope == "panel", paste(dataset, celltype, sep = " / "), dataset)]
tissue[, value := fifelse(is.finite(scent_mean) & scent_mean > 0, log10(scent_mean), NA_real_)]
for (sc in c("panel", "dataset_pool")) {
  z <- tissue[scope == sc]
  if (!nrow(z)) next
  z[, scope_label := factor(scope_label, levels = rev(unique(scope_label)))]
  p <- ggplot(z, aes(tissue, scope_label, fill = value)) +
    geom_tile(colour = "white", linewidth = .12) + facet_wrap(~method_label, ncol = 1) +
    scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426",
                         midpoint = 0, na.value = "grey92", name = "log10\nmean enrichment") +
    labs(title = paste("GTEx tissue-specific eQTL1 enrichment:", sc), x = NULL, y = NULL) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),
          panel.grid = element_blank(), legend.position = "bottom")
  name <- if (sc == "panel") "eQTL1_tissue_specific_panel_heatmap" else "eQTL1_tissue_specific_dataset_pool_heatmap"
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = 18, height = if (sc == "panel") 18 else 12, useDingbats = FALSE)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = 18, height = if (sc == "panel") 18 else 12, dpi = 220)
}

e2 <- fread(file.path(root, "eQTL2", "all_curves.tsv.gz"))[score_type == "primary"]
e2[, method_label := factor(ifelse(method == "signac", "Signac (corrected)", unname(labels[method])),
                            levels = levels)]
e2[, scope_label := fifelse(scope == "panel", paste(dataset, celltype, sep = " / "), dataset)]
rank_levels <- c("EpiLink-Poisson", "EpiLink-NB", "scMultiMap", "Signac (corrected)")
rank_cols <- cols[rank_levels]
for (sc in c("panel", "dataset_pool")) {
  z <- e2[scope == sc]
  if (!nrow(z)) next
  p <- ggplot(z, aes(recall, enrichment, colour = method_label, group = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_line(linewidth = .65) + geom_point(size = .7) +
    scale_colour_manual(values = rank_cols, breaks = rank_levels, drop = FALSE) +
    facet_wrap(~scope_label, ncol = if (sc == "panel") 4 else 3, scales = "free") +
    labs(title = paste("GTEx eQTL2 ranked enrichment:", sc),
         subtitle = "PIP > 0.5; EpiLink ranks use abs(beta) / abs(zval); Signac continuous score is shared by both native views",
         x = "Recall of fine-mapped variant-gene units", y = "Precision enrichment", colour = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  name <- if (sc == "panel") "eQTL2_ranked_panel" else "eQTL2_ranked_dataset_pool"
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = 16, height = if (sc == "panel") 14 else 9, useDingbats = FALSE)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = 16, height = if (sc == "panel") 14 else 9, dpi = 220)
}

native <- fread(file.path(root, "eQTL2", "all_native_significance.tsv"))
native[, method_label := factor(unname(labels[method]), levels = levels)]
native[, scope_label := fifelse(scope == "panel", paste(dataset, celltype, sep = " / "), dataset)]
for (sc in c("panel", "dataset_pool")) {
  z <- native[scope == sc]
  if (!nrow(z)) next
  p <- ggplot(z, aes(recall, enrichment, colour = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") + geom_point(size = 1.7) +
    scale_colour_manual(values = cols, drop = FALSE) +
    facet_wrap(~scope_label, ncol = if (sc == "panel") 4 else 3, scales = "free") +
    labs(title = paste("GTEx eQTL2 native significance:", sc),
         x = "Recall of PIP > 0.5 units", y = "Native significance enrichment", colour = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  name <- if (sc == "panel") "eQTL2_native_panel" else "eQTL2_native_dataset_pool"
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = 16, height = if (sc == "panel") 14 else 9, useDingbats = FALSE)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = 16, height = if (sc == "panel") 14 else 9, dpi = 220)
}

fwrite(e1, file.path(fig, "eQTL1_plot_data.tsv.gz"), sep = "\t")
fwrite(e2, file.path(fig, "eQTL2_ranked_plot_data.tsv.gz"), sep = "\t")
fwrite(native, file.path(fig, "eQTL2_native_plot_data.tsv"), sep = "\t")
message("GTEx eQTL figures written to ", fig)
