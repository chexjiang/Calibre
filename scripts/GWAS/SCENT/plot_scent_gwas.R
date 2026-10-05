#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
gwas_dir <- file.path(root, "results", "validation", "GWAS")
indir <- Sys.getenv("SCENT_RESULTS_DIR", file.path(gwas_dir, "SCENT_style", "testing", "results"))
outdir <- Sys.getenv("SCENT_FIGURE_DIR", file.path(gwas_dir, "SCENT_style", "testing", "figures"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- fread(file.path(indir, "scent_panel_summary.tsv"))
pt <- fread(file.path(indir, "scent_per_trait.tsv"))
methods_path <- Sys.getenv("GWAS_METHODS_FILE", file.path(gwas_dir, "SCENT_style", "config", "methods.tsv"))
method_order <- fread(methods_path)$label
exclude_labels <- trimws(strsplit(Sys.getenv("GWAS_PLOT_EXCLUDE_METHODS", "NB-GLM,Poisson-GLM"), ",", fixed = TRUE)[[1]])
method_order <- setdiff(method_order, exclude_labels)
x <- x[!method_label %in% exclude_labels]
pt <- pt[!method_label %in% exclude_labels]
x[, method_label := factor(method_label, levels = method_order)]
x[, dataset_celltype := paste(dataset, celltype, sep = " | ")]
x[, log2_enrichment := log2(mean_enrichment)]

theme_gwas <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
        axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")

for (pn in unique(x$panel)) {
  z <- x[panel == pn]
  z[, dataset_celltype := factor(dataset_celltype, levels = rev(unique(dataset_celltype)))]
  p <- ggplot(z, aes(method_label, dataset_celltype, fill = log2_enrichment)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_point(aes(size = pooled_causal_in), shape = 21, fill = NA, color = "black", stroke = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         na.value = "grey85", name = expression(log[2](mean~enrichment))) +
    scale_size_continuous(range = c(0.5, 4), name = "Causal hits") +
    labs(title = paste("SCENT-style GWAS enrichment:", pn), x = NULL, y = NULL,
         subtitle = "Fill: arithmetic mean across matched traits; circle: fine-mapped variants overlapping linked peaks") +
    theme_gwas
  h <- max(4.5, 0.31 * nrow(z) + 2)
  ggsave(file.path(outdir, paste0("scent_", pn, "_heatmap.pdf")), p, width = 8.2, height = h)
}

immune <- x[panel == "immune" & interpretation == "primary"]
immune[, dataset := factor(dataset, levels = unique(dataset))]
p_forest <- ggplot(immune, aes(method_label, mean_enrichment, color = celltype, group = celltype)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey45") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, position = position_dodge(width = 0.55), na.rm = TRUE) +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  scale_y_log10() +
  facet_wrap(~dataset, scales = "free_x") +
  labs(title = "Immune GWAS enrichment in linked peaks", subtitle = "Mean of RA, T1D, IBD and MS; 95% CI from 1,000 trait bootstraps",
       x = NULL, y = "Mean enrichment (log scale)", color = "Cell type") +
  theme_gwas
ggsave(file.path(outdir, "scent_immune_primary_forest.pdf"), p_forest, width = 12, height = 7)

pt[, dataset_celltype := paste(dataset, celltype, sep = " | ")]
p_counts <- ggplot(pt, aes(n_common_in_annotation, n_causal_in_annotation, color = method_label)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  facet_wrap(~trait_label, scales = "free") +
  labs(title = "Annotation size versus fine-mapped GWAS overlap", x = "Common variants in linked peaks (log scale)",
       y = "Fine-mapped variants in linked peaks", color = "Method") + theme_gwas
ggsave(file.path(outdir, "scent_overlap_counts.pdf"), p_counts, width = 11, height = 7)
