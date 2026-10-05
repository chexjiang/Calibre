#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
base <- file.path(root, "results", "validation", "GWAS", "sLDSC")
results_dir <- Sys.getenv("SLDSC_RESULTS_DIR", file.path(base, "results"))
out <- Sys.getenv("SLDSC_FIGURE_DIR", file.path(base, "figures"))
f <- file.path(results_dir, "sldsc_all_results.tsv")
if (!file.exists(f)) stop("Missing ", f, "; run run_sldsc.sh first")
x <- fread(f)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
methods_path <- Sys.getenv("GWAS_METHODS_FILE", file.path(root, "results", "validation", "GWAS", "SCENT_style", "config", "methods.tsv"))
methods <- fread(methods_path)
exclude_labels <- trimws(strsplit(Sys.getenv("GWAS_PLOT_EXCLUDE_METHODS", "NB-GLM,Poisson-GLM"), ",", fixed = TRUE)[[1]])
x <- x[!method_label %in% exclude_labels]
methods <- methods[!label %in% exclude_labels]
x[, method_label := factor(method_label, levels = methods$label)]
x[, dataset_celltype := paste(dataset, celltype, sep = " | ")]
x[, log2_enrichment := log2(Enrichment)]

theme_gwas <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1))
for (tr in unique(x$trait)) {
  z <- x[trait == tr]
  z[, dataset_celltype := factor(dataset_celltype, levels = rev(unique(dataset_celltype)))]
  p <- ggplot(z, aes(method_label, dataset_celltype, fill = log2_enrichment)) +
    geom_tile(color = "white", linewidth = .25) +
    geom_point(aes(size = -log10(Enrichment_p)), shape = 21, fill = NA, stroke = .25) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
    scale_size_continuous(range = c(.5, 4)) +
    labs(title = paste("S-LDSC heritability enrichment:", z$trait_label[1]),
         subtitle = "baseline-LD v2.2 + one linked-peak annotation; circle size is -log10(P)",
         x = NULL, y = NULL, fill = expression(log[2](enrichment)), size = expression(-log[10](P))) + theme_gwas
  h <- max(4.5, .31 * nrow(z) + 2)
  ggsave(file.path(out, paste0("sldsc_", tr, "_heatmap.pdf")), p, width = 8.3, height = h)
}

primary <- x[interpretation == "primary"]
p <- ggplot(primary, aes(method_label, Enrichment, color = celltype)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") + geom_point(position = position_jitter(width = .12), alpha = .75) +
  scale_y_log10() + facet_grid(trait_label ~ dataset, scales = "free_x", space = "free_x") +
  labs(title = "S-LDSC enrichment across primary matched cell types", x = NULL, y = "Heritability enrichment (log scale)", color = "Cell type") + theme_gwas
ggsave(file.path(out, "sldsc_primary_overview.pdf"), p, width = 14, height = 10)
