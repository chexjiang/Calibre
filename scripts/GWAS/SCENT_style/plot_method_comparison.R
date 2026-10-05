#!/usr/bin/env Rscript

# One PDF comparing the five requested methods. Each panel has one method and
# one GWAS source, so small enrichments remain readable even when one dataset
# has a much larger value.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
base <- file.path(root, "results", "validation", "GWAS", "SCENT_style", "UKBB", "sakaue_style", "method_comparison")
out_dir <- file.path(base, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- fread(file.path(base, "results", "scent_sakaue_summary.tsv"))
x <- x[is.finite(mean_enrichment) & is.finite(ci_low) & is.finite(ci_high)]
x[, dataset_label := fifelse(dataset == "Dogma_CTRL", "Dogma CTRL",
                             fifelse(dataset == "Dogma_STIM", "Dogma STIM", dataset))]
method_levels <- c("EpiLink-NB (correction = FALSE)",
                   "EpiLink-Poisson (correction = FALSE)",
                   "scMultiMap", "Signac (corrected)",
                   "Signac (uncorrected)")
x[, method_label := factor(method_label, levels = method_levels)]
x[, dataset_label := factor(dataset_label,
                            levels = c("pbmc", "bmmc", "Dogma CTRL", "Dogma STIM",
                                       "NEAT", "Arthritis", "Brain"))]
x[, source := factor(source, levels = c("FinnGen", "UKBB"))]

p <- ggplot(x, aes(x = dataset_label, y = mean_enrichment)) +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = pmax(ci_low, 0), ymax = ci_high),
                width = 0.18, linewidth = 0.55, color = "#2166ac") +
  geom_point(size = 2.5, color = "#2166ac") +
  facet_grid(source ~ method_label, scales = "free_y", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "Sakaue-style GWAS enrichment: EpiLink vs scMultiMap and Signac",
    subtitle = "Each point is the mean enrichment across GWAS traits; whisker = 95% trait-bootstrap CI",
    x = "Dataset (all primary cell types pooled within dataset)",
    y = "Mean causal-variant enrichment"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
    axis.title = element_text(face = "plain"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold"),
    panel.spacing = unit(0.9, "lines")
  )

ggsave(file.path(out_dir, "scent_sakaue_epilink_scmultimap_signac.pdf"),
       p, width = 18, height = 8.5, units = "in", device = cairo_pdf)

setorder(x, source, method_label, dataset_label)
fwrite(x, file.path(out_dir, "scent_sakaue_method_comparison_summary_readable.tsv"), sep = "\t")
message("Wrote ", file.path(out_dir, "scent_sakaue_epilink_scmultimap_signac.pdf"))
