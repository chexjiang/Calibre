#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
base <- file.path(root, "results", "validation", "GWAS", "SCENT", "sakaue_style",
                  "brain2_pituitary_ukbb")
out <- file.path(base, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- fread(file.path(base, "results", "scent_sakaue_summary.tsv"))
x <- x[is.finite(mean_enrichment) & is.finite(ci_low) & is.finite(ci_high)]
x[, dataset_label := factor(dataset,
                            levels = c("Brain2_ALL", "Brain2_CTRL", "Brain2_AD", "Pituitary"))]
x[, method_label := factor(method_label, levels = c(
  "EpiLink-NB", "EpiLink-Poisson",
  "scMultiMap", "Signac (corrected)", "Signac (uncorrected)"))]

theme_bench <- theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "top", panel.spacing = unit(1.1, "lines"))

z <- x[celltype_scope == "ALL_REQUESTED"]
if (nrow(z)) {
  p <- ggplot(z, aes(method_label, mean_enrichment, color = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = .35, color = "grey45") +
    geom_errorbar(aes(ymin = pmax(ci_low, 0), ymax = ci_high), width = .16,
                  linewidth = .55, na.rm = TRUE) +
    geom_point(size = 2.8, na.rm = TRUE) +
    facet_wrap(~dataset_label, ncol = 2, scales = "free_y") +
    labs(title = "UK Biobank SCENT-style enrichment: pooled cell types",
         subtitle = "Brain2 ALL/CTRL/AD are separate analyses; Pituitary uses six major cell types",
         x = NULL, y = "Mean causal-variant enrichment", color = NULL) + theme_bench +
    guides(color = "none")
  ggsave(file.path(out, "ukbb_pooled_brain2_pituitary.pdf"), p,
         width = 12, height = 8, device = cairo_pdf)
}

z <- x[celltype_scope != "ALL_REQUESTED"]
if (nrow(z)) {
  z[, panel := paste(dataset, celltype_scope, sep = " | ")]
  z[, panel := factor(panel, levels = unique(panel))]
  p <- ggplot(z, aes(method_label, mean_enrichment, color = method_label)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = .3, color = "grey45") +
    geom_errorbar(aes(ymin = pmax(ci_low, 0), ymax = ci_high), width = .14,
                  linewidth = .45, na.rm = TRUE) +
    geom_point(size = 2.1, na.rm = TRUE) +
    facet_wrap(~panel, ncol = 2, scales = "free_y") +
    labs(title = "UK Biobank SCENT-style enrichment: cell-type-specific",
         subtitle = "Point = mean across UKBB traits; whisker = 95% trait bootstrap CI",
         x = NULL, y = "Mean causal-variant enrichment", color = NULL) + theme_bench +
    guides(color = "none")
  ggsave(file.path(out, "ukbb_celltype_brain2_pituitary.pdf"), p,
         width = 13, height = 16, device = cairo_pdf)
}

fwrite(x, file.path(out, "ukbb_scent_summary_readable.tsv"), sep = "\t")
message("Wrote UKBB Brain2/Pituitary plots and summary to ", out)
