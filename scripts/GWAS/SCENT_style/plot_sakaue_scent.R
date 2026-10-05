#!/usr/bin/env Rscript

# Plot the Sakaue Fig. 3-style dataset-level summaries. FinnGen and UKBB are
# written as separate PDFs; each point is a dataset-level union across primary
# cell types, with a trait-bootstrap 95% CI.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
base <- file.path(root, "results", "validation", "GWAS", "SCENT_style", "UKBB", "sakaue_style")
out_dir <- file.path(base, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- fread(file.path(base, "results", "scent_sakaue_summary.tsv"))
x <- x[is.finite(mean_enrichment) & is.finite(ci_low) & is.finite(ci_high)]
x[, `:=`(
  dataset_label = fifelse(dataset == "Dogma_CTRL", "Dogma CTRL", fifelse(dataset == "Dogma_STIM", "Dogma STIM", dataset)),
  family = fifelse(method == "epilink_nb", "EpiLink-NB", "EpiLink-Poisson"),
  correction_label = fifelse(correction, "correction = TRUE", "correction = FALSE")
)]
x[, dataset_label := factor(dataset_label, levels = unique(dataset_label))]

plot_one <- function(source_name) {
  z <- x[source == source_name]
  if (!nrow(z)) return(invisible(NULL))
  dodge <- position_dodge(width = 0.58)
  p <- ggplot(z, aes(x = dataset_label, y = mean_enrichment,
                     color = correction_label, shape = correction_label,
                     group = correction_label)) +
    geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dashed", color = "grey45") +
    geom_errorbar(aes(ymin = pmax(ci_low, 0), ymax = ci_high),
                  width = 0.16, linewidth = 0.55, position = dodge) +
    geom_point(size = 2.7, stroke = 0.8, position = dodge) +
    facet_wrap(~family, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("correction = TRUE" = "#1b6ca8", "correction = FALSE" = "#d97706")) +
    scale_shape_manual(values = c("correction = TRUE" = 16, "correction = FALSE" = 17)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = paste0("SCENT-style GWAS enrichment: ", source_name),
      subtitle = "Dataset-level union of primary cell types; point = mean across traits; whisker = 95% trait bootstrap CI",
      x = "Dataset (all primary cell types pooled within dataset)",
      y = "Mean causal-variant enrichment",
      color = NULL, shape = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.title = element_text(face = "plain"),
      strip.background = element_blank(), strip.text = element_text(face = "bold"),
      legend.position = "top", plot.title = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines")
    )
  ggsave(file.path(out_dir, paste0("scent_sakaue_", tolower(source_name), ".pdf")),
         p, width = 12, height = 6.5, units = "in", device = cairo_pdf)
  invisible(p)
}

plot_one("FinnGen")
plot_one("UKBB")

setorder(x, source, family, dataset_label, correction)
fwrite(x, file.path(out_dir, "scent_sakaue_summary_readable.tsv"), sep = "\t")
message("Wrote FinnGen and UKBB PDFs to ", out_dir)
