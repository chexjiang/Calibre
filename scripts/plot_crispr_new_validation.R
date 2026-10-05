#!/usr/bin/env Rscript

# Plot pgBoost-style CRISPR SNP-gene benchmark results. PDF output only.

raw_args <- commandArgs(trailingOnly = TRUE)
dataset_args <- raw_args[!grepl("^--", raw_args)]

file_arguments <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arguments)) {
  script_path <- normalizePath(sub("^--file=", "", file_arguments[[1L]]))
  project_root <- normalizePath(file.path(dirname(script_path), ".."),
                               mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
  if (!file.exists(file.path(project_root, "R", "validation_plots.R"))) {
    project_root <- normalizePath(file.path(project_root, ".."),
                                  mustWork = TRUE)
  }
}
setwd(project_root)

suppressPackageStartupMessages(library(ggplot2))

method_order <- c(
  "EpiLink-Poisson", "EpiLink-NB", "Poisson-GLM", "NB-GLM",
  "scMultiMap", "scMultiMap-noBatch", "Signac (corrected)",
  "Signac (uncorrected)"
)
method_colours <- c(
  "EpiLink-Poisson" = "#D55E5E", "EpiLink-NB" = "#8E6BBE",
  "Poisson-GLM" = "#CC79A7", "NB-GLM" = "#377EB8",
  "scMultiMap" = "#E69F00", "scMultiMap-noBatch" = "#009E73",
  "Signac (corrected)" = "#6E6E6E", "Signac (uncorrected)" = "#B07AA1"
)
preferred_celltypes <- c(
  "Bcell", "Dendritic", "Erythroid", "Myeloid", "T", "Tnk", "All"
)
distance_order <- c(1000L, 10000L, 100000L)
distance_labels <- c(">1 kb", ">10 kb", ">100 kb")

real_root <- file.path(project_root, "results", "real")
validation_root <- file.path(project_root, "results", "validation", "CRISPR",
                             "CRISPR_pgBoost")
figure_dir <- file.path(validation_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
available <- list.dirs(real_root, recursive = FALSE, full.names = FALSE)
available <- setdiff(available, c("data_explore", "SatpathyK562"))
if (!length(dataset_args)) {
  datasets <- available
} else {
  datasets <- vapply(dataset_args, function(dataset) {
    hit <- available[tolower(available) == tolower(dataset)]
    if (length(hit) != 1L) {
      stop("No results/real dataset named ", dataset, call. = FALSE)
    }
    hit
  }, character(1L))
}

theme_validation <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = element_rect(fill = "grey95", colour = "grey70")
    )
}

for (dataset in unique(datasets)) {
  metrics_file <- file.path(
    validation_root, paste0(dataset, "_enrichment_metrics.tsv")
  )
  curves_file <- file.path(
    validation_root, paste0(dataset, "_enrichment_recall_curves.tsv.gz")
  )
  if (!file.exists(metrics_file) || !file.exists(curves_file)) {
    message("Skipping ", dataset, ": missing CRISPR_pgBoost validation output")
    next
  }

  metrics <- read.delim(metrics_file, stringsAsFactors = FALSE,
                        check.names = FALSE)
  curves <- read.delim(gzfile(curves_file), stringsAsFactors = FALSE,
                       check.names = FALSE)
  # Bootstrap percentile columns receive names such as
  # average_enrichment_low.2.5% when written by data.table.  Normalize those
  # names so both bootstrap and --n-bootstrap=0 outputs use the same plotting
  # code (the latter simply has NA intervals).
  names(metrics) <- sub("^average_enrichment_low.*$",
                         "average_enrichment_low", names(metrics))
  names(metrics) <- sub("^average_enrichment_high.*$",
                         "average_enrichment_high", names(metrics))
  names(metrics) <- sub("^auprc_low.*$", "auprc_low", names(metrics))
  names(metrics) <- sub("^auprc_high.*$", "auprc_high", names(metrics))
  observed_methods <- method_order[method_order %in% unique(metrics$method)]
  observed_celltypes <- unique(as.character(metrics$celltype))
  celltype_order <- c(
    preferred_celltypes[preferred_celltypes %in% observed_celltypes],
    sort(setdiff(observed_celltypes, preferred_celltypes))
  )
  format_data <- function(x) {
    x$method <- factor(x$method, levels = observed_methods)
    x$celltype <- factor(x$celltype, levels = celltype_order)
    x$distance <- factor(
      as.integer(x$distance_threshold), levels = distance_order,
      labels = distance_labels
    )
    x
  }
  metrics <- format_data(metrics)
  curves <- format_data(curves)

  facet_height <- max(6, 2.4 * ceiling(length(celltype_order) / 2))

  p_enrichment <- ggplot(
    metrics,
    aes(x = distance, y = average_enrichment, colour = method, group = method)
  ) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_line(position = position_dodge(width = 0.35), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.35), size = 1.8) +
    facet_wrap(~celltype, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    labs(
      title = paste(dataset, "CRISPR SNP-gene benchmark"),
      subtitle = "Area under enrichment-recall curve (random ranking = 1)",
      x = "Minimum SNP-to-TSS distance", y = "Average enrichment",
      colour = NULL
    ) + theme_validation()

  p_auprc <- ggplot(
    metrics, aes(x = distance, y = auprc, colour = method, group = method)
  ) +
    geom_line(position = position_dodge(width = 0.35), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.35), size = 1.8) +
    facet_wrap(~celltype, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    labs(
      title = paste(dataset, "CRISPR SNP-gene benchmark"),
      x = "Minimum SNP-to-TSS distance", y = "Area under precision-recall curve",
      colour = NULL
    ) + theme_validation()

  p_curve <- ggplot(
    curves, aes(x = recall, y = enrichment, colour = method, group = method)
  ) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_step(linewidth = 0.55, alpha = 0.9) +
    facet_grid(celltype ~ distance, scales = "free_y") +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = paste(dataset, "CRISPR enrichment-recall curves"),
      x = "Recall of validated SNP-gene links", y = "Fold enrichment",
      colour = NULL
    ) + theme_validation()

  p_or <- ggplot(
    metrics, aes(x = distance, y = odds_ratio, colour = method, group = method)
  ) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_line(position = position_dodge(width = 0.35), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.35), size = 1.8) +
    facet_wrap(~celltype, ncol = 2) +
    scale_y_log10() +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    labs(
      title = paste(dataset, "CRISPR enrichment among top 5% links"),
      x = "Minimum SNP-to-TSS distance", y = "Fisher odds ratio",
      colour = NULL
    ) + theme_validation()

  p_pvalue <- ggplot(
    metrics,
    aes(x = distance, y = -log10(pmax(pval, .Machine$double.xmin)),
        colour = method, group = method)
  ) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey55") +
    geom_line(position = position_dodge(width = 0.35), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.35), size = 1.8) +
    facet_wrap(~celltype, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    labs(
      title = paste(dataset, "CRISPR top-5% enrichment significance"),
      x = "Minimum SNP-to-TSS distance", y = expression(-log[10](italic(P))),
      colour = NULL
    ) + theme_validation()

  ggplot2::ggsave(file.path(figure_dir,
                             paste0(dataset, "_average_enrichment_crispr_pgboost.pdf")),
                  p_enrichment, width = 12, height = facet_height)
  ggplot2::ggsave(file.path(figure_dir,
                             paste0(dataset, "_auprc_crispr_pgboost.pdf")),
                  p_auprc, width = 12, height = facet_height)
  ggplot2::ggsave(file.path(figure_dir,
                             paste0(dataset, "_enrichment_recall_crispr_pgboost.pdf")),
                  p_curve, width = 14, height = max(7, 2.2 * length(celltype_order)))
  ggplot2::ggsave(file.path(figure_dir,
                             paste0(dataset, "_odds_ratio_crispr_pgboost.pdf")),
                  p_or, width = 12, height = facet_height)
  ggplot2::ggsave(file.path(figure_dir,
                             paste0(dataset, "_pvalue_crispr_pgboost.pdf")),
                  p_pvalue, width = 12, height = facet_height)
  message("Wrote CRISPR_pgBoost PDFs for ", dataset, " to ", figure_dir)
}
