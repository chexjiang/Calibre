#!/usr/bin/env Rscript

# Collect the current GWAS validation results into one stable report directory
# and generate pooled/cell-type comparison figures.  This script deliberately
# does not recompute SCENT or sLDSC; it makes a reproducible presentation layer
# from the completed result tables.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
out <- file.path(gwas, "SCENT_style", "benchmark")
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
source_results <- source_results_paths(root)
dir.create(file.path(out, "config"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

method_colors <- c(
  "EpiLink-Poisson" = "#D65F5F",
  "EpiLink-NB" = "#E39C37",
  "scMultiMap" = "#009E73",
  "Signac (corrected)" = "#7F7F7F",
  "Signac (uncorrected)" = "#2b8cbe"
)
methods <- names(method_colors)
method_levels <- methods

canonical_method <- function(x) {
  y <- trimws(as.character(x))
  y[grepl("^EpiLink-NB", y, ignore.case = TRUE)] <- "EpiLink-NB"
  y[grepl("^EpiLink-Poisson", y, ignore.case = TRUE)] <- "EpiLink-Poisson"
  y[grepl("^scMultiMap$", y, ignore.case = TRUE)] <- "scMultiMap"
  # Match uncorrected first and exclude it from the corrected replacement:
  # the string "uncorrected" contains "corrected".
  is_uncorrected <- grepl("^Signac.*uncorrected", y, ignore.case = TRUE)
  y[is_uncorrected] <- "Signac (uncorrected)"
  y[!is_uncorrected & grepl("^Signac.*corrected", y, ignore.case = TRUE)] <- "Signac (corrected)"
  y
}

read_required <- function(path) {
  if (!file.exists(path)) stop("Required result table is missing: ", path)
  fread(path)
}

clean_scent <- function(path, analysis, view, source_label) {
  x <- read_required(path)
  x[, method_label := canonical_method(method_label)]
  x <- x[method_label %in% methods]
  if (!nrow(x)) stop("No comparison methods found in ", path)
  x[, `:=`(
    analysis = analysis,
    view = view,
    source = source_label,
    celltype = as.character(celltype),
    method_label = factor(method_label, levels = method_levels)
  )]
  x[, metric := as.numeric(mean_enrichment)]
  x[, metric_log2 := log2(1 + pmax(metric, 0))]
  x[, source_file := path]
  x[]
}

clean_ukbb <- function(path, view, source_label) {
  # Sakaue-style summaries call this field celltype_scope rather than
  # celltype; normalize it before using the common collector.
  raw <- read_required(path)
  raw[, celltype := as.character(celltype_scope)]
  tmp <- tempfile(fileext = ".tsv")
  fwrite(raw, tmp, sep = "\t")
  on.exit(unlink(tmp), add = TRUE)
  x <- clean_scent(tmp, "SCENT (UKBB 94 traits)", view, source_label)
  # The older custom table uses ALL_REQUESTED for pooled annotations.  Keep a
  # readable cell-type label in the unified table.
  x[, celltype := ifelse(celltype == "ALL_REQUESTED", "ALL", celltype)]
  x
}

clean_sldsc <- function(path, view, source_label) {
  x <- read_required(path)
  x[, method_label := canonical_method(method_label)]
  x <- x[method_label %in% methods & is.finite(Enrichment) & Enrichment > 0]
  if (!nrow(x)) stop("No usable sLDSC comparison methods found in ", path)
  # sLDSC has one row per trait.  The summary metric is the arithmetic mean of
  # log2(1 + enrichment) across the matched traits; the raw trait rows remain
  # available in the original sLDSC result table.
  z <- x[, .(
    n_traits = .N,
    mean_enrichment = mean(Enrichment),
    median_enrichment = median(Enrichment),
    mean_log2_enrichment = mean(log2(1 + Enrichment)),
    n_nominal_p05 = sum(is.finite(Enrichment_p) & Enrichment_p < 0.05),
    n_fdr05 = sum(is.finite(enrichment_fdr) & enrichment_fdr < 0.05),
    mean_enrichment_p = mean(Enrichment_p, na.rm = TRUE),
    source_file = path
  ), by = .(dataset, celltype, panel, interpretation, method, method_label)]
  z[, `:=`(
    analysis = "sLDSC (FinnGen autoimmune traits)",
    view = view,
    source = source_label,
    metric = mean_enrichment,
    metric_log2 = mean_log2_enrichment
  )]
  z[, method_label := factor(method_label, levels = method_levels)]
  z[]
}

# Current primary analysis: restricted_v2 is the latest result set with the
# requested dataset/cell-type scope.  The pooled sLDSC run was completed in
# the sLDSC/pooled directory; its trait-level table is used directly.
scent_cell <- clean_scent(
  file.path(gwas, "SCENT_style", "testing", "restricted_v2", "celltype", "SCENT", "results", "scent_panel_summary.tsv"),
  "SCENT (FinnGen autoimmune traits)", "celltype", "FinnGen R12"
)
scent_pool <- clean_scent(
  file.path(gwas, "SCENT_style", "testing", "restricted_v2", "pooled", "SCENT", "results", "scent_panel_summary.tsv"),
  "SCENT (FinnGen autoimmune traits)", "pooled", "FinnGen R12"
)
sldsc_cell <- clean_sldsc(
  file.path(gwas, "sLDSC", "testing", "restricted_v2", "celltype", "sLDSC", "results", "sldsc_all_results.tsv"),
  "celltype", "FinnGen R12"
)
sldsc_pool <- clean_sldsc(
  file.path(gwas, "sLDSC", "pooled", "results", "sldsc_all_results.tsv"),
  "pooled", "FinnGen R12"
)

# Sakaue-style UKBB analysis: pooled and cell-type source results already
# include the Brain2/Pituitary extension.
ukbb_cell <- clean_ukbb(
  file.path(source_results$celltype_results, "scent_sakaue_summary.tsv"),
  "celltype", "UKBB 94 traits: merged source results"
)
ukbb_pool <- clean_ukbb(
  file.path(source_results$pooled_results, "scent_sakaue_summary.tsv"),
  "pooled", "UKBB 94 traits: merged source results"
)
# Current figure scope: do not report NEAT or Brain2_CTRL in the UKBB
# comparison.  The original Brain dataset remains available.
ukbb_cell <- ukbb_cell[dataset != "NEAT"]
ukbb_pool <- ukbb_pool[dataset != "NEAT"]
ukbb_cell <- ukbb_cell[dataset != "Brain2_CTRL"]
ukbb_pool <- ukbb_pool[dataset != "Brain2_CTRL"]

write_table <- function(x, path) {
  y <- copy(x)
  if ("method_label" %in% names(y)) y[, method_label := as.character(method_label)]
  fwrite(y, path, sep = "\t", na = "NA")
}

write_table(scent_pool, file.path(out, "tables", "scent_pooled.tsv"))
write_table(scent_cell, file.path(out, "tables", "scent_celltype.tsv"))
write_table(sldsc_pool, file.path(out, "tables", "sldsc_pooled.tsv"))
write_table(sldsc_cell, file.path(out, "tables", "sldsc_celltype.tsv"))
write_table(ukbb_pool, file.path(out, "tables", "scent_ukbb_pooled.tsv"))
write_table(ukbb_cell, file.path(out, "tables", "scent_ukbb_celltype.tsv"))

# The master tables use one common metric definition.  For sLDSC the metric is
# a trait-level summary; for SCENT it is the trait-bootstrap mean enrichment.
master_pool <- rbindlist(list(
  scent_pool[, .(analysis, view, source, dataset, celltype, method_label,
                 metric, metric_log2, n_traits, source_file)],
  sldsc_pool[, .(analysis, view, source, dataset, celltype, method_label,
                 metric, metric_log2, n_traits, source_file)],
  ukbb_pool[, .(analysis, view, source, dataset, celltype, method_label,
                metric, metric_log2, n_traits, source_file)]
), fill = TRUE)
master_cell <- rbindlist(list(
  scent_cell[, .(analysis, view, source, dataset, celltype, method_label,
                 metric, metric_log2, n_traits, source_file)],
  sldsc_cell[, .(analysis, view, source, dataset, celltype, method_label,
                 metric, metric_log2, n_traits, source_file)],
  ukbb_cell[, .(analysis, view, source, dataset, celltype, method_label,
                metric, metric_log2, n_traits, source_file)]
), fill = TRUE)
write_table(master_pool, file.path(out, "tables", "gwas_master_pooled.tsv"))
write_table(master_cell, file.path(out, "tables", "gwas_master_celltype.tsv"))

manifest <- data.table(
  analysis = c("SCENT (FinnGen autoimmune traits)", "sLDSC (FinnGen autoimmune traits)",
               "SCENT (UKBB 94 traits)"),
  pooled_source = c(
    file.path(gwas, "SCENT_style", "testing", "restricted_v2", "pooled", "SCENT", "results", "scent_panel_summary.tsv"),
    file.path(gwas, "sLDSC", "pooled", "results", "sldsc_all_results.tsv"),
    file.path(source_results$pooled_results, "scent_sakaue_summary.tsv")
  ),
  celltype_source = c(
    file.path(gwas, "SCENT_style", "testing", "restricted_v2", "celltype", "SCENT", "results", "scent_panel_summary.tsv"),
    file.path(gwas, "sLDSC", "testing", "restricted_v2", "celltype", "sLDSC", "results", "sldsc_all_results.tsv"),
    file.path(source_results$celltype_results, "scent_sakaue_summary.tsv")
  ),
  metric = c("trait-bootstrap mean enrichment", "mean across matched traits of log2(1 + enrichment)",
             "trait-bootstrap mean enrichment"),
  methods = paste(methods, collapse = "; ")
)
fwrite(manifest, file.path(out, "config", "analysis_manifest.tsv"), sep = "\t")
fwrite(data.table(method_order = seq_along(methods), method_label = methods),
       file.path(out, "config", "method_manifest.tsv"), sep = "\t")

theme_benchmark <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

plot_heatmap <- function(x, file, title, subtitle) {
  if (!nrow(x)) return(invisible(NULL))
  x <- copy(x)
  x[, method_label := factor(as.character(method_label), levels = method_levels)]
  x[, row_label := ifelse(view == "pooled" | celltype %in% c("ALL", "ALL_REQUESTED"),
                          as.character(dataset), paste(dataset, celltype, sep = " | "))]
  x[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  x[, analysis := factor(analysis, levels = c("SCENT (FinnGen autoimmune traits)",
                                               "sLDSC (FinnGen autoimmune traits)",
                                               "SCENT (UKBB 94 traits)"))]
  p <- ggplot(x, aes(method_label, row_label, fill = metric_log2)) +
    geom_tile(color = "white", linewidth = 0.25, na.rm = TRUE) +
    geom_point(aes(color = method_label), shape = 16, size = 1.8, na.rm = TRUE) +
    geom_text(aes(label = ifelse(is.finite(metric), sprintf("%.1f", metric), "NA")),
              size = 2.5, na.rm = TRUE) +
    facet_wrap(~analysis, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         name = "log2(1 + enrichment)", na.value = "grey85") +
    scale_color_manual(values = method_colors, drop = FALSE, name = "Method") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) + theme_benchmark
  h <- max(6.2, 0.23 * length(unique(x$row_label)) * length(unique(x$analysis)) + 3)
  ggsave(file, p, width = 10.5, height = h, units = "in", device = cairo_pdf)
}

plot_heatmap(master_pool,
             file.path(out, "figures", "gwas_benchmark_pooled.pdf"),
             "GWAS validation benchmark: pooled cell types",
             "Each tile is the method-level enrichment summary; labels show the raw enrichment")
plot_heatmap(master_cell,
             file.path(out, "figures", "gwas_benchmark_celltype.pdf"),
             "GWAS validation benchmark: cell-type-specific",
             "Each tile is evaluated within one dataset × cell type; labels show the raw enrichment")

# UKBB-only figures retain the bootstrap intervals and are easier to use for
# the Sakaue-style comparison than the compact master heatmaps.
plot_ukbb_forest <- function(x, file, title, view_label) {
  x <- copy(x)
  x <- x[is.finite(mean_enrichment) & is.finite(ci_low) & is.finite(ci_high)]
  if (!nrow(x)) return(invisible(NULL))
  x[, method_label := factor(as.character(method_label), levels = method_levels)]
  x[, row_label := ifelse(view == "pooled" | celltype == "ALL", as.character(dataset),
                          paste(dataset, celltype, sep = " | "))]
  x[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  p <- ggplot(x, aes(mean_enrichment, row_label, color = method_label)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50", linewidth = 0.35) +
    geom_errorbar(aes(y = row_label, xmin = pmax(ci_low, 0), xmax = ci_high),
                  orientation = "y", height = 0.18,
                  position = position_dodge(width = 0.6), linewidth = 0.45,
                  na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.6), size = 1.9, na.rm = TRUE) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 20, 50, 100),
                       labels = c("0", "1", "2", "5", "10", "20", "50", "100"),
                       expand = expansion(mult = c(0.02, 0.08))) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    labs(title = title, subtitle = paste0("UKBB 94 traits; ", view_label,
                                          "; point = mean across traits, line = 95% trait-bootstrap CI"),
         x = "Mean causal-variant enrichment (log1p scale)", y = NULL, color = "Method") +
    theme_benchmark
  h <- max(5.5, 0.27 * length(unique(x$row_label)) + 2.3)
  ggsave(file, p, width = 11.5, height = h, units = "in", device = cairo_pdf)
}

plot_ukbb_forest(ukbb_pool, file.path(out, "figures", "scent_ukbb_pooled.pdf"),
                 "SCENT-style UKBB benchmark", "pooled")
plot_ukbb_forest(ukbb_cell, file.path(out, "figures", "scent_ukbb_celltype.pdf"),
                 "SCENT-style UKBB benchmark", "cell-type-specific")

message("Wrote unified GWAS benchmark to: ", out)
message("  pooled rows: ", nrow(master_pool), "; celltype rows: ", nrow(master_cell))
