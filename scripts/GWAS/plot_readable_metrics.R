#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")

scent_results <- Sys.getenv("SCENT_RESULTS_DIR", file.path(gwas, "SCENT_style", "testing", "results"))
sldsc_results <- Sys.getenv("SLDSC_RESULTS_DIR", file.path(gwas, "sLDSC", "results"))
scent_metric_out <- Sys.getenv("SCENT_METRIC_OUT_DIR", file.path(gwas, "SCENT_style", "UKBB", "metrics"))
sldsc_metric_out <- Sys.getenv("SLDSC_METRIC_OUT_DIR", file.path(gwas, "sLDSC", "metrics"))
scent_fig <- Sys.getenv("SCENT_READABLE_FIGURE_DIR", file.path(gwas, "SCENT_style", "testing", "figures", "readable"))
sldsc_fig <- Sys.getenv("SLDSC_READABLE_FIGURE_DIR", file.path(gwas, "sLDSC", "figures", "readable"))
dir.create(scent_metric_out, recursive = TRUE, showWarnings = FALSE)
dir.create(sldsc_metric_out, recursive = TRUE, showWarnings = FALSE)
dir.create(scent_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(sldsc_fig, recursive = TRUE, showWarnings = FALSE)

methods <- fread(file.path(gwas, "SCENT_style", "config", "methods.tsv"))
method_order <- methods$label
exclude_labels <- trimws(strsplit(Sys.getenv("GWAS_PLOT_EXCLUDE_METHODS", "NB-GLM,Poisson-GLM"), ",", fixed = TRUE)[[1]])
method_order <- setdiff(method_order, exclude_labels)

safe_mean <- function(x) if (length(x) && any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
safe_median <- function(x) if (length(x) && any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
safe_rate <- function(x) if (length(x) && any(!is.na(x))) mean(x, na.rm = TRUE) else NA_real_
safe_file <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

theme_readable <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey94", color = "grey40"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

save_both <- function(p, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, limitsize = FALSE)
}

## -------------------------------
## SCENT metrics and plots
## -------------------------------
scent_file <- file.path(scent_results, "scent_panel_summary.tsv")
if (file.exists(scent_file)) {
  scent <- fread(scent_file)
  scent <- scent[!method_label %in% exclude_labels]
  ## Keep true zero-overlap annotations. They are valid SCENT results, not missing data.
  scent <- scent[interpretation == "primary" & is.finite(mean_enrichment) & mean_enrichment >= 0]
  scent[, `:=`(
    method_label = factor(method_label, levels = method_order),
    log2_enrichment = log2(1 + mean_enrichment),
    ci_excludes_1 = (ci_low > 1 | ci_high < 1)
  )]

  scent_method <- scent[, .(
    n_annotations = .N,
    mean_enrichment = safe_mean(mean_enrichment),
    median_enrichment = safe_median(mean_enrichment),
    median_log2_enrichment = safe_median(log2_enrichment),
    fraction_enrichment_gt1 = safe_rate(mean_enrichment > 1),
    n_ci_excludes_1 = sum(ci_excludes_1, na.rm = TRUE),
    fraction_ci_excludes_1 = safe_rate(ci_excludes_1)
  ), by = method_label]
  scent_method <- merge(data.table(method_label = method_order), scent_method,
                        by = "method_label", all.x = TRUE, sort = FALSE)
  fwrite(scent_method, file.path(scent_metric_out, "scent_primary_metrics_by_method.tsv"), sep = "\t", na = "NA")

  scent_detail <- scent[, .(
    n_traits = first(n_traits),
    mean_enrichment = first(mean_enrichment),
    median_enrichment = first(median_enrichment),
    ci_low = first(ci_low),
    ci_high = first(ci_high),
    pooled_causal_in = first(pooled_causal_in),
    pooled_causal_total = first(pooled_causal_total)
  ), by = .(dataset, celltype, method_label)]
  fwrite(scent_detail, file.path(scent_metric_out, "scent_primary_metrics_by_dataset_celltype.tsv"), sep = "\t", na = "NA")

  ## One readable panel per dataset; each dataset has its own y scale.
  p_scent <- ggplot(scent, aes(method_label, log2_enrichment, color = celltype)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
    geom_point(position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.68), size = 2.4, alpha = 0.9) +
    facet_wrap(~dataset, scales = "free_y", ncol = 3) +
    labs(
      title = "SCENT: primary immune GWAS enrichment by method",
      subtitle = "Y-axis is log2(1 + mean enrichment); true zero-overlap annotations are retained; each dataset has an independent y-axis",
      x = NULL, y = "log2(1 + mean enrichment)", color = "Cell type"
    ) + theme_readable
  save_both(p_scent, file.path(scent_fig, "scent_primary_method_comparison_log2"), 16, 10)
}

## -------------------------------
## sLDSC metrics and plots
## -------------------------------
sldsc_file <- file.path(sldsc_results, "sldsc_all_results.tsv")
if (file.exists(sldsc_file)) {
  sldsc <- fread(sldsc_file)
  sldsc <- sldsc[!method_label %in% exclude_labels]
  sldsc <- sldsc[interpretation == "primary" & is.finite(Enrichment) & Enrichment > 0]
  sldsc[, `:=`(
    method_label = factor(method_label, levels = method_order),
    log2_enrichment = log2(Enrichment),
    nominal_p05 = Enrichment_p < 0.05,
    fdr05 = enrichment_fdr < 0.05
  )]

  sldsc_method <- sldsc[, .(
    n_trait_annotations = .N,
    mean_enrichment = safe_mean(Enrichment),
    median_enrichment = safe_median(Enrichment),
    median_log2_enrichment = safe_median(log2_enrichment),
    fraction_enrichment_gt1 = safe_rate(Enrichment > 1),
    n_nominal_p05 = sum(nominal_p05, na.rm = TRUE),
    fraction_nominal_p05 = safe_rate(nominal_p05),
    n_fdr05 = sum(fdr05, na.rm = TRUE),
    fraction_fdr05 = safe_rate(fdr05)
  ), by = method_label]
  sldsc_method <- merge(data.table(method_label = method_order), sldsc_method,
                        by = "method_label", all.x = TRUE, sort = FALSE)
  fwrite(sldsc_method, file.path(sldsc_metric_out, "sldsc_primary_metrics_by_method.tsv"), sep = "\t", na = "NA")

  sldsc_detail <- sldsc[, .(
    Enrichment = first(Enrichment),
    Enrichment_std_error = first(Enrichment_std_error),
    Enrichment_p = first(Enrichment_p),
    enrichment_fdr = first(enrichment_fdr),
    log2_enrichment = first(log2_enrichment)
  ), by = .(trait, trait_label, dataset, celltype, method_label)]
  fwrite(sldsc_detail, file.path(sldsc_metric_out, "sldsc_primary_metrics_by_trait_dataset_celltype.tsv"), sep = "\t", na = "NA")

  ## One plot per trait keeps the y axis tall enough for direct method comparison.
  for (tr in unique(sldsc$trait)) {
    z <- sldsc[trait == tr]
    p <- ggplot(z, aes(method_label, log2_enrichment, color = celltype, size = -log10(pmax(Enrichment_p, 1e-300)))) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
      geom_point(position = position_jitter(width = 0.07, height = 0), alpha = 0.9) +
      facet_wrap(~dataset, scales = "free_y", ncol = 3) +
      scale_size_continuous(range = c(1.8, 4.5), name = expression(-log[10](P))) +
      labs(
        title = paste0("sLDSC: ", z$trait_label[1], " primary GWAS enrichment"),
        subtitle = "Y-axis is log2(heritability enrichment); each dataset has an independent y-axis",
        x = NULL, y = "log2(heritability enrichment)", color = "Cell type"
    ) + theme_readable
    save_both(p, file.path(sldsc_fig, paste0("sldsc_primary_", safe_file(tr), "_method_comparison_log2")), 16, 9)
  }
}

message("Wrote readable SCENT metric tables to ", scent_metric_out)
message("Wrote readable sLDSC metric tables to ", sldsc_metric_out)
message("Wrote readable SCENT figures to ", scent_fig)
message("Wrote readable sLDSC figures to ", sldsc_fig)
