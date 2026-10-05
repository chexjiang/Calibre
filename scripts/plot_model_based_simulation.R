#!/usr/bin/env Rscript

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
if (!requireNamespace("yaml", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("yaml and ggplot2 are required", call. = FALSE)
}
configuration <- yaml::read_yaml("config/simulation.yml")
settings <- utils::modifyList(configuration$defaults, configuration$profiles$model_based)
output_dir <- settings$output_dir
metrics <- readRDS(file.path(output_dir, "replicate_metrics.rds"))
scent_file <- file.path(output_dir, "scent_replicate_metrics.rds")
if (file.exists(scent_file) && !any(metrics$method_id == "scent")) {
  metrics <- rbind(metrics, readRDS(scent_file))
}
archr_file <- file.path(output_dir, "archr_replicate_metrics.rds")
if (file.exists(archr_file) && !any(metrics$method_id == "archr")) {
  metrics <- rbind(metrics, readRDS(archr_file))
}
reference_file <- settings$reference_500_file
if (!is.null(reference_file) && file.exists(reference_file)) {
  reference <- readRDS(reference_file)
  scent_reference_file <- file.path(output_dir, "scent_reference_500_metrics.rds")
  if (file.exists(scent_reference_file) && !any(reference$method_id == "scent")) {
    reference <- rbind(reference, readRDS(scent_reference_file))
  }
  metrics <- rbind(metrics, reference)
}
metrics <- metrics[order(metrics$n_cells, metrics$dgp, metrics$effect_size,
                         metrics$replicate, metrics$method_id), ]
saveRDS(metrics, file.path(output_dir, "replicate_metrics_cells_100_500.rds"))
utils::write.table(
  metrics, gzfile(file.path(output_dir, "replicate_metrics_cells_100_500.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

long <- rbind(
  transform(metrics, metric = "FDR", value = fdr),
  transform(metrics, metric = "Power", value = power)
)
group <- interaction(long$n_cells, long$dgp, long$effect_size,
                     long$method_id, long$metric, drop = TRUE)
summary_rows <- lapply(split(long, group), function(x) {
  values <- x$value[is.finite(x$value)]
  data.frame(
    n_cells = x$n_cells[1], dgp = x$dgp[1], effect_size = x$effect_size[1],
    method = x$method[1], method_id = x$method_id[1], metric = x$metric[1],
    mean = if (length(values)) mean(values) else NA_real_,
    se = if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_,
    n = length(values)
  )
})
summary_metrics <- do.call(rbind, summary_rows)
summary_metrics$dgp <- factor(
  summary_metrics$dgp, levels = c("poisson", "nb", "zinb"),
  labels = c("Poisson DGP", "NB DGP", "ZINB DGP")
)
summary_metrics$metric <- factor(summary_metrics$metric, levels = c("FDR", "Power"))
utils::write.table(
  summary_metrics, file.path(output_dir, "summary_metrics_cells_100_500.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

runtime_group <- interaction(metrics$n_cells, metrics$method_id, drop = TRUE)
runtime_rows <- lapply(split(metrics, runtime_group), function(x) {
  values <- x$dataset_runtime_seconds[is.finite(x$dataset_runtime_seconds)]
  data.frame(
    n_cells = x$n_cells[1], method = x$method[1], method_id = x$method_id[1],
    mean_seconds = if (length(values)) mean(values) else NA_real_,
    se_seconds = if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_,
    n = length(values), stringsAsFactors = FALSE
  )
})
runtime_summary <- do.call(rbind, runtime_rows)
runtime_summary <- runtime_summary[
  order(runtime_summary$n_cells, runtime_summary$method_id), , drop = FALSE
]
utils::write.table(
  runtime_summary, file.path(output_dir, "runtime_summary_cells_100_400.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

palette <- c(
  "glm_poisson_bh" = "#9ecae1", "glm_nb_bh" = "#fa9fb5",
  "epilink_poisson_residual" = "#08519c",
  "epilink_nb_residual" = "#ce1256", "scmultimap" = "#E39C37",
  "signac" = "#525252", "scent" = "#41ae76", "archr" = "#756bb1"
)
method_lookup <- unique(summary_metrics[, c("method_id", "method")])
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

make_plot <- function(cell_count) {
  plot_data <- summary_metrics[summary_metrics$n_cells == cell_count &
                                 is.finite(summary_metrics$mean), , drop = FALSE]
  error_data <- plot_data[is.finite(plot_data$se), , drop = FALSE]
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(effect_size, mean, color = method_id, group = method_id)
  ) +
    ggplot2::geom_hline(
      data = data.frame(metric = factor("FDR", levels = c("FDR", "Power")),
                        y = settings$fdr),
      ggplot2::aes(yintercept = y), linetype = "dashed", color = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_errorbar(
      data = error_data,
      ggplot2::aes(ymin = pmax(mean - se, 0), ymax = pmin(mean + se, 1)),
      width = 0.02, linewidth = 0.3, alpha = 0.55
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_grid(metric ~ dgp) +
    ggplot2::scale_color_manual(
      values = palette, breaks = method_lookup$method_id,
      labels = method_lookup$method
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title = paste0("Simulation benchmark: ", cell_count, " cells"),
      subtitle = paste0(settings$n_pairs, " pairs, ", settings$n_replicates,
                        " replicates; target FDR = ", settings$fdr),
      x = "Enhancer effect size (log scale coefficient)", y = NULL, color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "grey45", fill = NA),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey45"),
      legend.position = "bottom", legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
}

for (cell_count in sort(unique(summary_metrics$n_cells))) {
  plot <- make_plot(cell_count)
  stem <- file.path(figure_dir, paste0("fdr_power_cells_", cell_count))
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = 11, height = 6.5)
}

selected_effects <- c(0.2, 0.4, 0.6)
make_cell_count_plot <- function(metric_name) {
  plot_data <- summary_metrics[
    summary_metrics$metric == metric_name &
      summary_metrics$effect_size %in% selected_effects &
      is.finite(summary_metrics$mean), , drop = FALSE
  ]
  plot_data$effect_label <- factor(
    paste0("Effect = ", plot_data$effect_size),
    levels = paste0("Effect = ", selected_effects)
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(n_cells, mean, color = method_id, group = method_id)
  ) +
    {if (metric_name == "FDR") ggplot2::geom_hline(
      yintercept = settings$fdr, linetype = "dashed", color = "black"
    )} +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(mean - se, 0), ymax = pmin(mean + se, 1)),
      width = 8, linewidth = 0.3, alpha = 0.5
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_grid(effect_label ~ dgp) +
    ggplot2::scale_x_continuous(breaks = sort(unique(plot_data$n_cells))) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_color_manual(
      values = palette, breaks = method_lookup$method_id,
      labels = method_lookup$method
    ) +
    ggplot2::labs(
      title = paste0(metric_name, " by cell number and enhancer effect size"),
      subtitle = paste0(settings$n_pairs, " pairs, ", settings$n_replicates,
                        " replicates; target FDR = ", settings$fdr),
      x = "Number of cells", y = metric_name, color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "grey45", fill = NA),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey45"),
      legend.position = "bottom", legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
}

for (metric_name in c("Power", "FDR")) {
  plot <- make_cell_count_plot(metric_name)
  stem <- file.path(figure_dir, paste0(tolower(metric_name), "_by_cell_size"))
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = 11, height = 9)
}

make_effect_size_plot <- function(metric_name) {
  plot_data <- summary_metrics[
    summary_metrics$metric == metric_name & is.finite(summary_metrics$mean),
    , drop = FALSE
  ]
  error_data <- plot_data[is.finite(plot_data$se), , drop = FALSE]
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(effect_size, mean, color = method_id, group = method_id)
  ) +
    {if (metric_name == "FDR") ggplot2::geom_hline(
      yintercept = settings$fdr, linetype = "dashed", color = "black"
    )} +
    ggplot2::geom_errorbar(
      data = error_data,
      ggplot2::aes(ymin = pmax(mean - se, 0), ymax = pmin(mean + se, 1)),
      width = 0.02, linewidth = 0.25, alpha = 0.45
    ) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 1.25) +
    ggplot2::facet_grid(n_cells ~ dgp) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_color_manual(
      values = palette, breaks = method_lookup$method_id,
      labels = method_lookup$method
    ) +
    ggplot2::labs(
      title = paste0(metric_name, " across enhancer effect size and cell number"),
      subtitle = paste0(settings$n_pairs, " pairs, ", settings$n_replicates,
                        " replicates; target FDR = ", settings$fdr),
      x = "Enhancer effect size (log scale coefficient)", y = metric_name,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "grey45", fill = NA),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey45"),
      legend.position = "bottom", legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
}

for (metric_name in c("Power", "FDR")) {
  plot <- make_effect_size_plot(metric_name)
  stem <- file.path(figure_dir, paste0(tolower(metric_name), "_by_effect_size_and_cell_size"))
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = 11, height = 13)
}
message("Sample-size summaries and plots written to ", output_dir)
