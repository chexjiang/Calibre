#!/usr/bin/env Rscript

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
if (!requireNamespace("yaml", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("yaml and ggplot2 are required", call. = FALSE)
}
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(source_file)

settings <- yaml::read_yaml("config/empirical_null.yml")
output_dir <- settings$output_dir
summary_file <- file.path(output_dir, "method_summary.tsv")
replicate_file <- file.path(output_dir, "replicate_metrics.tsv")
if (!file.exists(summary_file) || !file.exists(replicate_file)) {
  stop("Run scripts/summarize_empirical_null.R before plotting", call. = FALSE)
}
method_summary <- utils::read.delim(summary_file, stringsAsFactors = FALSE)
replicate_metrics <- utils::read.delim(replicate_file, stringsAsFactors = FALSE)
registry <- empirical_null_methods()

# The complete-null benchmark has no true positives: every selected pair is a
# false positive. These are the only two requested summary quantities.
complete_replicates <- replicate_metrics[
  replicate_metrics$status == "complete", , drop = FALSE
]
replicate_counts <- stats::aggregate(
  operational_any_rejection ~ method_id, data = complete_replicates,
  FUN = function(x) sum(as.logical(x), na.rm = TRUE)
)
names(replicate_counts)[names(replicate_counts) == "operational_any_rejection"] <-
  "null_replicates_with_discovery"
discovery_summary <- merge(
  method_summary[, c("method_id", "method", "n_complete_replicates",
                     "mean_operational_false_discoveries")],
  replicate_counts, by = "method_id", all.x = TRUE, sort = FALSE
)
discovery_summary <- discovery_summary[match(registry$method_id,
                                             discovery_summary$method_id), ]
utils::write.table(
  discovery_summary,
  file.path(output_dir, "false_positive_discovery_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

discovery_data <- rbind(
  data.frame(
    method_id = discovery_summary$method_id,
    method = discovery_summary$method,
    metric = "Mean false-positive discoveries per null replicate",
    value = discovery_summary$mean_operational_false_discoveries,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method_id = discovery_summary$method_id,
    method = discovery_summary$method,
    metric = "Null replicates with >=1 significant pair",
    value = discovery_summary$null_replicates_with_discovery,
    stringsAsFactors = FALSE
  )
)
discovery_data$method <- factor(discovery_data$method, levels = registry$method)
discovery_data$metric <- factor(
  discovery_data$metric,
  levels = c("Mean false-positive discoveries per null replicate",
             "Null replicates with >=1 significant pair")
)

palette <- c(
  "EpiLink-Poisson" = "#D65F5F", "EpiLink-NB" = "#8E6BBE",
  "Poisson-GLM" = "#2E6F9E", "NB-GLM" = "#4E9F6E",
  "scMultiMap" = "#E39C37", "Signac" = "#8A8A8A", "SCENT" = "#333333"
)
p_discoveries <- ggplot2::ggplot(
  discovery_data,
  ggplot2::aes(x = method, y = value, fill = method)
) +
  ggplot2::geom_col(width = 0.68, na.rm = TRUE) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(is.finite(value), format(round(value, 2),
                                                       trim = TRUE), "NA")),
    vjust = -0.3, size = 3.2, na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
  ggplot2::labs(
    title = "False-positive discoveries in PBMC/Tnk empirical-null data",
    subtitle = paste0(
      settings$n_replicates, " null replicates; all selected pairs are false positives; ",
      "final method-specific selection at ", settings$alpha
    ),
    x = NULL, y = NULL, fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.border = ggplot2::element_rect(
      color = "grey45", fill = NA, linewidth = 0.6
    ),
    strip.background = ggplot2::element_rect(
      fill = "grey92", color = "grey45"
    ),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  file.path(figure_dir, "false_positive_discoveries.pdf"), p_discoveries,
  width = 12, height = 5.8
)
message("False-positive discovery figure written to ", figure_dir)
