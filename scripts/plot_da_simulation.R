#!/usr/bin/env Rscript

file_argument <- sub(
  "^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
)
project_root <- normalizePath(file.path(dirname(file_argument), ".."),
                              mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("yaml", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("yaml and ggplot2 are required", call. = FALSE)
}
for (source_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  source(source_file)
}

settings <- yaml::read_yaml("config/simulation_da.yml")
summary_file <- file.path(settings$output_dir, "summary_metrics.tsv")
if (!file.exists(summary_file)) {
  stop("Missing ", summary_file, "; run run_da_simulation.R first",
       call. = FALSE)
}

summary <- utils::read.delim(summary_file, stringsAsFactors = FALSE)
required <- c("dgp", "n_cells", "effect_size", "method_id", "metric",
              "mean", "se")
missing <- setdiff(required, names(summary))
if (length(missing) > 0L) {
  stop("summary_metrics.tsv is missing: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

method_order <- c(
  "glm_poisson", "glm_nb", "epilink_poisson", "epilink_nb", "scmultimap"
)
method_labels <- c(
  glm_poisson = "Poisson-GLM", glm_nb = "NB-GLM",
  epilink_poisson = "EpiLink-Poisson", epilink_nb = "EpiLink-NB",
  scmultimap = "scMultiMap"
)
method_colors <- c(
  glm_poisson = "#9ecae1", glm_nb = "#fa9fb5",
  epilink_poisson = "#08519c", epilink_nb = "#ce1256",
  scmultimap = "#E39C37"
)

summary$method_id <- factor(summary$method_id, levels = method_order)
summary$dgp <- factor(
  tolower(summary$dgp), levels = c("poisson", "nb", "zinb"),
  labels = c("Poisson", "NB", "ZINB")
)
summary$n_cells <- factor(
  summary$n_cells, levels = sort(unique(summary$n_cells)), ordered = TRUE
)
summary$metric <- as.character(summary$metric)
summary$mean <- as.numeric(summary$mean)
summary$se <- as.numeric(summary$se)
summary$effect_size <- as.numeric(summary$effect_size)

figure_dir <- file.path(settings$output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

make_plot <- function(metric_name, title, output_file, add_target = FALSE) {
  dat <- summary[summary$metric == metric_name &
                   is.finite(summary$mean), , drop = FALSE]
  if (nrow(dat) == 0L) {
    warning("No rows available for ", metric_name, "; skipping plot")
    return(invisible(NULL))
  }

  # Keep lines within each cell-size/DGP/method facet.
  dat$line_group <- interaction(dat$n_cells, dat$dgp, dat$method_id,
                                drop = TRUE)
  line_dat <- dat[
    ave(dat$mean, dat$line_group, FUN = length) >= 2L, , drop = FALSE
  ]

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(effect_size, mean, color = method_id, group = line_group)
  )
  if (add_target) {
    target_dat <- expand.grid(
      n_cells = levels(summary$n_cells),
      dgp = levels(summary$dgp),
      target = as.numeric(settings$fdr),
      stringsAsFactors = FALSE
    )
    target_dat$n_cells <- factor(target_dat$n_cells,
                                 levels = levels(summary$n_cells), ordered = TRUE)
    target_dat$dgp <- factor(target_dat$dgp, levels = levels(summary$dgp))
    p <- p + ggplot2::geom_hline(
      data = target_dat, ggplot2::aes(yintercept = target),
      inherit.aes = FALSE, linetype = "dashed", linewidth = 0.35,
      color = "grey45"
    )
  }
  p <- p +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(mean - se, 0), ymax = pmin(mean + se, 1)),
      width = 0.02, linewidth = 0.3, alpha = 0.5, na.rm = TRUE
    ) +
    ggplot2::geom_line(data = line_dat, linewidth = 0.7, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.7, na.rm = TRUE) +
    ggplot2::facet_grid(n_cells ~ dgp, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = method_colors, breaks = method_order,
      labels = unname(method_labels[method_order]), drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(dat$effect_size)),
      expand = ggplot2::expansion(mult = c(0.04, 0.06))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.03, 0.05))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Differential peak-gene association simulation",
      x = "Effect size", y = metric_name, color = NULL
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    file.path(figure_dir, output_file), p,
    width = 9, height = 8, units = "in"
  )
  message("Wrote ", file.path(figure_dir, output_file))
  invisible(p)
}

make_plot("FDR", "Differential association FDR",
          "fdr_by_effect_size_and_cell_number.pdf", add_target = TRUE)
make_plot("Power", "Differential association power",
          "power_by_effect_size_and_cell_number.pdf", add_target = FALSE)
