# Shared odds-ratio, confidence-interval, and p-value plots for validation.

validation_method_levels <- c(
  "Calibre-Poisson", "Calibre-NB", "scMultiMap", "scMultiMap-noBatch",
  "Signac (corrected)", "Signac (uncorrected)", "ArchR"
)

validation_palette <- c(
  "Calibre-Poisson" = "#F94144", "Calibre-NB" = "#E7A2C2",
  "scMultiMap" = "#54b2a2", 
  "Signac (corrected)" = "#08519c", "Signac (uncorrected)" = "grey",
  "ArchR" = "#F8961E"
)

validation_celltypes <- c("Bcell", "Myeloid", "Tnk", "All")
validation_cutoffs <- c(0.05, 0.10, 0.15, 0.20)

fisher_metric_row <- function(template, tp, fp, fn, tn, celltype = "All") {
  table <- matrix(c(tp, fp, fn, tn), nrow = 2L, byrow = TRUE,
                  dimnames = list(selection = c("selected", "not_selected"),
                                  validation = c("positive", "negative")))
  defined <- all(rowSums(table) > 0L) && all(colSums(table) > 0L)
  test <- if (defined) stats::fisher.test(table, alternative = "greater") else NULL
  out <- template[1L, , drop = FALSE]
  out$celltype <- celltype
  out$TP <- tp
  out$FP <- fp
  out$FN <- fn
  out$TN <- tn
  out$odds_ratio <- if (defined) unname(test$estimate) else NA_real_
  out$CI_low <- if (defined) test$conf.int[[1L]] else NA_real_
  out$CI_high <- if (defined) test$conf.int[[1L + 1L]] else NA_real_
  out$pval <- if (defined) test$p.value else NA_real_
  if ("n_universe" %in% names(out)) out$n_universe <- tp + fp + fn + tn
  if ("n_selected" %in% names(out)) out$n_selected <- tp + fp
  if ("n_pchic_links" %in% names(out)) out$n_pchic_links <- tp + fn
  out
}

# PCHiC metrics are stored per cell type. Pool their 2-by-2 tables to make the
# requested All panel without rerunning pair selection or enrichment tests.
append_validation_all <- function(metrics, celltypes = validation_celltypes[1:3]) {
  metrics$cutoff <- as.numeric(metrics$cutoff)
  if ("All" %in% as.character(metrics$celltype)) return(metrics)
  keep <- metrics[as.character(metrics$celltype) %in% celltypes, , drop = FALSE]
  if (!nrow(keep)) return(metrics)
  groups <- split(seq_len(nrow(keep)),
                  interaction(keep$method, keep$cutoff, drop = TRUE, lex.order = TRUE))
  all_rows <- lapply(groups, function(index) {
    x <- keep[index, , drop = FALSE]
    fisher_metric_row(
      x, sum(x$TP), sum(x$FP), sum(x$FN), sum(x$TN), celltype = "All"
    )
  })
  rbind(metrics, do.call(rbind, all_rows))
}

prepare_validation_plot_data <- function(metrics, use_haldane = FALSE) {
  metrics$cutoff <- as.numeric(metrics$cutoff)
  # Plot only the requested validation methods. The GLM-derived Poisson/NB rows
  # remain in enrichment_metrics.tsv but are intentionally excluded from all
  # validation figures.
  metrics <- metrics[as.character(metrics$method) %in% validation_method_levels,
                     , drop = FALSE]
  # An infinite odds ratio is not a finite point estimate. Keep its row so
  # the cutoff/method layout stays aligned, but make the OR and CI missing so
  # ggplot leaves a gap instead of drawing a clipped artificial value.
  infinite_or <- is.infinite(metrics$odds_ratio)
  observed_celltypes <- unique(as.character(metrics$celltype))
  ordered_celltypes <- c(
    validation_celltypes[validation_celltypes != "All" &
                          validation_celltypes %in% observed_celltypes],
    sort(setdiff(observed_celltypes, c(validation_celltypes, "All"))),
    if ("All" %in% observed_celltypes) "All"
  )
  metrics$celltype <- factor(as.character(metrics$celltype), levels = ordered_celltypes)
  metrics$method <- factor(as.character(metrics$method), levels = validation_method_levels)
  finite_values <- c(
    metrics$odds_ratio[is.finite(metrics$odds_ratio) & metrics$odds_ratio > 0],
    metrics$CI_high[is.finite(metrics$CI_high) & metrics$CI_high > 0]
  )
  max_finite <- if (length(finite_values)) max(finite_values, na.rm = TRUE) else 1
  cap <- max(10, max_finite * 1.5)
  metrics$OR_inf <- infinite_or
  if (use_haldane && all(c("TP", "FP", "FN", "TN") %in% names(metrics))) {
    # Recalculate the Haldane--Anscombe display values from the four counts so
    # this helper also works with older metric tables.  Raw Fisher estimates
    # and exact confidence intervals remain separate columns.
    tp <- pmax(as.numeric(metrics$TP), 0)
    fp <- pmax(as.numeric(metrics$FP), 0)
    fn <- pmax(as.numeric(metrics$FN), 0)
    tn <- pmax(as.numeric(metrics$TN), 0)
    metrics$odds_ratio_ha <- (tp + 0.5) * (tn + 0.5) /
      ((fp + 0.5) * (fn + 0.5))
    se_log_or <- sqrt(1 / (tp + 0.5) + 1 / (fp + 0.5) +
                        1 / (fn + 0.5) + 1 / (tn + 0.5))
    metrics$CI_low_ha <- exp(log(metrics$odds_ratio_ha) - 1.96 * se_log_or)
    metrics$CI_high_ha <- exp(log(metrics$odds_ratio_ha) + 1.96 * se_log_or)
    valid_margins <- (tp + fp > 0) & (fn + tn > 0) &
      (tp + fn > 0) & (fp + tn > 0)
    metrics$odds_ratio_ha[!valid_margins] <- NA_real_
    metrics$CI_low_ha[!valid_margins] <- NA_real_
    metrics$CI_high_ha[!valid_margins] <- NA_real_
    metrics$odds_ratio_plot <- pmax(metrics$odds_ratio_ha, 0.05)
    metrics$CI_low_plot <- pmax(metrics$CI_low_ha, 0.05)
    metrics$CI_high_plot <- pmax(metrics$CI_high_ha, 0.05)
  } else {
    metrics$odds_ratio_plot <- ifelse(
      metrics$OR_inf, NA_real_, pmax(metrics$odds_ratio, 0.05)
    )
    metrics$CI_low_plot <- ifelse(metrics$OR_inf, NA_real_,
                                  pmax(metrics$CI_low, 0.05))
    metrics$CI_high_plot <- ifelse(metrics$OR_inf, NA_real_,
                                   ifelse(is.infinite(metrics$CI_high), cap,
                                          pmax(metrics$CI_high, 0.05)))
  }
  metrics$cutoff_label <- factor(
    paste0("FDR = ", format(metrics$cutoff, trim = TRUE, nsmall = 2)),
    levels = paste0("FDR = ", format(validation_cutoffs, trim = TRUE, nsmall = 2))
  )
  attr(metrics, "plot_cap") <- cap
  attr(metrics, "use_haldane") <- isTRUE(use_haldane)
  metrics
}

plot_validation_odds_ratio <- function(metrics, dataset_name, benchmark_label,
                                       use_haldane = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot validation results", call. = FALSE)
  }
  metrics <- prepare_validation_plot_data(metrics, use_haldane = use_haldane)
  subtitle <- if (use_haldane) {
    "Haldane--Anscombe 0.5-corrected OR for plotting; triangles mark raw Fisher OR = Inf"
  } else {
    "Odds ratio versus selection FDR cutoff; All = pooled cell-type tables"
  }
  ggplot2::ggplot(metrics, ggplot2::aes(
    x = cutoff, y = odds_ratio_plot, color = method, group = method
  )) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.1, na.rm = TRUE) +
    {if (use_haldane && any(metrics$OR_inf, na.rm = TRUE))
      ggplot2::geom_point(
        data = metrics[metrics$OR_inf, , drop = FALSE],
        ggplot2::aes(x = cutoff, y = odds_ratio_plot, color = method),
        shape = 24, fill = "white", size = 3.1, inherit.aes = FALSE
      )} +
    ggplot2::scale_x_continuous(breaks = validation_cutoffs) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = validation_palette, drop = FALSE) +
    ggplot2::facet_wrap(~ celltype, ncol = 2, drop = FALSE) +
    ggplot2::labs(
      title = paste0(toupper(dataset_name), " ", benchmark_label, " enrichment"),
      subtitle = subtitle,
      x = "Selection FDR cutoff",
      y = if (use_haldane) "Haldane--Anscombe odds ratio (log scale)" else
        "Fisher odds ratio (log scale)", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey60")
    )
}

plot_validation_pvalue <- function(metrics, dataset_name, benchmark_label) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot validation results", call. = FALSE)
  }
  metrics <- prepare_validation_plot_data(metrics)
  # Fisher p-values can be exactly zero after numerical underflow. Cap them at
  # the smallest representable positive value before transforming to -log10.
  pval <- as.numeric(metrics$pval)
  pval <- pmax(pval, .Machine$double.xmin)
  pval[!is.finite(metrics$pval)] <- NA_real_
  metrics$pval_plot <- -log10(pval)
  ggplot2::ggplot(metrics, ggplot2::aes(
    x = cutoff, y = pval_plot, color = method, group = method
  )) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45") +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.1, na.rm = TRUE) +
    ggplot2::scale_x_continuous(breaks = validation_cutoffs) +
    ggplot2::scale_color_manual(values = validation_palette, drop = FALSE) +
    ggplot2::facet_wrap(~ celltype, ncol = 2, drop = FALSE) +
    ggplot2::labs(
      title = paste0(toupper(dataset_name), " ", benchmark_label, " enrichment"),
      subtitle = "Fisher exact-test p-value versus selection FDR cutoff; All = pooled cell-type tables",
      x = "Selection FDR cutoff", y = expression(-log[10](Fisher~p-value)), color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey60")
  )
}

plot_validation_true_positive_count <- function(metrics, dataset_name,
                                                benchmark_label) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot validation results", call. = FALSE)
  }
  if (!"TP" %in% names(metrics)) {
    stop("Validation metrics must contain a TP column", call. = FALSE)
  }
  metrics <- prepare_validation_plot_data(metrics)
  metrics$TP <- as.numeric(metrics$TP)
  ggplot2::ggplot(metrics, ggplot2::aes(
    x = cutoff, y = TP, color = method, group = method
  )) +
    ggplot2::geom_vline(xintercept = 0.10, linetype = 3, color = "grey75") +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.1, na.rm = TRUE) +
    ggplot2::scale_x_continuous(breaks = validation_cutoffs) +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                expand = ggplot2::expansion(mult = c(0.05, 0.10))) +
    ggplot2::scale_color_manual(values = validation_palette, drop = FALSE) +
    ggplot2::facet_wrap(~ celltype, ncol = 2, scales = "free_y", drop = FALSE) +
    ggplot2::labs(
      title = paste0(toupper(dataset_name), " ", benchmark_label,
                     " validated positive links"),
      subtitle = paste0("True-positive count versus selection FDR cutoff; ",
                        "the vertical line marks FDR = 0.10"),
      x = "Selection FDR cutoff", y = "True positives (count)", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey60")
    )
}
