#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
})

project_dir <- if (length(commandArgs(trailingOnly = TRUE))) {
  normalizePath(commandArgs(trailingOnly = TRUE)[[1L]])
} else {
  normalizePath(".")
}

out_dir <- file.path(
  project_dir, "results", "validation", "functional", "Brain2_diff", "figures"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pair_gene <- "NEAT1"
pair_peak <- "chr11-65498624-65499539"

diff_file <- file.path(
  project_dir, "results", "real", "Brain2_diff", "Astrocytes",
  "epilink_diff_poisson_no_donor_effect_no_donor_robust", "all_pairs.tsv.gz"
)
ordinary_file <- file.path(
  project_dir, "results", "real", "Brain2", "Astrocytes",
  "epilink_poisson", "all_pairs.tsv.gz"
)
scmultimap_file <- file.path(
  project_dir, "results", "real", "Brain2_diff", "Astrocytes",
  "scmultimap", "all_pairs.tsv.gz"
)

diff <- fread(diff_file)[gene == pair_gene & peak == pair_peak]
ordinary <- fread(ordinary_file)[gene == pair_gene & peak == pair_peak]
scmultimap <- fread(scmultimap_file)[gene == pair_gene & peak == pair_peak]
if (nrow(diff) != 1L || nrow(ordinary) != 1L || nrow(scmultimap) != 1L) {
  stop("The requested NEAT1 pair was not found exactly once in all three result tables.")
}

gene_gr <- genes(
  EnsDb.Hsapiens.v86,
  filter = GeneNameFilter(pair_gene)
)
exon_gr <- exons(
  EnsDb.Hsapiens.v86,
  filter = GeneNameFilter(pair_gene)
)
gene_start <- start(gene_gr)[[1L]]
gene_end <- end(gene_gr)[[1L]]
peak_start <- as.integer(sub(".*-(\\d+)-\\d+$", "\\1", pair_peak))
peak_end <- as.integer(sub(".*-(\\d+)$", "\\1", pair_peak))
peak_mid <- (peak_start + peak_end) / 2
plot_start <- min(gene_start, peak_start) - 7000L
plot_end <- max(gene_end, peak_end) + 7000L

exon_df <- data.frame(
  start = start(exon_gr), end = end(exon_gr), y = 0.52
)
gene_df <- data.frame(start = gene_start, end = gene_end, y = 0.52)
peak_df <- data.frame(start = peak_start, end = peak_end, y = 1.04)
link_df <- data.frame(
  x = gene_start, xend = peak_mid, y = 0.64, yend = 0.96
)

panel_a <- ggplot() +
  geom_segment(
    data = gene_df, aes(x = start, xend = end, y = y, yend = y),
    linewidth = 0.55, colour = "#2B5C8A"
  ) +
  geom_rect(
    data = exon_df,
    aes(xmin = start, xmax = end, ymin = y - 0.09, ymax = y + 0.09),
    fill = "#2B5C8A", colour = NA
  ) +
  geom_curve(
    data = link_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    curvature = -0.35, linewidth = 0.55, colour = "#6B7280",
    arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")
  ) +
  geom_rect(
    data = peak_df,
    aes(xmin = start, xmax = end, ymin = y - 0.10, ymax = y + 0.10),
    fill = "#D97706", colour = NA
  ) +
  annotate("text", x = (gene_start + gene_end) / 2, y = 0.27,
           label = "NEAT1", colour = "#2B5C8A", fontface = "bold", size = 3.7) +
  annotate("text", x = peak_mid, y = 1.30,
           label = "candidate ATAC peak", colour = "#9A5B00", size = 3.2) +
  annotate("text", x = peak_mid, y = 1.17,
           label = paste0("chr11:", format(peak_start, big.mark = ","), "–",
                          format(peak_end, big.mark = ",")),
           colour = "#9A5B00", size = 2.9) +
  scale_x_continuous(
    limits = c(plot_start, plot_end),
    labels = function(x) format(round(x / 1e6, 3), trim = TRUE)
  ) +
  scale_y_continuous(NULL, breaks = NULL, limits = c(0.05, 1.45)) +
  labs(title = "A  NEAT1 locus and candidate peak", x = "chr11 position (Mb)") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

beta_df <- data.frame(
  condition = factor(c("Ctrl", "AD"), levels = c("Ctrl", "AD")),
  beta = c(diff$beta_control, diff$beta_disease),
  colour = c("Ctrl", "AD")
)
panel_b <- ggplot(beta_df, aes(y = condition, x = beta, colour = colour)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#6B7280") +
  geom_segment(aes(x = 0, xend = beta, y = condition, yend = condition), linewidth = 1.3) +
  geom_point(size = 3.4) +
  geom_text(
    aes(label = sprintf("%+.3f", beta)), hjust = ifelse(beta_df$beta >= 0, -0.22, 1.22),
    colour = "#374151", size = 3.4
  ) +
  scale_colour_manual(values = c(Ctrl = "#2563EB", AD = "#DC2626"), guide = "none") +
  scale_x_continuous(limits = c(-0.16, 0.16), breaks = c(-0.15, -0.1, 0, 0.1, 0.15)) +
  labs(
    title = "B  EpiLink condition-specific association",
    subtitle = sprintf("βdiff (AD − Ctrl) = %+.3f", diff$beta_diff),
    x = "EpiLink association coefficient (β)", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(colour = "#4B5563"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

method_df <- data.frame(
  method = factor(
    c("Ordinary EpiLink", "scMultiMap-Diff", "EpiLink-Diff"),
    levels = c("EpiLink-Diff", "scMultiMap-Diff", "Ordinary EpiLink")
  ),
  p = c(ordinary$pval, scmultimap$pval, diff$wald_p),
  status = c("not selected", "not selected", "selected")
)
method_df$minus_log10_p <- -log10(pmax(method_df$p, .Machine$double.xmin))
panel_c <- ggplot(method_df, aes(x = minus_log10_p, y = method, colour = status)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = "#6B7280") +
  geom_segment(aes(x = 0, xend = minus_log10_p, y = method, yend = method), linewidth = 1.1) +
  geom_point(size = 3.2) +
  geom_text(
    aes(label = ifelse(p < 1e-10, "p < 1e−10", sprintf("p = %.3g", p))),
    hjust = ifelse(method_df$method == "EpiLink-Diff", 1.05, -0.08),
    colour = "#374151", size = 3.2
  ) +
  scale_colour_manual(values = c(selected = "#B91C1C", `not selected` = "#6B7280"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "C  Differential signal is specific to EpiLink-Diff",
    subtitle = "Dashed line: nominal p = 0.05",
    x = expression(-log[10](p)), y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(colour = "#4B5563"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

figure <- panel_a / (panel_b | panel_c) +
  plot_annotation(
    title = "Astrocyte NEAT1 peak–gene association in Brain2",
    subtitle = "A condition-dependent association detected by EpiLink-Diff",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(colour = "#4B5563", size = 10)
    )
  )

ggsave(file.path(out_dir, "NEAT1_epilink_diff_association.pdf"), figure,
       width = 10, height = 8.2, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "NEAT1_epilink_diff_association.png"), figure,
       width = 10, height = 8.2, units = "in", dpi = 320)

writeLines(
  c(
    "NEAT1 Astrocyte EpiLink-Diff visualization",
    paste0("Pair: ", pair_gene, " / ", pair_peak),
    paste0("EpiLink-Diff beta_diff: ", signif(diff$beta_diff, 5)),
    paste0("EpiLink-Diff wald p: ", format(diff$wald_p, scientific = TRUE)),
    paste0("Ordinary EpiLink p: ", format(ordinary$pval, scientific = TRUE)),
    paste0("scMultiMap-Diff p: ", format(scmultimap$pval, scientific = TRUE)),
    "Interpretation: positive beta_diff means stronger peak-gene association in AD.",
    "This is a condition-dependent association, not a gene-expression fold change."
  ),
  file.path(out_dir, "NEAT1_epilink_diff_association.txt")
)

message("Wrote NEAT1 figure to ", normalizePath(out_dir))
