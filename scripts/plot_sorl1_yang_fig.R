#!/usr/bin/env Rscript

# Plot a Yang-supported SORL1 peak-gene pair in a Fig. 4d/e-inspired layout.
# The locus/peak/variant match comes from the non-differential EpiLink
# condition-effect run. Condition-specific slopes are read from the saved
# EpiLinkDiff fit because ordinary EpiLink stores only the pooled beta.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_dir)

pair_gene <- "SORL1"
pair_peak <- "chr11-121589776-121590720"
yang_match_file <- file.path(
  project_dir, "results", "validation", "functional", "Brain2_diff",
  "Microglia", "yang_overlap", "all_current_methods",
  "epilink_poisson_condition_effect_yang_matches.tsv.gz"
)
epilink_file <- file.path(
  project_dir, "results", "real", "Brain2", "Microglia",
  "epilink_poisson_condition_effect", "all_pairs.tsv.gz"
)
diff_file <- file.path(
  project_dir, "results", "real", "Brain2_diff", "Microglia",
  "epilink_diff_poisson_no_donor_effect_donor_robust", "all_pairs.tsv.gz"
)
out_dir <- file.path(
  project_dir, "results", "validation", "functional", "Brain2_diff",
  "Microglia", "figures"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

yang <- fread(yang_match_file, showProgress = FALSE)
yang <- yang[gene == pair_gene & peak == pair_peak]
if (nrow(yang) != 1L) stop("Expected one SORL1 Yang match, found ", nrow(yang), call. = FALSE)

pooled <- fread(epilink_file, showProgress = FALSE)[gene == pair_gene & peak == pair_peak]
diff <- fread(diff_file, showProgress = FALSE)[gene == pair_gene & peak == pair_peak]
if (nrow(pooled) != 1L || nrow(diff) != 1L) stop("SORL1 pair missing from EpiLink output", call. = FALSE)

peak_start <- 121589776
peak_end <- 121590720
variant_pos <- 121590137
# GRCh38 SORL1 gene span from ClinGen/NCBI-Ensembl-UCSC.
gene_start <- 121452314
gene_end <- 121633763
x0 <- 121575000
x1 <- 121605000
gene_start_plot <- max(gene_start, x0)
gene_end_plot <- min(gene_end, x1)

locus <- data.table(
  feature = c("SORL1", "Brain2 peak"),
  start = c(gene_start_plot, peak_start),
  end = c(gene_end_plot, peak_end),
  y = c(2, 1)
)

p_locus <- ggplot() +
  geom_segment(data = locus[locus$feature == "SORL1"],
               aes(x = start, xend = end, y = y, yend = y),
               linewidth = 1.2, colour = "#2C7FB8") +
  geom_rect(data = locus[locus$feature == "Brain2 peak"],
            aes(xmin = start, xmax = end, ymin = y - 0.18, ymax = y + 0.18),
            fill = "#F28E2B", colour = NA) +
  geom_vline(xintercept = variant_pos, linewidth = 0.8, colour = "#2CA25F") +
  annotate("text", x = gene_start_plot, y = 2.25, label = "SORL1", hjust = 0,
           fontface = "bold", size = 4) +
  annotate("text", x = mean(c(peak_start, peak_end)), y = 0.68,
           label = "Brain2 peak", hjust = 0.5, size = 3.3) +
  annotate("text", x = variant_pos, y = 1.43, label = "rs2276412",
           angle = 90, vjust = -0.35, colour = "#2CA25F", size = 3.3) +
  scale_x_continuous(limits = c(x0, x1), labels = function(x) sprintf("%.3f", x / 1e6)) +
  scale_y_continuous(NULL, breaks = c(1, 2), labels = c("peak", "gene"), limits = c(0.3, 2.55)) +
  labs(title = "A  Yang-supported SORL1 locus", x = "chr11 position (Mb)",
       subtitle = "GRCh38; rs2276412 overlaps the Brain2 peak and Yang target gene SORL1") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

strength <- data.table(
  condition = factor(c("CTRL", "AD"), levels = c("CTRL", "AD")),
  beta = c(diff$beta_control, diff$beta_disease)
)

p_strength <- ggplot(strength, aes(condition, beta, group = 1)) +
  geom_hline(yintercept = 0, colour = "grey65") +
  geom_line(linewidth = 0.8, colour = "#7A0177") +
  geom_point(aes(colour = condition), size = 4) +
  geom_text(aes(label = sprintf("%+.3f", beta)), vjust = -1.1, size = 3.8) +
  scale_colour_manual(values = c(CTRL = "#2B8CBE", AD = "#D7301F"), guide = "none") +
  labs(title = "B  EpiLink condition-specific link strength",
       x = NULL, y = "ATAC slope (beta)",
       subtitle = sprintf("beta_diff = %+.3f; donor-robust Wald p = %.3g",
                          diff$beta_diff, diff$wald_p)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

fig <- p_locus / p_strength + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(title = "SORL1–rs2276412 peak-gene association in Brain2 Microglia",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

png_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink.png")
pdf_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink.pdf")
ggsave(png_file, fig, width = 9, height = 7, dpi = 320)
ggsave(pdf_file, fig, width = 9, height = 7, device = cairo_pdf)

message("Wrote ", normalizePath(png_file))
message("Wrote ", normalizePath(pdf_file))
