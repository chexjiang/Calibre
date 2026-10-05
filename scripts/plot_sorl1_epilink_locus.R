#!/usr/bin/env Rscript

# Fig. 4d/e-style visualization for the Yang-supported SORL1 locus.
# The top track shows annotated genes, the middle track Yang Table 6d
# variants, and the bottom track strong local EpiLink links. Link colour is
# the pooled beta from exprs ~ atac + log(nCount_RNA) + condition.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Seurat)
  library(Signac)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_dir)

pair_gene <- "SORL1"
pair_peak <- "chr11-121589776-121590720"
pair_rsid <- "rs2276412"
peak_start <- 121589776L
peak_end <- 121590720L
variant_pos <- 121590137L
beta_cutoff <- 0.1

# Use a compromise window: wide enough to show the local neighbouring-gene
# context around SORL1, but without truncating the next gene at the right edge.
# This is chr11:121.00-121.90 Mb (900 kb).
x0 <- 121000000L
x1 <- 121900000L

result_dir <- file.path(project_dir, "results", "real", "Brain2", "Microglia",
                        "ALL", "epilink_poisson_condition_effect")
signac_object_file <- file.path(project_dir, "results", "real", "Brain2",
                                "Microglia", "ALL", "signac", "signac_object.rds")
yang_file <- file.path(project_dir, "results", "validation", "functional",
                       "Brain2_diff", "Microglia", "yang_overlap",
                       "yang_supp_table6d_grch38.tsv.gz")
out_dir <- file.path(project_dir, "results", "validation", "functional", "Brain2",
                     "Microglia", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

epilink <- fread(file.path(result_dir, "all_pairs.tsv.gz"), showProgress = FALSE)
epilink[, `:=`(gene = as.character(gene), peak = as.character(peak),
               beta = as.numeric(beta))]
peak_parts <- tstrsplit(epilink$peak, "-", fixed = TRUE)
epilink[, `:=`(peak_chr = peak_parts[[1L]], peak_start = as.integer(peak_parts[[2L]]),
               peak_end = as.integer(peak_parts[[3L]]))]
target_row <- epilink[gene == pair_gene & peak == pair_peak]
if (nrow(target_row) != 1L) {
  stop("Expected one SORL1 EpiLink row, found ", nrow(target_row), call. = FALSE)
}

# Reuse the exact gene annotation attached to the Brain2 Signac object.
if (!file.exists(signac_object_file)) stop("Missing ", signac_object_file, call. = FALSE)
signac_object <- readRDS(signac_object_file)
annotation <- as.data.table(as.data.frame(Signac::Annotation(signac_object[["ATAC"]])))
annotation[, `:=`(start = as.integer(start), end = as.integer(end),
                  gene_name = as.character(gene_name), strand = as.character(strand),
                  type = as.character(type))]
annotation <- annotation[seqnames == "chr11" & end >= x0 & start <= x1 &
                         !is.na(gene_name) & nzchar(gene_name)]
if (!nrow(annotation)) stop("No gene annotation in plotting interval", call. = FALSE)

gene_summary <- annotation[, .(
  gene_start = min(start), gene_end = max(end), strand = strand[1L]
), by = gene_name][order(gene_start, gene_end)]
gene_summary[, y := rev(seq_len(.N))]
gene_levels <- gene_summary$gene_name
exons <- annotation[type == "exon" & gene_name %in% gene_levels]
exons <- merge(exons, gene_summary[, .(gene_name, y)], by = "gene_name", all = FALSE)
exons[, `:=`(start = pmax(start, x0), end = pmin(end, x1))]
exons <- exons[start <= end]
exon_half_height <- if (nrow(gene_summary) == 1L) 0.045 else 0.13
exons[, `:=`(ymin = y - exon_half_height, ymax = y + exon_half_height)]
gene_y_limits <- if (nrow(gene_summary) == 1L) c(0.5, 1.5) else NULL
gene_y_expand <- if (nrow(gene_summary) == 1L) {
  expansion(mult = c(0, 0))
} else {
  expansion(mult = c(0.08, 0.18))
}

# Yang Table 6d variants in the interval; rs2276412 is highlighted.
yang <- fread(yang_file, showProgress = FALSE)
yang[, `:=`(pos_hg38 = as.integer(pos_hg38), chr_grch38 = as.character(chr_grch38),
            rsid = as.character(rsid))]
variants <- unique(yang[chr_grch38 == "chr11" & pos_hg38 >= x0 & pos_hg38 <= x1,
                        .(pos = pos_hg38, rsid)])
if (!nrow(variants)) stop("No Yang variants in plotting interval", call. = FALSE)
variants[, highlighted := rsid == pair_rsid]

# Match the paper-style sparse link landscape: retain strong local links and
# explicitly retain the Yang-supported target even though its beta is 0.134.
local_links <- epilink[peak_chr == "chr11" & peak_end >= x0 &
                       peak_start <= x1 & is.finite(beta) &
                       (abs(beta) > beta_cutoff |
                          (gene == pair_gene & peak == pair_peak))]
gene_tss <- gene_summary[, .(
  gene_tss = ifelse(strand == "+", gene_start, gene_end)
), by = .(gene = gene_name)]
local_links <- merge(local_links, gene_tss, by = "gene", all = FALSE)
local_links[, peak_mid := (peak_start + peak_end) / 2]
local_links[, highlighted := gene == pair_gene & peak == pair_peak]
if (!nrow(local_links)) stop("No local EpiLink links in plotting interval", call. = FALSE)

make_arc <- function(row, n = 50L) {
  t <- seq(0, 1, length.out = n)
  distance <- abs(row$gene_tss - row$peak_mid)
  height <- max(0.04, min(0.52, distance / (x1 - x0) * 0.9))
  data.table(link_id = row$link_id,
             x = (1 - t) * row$peak_mid + t * row$gene_tss,
             y = -4 * height * t * (1 - t), beta = row$beta,
             highlighted = row$highlighted)
}
local_links[, link_id := .I]
arc_data <- rbindlist(lapply(seq_len(nrow(local_links)), function(i) make_arc(local_links[i])))

gene_plus <- gene_summary[strand == "+"]
gene_minus <- gene_summary[strand == "-"]

p_gene <- ggplot() +
  geom_segment(data = gene_plus,
               aes(x = gene_start, xend = gene_end, y = y, yend = y),
               linewidth = 0.35, colour = "#2C7FB8",
               arrow = arrow(length = unit(0.07, "inches"), type = "closed")) +
  geom_segment(data = gene_minus,
               aes(x = gene_start, xend = gene_end, y = y, yend = y),
               linewidth = 0.35, colour = "#2C7FB8",
               arrow = arrow(length = unit(0.07, "inches"), type = "closed",
                             ends = "first")) +
  geom_rect(data = exons,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
            fill = "#2C7FB8", colour = NA) +
  geom_vline(xintercept = variant_pos, colour = "#74C476", alpha = 0.30,
             linewidth = 5) +
  scale_x_continuous(limits = c(x0, x1), labels = function(x) sprintf("%.2f", x / 1e6)) +
  scale_y_continuous("Genes", breaks = gene_summary$y, labels = gene_summary$gene_name,
                     limits = gene_y_limits, expand = gene_y_expand) +
  labs(title = "A  Gene annotation") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "none")

p_snp <- ggplot() +
  geom_segment(aes(x = x0, xend = x1, y = 0, yend = 0),
               linewidth = 0.3, colour = "grey50") +
  geom_point(data = variants[highlighted == FALSE], aes(x = pos, y = 0),
             shape = 124, size = 3, colour = "grey25") +
  geom_point(data = variants[highlighted == TRUE], aes(x = pos, y = 0),
             shape = 124, size = 5, colour = "#238B45") +
  geom_text(data = variants[highlighted == TRUE], aes(x = pos, y = 0.12, label = rsid),
            colour = "#238B45", hjust = 0, size = 3.4) +
  scale_x_continuous(limits = c(x0, x1), labels = function(x) sprintf("%.2f", x / 1e6)) +
  scale_y_continuous("SNPs", breaks = 0, labels = "Yang GWAS", limits = c(-0.12, 0.34)) +
  labs(title = "B  Yang GWAS variant", x = "chr11 position (Mb)") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"), legend.position = "none")

p_links <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_path(data = arc_data[highlighted == FALSE],
            aes(x = x, y = y, group = link_id, colour = beta),
            linewidth = 0.35, alpha = 0.24) +
  # A neutral halo makes the target arc visible without changing beta's
  # colour encoding.
  geom_path(data = arc_data[highlighted == TRUE],
            aes(x = x, y = y, group = link_id),
            linewidth = 2.8, colour = "black", alpha = 0.9) +
  geom_path(data = arc_data[highlighted == TRUE],
            aes(x = x, y = y, group = link_id, colour = beta),
            linewidth = 1.7, alpha = 1) +
  geom_point(data = local_links, aes(x = peak_mid, y = 0, colour = beta),
             size = 0.9, alpha = 0.55) +
  geom_point(data = local_links[highlighted == TRUE], aes(x = peak_mid, y = 0),
             size = 3.4, colour = "black") +
  geom_point(data = local_links[highlighted == TRUE], aes(x = peak_mid, y = 0, colour = beta),
             size = 2.1) +
  geom_text(data = local_links[highlighted == TRUE],
            aes(x = (peak_mid + gene_tss) / 2, y = -0.12,
                label = paste0(pair_gene, "–", pair_rsid)),
            fontface = "bold", size = 3.2, colour = "black") +
  scale_x_continuous(limits = c(x0, x1), labels = function(x) sprintf("%.2f", x / 1e6)) +
  scale_y_continuous("EpiLink links", breaks = NULL, limits = c(-0.58, 0.08)) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         midpoint = 0, name = "EpiLink beta") +
  labs(title = "C  Local EpiLink link landscape",
       subtitle = sprintf("%d links with |beta| > %.1f plus highlighted target: %s–%s (beta = %+.3f)",
                          nrow(local_links), beta_cutoff, pair_gene, pair_rsid,
                          target_row$beta),
       x = "chr11 position (Mb)") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "right")

fig <- p_gene / p_snp / p_links +
  plot_layout(heights = c(1.55, 0.75, 2.6)) +
  plot_annotation(
    title = "SORL1–rs2276412 peak–gene association in Brain2 Microglia",
    subtitle = "Pooled EpiLink model: exprs ~ atac + log(nCount_RNA) + condition; GRCh38",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10))
  )

png_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink_ALL.png")
pdf_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink_ALL.pdf")
links_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink_ALL_links.tsv.gz")
fwrite(local_links, links_file, sep = "\t")
ggsave(png_file, fig, width = 10.5, height = 9.5, dpi = 320)
ggsave(pdf_file, fig, width = 10.5, height = 9.5, device = cairo_pdf)
message("Wrote ", normalizePath(png_file))
message("Wrote ", normalizePath(pdf_file))
message("Wrote ", normalizePath(links_file))
