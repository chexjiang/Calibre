#!/usr/bin/env Rscript

# Compact single-track Fig. 4d/e-style plot for the Yang-supported SORL1 pair.
# The coordinate system contains gene bodies, Yang variants, Brain2 peaks, and
# EpiLink arcs without separate boxed panels.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Seurat)
  library(Signac)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_dir)

pair_gene <- "SORL1"
pair_peak <- "chr11-121589776-121590720"
pair_rsid <- "rs2276412"
variant_pos <- 121590137L
beta_cutoff <- 0.1

# A compact interval retaining the relevant neighbouring-gene context.
x0 <- 121270000L
x1 <- 121750000L

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

signac_object <- readRDS(signac_object_file)
annotation <- as.data.table(as.data.frame(Signac::Annotation(signac_object[["ATAC"]])))
annotation[, `:=`(start = as.integer(start), end = as.integer(end),
                  gene_name = as.character(gene_name), strand = as.character(strand),
                  type = as.character(type))]
annotation <- annotation[seqnames == "chr11" & end >= x0 & start <= x1 &
                         !is.na(gene_name) & nzchar(gene_name)]
gene_summary <- annotation[, .(
  gene_start = min(start), gene_end = max(end), strand = strand[1L]
), by = gene_name][order(gene_start, gene_end)]
gene_summary[, y := 2.05]
exons <- annotation[type == "exon" & gene_name %in% gene_summary$gene_name]
exons <- merge(exons, gene_summary[, .(gene_name, y)], by = "gene_name", all = FALSE)
exons[, `:=`(start = pmax(start, x0), end = pmin(end, x1),
             ymin = y - 0.075, ymax = y + 0.075)]
exons <- exons[start <= end]

yang <- fread(yang_file, showProgress = FALSE)
yang[, `:=`(pos_hg38 = as.integer(pos_hg38), chr_grch38 = as.character(chr_grch38),
            rsid = as.character(rsid))]
variants <- unique(yang[chr_grch38 == "chr11" & pos_hg38 >= x0 & pos_hg38 <= x1,
                        .(pos = pos_hg38, rsid)])
variants[, highlighted := rsid == pair_rsid]

local_links <- epilink[peak_chr == "chr11" & peak_end >= x0 & peak_start <= x1 &
                       is.finite(beta) &
                       (abs(beta) > beta_cutoff |
                          (gene == pair_gene & peak == pair_peak))]
gene_tss <- gene_summary[, .(
  gene_tss = ifelse(strand == "+", gene_start, gene_end)
), by = .(gene = gene_name)]
local_links <- merge(local_links, gene_tss, by = "gene", all = FALSE)
local_links[, `:=`(peak_mid = (peak_start + peak_end) / 2,
                   highlighted = gene == pair_gene & peak == pair_peak)]
if (!nrow(local_links)) stop("No local links in plotting interval", call. = FALSE)

peaks <- unique(local_links[, .(
  peak, peak_start, peak_end,
  highlighted = any(highlighted)
)])

make_arc <- function(row, n = 60L) {
  t <- seq(0, 1, length.out = n)
  distance <- abs(row$gene_tss - row$peak_mid)
  height <- max(0.12, min(0.85, distance / (x1 - x0) * 1.25))
  data.table(
    link_id = row$link_id,
    x = (1 - t) * row$peak_mid + t * row$gene_tss,
    y = 1.00 - 4 * height * t * (1 - t),
    beta = row$beta,
    highlighted = row$highlighted
  )
}
local_links[, link_id := .I]
arc_data <- rbindlist(lapply(seq_len(nrow(local_links)), function(i) make_arc(local_links[i])))

gene_plus <- gene_summary[strand == "+"]
gene_minus <- gene_summary[strand == "-"]
target_arc <- arc_data[highlighted == TRUE]
target_label <- local_links[highlighted == TRUE]

p <- ggplot() +
  # Yang variant highlight behind the gene/peak tracks.
  geom_rect(aes(xmin = variant_pos - 900, xmax = variant_pos + 900,
                ymin = 1.15, ymax = 2.62),
            fill = "#A1D99B", alpha = 0.48, colour = NA) +
  geom_segment(data = gene_plus,
               aes(x = gene_start, xend = gene_end, y = y, yend = y),
               linewidth = 0.45, colour = "#253494",
               arrow = arrow(length = unit(0.075, "inches"), type = "closed")) +
  geom_segment(data = gene_minus,
               aes(x = gene_start, xend = gene_end, y = y, yend = y),
               linewidth = 0.45, colour = "#253494",
               arrow = arrow(length = unit(0.075, "inches"), type = "closed",
                             ends = "first")) +
  geom_rect(data = exons,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
            fill = "#253494", colour = NA) +
  geom_text(data = gene_summary,
            aes(x = (gene_start + gene_end) / 2, y = 2.25, label = gene_name),
            size = 3.3, colour = "#1F1F1F", fontface = "bold") +
  geom_text(data = variants[highlighted == TRUE],
            aes(x = pos, y = 2.64, label = rsid),
            hjust = 0, size = 3.4, colour = "#238B45", fontface = "bold") +
  geom_rect(data = peaks[highlighted == FALSE],
            aes(xmin = peak_start, xmax = peak_end, ymin = 1.22, ymax = 1.38),
            fill = "#737373", colour = NA) +
  geom_rect(data = peaks[highlighted == TRUE],
            aes(xmin = peak_start, xmax = peak_end, ymin = 1.20, ymax = 1.40),
            fill = "#74C476", colour = "#238B45", linewidth = 0.35) +
  geom_segment(aes(x = x0, xend = x1, y = 1.00, yend = 1.00),
               colour = "#BDBDBD", linewidth = 0.45) +
  geom_path(data = arc_data[highlighted == FALSE],
            aes(x = x, y = y, group = link_id, colour = beta),
            linewidth = 0.4, alpha = 0.22) +
  geom_path(data = target_arc,
            aes(x = x, y = y, group = link_id),
            linewidth = 2.8, colour = "black", alpha = 0.9) +
  geom_path(data = target_arc,
            aes(x = x, y = y, group = link_id, colour = beta),
            linewidth = 1.75, alpha = 1) +
  geom_text(data = target_label,
            aes(x = (peak_mid + gene_tss) / 2, y = 0.78,
                label = paste0(pair_gene, "–", pair_rsid)),
            size = 3.2, colour = "#1F1F1F", fontface = "bold") +
  annotate("text", x = x0 - 18000, y = 2.05, label = "Genes", hjust = 1,
           size = 3.4, fontface = "bold") +
  annotate("text", x = x0 - 18000, y = 1.30, label = "Peaks", hjust = 1,
           size = 3.4, fontface = "bold") +
  annotate("text", x = x0 - 18000, y = 0.98, label = "Links", hjust = 1,
           size = 3.4, fontface = "bold") +
  scale_x_continuous(limits = c(x0 - 90000, x1),
                     labels = function(x) format(x, big.mark = ",",
                                                 scientific = FALSE, trim = TRUE)) +
  scale_y_continuous(limits = c(0.05, 2.75), breaks = NULL, expand = c(0, 0)) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         midpoint = 0, name = "EpiLink beta") +
  labs(x = "chr11 position (bp)") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), legend.position = "right",
        plot.margin = margin(5, 18, 5, 12))

png_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink_ALL_compact.png")
pdf_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink_ALL_compact.pdf")
links_file <- file.path(out_dir, "Fig4_style_SORL1_rs2276412_EpiLink_ALL_compact_links.tsv.gz")
fwrite(local_links, links_file, sep = "\t")
ggsave(png_file, p, width = 11, height = 3.6, dpi = 320)
ggsave(pdf_file, p, width = 11, height = 3.6, device = cairo_pdf)
message("Wrote ", normalizePath(png_file))
message("Wrote ", normalizePath(pdf_file))
message("Wrote ", normalizePath(links_file))
