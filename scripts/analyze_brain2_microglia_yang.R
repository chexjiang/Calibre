#!/usr/bin/env Rscript

# Integrate Yang et al. (Nature Genetics 2023) Supplementary Table 6 with
# Brain2 Microglia candidate pairs. Table 6d already contains GRCh38 variant
# positions (pos hg38), so no liftover is applied. A valid match requires both
# (i) the GRCh38 variant to fall inside the Brain2 peak and (ii) the Brain2
# gene to equal a Yang pcHi-C interacting target gene. EpiLink-Diff hits are
# then selected from the exact matches.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_dir)

out_dir <- file.path(project_dir, "results", "validation", "functional", "Brain2_diff", "Microglia")
overlap_dir <- file.path(out_dir, "yang_overlap")
figure_dir <- file.path(out_dir, "figures")
dir.create(overlap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

yang_xlsx <- file.path(project_dir, "reference", "GWAS", "yang_2023_MOESM3.xlsx")
epilink_dir <- file.path(project_dir, "results", "real", "Brain2_diff", "Microglia",
                         "epilink_diff_poisson_no_donor_effect_no_donor_robust")
sc_dir <- file.path(project_dir, "results", "real", "Brain2_diff", "Microglia", "scmultimap")
if (!file.exists(yang_xlsx)) stop("Missing ", yang_xlsx, call. = FALSE)

write_tsv_gz <- function(x, path) {
  if (nrow(x)) {
    fwrite(x, path, sep = "\t")
  } else {
    con <- gzfile(path, open = "wt")
    if (ncol(x)) writeLines(paste(names(x), collapse = "\t"), con)
    close(con)
  }
}

read_yang <- function(sheet, names_out) {
  x <- as.data.table(read_excel(yang_xlsx, sheet = sheet, skip = 1))
  if (ncol(x) != length(names_out)) stop("Unexpected columns in ", sheet, call. = FALSE)
  setnames(x, names_out)
  x
}

yang6c <- read_yang(
  "Supplementary Table 6c",
  c("gwas_locus", "rsid", "chr", "pos_hg19", "pos_hg38", "gwas_log10p",
    "paintor_pip", "paintor_rank", "fgwas_ppa", "fgwas_rank", "variant_type",
    "prioritized_schwartzentruber", "new_ccre", "eqtl_rasqual", "eqtl_lmm",
    "eqtl_kosoy", "atac_start_hg19", "atac_end_hg19", "atac_score",
    "max_pchic_score", "aggregate_pchic_expressed", "aggregate_pchic_eqtl")
)
yang6d <- read_yang(
  "Supplementary Table 6d",
  c("gwas_locus", "rsid", "chr", "pos_hg19", "pos_hg38", "pchic_residing_start_hg19",
    "pchic_residing_end_hg19", "residing_genes", "pchic_interacting_start_hg19",
    "pchic_interacting_end_hg19", "target_genes", "target_rpkm")
)
yang6c[, `:=`(chr = as.integer(chr), pos_hg38 = as.integer(pos_hg38),
              chr_grch38 = paste0("chr", as.integer(chr)), genome_build = "GRCh38")]
yang6d[, `:=`(chr = as.integer(chr), pos_hg38 = as.integer(pos_hg38),
              chr_grch38 = paste0("chr", as.integer(chr)), genome_build = "GRCh38")]
fwrite(yang6c, file.path(overlap_dir, "yang_supp_table6c_grch38.tsv.gz"), sep = "\t")
fwrite(yang6d, file.path(overlap_dir, "yang_supp_table6d_grch38.tsv.gz"), sep = "\t")

target_long <- yang6d[!is.na(pos_hg38) & !is.na(target_genes),
  .(target_gene = unlist(strsplit(gsub("\\s+", "", target_genes), ";", fixed = TRUE))),
  by = .(gwas_locus, rsid, chr_grch38, pos_hg38, target_rpkm)]
target_long <- target_long[target_gene != "." & nzchar(target_gene)]

candidate <- as.data.table(readRDS(file.path(epilink_dir, "candidate_pairs.rds")))
candidate <- candidate[, .(gene = as.character(gene), peak = as.character(peak))]
peak_parts <- tstrsplit(candidate$peak, "-", fixed = TRUE)
candidate[, `:=`(chr_grch38 = peak_parts[[1L]], peak_start = as.integer(peak_parts[[2L]]),
                 peak_end = as.integer(peak_parts[[3L]]))]
setkey(candidate, chr_grch38, peak_start, peak_end)

variant_intervals <- target_long[, .(
  chr_grch38, variant_start = pos_hg38, variant_end = pos_hg38,
  gwas_locus, rsid, pos_hg38, target_gene, target_rpkm
)]
overlap_all <- foverlaps(
  variant_intervals, candidate,
  by.x = c("chr_grch38", "variant_start", "variant_end"),
  by.y = c("chr_grch38", "peak_start", "peak_end"),
  type = "within", nomatch = 0L
)
overlap_all[, gene_match := target_gene == gene]
overlap_gene <- unique(overlap_all[gene_match == TRUE])
fwrite(overlap_all, file.path(overlap_dir, "yang_brain2_peak_overlaps.tsv.gz"), sep = "\t")
fwrite(overlap_gene, file.path(overlap_dir, "yang_brain2_peak_gene_matches.tsv.gz"), sep = "\t")

ep_all <- fread(file.path(epilink_dir, "all_pairs.tsv.gz"), showProgress = FALSE)
ep_sig <- fread(file.path(epilink_dir, "significant_pairs.tsv.gz"), showProgress = FALSE)
sc_all <- fread(file.path(sc_dir, "all_pairs.tsv.gz"), showProgress = FALSE)
sc_sig <- fread(file.path(sc_dir, "significant_pairs.tsv.gz"), showProgress = FALSE)

ep_cols <- intersect(c("gene", "peak", "beta_control", "beta_disease", "beta_diff", "wald_p", "padj", "calibrated_q", "significant"), names(ep_all))
sc_cols <- intersect(c("gene", "peak", "covar_control", "covar_disease", "cor_control", "cor_disease", "pval", "padj"), names(sc_all))
ep_hit <- merge(overlap_gene, ep_all[, ..ep_cols], by = c("gene", "peak"), all.x = TRUE)
sc_hit <- merge(overlap_gene, sc_all[, ..sc_cols], by = c("gene", "peak"), all.x = TRUE)
ep_sig_hit <- merge(overlap_gene, ep_sig, by = c("gene", "peak"), all = FALSE)
sc_sig_hit <- merge(overlap_gene, sc_sig, by = c("gene", "peak"), all = FALSE)
write_tsv_gz(ep_hit, file.path(overlap_dir, "yang_epilink_all_peak_gene_matches.tsv.gz"))
write_tsv_gz(ep_sig_hit, file.path(overlap_dir, "yang_epilinkdiff_significant_matches.tsv.gz"))
write_tsv_gz(sc_hit, file.path(overlap_dir, "yang_scmultimap_all_peak_gene_matches.tsv.gz"))
write_tsv_gz(sc_sig_hit, file.path(overlap_dir, "yang_scmultimap_significant_matches.tsv.gz"))

plot_match <- function(hit, index) {
  pair <- unique(hit[, .(gene, peak, chr_grch38, peak_start, peak_end, rsid, pos_hg38, gwas_locus)])[1L]
  ep <- ep_all[gene == pair$gene & peak == pair$peak][1L]
  sc <- sc_all[gene == pair$gene & peak == pair$peak][1L]
  x0 <- min(pair$peak_start, pair$pos_hg38) - 5000
  x1 <- max(pair$peak_end, pair$pos_hg38) + 5000
  p1 <- ggplot() +
    geom_segment(aes(x = x0, xend = x1, y = 0.5, yend = 0.5), linewidth = 0.5) +
    geom_rect(aes(xmin = pair$peak_start, xmax = pair$peak_end, ymin = 0.4, ymax = 0.6), fill = "#D97706") +
    geom_vline(aes(xintercept = pair$pos_hg38), colour = "#2B8CBE", linewidth = 0.8) +
    annotate("text", x = pair$pos_hg38, y = 0.85, label = pair$rsid, colour = "#2B8CBE", hjust = 0) +
    annotate("text", x = mean(c(pair$peak_start, pair$peak_end)), y = 0.25, label = pair$gene, fontface = "bold") +
    scale_x_continuous(limits = c(x0, x1), labels = function(z) round(z / 1e6, 3)) +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0, 1)) +
    labs(title = paste0("A  ", pair$gene, " locus and GWAS variant"), x = paste0(pair$chr_grch38, " position (Mb)")) +
    theme_bw(base_size = 10) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  beta <- data.table(condition = factor(c("Ctrl", "AD"), levels = c("Ctrl", "AD")), beta = c(ep$beta_control, ep$beta_disease))
  p2 <- ggplot(beta, aes(condition, beta, colour = condition)) +
    geom_hline(yintercept = 0, colour = "grey60") + geom_point(size = 3) + geom_line(aes(group = 1), linewidth = 0.7) +
    scale_colour_manual(values = c(Ctrl = "#2B8CBE", AD = "#D7301F"), guide = "none") +
    labs(title = "B  EpiLink condition-specific association", y = "ATAC slope (GLM linear-predictor scale)", x = NULL,
         subtitle = sprintf("beta_diff = %+.3f; Wald p = %.3g", ep$beta_diff, ep$wald_p)) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  p3 <- ggplot(data.table(condition = factor(c("Ctrl", "AD"), levels = c("Ctrl", "AD")), correlation = c(sc$cor_control, sc$cor_disease)), aes(condition, correlation, colour = condition)) +
    geom_hline(yintercept = 0, colour = "grey60") + geom_point(size = 3) + geom_line(aes(group = 1), linewidth = 0.7) +
    scale_colour_manual(values = c(Ctrl = "#2B8CBE", AD = "#D7301F"), guide = "none") +
    labs(title = "C  scMultiMap association", y = "Latent correlation", x = NULL, subtitle = sprintf("Differential p = %.3g", sc$pval)) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  fig <- p1 / (p2 | p3) + plot_annotation(title = paste0(pair$gene, " GWAS-supported peak–gene association"))
  stem <- sprintf("GWAS_match_%02d_%s_%s", index, pair$gene, gsub("[^A-Za-z0-9]+", "_", pair$rsid))
  ggsave(file.path(figure_dir, paste0(stem, ".png")), fig, width = 9, height = 7, dpi = 320)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), fig, width = 9, height = 7, device = cairo_pdf)
}

if (nrow(ep_sig_hit)) {
  hits <- unique(ep_sig_hit[, .(gene, peak, chr_grch38, peak_start, peak_end, rsid, pos_hg38, gwas_locus)])
  for (i in seq_len(min(nrow(hits), 10L))) plot_match(hits[i], i)
  figure_status <- paste0("Generated ", min(nrow(hits), 10L), " Fig. 4d/e-style figure(s).")
} else {
  figure_status <- "No exact Yang-variant/Brain2-peak/Yang-target/EpiLinkDiff-significant match; no positive Fig. 4d/e-style figure was generated."
  writeLines(figure_status, file.path(figure_dir, "NO_FIGURE_EXACT_EpiLinkDiff_MATCH.txt"))
}

summary <- data.table(
  quantity = c("Yang Table 6d rows", "Yang Table 6d unique variants", "Yang Table 6d target genes",
               "Brain2 candidate pairs", "variant-containing Brain2 pairs",
               "variant-containing pairs with Yang target-gene match",
               "EpiLinkDiff significant exact matches", "scMultiMap significant exact matches"),
  value = c(nrow(yang6d), uniqueN(yang6d$rsid), uniqueN(target_long$target_gene), nrow(candidate),
            uniqueN(overlap_all[, .(gene, peak)]), uniqueN(overlap_gene[, .(gene, peak)]),
            uniqueN(ep_sig_hit[, .(gene, peak)]), uniqueN(sc_sig_hit[, .(gene, peak)]))
)
fwrite(summary, file.path(out_dir, "yang_overlap_summary.tsv"), sep = "\t")
writeLines(c(
  "Brain2 Microglia Yang Supplementary Table 6 integration",
  "Yang Table 6c/6d pos hg38 was used directly; genome_build=GRCh38.",
  paste0("Figure status: ", figure_status)
), file.path(out_dir, "YANG_overlap_README.txt"))
message(figure_status)
message("Wrote overlap results to ", normalizePath(overlap_dir))
