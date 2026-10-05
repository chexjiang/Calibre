#!/usr/bin/env Rscript

# Direct GO enrichment and Yang et al. Supplementary Table 6 integration for
# Brain2 Microglia.  The Yang table is already reported in GRCh38; this script
# makes the build explicit, intersects variant positions with Brain2 peaks,
# requires a matching Yang target gene, and then filters to EpiLink-Diff hits.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_dir)

out_dir <- file.path(
  project_dir, "results", "validation", "functional", "Brain2_diff", "Microglia"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
go_dir <- file.path(out_dir, "go")
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
overlap_dir <- file.path(out_dir, "yang_overlap")
dir.create(overlap_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(out_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

epilink_dir <- file.path(
  project_dir, "results", "real", "Brain2_diff", "Microglia",
  "epilink_diff_poisson_no_donor_effect_no_donor_robust"
)
scmultimap_dir <- file.path(
  project_dir, "results", "real", "Brain2_diff", "Microglia", "scmultimap"
)
yang_xlsx <- file.path(project_dir, "reference", "GWAS", "yang_2023_MOESM3.xlsx")
if (!file.exists(yang_xlsx)) {
  stop("Missing Yang supplementary workbook: ", yang_xlsx, call. = FALSE)
}

read_pairs <- function(path) {
  x <- fread(path, showProgress = FALSE)
  x[, gene := as.character(gene)]
  x[, peak := as.character(peak)]
  x
}

epilink_all <- read_pairs(file.path(epilink_dir, "all_pairs.tsv.gz"))
epilink_sig <- read_pairs(file.path(epilink_dir, "significant_pairs.tsv.gz"))
sc_all <- read_pairs(file.path(scmultimap_dir, "all_pairs.tsv.gz"))
sc_sig <- read_pairs(file.path(scmultimap_dir, "significant_pairs.tsv.gz"))

# -------------------------------------------------------------------------
# 1. Direct GO enrichment
# -------------------------------------------------------------------------

run_go <- function(label, significant, tested, prefix) {
  sig_genes <- sort(unique(significant$gene[!is.na(significant$gene)]))
  universe <- sort(unique(tested$gene[!is.na(tested$gene)]))
  fwrite(data.table(gene = sig_genes), file.path(go_dir, paste0(prefix, "_input_genes.tsv")), sep = "\t")
  fwrite(data.table(gene = universe), file.path(go_dir, paste0(prefix, "_universe_genes.tsv")), sep = "\t")

  mapped_sig <- AnnotationDbi::select(
    org.Hs.eg.db, keys = sig_genes, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")
  )
  mapped_universe <- AnnotationDbi::select(
    org.Hs.eg.db, keys = universe, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")
  )
  mapped_sig <- unique(mapped_sig[!is.na(mapped_sig$ENTREZID), c("SYMBOL", "ENTREZID")])
  mapped_universe <- unique(mapped_universe[!is.na(mapped_universe$ENTREZID), c("SYMBOL", "ENTREZID")])
  fwrite(mapped_sig, file.path(go_dir, paste0(prefix, "_mapped_genes.tsv")), sep = "\t")
  fwrite(mapped_universe, file.path(go_dir, paste0(prefix, "_mapped_universe.tsv")), sep = "\t")

  if (!nrow(mapped_sig) || !nrow(mapped_universe)) {
    fwrite(data.table(), file.path(go_dir, paste0(prefix, "_GO_BP.tsv.gz")), sep = "\t")
    return(data.table(method = label, n_pairs = nrow(significant), n_genes = length(sig_genes),
                      mapped_genes = nrow(mapped_sig), n_terms = 0L))
  }

  ego <- enrichGO(
    gene = unique(mapped_sig$ENTREZID),
    universe = unique(mapped_universe$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  go_tab <- as.data.table(as.data.frame(ego))
  if (nrow(go_tab)) {
    setorder(go_tab, p.adjust, pvalue)
    fwrite(go_tab, file.path(go_dir, paste0(prefix, "_GO_BP.tsv.gz")), sep = "\t")
    top <- go_tab[order(p.adjust, pvalue)][seq_len(min(20L, .N))]
    top[, Description := factor(Description, levels = rev(Description))]
    p <- ggplot(top, aes(x = -log10(p.adjust), y = Description, size = Count, colour = -log10(p.adjust))) +
      geom_point() +
      scale_colour_viridis_c() +
      labs(title = paste0(label, " direct GO biological-process enrichment"),
           x = expression(-log[10](BH~adjusted~p)), y = NULL, size = "Genes") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
    ggsave(file.path(go_dir, paste0(prefix, "_GO_BP.png")), p, width = 8, height = 6, dpi = 320)
    ggsave(file.path(go_dir, paste0(prefix, "_GO_BP.pdf")), p, width = 8, height = 6, device = cairo_pdf)
  } else {
    fwrite(data.table(), file.path(go_dir, paste0(prefix, "_GO_BP.tsv.gz")), sep = "\t")
  }
  data.table(method = label, n_pairs = nrow(significant), n_genes = length(sig_genes),
            mapped_genes = nrow(mapped_sig), n_terms = nrow(go_tab))
}

go_summary <- rbind(
  run_go("EpiLink-Diff Poisson", epilink_sig, epilink_all, "epilink_diff_poisson"),
  run_go("scMultiMap-Diff", sc_sig, sc_all, "scmultimap_diff")
)
fwrite(go_summary, file.path(go_dir, "GO_analysis_summary.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# 2. Yang Supplementary Table 6, GRCh38, and Brain2 peak/gene overlap
# -------------------------------------------------------------------------

read_yang <- function(sheet, names_out) {
  x <- as.data.table(read_excel(yang_xlsx, sheet = sheet, skip = 1))
  if (length(names_out) != ncol(x)) stop("Unexpected columns in ", sheet)
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

yang6c[, chr := as.integer(chr)]
yang6c[, pos_hg38 := as.integer(pos_hg38)]
yang6d[, chr := as.integer(chr)]
yang6d[, pos_hg38 := as.integer(pos_hg38)]
yang6c[, chr_grch38 := paste0("chr", chr)]
yang6c[, genome_build := "GRCh38"]
yang6d[, chr_grch38 := paste0("chr", chr)]
yang6d[, genome_build := "GRCh38"]
fwrite(yang6c, file.path(overlap_dir, "yang_supp_table6c_grch38.tsv.gz"), sep = "\t")
fwrite(yang6d, file.path(overlap_dir, "yang_supp_table6d_grch38.tsv.gz"), sep = "\t")

target_long <- yang6d[!is.na(pos_hg38) & !is.na(target_genes), .(
  gwas_locus, rsid, chr_grch38, pos_hg38, target_gene = unlist(strsplit(gsub("\\s+", "", target_genes), ";", fixed = TRUE)),
  target_rpkm
)]
target_long <- target_long[target_gene != "." & nzchar(target_gene)]

candidate <- readRDS(file.path(epilink_dir, "candidate_pairs.rds"))
candidate <- as.data.table(candidate)[, .(gene = as.character(gene), peak = as.character(peak))]
peak_parts <- tstrsplit(candidate$peak, "-", fixed = TRUE)
candidate[, `:=`(chr_grch38 = peak_parts[[1L]], peak_start = as.integer(peak_parts[[2L]]), peak_end = as.integer(peak_parts[[3L]]))]
setkey(candidate, chr_grch38, peak_start, peak_end)
variant_intervals <- target_long[, .(chr_grch38, variant_start = pos_hg38, variant_end = pos_hg38,
                                     gwas_locus, rsid, pos_hg38, target_gene, target_rpkm)]

# foverlaps avoids a large chromosome-wise Cartesian merge.
overlap_all <- foverlaps(
  variant_intervals, candidate,
  by.x = c("chr_grch38", "variant_start", "variant_end"),
  by.y = c("chr_grch38", "peak_start", "peak_end"),
  type = "within", nomatch = 0L, allow.cartesian = TRUE
)
setnames(overlap_all, c("start", "end"), c("peak_start", "peak_end"), skip_absent = TRUE)
overlap_all[, gene_match := target_gene == gene]
overlap_all[, peak_gene_match := gene_match]
overlap_gene <- unique(overlap_all[gene_match == TRUE])
fwrite(overlap_all, file.path(overlap_dir, "yang_brain2_peak_overlaps.tsv.gz"), sep = "\t")
fwrite(overlap_gene, file.path(overlap_dir, "yang_brain2_peak_gene_matches.tsv.gz"), sep = "\t")

ep_cols <- intersect(c("gene", "peak", "beta_control", "beta_disease", "beta_diff", "wald_p", "padj", "calibrated_q", "significant"), names(epilink_all))
sc_cols <- intersect(c("gene", "peak", "covar_control", "covar_disease", "cor_control", "cor_disease", "pval", "padj"), names(sc_all))
ep_hit <- merge(overlap_gene, epilink_all[, ..ep_cols], by = c("gene", "peak"), all.x = TRUE)
sc_hit <- merge(overlap_gene, sc_all[, ..sc_cols], by = c("gene", "peak"), all.x = TRUE)
ep_sig_hit <- ep_hit[(isTRUE(significant) | (!is.na(significant) & significant == TRUE)) |
                       (!is.na(calibrated_q) & calibrated_q <= 0.1)]
fwrite(ep_hit, file.path(overlap_dir, "yang_epilink_all_peak_gene_matches.tsv.gz"), sep = "\t")
fwrite(ep_sig_hit, file.path(overlap_dir, "yang_epilinkdiff_significant_matches.tsv.gz"), sep = "\t")
fwrite(sc_hit, file.path(overlap_dir, "yang_scmultimap_all_peak_gene_matches.tsv.gz"), sep = "\t")

# -------------------------------------------------------------------------
# 3. Fig. 4d/e-style plotting for genuine matches
# -------------------------------------------------------------------------

plot_match <- function(hit, index) {
  pair <- unique(hit[, .(gene, peak, chr_grch38, peak_start, peak_end, rsid, pos_hg38, gwas_locus)])
  pair <- pair[1L]
  ep <- epilink_all[gene == pair$gene & peak == pair$peak][1L]
  sc <- sc_all[gene == pair$gene & peak == pair$peak][1L]
  x0 <- min(pair$peak_start, pair$pos_hg38) - 5000
  x1 <- max(pair$peak_end, pair$pos_hg38) + 5000
  locus <- data.table(x = c(pair$peak_start, pair$peak_end), y = 1)
  variant <- data.table(x = pair$pos_hg38, label = pair$rsid)
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
  beta <- data.table(condition = factor(c("Ctrl", "AD"), levels = c("Ctrl", "AD")),
                     beta = c(ep$beta_control, ep$beta_disease))
  p2 <- ggplot(beta, aes(x = condition, y = beta, colour = condition)) +
    geom_hline(yintercept = 0, colour = "grey60") + geom_point(size = 3) + geom_line(aes(group = 1), linewidth = 0.7) +
    scale_colour_manual(values = c(Ctrl = "#2B8CBE", AD = "#D7301F"), guide = "none") +
    labs(title = "B  EpiLink condition-specific association", y = "ATAC slope (GLM linear-predictor scale)", x = NULL,
         subtitle = sprintf("beta_diff = %+.3f; Wald p = %.3g", ep$beta_diff, ep$wald_p)) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  p3 <- ggplot(data.table(condition = factor(c("Ctrl", "AD"), levels = c("Ctrl", "AD")),
                          correlation = c(sc$cor_control, sc$cor_disease)), aes(condition, correlation, colour = condition)) +
    geom_hline(yintercept = 0, colour = "grey60") + geom_point(size = 3) + geom_line(aes(group = 1), linewidth = 0.7) +
    scale_colour_manual(values = c(Ctrl = "#2B8CBE", AD = "#D7301F"), guide = "none") +
    labs(title = "C  scMultiMap association", y = "Latent correlation", x = NULL,
         subtitle = sprintf("Differential p = %.3g", sc$pval)) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  fig <- p1 / (p2 | p3) + plot_annotation(title = paste0(pair$gene, " GWAS-supported peak–gene association"))
  stem <- sprintf("GWAS_match_%02d_%s_%s", index, pair$gene, gsub("[^A-Za-z0-9]+", "_", pair$rsid))
  ggsave(file.path(figure_dir, paste0(stem, ".png")), fig, width = 9, height = 7, dpi = 320)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), fig, width = 9, height = 7, device = cairo_pdf)
}

if (nrow(ep_sig_hit)) {
  unique_hits <- unique(ep_sig_hit[, .(gene, peak, chr_grch38, peak_start, peak_end, rsid, pos_hg38, gwas_locus)])
  for (i in seq_len(min(nrow(unique_hits), 10L))) plot_match(unique_hits[i], i)
  figure_status <- paste0("Generated ", min(nrow(unique_hits), 10L), " Fig. 4d/e-style figure(s).")
} else {
  figure_status <- "No exact Yang-variant/Brain2-peak/Yang-target/EpiLinkDiff-significant match; no positive Fig. 4d/e-style figure was generated."
  writeLines(figure_status, file.path(figure_dir, "NO_FIGURE_EXACT_EpiLinkDiff_MATCH.txt"))
}

summary <- data.table(
  quantity = c("EpiLinkDiff significant pairs", "EpiLinkDiff significant genes",
               "scMultiMap-Diff significant pairs", "scMultiMap-Diff significant genes",
               "Yang Table 6d rows", "Yang Table 6d unique variants", "Yang Table 6d target genes",
               "Brain2 candidate pairs", "variant-containing Brain2 pairs",
               "variant-containing pairs with Yang target-gene match",
               "EpiLinkDiff significant exact matches"),
  value = c(nrow(epilink_sig), uniqueN(epilink_sig$gene), nrow(sc_sig), uniqueN(sc_sig$gene),
            nrow(yang6d), uniqueN(yang6d$rsid), uniqueN(target_long$target_gene), nrow(candidate),
            uniqueN(overlap_all[, .(gene, peak)]), uniqueN(overlap_gene[, .(gene, peak)]), uniqueN(ep_sig_hit[, .(gene, peak)]))
)
fwrite(summary, file.path(out_dir, "analysis_summary.tsv"), sep = "\t")
writeLines(c(
  "Brain2 Microglia direct GO and Yang Table 6 analysis",
  paste0("Output directory: ", normalizePath(out_dir)),
  paste0("Figure status: ", figure_status),
  "Yang Table 6d pos hg38 was used directly and marked GRCh38; no liftover was applied."
), file.path(out_dir, "README_analysis.txt"))
message(figure_status)
message("Wrote GO results to ", normalizePath(go_dir))
message("Wrote Yang overlap results to ", normalizePath(overlap_dir))
