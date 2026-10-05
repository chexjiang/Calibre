#!/usr/bin/env Rscript

# Direct GO-BP enrichment for Brain2 Microglia EpiLink-Diff and scMultiMap-Diff
# significant genes. The universe is all genes in the corresponding candidate
# pair result table, so enrichment is conditional on the tested gene set.

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_dir)

out_dir <- file.path(project_dir, "results", "validation", "functional", "Brain2_diff", "Microglia", "go")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ep_dir <- file.path(project_dir, "results", "real", "Brain2_diff", "Microglia", "epilink_diff_poisson_no_donor_effect_no_donor_robust")
sc_dir <- file.path(project_dir, "results", "real", "Brain2_diff", "Microglia", "scmultimap")

run_go <- function(label, sig_path, all_path, prefix) {
  sig <- fread(sig_path, showProgress = FALSE)
  tested <- fread(all_path, showProgress = FALSE)
  sig_genes <- sort(unique(sig$gene[!is.na(sig$gene)]))
  universe <- sort(unique(tested$gene[!is.na(tested$gene)]))
  fwrite(data.table(gene = sig_genes), file.path(out_dir, paste0(prefix, "_input_genes.tsv")), sep = "\t")
  fwrite(data.table(gene = universe), file.path(out_dir, paste0(prefix, "_universe_genes.tsv")), sep = "\t")

  mapped_sig <- AnnotationDbi::select(org.Hs.eg.db, keys = sig_genes, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
  mapped_universe <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
  mapped_sig <- unique(mapped_sig[!is.na(mapped_sig$ENTREZID), c("SYMBOL", "ENTREZID")])
  mapped_universe <- unique(mapped_universe[!is.na(mapped_universe$ENTREZID), c("SYMBOL", "ENTREZID")])
  fwrite(mapped_sig, file.path(out_dir, paste0(prefix, "_mapped_genes.tsv")), sep = "\t")
  fwrite(mapped_universe, file.path(out_dir, paste0(prefix, "_mapped_universe.tsv")), sep = "\t")

  ego <- enrichGO(
    gene = unique(mapped_sig$ENTREZID), universe = unique(mapped_universe$ENTREZID),
    OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
  )
  tab <- as.data.table(as.data.frame(ego))
  if (nrow(tab)) {
    setorder(tab, p.adjust, pvalue)
    fwrite(tab, file.path(out_dir, paste0(prefix, "_GO_BP.tsv.gz")), sep = "\t")
    top <- tab[seq_len(min(20L, .N))]
    top[, Description := factor(Description, levels = rev(Description))]
    p <- ggplot(top, aes(-log10(p.adjust), Description, size = Count, colour = -log10(p.adjust))) +
      geom_point() + scale_colour_viridis_c() +
      labs(title = paste0(label, " direct GO biological-process enrichment"),
           x = expression(-log[10](BH~adjusted~p)), y = NULL, size = "Genes") +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
    ggsave(file.path(out_dir, paste0(prefix, "_GO_BP.png")), p, width = 8, height = 6, dpi = 320)
    ggsave(file.path(out_dir, paste0(prefix, "_GO_BP.pdf")), p, width = 8, height = 6, device = cairo_pdf)
  } else {
    con <- gzfile(file.path(out_dir, paste0(prefix, "_GO_BP.tsv.gz")), "wt"); close(con)
  }
  data.table(method = label, n_pairs = nrow(sig), n_genes = length(sig_genes), mapped_genes = nrow(mapped_sig), n_terms = nrow(tab))
}

summary <- rbind(
  run_go("EpiLink-Diff Poisson", file.path(ep_dir, "significant_pairs.tsv.gz"), file.path(ep_dir, "all_pairs.tsv.gz"), "epilink_diff_poisson"),
  run_go("scMultiMap-Diff", file.path(sc_dir, "significant_pairs.tsv.gz"), file.path(sc_dir, "all_pairs.tsv.gz"), "scmultimap_diff")
)
fwrite(summary, file.path(out_dir, "GO_analysis_summary.tsv"), sep = "\t")
message("Wrote direct GO results to ", normalizePath(out_dir))
