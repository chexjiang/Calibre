#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
base <- file.path(root, "results", "validation", "GWAS", "sLDSC")
ldscore_dir <- Sys.getenv("SLDSC_LDSCORE_DIR", file.path(base, "ldscores"))
results_dir <- Sys.getenv("SLDSC_RESULTS_DIR", file.path(base, "results"))
manifest <- fread(file.path(ldscore_dir, "sldsc_annotation_manifest.tsv"))
traits <- fread(file.path(root, "results", "validation", "GWAS", "SCENT_style", "config", "traits.tsv"))

rows <- list()
k <- 0L
for (i in seq_len(nrow(manifest))) {
  a <- manifest[i]
  for (tr in unique(traits[panel == a$panel, trait])) {
    f <- file.path(results_dir, tr, paste0(a$annotation_id, ".results"))
    if (!file.exists(f)) next
    x <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(x) || !nrow(x)) next
    target <- x[grepl(a$annotation_id, Category, fixed = TRUE)]
    if (!nrow(target)) target <- x[.N]
    target <- target[1]
    k <- k + 1L
    rows[[k]] <- cbind(
      a[, .(annotation_id, dataset, celltype, panel, interpretation, method, method_label)],
      traits[trait == tr][1, .(trait, trait_label = label, phenotype, cases, controls)],
      target
    )
  }
}
ans <- rbindlist(rows, fill = TRUE)
if (nrow(ans)) {
  pcol <- intersect(names(ans), c("Enrichment_p", "Enrichment.P", "Coefficient_p"))
  if (length(pcol)) ans[, enrichment_fdr := p.adjust(get(pcol[[1]]), method = "BH")]
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(ans, file.path(results_dir, "sldsc_all_results.tsv"), sep = "\t")
message("Collected ", nrow(ans), " sLDSC annotation-trait results")
