#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork=TRUE)
view <- if (length(args) >= 2L) args[[2]] else "celltype"
source(file.path(root, "scripts", "GWAS", "sLDSC", "source_results_paths.R"))
paths <- sldsc_source_results_paths(root)
base <- paths[[view]]
manifest <- fread(file.path(base, "manifest.tsv"))
jobs <- fread(file.path(base, "ldscores", "sldsc_jobs.tsv"))
traits <- unique(fread(file.path(root, "results", "validation", "GWAS", "SCENT_style", "config", "traits.tsv"))[, .(trait, trait_label=label)])

rows <- list(); audit <- list(); k <- 0L; ka <- 0L
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i]
  a <- manifest[annotation_id == j$annotation_id][1]
  f <- file.path(base, "results", "raw", j$trait, paste0(j$annotation_id, ".results"))
  status <- "ok"; reason <- ""
  x <- if (file.exists(f)) tryCatch(fread(f), error=function(e) NULL) else NULL
  if (is.null(x) || !nrow(x)) {
    status <- "missing"; reason <- "LDSC result absent or unreadable"
  } else {
    target <- x[Category == j$annotation_id]
    if (!nrow(target)) target <- x[grepl(j$annotation_id, Category, fixed=TRUE)]
    if (!nrow(target)) {
      status <- "missing"; reason <- "target annotation row absent"
    } else {
      bg <- x[Category == j$background_id]
      if (!nrow(bg)) bg <- x[grepl(j$background_id, Category, fixed=TRUE)]
      k <- k + 1L
      rows[[k]] <- cbind(
        a,
        traits[trait == j$trait][1],
        target[1],
        data.table(background_coefficient=if (nrow(bg)) bg$Coefficient[1] else NA_real_,
                   background_coefficient_z=if (nrow(bg)) bg[["Coefficient_z-score"]][1] else NA_real_,
                   conditioning_model="baseline-LD v2.2 + all tested peaks + linked peaks",
                   gwas_reference="FinnGen R12")
      )
    }
  }
  ka <- ka + 1L
  audit[[ka]] <- data.table(annotation_id=j$annotation_id, background_id=j$background_id,
                            trait=j$trait, status=status, reason=reason, result_file=f)
}
ans <- rbindlist(rows, fill=TRUE)
if (nrow(ans)) {
  ans[, enrichment_fdr_global := p.adjust(Enrichment_p, method="BH")]
  ans[, enrichment_fdr_within_trait := p.adjust(Enrichment_p, method="BH"), by=trait]
}
dir.create(file.path(base, "results"), recursive=TRUE, showWarnings=FALSE)
fwrite(ans, file.path(base, "results", "sldsc_all_results.tsv"), sep="\t")
fwrite(rbindlist(audit, fill=TRUE), file.path(base, "results", "result_audit.tsv"), sep="\t")
message(view, ": collected ", nrow(ans), "/", nrow(jobs), " results")
