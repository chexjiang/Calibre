#!/usr/bin/env Rscript

# Export source-annotation counts together with UKBB overlap counts for the
# requested pooled and cell-type-specific analyses.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork=TRUE)
scope_name <- if (length(args) >= 2L) args[[2]] else "custom"
gwas <- file.path(root, "results", "validation", "GWAS")
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
source_results <- source_results_paths(root)
pooled_base <- source_results$pooled
cell_base <- source_results$celltype
out <- file.path(source_results$shared, paste0("ukbb_method_count_table_", scope_name, ".tsv"))

one <- function(base, analysis) {
  m <- fread(file.path(base, "manifest.tsv"))
  s <- fread(file.path(base, "results", "scent_sakaue_summary.tsv"))[source == "UKBB"]
  p <- fread(file.path(base, "results", "scent_sakaue_per_trait.tsv.gz"))[source == "UKBB"]
  counts <- p[, .(
    n_causal_in_sum = sum(n_causal_in_annotation, na.rm=TRUE),
    n_causal_total_sum = sum(n_causal_total, na.rm=TRUE),
    n_common_in_sum = sum(n_common_in_annotation, na.rm=TRUE),
    n_common_locus_total_sum = sum(n_common_locus_total, na.rm=TRUE)
  ), by=annotation_id]
  ans <- merge(m[, .(annotation_id, dataset, celltype, celltypes, method, method_label,
                     n_pairs, n_hg38_intervals)],
               s[, .(annotation_id, mean_enrichment, ci_low, ci_high, pooled_enrichment)],
               by="annotation_id", all.x=TRUE)
  ans <- merge(ans, counts, by="annotation_id", all.x=TRUE)
  ans[, analysis := analysis]
  setcolorder(ans, c("analysis", "dataset", "celltype", "celltypes", "method",
                     "method_label", "n_pairs", "n_hg38_intervals",
                     "n_causal_in_sum", "n_causal_total_sum", "n_common_in_sum",
                     "n_common_locus_total_sum", "mean_enrichment", "ci_low",
                     "ci_high", "pooled_enrichment", "annotation_id"))
  ans
}

out_dt <- rbindlist(list(one(pooled_base, "pooled"), one(cell_base, "celltype")), fill=TRUE)
fwrite(out_dt, out, sep="\t")
message("Wrote ", nrow(out_dt), " count rows to ", out)
