#!/usr/bin/env Rscript

# One-time migration helper for merging the main-dataset and Brain2/Pituitary
# UKBB SCENT-style results within each analysis view. The current canonical
# layout has already completed this migration: BEDs live in the view-level
# annotations/ directory and provenance is retained in annotation_source.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
paths <- source_results_paths(root)

read_tagged <- function(path, tag) {
  if (!file.exists(path)) stop("Missing source result: ", path)
  x <- fread(path)
  x[, annotation_source := tag]
  x
}

merge_view <- function(view, main_dir, brain_dir, output_dir,
                       brain_result_name, brain_summary_name) {
  main_results <- file.path(main_dir, "results")
  brain_results <- file.path(brain_dir, "results")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "results"), recursive = TRUE, showWarnings = FALSE)

  main_per <- read_tagged(file.path(main_results, "scent_sakaue_per_trait.tsv.gz"),
                          "main_datasets")
  brain_per <- read_tagged(file.path(brain_results, brain_result_name),
                           "brain2_pituitary")
  per <- rbindlist(list(main_per, brain_per), fill = TRUE, use.names = TRUE)
  setorder(per, dataset, celltype_scope, method_label, trait, annotation_id)
  fwrite(per, file.path(output_dir, "results", "scent_sakaue_per_trait.tsv.gz"),
         sep = "\t")

  main_summary <- read_tagged(file.path(main_results, "scent_sakaue_summary.tsv"),
                              "main_datasets")
  brain_summary <- read_tagged(file.path(brain_results, brain_summary_name),
                               "brain2_pituitary")
  summary <- rbindlist(list(main_summary, brain_summary), fill = TRUE, use.names = TRUE)
  setorder(summary, dataset, celltype_scope, method_label, annotation_id)
  fwrite(summary, file.path(output_dir, "results", "scent_sakaue_summary.tsv"),
         sep = "\t")

  main_manifest <- fread(file.path(main_dir, "manifest.tsv"))
  main_manifest[, annotation_source := "main_datasets"]
  brain_manifest <- fread(file.path(brain_dir, "manifest.tsv"))
  brain_manifest[, annotation_source := "brain2_pituitary"]
  manifest <- rbindlist(list(main_manifest, brain_manifest), fill = TRUE, use.names = TRUE)
  setorder(manifest, dataset, celltype, method, annotation_id)
  fwrite(manifest, file.path(output_dir, "manifest.tsv"), sep = "\t")

  # Recompute the audit from the merged per-trait rows.  The two source audits
  # cannot simply be concatenated because they summarize different subsets.
  audit <- per[, .(
    n_annotations = uniqueN(annotation_id),
    n_traits = uniqueN(trait),
    zero_causal_hits = sum(n_causal_in_annotation == 0L, na.rm = TRUE),
    zero_annotation_hits = sum(n_common_in_annotation == 0L, na.rm = TRUE),
    min_common_locus = min(n_common_locus_total, na.rm = TRUE),
    max_common_locus = max(n_common_locus_total, na.rm = TRUE)
  ), by = source]
  fwrite(audit, file.path(output_dir, "results", "scent_sakaue_audit.tsv"),
         sep = "\t")

  message(view, ": merged per-trait rows=", nrow(per),
          "; annotations=", nrow(manifest),
          "; output=", output_dir)
}

merge_view(
  "pooled", paths$pooled_main, paths$pooled_brain2, paths$pooled,
  "brain2_pituitary_pooled_per_trait.tsv.gz",
  "brain2_pituitary_pooled_summary.tsv"
)
merge_view(
  "celltype", paths$celltype_main, paths$celltype_brain2, paths$celltype,
  "brain2_pituitary_celltype_per_trait.tsv.gz",
  "brain2_pituitary_celltype_summary.tsv"
)
