#!/usr/bin/env Rscript

# Create the canonical UKBB SCENT-style source-result layout.  This migration
# is copy-based: the previous sakaue_style tree is left intact until the new
# paths have been checked by the user.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
paths <- source_results_paths(root)

old <- file.path(root, "results", "validation", "GWAS", "SCENT_style", "UKBB", "sakaue_style")
old_pooled <- file.path(old, "method_comparison_custom")
old_celltype <- file.path(old, "celltype_comparison_custom")
old_brain <- file.path(old, "brain2_pituitary_ukbb")

required <- c(
  file.path(old_pooled, "results", "scent_sakaue_per_trait.tsv.gz"),
  file.path(old_celltype, "results", "scent_sakaue_per_trait.tsv.gz"),
  file.path(old_brain, "results", "scent_sakaue_per_trait.tsv.gz"),
  file.path(old_brain, "sakaue_annotation_manifest.tsv"),
  file.path(old_brain, "dataset_scope.tsv")
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing source files:\n", paste(missing, collapse = "\n"))

copy_one <- function(src, dst) {
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(src, dst, overwrite = TRUE)) stop("Could not copy: ", src)
}

clear_files <- function(dir) {
  if (!dir.exists(dir)) return(invisible(NULL))
  old_files <- list.files(dir, full.names = TRUE, recursive = FALSE)
  if (length(old_files)) unlink(old_files, recursive = TRUE, force = TRUE)
}

copy_dir_files <- function(src_dir, dst_dir) {
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src_dir, full.names = TRUE, recursive = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  if (length(files)) for (src in files) copy_one(src, file.path(dst_dir, basename(src)))
}

copy_analysis <- function(src, dst) {
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  clear_files(file.path(dst, "annotations"))
  clear_files(file.path(dst, "results"))
  copy_dir_files(file.path(src, "annotations"), file.path(dst, "annotations"))
  copy_dir_files(file.path(src, "results"), file.path(dst, "results"))
  copy_one(file.path(src, "sakaue_annotation_manifest.tsv"), file.path(dst, "manifest.tsv"))
  if (file.exists(file.path(src, "dataset_scope.tsv")))
    copy_one(file.path(src, "dataset_scope.tsv"), file.path(dst, "dataset_scope.tsv"))
}

copy_analysis(old_pooled, paths$pooled_main)
copy_analysis(old_celltype, paths$celltype_main)

# Keep a shared copy of the original Brain2/Pituitary inputs before splitting.
dir.create(paths$shared_brain2, recursive = TRUE, showWarnings = FALSE)
copy_one(file.path(old_brain, "dataset_scope.tsv"),
         file.path(paths$shared_brain2, "dataset_scope.tsv"))
copy_one(file.path(old_brain, "sakaue_annotation_manifest.tsv"),
         file.path(paths$shared_brain2, "original_manifest.tsv"))
copy_one(file.path(old_brain, "results", "scent_sakaue_per_trait.tsv.gz"),
         file.path(paths$shared_brain2, "original_per_trait.tsv.gz"))
copy_one(file.path(old_brain, "results", "scent_sakaue_summary.tsv"),
         file.path(paths$shared_brain2, "original_summary.tsv"))
copy_one(file.path(old_brain, "results", "scent_sakaue_audit.tsv"),
         file.path(paths$shared_brain2, "original_audit.tsv"))

# The custom dataset scope is shared by pooled and celltype preparation.
copy_one(file.path(old, "custom_scope.tsv"),
         file.path(paths$shared, "custom_scope.tsv"))

brain_manifest <- fread(file.path(old_brain, "sakaue_annotation_manifest.tsv"))
brain_per_trait <- fread(file.path(old_brain, "results", "scent_sakaue_per_trait.tsv.gz"))
brain_summary <- fread(file.path(old_brain, "results", "scent_sakaue_summary.tsv"))

# The historical Brain2/Pituitary manifest used celltype_scope=
# "celltype_specific" even for the pooled ALL_REQUESTED rows.  Use the
# explicit celltype/annotation identifier instead of that inconsistent field.
is_pooled <- brain_manifest$celltype == "ALL_REQUESTED" |
  grepl("__ALL_REQUESTED__", brain_manifest$annotation_id, fixed = TRUE)
pooled_manifest <- brain_manifest[is_pooled]
celltype_manifest <- brain_manifest[!is_pooled]
pooled_ids <- pooled_manifest$annotation_id
celltype_ids <- celltype_manifest$annotation_id

copy_brain_annotations <- function(ids, dst) {
  clear_files(file.path(dst, "annotations"))
  dir.create(file.path(dst, "annotations"), recursive = TRUE, showWarnings = FALSE)
  for (id in ids) {
    src <- file.path(old_brain, "annotations", paste0(id, ".hg38.bed"))
    if (!file.exists(src)) stop("Missing Brain2/Pituitary annotation: ", src)
    copy_one(src, file.path(dst, "annotations", basename(src)))
  }
}

write_view <- function(manifest, per_trait, summary, dst, per_trait_name, summary_name) {
  clear_files(file.path(dst, "results"))
  dir.create(file.path(dst, "results"), recursive = TRUE, showWarnings = FALSE)
  fwrite(manifest, file.path(dst, "manifest.tsv"), sep = "\t")
  fwrite(per_trait, file.path(dst, "results", per_trait_name), sep = "\t")
  fwrite(summary, file.path(dst, "results", summary_name), sep = "\t")
}

pooled_per_trait <- brain_per_trait[annotation_id %chin% pooled_ids]
celltype_per_trait <- brain_per_trait[annotation_id %chin% celltype_ids]
pooled_summary <- brain_summary[annotation_id %chin% pooled_ids]
celltype_summary <- brain_summary[annotation_id %chin% celltype_ids]

write_view(
  pooled_manifest, pooled_per_trait, pooled_summary, paths$pooled_brain2,
  "brain2_pituitary_pooled_per_trait.tsv.gz",
  "brain2_pituitary_pooled_summary.tsv"
)
write_view(
  celltype_manifest, celltype_per_trait, celltype_summary, paths$celltype_brain2,
  "brain2_pituitary_celltype_per_trait.tsv.gz",
  "brain2_pituitary_celltype_summary.tsv"
)
copy_brain_annotations(pooled_ids, paths$pooled_brain2)
copy_brain_annotations(celltype_ids, paths$celltype_brain2)

message("Created canonical UKBB source-result layout at: ", paths$root)
message("Pooled Brain2/Pituitary rows: ", nrow(pooled_per_trait),
        "; celltype rows: ", nrow(celltype_per_trait))
