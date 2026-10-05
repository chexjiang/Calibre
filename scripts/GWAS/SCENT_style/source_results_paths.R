# Canonical paths for the UKBB SCENT-style source-result layer.
# Keep all path construction in one place so benchmark and preparation scripts
# agree on the pooled/celltype layout.

source_results_paths <- function(root) {
  gwas <- file.path(root, "results", "validation", "GWAS")
  sr <- file.path(gwas, "SCENT_style", "UKBB", "source_results")
  list(
    root = sr,
    pooled_main = file.path(sr, "pooled", "main_datasets"),
    celltype_main = file.path(sr, "celltype", "main_datasets"),
    pooled_brain2 = file.path(sr, "pooled", "brain2_pituitary"),
    celltype_brain2 = file.path(sr, "celltype", "brain2_pituitary"),
    pooled = file.path(sr, "pooled"),
    celltype = file.path(sr, "celltype"),
    pooled_annotations = file.path(sr, "pooled", "annotations"),
    celltype_annotations = file.path(sr, "celltype", "annotations"),
    pooled_results = file.path(sr, "pooled", "results"),
    celltype_results = file.path(sr, "celltype", "results"),
    shared = file.path(sr, "shared"),
    shared_brain2 = file.path(sr, "shared", "brain2_pituitary")
  )
}
