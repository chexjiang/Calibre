# Canonical paths for the scMultiMap-style conditional S-LDSC benchmark.

sldsc_source_results_paths <- function(root) {
  base <- file.path(root, "results", "validation", "GWAS", "sLDSC")
  sr <- file.path(base, "source_results")
  list(
    base = base,
    root = sr,
    celltype = file.path(sr, "celltype"),
    pooled = file.path(sr, "pooled"),
    celltype_annotations = file.path(sr, "celltype", "annotations"),
    pooled_annotations = file.path(sr, "pooled", "annotations"),
    celltype_backgrounds = file.path(sr, "celltype", "backgrounds"),
    pooled_backgrounds = file.path(sr, "pooled", "backgrounds"),
    celltype_ldscores = file.path(sr, "celltype", "ldscores"),
    pooled_ldscores = file.path(sr, "pooled", "ldscores"),
    celltype_results = file.path(sr, "celltype", "results"),
    pooled_results = file.path(sr, "pooled", "results"),
    figures = file.path(base, "figures", "scmultimap_style"),
    tables = file.path(base, "tables", "scmultimap_style")
  )
}
