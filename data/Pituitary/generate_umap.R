#!/usr/bin/env Rscript

# Generate an RNA PCA/UMAP plot from the current Pituitary benchmark inputs.
# This script intentionally does not rebuild the RNA/ATAC matrices or candidate
# pairs; it only reads the existing RDS inputs and writes exploration artifacts.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

file_arg <- commandArgs(FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
project_root <- if (length(file_arg)) {
  script_path <- sub("^--file=", "", file_arg[[1L]])
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
if (!dir.exists(file.path(project_root, "data", "Pituitary"))) {
  stop("Could not locate data/Pituitary under project root: ", project_root,
       call. = FALSE)
}
setwd(project_root)

out_dir <- file.path(project_root, "data", "Pituitary")
rna <- readRDS(file.path(out_dir, "Pituitary.rna.rds"))
meta <- readRDS(file.path(out_dir, "Pituitary.meta.rds"))
stopifnot(identical(colnames(rna), rownames(meta)))
stopifnot("newCT" %in% names(meta))

celltype <- as.character(meta$newCT)
celltype_counts <- table(celltype, useNA = "no")
keep_types <- names(celltype_counts)[celltype_counts > 500]
keep_cells <- !is.na(celltype) & celltype %in% keep_types
rna <- rna[, keep_cells, drop = FALSE]
meta <- meta[keep_cells, , drop = FALSE]
meta$plot_celltype <- as.character(meta$newCT)
stopifnot(identical(colnames(rna), rownames(meta)))

object <- CreateSeuratObject(counts = rna, meta.data = meta, assay = "RNA")
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, nfeatures = 2000, verbose = FALSE)
object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
object <- RunPCA(object, features = VariableFeatures(object), npcs = 50,
                 verbose = FALSE)
object <- RunUMAP(object, reduction = "pca", dims = 1:30,
                  reduction.name = "umap_rna", n.neighbors = 30,
                  min.dist = 0.3, verbose = FALSE)

embedding <- as.data.frame(Embeddings(object, "umap_rna")[, 1:2, drop = FALSE])
colnames(embedding) <- c("UMAP1", "UMAP2")
embedding$celltype <- meta$plot_celltype[match(rownames(embedding), rownames(meta))]
label_counts <- table(embedding$celltype)
embedding$label <- paste0(embedding$celltype, " (",
                          label_counts[embedding$celltype], ")")

p <- ggplot(embedding, aes(UMAP1, UMAP2, colour = label)) +
  geom_point(size = 0.25, alpha = 0.55) +
  theme_classic(base_size = 11) +
  labs(title = "Pituitary",
       subtitle = "RNA PCA/UMAP from current benchmark input",
       colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(out_dir, "Pituitary.celltype_umap.png"), p,
       width = 8.5, height = 6.5, dpi = 220)

saveRDS(list(
  pca = Embeddings(object, "pca"),
  umap = embedding
), file.path(out_dir, "Pituitary.exploration_embedding.rds"), compress = "gzip")

cat("Wrote Pituitary UMAP for", nrow(embedding),
    "cells across", length(keep_types), "celltypes (n > 500).\n")
print(label_counts)
