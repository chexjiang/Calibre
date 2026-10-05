#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(harmony)
  library(GenomicRanges)
  library(reticulate)
  library(scMultiMap)
  library(ggplot2)
})

# Run this script from the repository root.
root <- normalizePath(getwd())
out <- file.path(root, "data", "bmmc")
annotation <- readRDS(file.path(root, "data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
source(file.path(root, "R", "helper.R"))

rna_file <- file.path(out, "NeurIPS_multiome.rna.rds")
atac_file <- file.path(out, "NeurIPS_multiome.atac.rds")
meta_file <- file.path(out, "NeurIPS_multiome.meta.rds")

if (all(file.exists(rna_file, atac_file, meta_file))) {
  rna <- readRDS(rna_file)
  atac <- readRDS(atac_file)
  meta <- readRDS(meta_file)
} else {
  use_python(Sys.getenv("RETICULATE_PYTHON", "/home/chenxin/miniforge3/bin/python"),
             required = TRUE)
  anndata <- import("anndata", convert = FALSE)
  adata <- anndata$read_h5ad(file.path(root, "data", "bmmc", "GSE194122",
                                       "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"))
  meta <- py_to_r(adata$obs)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  cells <- rownames(meta)
  feature_types <- as.character(py_to_r(adata$var[["feature_types"]]))
  feature_names <- rownames(py_to_r(adata$var))
  counts <- adata$layers[["counts"]]

  read_assay <- function(assay) {
    keep <- which(feature_types == assay)
    py_keep <- r_to_py(as.integer(keep - 1L))
    xpy <- counts[, py_keep]$T$tocsc()
    x <- new("dgCMatrix",
             i = as.integer(py_to_r(xpy$indices)),
             p = as.integer(py_to_r(xpy$indptr)),
             x = as.numeric(py_to_r(xpy$data)),
             Dim = as.integer(py_to_r(xpy$shape)))
    rownames(x) <- feature_names[keep]
    colnames(x) <- cells
    if (assay == "ATAC") {
      rownames(x) <- sub("^([^:]+)-([0-9]+)-([0-9]+)$", "\\1:\\2-\\3", rownames(x))
    }
    x
  }
  rna <- read_assay("GEX")
  atac <- read_assay("ATAC")
  saveRDS(rna, rna_file)
  saveRDS(atac, atac_file)
  saveRDS(meta, meta_file)
}

stopifnot(identical(colnames(rna), colnames(atac)))
stopifnot(identical(colnames(rna), rownames(meta)))

# Save the matrices and metadata used by EpiLink.

# The h5ad already contains a fixed peak set; only derive genomic candidate
# gene-peak pairs on that set.  No CallPeaks step is performed here.
object <- CreateSeuratObject(counts = rna, meta.data = meta, assay = "RNA")
object[["ATAC"]] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  annotation = annotation
)
gene_peak_df <- GetGenePeakPairs(object, annotation)
# Keep the same chr-start-end convention as the Open Problems feature IDs and
# the existing EpiLink inputs.
gene_peak_df$peak <- sub(":", "-", gene_peak_df$peak, fixed = TRUE)
dimnames(atac)[[1]] <- sub(":", "-", dimnames(atac)[[1]], fixed = TRUE)
saveRDS(atac, atac_file)
saveRDS(gene_peak_df, file.path(out, "NeurIPS_multiome.gene_peak_df.rds"))

cat("Saved", ncol(rna), "cells;", nrow(rna), "RNA features;",
    nrow(atac), "ATAC peaks;", nrow(gene_peak_df), "gene-peak pairs\n")



# RNA preprocessing + UMAP
DefaultAssay(object) <- "RNA"

object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 50)
# UMAP before Harmony
object <- RunUMAP(object, reduction = "pca", dims = 1:30, reduction.name = "umap_pca")
DimPlot(object, reduction = "umap_pca", group.by = "batch") + ggtitle("Luecken BMMC with batch")

# Run Harmony
object <- RunHarmony(object, group.by.vars = "batch", reduction = "pca", dims.use = 1:30, reduction.save = "harmony")
object <- RunUMAP(object, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
table(object$cell_type)
# Cell type annotation
DimPlot(object, reduction = "umap_harmony", group.by = "cell_type")
# Rename newCT
object$newCT <- "Others"
object$newCT[object$cell_type %in% c("B1 B", "Naive CD20+ B", "Transitional B")] <- "Bcell"
object$newCT[object$cell_type %in% c("CD4+ T activated", "CD4+ T naive", "CD8+ T", "CD8+ T naive")] <- "T"
object$newCT[object$cell_type %in% c("CD14+ Mono", "CD16+ Mono", "ID2-hi myeloid prog")] <- "Myeloid"
object$newCT[object$cell_type %in% c("Erythroblast", "Proerythroblast", "Normoblast")] <- "Erythroid"
object$newCT[object$cell_type %in% c("cDC2", "pDC")] <- "Dendritic"
table(object$newCT, useNA = "ifany")

# Bcell Dendritic Erythroid   Myeloid    Others         T 
# 9752      2050      8996     12845     13081     22525 

cell_counts <- table(object$newCT)
legend_labels <- paste0(names(cell_counts), " (", as.integer(cell_counts), ")")

DimPlot(object, reduction = "umap_harmony", group.by = "newCT") +
  scale_color_discrete(labels = legend_labels) + ggtitle("Luecken BMMC")

meta <- meta[colnames(object), , drop = FALSE]
meta$newCT <- as.character(object$newCT)
saveRDS(meta, meta_file)

DimPlot(object, reduction = "umap_harmony", group.by = "batch") + ggtitle("Luecken BMMC")
saveRDS(object, file = file.path(out, "NeurIPS_object.rds"))
