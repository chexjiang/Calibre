# Finish the arthritis outputs after the large H5AD-to-RDS conversion.
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(Signac); library(ggplot2); library(rhdf5)
})
source(file.path("R", "helper.R"))
out_dir <- file.path("data", "Arthritis")
rna <- readRDS(file.path(out_dir, "Arthritis.rna.rds"))
meta <- readRDS(file.path(out_dir, "Arthritis.meta.rds"))
annotation <- readRDS(file.path("data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
obj <- SeuratObject::CreateSeuratObject(counts=rna, assay="RNA", meta.data=meta)

# RNA preprocessing + UMAP
DefaultAssay(object) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)
# UMAP before Harmony
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap_pca")
DimPlot(obj, reduction = "umap_pca", group.by = "donor") + ggtitle("Brain with batch")

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
