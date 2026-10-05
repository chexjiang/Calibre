### NEAT-seq

library(Signac)
library(Seurat)
library(harmony)
library(scMultiMap)
library(data.table)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(Signac)
library(GenomeInfoDb)


# Annotation
annotation <- readRDS("~/proj/EpiLink/data/hg38_EnsDb_Hsapiens_v86_annotation.rds")

# Read sc-multiome data
input_atac <- "data/NEAT/NEAT.atac.rds"
input_mrna <- "data/NEAT/NEAT.rna.rds"
input_meta <- "data/NEAT/NEAT.meta.rds"

atac <- readRDS(input_atac)
mrna <- readRDS(input_mrna)
meta <- readRDS(input_meta)

# Ensure metadata matches RNA/ATAC
stopifnot(all(colnames(mrna) == colnames(atac)))   # Same cell order
stopifnot(all(rownames(meta) %in% colnames(mrna))) # Metadata cells exist
meta <- meta[colnames(mrna), , drop = FALSE]

# Construct Seurat object
object <- CreateSeuratObject(counts = mrna, meta.data = meta, assay = "RNA")
object[["ATAC"]] <- CreateChromatinAssay(counts = atac)

# Get candidate gene-peak pairs
gene_peak_file <- "data/NEAT/NEAT.gene_peak_df.rds"
if (!file.exists(gene_peak_file)) {
  gene_peak_df <- GetGenePeakPairs(object, annotation)
  saveRDS(gene_peak_df, file = gene_peak_file)
}




# RNA preprocessing + UMAP
DefaultAssay(object) <- "RNA"

object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 50)
object <- RunUMAP(object, reduction = "pca", dims = 1:30, reduction.name = "umap_rna")

object$newCT

# Example: color by cell type annotation
DimPlot(object, reduction = "umap_rna", group.by = "newCT")
table(object$newCT)
# CD4T 
# 8472 

# Example: color by sample lane
DimPlot(object, reduction = "umap_rna", group.by = "Sample")
table(object$Sample)
# lane1 lane2 
# 4511  3961 
