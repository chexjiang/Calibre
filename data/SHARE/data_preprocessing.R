### SHARE-seq

library(Signac)
library(Seurat)
library(harmony)
library(scMultiMap)
library(data.table)
library(stringr)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Signac)
library(GenomeInfoDb)

# annotation <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75))
# seqlevelsStyle(annotation) <- "UCSC"
# saveRDS(annotation, "data/hg19_EnsDb_Hsapiens_v75_annotation.rds")

# Annotation
annotation <- readRDS("~/proj/EpiLink/data/hg19_EnsDb_Hsapiens_v75_annotation.rds")


# Read sc-multiome data
input_atac <- "data/SHARE/SHARE.rep3.atac.rds"
input_mrna <- "data/SHARE/SHARE.rep3.rna.rds"
input_meta <- "data/SHARE/SHARE.rep3.meta.rds"


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
gene_peak_file <- "data/SHARE/SHARE.rep3.gene_peak_df.rds"
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


# Example: color by cell type annotation
DimPlot(object, reduction = "umap_rna", group.by = "predicted.id")
DimPlot(object, reduction = "umap_rna", group.by = "newCT")
