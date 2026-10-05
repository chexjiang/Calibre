### Pituitary

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
input_atac <- "data/Pituitary/Pituitary.atac.rds"
input_mrna <- "data/Pituitary/Pituitary.rna.rds"
input_meta <- "data/Pituitary/Pituitary.meta.rds"

atac <- readRDS(input_atac)
mrna <- readRDS(input_mrna)
meta <- readRDS(input_meta)
head(row.names(atac))

# Ensure metadata matches RNA/ATAC
stopifnot(all(colnames(mrna) == colnames(atac)))   # Same cell order
stopifnot(all(rownames(meta) %in% colnames(mrna))) # Metadata cells exist
meta <- meta[colnames(mrna), , drop = FALSE]

# Construct Seurat object
object <- CreateSeuratObject(counts = mrna, meta.data = meta, assay = "RNA")
object[["ATAC"]] <- CreateChromatinAssay(counts = atac)

# Get candidate gene-peak pairs
gene_peak_df <- GetGenePeakPairs(object, annotation)
gene_peak_df <- FilterGenePeakPairs(object, gene_peak_df, selection = "rate")
saveRDS(gene_peak_df, file = "data/Pituitary/Pituitary.gene_peak_df.rds")



# RNA preprocessing + UMAP
DefaultAssay(object) <- "RNA"

object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, npcs = 50)
object <- RunUMAP(object, reduction = "pca", dims = 1:30, reduction.name = "umap_rna")


# Example: color by cell type annotation
DimPlot(object, reduction = "umap_rna", group.by = "newCT")
table(object$newCT)
# Corticotropes     Endothelial    Gonadotropes    Immune cells     Lactotropes       Pericytes 
# 898              94            4767             120            1222              94 
# Pituicytes    Somatotropes      Stem cells Stem committing     Thyrotropes 
# 245            2617            1799             407             988 
names(which(table(object$newCT) > 500))
# [1] "Corticotropes" "Gonadotropes"  "Lactotropes"   "Somatotropes"  "Stem cells"    "Thyrotropes"  