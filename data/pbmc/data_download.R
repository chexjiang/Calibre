# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi

library(Signac)
library(Seurat)
library(SeuratDisk)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("~/proj/EpiLink/data/pbmc")
# load the RNA and ATAC data
counts <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation)
pbmc

#------ QC
DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)
DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1)
pbmc
 

#------ Data processing    
DefaultAssay(pbmc) <- "RNA"
options(future.globals.maxSize = 2 * 1024^3)
pbmc <- SCTransform(pbmc, vst.flavor = "v1")
pbmc <- RunPCA(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")


DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
#saveRDS(pbmc, file = "pbmc.rds")

#------ Annotating cell types

# load PBMC reference from https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
#reference <- LoadH5Seurat("pbmc_multimodal.h5seurat", assays = list("SCT" = "counts"), reductions = 'spca')
reference <- readRDS("~/proj/EpiLink/data/pbmc/reference.rds")
reference <- UpdateSeuratObject(reference)

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"
table(pbmc@meta.data[["predicted.id"]])

# remove low-quality predictions
pbmc <- pbmc[, pbmc$prediction.score.max > 0.5]


#------ Joint UMAP visualization
# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend() + ggtitle("PBMC Celltype")


#------ Linking peaks to genes
DefaultAssay(pbmc) <- "ATAC"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = c("LYZ", "MS4A1")
)


links <- Links(pbmc)
links

idents.plot <- c("B naive", "B intermediate", "B memory",
                 "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive")

pbmc <- SortIdents(pbmc)

p1 <- CoveragePlot(
  object = pbmc,
  region = "MS4A1",
  features = "MS4A1",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = pbmc,
  region = "LYZ",
  features = "LYZ",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)


### Map cell type
predicted <- as.character(pbmc@meta.data[["predicted.id"]])

mapping_3class <- c(
  # B cells
  "B intermediate" = "Bcell", "B memory" = "Bcell",
  "B naive" = "Bcell", "Plasmablast" = "Bcell",
  
  # T / NK / innate lymphoid cells
  "CD4 CTL" = "T", "CD4 Naive" = "T", "CD4 Proliferating" = "T",
  "CD4 TCM" = "T", "CD4 TEM" = "T",
  "CD8 Naive" = "T", "CD8 TCM" = "T", "CD8 TEM" = "T",
  "dnT" = "T", "gdT" = "T", "MAIT" = "T", "Treg" = "T",
  "NK" = "T", "NK Proliferating" = "T", "NK_CD56bright" = "T", "ILC" = "T",
  
  # Myeloid / dendritic
  "CD14 Mono" = "Myeloid", "CD16 Mono" = "Myeloid",
  "ASDC" = "Myeloid", "cDC1" = "Myeloid",
  "cDC2" = "Myeloid", "pDC" = "Myeloid",
  
  # Excluded
  "Eryth" = NA, "HSPC" = NA
)

pbmc@meta.data$newCT <- unname(mapping_3class[predicted])
pbmc@meta.data$newCT[is.na(pbmc@meta.data$newCT)] <- "Others"
table(pbmc@meta.data$newCT, useNA = "ifany")
saveRDS(pbmc, "processed_pbmc.rds")


cell_counts <- table(pbmc$newCT)
legend_labels <- paste0(names(cell_counts), " (", as.integer(cell_counts), ")")

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap", group.by = "newCT") +  
  scale_color_discrete(labels = legend_labels) + ggtitle("PBMC")



# -
# Call peaks by cell type newCT
# -
print('Call peaks by cell type')

DefaultAssay(pbmc) <- "ATAC"
# call peaks using MACS2 by cell type
peaks <- CallPeaks(pbmc, group.by = 'newCT', macs2.path = "/home/chenxin/miniforge3/envs/macs3_env/bin/macs3")
warnings() %>% print

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
saveRDS(peaks, "data/MACS2_peaks_newCT.rds")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["ATAC_ct"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# get DNA sequence information
pbmc[['ATAC_ct']] <- RegionStats(
  object = pbmc[['ATAC_ct']],
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

print('peak calling by cell types ended')
saveRDS(pbmc, "data/processed_pbmc_peak_call.rds")


# Extract information
# RNA count matrix (genes x cells)
mrna <- GetAssayData(pbmc, assay = "RNA", layer = "counts")
# ATAC count matrix (peaks x cells)
atac <- GetAssayData(pbmc, assay = "ATAC", layer = "counts")
atac_pc <- GetAssayData(pbmc, assay = "ATAC_ct", layer = "counts")
# Pair
gene_peak_df <- GetGenePeakPairs(pbmc,
                                 expression.assay = "RNA",
                                 peak.assay = "ATAC_ct")
# Metadata (cell-level)
meta <- pbmc@meta.data
meta <- meta[, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                 "nCount_ATAC" , "nFeature_ATAC", 
                 "nFeature_ATAC_ct", "nCount_ATAC_ct", 
                 "percent.mt", "predicted.id", "newCT")]

saveRDS(mrna, "data/pbmc_multimodal.rna.rds")
saveRDS(atac, "data/pbmc_multimodal.atac.rds")
saveRDS(atac_pc, "data/pbmc_multimodal.atac_pc.rds")
saveRDS(meta, "data/pbmc_multimodal.meta.rds")
saveRDS(gene_peak_df, "data/pbmc_gene_peak_df.rds")




# -
# Call peaks by cell type predicted.id
# -
print('Call peaks by cell type')

DefaultAssay(pbmc) <- "ATAC"
# call peaks using MACS2 by cell type
peaks <- CallPeaks(pbmc, group.by = 'predicted.id', 
                   macs2.path = "/home/chenxin/miniforge3/envs/macs3_env/bin/macs3")
warnings() %>% print

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
saveRDS(peaks, "data/MACS2_peaks_predicted.id.rds")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["ATAC_ct2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# get DNA sequence information
pbmc[['ATAC_ct2']] <- RegionStats(
  object = pbmc[['ATAC_ct2']],
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

print('peak calling by cell types ended')
saveRDS(pbmc, "data/processed_pbmc_peak_call2.rds")


atac_pc2 <- GetAssayData(pbmc, assay = "ATAC_ct2", layer = "counts")
saveRDS(atac_pc2, "data/pbmc_multimodal.atac_pc2.rds")




## PCHiC
mapping_vec2 <- c(
  # 🧫 B cell lineage → Naïve / Total B (nB / tB)
  "B intermediate"   = "tB",
  "B memory"         = "tB",
  "B naive"          = "nB",
  "Plasmablast"      = "tB",
  
  # 🧬 Monocyte lineage → Monocytes (Mon)
  "CD14 Mono"        = "Mon",
  "CD16 Mono"        = "Mon",
  
  # 🩸 CD4 T cell lineage → Naïve / Total / Activated CD4 (nCD4 / tCD4 / aCD4)
  "CD4 Naive"        = "nCD4",
  "CD4 TCM"          = "tCD4",
  "CD4 TEM"          = "tCD4",
  "CD4 CTL"          = "aCD4",
  "CD4 Proliferating"= "aCD4",
  "Treg"             = "tCD4",
  
  # 🔥 CD8 T cell lineage → Naïve / Total CD8 (nCD8 / tCD8)
  "CD8 Naive"        = "nCD8",
  "CD8 TCM"          = "tCD8",
  "CD8 TEM"          = "tCD8",
  
  # ❤️ Erythroid lineage → Erythroblasts (Ery)
  "Eryth"            = "Ery",
  
  # ⚠️ Unmatched / not present in PCHi-C
  "ASDC"             = NA,
  "cDC1"             = NA,
  "cDC2"             = NA,
  "pDC"              = NA,
  "HSPC"             = NA,
  "dnT"              = NA,
  "gdT"              = NA,
  "MAIT"             = NA,
  "ILC"              = NA,
  "NK"               = NA,
  "NK Proliferating" = NA,
  "NK_CD56bright"    = NA
)



pbmc@meta.data$PCHiC_celltype <- mapping_vec2[pbmc@meta.data$predicted.id]
table(pbmc$PCHiC_celltype, useNA = "ifany")
table(pbmc$predicted.id, useNA = "ifany")
