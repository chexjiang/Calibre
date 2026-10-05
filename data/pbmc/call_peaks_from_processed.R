#!/usr/bin/env Rscript

# Call cell-type-specific ATAC peaks from the QC-filtered PBMC object.
# The input object already contains the 10x peak matrix and the linked
# fragment file; MACS2 is used here to create a new B/T/Myeloid peak set.

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

data_dir <- "/home/chenxin/proj/EpiLink/data/pbmc"
setwd(data_dir)

input_file <- "processed_pbmc.rds"
output_object <- "processed_pbmc_peak_call.rds"
output_peaks <- "MACS2_peaks_BTM.rds"
output_rna <- "pbmc_multimodal.rna.rds"
output_atac <- "pbmc_multimodal.atac.rds"
output_meta <- "pbmc_multimodal.meta.rds"
output_gene_peak <- "pbmc_gene_peak_df.rds"

macs2_path <- "/home/chenxin/miniforge3/envs/macs3_env/bin/macs3"
if (!file.exists(macs2_path)) {
  stop("MACS3 executable not found: ", macs2_path, call. = FALSE)
}

pbmc <- readRDS(input_file)
if (!"newCT" %in% colnames(pbmc[[]])) {
  stop("newCT is not present in ", input_file, call. = FALSE)
}

target_celltypes <- c("Bcell", "Myeloid", "T")
pbmc <- subset(pbmc, subset = newCT %in% target_celltypes)
pbmc$newCT <- factor(as.character(pbmc$newCT), levels = target_celltypes)
print(table(pbmc$newCT, useNA = "ifany"))

DefaultAssay(pbmc) <- "ATAC"
fragments <- Fragments(pbmc)
if (!length(fragments)) {
  stop("No fragment file is linked to the ATAC assay", call. = FALSE)
}

peaks <- CallPeaks(
  object = pbmc,
  group.by = "newCT",
  macs2.path = macs2_path
)

# Keep standard hg38 chromosomes. Apply the project blacklist here as well if
# blacklist_hg38_unified is loaded in the calling environment.
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
if (exists("blacklist_hg38_unified", inherits = TRUE)) {
  peaks <- subsetByOverlaps(
    x = peaks, ranges = get("blacklist_hg38_unified"), invert = TRUE
  )
}
saveRDS(peaks, output_peaks)

peak_counts <- FeatureMatrix(
  fragments = fragments,
  features = peaks,
  cells = colnames(pbmc)
)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

pbmc[["ATAC_ct"]] <- CreateChromatinAssay(
  counts = peak_counts,
  fragments = fragments,
  annotation = annotation
)
pbmc[["ATAC_ct"]] <- RegionStats(
  object = pbmc[["ATAC_ct"]],
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)

saveRDS(pbmc, output_object)
saveRDS(GetAssayData(pbmc, assay = "RNA", layer = "counts"), output_rna)
saveRDS(GetAssayData(pbmc, assay = "ATAC_ct", layer = "counts"), output_atac)
saveRDS(pbmc@meta.data, output_meta)
saveRDS(
  GetGenePeakPairs(
    pbmc, expression.assay = "RNA", peak.assay = "ATAC_ct"
  ),
  output_gene_peak
)

message("Finished B/T/Myeloid peak calling")
