suppressPackageStartupMessages({library(Matrix); library(Seurat); library(Signac)})
source(file.path("R", "helper.R"))
rna <- readRDS("data/Arthritis/Arthritis.rna.rds")
atac <- readRDS("data/Arthritis/Arthritis.atac.rds")
meta <- readRDS("data/Arthritis/Arthritis.meta.rds")
pairs <- readRDS("data/Arthritis/Arthritis.gene_peak_df.rds")
keep <- rownames(meta)[meta$newCT == "Fibroblast"]
atac <- atac[, keep, drop = FALSE]
object <- Seurat::CreateSeuratObject(counts = rna[, keep, drop = FALSE], meta.data = meta[keep, , drop = FALSE])
object[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac, sep = InferPeakSeparator(rownames(atac))
)
pairs$peak <- CanonicalPeakNames(pairs$peak)
atac_counts <- Seurat::GetAssayData(object, assay = "ATAC", layer = "counts")
peak_names <- rownames(atac_counts)
overlap <- sum(pairs$peak %in% peak_names)
detect <- Matrix::rowMeans(atac_counts > 0)
rate <- detect[pairs$peak]
keep_pairs <- pairs$gene %in% rownames(rna) & is.finite(rate) & rate >= 0.05
cat("assay peaks:", length(peak_names), "\n")
cat("candidate pairs:", nrow(pairs), "\n")
cat("pair peak overlap:", overlap, "\n")
cat("pairs with RNA gene and ATAC detection >= 5%:", sum(keep_pairs), "\n")
stopifnot(overlap > 0L, sum(keep_pairs) > 0L)
