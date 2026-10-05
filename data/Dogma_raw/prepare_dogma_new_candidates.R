#!/usr/bin/env Rscript

# Generate the common +/-500 kb gene-peak candidate universe for the two new
# DOGMA objects. The matrices are already QC filtered by process_dogma_new.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
})

file_argument <- sub("^--file=", "", commandArgs(FALSE)[
  grep("^--file=", commandArgs(FALSE))[[1L]]
])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)
source("R/helper.R")
annotation <- readRDS(file.path(project_root, "data",
                                "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
names(annotation) <- NULL

for (condition in c("CTRL", "STIM")) {
  directory <- file.path(project_root, "data", paste0("Dogma_", condition, "_new"))
  prefix <- paste0("Dogma_", condition, "_new")
  rna <- readRDS(file.path(directory, paste0(prefix, ".rna.rds")))
  atac <- readRDS(file.path(directory, paste0(prefix, ".atac.rds")))
  meta <- readRDS(file.path(directory, paste0(prefix, ".meta.rds")))
  stopifnot(identical(colnames(rna), colnames(atac)))
  meta <- meta[colnames(rna), , drop = FALSE]
  object <- CreateSeuratObject(counts = rna, meta.data = meta, assay = "RNA")
  object[["ATAC"]] <- CreateChromatinAssay(
    counts = atac, annotation = annotation, sep = c("-", "-")
  )
  pairs <- GetGenePeakPairs(object, annotation = annotation, distance = 5e5)
  saveRDS(pairs, file.path(directory, paste0(prefix, ".gene_peak_df.rds")))
  message("Wrote ", prefix, ".gene_peak_df.rds: ", nrow(pairs), " pairs")
}
