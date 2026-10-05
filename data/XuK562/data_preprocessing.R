# Build candidate gene-peak pairs for Xu K562.

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "data/XuK562/data_preprocessing.R"
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(project_root, "R", "helper.R"))

data_dir <- file.path(project_root, "data", "XuK562")
annotation <- readRDS(file.path(project_root, "data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
rna <- readRDS(file.path(data_dir, "XuK562.rna.rds"))
atac <- readRDS(file.path(data_dir, "XuK562.atac.rds"))
meta <- readRDS(file.path(data_dir, "XuK562.meta.rds"))

stopifnot(identical(colnames(rna), colnames(atac)))
meta <- meta[colnames(rna), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(rna)))

object <- CreateSeuratObject(counts = rna, meta.data = meta, assay = "RNA")
object[["ATAC"]] <- CreateChromatinAssay(counts = atac)
gene_peak_df <- GetGenePeakPairs(object, annotation)
saveRDS(gene_peak_df, file.path(data_dir, "XuK562.gene_peak_df.rds"), compress = TRUE)

message("Saved ", nrow(gene_peak_df), " candidate gene-peak pairs.")
