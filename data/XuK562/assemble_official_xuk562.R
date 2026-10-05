#!/usr/bin/env Rscript

# Assemble the official Xu et al. 10x Multiome outputs after RNA/ATAC
# alignment.  This script keeps only paired barcodes, calculates the same
# per-cell QC quantities listed in Dorans et al. Supplementary Table 14, and
# writes the EpiLink-format matrices/metadata.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "data/XuK562/assemble_official_xuk562.R"
root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
data_dir <- file.path(root, "data", "XuK562")
rna_dir <- file.path(data_dir, "official_pipeline", "rna", "outs", "Solo.out", "Gene", "filtered")
fragment_file <- file.path(data_dir, "official_pipeline", "atac", "outs", "fragments.tsv.gz")
peak_file <- file.path(data_dir, "official_pipeline", "atac", "macs2_full", "aggregate_peaks.narrowPeak")

need <- c(
  file.path(rna_dir, c("matrix.mtx", "features.tsv", "barcodes.tsv")),
  fragment_file, peak_file
)
if (any(!file.exists(need))) {
  stop("Missing official pipeline output(s): ", paste(need[!file.exists(need)], collapse = ", "))
}

read_mex <- function(dir, feature_names = c("id", "name")) {
  mat <- readMM(file.path(dir, "matrix.mtx"))
  features <- read.delim(file.path(dir, "features.tsv"), header = FALSE,
                         stringsAsFactors = FALSE, check.names = FALSE)
  barcodes <- readLines(file.path(dir, "barcodes.tsv"))
  if (nrow(features) < nrow(mat)) stop("Feature file is shorter than matrix rows: ", dir)
  rownames(mat) <- make.unique(as.character(features[[min(2L, ncol(features))]]))
  colnames(mat) <- barcodes
  as(mat, "dgCMatrix")
}

rna <- read_mex(rna_dir)
common <- colnames(rna)
if (!length(common)) stop("RNA matrix has no cells")

annotation <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
annotation_levels <- GenomeInfoDb::seqlevels(annotation)
GenomeInfoDb::seqlevels(annotation) <- ifelse(annotation_levels == "MT", "chrM", paste0("chr", annotation_levels))

meta <- data.frame(
  barcodes = common,
  newCT = "K562",
  nCount_RNA = Matrix::colSums(rna),
  nCount_ATAC = NA_real_,
  row.names = common,
  stringsAsFactors = FALSE
)
obj <- CreateSeuratObject(counts = rna, assay = "RNA", meta.data = meta)
frags <- CreateFragmentObject(path = fragment_file, cells = common, validate.fragments = FALSE)
peaks <- rtracklayer::import(peak_file)
peaks <- peaks[as.character(GenomeInfoDb::seqnames(peaks)) != "chrM"]
blacklist_file <- file.path(root, "reference", "XuK562", "hg38-blacklist.bed.gz")
if (file.exists(blacklist_file)) {
  bl <- rtracklayer::import(blacklist_file)
  peaks <- peaks[!overlapsAny(peaks, bl)]
}
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- sort(unique(peaks))
atac <- FeatureMatrix(fragments = frags, features = peaks, cells = common,
                      process_n = 2000, verbose = TRUE)
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = atac, fragments = frags, annotation = annotation,
  sep = c("-", "-")
)
obj[["nCount_ATAC"]] <- Matrix::colSums(atac)
obj <- NucleosomeSignal(obj, assay = "ATAC")
obj <- TSSEnrichment(obj, assay = "ATAC", fast = TRUE)
meta <- obj[[]]
meta$barcodes <- rownames(meta)
meta$NucleosomeRatio <- as.numeric(meta$nucleosome_signal)
meta$TSSEnrichment <- as.numeric(meta$TSS.enrichment)
meta$qc_pass <- with(meta,
  is.finite(nCount_RNA) & nCount_RNA > 1000 &
    is.finite(nCount_ATAC) & nCount_ATAC > 200 &
    is.finite(NucleosomeRatio) & NucleosomeRatio < 2 &
    is.finite(TSSEnrichment) & TSSEnrichment > 2
)

keep <- rownames(meta)[meta$qc_pass]
rna <- rna[, keep, drop = FALSE]
atac <- atac[, keep, drop = FALSE]
meta <- meta[keep, c("barcodes", "newCT", "nCount_RNA", "nCount_ATAC",
                      "NucleosomeRatio", "TSSEnrichment", "qc_pass"), drop = FALSE]
stopifnot(identical(colnames(rna), colnames(atac)), identical(colnames(rna), rownames(meta)))

# STARsolo feature names are normally Ensembl IDs. Map them to symbols once,
# preserving already-symbol-named rows and removing duplicate symbols.
ensdb_annotation <- readRDS(file.path(root, "data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
feature_ids <- rownames(rna)
symbols <- as.character(ensdb_annotation$gene_name)[match(feature_ids, as.character(ensdb_annotation$gene_id))]
symbols[!grepl("^ENSG", feature_ids)] <- feature_ids[!grepl("^ENSG", feature_ids)]
gene_keep <- !is.na(symbols) & nzchar(symbols) & !duplicated(symbols)
rna <- rna[gene_keep, , drop = FALSE]
rownames(rna) <- symbols[gene_keep]

saveRDS(rna, file.path(data_dir, "XuK562.rna.rds"), compress = FALSE)
saveRDS(atac, file.path(data_dir, "XuK562.atac.rds"), compress = FALSE)
saveRDS(meta, file.path(data_dir, "XuK562.meta.rds"), compress = FALSE)
write.table(
  data.frame(metric = c("paired_cells_before_qc", "cells_passing_qc", "mean_rna_counts",
                        "mean_atac_counts"),
             value = c(length(common), nrow(meta), mean(meta$nCount_RNA), mean(meta$nCount_ATAC))),
  file.path(data_dir, "XuK562.qc_summary.tsv"), sep = "\t", quote = FALSE,
  row.names = FALSE
)
message("Saved official XuK562 matrices with ", nrow(meta), " QC-passing paired cells.")
