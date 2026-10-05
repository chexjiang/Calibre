#!/usr/bin/env Rscript

# Build the Satpathy K562 benchmark input from the released IGVF analysis set
# IGVFDS2018URDY. The analysis set is a MULTI-seq multiplex of Jurkat, K562,
# THP-1, and TeloHAEC. We use the RNA matrix and MULTI-seq calls to identify
# the K562 barcode groups, then build the ATAC peak matrix from the official
# GRCh38 fragments. Only cells with both modalities are retained.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(EnsDb.Hsapiens.v86)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "data/SatpathyK562/data_preprocessing.R"
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
data_dir <- file.path(project_root, "data", "SatpathyK562")
raw_dir <- file.path(data_dir, "raw")

rna_h5 <- file.path(raw_dir, "IGVFFI2641LMFU.h5")
hash_file <- file.path(raw_dir, "IGVFFI5156ENLL.tsv.gz")
fragment_file <- file.path(raw_dir, "IGVFFI0627MDWY.tsv.gz")
fragment_index <- file.path(raw_dir, "IGVFFI3930OIYP.tbi")
peak_file <- file.path(raw_dir, "IGVFFI8318PQLH.bed.gz")
needed <- c(rna_h5, hash_file, fragment_file, fragment_index, peak_file)
if (any(!file.exists(needed))) {
  stop("Missing official Satpathy input(s): ", paste(needed[!file.exists(needed)], collapse = ", "))
}
fragment_index_expected <- paste0(fragment_file, ".tbi")
if (!file.exists(fragment_index_expected)) {
  file.copy(fragment_index, fragment_index_expected, overwrite = FALSE)
}

message("Reading official RNA matrix")
selected_rna_dir <- file.path(raw_dir, "selected_rna")
if (file.exists(file.path(selected_rna_dir, "matrix.mtx"))) {
  message("Using cached RNA MEX for MULTI-seq barcodes")
  rna_all <- readMM(file.path(selected_rna_dir, "matrix.mtx"))
  feature_table <- read.delim(file.path(selected_rna_dir, "features.tsv"), header = FALSE,
                              stringsAsFactors = FALSE, check.names = FALSE)
  rna_barcodes <- readLines(file.path(selected_rna_dir, "barcodes.tsv"))
  rownames(rna_all) <- make.unique(as.character(feature_table[[min(2L, ncol(feature_table))]]))
  colnames(rna_all) <- rna_barcodes
  rna_all <- as(rna_all, "dgCMatrix")
} else {
  rna_all <- Read10X_h5(rna_h5, use.names = TRUE, unique.features = TRUE)
  if (is.list(rna_all)) {
    if (!"Gene Expression" %in% names(rna_all)) stop("RNA H5 has no Gene Expression matrix")
    rna_all <- rna_all[["Gene Expression"]]
  }
  rna_all <- as(rna_all, "dgCMatrix")
}

message("Reading MULTI-seq barcode calls")
hash <- fread(hash_file, data.table = FALSE, check.names = FALSE)
if (!"cellBC" %in% names(hash)) stop("Unexpected MULTI-seq table: no cellBC column")
hash$cellBC <- paste0(hash$cellBC, "-1")
tag_columns <- setdiff(names(hash), "cellBC")
hash_counts <- as.matrix(hash[, tag_columns, drop = FALSE])
storage.mode(hash_counts) <- "numeric"
hash_tag_index <- max.col(hash_counts, ties.method = "first")
hash_total <- rowSums(hash_counts)
hash_max <- apply(hash_counts, 1L, max)
hash_tag <- tag_columns[hash_tag_index]
names(hash_tag) <- hash$cellBC
hash_total <- setNames(hash_total, hash$cellBC)
hash_max <- setNames(hash_max, hash$cellBC)

common <- intersect(colnames(rna_all), hash$cellBC)
if (!length(common)) stop("No RNA barcodes overlap the MULTI-seq table")
rna_hash <- rna_all[, common, drop = FALSE]
hash_tag <- hash_tag[common]
hash_total <- hash_total[common]
hash_max <- hash_max[common]

# Use the erythroid/K562 expression program to identify the two K562 MULTI-seq
# tag groups. This is deliberately data-driven: the released file contains
# several multiplexed cell lines, while the IGVF metadata does not provide a
# cell-type label for each barcode in the count matrix. The two top-scoring
# tags are expected to be the K562 technical/sample groups. The selected tags
# and their scores are written to the QC summary for auditability.
marker_sets <- list(
  K562 = c("HBA1", "HBA2", "HBB", "GATA1", "KLF1", "TAL1", "NFE2", "EPOR", "SLC25A37", "AHSP", "ALAS2", "SLC4A1"),
  Jurkat = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "LTB"),
  THP1 = c("LST1", "S100A8", "S100A9", "CTSS", "FCN1", "LYZ"),
  TeloHAEC = c("VWF", "PECAM1", "KDR", "ENG", "EMCN", "RAMP2")
)
score_markers <- function(m, genes) {
  genes <- intersect(genes, rownames(m))
  if (!length(genes)) return(setNames(rep(0, ncol(m)), colnames(m)))
  Matrix::colSums(m[genes, , drop = FALSE])
}
ncount_rna_hash <- Matrix::colSums(rna_hash)
marker_scores <- lapply(marker_sets, score_markers, m = rna_hash)
marker_scores <- lapply(marker_scores, function(x) 1e4 * x / pmax(ncount_rna_hash, 1))
tag_score <- tapply(marker_scores$K562, hash_tag, mean, na.rm = TRUE)
tag_n <- table(hash_tag)
tag_score <- sort(tag_score[names(tag_n) %in% names(tag_score)], decreasing = TRUE)
if (length(tag_score) < 2L) stop("Fewer than two usable MULTI-seq tag groups")
k562_tags <- names(tag_score)[seq_len(min(2L, length(tag_score)))]
message("Selected K562 MULTI-seq tags: ", paste(k562_tags, collapse = ", "))

selected <- common[hash_tag %in% k562_tags]
if (!length(selected)) stop("No K562 cells selected from MULTI-seq tags")
rna <- rna_all[, selected, drop = FALSE]
rm(rna_all, rna_hash, hash, hash_counts)
gc()

# 10x Multiome uses distinct GEX and ATAC barcode sequences. The two
# 737K-arc-v1 whitelist files are line-aligned by gel bead; translate the
# selected GEX barcodes to their corresponding ATAC barcodes before reading
# fragments. Without this translation the modalities have zero exact overlap.
gex_whitelist_file <- file.path(data_dir, "reference", "gex_737K-arc-v1.txt.gz")
atac_whitelist_file <- file.path(data_dir, "reference", "atac_737K-arc-v1.txt.gz")
if (!all(file.exists(c(gex_whitelist_file, atac_whitelist_file)))) {
  stop("Missing 10x Multiome GEX/ATAC barcode whitelists in ",
       file.path(data_dir, "reference"))
}
gex_whitelist <- readLines(gzfile(gex_whitelist_file))
atac_whitelist <- readLines(gzfile(atac_whitelist_file))
gex_selected <- sub("-1$", "", selected)
barcode_index <- match(gex_selected, sub("-1$", "", gex_whitelist))
if (anyNA(barcode_index)) {
  stop("Some selected GEX barcodes are absent from the 737K-arc-v1 whitelist")
}
selected_atac <- paste0(atac_whitelist[barcode_index], "-1")
names(selected_atac) <- selected

# Keep gene symbols in the same format as the other EpiLink datasets.
annotation <- readRDS(file.path(project_root, "data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
feature_ids <- rownames(rna)
symbols <- as.character(annotation$gene_name)[match(feature_ids, as.character(annotation$gene_id))]
symbols[!grepl("^ENSG", feature_ids)] <- feature_ids[!grepl("^ENSG", feature_ids)]
gene_keep <- !is.na(symbols) & nzchar(symbols) & !duplicated(symbols)
rna <- rna[gene_keep, , drop = FALSE]
rownames(rna) <- symbols[gene_keep]

# Use the published GRCh38 peaks, excluding mitochondrial and blacklist peaks.
peaks <- rtracklayer::import(peak_file)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- peaks[!as.character(seqnames(peaks)) %in% c("chrM", "MT", "M")]
blacklist_file <- file.path(project_root, "reference", "XuK562", "hg38-blacklist.bed.gz")
if (file.exists(blacklist_file)) {
  blacklist <- rtracklayer::import(blacklist_file)
  peaks <- peaks[!overlapsAny(peaks, blacklist)]
}
peaks <- sort(unique(peaks))
peak_limit <- as.integer(Sys.getenv("SATPATHY_PEAK_LIMIT", "0"))
if (is.finite(peak_limit) && peak_limit > 0L && peak_limit < length(peaks)) {
  # Optional practical mode for shared nodes: retain an evenly spaced subset
  # of the released peaks. The default (0) uses the complete peak set.
  peaks <- peaks[unique(round(seq(1, length(peaks), length.out = peak_limit)))]
  message("SATPATHY_PEAK_LIMIT selected ", length(peaks), " peaks")
}
peak_names <- paste0(as.character(seqnames(peaks)), "-", start(peaks), "-", end(peaks))
peak_names <- make.unique(peak_names)
names(peaks) <- peak_names

meta <- data.frame(
  barcodes = selected,
  atac_barcodes = unname(selected_atac),
  newCT = "K562",
  hash_tag = unname(hash_tag[selected]),
  hash_count = unname(hash_total[selected]),
  hash_max = unname(hash_max[selected]),
  nCount_RNA = Matrix::colSums(rna),
  row.names = selected,
  stringsAsFactors = FALSE
)
obj <- CreateSeuratObject(counts = rna, assay = "RNA", meta.data = meta)
frags <- CreateFragmentObject(path = fragment_file, cells = selected_atac,
                               validate.fragments = FALSE)
message("Building ATAC matrix for ", length(selected), " K562 cells and ", length(peaks), " peaks")
atac_cache <- file.path(data_dir, "SatpathyK562.atac.pre_qc.rds")
matrix_helper <- file.path(data_dir, "atac_matrix_cpp")
matrix_dir <- file.path(data_dir, "matrix_cache")
dir.create(matrix_dir, showWarnings = FALSE)
selected_file <- file.path(matrix_dir, "selected_barcodes.txt")
peak_bed <- file.path(matrix_dir, "peaks.bed")
matrix_file <- file.path(matrix_dir, "atac.mtx")
if (file.exists(atac_cache)) {
  message("Reading cached pre-QC ATAC matrix")
  atac <- readRDS(atac_cache)
} else if (file.exists(matrix_file)) {
  message("Reading cached Matrix Market ATAC matrix")
  atac <- Matrix::readMM(matrix_file)
  atac <- as(atac, "dgCMatrix")
  saveRDS(atac, atac_cache, compress = FALSE)
} else {
  # FeatureMatrix can allocate a very large intermediate object for the
  # released 189k-peak file.  The bundled streaming helper makes the same
  # overlap/count matrix without materialising all fragment-feature pairs.
  if (!file.exists(matrix_helper)) {
    source_helper <- file.path(data_dir, "atac_matrix.cpp")
    status <- system2("g++", c("-O3", "-std=c++17", source_helper,
                               "-lz", "-o", matrix_helper))
    if (!identical(status, 0L)) stop("Failed to compile atac_matrix.cpp")
  }
  writeLines(unname(selected_atac), selected_file)
  write.table(data.frame(chrom = as.character(seqnames(peaks)),
                         start = start(peaks) - 1L,
                         end = end(peaks)),
              peak_bed, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE)
  status <- system2(matrix_helper, c(fragment_file, selected_file,
                                     peak_bed, matrix_file))
  if (!identical(status, 0L)) stop("atac_matrix_cpp failed with status ", status)
  atac <- Matrix::readMM(matrix_file)
  atac <- as(atac, "dgCMatrix")
  saveRDS(atac, atac_cache, compress = FALSE)
}
colnames(atac) <- selected
rownames(atac) <- peak_names
chromatin_annotation <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
annotation_levels <- seqlevels(chromatin_annotation)
seqlevels(chromatin_annotation) <- ifelse(annotation_levels == "MT", "chrM", paste0("chr", annotation_levels))
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = atac, fragments = frags,
  annotation = chromatin_annotation,
  sep = c("-", "-")
)

# The official fragment file is several hundred million rows.  Signac's
# CountFragments(), NucleosomeSignal(), and TSSEnrichment() each scan it, which
# is unnecessarily expensive and can exhaust memory on shared nodes.  Compute
# the same QC inputs in one streaming pass with the bundled helper instead.
qc_dir <- file.path(data_dir, "qc_cache")
dir.create(qc_dir, showWarnings = FALSE)
selected_file <- file.path(qc_dir, "selected_barcodes.txt")
tss_file <- file.path(qc_dir, "tss.tsv")
stats_file <- file.path(qc_dir, "fragment_stats.tsv")
writeLines(unname(selected_atac), selected_file)
tss_pos <- ifelse(as.character(strand(chromatin_annotation)) == "+",
                  start(chromatin_annotation), end(chromatin_annotation))
tss_tbl <- unique(data.frame(chrom = as.character(seqnames(chromatin_annotation)),
                             pos = tss_pos, stringsAsFactors = FALSE))
write.table(tss_tbl, tss_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
skip_fragment_qc <- identical(Sys.getenv("SATPATHY_SKIP_FRAGMENT_QC", "0"), "1")
if (!skip_fragment_qc && (!file.exists(stats_file) || file.info(stats_file)$size == 0)) {
  helper_cpp <- file.path(data_dir, "fragment_qc_cpp")
  if (file.exists(helper_cpp)) {
    status <- system2(helper_cpp, c(fragment_file, selected_file,
                                    tss_file, stats_file))
  } else {
    helper <- file.path(data_dir, "fragment_qc.py")
    status <- system2("python3", c(helper, fragment_file, selected_file,
                                    tss_file, stats_file))
  }
  if (!identical(status, 0L)) stop("fragment_qc.py failed with status ", status)
}
if (skip_fragment_qc) {
  warning("SATPATHY_SKIP_FRAGMENT_QC=1: retaining cells using RNA/peak counts only; " ,
          "nucleosome/TSS metrics are intentionally left NA")
  fragment_stats <- data.frame(barcode = selected_atac,
    fragment_count = Matrix::colSums(atac), nucleosome_fragments = NA_real_,
    nfr_fragments = NA_real_, cut_core = NA_real_, cut_window = NA_real_)
} else {
  fragment_stats <- read.delim(stats_file, stringsAsFactors = FALSE,
                               check.names = FALSE)
  fragment_stats <- fragment_stats[match(selected_atac, fragment_stats$barcode), , drop = FALSE]
}
obj[["nFragments"]] <- fragment_stats$fragment_count
peak_counts <- as.numeric(Matrix::colSums(atac))
names(peak_counts) <- selected
obj[["nCount_ATAC"]] <- peak_counts
obj[["nucleosome_signal"]] <- with(fragment_stats,
  nucleosome_fragments / pmax(nfr_fragments, 1))
# This is a cut-site enrichment ratio: +/-100 bp around a TSS divided by the
# surrounding +/-1000 bp window, with the window-size correction used by the
# standard TSS-enrichment definition.
obj[["TSS.enrichment"]] <- with(fragment_stats,
  (cut_core / 200) / pmax((cut_window - cut_core) / 1800, 1e-12))

meta <- obj[[]]
meta$barcodes <- rownames(meta)
meta$NucleosomeRatio <- as.numeric(meta$nucleosome_signal)
meta$TSSEnrichment <- as.numeric(meta$TSS.enrichment)
meta$qc_pass <- with(meta,
  is.finite(nCount_RNA) & nCount_RNA > 1000 &
    is.finite(nCount_ATAC) & nCount_ATAC > 200 &
    if (skip_fragment_qc) TRUE else
      is.finite(NucleosomeRatio) & NucleosomeRatio < 2 &
      is.finite(TSSEnrichment) & TSSEnrichment > 2
)

keep <- rownames(meta)[meta$qc_pass]
# FeatureMatrix can return cells in a different order from the RNA H5. Align
# by barcode explicitly before writing the paired matrices.
common_cells <- Reduce(intersect, list(colnames(rna), colnames(atac), rownames(meta)))
keep <- intersect(keep, common_cells)
rna <- rna[, keep, drop = FALSE]
atac <- atac[, keep, drop = FALSE]
meta <- meta[keep, c("barcodes", "atac_barcodes", "newCT", "hash_tag", "hash_count", "hash_max",
                     "nCount_RNA", "nCount_ATAC", "nFragments", "NucleosomeRatio",
                     "TSSEnrichment", "qc_pass"), drop = FALSE]
stopifnot(identical(colnames(rna), colnames(atac)), identical(colnames(rna), rownames(meta)))

saveRDS(rna, file.path(data_dir, "SatpathyK562.rna.rds"), compress = FALSE)
saveRDS(atac, file.path(data_dir, "SatpathyK562.atac.rds"), compress = FALSE)
saveRDS(meta, file.path(data_dir, "SatpathyK562.meta.rds"), compress = FALSE)

summary_rows <- rbind(
  data.frame(metric = "paired_k562_cells_before_qc", value = length(selected)),
  data.frame(metric = "cells_passing_qc", value = nrow(meta)),
  data.frame(metric = "mean_rna_counts", value = mean(meta$nCount_RNA)),
  data.frame(metric = "mean_atac_peak_counts", value = mean(meta$nCount_ATAC)),
  data.frame(metric = "mean_atac_fragments", value = mean(meta$nFragments)),
  data.frame(metric = "k562_tags", value = paste(k562_tags, collapse = ",")),
  data.frame(metric = "qc_rna_threshold", value = ">1000"),
  data.frame(metric = "qc_atac_threshold", value = ">200"),
  data.frame(metric = "qc_nucleosome_threshold", value = "<2"),
  data.frame(metric = "qc_tss_threshold", value = ">2")
)
write.table(summary_rows, file.path(data_dir, "SatpathyK562.qc_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
tag_summary <- data.frame(tag = names(tag_n), cells = as.integer(tag_n[names(tag_n)]),
                          k562_score = unname(tag_score[names(tag_n)]),
                          selected_k562 = names(tag_n) %in% k562_tags)
write.table(tag_summary, file.path(data_dir, "SatpathyK562.multiseq_tag_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved SatpathyK562 matrices with ", nrow(meta), " QC-passing paired cells")
