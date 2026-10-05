suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
  library(Signac)
  library(GenomicRanges)
})
source(file.path("R", "helper.R"))

out_dir <- file.path("data", "Brain2")
raw_dir <- file.path(out_dir, "raw")
h5_file <- file.path(raw_dir, "GSE214979_filtered_feature_bc_matrix.h5")
metadata_file <- file.path(raw_dir, "GSE214979_cell_metadata.csv.gz")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

barcodes <- as.character(h5read(h5_file, "matrix/barcodes"))
feature_type <- as.character(h5read(h5_file, "matrix/features/feature_type"))
feature_id <- as.character(h5read(h5_file, "matrix/features/id"))
feature_name <- as.character(h5read(h5_file, "matrix/features/name"))
gene_mask <- feature_type == "Gene Expression"
peak_mask <- feature_type == "Peaks"
n_gene <- sum(gene_mask)
n_peak <- sum(peak_mask)
gene_names <- make.unique(feature_name[gene_mask])
peak_names <- CanonicalPeakNames(feature_id[peak_mask])

# Read the 10x CSC matrix in column chunks instead of materializing the
# combined 187k-feature matrix before splitting RNA and ATAC.
read_modalities <- function(path, n_gene, n_peak, barcodes, chunk_cells = 2000L) {
  indptr <- as.numeric(h5read(path, "matrix/indptr"))
  n_cells <- length(barcodes)
  counts_rna <- integer(n_cells)
  for (chunk_i in seq_len(ceiling(n_cells / chunk_cells))) {
    first <- (chunk_i - 1L) * chunk_cells + 1L
    last <- min(chunk_i * chunk_cells, n_cells)
    offset <- indptr[first]
    nnz <- indptr[last + 1L] - offset
    if (nnz == 0) next
    idx <- as.integer(h5read(path, "matrix/indices", start = offset + 1, count = nnz)) + 1L
    column <- rep(seq_len(last - first + 1L), diff(indptr[first:(last + 1L)]))
    counts_rna[first:last] <- tabulate(column[idx <= n_gene], nbins = last - first + 1L)
  }
  counts_atac <- diff(indptr) - counts_rna
  p_rna <- as.integer(c(0, cumsum(counts_rna)))
  p_atac <- as.integer(c(0, cumsum(counts_atac)))
  i_rna <- integer(p_rna[n_cells + 1L]); x_rna <- numeric(length(i_rna))
  i_atac <- integer(p_atac[n_cells + 1L]); x_atac <- numeric(length(i_atac))

  # Fill the CSC slots directly. This avoids the large temporary copies made
  # by cbind() when concatenating dozens of sparse matrix chunks.
  for (chunk_i in seq_len(ceiling(n_cells / chunk_cells))) {
    first <- (chunk_i - 1L) * chunk_cells + 1L
    last <- min(chunk_i * chunk_cells, n_cells)
    offset <- indptr[first]
    nnz <- indptr[last + 1L] - offset
    if (nnz == 0) next
    idx <- as.integer(h5read(path, "matrix/indices", start = offset + 1, count = nnz)) + 1L
    values <- as.numeric(h5read(path, "matrix/data", start = offset + 1, count = nnz))
    local_ptr <- indptr[first:(last + 1L)] - offset
    for (jj in seq_len(last - first + 1L)) {
      lo <- local_ptr[jj] + 1L
      hi <- local_ptr[jj + 1L]
      if (hi < lo) next
      col <- first + jj - 1L
      ix <- idx[lo:hi]
      vv <- values[lo:hi]
      is_rna <- ix <= n_gene
      if (any(is_rna)) {
        dest <- (p_rna[col] + 1L):(p_rna[col + 1L])
        i_rna[dest] <- ix[is_rna] - 1L
        x_rna[dest] <- vv[is_rna]
      }
      if (any(!is_rna)) {
        dest <- (p_atac[col] + 1L):(p_atac[col + 1L])
        i_atac[dest] <- ix[!is_rna] - n_gene - 1L
        x_atac[dest] <- vv[!is_rna]
      }
    }
    if (chunk_i %% 5L == 0L || last == n_cells) message("matrix columns ", last, "/", n_cells)
  }
  rna <- new("dgCMatrix", i = i_rna, p = p_rna, x = x_rna,
             Dim = as.integer(c(n_gene, n_cells)),
             Dimnames = list(gene_names, barcodes))
  atac <- new("dgCMatrix", i = i_atac, p = p_atac, x = x_atac,
              Dim = as.integer(c(n_peak, n_cells)),
              Dimnames = list(peak_names, barcodes))
  dim(rna) <- c(n_gene, n_cells)
  dim(atac) <- c(n_peak, n_cells)
  list(rna = rna, atac = atac)
}

matrices <- read_modalities(h5_file, n_gene, n_peak, barcodes)
rna <- matrices$rna
atac <- matrices$atac
rm(matrices)
gc()

meta_raw <- read.csv(gzfile(metadata_file), check.names = FALSE,
                     stringsAsFactors = FALSE)
cell_ids <- as.character(meta_raw[[1L]])
rownames(meta_raw) <- cell_ids
meta_raw <- meta_raw[, -1L, drop = FALSE]
meta_raw <- meta_raw[barcodes, , drop = FALSE]
stopifnot(identical(rownames(meta_raw), colnames(rna)))

meta_raw$newCT <- unname(c(
  Excitatory = "Excitatory neurons",
  Inhibitory = "Inhibitory neurons",
  Oligodendrocytes = "Oligodendrocytes",
  Astrocytes = "Astrocytes",
  Microglia = "Microglia",
  OPCs = "OPCs",
  Endothelial = "Endothelial",
  Pericytes = "Pericytes"
)[as.character(meta_raw$predicted.id)])
meta_raw$celltype <- meta_raw$newCT
meta_raw$cell <- rownames(meta_raw)
meta_raw$donor <- as.character(meta_raw$id)
meta_raw$sample <- as.character(meta_raw$id)
meta_raw$condition <- as.character(meta_raw$Status)
meta_raw$disease <- as.character(meta_raw$Status)
meta_raw$batch <- as.character(meta_raw$`Sub.batch`)

# Keep the project-wide benchmark rule n_cells > 500. The source-level
# eight-celltype table remains available in Brain2.source_celltype_counts.tsv.
source_counts <- as.data.frame(table(meta_raw$newCT), stringsAsFactors = FALSE)
colnames(source_counts) <- c("celltype", "n_cells")
write.table(source_counts, file.path(out_dir, "Brain2.source_celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
keep_types <- source_counts$celltype[source_counts$n_cells > 500]
keep <- meta_raw$newCT %in% keep_types
rna <- rna[, keep, drop = FALSE]
atac <- atac[, keep, drop = FALSE]
meta <- meta_raw[keep, , drop = FALSE]

# The GEO metadata already contains the QC values used by the source study;
# reusing them avoids materializing rna > 0 and atac > 0 for this large matrix.
meta$log_nCount_RNA <- log(pmax(meta$nCount_RNA, 1))
meta$log_nCount_ATAC <- log(pmax(meta$nCount_ATAC, 1))
rm(meta_raw)
gc()
stopifnot(identical(colnames(rna), colnames(atac)))
stopifnot(identical(rownames(meta), colnames(rna)))
stopifnot(all(table(meta$newCT) > 500))

saveRDS(rna, file.path(out_dir, "Brain2.rna.rds"), compress = FALSE)
saveRDS(atac, file.path(out_dir, "Brain2.atac.rds"), compress = FALSE)
saveRDS(meta, file.path(out_dir, "Brain2.meta.rds"), compress = "gzip")

counts <- as.data.frame(table(meta$newCT), stringsAsFactors = FALSE)
colnames(counts) <- c("celltype", "n_cells")
write.table(counts, file.path(out_dir, "Brain2.celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

qc_metrics <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC",
                "nucleosome_signal", "TSS.enrichment", "blacklist_fraction")
qc_summary <- do.call(rbind, lapply(qc_metrics, function(nm) {
  x <- as.numeric(meta[[nm]])
  data.frame(metric = nm, n = sum(is.finite(x)), min = min(x, na.rm = TRUE),
             median = median(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE),
             max = max(x, na.rm = TRUE))
}))
write.table(qc_summary, file.path(out_dir, "Brain2.qc_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sample_counts <- as.data.frame(table(meta$donor, meta$condition, meta$newCT),
                               stringsAsFactors = FALSE)
colnames(sample_counts) <- c("donor", "condition", "celltype", "n_cells")
sample_counts <- sample_counts[sample_counts$n_cells > 0, , drop = FALSE]
write.table(sample_counts, file.path(out_dir, "Brain2.sample_celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Candidate pairs: hg38 gene body +/-500 kb against the filtered peak set.
annotation <- readRDS(file.path("data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
annotation_levels <- GenomeInfoDb::seqlevels(annotation)
GenomeInfoDb::seqlevels(annotation) <- ifelse(
  annotation_levels == "MT", "chrM",
  ifelse(grepl("^chr", annotation_levels), annotation_levels,
         paste0("chr", annotation_levels))
)
gene_range <- CollapseToLongestTranscript(annotation)
gene_window <- Signac::Extend(gene_range, upstream = 5e5, downstream = 5e5)
peaks <- Signac::StringToGRanges(rownames(atac), sep = c("-", "-"))
common_chr <- intersect(GenomeInfoDb::seqlevels(peaks), GenomeInfoDb::seqlevels(gene_window))
peaks <- GenomeInfoDb::keepSeqlevels(peaks, common_chr, pruning.mode = "coarse")
gene_window <- GenomeInfoDb::keepSeqlevels(gene_window, common_chr, pruning.mode = "coarse")
hits <- GenomicRanges::findOverlaps(peaks, gene_window)
gene_peak_df <- data.frame(
  gene = gene_window$gene_name[S4Vectors::subjectHits(hits)],
  peak = rownames(atac)[S4Vectors::queryHits(hits)],
  stringsAsFactors = FALSE
)
gene_peak_df <- gene_peak_df[gene_peak_df$gene %in% rownames(rna), , drop = FALSE]
gene_peak_df <- unique(gene_peak_df[, c("gene", "peak"), drop = FALSE])
rownames(gene_peak_df) <- NULL
saveRDS(gene_peak_df, file.path(out_dir, "Brain2.gene_peak_df.rds"), compress = "gzip")

qc_thresholds <- data.frame(
  dataset = "GSE214979", source_cells = 105332L, source_celltypes = 8L,
  retained_cells = ncol(rna), retained_celltypes = length(keep_types),
  rna_genes = nrow(rna), atac_peaks = nrow(atac),
  source_qc = "GEO filtered feature-barcode matrix + supplied metadata",
  additional_filter = "celltype n_cells > 500",
  stringsAsFactors = FALSE
)
write.table(qc_thresholds, file.path(out_dir, "Brain2.qc_thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "Brain2 benchmark preparation", "",
  "Source: GEO GSE214979, Anderson et al. (2023), postmortem DLPFC snRNA+snATAC multiome.",
  "The filtered 10x H5 contains 105,332 nuclei, 36,601 genes and 150,614 peaks.",
  "The supplied cell metadata contains RNA/ATAC depth, TSS and nucleosome QC metrics.",
  "The source has 8 annotated cell types; this benchmark export retains cell types with >500 cells, leaving 6 types.",
  "The source-level 8-type counts are in Brain2.source_celltype_counts.tsv.",
  "Candidate gene_peak_df uses hg38 gene bodies +/-500 kb and the filtered peak set.",
  "The scMultiMap paper reports 380,154 peaks in its summary table; this differs from the 150,614 peak features in the GEO filtered H5.",
  "The 59 GB GEO fragment file was not downloaded; the supplied metadata is used for QC and the filtered count matrix for benchmarking."
), file.path(out_dir, "README.md"))

cat("Brain2 complete:", ncol(rna), "cells;", length(keep_types),
    "celltypes;", nrow(rna), "genes;", nrow(atac), "peaks;",
    nrow(gene_peak_df), "gene-peak pairs.\n")
