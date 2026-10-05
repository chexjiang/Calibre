suppressPackageStartupMessages({library(rhdf5); library(Signac); library(GenomicRanges)})
source(file.path("R", "helper.R"))
out_dir <- file.path("data", "Brain2")
h5_file <- file.path(out_dir, "raw", "GSE214979_filtered_feature_bc_matrix.h5")
feature_type <- as.character(h5read(h5_file, "matrix/features/feature_type"))
feature_id <- as.character(h5read(h5_file, "matrix/features/id"))
feature_name <- as.character(h5read(h5_file, "matrix/features/name"))
gene_names <- make.unique(feature_name[feature_type == "Gene Expression"])
peak_names <- CanonicalPeakNames(feature_id[feature_type == "Peaks"])
annotation <- readRDS(file.path("data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
annotation_levels <- GenomeInfoDb::seqlevels(annotation)
GenomeInfoDb::seqlevels(annotation) <- ifelse(
  annotation_levels == "MT", "chrM",
  ifelse(grepl("^chr", annotation_levels), annotation_levels,
         paste0("chr", annotation_levels))
)
gene_range <- CollapseToLongestTranscript(annotation)
gene_window <- Signac::Extend(gene_range, upstream = 5e5, downstream = 5e5)
peaks <- Signac::StringToGRanges(peak_names, sep = c("-", "-"))
common_chr <- intersect(GenomeInfoDb::seqlevels(peaks), GenomeInfoDb::seqlevels(gene_window))
peaks <- GenomeInfoDb::keepSeqlevels(peaks, common_chr, pruning.mode = "coarse")
gene_window <- GenomeInfoDb::keepSeqlevels(gene_window, common_chr, pruning.mode = "coarse")
hits <- GenomicRanges::findOverlaps(peaks, gene_window)
gene_peak_df <- data.frame(
  gene = gene_window$gene_name[S4Vectors::subjectHits(hits)],
  peak = peak_names[S4Vectors::queryHits(hits)], stringsAsFactors = FALSE
)
gene_peak_df <- gene_peak_df[gene_peak_df$gene %in% gene_names, , drop = FALSE]
gene_peak_df <- unique(gene_peak_df[, c("gene", "peak"), drop = FALSE])
rownames(gene_peak_df) <- NULL
saveRDS(gene_peak_df, file.path(out_dir, "Brain2.gene_peak_df.rds"), compress = "gzip")
cat("Brain2 gene-peak pairs:", nrow(gene_peak_df), "\n")
