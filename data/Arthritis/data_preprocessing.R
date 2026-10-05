# Prepare the arthritis-tissue benchmark inputs used by Sakaue et al. (2024).
#
# GEO provides both complete raw-count H5AD matrices and a processed ATAC H5AD.
# We use the complete raw gene/peak matrices here.  The source matrices are
# already post-QC; this script documents the paper's QC thresholds but does
# not apply a second, incompatible QC pass.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(rhdf5)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(ggplot2)
})

source(file.path("R", "helper.R"))

out_dir <- file.path("data", "Arthritis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_counts_h5ad <- function(path) {
  # AnnData stores X as CSR (cells x features).  Reading the three CSR arrays
  # directly avoids a dense DelayedArray coercion and is substantially faster
  # for the 295-million-entry arthritis ATAC matrix.
  feature_names <- as.character(rhdf5::h5read(path, "var/_index"))
  cell_names <- as.character(rhdf5::h5read(path, "obs/_index"))
  sce_meta <- zellkonverter::readH5AD(path, use_hdf5 = TRUE,
                                      reader = "R", verbose = FALSE)
  metadata <- as.data.frame(SummarizedExperiment::colData(sce_meta),
                            stringsAsFactors = FALSE)
  rm(sce_meta)
  shape <- as.integer(c(length(cell_names), length(feature_names)))
  data <- rhdf5::h5read(path, "X/data")
  indices <- rhdf5::h5read(path, "X/indices")
  indptr <- rhdf5::h5read(path, "X/indptr")
  csr <- Matrix::sparseMatrix(
    p = as.integer(indptr), j = as.integer(indices),
    x = as.numeric(data), dims = shape, repr = "R", index1 = FALSE
  )
  counts <- Matrix::t(csr)
  rownames(counts) <- feature_names
  colnames(counts) <- cell_names
  rm(csr, data, indices, indptr)
  gc(verbose = FALSE)
  list(
    counts = counts,
    feature_names = feature_names,
    cell_names = cell_names,
    meta = metadata
  )
}

rna_src <- read_counts_h5ad(
  file.path(out_dir, "GSE243917_genesXcells_ALL_raw.h5ad")
)
atac_src <- read_counts_h5ad(
  file.path(out_dir, "GSE243917_peaksXcells_ALL_raw.h5ad")
)

rna <- rna_src$counts
atac <- atac_src$counts
rownames(rna) <- make.unique(rna_src$feature_names)
rownames(atac) <- make.unique(atac_src$feature_names)
colnames(rna) <- rna_src$cell_names
colnames(atac) <- atac_src$cell_names

common_cells <- intersect(colnames(rna), colnames(atac))
stopifnot(length(common_cells) == ncol(rna), length(common_cells) == ncol(atac))
rna <- rna[, common_cells, drop = FALSE]
atac <- atac[, common_cells, drop = FALSE]
meta <- rna_src$meta[common_cells, , drop = FALSE]

# Match the SCENT arthritis panel.  The complete GEO matrices also contain
# 654 mural cells; these are the gray/other population omitted from the
# published 30,893-cell arthritis panel.  B/plasma and T/NK are merged into
# the broad labels used in the figure.
meta$newCT <- unname(c(
  Bcell = "B", plasma = "B",
  Tcell = "T", NK = "T",
  myeloid = "Myeloid", fibroblast = "Fibroblast",
  endothelial = "Endothelial"
)[as.character(meta$ct)])
keep <- !is.na(meta$newCT)
rna <- rna[, keep, drop = FALSE]
atac <- atac[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
meta$celltype <- meta$newCT
meta$cell <- rownames(meta)
meta$donor <- sub("_.*$", "", rownames(meta))

# These matrices are raw-count inputs.  Keep integer-like accessibility
# counts, since SCENT binarizes ATAC method-side rather than requiring a
# pre-binarized RDS.
rna@x <- as.numeric(round(rna@x))
atac@x <- as.numeric(round(atac@x))
meta$nCount_RNA <- as.numeric(Matrix::colSums(rna))
meta$nFeature_RNA <- Matrix::colSums(rna > 0)
meta$nCount_ATAC <- as.numeric(Matrix::colSums(atac))
meta$nFeature_ATAC <- Matrix::colSums(atac > 0)
meta$log_nCount_RNA <- log(pmax(meta$nCount_RNA, 1))
meta$log_nCount_ATAC <- log(pmax(meta$nCount_ATAC, 1))

desired_order <- c("T", "Myeloid", "B",
                   "Fibroblast", "Endothelial")
meta$newCT <- factor(meta$newCT, levels = desired_order)
meta$celltype <- as.character(meta$newCT)
stopifnot(identical(colnames(rna), colnames(atac)))
stopifnot(identical(rownames(meta), colnames(rna)))
stopifnot(ncol(rna) == 30893L)
stopifnot(all(table(meta$newCT) > 500))

saveRDS(rna, file.path(out_dir, "Arthritis.rna.rds"), compress = "gzip")
# The 295-million-entry ATAC matrix is intentionally uncompressed: gzip
# serialization is prohibitively slow for this sparse object and the RDS is
# still much smaller than the downloaded H5AD source.
saveRDS(atac, file.path(out_dir, "Arthritis.atac.rds"), compress = FALSE)
saveRDS(meta, file.path(out_dir, "Arthritis.meta.rds"), compress = "gzip")

# Candidate cis pairs: all annotated hg38 genes whose gene body +/-500 kb
# overlaps a shared ATAC peak.  Detection-rate filtering remains method-side.
annotation <- readRDS(file.path("data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac, sep = c(":", "-"), genome = NULL
)
Signac::Annotation(chrom_assay) <- annotation
obj <- SeuratObject::CreateSeuratObject(
  counts = rna, assay = "RNA", meta.data = meta
)
obj[["ATAC"]] <- chrom_assay
gene_peak_df <- GetGenePeakPairs(
  obj, annotation = annotation, expression.assay = "RNA",
  peak.assay = "ATAC", distance = 5e5
)
gene_peak_df <- gene_peak_df[!duplicated(gene_peak_df), c("gene", "peak"),
                             drop = FALSE]
gene_peak_df$peak <- CanonicalPeakNames(gene_peak_df$peak)
saveRDS(gene_peak_df, file.path(out_dir, "Arthritis.gene_peak_df.rds"),
        compress = "gzip")

counts <- as.data.frame(table(meta$newCT), stringsAsFactors = FALSE)
colnames(counts) <- c("celltype", "n_cells")
counts <- counts[order(match(counts$celltype, desired_order)), , drop = FALSE]
write.table(counts, file.path(out_dir, "Arthritis.celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

qc_metrics <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC")
qc_summary <- do.call(rbind, lapply(qc_metrics, function(nm) {
  x <- as.numeric(meta[[nm]])
  data.frame(metric = nm, n = length(x), min = min(x),
             median = median(x), mean = mean(x), max = max(x))
}))
write.table(qc_summary, file.path(out_dir, "Arthritis.qc_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

qc_thresholds <- data.frame(
  dataset = "arthritis-tissue", tissue = "joint",
  RNA_count = "> 500 genes", ATAC_count = "> 5000",
  nucleosome = "10,000 reads (paper table)",
  TSS_enrichment = ">= 10% of reads in promoters",
  source_matrix = "GEO post-QC matrices; no second QC pass",
  stringsAsFactors = FALSE
)
write.table(qc_thresholds, file.path(out_dir, "Arthritis.qc_thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Data exploration: PCA/UMAP on RNA, colored and counted by the broad celltype.
explore <- obj
explore <- Seurat::NormalizeData(explore, verbose = FALSE)
explore <- Seurat::FindVariableFeatures(explore, nfeatures = 2000,
                                        verbose = FALSE)
explore <- Seurat::ScaleData(explore, features = Seurat::VariableFeatures(explore),
                             verbose = FALSE)
explore <- Seurat::RunPCA(explore, features = Seurat::VariableFeatures(explore),
                          npcs = 30, verbose = FALSE)
explore <- Seurat::RunUMAP(explore, dims = 1:20, n.neighbors = 30,
                           min.dist = 0.3, verbose = FALSE)
emb <- as.data.frame(SeuratObject::Embeddings(explore, "umap")[, 1:2,
                                                              drop = FALSE])
colnames(emb) <- c("UMAP1", "UMAP2")
emb$celltype <- meta$celltype[match(rownames(emb), rownames(meta))]
n_by_type <- table(emb$celltype)
emb$label <- paste0(emb$celltype, " (", n_by_type[emb$celltype], ")")
p <- ggplot(emb, aes(UMAP1, UMAP2, colour = label)) +
  geom_point(size = 0.22, alpha = 0.45) +
  theme_classic(base_size = 11) +
  labs(title = "Arthritis tissue (joint)",
       subtitle = "RNA PCA/UMAP; GEO post-QC shared cells", colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(out_dir, "Arthritis.celltype_umap.png"), p,
       width = 6.5, height = 5.0, dpi = 220)

qc_long <- do.call(rbind, lapply(qc_metrics, function(nm) {
  data.frame(celltype = meta$celltype, metric = nm,
             value = as.numeric(meta[[nm]]))
}))
q <- ggplot(qc_long, aes(celltype, value, fill = celltype)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                                        legend.position = "none") +
  labs(title = "Arthritis data exploration: count distributions", x = NULL, y = NULL)
ggsave(file.path(out_dir, "Arthritis.qc_violin.png"), q,
       width = 9, height = 6, dpi = 220)

saveRDS(explore, file.path(out_dir, "Arthritis.exploration_seurat.rds"),
        compress = "gzip")
writeLines(c(
  "Arthritis tissue benchmark preparation",
  "",
  "Source: GEO GSE243917. The complete raw gene/peak H5AD matrices are already post-QC.",
  "The paper's arthritis QC thresholds are recorded in Arthritis.qc_thresholds.tsv.",
  "The complete GEO object has 654 mural cells; these are excluded to reproduce the",
  "published 30,893-cell SCENT arthritis panel. Bcell+plasma and Tcell+NK are merged.",
  "ATAC.rds is a peak-by-cell count matrix; SCENT should binarize it method-side.",
  "Candidate gene_peak_df uses hg38 gene bodies +/-500 kb and the shared ATAC peak set.",
  "The UMAP is an RNA PCA/UMAP exploration plot, not the paper's unpublished embedding."
), file.path(out_dir, "README.md"))

cat("Prepared", ncol(rna), "cells,", nrow(rna), "RNA features,", nrow(atac),
    "ATAC peaks, and", nrow(gene_peak_df), "candidate gene-peak pairs.\n")
