# Prepare the pituitary benchmark inputs used by Sakaue et al. (2024).
#
# GSE178453 is a 10x multiome sample.  Its filtered H5 contains both the
# gene-expression and Cell Ranger peak-by-cell matrices; the companion
# fragments file is retained for QC/provenance.  The supplied Peaks matrix is
# therefore already peak-called, so this benchmark export does not call peaks
# a second time.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(rhdf5)
})

source(file.path("R", "helper.R"))

out_dir <- file.path("data", "Pituitary")
raw_dir <- file.path(out_dir, "raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(raw_dir,
  "GSM5391776_HF1M-Cryo_filtered_feature_bc_matrix.h5")
fragment_file <- file.path(raw_dir,
  "GSM5391776_HF1M-Cryo_atac_fragments.tsv.gz")
stopifnot(file.exists(matrix_file), file.exists(fragment_file))

read_10x_modalities <- function(path) {
  feature_type <- as.character(rhdf5::h5read(path, "matrix/features/feature_type"))
  feature_name <- as.character(rhdf5::h5read(path, "matrix/features/name"))
  barcodes <- as.character(rhdf5::h5read(path, "matrix/barcodes"))
  shape <- as.integer(rhdf5::h5read(path, "matrix/shape"))
  data <- rhdf5::h5read(path, "matrix/data")
  indices <- rhdf5::h5read(path, "matrix/indices")
  indptr <- rhdf5::h5read(path, "matrix/indptr")
  x <- Matrix::sparseMatrix(i=as.integer(indices),
                            p=as.integer(indptr),
                            x=as.numeric(data), dims=shape,
                            repr="C", index1=FALSE)
  list(
    `Gene Expression` = x[feature_type == "Gene Expression",,drop=FALSE],
    Peaks = x[feature_type == "Peaks",,drop=FALSE],
    gene_names = feature_name[feature_type == "Gene Expression"],
    peak_names = feature_name[feature_type == "Peaks"],
    barcodes = barcodes
  )
}
tenx <- read_10x_modalities(matrix_file)
rna <- tenx[["Gene Expression"]]
atac <- tenx[["Peaks"]]
rownames(rna) <- tenx$gene_names; rownames(atac) <- tenx$peak_names
colnames(rna) <- tenx$barcodes; colnames(atac) <- tenx$barcodes
rownames(rna) <- make.unique(rownames(rna))
rownames(atac) <- make.unique(rownames(atac))

common_cells <- intersect(colnames(rna), colnames(atac))
rna <- rna[, common_cells, drop = FALSE]
atac <- atac[, common_cells, drop = FALSE]

annotation <- readRDS(file.path("data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
# Create the tabix index before CreateChromatinAssay validates fragments.
if (!file.exists(paste0(fragment_file, ".tbi"))) {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Rsamtools is required to index the pituitary fragment file")
  }
  Rsamtools::indexTabix(fragment_file, format = "bed")
}
chrom_assay <- Signac::CreateChromatinAssay(
  counts = atac, fragments = NULL, genome = NULL,
  annotation = annotation, sep = c(":", "-")
)
obj <- SeuratObject::CreateSeuratObject(counts = rna, assay = "RNA")
obj[["ATAC"]] <- chrom_assay

# Apply the count thresholds from Supplementary Table 9.  The fragments file
# is retained under raw/ for optional full fragment-level QC; this export uses
# the public filtered peak matrix and does not rescan 3 GB of fragments.
md <- obj[[]]
md$qc_pass <- with(md,
  nCount_RNA > 1000 & nCount_RNA < 25000 &
  nCount_ATAC > 1000 & nCount_ATAC < 100000)
keep <- rownames(md)[md$qc_pass]
obj <- subset(obj, cells = keep)

# Assign the broad labels used in the Sakaue panel from canonical pituitary
# hormone markers.  The supplied multiome sample is unlabeled in GEO.  Scores
# are computed on normalized RNA and the strongest lineage score is retained.
obj <- Seurat::NormalizeData(obj, assay = "RNA", verbose = FALSE)
marker_sets <- list(
  Somatotropes = c("GH1", "GH2", "GHRHR"),
  Lactotropes = c("PRL", "PIT1", "POU1F1"),
  Corticotropes = c("POMC", "TBX19"),
  Gonadotropes = c("LHB", "FSHB", "GNRHR"),
  Thyrotropes = c("TSHB", "TRHR")
)
data_norm <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
score_mat <- sapply(marker_sets, function(gs) {
  present <- intersect(gs, rownames(data_norm))
  if (!length(present)) return(rep(0, ncol(data_norm)))
  Matrix::colMeans(data_norm[present, , drop = FALSE])
})
score_mat <- as.matrix(score_mat)
obj$celltype <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
obj$newCT <- obj$celltype

desired_order <- c("Somatotropes", "Corticotropes", "Lactotropes",
                   "Thyrotropes", "Gonadotropes")
obj$celltype <- factor(obj$celltype, levels = desired_order)
obj$newCT <- as.character(obj$celltype)
counts <- table(obj$celltype)
keep_types <- names(counts)[counts > 500]
keep_types <- intersect(desired_order, keep_types)
obj <- subset(obj, cells = colnames(obj)[obj$celltype %in% keep_types])
obj$celltype <- droplevels(obj$celltype)
obj$newCT <- as.character(obj$celltype)
stopifnot(all(table(obj$celltype) > 500))

rna <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
atac <- SeuratObject::GetAssayData(obj, assay = "ATAC", layer = "counts")
meta <- obj[[]]
meta$cell <- rownames(meta)
meta$donor <- "HF1M-Cryo"
meta$log_nCount_RNA <- log(pmax(meta$nCount_RNA, 1))
meta$log_nCount_ATAC <- log(pmax(meta$nCount_ATAC, 1))
meta$celltype <- as.character(meta$celltype)
meta$newCT <- meta$celltype
rna <- rna[, rownames(meta), drop = FALSE]
atac <- atac[, rownames(meta), drop = FALSE]
meta <- meta[colnames(rna), , drop = FALSE]

saveRDS(rna, file.path(out_dir, "Pituitary.rna.rds"), compress = "gzip")
# Avoid very slow gzip serialization for the 95k-peak sparse matrix.
saveRDS(atac, file.path(out_dir, "Pituitary.atac.rds"), compress = FALSE)
saveRDS(meta, file.path(out_dir, "Pituitary.meta.rds"), compress = "gzip")

obj_for_pairs <- SeuratObject::CreateSeuratObject(
  counts = rna, assay = "RNA", meta.data = meta
)
obj_for_pairs[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac, fragments = NULL, genome = NULL, annotation = annotation,
  sep = c(":", "-")
)
gene_peak_df <- GetGenePeakPairs(
  obj_for_pairs, annotation = annotation, expression.assay = "RNA",
  peak.assay = "ATAC", distance = 5e5
)
gene_peak_df <- gene_peak_df[!duplicated(gene_peak_df), c("gene", "peak"),
                             drop = FALSE]
saveRDS(gene_peak_df, file.path(out_dir, "Pituitary.gene_peak_df.rds"),
        compress = "gzip")

counts_df <- as.data.frame(table(meta$newCT), stringsAsFactors = FALSE)
colnames(counts_df) <- c("celltype", "n_cells")
counts_df <- counts_df[order(match(counts_df$celltype, desired_order)), , drop = FALSE]
write.table(counts_df, file.path(out_dir, "Pituitary.celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

qc_metrics <- intersect(c("nCount_RNA", "nFeature_RNA", "nCount_ATAC",
                          "nFeature_ATAC"),
                        colnames(meta))
qc_summary <- do.call(rbind, lapply(qc_metrics, function(nm) {
  x <- as.numeric(meta[[nm]])
  data.frame(metric = nm, n = length(x), min = min(x),
             median = median(x), mean = mean(x), max = max(x))
}))
write.table(qc_summary, file.path(out_dir, "Pituitary.qc_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
qc_thresholds <- data.frame(
  dataset = "Pituitary", tissue = "pituitary",
  RNA_count = ">1000 & <25000", ATAC_count = ">1000 & <100000",
  nucleosome = "<2", TSS_enrichment = ">1",
  source_matrix = "10x filtered gene/peak matrix; count QC applied; fragments retained",
  stringsAsFactors = FALSE
)
write.table(qc_thresholds, file.path(out_dir, "Pituitary.qc_thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# RNA exploration UMAP and QC violin plot.
explore <- Seurat::NormalizeData(obj_for_pairs, verbose = FALSE)
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
  geom_point(size = 0.28, alpha = 0.5) + theme_classic(base_size = 11) +
  labs(title = "Pituitary", subtitle = "RNA PCA/UMAP after QC",
       colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(out_dir, "Pituitary.celltype_umap.png"), p,
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
  labs(title = "Pituitary data exploration: QC distributions", x = NULL, y = NULL)
ggsave(file.path(out_dir, "Pituitary.qc_violin.png"), q,
       width = 9, height = 7, dpi = 220)
saveRDS(explore, file.path(out_dir, "Pituitary.exploration_seurat.rds"),
        compress = "gzip")
explore_embedding <- list(
  pca = SeuratObject::Embeddings(explore, "pca"),
  umap = data.frame(
    SeuratObject::Embeddings(explore, "umap")[, 1:2, drop = FALSE],
    celltype = meta$celltype[match(
      rownames(SeuratObject::Embeddings(explore, "umap")),
      rownames(meta)
    )],
    row.names = rownames(SeuratObject::Embeddings(explore, "umap")),
    check.names = FALSE
  )
)
saveRDS(explore_embedding,
        file.path(out_dir, "Pituitary.exploration_embedding.rds"),
        compress = "gzip")

writeLines(c(
  "Pituitary benchmark preparation",
  "",
  "Source: GEO GSE178453 / GSM5391776 (HF1M-Cryo human pituitary multiome).",
  "The 10x H5 contains Gene Expression and a filtered peak-by-cell matrix.",
  "The peak matrix is already peak-called by Cell Ranger; the original fragments",
  "file is retained under raw/ for optional fragment-level nucleosome/TSS QC.",
  "Because GEO does not provide broad labels for this sample, celltype/newCT are",
  "marker-score annotations based on canonical hormone/lineage genes.",
  "Only cell types with >500 cells are exported. Candidate pairs use hg38 +/-500 kb."
), file.path(out_dir, "README.md"))

cat("Prepared", ncol(rna), "cells,", nrow(rna), "RNA features,", nrow(atac),
    "ATAC peaks, and", nrow(gene_peak_df), "candidate gene-peak pairs.\n")
