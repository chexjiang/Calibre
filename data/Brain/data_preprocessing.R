# Prepare the Brain multiome benchmark inputs used by Sakaue et al. (2024).
#
# The GEO RDS is a doubly-compressed RDS: gzcon(gzfile(...)) is intentional.
# The source object contains the original QC metrics, cell labels, reductions,
# RNA counts, and two ATAC peak spaces.  We use the source ATAC assay as the
# shared peak-by-cell matrix for the benchmark.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(ggplot2)
})

source(file.path("R", "helper.R"))

input <- file.path("data", "Brain", "GSE193240_HS_GM_Ctrl_multiomics_UCSC.rds.gz")
out_dir <- file.path("data", "Brain")

read_nested_gz_rds <- function(path) {
  con <- gzcon(gzfile(path, open = "rb"))
  on.exit(close(con), add = TRUE)
  readRDS(con)
}

source_object <- read_nested_gz_rds(input)
stopifnot(inherits(source_object, "Seurat"))
stopifnot(all(c("RNA", "ATAC") %in% SeuratObject::Assays(source_object)))

source_labels <- c("EXCNEU_CTR", "ASTRO_CTR", "OLIGO_CTR")
label_map <- c(
  EXCNEU_CTR = "Excitatory neurons",
  ASTRO_CTR = "Astrocytes",
  OLIGO_CTR = "Oligodendroglia"
)

source_meta <- source_object@meta.data
keep <- rownames(source_meta)[source_meta$broadcelltypes %in% source_labels]
stopifnot(length(keep) == 4272L)

rna <- SeuratObject::GetAssayData(
  source_object, assay = "RNA", layer = "counts"
)[, keep, drop = FALSE]
atac <- SeuratObject::GetAssayData(
  source_object, assay = "ATAC", layer = "counts"
)[, keep, drop = FALSE]
meta <- source_meta[keep, , drop = FALSE]

stopifnot(identical(colnames(rna), colnames(atac)))
stopifnot(identical(rownames(meta), colnames(rna)))
stopifnot(identical(as.character(meta$broadcelltypes),
                    as.character(source_labels[match(meta$broadcelltypes, source_labels)])))

meta$newCT <- unname(label_map[as.character(meta$broadcelltypes)])
meta$celltype <- meta$newCT
meta$cell <- rownames(meta)
meta$log_nCount_RNA <- log(pmax(as.numeric(meta$nCount_RNA), 1))
meta$log_nCount_ATAC <- log(pmax(as.numeric(meta$nCount_ATAC), 1))

# The source metadata are retained without a second QC pass.  These checks
# document that its count summaries agree with the matrices exported here.
stopifnot(all(as.numeric(Matrix::colSums(rna)) == as.numeric(meta$nCount_RNA)))
stopifnot(all(as.numeric(Matrix::colSums(atac)) == as.numeric(meta$nCount_ATAC)))

saveRDS(rna, file.path(out_dir, "Brain.rna.rds"), compress = "gzip")
saveRDS(atac, file.path(out_dir, "Brain.atac.rds"), compress = "gzip")
saveRDS(meta, file.path(out_dir, "Brain.meta.rds"), compress = "gzip")


# Candidate cis pairs are generated from the common ATAC peak space and the
# source hg38 gene annotation.  Detection-rate filtering remains method-side
# so EpiLink and SCENT can apply their documented, method-specific filters.
gene_peak_df <- GetGenePeakPairs(
  source_object,
  expression.assay = "RNA",
  peak.assay = "ATAC",
  distance = 5e5
)
gene_peak_df <- gene_peak_df[!duplicated(gene_peak_df), c("gene", "peak"), drop = FALSE]
saveRDS(gene_peak_df,
        file.path(out_dir, "Brain.gene_peak_df.rds"),
        compress = "gzip")

counts <- as.data.frame(table(meta$newCT), stringsAsFactors = FALSE)
colnames(counts) <- c("celltype", "n_cells")
counts <- counts[order(match(counts$celltype, unname(label_map))), , drop = FALSE]
write.table(counts, file.path(out_dir, "Brain.celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

qc_names <- intersect(
  c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC",
    "nucleosome_signal", "TSS.enrichment", "FRiP", "blacklist_fraction"),
  colnames(meta)
)
qc_summary <- do.call(rbind, lapply(qc_names, function(nm) {
  x <- as.numeric(meta[[nm]])
  finite_x <- x[is.finite(x)]
  if (!length(finite_x)) {
    return(data.frame(
      metric = nm, n = 0L, min = NA_real_, median = NA_real_,
      mean = NA_real_, max = NA_real_, stringsAsFactors = FALSE
    ))
  }
  data.frame(
    metric = nm,
    n = length(finite_x),
    min = min(finite_x),
    median = median(finite_x),
    mean = mean(finite_x),
    max = max(finite_x),
    stringsAsFactors = FALSE
  )
}))
write.table(qc_summary, file.path(out_dir, "Brain.qc_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Reproduce the Brain panel using the precomputed source WNN UMAP.
if ("wnn.umap" %in% names(source_object@reductions)) {
  umap <- as.data.frame(SeuratObject::Embeddings(source_object, "wnn.umap")[keep, 1:2, drop = FALSE])
  colnames(umap) <- c("UMAP1", "UMAP2")
  umap$celltype <- meta$newCT[match(rownames(umap), rownames(meta))]
  n_by_type <- table(umap$celltype)
  umap$label <- paste0(umap$celltype, " (", n_by_type[umap$celltype], ")")
  p <- ggplot(umap, aes(UMAP1, UMAP2, colour = label)) +
    geom_point(size = 0.35, alpha = 0.55) +
    theme_classic(base_size = 11) +
    labs(title = "Brain multiome: Sakaue et al. Brain subset",
         subtitle = "Precomputed WNN UMAP; source QC retained",
         colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ggsave(file.path(out_dir, "Brain.wnn_umap.png"), p,
         width = 6, height = 4.5, dpi = 220)
}

cat("Prepared", ncol(rna), "cells,", nrow(rna), "RNA features,",
    nrow(atac), "ATAC peaks, and", nrow(gene_peak_df),
    "candidate gene-peak pairs.\n")
