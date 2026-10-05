suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(irlba); library(uwot)})
out_dir <- file.path("data", "Pituitary")
rna <- readRDS(file.path(out_dir, "Pituitary.rna.rds"))
meta <- readRDS(file.path(out_dir, "Pituitary.meta.rds"))
stopifnot(identical(colnames(rna), rownames(meta)))
lib <- Matrix::colSums(rna)
scale_factor <- 1e4 / pmax(lib, 1)
for (j in seq_len(ncol(rna))) {
  if (rna@p[j] < rna@p[j + 1]) {
    idx <- (rna@p[j] + 1L):rna@p[j + 1L]
    rna@x[idx] <- log1p(rna@x[idx] * scale_factor[j])
  }
}
variance <- Matrix::rowMeans(rna^2) - Matrix::rowMeans(rna)^2
hv <- order(variance, decreasing = TRUE)[seq_len(min(2000, nrow(rna)))]
pca <- irlba::prcomp_irlba(t(rna[hv, , drop = FALSE]), n = 20, center = TRUE, scale. = FALSE)
um <- uwot::umap(pca$x, n_neighbors = 30, min_dist = .3, metric = "cosine", n_threads = 2, verbose = FALSE)
emb <- as.data.frame(um); colnames(emb) <- c("UMAP1", "UMAP2")
emb$celltype <- meta$newCT
n_by_type <- table(emb$celltype)
emb$label <- paste0(emb$celltype, " (", n_by_type[emb$celltype], ")")
p <- ggplot(emb, aes(UMAP1, UMAP2, colour = label)) +
  geom_point(size = .28, alpha = .5) + theme_classic(base_size = 11) +
  labs(title = "Pituitary", subtitle = "RNA PCA/UMAP after count QC", colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(out_dir, "Pituitary.celltype_umap.png"), p, width = 6.5, height = 5, dpi = 220)
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC")
qlong <- do.call(rbind, lapply(qc_metrics, function(nm)
  data.frame(celltype = meta$newCT, metric = nm, value = as.numeric(meta[[nm]]))))
q <- ggplot(qlong, aes(celltype, value, fill = celltype)) +
  geom_violin(scale = "width", trim = TRUE) + facet_wrap(~metric, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
  labs(title = "Pituitary data exploration: count distributions", x = NULL, y = NULL)
ggsave(file.path(out_dir, "Pituitary.qc_violin.png"), q, width = 9, height = 6, dpi = 220)
saveRDS(list(pca = pca$x, umap = emb), file.path(out_dir, "Pituitary.exploration_embedding.rds"), compress = "gzip")
cat("Wrote pituitary exploration for", nrow(emb), "cells.\n")
