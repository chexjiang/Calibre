suppressPackageStartupMessages(library(ggplot2))
out_dir <- file.path("data", "Arthritis")
x <- readRDS(file.path(out_dir, "Arthritis.exploration_embedding.rds"))
meta <- readRDS(file.path(out_dir, "Arthritis.meta.rds"))
emb <- as.data.frame(x$umap[, 1:2, drop = FALSE])
colnames(emb) <- c("UMAP1", "UMAP2")
stopifnot(nrow(emb) == nrow(meta))
emb$celltype <- meta$celltype
n_by_type <- table(emb$celltype)
emb$label <- paste0(emb$celltype, " (", n_by_type[emb$celltype], ")")
p <- ggplot(emb, aes(UMAP1, UMAP2, colour = label)) +
  geom_point(size = .22, alpha = .45) + theme_classic(base_size = 11) +
  labs(title = "Arthritis tissue (joint)", subtitle = "Fast RNA PCA/UMAP exploration; GEO post-QC cells", colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(out_dir, "Arthritis.celltype_umap.png"), p, width = 6.5, height = 5, dpi = 220)
cat("Replotted arthritis UMAP for", nrow(emb), "cells.\n")
