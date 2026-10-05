suppressPackageStartupMessages(library(ggplot2))
out_dir <- file.path("data", "Arthritis")
meta <- readRDS(file.path(out_dir, "Arthritis.meta.rds"))
map <- c(T = "T cells/NK cells", Myeloid = "Myeloid cells", B = "B cells")
old <- as.character(meta$newCT)
renamed <- unname(map[old])
meta$newCT <- ifelse(is.na(renamed), old, renamed)
meta$celltype <- meta$newCT
saveRDS(meta, file.path(out_dir, "Arthritis.meta.rds"), compress = "gzip")
desired_order <- c("T cells/NK cells", "Myeloid cells", "B cells", "Fibroblast", "Endothelial")
counts <- as.data.frame(table(meta$newCT), stringsAsFactors = FALSE)
colnames(counts) <- c("celltype", "n_cells")
counts <- counts[order(match(counts$celltype, desired_order)), , drop = FALSE]
write.table(counts, file.path(out_dir, "Arthritis.celltype_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
x <- readRDS(file.path(out_dir, "Arthritis.exploration_embedding.rds"))
emb <- as.data.frame(x$umap[, 1:2, drop = FALSE]); colnames(emb) <- c("UMAP1", "UMAP2")
emb$celltype <- meta$celltype
n_by_type <- table(emb$celltype); emb$label <- paste0(emb$celltype, " (", n_by_type[emb$celltype], ")")
p <- ggplot(emb, aes(UMAP1, UMAP2, colour = label)) + geom_point(size = .22, alpha = .45) +
  theme_classic(base_size = 11) + labs(title = "Arthritis tissue (joint)",
  subtitle = "Fast RNA PCA/UMAP exploration; GEO post-QC cells", colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(out_dir, "Arthritis.celltype_umap.png"), p, width = 6.5, height = 5, dpi = 220)
