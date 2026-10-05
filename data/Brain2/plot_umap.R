suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(scales)
})

out_dir <- "data/Brain2"
rna <- readRDS(file.path(out_dir, "Brain2.rna.rds"))
meta <- readRDS(file.path(out_dir, "Brain2.meta.rds"))

stopifnot(identical(colnames(rna), rownames(meta)))
celltype <- as.character(meta$newCT)
stopifnot(!anyNA(celltype), all(nchar(celltype) > 0))

obj <- CreateSeuratObject(
  counts = rna,
  meta.data = meta,
  assay = "RNA",
  project = "Brain2",
  min.cells = 0,
  min.features = 0
)
rm(rna, meta)
gc()

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000,
                            verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, n.neighbors = 30,
               min.dist = 0.3, metric = "cosine", seed.use = 1,
               return.model = FALSE, verbose = FALSE)

obj$plot_celltype <- factor(as.character(obj$newCT), levels = sort(unique(obj$newCT)))
emb <- as.data.frame(Embeddings(obj, "umap"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
emb$cell <- rownames(emb)
emb$celltype <- as.character(obj$newCT)
emb$celltype <- factor(emb$celltype, levels = levels(obj$plot_celltype))
write.table(emb, file.path(out_dir, "Brain2.rna_umap.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

counts <- as.data.frame(table(emb$celltype), stringsAsFactors = FALSE)
names(counts) <- c("celltype", "n_cells")
counts$label <- paste0(counts$celltype, "\n(n=", format(counts$n_cells, big.mark = ","), ")")

# Label at the median UMAP position. For elongated populations this remains
# readable while preserving all cells in the plot.
centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype, data = emb, FUN = median)
centers <- merge(centers, counts, by = "celltype", sort = FALSE)

palette <- c(
  "Astrocytes" = "#E69F00",
  "Excitatory neurons" = "#3B75AF",
  "Inhibitory neurons" = "#78B7D0",
  "Microglia" = "#5B3A8E",
  "Oligodendrocytes" = "#4E9F50",
  "OPCs" = "#2F7F4F"
)
palette <- palette[levels(emb$celltype)]

p <- ggplot(emb, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 0.12, alpha = 0.55, stroke = 0) +
  geom_label(data = centers,
             aes(x = UMAP_1, y = UMAP_2, label = label), inherit.aes = FALSE,
             color = "black", fill = "white", alpha = 0.88, size = 3.0,
             linewidth = 0.15, lineheight = 0.9) +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(title = "Brain2 RNA UMAP", subtitle = "GSE214979; cell types with n > 500",
       x = "UMAP1", y = "UMAP2", color = "Cell type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Brain2.rna_umap_celltype.png"), p,
       width = 10, height = 7.5, dpi = 300, limitsize = FALSE)
ggsave(file.path(out_dir, "Brain2.rna_umap_celltype.pdf"), p,
       width = 10, height = 7.5, limitsize = FALSE)
saveRDS(obj, file.path(out_dir, "Brain2.rna_umap_seurat.rds"), compress = FALSE)
write.table(counts[, c("celltype", "n_cells")],
            file.path(out_dir, "Brain2.umap_celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Saved Brain2 UMAP plots and embedding to ", out_dir)
