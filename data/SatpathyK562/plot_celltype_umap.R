#!/usr/bin/env Rscript

# Plot the RNA UMAP of QC-passing Satpathy K562 cells, coloured by cell type.
# SatpathyK562 is a K562-only subset, so the plot intentionally has one
# cell-type label; the UMAP still documents the within-line structure used by
# downstream benchmark inputs.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else
  "data/SatpathyK562/plot_celltype_umap.R"
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
data_dir <- file.path(project_root, "data", "SatpathyK562")

rna <- readRDS(file.path(data_dir, "SatpathyK562.rna.rds"))
meta <- readRDS(file.path(data_dir, "SatpathyK562.meta.rds"))
stopifnot(identical(colnames(rna), rownames(meta)))

object <- CreateSeuratObject(counts = rna, meta.data = meta)
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst",
                               nfeatures = min(3000L, nrow(object)),
                               verbose = FALSE)
object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
object <- RunPCA(object, features = VariableFeatures(object), npcs = 30,
                 verbose = FALSE)
n_pcs <- min(30L, ncol(Embeddings(object, "pca")))
object <- RunUMAP(object, reduction = "pca", dims = seq_len(n_pcs),
                  reduction.name = "umap_rna", reduction.key = "RNAUMAP_",
                  verbose = FALSE)

plot_data <- as.data.frame(Embeddings(object, "umap_rna")[, 1:2, drop = FALSE])
colnames(plot_data) <- c("UMAP1", "UMAP2")
plot_data$celltype <- meta[rownames(plot_data), "newCT"]
n_by_type <- table(plot_data$celltype)
plot_data$celltype_label <- paste0(plot_data$celltype, " (",
                                   n_by_type[plot_data$celltype], ")")

p <- ggplot(plot_data, aes(UMAP1, UMAP2, colour = celltype_label)) +
  geom_point(size = 0.35, alpha = 0.65, stroke = 0) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Satpathy K562: QC-passing cells",
       subtitle = "RNA PCA/UMAP; colour indicates cell type",
       x = "UMAP 1", y = "UMAP 2", colour = "Cell type") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave(file.path(data_dir, "SatpathyK562.celltype_umap.pdf"), p,
       width = 7, height = 5.5, useDingbats = FALSE)
ggsave(file.path(data_dir, "SatpathyK562.celltype_umap.png"), p,
       width = 7, height = 5.5, dpi = 300)
write.table(as.data.frame(table(meta$newCT)),
            file.path(data_dir, "SatpathyK562.celltype_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved SatpathyK562.celltype_umap.pdf/.png")
