# Import the public K562 10x Multiome object used as a practical processed
# starting point for the Xu K562 dataset. The source object is an h5mu file
# and stores each modality as a cell-by-feature CSC matrix; EpiLink stores
# feature-by-cell dgCMatrix objects instead.

suppressPackageStartupMessages({
  library(Matrix)
  library(reticulate)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "data/XuK562/import_h5mu.R"
project_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
data_dir <- file.path(project_root, "data", "XuK562")

source_h5mu <- Sys.getenv("XUK562_H5MU", unset = "/tmp/K562.h5mu")
if (!file.exists(source_h5mu)) {
  stop("Missing h5mu source. Set XUK562_H5MU to the downloaded K562.h5mu file.")
}

python_bin <- Sys.getenv("EPILINK_PYTHON", unset = "/home/chenxin/miniforge3/bin/python")
use_python(python_bin, required = TRUE)
h5py <- import("h5py")
numpy <- import("numpy")

decode_dataset <- py_eval("lambda ds: [x.decode() for x in ds[:]]")

read_counts <- function(mod) {
  x_group <- mod[["layers"]][["counts"]]
  shape <- as.integer(py_to_r(numpy$asarray(x_group$attrs[["shape"]], dtype = "int32")))
  values <- as.numeric(py_to_r(numpy$asarray(x_group[["data"]], dtype = "float64")))
  indices <- as.integer(py_to_r(numpy$asarray(x_group[["indices"]], dtype = "int32")))
  indptr <- as.integer(py_to_r(numpy$asarray(x_group[["indptr"]], dtype = "int32")))
  matrix_cell_feature <- new(
    "dgCMatrix", Dim = shape, p = indptr, i = indices, x = values
  )
  matrix_feature_cell <- Matrix::t(matrix_cell_feature)
  rownames(matrix_feature_cell) <- py_to_r(decode_dataset(mod[["var"]][["_index"]]))
  colnames(matrix_feature_cell) <- py_to_r(decode_dataset(mod[["obs"]][["_index"]]))
  matrix_feature_cell
}

h5 <- h5py$File(source_h5mu, "r")
message("Reading RNA counts from ", source_h5mu)
rna <- read_counts(h5[["mod"]][["rna"]])
message("RNA loaded: ", nrow(rna), " x ", ncol(rna))
message("Reading ATAC counts")
atac <- read_counts(h5[["mod"]][["atac"]])
message("ATAC loaded: ", nrow(atac), " x ", ncol(atac))

# Apply the Xu K562 QC thresholds used in the benchmarking table:
# RNA > 1,000; ATAC > 200; nucleosome signal < 2; TSS enrichment > 2.
obs <- h5[["mod"]][["rna"]][["obs"]]
barcodes <- py_to_r(decode_dataset(obs[["_index"]]))
read_obs_numeric <- function(field) {
  as.numeric(py_to_r(numpy$asarray(obs[[field]], dtype = "float64")))
}
meta <- data.frame(
  barcodes = barcodes,
  newCT = "K562",
  nCount_RNA = read_obs_numeric("nCount_RNA"),
  nCount_ATAC = read_obs_numeric("nCount_ATAC"),
  NucleosomeRatio = read_obs_numeric("NucleosomeRatio"),
  TSSEnrichment = read_obs_numeric("TSSEnrichment"),
  row.names = barcodes,
  stringsAsFactors = FALSE
)
meta$qc_pass <- with(meta,
  is.finite(nCount_RNA) & nCount_RNA > 1000 &
    is.finite(nCount_ATAC) & nCount_ATAC > 200 &
    is.finite(NucleosomeRatio) & NucleosomeRatio < 2 &
    is.finite(TSSEnrichment) & TSSEnrichment > 2
)
message("Xu K562 QC retained ", sum(meta$qc_pass), " / ", nrow(meta), " cells")

if (!identical(colnames(rna), colnames(atac))) {
  common <- intersect(colnames(rna), colnames(atac))
  rna <- rna[, common, drop = FALSE]
  atac <- atac[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
}

keep_cells <- intersect(colnames(rna), rownames(meta)[meta$qc_pass])
rna <- rna[, keep_cells, drop = FALSE]
atac <- atac[, keep_cells, drop = FALSE]
meta <- meta[keep_cells, , drop = FALSE]

# The h5mu object uses Ensembl IDs. EpiLink candidate pairs use the gene
# symbols in the hg38 annotation, so map IDs once here and retain one row per
# symbol to avoid duplicated feature names.
annotation <- readRDS(file.path(project_root, "data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
message("Mapping RNA Ensembl IDs to hg38 gene symbols")
feature_ids <- rownames(rna)
symbols <- as.character(annotation$gene_name)[match(feature_ids, as.character(annotation$gene_id))]
# Cell Ranger can emit Ensembl IDs for some features and gene symbols for
# others. Keep already-symbol-named features and map only ENSG identifiers.
symbols[!grepl("^ENSG", feature_ids)] <- feature_ids[!grepl("^ENSG", feature_ids)]
keep <- !is.na(symbols) & nzchar(symbols)
keep <- keep & !duplicated(symbols)
rna <- rna[keep, , drop = FALSE]
rownames(rna) <- symbols[keep]
rm(annotation, feature_ids, symbols, keep)
gc()

meta$barcodes <- rownames(meta)

saveRDS(rna, file.path(data_dir, "XuK562.rna.rds"), compress = FALSE)
saveRDS(atac, file.path(data_dir, "XuK562.atac.rds"), compress = FALSE)
saveRDS(meta, file.path(data_dir, "XuK562.meta.rds"), compress = FALSE)

message("Saved XuK562 RNA ", nrow(rna), " x ", ncol(rna),
        "; ATAC ", nrow(atac), " x ", ncol(atac),
        "; cells with paired modalities.")
