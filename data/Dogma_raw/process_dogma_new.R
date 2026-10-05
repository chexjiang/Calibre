#!/usr/bin/env Rscript

# Rebuild DOGMA-seq CTRL/STIM objects from the original GSE156478 files.
# QC thresholds follow the supplied Sakaue/SCENT table:
# RNA 1,000--25,000; ATAC 1,000--100,000; nucleosome <2; TSS >1.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
})

file_argument <- sub("^--file=", "", commandArgs(FALSE)[
  grep("^--file=", commandArgs(FALSE))[[1L]]
])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

raw_dir <- file.path(project_root, "data", "Dogma_new_raw")
if (!dir.exists(raw_dir)) {
  stop("Missing ", raw_dir, "; run bash data/download_dogma_new.sh first",
       call. = FALSE)
}

annotation_file <- file.path(project_root, "data", "hg38_EnsDb_Hsapiens_v86_annotation.rds")
if (!file.exists(annotation_file)) {
  stop("Missing annotation: ", annotation_file, call. = FALSE)
}
annotation <- readRDS(annotation_file)
names(annotation) <- NULL

library_specs <- data.frame(
  condition = c("CTRL", "CTRL", "STIM", "STIM"),
  lysis = c("LLL", "DIG", "LLL", "DIG"),
  matrix_prefix = c(
    "GSM5065525_LLL_CTRL_GExp_ATAC_filtered",
    "GSM5065531_DIG_CTRL_GExp_ATAC_filtered",
    "GSM5065528_LLL_STIM_GExp_ATAC_filtered",
    "GSM5065534_DIG_STIM_GExp_ATAC_filtered"
  ),
  metrics_file = c(
    "GSM5065524_LLL_CTRL_per_barcode_metrics.csv.gz",
    "GSM5065530_DIG_CTRL_per_barcode_metrics.csv.gz",
    "GSM5065527_LLL_STIM_per_barcode_metrics.csv.gz",
    "GSM5065533_DIG_STIM_per_barcode_metrics.csv.gz"
  ),
  fragment_file = c(
    "GSM5065524_LLL_CTRL_fragments.tsv.gz",
    "GSM5065530_DIG_CTRL_fragments.tsv.gz",
    "GSM5065527_LLL_STIM_fragments.tsv.gz",
    "GSM5065533_DIG_STIM_fragments.tsv.gz"
  ),
  stringsAsFactors = FALSE
)
# The SCENT Dogma objects are the LLL libraries; keep the DIG downloads in
# Dogma_new_raw for provenance, but do not use them in these benchmark inputs.
library_specs <- library_specs[library_specs$lysis == "LLL", , drop = FALSE]

required_raw <- unlist(lapply(seq_len(nrow(library_specs)), function(i) {
  c(paste0(library_specs$matrix_prefix[[i]], c("_barcodes.tsv.gz",
    "_features.tsv.gz", "_matrix.mtx.gz")), library_specs$metrics_file[[i]],
    library_specs$fragment_file[[i]])
}))
missing_raw <- required_raw[!file.exists(file.path(raw_dir, required_raw))]
if (length(missing_raw)) {
  stop("Missing raw DOGMA files:\n", paste(missing_raw, collapse = "\n"),
       call. = FALSE)
}

read_10x <- function(prefix) {
  barcode_file <- file.path(raw_dir, paste0(prefix, "_barcodes.tsv.gz"))
  feature_file <- file.path(raw_dir, paste0(prefix, "_features.tsv.gz"))
  matrix_file <- file.path(raw_dir, paste0(prefix, "_matrix.mtx.gz"))
  barcodes <- data.table::fread(barcode_file, header = FALSE,
                                data.table = FALSE)[[1L]]
  features <- data.table::fread(feature_file, header = FALSE,
                                 data.table = FALSE)
  if (ncol(features) < 3L) {
    stop("Expected 10x feature type column in ", feature_file, call. = FALSE)
  }
  counts <- Matrix::readMM(gzfile(matrix_file))
  counts <- as(counts, "dgCMatrix")
  stopifnot(nrow(counts) == nrow(features), ncol(counts) == length(barcodes))
  feature_type <- as.character(features[[3L]])
  feature_name <- as.character(features[[2L]])
  feature_id <- as.character(features[[1L]])
  gene <- feature_type %in% c("Gene Expression", "Gene Expression")
  peak <- feature_type %in% c("Peaks", "Peaks")
  if (!any(gene) || !any(peak)) {
    stop("Could not find both Gene Expression and Peaks in ", feature_file,
         "; feature types: ", paste(unique(feature_type), collapse = ", "),
         call. = FALSE)
  }
  feature_labels <- ifelse(gene | peak, feature_name, feature_id)
  feature_labels[peak] <- gsub(":", "-", feature_labels[peak], fixed = TRUE)
  rownames(counts) <- make.unique(feature_labels)
  colnames(counts) <- barcodes
  list(
    rna = counts[gene, , drop = FALSE],
    atac = counts[peak, , drop = FALSE],
    feature_id = feature_id[peak],
    feature_name = feature_name[peak]
  )
}

read_metrics <- function(file) {
  metrics <- data.table::fread(file, showProgress = FALSE)
  names(metrics) <- make.names(names(metrics), unique = TRUE)
  metrics
}

find_column <- function(x, candidates, label) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) {
    stop("Could not find ", label, " column. Available columns: ",
         paste(names(x), collapse = ", "), call. = FALSE)
  }
  hit[[1L]]
}

parse_qc <- function(metrics, barcodes) {
  barcode_col <- find_column(metrics,
    c("barcode", "Barcode", "cell_barcode", "cell"), "barcode")
  metric_index <- match(barcodes, as.character(metrics[[barcode_col]]))
  if (anyNA(metric_index)) stop("Barcode matching failed in metrics", call. = FALSE)
  metrics <- metrics[metric_index, , drop = FALSE]
  gex_col <- find_column(metrics,
    c("gex_barcode", "GEX_barcode"), "GEX barcode")
  atac_col <- find_column(metrics,
    c("atac_barcode", "ATAC_barcode"), "ATAC barcode")
  n_atac_col <- find_column(metrics,
    c("atac_fragments", "ATAC_fragments", "atac_fragments_count",
      "atac_fragments_count_per_barcode", "nFrags"), "ATAC fragments")
  out <- data.frame(
    nCount_ATAC = as.numeric(metrics[[n_atac_col]]),
    gex_barcode = as.character(metrics[[gex_col]]),
    atac_barcode = as.character(metrics[[atac_col]]),
    row.names = barcodes,
    check.names = FALSE
  )
  out
}

processed <- lapply(seq_len(nrow(library_specs)), function(i) {
  spec <- library_specs[i, ]
  message("Reading ", spec$condition, "/", spec$lysis)
  x <- read_10x(spec$matrix_prefix)
  metrics <- read_metrics(file.path(raw_dir, spec$metrics_file))
  qc <- parse_qc(metrics, colnames(x$rna))
  if (anyDuplicated(qc$atac_barcode)) {
    stop("ATAC barcodes are not unique in ", spec$metrics_file, call. = FALSE)
  }
  # The filtered matrix and fragment file use the GEX barcode. The separate
  # ATAC barcode is retained in metadata for provenance and is not used as
  # the Seurat cell key.
  rownames(qc) <- colnames(x$rna)

  fragment_path <- file.path(raw_dir, spec$fragment_file)
  fragment_index <- paste0(fragment_path, ".tbi")
  if (!file.exists(fragment_index)) {
    if (!nzchar(Sys.which("tabix"))) stop("tabix is required to index ", fragment_path,
                                           call. = FALSE)
    status <- system2(Sys.which("tabix"), c("-p", "bed", fragment_path))
    if (!identical(status, 0L) || !file.exists(fragment_index)) {
      stop("Could not create tabix index for ", fragment_path, call. = FALSE)
    }
  }
  fragments <- Signac::CreateFragmentObject(
    path = fragment_path, cells = colnames(x$atac), validate.fragments = TRUE,
    verbose = FALSE
  )
  atac_assay <- Signac::CreateChromatinAssay(
    counts = x$atac, fragments = fragments, annotation = annotation,
    sep = c("-", "-"), validate.fragments = TRUE,
    verbose = FALSE
  )
  qc_object <- Seurat::CreateSeuratObject(
    counts = x$rna, assay = "RNA", meta.data = qc
  )
  qc_object[["ATAC"]] <- atac_assay
  DefaultAssay(qc_object) <- "ATAC"
  qc_object <- Signac::NucleosomeSignal(qc_object, assay = "ATAC", verbose = FALSE)
  qc_object <- Signac::TSSEnrichment(
    qc_object, assay = "ATAC", fast = TRUE, verbose = FALSE
  )
  qc$nucleosome_signal <- qc_object$nucleosome_signal[colnames(x$atac)]
  qc$TSS.enrichment <- qc_object$TSS.enrichment[colnames(x$atac)]
  rna_counts <- Matrix::colSums(x$rna)
  qc$nCount_RNA <- as.numeric(rna_counts)
  qc$percent.mito <- if (any(grepl("^MT-", rownames(x$rna)))) {
    100 * Matrix::colSums(x$rna[grepl("^MT-", rownames(x$rna)), , drop = FALSE]) /
      pmax(qc$nCount_RNA, 1)
  } else 0
  qc$pass_qc <- with(qc,
    nCount_RNA > 1000 & nCount_RNA < 25000 &
      nCount_ATAC > 1000 & nCount_ATAC < 100000 &
      nucleosome_signal < 2 & TSS.enrichment > 1
  )
  qc$sample <- sub("_per_barcode_metrics.*$", "", spec$metrics_file)
  qc$lysis <- spec$lysis
  qc$condition <- spec$condition
  keep <- which(qc$pass_qc)
  message("  cells ", ncol(x$rna), "; pass QC ", length(keep))
  original_cells <- colnames(x$rna)
  prefixed_cells <- paste0(spec$condition, "_", spec$lysis, "_", original_cells)
  colnames(x$rna) <- prefixed_cells
  colnames(x$atac) <- prefixed_cells
  rownames(qc) <- prefixed_cells
  list(
    rna = x$rna[, keep, drop = FALSE],
    atac = x$atac[, keep, drop = FALSE],
    meta = qc[keep, , drop = FALSE],
    feature_id = x$feature_id,
    feature_name = x$feature_name,
    all_qc = qc
  )
})

combine <- function(items, slot) {
  matrices <- lapply(items, `[[`, slot)
  all_features <- unique(unlist(lapply(matrices, rownames), use.names = FALSE))
  aligned <- lapply(matrices, function(m) {
    if (identical(rownames(m), all_features)) return(m)
    nnz_by_column <- diff(m@p)
    Matrix::sparseMatrix(
      i = match(rownames(m), all_features)[rep(seq_len(ncol(m)), nnz_by_column)],
      j = rep(seq_len(ncol(m)), nnz_by_column),
      x = m@x,
      dims = c(length(all_features), ncol(m)),
      dimnames = list(all_features, colnames(m))
    )
  })
  do.call(cbind, aligned)
}

for (condition in c("CTRL", "STIM")) {
  # SCENT's Dogma control/stim inputs correspond to the low-loss lysis (LLL)
  # libraries. DIG is downloaded and QC-computed for provenance, but is not
  # mixed into the benchmark object.
  items <- processed[vapply(processed, function(x) {
    x$meta$condition[[1L]] == condition && x$meta$lysis[[1L]] == "LLL"
  }, logical(1L))]
  rna <- combine(items, "rna")
  atac <- combine(items, "atac")
  meta <- do.call(rbind, lapply(items, `[[`, "meta"))
  rownames(meta) <- colnames(rna)
  old_meta_file <- file.path(
    project_root, "data", paste0("Dogma_", condition),
    paste0("Dogma_", condition, ".meta.rds")
  )
  if (file.exists(old_meta_file)) {
    old_meta <- readRDS(old_meta_file)
    old_barcode <- if ("cell" %in% names(old_meta)) {
      as.character(old_meta$cell)
    } else rownames(old_meta)
    old_index <- match(meta$gex_barcode, old_barcode)
    meta$predicted.id <- as.character(old_meta$predicted.id[old_index])
    meta$newCT <- as.character(old_meta$newCT[old_index])
  } else {
    meta$predicted.id <- NA_character_
    meta$newCT <- NA_character_
  }
  meta$cell <- meta$gex_barcode
  meta$batch <- meta$sample
  output_dir <- file.path(project_root, "data", paste0("Dogma_", condition, "_new"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- paste0("Dogma_", condition, "_new")
  saveRDS(rna, file.path(output_dir, paste0(prefix, ".rna.rds")))
  saveRDS(atac, file.path(output_dir, paste0(prefix, ".atac.rds")))
  saveRDS(meta, file.path(output_dir, paste0(prefix, ".meta.rds")))
  saveRDS(data.frame(
    peak = rownames(atac),
    stringsAsFactors = FALSE
  ), file.path(output_dir, paste0(prefix, ".peaks.rds")))
  data.table::fwrite(meta, file.path(output_dir, paste0(prefix, ".qc.tsv")),
                     sep = "\t", row.names = TRUE)
  data.table::fwrite(
    library_specs[library_specs$condition == condition & library_specs$lysis == "LLL",
                  , drop = FALSE],
    file.path(output_dir, paste0(prefix, ".source_manifest.tsv")), sep = "\t"
  )
  message("Wrote ", prefix, ": ", ncol(rna), " cells, ", nrow(rna),
          " genes, ", nrow(atac), " peaks")
}
