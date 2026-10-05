#!/usr/bin/env Rscript

# Build one dataset-level annotation per method for a Sakaue-style comparison.
# The dataset scope is inherited from the existing Sakaue EpiLink manifest, so
# every method uses the same primary cell-type union.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
sakaue <- file.path(gwas, "SCENT_style", "UKBB", "sakaue_style")
base <- file.path(sakaue, "method_comparison")
ann_dir <- file.path(base, "annotations")
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)

source_manifest <- fread(file.path(gwas, "SCENT_style", "config", "annotation_manifest.tsv"))
sakaue_manifest <- fread(file.path(sakaue, "sakaue_annotation_manifest.tsv"))
scope <- sakaue_manifest[correction == FALSE,
                         .(dataset, celltypes, n_celltypes)]
scope <- unique(scope, by = "dataset")
scope[, celltype_vec := strsplit(celltypes, ";", fixed = TRUE)]

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L)
    return(data.table(chr=character(), start=integer(), end=integer()))
  x <- fread(path, header = FALSE, col.names = c("chr", "start0", "end"))
  x[, .(chr, start = as.integer(start0) + 1L, end = as.integer(end))]
}

reduce_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand = TRUE)
  data.table(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
}

write_bed <- function(x, path) {
  y <- copy(x)
  y[, start0 := pmax(start - 1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep = "\t", col.names = FALSE)
}

method_info <- data.table(
  method = c("epilink_nb", "epilink_poisson", "scmultimap",
             "signac_corrected", "signac_uncorrected"),
  method_label = c("EpiLink-NB (correction = FALSE)",
                   "EpiLink-Poisson (correction = FALSE)",
                   "scMultiMap", "Signac (corrected)",
                   "Signac (uncorrected)"),
  correction = c(FALSE, FALSE, NA, TRUE, FALSE)
)

rows <- list(); k <- 0L
for (i in seq_len(nrow(scope))) {
  ds <- scope$dataset[i]
  cells <- scope$celltype_vec[[i]]
  for (j in seq_len(nrow(method_info))) {
    mi <- method_info[j]
    id <- if (mi$method %chin% c("epilink_nb", "epilink_poisson"))
      paste0(ds, "__", mi$method, "__correction_F") else
      paste0(ds, "__", mi$method, "__primary")

    if (mi$method %chin% c("epilink_nb", "epilink_poisson")) {
      src_id <- paste0(ds, "__", mi$method, "__correction_F")
      bed <- read_bed(file.path(sakaue, "annotations", paste0(src_id, ".hg38.bed")))
      n_pairs <- sakaue_manifest[annotation_id == src_id, n_pairs][1]
      n_native <- sakaue_manifest[annotation_id == src_id, n_native_intervals][1]
    } else {
      z <- source_manifest[dataset == ds & celltype %chin% cells & method == mi$method]
      bed <- rbindlist(lapply(z$annotation_id, function(a)
        read_bed(file.path(gwas, "SCENT_style", "testing", "annotations", paste0(a, ".hg38.bed")))),
        fill = TRUE)
      bed <- reduce_bed(bed)
      n_pairs <- sum(z$n_pairs, na.rm = TRUE)
      n_native <- sum(z$n_native_intervals, na.rm = TRUE)
    }
    bed <- reduce_bed(bed)
    write_bed(bed, file.path(ann_dir, paste0(id, ".hg38.bed")))
    k <- k + 1L
    rows[[k]] <- data.table(
      annotation_id = id, dataset = ds, celltype = "ALL_PRIMARY",
      celltype_scope = "ALL_PRIMARY",
      n_celltypes = length(cells), celltypes = paste(cells, collapse = ";"),
      method = mi$method, method_label = mi$method_label,
      correction = mi$correction, fdr = 0.1, source_build = "hg38",
      n_pairs = n_pairs, n_native_intervals = n_native,
      n_hg38_intervals = nrow(bed), n_hg19_intervals = NA_integer_,
      lift_fraction_hg38 = 1, lift_fraction_hg19 = NA_real_
    )
  }
}

manifest <- rbindlist(rows, fill = TRUE)
fwrite(manifest, file.path(base, "sakaue_annotation_manifest.tsv"), sep = "\t")
fwrite(data.table(
  dataset = scope$dataset,
  celltypes = scope$celltypes,
  n_celltypes = scope$n_celltypes,
  note = "All methods are pooled over the same Sakaue-style primary cell-type scope"
), file.path(base, "dataset_scope.tsv"), sep = "\t")
message("Wrote ", nrow(manifest), " annotations for ", uniqueN(manifest$dataset),
        " datasets to ", base)
