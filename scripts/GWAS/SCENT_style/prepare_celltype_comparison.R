#!/usr/bin/env Rscript

# Build cell-type-specific hg38 annotations for the same five-method comparison
# used in method_comparison. EpiLink is re-filtered from its stored object with
# correction=FALSE so it is directly comparable to the pooled analysis.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
source_results <- source_results_paths(root)
analysis_arg <- if (length(args) >= 2L) args[[2]] else "celltype"
analysis_name <- if (analysis_arg %in% c("celltype_comparison", "celltype_comparison_custom", "celltype_main", "celltype"))
  "celltype" else analysis_arg
scope_path <- if (length(args) >= 3L) normalizePath(args[[3]], mustWork = TRUE) else
  file.path(source_results$shared, "custom_scope.tsv")
gwas <- file.path(root, "results", "validation", "GWAS")
sakaue <- source_results$root
base <- if (analysis_name == "celltype") source_results$celltype else
  file.path(sakaue, "celltype", analysis_name)
ann_dir <- file.path(base, "annotations")
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "R", "EpiLink.R"))

scope <- fread(scope_path)
scope[, celltype_vec := strsplit(celltypes, ";", fixed = TRUE)]
source_manifest <- fread(file.path(gwas, "SCENT_style", "config", "annotation_manifest.tsv"))
dataset_manifest <- fread(file.path(gwas, "SCENT_style", "config", "datasets.tsv"))

import_chain_gz <- function(path) {
  tmp <- tempfile(fileext = ".over.chain")
  con <- gzfile(path, open = "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  writeLines(lines, tmp)
  import.chain(tmp)
}
chain_19_38 <- import_chain_gz(file.path(gwas, "reference", "liftover", "hg19ToHg38.over.chain.gz"))

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L)
    return(data.table(chr=character(), start=integer(), end=integer()))
  x <- fread(path, header = FALSE, col.names = c("chr", "start0", "end"))
  x[, .(chr, start = as.integer(start0) + 1L, end = as.integer(end))]
}
parse_peaks <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(data.table(chr=character(), start=integer(), end=integer()))
  z <- tstrsplit(x, "-", fixed = TRUE)
  if (length(z) != 3L) stop("Peak IDs must have chr-start-end form")
  ans <- data.table(chr=z[[1]], start=as.integer(z[[2]]), end=as.integer(z[[3]]))
  ans <- ans[grepl("^chr([1-9]|1[0-9]|2[0-2]|X)$", chr) & !is.na(start) & !is.na(end) & start <= end]
  unique(ans)
}
reduce_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand=TRUE)
  data.table(chr=as.character(seqnames(gr)), start=start(gr), end=end(gr))
}
lift_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- GRanges(x$chr, IRanges(x$start, x$end))
  lifted <- unlist(liftOver(gr, chain_19_38), use.names=FALSE)
  if (!length(lifted)) return(data.table(chr=character(), start=integer(), end=integer()))
  reduce_bed(data.table(chr=as.character(seqnames(lifted)), start=start(lifted), end=end(lifted)))
}
write_bed <- function(x, path) {
  y <- copy(x); y[, start0 := pmax(start - 1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep="\t", col.names=FALSE)
}
safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

method_info <- data.table(
  method = c("epilink_nb", "epilink_poisson", "scmultimap",
             "signac_corrected", "signac_uncorrected"),
  method_label = c("EpiLink-NB", "EpiLink-Poisson", "scMultiMap",
                   "Signac (corrected)", "Signac (uncorrected)"),
  correction = c(FALSE, FALSE, NA, TRUE, FALSE)
)

rows <- list(); k <- 0L
for (i in seq_len(nrow(scope))) {
  ds <- scope$dataset[i]
  cells <- scope$celltype_vec[[i]]
  for (ct in cells) for (j in seq_len(nrow(method_info))) {
    mi <- method_info[j]
    id <- if (mi$method %chin% c("epilink_nb", "epilink_poisson"))
      paste(ds, ct, mi$method, "correction_F", sep="__") else
      paste(ds, ct, mi$method, "primary", sep="__")

    if (mi$method %chin% c("epilink_nb", "epilink_poisson")) {
      obj_path <- file.path(root, "results", "real", ds, ct, mi$method, "epilink_object.rds")
      obj <- readRDS(obj_path)
      sig <- EpiLinkFilter(obj, stat="zval", fdrcutoff=0.1,
                           correction=FALSE, conservative=FALSE, verbose=FALSE)
      native <- reduce_bed(parse_peaks(if (nrow(sig) && "peak" %in% names(sig)) sig$peak else character()))
      build <- dataset_manifest[dataset == ds & celltype == ct, build][1]
      bed <- if (identical(build, "hg19")) lift_bed(native) else native
      n_pairs <- nrow(sig); rm(obj); gc(verbose=FALSE)
    } else {
      src_id <- paste(safe_id(ds), safe_id(ct), mi$method, sep="__")
      z <- source_manifest[annotation_id == src_id]
      bed <- read_bed(file.path(gwas, "SCENT_style", "testing", "annotations", paste0(src_id, ".hg38.bed")))
      n_pairs <- if (nrow(z)) z$n_pairs[1] else NA_integer_
      build <- if (nrow(z)) z$source_build[1] else "hg38"
    }
    write_bed(bed, file.path(ann_dir, paste0(id, ".hg38.bed")))
    k <- k + 1L
    rows[[k]] <- data.table(
      annotation_id=id, dataset=ds, celltype=ct, celltype_scope=ct,
      n_celltypes=1L, celltypes=ct, method=mi$method,
      method_label=mi$method_label, correction=mi$correction, fdr=0.1,
      source_build=build, n_pairs=n_pairs, n_native_intervals=nrow(bed),
      n_hg38_intervals=nrow(bed), n_hg19_intervals=NA_integer_,
      lift_fraction_hg38=1, lift_fraction_hg19=NA_real_
    )
  }
}

manifest <- rbindlist(rows, fill=TRUE)
fwrite(manifest, file.path(base, "manifest.tsv"), sep="\t")
fwrite(scope[, .(dataset, celltypes, n_celltypes)], file.path(base, "dataset_scope.tsv"), sep="\t")
message("Wrote ", nrow(manifest), " cell-type-specific annotations to ", base)
