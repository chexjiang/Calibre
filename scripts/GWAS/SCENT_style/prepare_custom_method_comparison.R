#!/usr/bin/env Rscript

# Pool the cell-type-specific annotations created by
# prepare_celltype_comparison.R using a user-specified dataset scope.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork=TRUE)
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
source_results <- source_results_paths(root)
analysis_arg <- if (length(args) >= 2L) args[[2]] else "pooled"
scope_path <- if (length(args) >= 3L) normalizePath(args[[3]], mustWork=TRUE) else
  file.path(source_results$shared, "custom_scope.tsv")
celltype_arg <- if (length(args) >= 4L) args[[4]] else "celltype"
gwas <- file.path(root, "results", "validation", "GWAS")
sakaue <- source_results$root
analysis_name <- if (analysis_arg %in% c("method_comparison", "method_comparison_custom", "pooled_main", "pooled"))
  "pooled" else analysis_arg
celltype_name <- if (celltype_arg %in% c("celltype_comparison", "celltype_comparison_custom", "celltype_main", "celltype"))
  "celltype" else celltype_arg
base <- if (analysis_name == "pooled") source_results$pooled else
  file.path(sakaue, "pooled", analysis_name)
src_base <- if (celltype_name == "celltype") source_results$celltype else
  file.path(sakaue, "celltype", celltype_name)
ann_dir <- file.path(base, "annotations")
dir.create(ann_dir, recursive=TRUE, showWarnings=FALSE)

scope <- fread(scope_path)
scope[, celltype_vec := strsplit(celltypes, ";", fixed=TRUE)]
src_manifest <- fread(file.path(src_base, "manifest.tsv"))

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L)
    return(data.table(chr=character(), start=integer(), end=integer()))
  x <- fread(path, header=FALSE, col.names=c("chr", "start0", "end"))
  x[, .(chr, start=as.integer(start0)+1L, end=as.integer(end))]
}
reduce_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand=TRUE)
  data.table(chr=as.character(seqnames(gr)), start=start(gr), end=end(gr))
}
write_bed <- function(x, path) {
  y <- copy(x); y[, start0 := pmax(start-1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep="\t", col.names=FALSE)
}

method_info <- data.table(
  method=c("epilink_nb", "epilink_poisson", "scmultimap", "signac_corrected", "signac_uncorrected"),
  method_label=c("EpiLink-NB (correction = FALSE)",
                 "EpiLink-Poisson (correction = FALSE)", "scMultiMap",
                 "Signac (corrected)", "Signac (uncorrected)"),
  correction=c(FALSE, FALSE, NA, TRUE, FALSE)
)

rows <- list(); k <- 0L
for (i in seq_len(nrow(scope))) {
  ds <- scope$dataset[i]; cells <- scope$celltype_vec[[i]]
  for (j in seq_len(nrow(method_info))) {
    mi <- method_info[j]
    z <- src_manifest[dataset == ds & celltype %chin% cells & method == mi$method]
    if (!nrow(z)) stop("Missing source annotations for ", ds, " / ", mi$method)
    bed <- rbindlist(lapply(z$annotation_id, function(id)
      read_bed(file.path(src_base, "annotations", paste0(id, ".hg38.bed")))), fill=TRUE)
    bed <- reduce_bed(bed)
    id <- if (mi$method %chin% c("epilink_nb", "epilink_poisson"))
      paste(ds, mi$method, "correction_F", sep="__") else
      paste(ds, mi$method, "primary", sep="__")
    write_bed(bed, file.path(ann_dir, paste0(id, ".hg38.bed")))
    k <- k + 1L
    rows[[k]] <- data.table(
      annotation_id=id, dataset=ds, celltype="ALL_REQUESTED",
      celltype_scope="ALL_REQUESTED", n_celltypes=length(cells),
      celltypes=paste(cells, collapse=";"), method=mi$method,
      method_label=mi$method_label, correction=mi$correction, fdr=0.1,
      source_build="hg38", n_pairs=sum(z$n_pairs, na.rm=TRUE),
      n_native_intervals=sum(z$n_native_intervals, na.rm=TRUE),
      n_hg38_intervals=nrow(bed), n_hg19_intervals=NA_integer_,
      lift_fraction_hg38=1, lift_fraction_hg19=NA_real_
    )
  }
}

manifest <- rbindlist(rows, fill=TRUE)
fwrite(manifest, file.path(base, "manifest.tsv"), sep="\t")
fwrite(scope[, .(dataset, celltypes, n_celltypes, note)], file.path(base, "dataset_scope.tsv"), sep="\t")
message("Wrote ", nrow(manifest), " pooled annotations to ", base)
