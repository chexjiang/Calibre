#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
K <- as.integer(if (length(args) >= 2) args[[2]] else 1000L)
tag <- paste0("K", K)
gwas <- file.path(root, "results", "validation", "GWAS")
out_dir <- file.path(gwas, "SCENT_style", "testing", "size_matched", tag)
ann_dir <- file.path(out_dir, "annotations")
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- fread(file.path(gwas, "SCENT_style", "config", "datasets.tsv"))
methods <- fread(file.path(gwas, "SCENT_style", "config", "methods.tsv"))
safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

method_input <- function(root, dataset, celltype, method) {
  signac_dir <- file.path(root, "results", "real", dataset, celltype, "signac")
  if (method == "signac_corrected") {
    return(file.path(signac_dir, "significant_pairs_corrected.tsv.gz"))
  }
  if (method == "signac_uncorrected") {
    return(file.path(signac_dir, "significant_pairs_uncorrected.tsv.gz"))
  }
  file.path(root, "results", "real", dataset, celltype, method,
            "significant_pairs.tsv.gz")
}

parse_peaks <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(data.table(chr=character(), start=integer(), end=integer()))
  z <- tstrsplit(x, "-", fixed=TRUE)
  if (length(z) != 3L) stop("Unexpected peak format")
  ans <- data.table(chr=z[[1]], start=as.integer(z[[2]]), end=as.integer(z[[3]]))
  unique(ans[grepl("^chr([1-9]|1[0-9]|2[0-2]|X)$", chr) & !is.na(start) & !is.na(end) & start <= end])
}

reduce_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand=TRUE)
  data.table(chr=as.character(seqnames(gr)), start=start(gr), end=end(gr))
}

lift_bed <- function(x, chain_path) {
  if (!nrow(x)) return(x)
  tmp <- tempfile(fileext=".chain")
  writeLines(readLines(gzfile(chain_path), warn=FALSE), tmp)
  chain <- import.chain(tmp)
  gr <- GRanges(x$chr, IRanges(x$start, x$end))
  y <- unlist(liftOver(gr, chain), use.names=FALSE)
  if (!length(y)) return(data.table(chr=character(), start=integer(), end=integer()))
  reduce_bed(data.table(chr=as.character(seqnames(y)), start=start(y), end=end(y)))
}

write_bed <- function(x, path) {
  y <- copy(x)
  y[, start0 := pmax(start - 1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep="\t", col.names=FALSE)
}

chain <- file.path(gwas, "reference", "liftover", "hg19ToHg38.over.chain.gz")
rows <- list(); k <- 0L
for (i in seq_len(nrow(datasets))) for (j in seq_len(nrow(methods))) {
  d <- datasets[i]; m <- methods[j]
  f <- method_input(root, d$dataset, d$celltype, m$method)
  if (!file.exists(f)) next
  cols <- names(fread(f, nrows=0L, showProgress=FALSE))
  use <- intersect(c("peak", "pval"), cols)
  z <- fread(f, select=use, showProgress=FALSE)
  if (!"pval" %in% names(z)) z[, pval := Inf]
  z[, pval := as.numeric(pval)]
  z[!is.finite(pval), pval := Inf]
  # Fixed K, identical to all methods: rank unique peaks by best pair p-value.
  z <- z[order(pval)][, .(pval=min(pval, na.rm=TRUE)), by=peak][order(pval)]
  z <- z[seq_len(min(K, nrow(z)))]
  native <- reduce_bed(parse_peaks(z$peak))
  hg38 <- if (d$build == "hg38") native else lift_bed(native, chain)
  id <- paste(safe_id(d$dataset), safe_id(d$celltype), safe_id(m$method), sep="__")
  write_bed(hg38, file.path(ann_dir, paste0(id, ".hg38.bed")))
  k <- k + 1L
  rows[[k]] <- data.table(
    annotation_id=id, dataset=d$dataset, celltype=d$celltype,
    panel=d$panel, interpretation=d$interpretation, method=m$method,
    method_label=m$label, source_build="hg38", n_pairs=K,
    n_native_intervals=nrow(native), n_hg38_intervals=nrow(hg38),
    n_hg19_intervals=NA_integer_, lift_fraction_hg38=1,
    lift_fraction_hg19=NA_real_, ranking="unique_peak_min_pval", top_k=K
  )
}
manifest <- rbindlist(rows, fill=TRUE)
fwrite(manifest, file.path(out_dir, "annotation_manifest.tsv"), sep="\t")
message("Wrote ", nrow(manifest), " size-matched annotations at K=", K)
