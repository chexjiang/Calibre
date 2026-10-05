#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)

gwas <- file.path(root, "results", "validation", "GWAS")
config <- file.path(gwas, "SCENT_style", "config")
datasets <- fread(file.path(config, "datasets.tsv"))
methods <- fread(file.path(config, "methods.tsv"))

import_chain_gz <- function(path) {
  tmp <- tempfile(fileext = ".over.chain")
  in_con <- gzfile(path, open = "rt")
  lines <- readLines(in_con, warn = FALSE)
  close(in_con)
  writeLines(lines, tmp)
  import.chain(tmp)
}

chain_19_38 <- import_chain_gz(file.path(gwas, "reference", "liftover", "hg19ToHg38.over.chain.gz"))
chain_38_19 <- import_chain_gz(file.path(gwas, "sLDSC", "reference", "raw", "hg38ToHg19.over.chain.gz"))

safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

signac_input <- function(root, dataset, celltype, method) {
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
  if (!length(x)) return(data.table(chr = character(), start = integer(), end = integer()))
  z <- tstrsplit(x, "-", fixed = TRUE)
  if (length(z) != 3L) stop("Peak IDs must have chr-start-end form")
  ans <- data.table(chr = z[[1]], start = as.integer(z[[2]]), end = as.integer(z[[3]]))
  ans <- ans[grepl("^chr([1-9]|1[0-9]|2[0-2]|X)$", chr) & !is.na(start) & !is.na(end) & start <= end]
  unique(ans)
}

reduce_bed <- function(x) {
  if (!nrow(x)) return(data.table(chr = character(), start = integer(), end = integer()))
  gr <- GRanges(x$chr, IRanges(x$start, x$end))
  gr <- reduce(gr, ignore.strand = TRUE)
  data.table(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
}

lift_bed <- function(x, chain) {
  if (!nrow(x)) return(x)
  gr <- GRanges(x$chr, IRanges(x$start, x$end))
  lifted <- unlist(liftOver(gr, chain), use.names = FALSE)
  if (!length(lifted)) return(data.table(chr = character(), start = integer(), end = integer()))
  reduce_bed(data.table(chr = as.character(seqnames(lifted)), start = start(lifted), end = end(lifted)))
}

write_bed <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # BED is zero-based, half-open.
  y <- copy(x)
  if (nrow(y)) y[, start0 := pmax(start - 1L, 0L)] else y[, start0 := integer()]
  fwrite(y[, .(chr, start0, end)], path, sep = "\t", col.names = FALSE)
}

manifest <- list()
k <- 0L
for (i in seq_len(nrow(datasets))) {
  d <- datasets[i]
  for (j in seq_len(nrow(methods))) {
    m <- methods[j]
    input <- signac_input(root, d$dataset, d$celltype, m$method)
    if (!file.exists(input)) next
    dt <- tryCatch(fread(input, select = "peak", showProgress = FALSE), error = function(e) data.table(peak = character()))
    native <- reduce_bed(parse_peaks(unique(dt$peak)))
    if (d$build == "hg38") {
      hg38 <- native
      hg19 <- lift_bed(native, chain_38_19)
    } else {
      hg19 <- native
      hg38 <- lift_bed(native, chain_19_38)
    }
    id <- paste(safe_id(d$dataset), safe_id(d$celltype), safe_id(m$method), sep = "__")
    write_bed(hg38, file.path(gwas, "SCENT_style", "testing", "annotations", paste0(id, ".hg38.bed")))
    write_bed(hg19, file.path(gwas, "sLDSC", "annotations", paste0(id, ".hg19.bed")))
    k <- k + 1L
    manifest[[k]] <- data.table(
      annotation_id = id, dataset = d$dataset, celltype = d$celltype,
      panel = d$panel, interpretation = d$interpretation,
      method = m$method, method_label = m$label, source_build = d$build,
      n_pairs = nrow(dt), n_native_intervals = nrow(native),
      n_hg38_intervals = nrow(hg38), n_hg19_intervals = nrow(hg19),
      lift_fraction_hg38 = if (nrow(native)) nrow(hg38) / nrow(native) else NA_real_,
      lift_fraction_hg19 = if (nrow(native)) nrow(hg19) / nrow(native) else NA_real_
    )
  }
}

manifest <- rbindlist(manifest, fill = TRUE)
fwrite(manifest, file.path(config, "annotation_manifest.tsv"), sep = "\t")
message("Wrote ", nrow(manifest), " annotations to both hg38 (SCENT) and hg19 (sLDSC).")
