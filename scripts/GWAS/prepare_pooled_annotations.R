#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
config <- file.path(gwas, "SCENT_style", "config")
src <- fread(file.path(config, "annotation_manifest.tsv"))

out_scent <- file.path(gwas, "SCENT_style", "testing", "pooled", "annotations")
out_sldsc <- file.path(gwas, "sLDSC", "pooled", "annotations")
dir.create(out_scent, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sldsc, recursive = TRUE, showWarnings = FALSE)

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L)
    return(data.table(chr=character(), start=integer(), end=integer()))
  fread(path, header=FALSE, col.names=c("chr", "start0", "end"))[, .(
    chr, start=start0 + 1L, end
  )]
}

reduce_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand=TRUE)
  data.table(chr=as.character(seqnames(gr)), start=start(gr), end=end(gr))
}

write_bed <- function(x, path) {
  y <- copy(x)
  y[, start0 := pmax(start - 1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep="\t", col.names=FALSE)
}

safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
rows <- list(); k <- 0L

# Keep dataset and GWAS panel separate. This avoids mixing, for example,
# Arthritis immune cells with Arthritis stromal exploratory annotations.
groups <- unique(src[, .(dataset, panel)])
for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  for (m in unique(src$method)) {
    z <- src[dataset == g$dataset & panel == g$panel & method == m]
    if (!nrow(z)) next
    hg38 <- rbindlist(lapply(z$annotation_id, function(id)
      read_bed(file.path(gwas, "SCENT_style", "testing", "annotations", paste0(id, ".hg38.bed")))), fill=TRUE)
    hg19 <- rbindlist(lapply(z$annotation_id, function(id)
      read_bed(file.path(gwas, "sLDSC", "annotations", paste0(id, ".hg19.bed")))), fill=TRUE)
    hg38 <- reduce_bed(hg38)
    hg19 <- reduce_bed(hg19)
    id <- paste(safe_id(g$dataset), "ALL", safe_id(g$panel), safe_id(m), sep="__")
    write_bed(hg38, file.path(out_scent, paste0(id, ".hg38.bed")))
    write_bed(hg19, file.path(out_sldsc, paste0(id, ".hg19.bed")))
    k <- k + 1L
    rows[[k]] <- data.table(
      annotation_id=id, dataset=g$dataset, celltype="ALL",
      panel=g$panel, interpretation=if (all(z$interpretation == "primary")) "primary" else "mixed",
      method=m, method_label=z$method_label[1], source_build="hg38",
      n_source_celltypes=uniqueN(z$celltype), n_native_intervals=NA_integer_,
      n_hg38_intervals=nrow(hg38), n_hg19_intervals=nrow(hg19),
      n_source_annotations=nrow(z)
    )
  }
}

manifest <- rbindlist(rows, fill=TRUE)
fwrite(manifest, file.path(gwas, "SCENT_style", "testing", "pooled", "annotation_manifest.tsv"), sep="\t")
fwrite(manifest, file.path(gwas, "sLDSC", "pooled", "annotation_manifest.tsv"), sep="\t")
message("Wrote ", nrow(manifest), " pooled annotations across ", uniqueN(manifest$dataset), " datasets.")
